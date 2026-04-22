//! All-pairs cosine similarity + alignment verification clustering.
//!
//! Links pairs of sequences whose cosine similarity exceeds a threshold,
//! verifies with alignment identity, and forms connected components via
//! a lock-free UnionFind.  The longest member of each component is the
//! representative.

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::align::DpAligner;
use crate::utils::human::CommaReadable;
use crate::utils::union_find_2::UnionFind;

use super::{alignment_stats, cosine_similarity, GroupStats};

/// Run all-pairs cosine + alignment clustering on a single group.
///
/// `per_member` is indexed by group-local index and contains
/// `(freq_vector, l2_norm, sequence)` for each member.
///
/// `genomic_lens` gives the genomic length of each member (for selecting
/// the longest representative per cluster).
///
/// Returns the list of group-local indices chosen as representatives,
/// plus accumulated statistics.
pub fn cluster_group(
    per_member: &[(Vec<u8>, f32, Vec<u8>)],
    genomic_lens: &[usize],
    group_label: &str,
    min_cosine: f32,
    min_identity: f64,
    threads: usize,
    collect_stats: bool,
) -> (Vec<usize>, GroupStats) {
    let group_n = per_member.len();
    let uf = UnionFind::new(group_n);
    let mut stats = GroupStats::default();
    if collect_stats {
        stats.cosine_hist = vec![0usize; 21];
    }

    if group_n >= 2 {
        // Filter to members with non-zero frequency vectors
        let active: Vec<usize> = (0..group_n)
            .filter(|&i| per_member[i].1 > 0.0)
            .collect();

        let num_pairs = active.len() * (active.len().saturating_sub(1)) / 2;
        if num_pairs > 1_000_000 {
            log::warn!(
                "    {group_label}: {} active instances \u{2192} {} pairs (large group!)",
                active.len().commas(),
                num_pairs.commas()
            );
        }

        // Parallel all-pairs: split the triangular pair space into
        // equal-sized blocks across threads.  Each thread gets its
        // own Aligner and local counters; the lock-free UnionFind is
        // shared by reference.
        let progress = AtomicUsize::new(0);
        let percentage_done = AtomicUsize::new(0);

        struct ThreadResult {
            candidates: usize,
            linked: usize,
            cosine_pass: usize,
            already_merged: usize,
            aligned: usize,
            cosine_hist: Vec<usize>,
        }

        let num_threads = threads.min(num_pairs.max(1));

        let results: Vec<ThreadResult> = std::thread::scope(|s| {
            let handles: Vec<_> = (0..num_threads)
                .map(|tid| {
                    let uf = &uf;
                    let active = &active;
                    let per_member = per_member;
                    let progress = &progress;
                    let percentage_done = &percentage_done;
                    s.spawn(move || {
                        let mut aligner = DpAligner::with_defaults();
                        let mut local = ThreadResult {
                            candidates: 0,
                            linked: 0,
                            cosine_pass: 0,
                            already_merged: 0,
                            aligned: 0,
                            cosine_hist: if collect_stats {
                                vec![0usize; 21]
                            } else {
                                Vec::new()
                            },
                        };

                        // Divide pairs: thread tid handles linear indices
                        // [tid * block_size .. (tid+1) * block_size)
                        let block_size =
                            (num_pairs + num_threads - 1) / num_threads;
                        let pair_start = tid * block_size;
                        let pair_end = ((tid + 1) * block_size).min(num_pairs);

                        // Convert linear pair index → (ai, bi) in the
                        // upper triangle.
                        let na = active.len();
                        let mut ai = 0usize;
                        let mut row_start = 0usize;

                        // Advance ai to the row containing pair_start
                        while ai < na {
                            let row_len = na - 1 - ai;
                            if row_len == 0 {
                                break;
                            }
                            if row_start + row_len > pair_start {
                                break;
                            }
                            row_start += row_len;
                            ai += 1;
                        }

                        let mut p = pair_start;
                        while p < pair_end && ai < na.saturating_sub(1) {
                            let row_len = na - 1 - ai;
                            let offset_in_row = p - row_start;
                            let bi_start = ai + 1 + offset_in_row;
                            let bi_end_this_row =
                                (ai + 1 + row_len).min(ai + 1 + (pair_end - row_start).min(row_len));

                            let a = active[ai];
                            let (ref freq_a, norm_a, _) = per_member[a];

                            for bi in bi_start..bi_end_this_row {
                                let b = active[bi];
                                let (ref freq_b, norm_b, _) = per_member[b];
                                local.candidates += 1;

                                let c = cosine_similarity(
                                    freq_a, freq_b, norm_a, norm_b,
                                );
                                if collect_stats {
                                    let bucket =
                                        ((c * 20.0) as usize).min(20);
                                    local.cosine_hist[bucket] += 1;
                                }
                                if c >= min_cosine {
                                    local.cosine_pass += 1;
                                    if uf.find(a) == uf.find(b) {
                                        local.already_merged += 1;
                                    } else {
                                        local.aligned += 1;
                                        let seq_a = &per_member[a].2;
                                        let seq_b = &per_member[b].2;
                                        if let Some(aln) =
                                            aligner.align(seq_a, seq_b)
                                        {
                                            let (identity, _) =
                                                alignment_stats(&aln.cigar);
                                            if identity >= min_identity {
                                                local.linked += 1;
                                                uf.union(a, b);
                                            }
                                        }
                                    }
                                }
                            }

                            let done_this_batch = bi_end_this_row - bi_start;
                            p += done_this_batch;

                            // Progress reporting (approximate, via atomic)
                            if num_pairs > 0 {
                                let prev = progress.fetch_add(
                                    done_this_batch,
                                    Ordering::Relaxed,
                                );
                                let new_pct =
                                    ((prev + done_this_batch) * 100) / num_pairs;
                                let old_pct = percentage_done
                                    .load(Ordering::Relaxed);
                                if new_pct > old_pct {
                                    let _ = percentage_done
                                        .compare_exchange(
                                            old_pct,
                                            new_pct,
                                            Ordering::Relaxed,
                                            Ordering::Relaxed,
                                        );
                                    log::info!(
                                        "    {group_label}\t{new_pct}%"
                                    );
                                }
                            }

                            // Advance to next row
                            row_start += row_len;
                            ai += 1;
                        }

                        local
                    })
                })
                .collect();

            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });

        // Merge per-thread results
        for r in &results {
            stats.candidates += r.candidates;
            stats.linked += r.linked;
            stats.cosine_pass += r.cosine_pass;
            stats.already_merged += r.already_merged;
            stats.aligned += r.aligned;
            if collect_stats {
                for (i, &c) in r.cosine_hist.iter().enumerate() {
                    stats.cosine_hist[i] += c;
                }
            }
        }

        log::info!(
            "    {} pairs, {} merged",
            stats.candidates.commas(),
            stats.linked.commas()
        );
    }

    // Select best (longest) representative per cluster
    let mut cluster_best: HashMap<usize, (usize, usize)> = HashMap::new();
    let mut cluster_size: HashMap<usize, usize> = HashMap::new();

    for local in 0..group_n {
        let root = uf.find(local);
        let len = genomic_lens[local];
        *cluster_size.entry(root).or_insert(0) += 1;
        let entry = cluster_best.entry(root).or_insert((local, 0));
        if len > entry.1 {
            *entry = (local, len);
        }
    }

    stats.clusters = cluster_best.len();
    stats.singletons = cluster_size.values().filter(|&&sz| sz == 1).count();
    stats.cluster_sizes = cluster_size.into_values().collect();

    let rep_locals: Vec<usize> = cluster_best.values().map(|&(l, _)| l).collect();
    (rep_locals, stats)
}
