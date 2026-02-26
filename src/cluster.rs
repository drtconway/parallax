//! Cluster genomic repeat element instances and produce representative sequences.
//!
//! This subcommand reads a BED file of repeat element annotations
//! (produced by `scripts/ucsc_repeats_to_bed.py`) together with a reference
//! FASTA, extracts the genomic sequence of each instance (reverse-complemented
//! for minus-strand entries), and clusters instances within each
//! `(repClass, repFamily, repName)` group by sequence similarity using
//! dense k-mer frequency vectors and cosine similarity.
//!
//! For each cluster the longest (most complete) instance is written to the
//! output FASTA so that the resulting catalogue can be used directly by the
//! `annotate` subcommand.
//!
//! ## Algorithm
//!
//! 1. For each strand-resolved instance, build a dense `u8` frequency vector
//!    over the 4^K k-mer space (K=5 → 1024 entries, 1 KiB per instance).
//! 2. Precompute L2 norms for each vector.
//! 3. Within each `(repClass, repFamily, repName)` group, compute all-pairs
//!    cosine similarity and link pairs above `--min-cosine`.
//! 4. Form connected-component clusters via UnionFind; emit the longest
//!    member of each cluster.

#![allow(dead_code)]

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use clap::Args;

use crate::align::{Aligner, Kind, Op};
use crate::error::Result;
use crate::kmers::Kmer;
use crate::reference::InMemoryReference;
use crate::utils::human::CommaReadable;
use crate::utils::sequence::reverse_complement_into;
use crate::utils::union_find_2::UnionFind;

// ---------------------------------------------------------------------------
// K-mer / MinHash parameters
// ---------------------------------------------------------------------------

/// K-mer length for frequency vectors.  Smaller K produces a denser,
/// more informative vector at the cost of k-mer specificity.  K=5
/// gives 4^5 = 1024 distinct k-mers – compact enough for u8 vectors
/// (1 KiB each) while still discriminating repeat families well.
const K: usize = 5;

/// Total number of distinct k-mers: 4^K.
const NUM_KMERS: usize = 1 << (2 * K);

const DEFAULT_MIN_COSINE: f32 = 0.75;
const DEFAULT_MIN_IDENTITY: f64 = 0.85;

// ---------------------------------------------------------------------------
// CLI args
// ---------------------------------------------------------------------------

/// Cluster repeat element instances and produce representative sequences.
#[derive(Args, Debug, Clone)]
pub struct ClusterArgs {
    /// Reference FASTA file (plain or compressed; must be readable by niffler)
    pub reference: PathBuf,

    /// Input BED file produced by scripts/ucsc_repeats_to_bed.py
    /// (columns: chrom, start, end, repName|repClass|repFamily, score, strand)
    pub bed: PathBuf,

    /// Output FASTA file (default: stdout)
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Number of threads
    #[arg(short = 't', long, default_value = "4")]
    pub threads: usize,

    /// Minimum cosine similarity to consider a pair for alignment verification
    #[arg(long, default_value_t = DEFAULT_MIN_COSINE)]
    pub min_cosine: f32,

    /// Minimum alignment identity fraction to merge two instances
    #[arg(long, default_value_t = DEFAULT_MIN_IDENTITY)]
    pub min_identity: f64,

    /// Emit diagnostic histograms (cosine, cluster-size, inter-cluster identity)
    #[arg(long)]
    pub stats: bool,
}

// ---------------------------------------------------------------------------
// BED record
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct BedRecord {
    chrom: String,
    start: usize, // 0-based, half-open
    end: usize,
    rep_name: String,
    rep_class: String,
    rep_family: String,
    strand: char,
}

impl BedRecord {
    fn genomic_len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    fn group_key(&self) -> (&str, &str, &str) {
        (&self.rep_class, &self.rep_family, &self.rep_name)
    }
}

fn parse_bed(path: &PathBuf) -> std::io::Result<Vec<BedRecord>> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.splitn(7, '\t').collect();
        if fields.len() < 6 {
            log::warn!(
                "BED line {}: expected ≥6 fields, got {}, skipping",
                lineno + 1,
                fields.len()
            );
            continue;
        }
        let chrom = fields[0].to_string();
        let start: usize = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => {
                log::warn!("BED line {}: bad start field, skipping", lineno + 1);
                continue;
            }
        };
        let end: usize = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => {
                log::warn!("BED line {}: bad end field, skipping", lineno + 1);
                continue;
            }
        };
        let name = fields[3];
        let strand = fields[5].chars().next().unwrap_or('+');

        // Name format from ucsc_repeats_to_bed.py: "repName|repClass|repFamily"
        let mut parts = name.splitn(3, '|');
        let rep_name = match parts.next() {
            Some(s) => s.to_string(),
            None => {
                log::warn!(
                    "BED line {}: could not parse name '{}', skipping",
                    lineno + 1,
                    name
                );
                continue;
            }
        };
        let rep_class = match parts.next() {
            Some(s) => s.to_string(),
            None => {
                log::warn!(
                    "BED line {}: could not parse name '{}', skipping",
                    lineno + 1,
                    name
                );
                continue;
            }
        };
        let rep_family = match parts.next() {
            Some(s) => s.to_string(),
            None => {
                log::warn!(
                    "BED line {}: could not parse name '{}', skipping",
                    lineno + 1,
                    name
                );
                continue;
            }
        };

        records.push(BedRecord {
            chrom,
            start,
            end,
            rep_name,
            rep_class,
            rep_family,
            strand,
        });
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// K-mer frequency vector
// ---------------------------------------------------------------------------

/// Build a dense frequency vector over the 4^K k-mer space.
///
/// Counts are accumulated in a caller-supplied `u16` buffer (one per
/// worker thread) and then scaled to `u8`:
///   - if max ≤ 255: direct truncation (lossless);
///   - otherwise: linear rescale so that max maps to 255.
fn kmer_freq_vector(seq: &[u8], counts: &mut [u16; NUM_KMERS]) -> Vec<u8> {
    counts.iter_mut().for_each(|c| *c = 0);
    Kmer::<K>::kmerize_fwd(seq, |_pos, kmer| {
        let idx = kmer.0 as usize;
        counts[idx] = counts[idx].saturating_add(1);
    });
    let max_count = counts.iter().copied().max().unwrap_or(0);
    if max_count == 0 {
        return vec![0u8; NUM_KMERS];
    }
    if max_count <= 255 {
        counts.iter().map(|&c| c as u8).collect()
    } else {
        let scale = max_count as f32 / 255.0;
        counts.iter().map(|&c| (c as f32 / scale) as u8).collect()
    }
}

// ---------------------------------------------------------------------------
// Cosine similarity on dense u8 frequency vectors
// ---------------------------------------------------------------------------

/// Euclidean (L2) norm of a `u8` frequency vector.
///
/// The inner loop accumulates into `u32` so that 1024 × 255² = 66 million
/// fits comfortably.  The simple sum-of-squares loop auto-vectorises well
/// on x86-64 and aarch64.
#[inline]
fn l2_norm(v: &[u8]) -> f32 {
    let sum_sq: u32 = v.iter().map(|&x| (x as u32) * (x as u32)).sum();
    (sum_sq as f32).sqrt()
}

/// Cosine similarity between two dense `u8` frequency vectors whose L2
/// norms have been precomputed.
///
/// Returns 0.0 when either vector is the zero vector.  The dot-product
/// loop is written to mirror `l2_norm` for auto-vectorisation.
#[inline]
fn cosine_similarity(a: &[u8], b: &[u8], norm_a: f32, norm_b: f32) -> f32 {
    debug_assert_eq!(a.len(), b.len());
    let denom = norm_a * norm_b;
    if denom == 0.0 {
        return 0.0;
    }
    let dot: u32 = a.iter().zip(b.iter()).map(|(&x, &y)| x as u32 * y as u32).sum();
    dot as f32 / denom
}

// ---------------------------------------------------------------------------
// Alignment identity & edit distance from extended CIGAR
// ---------------------------------------------------------------------------

/// Compute identity fraction and edit distance from an extended CIGAR.
///
/// Identity = (sequence-match columns) / (total alignment columns).
/// Edit distance = mismatches + insertions + deletions (by column count).
fn alignment_stats(cigar: &[Op]) -> (f64, usize) {
    let mut matches = 0usize;
    let mut columns = 0usize;
    let mut edits = 0usize;
    for &op in cigar {
        let n = op.len();
        match op.kind() {
            Kind::SequenceMatch => {
                matches += n;
                columns += n;
            }
            Kind::SequenceMismatch => {
                edits += n;
                columns += n;
            }
            Kind::Insertion | Kind::Deletion => {
                edits += n;
                columns += n;
            }
            Kind::Match => {
                // Ambiguous M — count as column but not match or edit
                columns += n;
            }
            _ => {}
        }
    }
    let identity = if columns == 0 {
        0.0
    } else {
        matches as f64 / columns as f64
    };
    (identity, edits)
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

pub fn run(args: ClusterArgs) -> Result<()> {
    // 1. Load reference
    log::info!("Loading reference from {} ...", args.reference.display());
    let reference = Arc::new(InMemoryReference::load(&args.reference, false)?);
    let chrom_index: HashMap<String, usize> = (0..reference.num_chroms())
        .map(|i| (reference.chrom_name(i).to_string(), i))
        .collect();

    // 2. Parse BED and group by (repClass, repFamily, repName)
    log::info!("Parsing BED from {} ...", args.bed.display());
    let records = parse_bed(&args.bed)?;
    let n = records.len();
    log::info!("  {} records", n.commas());

    log::info!("Grouping into element families ...");
    let mut groups: HashMap<(String, String, String), Vec<usize>> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        groups
            .entry((
                rec.rep_class.clone(),
                rec.rep_family.clone(),
                rec.rep_name.clone(),
            ))
            .or_default()
            .push(i);
    }
    let mut sorted_keys: Vec<_> = groups.keys().cloned().collect();
    sorted_keys.sort();
    let num_groups = sorted_keys.len();
    log::info!("  {} groups", num_groups.commas());

    // 3. Open output
    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    // Global accumulators
    let mut aligner = Aligner::with_defaults();
    let mut total_candidates = 0usize;
    let mut total_cosine_pass = 0usize;
    let mut total_already_merged = 0usize;
    let mut total_aligned = 0usize;
    let mut total_linked = 0usize;
    let mut total_clusters = 0usize;
    let mut total_singletons = 0usize;
    let mut written = 0usize;
    let mut cosine_hist = if args.stats { vec![0usize; 21] } else { Vec::new() };
    let mut all_cluster_sizes: Vec<usize> = Vec::new();
    // For inter-cluster identity we collect (group_key, rep seqs) to process at the end
    let mut inter_cluster_seqs: Vec<Vec<Vec<u8>>> = Vec::new();

    // 4. Process each group independently
    log::info!(
        "Processing groups (min_cosine={:.2}, min_identity={:.2}, {} threads) ...",
        args.min_cosine,
        args.min_identity,
        args.threads,
    );

    for (group_idx, key) in sorted_keys.iter().enumerate() {
        let indices = &groups[key];
        let (rc, rf, rn) = key;
        let group_n = indices.len();

        log::info!(
            "  [{}/{num_groups}] {rc}/{rf}/{rn}: {} instances",
            group_idx + 1,
            group_n.commas()
        );

        // ── Extract sequences & build frequency vectors (parallel) ──────
        // Map group-local index (0..group_n) ↔ global record index.
        // Work items: (local_idx, chrom_idx, start, end, strand)
        let (tx_work, rx_work) =
            crossbeam::channel::bounded::<(usize, usize, usize, usize, u8)>(args.threads * 8);
        let (tx_result, rx_result) =
            crossbeam::channel::unbounded::<(usize, Vec<u8>, f32, Vec<u8>)>();

        let ref_arc = Arc::clone(&reference);
        let workers: Vec<_> = (0..args.threads)
            .map(|_| {
                let rx = rx_work.clone();
                let tx = tx_result.clone();
                let reference = Arc::clone(&ref_arc);
                std::thread::spawn(move || {
                    let mut rc_buf = Vec::new();
                    let mut counts = [0u16; NUM_KMERS];
                    for (local_idx, chrom_idx, start, end, strand) in rx {
                        let raw = reference.get_seq(chrom_idx, start, end);
                        if raw.iter().any(|&b| b == b'N' || b == b'n') {
                            continue;
                        }
                        let (freq, seq) = if strand == b'-' {
                            reverse_complement_into(raw, &mut rc_buf);
                            (kmer_freq_vector(&rc_buf, &mut counts), rc_buf.clone())
                        } else {
                            (kmer_freq_vector(raw, &mut counts), raw.to_vec())
                        };
                        let norm = l2_norm(&freq);
                        tx.send((local_idx, freq, norm, seq)).ok();
                    }
                })
            })
            .collect();

        // Adjust channel type: we send tuples with resolved data
        drop(tx_result);
        drop(rx_work);

        for (local_idx, &global_idx) in indices.iter().enumerate() {
            let rec = &records[global_idx];
            if let Some(&ci) = chrom_index.get(&rec.chrom) {
                tx_work
                    .send((local_idx, ci, rec.start, rec.end, rec.strand as u8))
                    .unwrap();
            } else {
                log::warn!(
                    "Chromosome '{}' not found in reference; skipping",
                    rec.chrom
                );
            }
        }
        drop(tx_work);

        // Collect into group-local arrays
        let mut per_member: Vec<(Vec<u8>, f32, Vec<u8>)> = (0..group_n)
            .map(|_| (vec![0u8; NUM_KMERS], 0.0f32, Vec::new()))
            .collect();
        for (local_idx, freq, norm, seq) in rx_result {
            per_member[local_idx] = (freq, norm, seq);
        }
        for w in workers {
            w.join().unwrap();
        }

        // ── Clustering ──────────────────────────────────────────────────
        let uf = UnionFind::new(group_n);

        if group_n >= 2 {
            // Norm diagnostics
            if log::log_enabled!(log::Level::Debug) {
                let norms: Vec<f32> = per_member.iter().map(|(_, n, _)| *n).collect();
                let non_zero: Vec<f32> = norms.iter().copied().filter(|&n| n > 0.0).collect();
                let mean = if non_zero.is_empty() {
                    0.0
                } else {
                    non_zero.iter().sum::<f32>() / non_zero.len() as f32
                };
                log::debug!(
                    "Group {rc}/{rf}/{rn}: {} instances, {}/{} non-zero norms, mean_norm={mean:.1}",
                    group_n.commas(),
                    non_zero.len().commas(),
                    group_n.commas(),
                );
            }

            // Filter to members with non-zero frequency vectors
            let active: Vec<usize> = (0..group_n)
                .filter(|&i| per_member[i].1 > 0.0)
                .collect();

            let num_pairs = active.len() * (active.len().saturating_sub(1)) / 2;
            if num_pairs > 1_000_000 {
                log::warn!(
                    "    {rc}/{rf}/{rn}: {} active instances \u{2192} {} pairs (large group!)",
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

            let num_threads = args.threads.min(num_pairs.max(1));
            let collect_stats = args.stats;
            let min_cosine = args.min_cosine;
            let min_identity = args.min_identity;

            let results: Vec<ThreadResult> = std::thread::scope(|s| {
                let handles: Vec<_> = (0..num_threads)
                    .map(|tid| {
                        let uf = &uf;
                        let active = &active;
                        let per_member = &per_member;
                        let progress = &progress;
                        let percentage_done = &percentage_done;
                        s.spawn(move || {
                            let mut aligner = Aligner::with_defaults();
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
                            // upper triangle.  For pair index p, find ai such
                            // that sum_{k=0}^{ai-1} (n-1-k) <= p.
                            let na = active.len();
                            let mut ai = 0usize;
                            let mut row_start = 0usize; // first pair index in row ai

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
                                        // Best-effort update; ok if another
                                        // thread beats us.
                                        let _ = percentage_done
                                            .compare_exchange(
                                                old_pct,
                                                new_pct,
                                                Ordering::Relaxed,
                                                Ordering::Relaxed,
                                            );
                                        log::info!(
                                            "    {rc}/{rf}/{rn}\t{new_pct}%"
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
            let mut group_candidates = 0usize;
            let mut group_linked = 0usize;
            for r in &results {
                group_candidates += r.candidates;
                group_linked += r.linked;
                total_cosine_pass += r.cosine_pass;
                total_already_merged += r.already_merged;
                total_aligned += r.aligned;
                if args.stats {
                    for (i, &c) in r.cosine_hist.iter().enumerate() {
                        cosine_hist[i] += c;
                    }
                }
            }

            log::info!("    {} pairs, {} merged", group_candidates.commas(), group_linked.commas());
            total_candidates += group_candidates;
            total_linked += group_linked;
        }

        // ── Select best representative per cluster ──────────────────────
        let mut cluster_best: HashMap<usize, (usize, usize)> = HashMap::new();
        let mut cluster_size: HashMap<usize, usize> = HashMap::new();

        for local in 0..group_n {
            let root = uf.find(local);
            let len = records[indices[local]].genomic_len();
            *cluster_size.entry(root).or_insert(0) += 1;
            let entry = cluster_best.entry(root).or_insert((local, 0));
            if len > entry.1 {
                *entry = (local, len);
            }
        }

        let group_clusters = cluster_best.len();
        let group_singletons = cluster_size.values().filter(|&&sz| sz == 1).count();
        total_clusters += group_clusters;
        total_singletons += group_singletons;
        all_cluster_sizes.extend(cluster_size.values().copied());

        // ── Write representatives for this group (sorted by chrom, start) ──
        let mut rep_locals: Vec<usize> = cluster_best.values().map(|&(l, _)| l).collect();
        rep_locals.sort_unstable_by(|&a, &b| {
            let ra = &records[indices[a]];
            let rb = &records[indices[b]];
            ra.chrom.cmp(&rb.chrom).then(ra.start.cmp(&rb.start))
        });

        // Collect representative sequences for inter-cluster stats
        if args.stats && rep_locals.len() >= 2 {
            let seqs: Vec<Vec<u8>> = rep_locals
                .iter()
                .filter_map(|&l| {
                    let seq = &per_member[l].2;
                    if seq.is_empty() { None } else { Some(seq.clone()) }
                })
                .collect();
            if seqs.len() >= 2 {
                inter_cluster_seqs.push(seqs);
            }
        }

        for &local in &rep_locals {
            let rec = &records[indices[local]];
            let seq = &per_member[local].2;
            if seq.is_empty() {
                continue;
            }

            writeln!(
                out,
                ">{}|{}|{} {}:{}-{}({})",
                rec.rep_class, rec.rep_family, rec.rep_name,
                rec.chrom, rec.start, rec.end, rec.strand
            )?;
            for chunk in seq.chunks(60) {
                out.write_all(chunk)?;
                writeln!(out)?;
            }
            written += 1;
        }
        log::info!("{rc}/{rf}/{rn} wrote {} representatives", rep_locals.len().commas());

        // per_member is dropped here, freeing this group's memory
    }

    // 5. Summary statistics
    log::info!(
        "  {} pairs checked, {} cosine pass, {} already merged, {} aligned, {} linked (cosine \u{2265} {:.2}, identity \u{2265} {:.2})",
        total_candidates.commas(),
        total_cosine_pass.commas(),
        total_already_merged.commas(),
        total_aligned.commas(),
        total_linked.commas(),
        args.min_cosine,
        args.min_identity
    );
    log::info!("  {} clusters ({} singletons)", total_clusters.commas(), total_singletons.commas());

    if args.stats {
        // Cosine histogram
        log::info!("Pairwise cosine distribution:");
        log::info!("  {:>12}  {:>8}  {}", "range", "count", "bar");
        for (i, &count) in cosine_hist.iter().enumerate() {
            let lo = i as f32 * 0.05;
            let hi = if i == 20 { 1.0f32 } else { lo + 0.05 };
            let bar_len = if total_candidates > 0 {
                ((count as f64 / total_candidates as f64) * 200.0).ceil() as usize
            } else {
                0
            };
            let bar = "#".repeat(bar_len);
            log::info!("  [{lo:.2}, {hi:.2})  {count:>8}  {bar}");
        }

        // Cluster size histogram (log2-scale buckets)
        all_cluster_sizes.sort_unstable();
        let mut size_hist: Vec<(String, usize)> = Vec::new();
        let mut lo = 1usize;
        loop {
            let hi = if lo <= 2 { lo } else { (lo * 2) - 1 };
            let count = all_cluster_sizes.iter().filter(|&&s| s >= lo && s <= hi).count();
            if count > 0 {
                let label = if lo == hi {
                    format!("{lo}")
                } else {
                    format!("{lo}-{hi}")
                };
                size_hist.push((label, count));
            }
            if hi >= *all_cluster_sizes.last().unwrap_or(&1) {
                break;
            }
            lo = hi + 1;
        }
        let max_count = size_hist.iter().map(|(_, c)| *c).max().unwrap_or(1);
        log::info!("Cluster size distribution:");
        log::info!("  {:>12}  {:>8}  {}", "size", "clusters", "bar");
        for (label, count) in &size_hist {
            let bar_len = (*count as f64 / max_count as f64 * 40.0).ceil() as usize;
            let bar = "#".repeat(bar_len);
            log::info!("  {:>12}  {:>8}  {}", label, count, bar);
        }

        // Inter-cluster identity histogram
        log::info!("Computing inter-cluster pairwise identity for representative sequences ...");
        let mut ident_hist = vec![0usize; 21];
        let mut total_rep_pairs = 0usize;
        for seqs in &inter_cluster_seqs {
            for (ai, seq_a) in seqs.iter().enumerate() {
                for seq_b in &seqs[ai + 1..] {
                    total_rep_pairs += 1;
                    if let Some(aln) = aligner.align(seq_a, seq_b) {
                        let (identity, _) = alignment_stats(&aln.cigar);
                        let pct = identity * 100.0;
                        let bucket = ((pct / 5.0) as usize).min(20);
                        ident_hist[bucket] += 1;
                    }
                }
            }
        }
        log::info!("Inter-cluster representative identity ({} pairs):", total_rep_pairs.commas());
        log::info!("  {:>12}  {:>8}  {}", "identity%", "count", "bar");
        let max_count = ident_hist.iter().copied().max().unwrap_or(1);
        for (i, &count) in ident_hist.iter().enumerate() {
            let lo = i as f64 * 5.0;
            let hi = if i == 20 { 100.0 } else { lo + 5.0 };
            let bar_len = if max_count > 0 {
                (count as f64 / max_count as f64 * 40.0).ceil() as usize
            } else {
                0
            };
            let bar = "#".repeat(bar_len);
            log::info!("  [{lo:5.1},{hi:5.1})  {count:>8}  {bar}");
        }
    }

    match &args.output {
        Some(path) => log::info!("Done. Wrote {} sequences to {}.", written.commas(), path.display()),
        None => log::info!("Done. Wrote {} sequences.", written.commas()),
    }

    Ok(())
}
