use core::panic;
use std::{io::Write, sync::Arc, usize};

use ordered_float::OrderedFloat;

use crate::{
    align::{AlignParams, Aligner},
    config::{self},
    error::{ParallaxError, Result},
    index::Index,
    kmers::Kmer,
    reads::{
        builder::{Flag, SegmentBuilder},
        chains::write_clusters_debug,
        seeds::{Read, SeedCluster, SeedHit, SeedSaver, analyze_gap_fills},
    },
    reference::InMemoryReference,
    scores::compute_mapq_from_diff,
    utils::{
        GroupByTrait,
        debug::{self, DebugFile},
        hasher::FnvHasher,
        heap::{Heap, HeapOrdering, Heapable},
        range_set::RangeSet,
        sequence::reverse_complement_into,
    },
    writer::AlignmentWriter,
};

pub mod builder;
pub mod chains;
pub mod seeds;

enum AlignmentError {
    NoClusters,
    #[allow(dead_code)]
    LowQuality,
}

/// SAM flags
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
#[allow(dead_code)]
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Collector for seed clusters with reusable buffers.
///
/// This struct holds all the intermediate buffers needed for seeding,
/// merging, extension, and clustering. Reusing these buffers across
/// multiple calls avoids repeated allocation.
struct ClusterCollector {
    /// Seed hits collected from k-mer index lookups
    hits: Vec<SeedHit>,
    /// Temporary buffer for index lookups
    hit_vec: Vec<(usize, usize)>,
    /// Scratch space for merging/deduplication
    merge_scratch: Vec<SeedHit>,
    /// Batch buffer for prefetched lookups: (read_pos, kmer_value)
    kmer_batch: Vec<(usize, u64)>,
}

impl ClusterCollector {
    /// Create a new collector with empty buffers
    fn new() -> Self {
        ClusterCollector {
            hits: Vec::new(),
            hit_vec: Vec::new(),
            merge_scratch: Vec::new(),
            kmer_batch: Vec::new(),
        }
    }

    /// Collect seed clusters from a single strand.
    ///
    /// This performs seeding, merging, extension, and DBSCAN clustering, returning
    /// the resulting clusters without building alignments. This separation allows
    /// for cross-strand analysis before alignment construction.
    fn collect_from_strand<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        strand_qual: &[u8],
        is_reverse: bool,
        index: &Index<K, S>,
        reference: &InMemoryReference,
        read_name: &str,
    ) -> Vec<SeedCluster> {
        // Collect the seeds using the index.
        // Populates self.hits.
        if config::get().seeding.batch_prefetch {
            self.gather_seeds_batched::<K, S>(
                strand_seq,
                strand_qual,
                is_reverse,
                index,
                reference,
                read_name,
            );
        } else {
            self.gather_seeds::<K, S>(
                strand_seq,
                strand_qual,
                is_reverse,
                index,
                reference,
                read_name,
            );
        }

        // Phase 4 & 5: Cluster hits by diagonal using DBSCAN, then build LIS chains.
        // Important: We must partition by chromosome first, since hits from different
        // chromosomes should never be clustered together. Hits are already sorted by
        // (chrom_id, diagonal, ref_pos), so we find chromosome boundaries and process
        // each partition separately.

        let mut clusters = Vec::new();

        self.hits.sort_unstable_by(|a, b| {
            a.chrom_id
                .cmp(&b.chrom_id)
                .then(a.diagonal.cmp(&b.diagonal))
                .then(a.ref_pos.cmp(&b.ref_pos))
        });

        // Process each chromosome partition
        for (chrom_id, partition) in self.hits.group_by(|seed| seed.chrom_id) {
            if partition.is_empty() {
                continue;
            }
            let chrom_name = reference.chrom_name(chrom_id).to_string();
            log::debug!(
                "Processing chromosome {} with {} seeds",
                chrom_name,
                partition.len()
            );
            let mut seeds: Vec<SeedHit> = partition.to_vec();
            let chrom_clusters =
                chains::kruskal::collect_chains(&mut seeds, &chrom_name, is_reverse);
            //let chrom_clusters = if USE_AGGLOMERATIVE_CHAINING {
            //    chains::agglomerative::collect_chains(&mut seeds, &chrom_name, is_reverse)
            //} else {
            //    chains::rmq_dp::collect_chains(&mut seeds, is_reverse)
            //};
            write_clusters_debug(
                &chrom_clusters,
                read_name,
                &chrom_name,
                strand_seq,
                strand_qual,
                strand_seq.len(),
                is_reverse,
            );
            clusters.extend(chrom_clusters);
        }

        clusters
    }

    fn gather_seeds<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        strand_qual: &[u8],
        is_reverse: bool,
        index: &Index<K, S>,
        reference: &InMemoryReference,
        read_name: &str,
    ) {
        let cfg = config::get();
        let seq_len = strand_seq.len();
        self.hits.clear();

        // Phase 1: Collect seed hits using forward-only syncmers
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            strand_seq,
            [(); S],
            |pos, kmer| {
                self.hit_vec.clear();
                index.with(&kmer, |chrom_id, chrom_pos| {
                    self.hit_vec.push((chrom_id, chrom_pos));
                });
                let kmer_uniqueness = self.hit_vec.len() as u32;
                // Use seeds up to occurrence threshold
                if self.hit_vec.len() <= cfg.seeding.max_seed_occurrences {
                    for &(chrom_id, chrom_pos) in self.hit_vec.iter() {
                        self.hits.push(SeedHit::new(
                            chrom_id,
                            chrom_pos,
                            pos,
                            kmer.0,
                            kmer_uniqueness,
                            K,
                        ));
                    }
                }
            },
        );

        let strand_name = if is_reverse { "REV" } else { "FWD" };
        metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase()))
            .record(self.hits.len() as f64);

        // Phase 2: Sort hits - SeedHit's Ord gives us (chrom_id, diagonal, ref_pos) order
        self.hits.sort_unstable();

        // Phase 3: Merge overlapping/adjacent hits on same diagonal
        self.merge_scratch.clear();
        for hit in self.hits.drain(..) {
            if let Some(last) = self.merge_scratch.last_mut() {
                if last
                    .extend(
                        hit.chrom_id,
                        hit.ref_pos,
                        hit.read_pos,
                        hit.kmer,
                        hit.kmer_uniqueness,
                        K,
                    )
                    .is_none()
                {
                    continue; // Successfully merged
                }
            }
            self.merge_scratch.push(hit);
        }
        std::mem::swap(&mut self.hits, &mut self.merge_scratch);

        // Phase 3b: Extend each seed's exact match bidirectionally
        // This is the minimap2-style extension that maximizes anchor length
        for hit in self.hits.iter_mut() {
            let ref_seq = reference.get_seq(hit.chrom_id, 0, usize::MAX);
            hit.extend_exact(strand_seq, ref_seq);
        }

        // Phase 3c: Remove duplicates created by extension
        self.merge_scratch.clear();
        for hit in self.hits.drain(..) {
            if let Some(last) = self.merge_scratch.last() {
                if hit.chrom_id == last.chrom_id
                    && hit.diagonal == last.diagonal
                    && hit.ref_pos == last.ref_pos
                    && hit.match_len == last.match_len
                {
                    continue; // Duplicate, skip
                }
            }
            self.merge_scratch.push(hit);
        }
        std::mem::swap(&mut self.hits, &mut self.merge_scratch);

        // Write debug SAM output for seed hits (if debug file is initialized)
        if debug::is_enabled(DebugFile::Seeds) {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                debug::write(
                    DebugFile::Seeds,
                    &hit.to_sam_line(read_name, chrom_name, is_reverse, strand_seq, strand_qual),
                );
            }
        }
        // Write debug TSV output for seed hits (if debug file is initialized)
        if debug::is_enabled(DebugFile::SeedsTsv) {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                let strand = if is_reverse { "-" } else { "+" };
                // Convert strand coordinates to forward coordinates
                let (fwd_start, fwd_end) = if is_reverse {
                    (seq_len - hit.read_end(), seq_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                debug::write(
                    DebugFile::SeedsTsv,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        read_name,
                        fwd_start,
                        fwd_end,
                        seq_len,
                        chrom_name,
                        hit.ref_pos,
                        hit.ref_end(),
                        strand,
                        hit.match_len,
                    ),
                );
            }
        }

        if false {
            let out = std::fs::File::create("seeds.json").unwrap();
            let mut writer = std::io::BufWriter::new(out);
            let seed_saver = SeedSaver {
                seeds: self.hits.clone(),
                is_reverse,
                read: Read::new(read_name, strand_seq, strand_qual),
            };
            serde_json::to_writer(&mut writer, &seed_saver).unwrap();
            writer.flush().unwrap();
            panic!("Wrote seeds.jsonl");
        }
    }

    /// Batched version of `gather_seeds` that uses software-pipelined prefetching.
    ///
    /// Instead of interleaving k-mer generation with index lookups (which causes
    /// serial cache misses in the multi-GB hash tables), this method:
    /// 1. Generates all syncmer k-mers into a batch buffer (sequential writes, cache-friendly)
    /// 2. Looks up all k-mers using `Index::lookup_batch` which prefetches PIPE steps ahead
    ///
    /// The rest of the pipeline (sort, merge, extend, dedup) is identical.
    fn gather_seeds_batched<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        strand_qual: &[u8],
        is_reverse: bool,
        index: &Index<K, S>,
        reference: &InMemoryReference,
        read_name: &str,
    ) {
        let cfg = config::get();
        let seq_len = strand_seq.len();
        self.hits.clear();
        self.kmer_batch.clear();

        // Phase 1a: Generate all syncmers into the batch buffer
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            strand_seq,
            [(); S],
            |pos, kmer| {
                self.kmer_batch.push((pos, kmer.0));
            },
        );

        // Phase 1b: Batched lookup with prefetching
        let max_occ = cfg.seeding.max_seed_occurrences as u32;
        index.lookup_batch(&self.kmer_batch, |read_pos, chrom_id, chrom_pos, kmer_val, hit_count| {
            if hit_count <= max_occ {
                self.hits.push(SeedHit::new(
                    chrom_id,
                    chrom_pos,
                    read_pos,
                    kmer_val,
                    hit_count,
                    K,
                ));
            }
        });

        let strand_name = if is_reverse { "REV" } else { "FWD" };
        metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase()))
            .record(self.hits.len() as f64);

        // Phase 2: Sort hits - SeedHit's Ord gives us (chrom_id, diagonal, ref_pos) order
        self.hits.sort_unstable();

        // Phase 3: Merge overlapping/adjacent hits on same diagonal
        self.merge_scratch.clear();
        for hit in self.hits.drain(..) {
            if let Some(last) = self.merge_scratch.last_mut() {
                if last
                    .extend(
                        hit.chrom_id,
                        hit.ref_pos,
                        hit.read_pos,
                        hit.kmer,
                        hit.kmer_uniqueness,
                        K,
                    )
                    .is_none()
                {
                    continue; // Successfully merged
                }
            }
            self.merge_scratch.push(hit);
        }
        std::mem::swap(&mut self.hits, &mut self.merge_scratch);

        // Phase 3b: Extend each seed's exact match bidirectionally
        for hit in self.hits.iter_mut() {
            let ref_seq = reference.get_seq(hit.chrom_id, 0, usize::MAX);
            hit.extend_exact(strand_seq, ref_seq);
        }

        // Phase 3c: Remove duplicates created by extension
        self.merge_scratch.clear();
        for hit in self.hits.drain(..) {
            if let Some(last) = self.merge_scratch.last() {
                if hit.chrom_id == last.chrom_id
                    && hit.diagonal == last.diagonal
                    && hit.ref_pos == last.ref_pos
                    && hit.match_len == last.match_len
                {
                    continue; // Duplicate, skip
                }
            }
            self.merge_scratch.push(hit);
        }
        std::mem::swap(&mut self.hits, &mut self.merge_scratch);

        // Write debug SAM output for seed hits (if debug file is initialized)
        if debug::is_enabled(DebugFile::Seeds) {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                debug::write(
                    DebugFile::Seeds,
                    &hit.to_sam_line(read_name, chrom_name, is_reverse, strand_seq, strand_qual),
                );
            }
        }
        // Write debug TSV output for seed hits (if debug file is initialized)
        if debug::is_enabled(DebugFile::SeedsTsv) {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                let strand = if is_reverse { "-" } else { "+" };
                let (fwd_start, fwd_end) = if is_reverse {
                    (seq_len - hit.read_end(), seq_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                debug::write(
                    DebugFile::SeedsTsv,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fwd_start,
                        fwd_end,
                        seq_len,
                        chrom_name,
                        hit.ref_pos,
                        hit.ref_end(),
                        strand,
                        hit.match_len,
                    ),
                );
            }
        }
    }
}

/// Align a single read and output SAM record(s) using the provided writer.
///
/// This is the core alignment function that can be called from FASTQ or uBAM readers.
///
/// # Arguments
/// * `index` - The k-mer index for seed lookup
/// * `reference` - The reference genome
/// * `writer` - The AlignmentWriter to output SAM records
/// * `read_name` - Name of the read
/// * `seq` - Read sequence (forward strand)
/// * `qual` - Quality scores (same orientation as seq), or None if unavailable
    pub fn align_read<const K: usize, const S: usize, W: std::io::Write>(
        index: &Index<K, S>,
        reference: &InMemoryReference,
        writer: &AlignmentWriter<W>,
        read_name: &str,
        seq: &[u8],
        qual: &[u8],
        alignment_params: &AlignParams,
    ) {
        match align_read_inner(index, reference, writer, read_name, seq, qual, alignment_params) {
            Ok(()) => (),
            Err(AlignmentError::NoClusters) => {
                log::info!("Read {}: no seed clusters found, outputting unmapped", read_name);
                let _ = writer.write_alignment(
                    read_name,
                    FLAG_UNMAPPED,
                    "*",
                    0,
                    0,
                    "*",
                    "*",
                    0,
                    0,
                    std::str::from_utf8(seq).unwrap(),
                    std::str::from_utf8(qual).unwrap(),
                    "",
                );
            }
            Err(AlignmentError::LowQuality) => {
                log::info!("Read {}: all seed clusters filtered as low quality, outputting unmapped", read_name);
                let _ = writer.write_alignment(
                    read_name,
                    FLAG_UNMAPPED,
                    "*",
                    0,
                    0,
                    "*",
                    "*",
                    0,
                    0,
                    std::str::from_utf8(seq).unwrap(),
                    std::str::from_utf8(qual).unwrap(),
                    "",
                );
            }
        }
    }
    
fn align_read_inner<const K: usize, const S: usize, W: std::io::Write>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    writer: &AlignmentWriter<W>,
    read_name: &str,
    seq: &[u8],
    qual: &[u8],
    alignment_params: &AlignParams,
) -> std::result::Result<(), AlignmentError> {
    let alignment_start = std::time::Instant::now();

    let seq_len = seq.len();

    // Reusable cluster collector
    let mut collector = ClusterCollector::new();

    // Compute reverse complement for reverse strand processing
    let mut rc_seq = Vec::with_capacity(seq_len);
    reverse_complement_into(seq, &mut rc_seq);

    // Reverse quality scores for reverse strand (if available)
    let rc_qual: Vec<u8> = qual.iter().rev().copied().collect();

    // =========================================================================
    // PASS 1: Collect all seed clusters from both strands
    // =========================================================================

    let mut all_clusters: Vec<SeedCluster> = Vec::new();

    // Collect clusters from forward strand
    let fwd_clusters = collector.collect_from_strand(seq, qual, false, index, reference, read_name);
    all_clusters.extend(fwd_clusters);

    // Collect clusters from reverse strand
    let rev_clusters =
        collector.collect_from_strand(&rc_seq, &rc_qual, true, index, reference, read_name);
    all_clusters.extend(rev_clusters);

    all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));

    log::debug!(
        "Read {}: collected {} seed clusters from both strands (coverage {:.2}%)",
        read_name,
        all_clusters.len(),
        all_clusters
            .iter()
            .map(|c| c.read_coverage(seq_len))
            .sum::<f64>()
            * 100.0,
    );

    if false {
        for cluster in &all_clusters {
            let (qry_line, ref_line) = cluster.format_seed_diagram();
            log::info!("Cluster seed diagram:\n{}\n{}", qry_line, ref_line);
        }
    }
    // =========================================================================
    // PASS 1.5: Align gaps and split at failed alignments
    // =========================================================================
    // For each cluster, run WFA on all gaps. If any gap alignment fails (None),
    // split the cluster at that point. This must happen before gap-fill analysis
    // so we only consider clusters with valid internal alignments.

    let cfg = config::get();

    let mut aligner = Aligner::new();

    let mut new_clusters = Vec::new();
    for cluster in all_clusters.into_iter() {
        let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
        let chrom_len = reference.chrom_length(cluster.chrom_id) as usize;
        let ref_seq = reference.get_seq(cluster.chrom_id, 0, chrom_len);

        // Align all gaps in the cluster
        let min_seed_length = K / 2;
        let aligned_clusters = cluster.align_gaps(read_name, strand_seq, ref_seq, min_seed_length, &mut aligner);
        new_clusters.extend(aligned_clusters);
    }

    let mut all_clusters = new_clusters;
    all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));

    if all_clusters.is_empty() {
        return Err(AlignmentError::NoClusters);
    }

    // Write debug chains TSV output (if debug file is initialized)
    // Each seed and each gap between seeds gets its own row
    if debug::is_enabled(DebugFile::Chains) {
        for (i, cluster) in all_clusters.iter().enumerate() {
            let strand = if cluster.is_reverse { "-" } else { "+" };
            let chrom_name = reference.chrom_name(cluster.chrom_id);

            for (j, seed) in cluster.chain.iter().enumerate() {
                // Write gap row before this seed (if not the first seed)
                if j > 0 {
                    let prev = &cluster.chain[j - 1];
                    let gap_read_start = prev.read_end();
                    let gap_read_end = seed.read_pos;
                    let gap_ref_start = prev.ref_end();
                    let gap_ref_end = seed.ref_pos;
                    let read_width = gap_read_end.saturating_sub(gap_read_start) as i64;
                    let ref_width = gap_ref_end.saturating_sub(gap_ref_start) as i64;
                    debug::write(
                        DebugFile::Chains,
                        &format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0",
                            read_name,
                            i,
                            "gap",
                            gap_read_start,
                            gap_read_end,
                            read_width,
                            gap_ref_start,
                            gap_ref_end,
                            ref_width,
                            chrom_name,
                            strand,
                        ),
                    );
                }

                // Write seed row
                let read_width = seed.match_len;
                let ref_width = seed.match_len;
                debug::write(
                    DebugFile::Chains,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        read_name,
                        i,
                        "seed",
                        seed.read_pos,
                        seed.read_end(),
                        read_width,
                        seed.ref_pos,
                        seed.ref_end(),
                        ref_width,
                        chrom_name,
                        strand,
                        seed.kmer_uniqueness
                    ),
                );
            }
        }
    } // Write debug chain SAM with SA tags linking seeds
    if debug::is_enabled(DebugFile::ChainsSam) {
        for (i, cluster) in all_clusters.iter().enumerate() {
            let chrom_name = reference.chrom_name(cluster.chain[0].chrom_id);
            let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
            let strand_qual = if cluster.is_reverse { &rc_qual } else { qual };
            cluster.write_chain_sam(read_name, i, chrom_name, strand_seq, strand_qual);
        }
    }

    if log::log_enabled!(log::Level::Debug) {
        let mut all_gaps: Vec<(usize, usize, usize, usize, usize)> = Vec::new();
        for (i, cluster) in all_clusters.iter().enumerate() {
            let (read_start, read_end) = cluster.fwd_read_range(seq_len);
            log::debug!(
                "  Cluster {}: {}-{} {} seeds on {} strand (chrom: {},seed length: {},coverage {:.2}%, density {:.2})",
                i + 1,
                read_start,
                read_end,
                cluster.chain.len(),
                if cluster.is_reverse { "REV" } else { "FWD" },
                reference.chrom_name(cluster.chrom_id),
                cluster.total_seed_length(),
                cluster.read_coverage(seq_len) * 100.0,
                cluster.seed_density(),
            );

            for ((gap_start, gap_end), gap_index) in cluster.gaps(seq_len, 2 * K) {
                all_gaps.push((gap_start, gap_end, i, gap_index, cluster.chrom_id));
            }
        }

        all_gaps.sort_unstable();

        for (gap_start, gap_end, cluster_index, gap_index, chrom_id) in all_gaps.iter() {
            log::debug!(
                "  Gap: read {}-{} (length {}) in cluster {}/{} (chrom {})",
                gap_start,
                gap_end,
                gap_end - gap_start,
                cluster_index,
                gap_index,
                reference.chrom_name(*chrom_id),
            );
        }
    }

    // =========================================================================
    // PASS 1.6: Split clusters at gaps filled by other clusters
    // =========================================================================
    // Identify gaps where another cluster provides coverage, indicating a
    // potential chimeric breakpoint. Split the cluster at such gaps rather
    // than bridging them with WFA.

    let gap_fills = analyze_gap_fills(
        read_name,
        &all_clusters,
        seq_len,
        cfg.seeding.min_gap_for_split,
        2 * K,
        cfg.seeding.gap_fill_tolerance,
        alignment_params,
    );

    if !gap_fills.is_empty() {
        log::info!(
            "Read {}: found {} gap fills for potential splitting",
            read_name,
            gap_fills.len(),
        );

        // Group splits by cluster and sort by gap index descending
        // so we can apply splits from back to front without invalidating indices
        let mut splits_by_cluster: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for fill in &gap_fills {
            splits_by_cluster
                .entry(fill.cluster_idx)
                .or_default()
                .push(fill.gap_seed_idx);
        }

        // Sort each cluster's splits in descending order
        for indices in splits_by_cluster.values_mut() {
            indices.sort_unstable_by(|a, b| b.cmp(a));
            indices.dedup(); // Remove duplicate split points
        }

        // Apply splits (in descending cluster index order to preserve indices)
        let mut cluster_indices: Vec<_> = splits_by_cluster.keys().copied().collect();
        cluster_indices.sort_unstable_by(|a, b| b.cmp(a));

        for cluster_idx in cluster_indices {
            let split_indices = &splits_by_cluster[&cluster_idx];
            for &gap_seed_idx in split_indices {
                if let Some((new_cluster, _)) = all_clusters[cluster_idx].split_at_gap(gap_seed_idx)
                {
                    all_clusters.push(new_cluster);
                }
            }
        }

        // Re-sort after splitting
        all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));
    }

    let segment_sets = form_covering_sets(&all_clusters, read_name, seq_len);

    let set_scores: Vec<f64> = segment_sets
        .iter()
        .map(|set| {
            set.iter()
                .map(|cluster| cluster.quality(alignment_params).value())
                .sum()
        })
        .collect();

    let mut mapqs: Vec<Vec<f64>> = segment_sets
        .iter()
        .map(|set| {
            set.iter()
                .map(|cluster| cluster.quality(alignment_params).value())
                .collect()
        })
        .collect();

    for (i, mq) in mapqs.iter().enumerate() {
        let mut mq = mq.clone();
        mq.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        log::debug!(
            "Segment set {}: {} clusters, MQ range {:?}",
            i + 1,
            mq.len(),
            mq,
        );
    }

    if mapqs.len() > 0 {
        let n = mapqs[0].len();
        let ranges: Vec<(usize, usize)> = segment_sets[0]
            .iter()
            .map(|set| set.fwd_read_range(seq_len))
            .collect();
        let mut best_covering_score: Vec<f64> = vec![f64::NEG_INFINITY; n];
        for (i, set) in segment_sets.iter().enumerate().skip(1) {
            for cluster in set.iter() {
                let (start, end) = cluster.fwd_read_range(seq_len);
                for k in 0..n {
                    let (cov_start, cov_end) = ranges[k];
                    if end <= cov_start || start >= cov_end {
                        continue; // No overlap
                    }
                    if set_scores[i] > best_covering_score[k] {
                        best_covering_score[k] = set_scores[i];
                    }
                }
            }
        }
        log::debug!(
            "Best covering scores for secondary sets: {:?}",
            best_covering_score
        );
        let scale = 10.0; // MapQ scaling factor
        for k in 0..n {
            if best_covering_score[k] > set_scores[0] {
                mapqs[0][k] = 0.0; // Set MQ to 0 if covered by better secondary
            }
            if best_covering_score[k] < set_scores[0] {
                let num_seeds = segment_sets[0][k].total_seed_length() / K;
                let mq = compute_mapq_from_diff(
                    set_scores[0],
                    Some(best_covering_score[k]),
                    num_seeds,
                    scale,
                );
                mapqs[0][k] = mq as f64;
                log::info!(
                    "{}: Segment {} in primary set: score {:.2}, best covering score {:.2}, num seeds {}, assigned MQ {}",
                    read_name,
                    k + 1,
                    set_scores[0],
                    best_covering_score[k],
                    num_seeds,
                    mq,
                );
            }
        }

        log::info!("{}: After MQ adjustment, primary MQs: {:?}", read_name, mapqs[0],);

        // All secondary mappings are assigned MQ=0.
        for i in 1..mapqs.len() {
            for j in 0..mapqs[i].len() {
                mapqs[i][j] = 0.0;
            }
        }
    }

    for (i, set) in segment_sets.into_iter().enumerate() {
        let summaries = set
            .iter()
            .map(|cluster| cluster.summary(seq_len))
            .collect::<Vec<_>>();

        let leftmost_index = (0..set.len())
            .min_by_key(|&j| set[j].fwd_read_range(seq_len).0)
            .unwrap_or(0);
        let rightmost_index = (0..set.len())
            .max_by_key(|&j| set[j].fwd_read_range(seq_len).1)
            .unwrap_or(0);

        for (j, cluster) in set.iter().enumerate() {
            let mut flags = Vec::new();
            if cluster.is_reverse {
                flags.push(Flag::ReverseComplement);
            }
            if i > 0 {
                flags.push(Flag::SecondaryAlignment);
            }
            if j > 0 {
                flags.push(Flag::SupplementaryAlignment);
            }

            let primary = if i == 0 && j > 0 {
                let primary_cluster = &set[0];
                let rnext = reference.chrom_name(primary_cluster.chrom_id);
                let pnext = primary_cluster.ref_start() + 1;
                Some((rnext, pnext))
            } else {
                None
            };

            let summary: String = summaries
                .iter()
                .enumerate()
                .filter_map(|(k, s)| if k == j { Some(s) } else { None })
                .map(|(chrom_id, ref_pos, is_rc, cig, nm)| {
                    format!(
                        "{},{},{},{},{},{}",
                        reference.chrom_name(*chrom_id),
                        ref_pos,
                        if *is_rc { "-" } else { "+" },
                        cig,
                        mapqs[i][j],
                        nm,
                    )
                })
                .collect::<Vec<_>>()
                .join(";");

            let mc = cluster.chain.len();

            let soft_clip = i == 0 && j == 0;

            let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
            let strand_qual = if cluster.is_reverse { &rc_qual } else { qual };
            let chrom_len = reference.chrom_length(cluster.chrom_id) as usize;
            let ref_seq = reference.get_seq(cluster.chrom_id, 0, chrom_len);

            let fwd_extend_left = j == leftmost_index;
            let fwd_extend_right = j == rightmost_index;

            let (alignment, ref_start_adjustment, seq_start, seq_end) = cluster.clone().into_alignment(
                fwd_extend_left,
                fwd_extend_right,
                soft_clip,
                strand_seq,
                ref_seq,
                &mut aligner,
            );

            // Use the sequence range returned by into_alignment
            // This accounts for extensions (included) and hard clips (excluded)
            let seq_segment = &strand_seq[seq_start..seq_end];
            let qual_segment = &strand_qual[seq_start..seq_end];

            let cigar = alignment.cigar_string();

            // Adjust reference position: the left extension consumes reference bases
            // before the first seed, so subtract the adjustment
            let ref_pos = cluster.ref_start().saturating_sub(ref_start_adjustment) + 1;

            // Validate the alignment against sequences
            // For validation, we pass the full strand_seq and tell it where to start
            let query_start = alignment.leading_hard_clip();
            let ref_slice = &ref_seq[ref_pos - 1..]; // ref_pos is 1-based
            if let Err(err) = alignment.validate(ref_slice, strand_seq, query_start) {
                log::error!(
                    "Alignment validation failed for read {}: {} | \
                     cluster.read_range={:?}, fwd_extend_left={}, fwd_extend_right={}, \
                     ref_pos={}, strand_seq.len={}, seq_range={}..{}, cigar={}",
                    read_name,
                    err,
                    cluster.read_range(),
                    fwd_extend_left,
                    fwd_extend_right,
                    ref_pos,
                    strand_seq.len(),
                    seq_start,
                    seq_end,
                    cigar
                );
                panic!("Alignment validation failed");
            }

            let builder = SegmentBuilder::new(read_name)
                .with_flags(&flags)
                .with_reference(
                    reference.chrom_name(cluster.chrom_id),
                    ref_pos,
                )
                .with_mapping_quality(mapqs[i][j] as u8)
                .with_cigar(&cigar)
                .with_primary(primary)
                .with_sequence_and_quality(seq_segment, qual_segment)
                .with_tag_and_value("mc", mc)
                .with_tag_and_value("SA", summary);
            builder.write(writer).expect("write failed");
        }
    }

    metrics::histogram!("analysis_alignment").record(alignment_start.elapsed().as_secs_f64());

    Ok(())
}

type SegmentSet = (RangeSet, Vec<usize>); // (covered read segments, cluster indices)

struct SegmentSetHeap<'a> {
    clusters: &'a [SeedCluster],
    params: AlignParams,
}

impl<'a> Heapable for SegmentSetHeap<'a> {
    type Item = SegmentSet;

    const ORDERING: HeapOrdering = HeapOrdering::Max;

    fn cmp(&self, lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
        let l = lhs
            .1
            .iter()
            .map(|&i| self.clusters[i].quality(&self.params).0)
            .sum::<f64>();
        let r = rhs
            .1
            .iter()
            .map(|&i| self.clusters[i].quality(&self.params).0)
            .sum::<f64>();
        l.partial_cmp(&r).unwrap_or(std::cmp::Ordering::Equal)
    }
}

fn form_covering_sets(
    clusters: &[SeedCluster],
    read_name: &str,
    read_len: usize,
) -> Vec<Vec<SeedCluster>> {
    let mut order_by_quality: Vec<usize> = (0..clusters.len()).collect();
    let params = AlignParams::default();
    order_by_quality.sort_by_key(|i| OrderedFloat(-clusters[*i].quality(&params).0));

    let mut segment_set_heap = Heap::new(SegmentSetHeap { clusters, params });
    let mut wanted_segment_set: Option<SegmentSet> = None;
    let mut stack: Vec<SegmentSet> = vec![];

    for &i in order_by_quality.iter() {
        let cluster = &clusters[i];
        let (read_start, read_end) = cluster.fwd_read_range(read_len);

        // Find the highest quality segment set that does not overlap with this cluster's read range
        while let Some(segment_set) = segment_set_heap.pop() {
            if segment_set.0.overlaps(&(read_start, read_end)) {
                stack.push(segment_set);
            } else {
                wanted_segment_set = Some(segment_set);
                break;
            }
        }

        // If we found a non-overlapping segment set, add this cluster to it. Otherwise, create a new segment set for this cluster.
        if let Some((mut ranges, mut set)) = wanted_segment_set.take() {
            assert!(!ranges.overlaps(&(read_start, read_end)));
            ranges.add_range(read_start, read_end);
            set.push(i);
            segment_set_heap.push((ranges, set));
        } else {
            let mut ranges = RangeSet::new();
            ranges.add_range(read_start, read_end);
            let set = vec![i];
            segment_set_heap.push((ranges, set));
        }

        let best = stack.len() + 1;

        // Put back the segment sets we popped off
        while let Some(segment_set) = stack.pop() {
            segment_set_heap.push(segment_set);
        }

        log::debug!(
            "{}/{}: read {}-{} (length {}) assigned to segment set {}, quality {:.2}",
            read_name,
            i,
            read_start,
            read_end,
            read_end - read_start,
            best,
            cluster.quality(&params).0
        );
    }

    let n = segment_set_heap.len();
    log::debug!(
        "Read {}: assigned {} clusters to {} segment sets",
        read_name,
        clusters.len(),
        n
    );

    let mut segment_set_index: Vec<usize> = vec![0; clusters.len()];
    for (idx, (_, set)) in segment_set_heap.drain().enumerate() {
        for &i in set.iter() {
            segment_set_index[i] = idx;
        }
    }

    let mut segment_sets: Vec<Vec<SeedCluster>> = (0..n).map(|_| vec![]).collect();
    for (i, cluster) in clusters.into_iter().enumerate() {
        let idx = segment_set_index[i];
        segment_sets[idx].push(cluster.clone());
    }

    segment_sets
}

/// Process reads from a FASTQ file (handles gzip, bzip2, xz compression transparently)
#[allow(dead_code)]
pub fn process_reads_from_fastq<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    fastq: &str,
    command_line: &str,
    read_group_header: Option<&str>,
) -> Result<()> {
    log::info!("Processing reads from {}", fastq);

    let params = AlignParams::default();

    let stdout = std::io::stdout();
    let writer = AlignmentWriter::builder(stdout.lock())
        .add_contigs(reference.chromosomes())
        .read_group(read_group_header.map(String::from))
        .command_line(command_line)
        .build()?;

    let (decompressed_reader, format) = niffler::from_path(std::path::Path::new(fastq))
        .map_err(|e| ParallaxError::Other(Box::new(e)))?;
    if format != niffler::Format::No {
        log::info!("Detected {:?} compression", format);
    }
    let reader = std::io::BufReader::new(decompressed_reader);
    let mut reader = noodles::fastq::io::Reader::new(reader);

    for record in reader.records() {
        let record = record?;
        let read_name = std::str::from_utf8(record.name()).unwrap_or("?");
        let seq: &[u8] = record.sequence().as_ref();
        let qual: &[u8] = record.quality_scores().as_ref();

        align_read(index, reference, &writer, read_name, seq, qual, &params);
    }

    writer.flush()?;
    Ok(())
}

/// A read to be processed by a worker thread
struct ReadWork {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Process reads from a FASTQ file using multiple threads.
///
/// Reads are distributed to worker threads via a channel. The InMemoryReference
/// is shared across all threads via Arc (no per-thread cloning needed).
pub fn process_reads_parallel<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    fastq: &str,
    sam: Option<&str>,
    num_threads: usize,
    command_line: &str,
    read_group_header: Option<&str>,
) -> Result<()> {
    use crossbeam::channel::bounded;

    log::info!(
        "Processing reads from {} using {} threads",
        fastq,
        num_threads
    );

    let now = std::time::Instant::now();
    let mut num_records = 0;

    // Initialize debug files from config
    let cfg = config::get();
    debug::init(&cfg, reference)?;

    let params = AlignParams::default();

    // Create writer with headers - either to file or stdout
    let output: Box<dyn std::io::Write + Send> = match sam {
        Some(path) => {
            log::info!("Writing output to {}", path);
            Box::new(std::fs::File::create(path)?)
        }
        None => Box::new(std::io::stdout()),
    };
    let writer = Arc::new(
        AlignmentWriter::builder(output)
            .add_contigs(reference.chromosomes())
            .read_group(read_group_header.map(String::from))
            .command_line(command_line)
            .build()?,
    );

    // Create a bounded channel for backpressure
    let (sender, receiver) = bounded::<ReadWork>(num_threads * 100);

    // Use crossbeam's scoped threads to safely borrow index, reference, and writer
    crossbeam::scope(|scope| {
        // Spawn worker threads
        for _ in 0..num_threads {
            let receiver = receiver.clone();
            let writer = writer.clone();
            scope.spawn(move |_| {
                while let Ok(work) = receiver.recv() {
                    align_read(
                        index, reference, &writer, &work.name, &work.seq, &work.qual, &params,
                    );
                }
            });
        }

        // Read FASTQ and send to workers (in main thread)
        // Use niffler for transparent decompression (gzip, bzip2, xz)
        let (decompressed_reader, format) =
            niffler::from_path(std::path::Path::new(fastq)).expect("Failed to open FASTQ file");
        if format != niffler::Format::No {
            log::info!("Detected {:?} compression", format);
        }
        let reader = std::io::BufReader::new(decompressed_reader);
        let mut reader = noodles::fastq::io::Reader::new(reader);

        for record in reader.records() {
            let record = record.expect("Failed to read FASTQ record");
            let seq: &[u8] = record.sequence().as_ref();
            let qual: &[u8] = record.quality_scores().as_ref();
            let work = ReadWork {
                name: String::from_utf8_lossy(record.name()).into_owned(),
                seq: seq.to_vec(),
                qual: qual.to_vec(),
            };
            sender.send(work).expect("Failed to send work to thread");
            num_records += 1;
        }

        // Signal completion by dropping sender
        drop(sender);

        // Scoped threads automatically join when scope ends
    })
    .expect("Scoped thread panicked");

    writer.flush()?;

    // Flush all debug files
    debug::flush_all();

    let elapsed = now.elapsed();
    log::info!(
        "Completed processing reads {} from {} in {:.2?}",
        num_records,
        fastq,
        elapsed
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create a SeedHit with a dummy kmer value
    fn make_hit(chrom_id: usize, ref_pos: usize, read_pos: usize, match_len: usize) -> SeedHit {
        SeedHit::new(chrom_id, ref_pos, read_pos, 0, 1, match_len)
    }

    #[test]
    fn test_seed_hit_new() {
        let hit = SeedHit::new(1, 100, 50, 12345, 1, 20);
        assert_eq!(hit.chrom_id, 1);
        assert_eq!(hit.ref_pos, 100);
        assert_eq!(hit.read_pos, 50);
        assert_eq!(hit.kmer, 12345);
        assert_eq!(hit.kmer_uniqueness, 1);
        assert_eq!(hit.match_len, 20);
        // diagonal = ref_pos - read_pos = 100 - 50 = 50
        assert_eq!(hit.diagonal, 50);
    }

    #[test]
    fn test_seed_hit_read_end() {
        let hit = make_hit(0, 100, 50, 20);
        assert_eq!(hit.read_end(), 70); // 50 + 20
    }

    #[test]
    fn test_seed_hit_ref_end() {
        let hit = make_hit(0, 100, 50, 20);
        assert_eq!(hit.ref_end(), 120); // 100 + 20
    }

    #[test]
    fn test_seed_hit_diagonal_calculation() {
        // Forward diagonal (ref ahead of read)
        let hit1 = make_hit(0, 1000, 100, 20);
        assert_eq!(hit1.diagonal, 900);

        // Negative diagonal (read ahead of ref)
        let hit2 = make_hit(0, 100, 1000, 20);
        assert_eq!(hit2.diagonal, -900);

        // Zero diagonal (same position)
        let hit3 = make_hit(0, 500, 500, 20);
        assert_eq!(hit3.diagonal, 0);
    }

    #[test]
    fn test_extend_overlapping_same_diagonal() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);

        // Second seed at read pos 10, ref pos 110, length 20
        // This is on the same diagonal (110-10 = 100-0 = 100)
        // And overlaps: ref 110 < ref_end 120, gap is 10 < k=20
        let k = 20;
        let result = hit.extend(0, 110, 10, 0, 1, k);

        assert!(
            result.is_none(),
            "Should extend in place, not return new hit"
        );
        // New end should be: (110 - 100) + 20 = 30
        assert_eq!(hit.match_len, 30);
        assert_eq!(hit.ref_end(), 130);
        assert_eq!(hit.read_end(), 30);
    }

    #[test]
    fn test_extend_adjacent_same_diagonal() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);

        // Second seed starts exactly where first ends
        // read pos 20, ref pos 120, still same diagonal
        let k = 20;
        let result = hit.extend(0, 120, 20, 0, 1, k);

        // Gap is exactly 20, which equals k, so should NOT extend
        // because condition is: chrom_pos - self.ref_pos < self.match_len + k
        // 120 - 100 = 20 < 20 + 20 = 40, so it SHOULD extend
        assert!(result.is_none(), "Should extend in place");
        assert_eq!(hit.match_len, 40);
    }

    #[test]
    fn test_extend_gap_too_large() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Second seed with large gap (beyond match_len + k)
        // ref pos 200, read pos 100 (same diagonal = 100)
        // Gap check: 200 - 100 = 100 >= 20 + 20 = 40
        let k = 20;
        let result = hit.extend(0, 200, 100, 999, 1, k);

        assert!(result.is_some(), "Should return new hit due to large gap");
        assert_eq!(hit.match_len, original_len, "Original should be unchanged");

        let new_hit = result.unwrap();
        assert_eq!(new_hit.ref_pos, 200);
        assert_eq!(new_hit.read_pos, 100);
        assert_eq!(new_hit.kmer, 999);
    }

    #[test]
    fn test_extend_different_chromosome() {
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Same positions but different chromosome
        let k = 20;
        let result = hit.extend(1, 110, 10, 0, 1, k);

        assert!(
            result.is_some(),
            "Different chromosome should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
        assert_eq!(result.unwrap().chrom_id, 1);
    }

    #[test]
    fn test_extend_different_diagonal() {
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Different diagonal: ref_pos - read_pos = 111 - 10 = 101 != 100
        let k = 20;
        let result = hit.extend(0, 111, 10, 0, 1, k);

        assert!(result.is_some(), "Different diagonal should create new hit");
        assert_eq!(hit.match_len, original_len);

        let new_hit = result.unwrap();
        assert_eq!(new_hit.diagonal, 101);
    }

    #[test]
    fn test_extend_backwards_ref_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;

        // New ref_pos before current ref_pos
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, 1, k);

        assert!(
            result.is_some(),
            "Backwards ref position should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_backwards_read_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;

        // New read_pos before current read_pos (even if same diagonal)
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, 1, k);

        assert!(
            result.is_some(),
            "Backwards read position should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_fully_contained() {
        // Seed covering positions 0-20 in read, 100-120 in ref
        let mut hit = make_hit(0, 100, 0, 20);

        // New seed at pos 5-25 overlaps significantly
        // ref 105, read 5, same diagonal (100)
        let k = 20;
        let result = hit.extend(0, 105, 5, 0, 1, k);

        assert!(result.is_none(), "Overlapping hit should extend");
        // New end: (105 - 100) + 20 = 25
        assert_eq!(hit.match_len, 25);
    }

    #[test]
    fn test_extend_no_length_change_if_contained() {
        // Seed covering positions 0-30 in read
        let mut hit = make_hit(0, 100, 0, 30);

        // New seed fully contained within existing match
        // ref 110, read 10, k=20 means it ends at read 30, ref 130
        // That's exactly where the original ends, so no extension needed
        let k = 20;
        let result = hit.extend(0, 110, 10, 0, 1, k);

        assert!(result.is_none(), "Contained hit should not create new hit");
        // (110 - 100) + 20 = 30, which equals original, so no change
        assert_eq!(hit.match_len, 30);
    }

    #[test]
    fn test_extend_sequence_of_hits() {
        let k = 20;
        let mut hit = make_hit(0, 100, 0, k);

        // Simulate a sequence of overlapping syncmers ~6 bases apart
        // All on the same diagonal
        for i in 1..10 {
            let read_pos = i * 6;
            let ref_pos = 100 + i * 6;
            let result = hit.extend(0, ref_pos, read_pos, 0, 1, k);
            assert!(result.is_none(), "Hit {} should extend in place", i);
        }

        // Final length should cover from 0 to (9*6 + 20) = 74
        assert_eq!(hit.match_len, 9 * 6 + k);
        assert_eq!(hit.read_end(), 74);
        assert_eq!(hit.ref_end(), 174);
    }

    // =========================================================================
    // SeedCluster tests
    // =========================================================================

    #[test]
    fn test_seed_cluster_new_sorts_by_read_pos() {
        // Create seeds in reverse read_pos order
        let seeds = vec![
            make_hit(0, 300, 200, 20), // read_pos = 200
            make_hit(0, 100, 0, 20),   // read_pos = 0
            make_hit(0, 200, 100, 20), // read_pos = 100
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Should be sorted by read_pos
        assert_eq!(cluster.chain[0].read_pos, 0);
        assert_eq!(cluster.chain[1].read_pos, 100);
        assert_eq!(cluster.chain[2].read_pos, 200);

        assert_eq!(cluster.read_start, 0);
        assert_eq!(cluster.read_end, 220); // 200 + 20
    }

    #[test]
    fn test_seed_cluster_empty_returns_none() {
        let seeds: Vec<SeedHit> = vec![];
        assert!(SeedCluster::new(seeds, false, 1).is_none());
    }

    #[test]
    fn test_seed_cluster_fwd_read_range_forward_strand() {
        let seeds = vec![make_hit(0, 100, 50, 20)];
        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        let (start, end) = cluster.fwd_read_range(1000);
        assert_eq!(start, 50);
        assert_eq!(end, 70);
    }

    #[test]
    fn test_seed_cluster_fwd_read_range_reverse_strand() {
        // For reverse strand, coordinates need to be flipped
        let seeds = vec![make_hit(0, 100, 50, 20)];
        let cluster = SeedCluster::new(seeds, true, 1).unwrap();

        // read_start=50, read_end=70, read_len=1000
        // fwd_start = 1000 - 70 = 930, fwd_end = 1000 - 50 = 950
        let (start, end) = cluster.fwd_read_range(1000);
        assert_eq!(start, 930);
        assert_eq!(end, 950);
    }

    #[test]
    fn test_seed_cluster_split_at_gap() {
        use crate::align::{Alignment, CigarOp};
        use crate::scores::DivergenceScore;

        // Create a chain with a gap between seeds 1 and 2
        let seeds = vec![
            make_hit(0, 100, 0, 20),   // 0-20
            make_hit(0, 200, 50, 20),  // 50-70
            make_hit(0, 400, 200, 20), // 200-220 (gap here)
            make_hit(0, 500, 250, 20), // 250-270
        ];

        let mut cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Add dummy gap alignments (3 gaps for 4 seeds)
        cluster.gap_alignments = vec![
            Alignment {
                divergence: DivergenceScore::new(10.0),
                cigar: vec![CigarOp::Match(10)],
            },
            Alignment {
                divergence: DivergenceScore::new(20.0),
                cigar: vec![CigarOp::Match(20)],
            },
            Alignment {
                divergence: DivergenceScore::new(15.0),
                cigar: vec![CigarOp::Match(15)],
            },
        ];

        assert_eq!(cluster.chain.len(), 4);

        // Split at gap between index 1 and 2
        let (tail, _dropped_alignment) = cluster.split_at_gap(1).unwrap();

        // Original cluster should have seeds 0, 1
        assert_eq!(cluster.chain.len(), 2);
        assert_eq!(cluster.read_start, 0);
        assert_eq!(cluster.read_end, 70);

        // Tail cluster should have seeds 2, 3
        assert_eq!(tail.chain.len(), 2);
        assert_eq!(tail.read_start, 200);
        assert_eq!(tail.read_end, 270);
        assert_eq!(tail.is_reverse, cluster.is_reverse);
    }

    #[test]
    fn test_seed_cluster_split_preserves_strand() {
        use crate::align::{Alignment, CigarOp};
        use crate::scores::DivergenceScore;

        let seeds = vec![make_hit(0, 100, 0, 20), make_hit(0, 300, 100, 20)];

        let mut cluster = SeedCluster::new(seeds, true, 1).unwrap();

        // Add dummy gap alignment (1 gap for 2 seeds)
        cluster.gap_alignments = vec![Alignment {
            divergence: DivergenceScore::new(10.0),
            cigar: vec![CigarOp::Match(10)],
        }];

        let (tail, _dropped_alignment) = cluster.split_at_gap(0).unwrap();

        assert!(cluster.is_reverse);
        assert!(tail.is_reverse);
    }

    #[test]
    fn test_seed_cluster_gaps_forward_strand() {
        let seeds = vec![
            make_hit(0, 100, 0, 20),   // 0-20
            make_hit(0, 200, 50, 20),  // 50-70, gap of 30
            make_hit(0, 400, 200, 20), // 200-220, gap of 130
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();
        let gaps: Vec<_> = cluster.gaps(1000, 50).collect();

        // Only the gap of 130 should be returned (min_gap=50)
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].0, (70, 200)); // gap from 70 to 200
        assert_eq!(gaps[0].1, 1); // seed index 1 before this gap
    }

    #[test]
    fn test_seed_cluster_gaps_reverse_strand() {
        // For reverse strand, the chain is in RC coordinates
        // but gaps() should return forward-strand coordinates
        let seeds = vec![
            make_hit(0, 100, 0, 20),   // RC pos 0-20
            make_hit(0, 200, 50, 20),  // RC pos 50-70, gap of 30
            make_hit(0, 400, 200, 20), // RC pos 200-220, gap of 130
        ];

        let cluster = SeedCluster::new(seeds, true, 1).unwrap();
        let read_len = 1000;
        let gaps: Vec<_> = cluster.gaps(read_len, 50).collect();

        // Gap in RC coords: 70-200
        // In forward coords: (1000-200, 1000-70) = (800, 930)
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].0, (800, 930));
    }

    // =========================================================================
    // Chain colinearity tests
    // =========================================================================

    #[test]
    fn test_colinear_chain_both_dimensions_increasing() {
        // A proper colinear chain should have both ref_pos and read_pos increasing
        let seeds = vec![
            make_hit(0, 100, 10, 20),
            make_hit(0, 200, 50, 20),
            make_hit(0, 300, 100, 20),
            make_hit(0, 400, 150, 20),
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Verify both dimensions are strictly increasing
        for i in 1..cluster.chain.len() {
            assert!(
                cluster.chain[i].read_pos > cluster.chain[i - 1].read_end() - 1,
                "read_pos should be increasing: {} vs {}",
                cluster.chain[i].read_pos,
                cluster.chain[i - 1].read_end()
            );
            assert!(
                cluster.chain[i].ref_pos >= cluster.chain[i - 1].ref_end(),
                "ref_pos should be increasing: {} vs {}",
                cluster.chain[i].ref_pos,
                cluster.chain[i - 1].ref_end()
            );
        }
    }

    #[test]
    fn test_chain_ref_pos_monotonic_after_read_sort() {
        // This tests the invariant that should hold after SeedCluster::new
        // Even if seeds come in arbitrary order, after sorting by read_pos,
        // ref_pos should also be increasing for a proper colinear chain
        let seeds = vec![
            make_hit(0, 400, 150, 20), // Will be last after sort
            make_hit(0, 100, 10, 20),  // Will be first after sort
            make_hit(0, 300, 100, 20), // Will be third after sort
            make_hit(0, 200, 50, 20),  // Will be second after sort
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Verify ref_pos is monotonically increasing
        for i in 1..cluster.chain.len() {
            assert!(
                cluster.chain[i].ref_pos >= cluster.chain[i - 1].ref_pos,
                "ref_pos not monotonic at {}: {} < {}",
                i,
                cluster.chain[i].ref_pos,
                cluster.chain[i - 1].ref_pos
            );
        }
    }
}
