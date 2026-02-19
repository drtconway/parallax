use core::panic;
use std::{io::Write, sync::Arc, usize};

use noodles::sam::alignment::{
    record::Flags,
    record::data::field::Tag,
    record_buf::{Data, data::field::Value},
};
use ordered_float::OrderedFloat;

use crate::{
    align::{AlignParams, Aligner, Alignment, Kind, Op},
    config::{self},
    error::Result,
    index::Index,
    kmers::Kmer,
    reads::{
        builder::{build_record, build_unmapped_record},
        chains::write_clusters_debug,
        seeds::{Read, SeedCluster, SeedHit, SeedSaver, analyze_gap_fills},
    },
    reference::InMemoryReference,
    scores::compute_mapq_from_diff,
    utils::{
        GroupByTrait,
        debug::{self, DebugFile, DebugOutput, DebugTsvWriter, TsvRow},
        hasher::FnvHasher,
        heap::{Heap, HeapOrdering, Heapable},
        range_set::RangeSet,
        sequence::reverse_complement_into,
    },
    writer::{AlignmentWriter, OutputFormat},
};

pub mod builder;
pub mod chains;
pub mod seeds;

// ── Debug file statics ──────────────────────────────────────────────────────

/// Debug SAM file with extended seeds (before clustering).
static SEEDS_SAM: DebugFile<SeedsSamDebug> = DebugFile::new();

/// Debug TSV file with candidate seeds.
static SEEDS_TSV: DebugFile<SeedsTsvDebug> = DebugFile::new();

/// Debug TSV file with seed chains (after chaining, before alignment).
static CHAINS_TSV: DebugFile<ChainsTsvDebug> = DebugFile::new();

// ── Concrete debug types ─────────────────────────────────────────────────────

pub(crate) struct SeedsSamDebug(DebugTsvWriter);

impl DebugOutput for SeedsSamDebug {
    type Item<'a> = str;
    fn create() -> Option<Self> {
        let path = &config::get().seeding.debug_seeds_sam;
        if path.is_empty() { return None; }
        DebugTsvWriter::open(path, debug::sam_header().as_deref()).ok().map(Self)
    }
    fn append(&self, item: &str) { self.0.append(item); }
    fn finish(&self) { self.0.finish(); }
}

type SeedsTsvRow<'a> = (&'a str, usize, usize, usize, &'a str, usize, usize, &'a str, usize);

pub(crate) struct SeedsTsvDebug(DebugTsvWriter);

impl SeedsTsvDebug {
    const HEADERS: &[&str] = &[
        "read_name", "read_start", "read_end", "read_len",
        "chrom", "ref_start", "ref_end", "strand", "score",
    ];
    const _CHECK: () = assert!(Self::HEADERS.len() == <SeedsTsvRow<'static> as TsvRow>::NUM_FIELDS);
}

impl DebugOutput for SeedsTsvDebug {
    type Item<'a> = SeedsTsvRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_seeds_tsv;
        if path.is_empty() { return None; }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &SeedsTsvRow<'_>) { self.0.append_row(item); }
    fn finish(&self) { self.0.finish(); }
}

type ChainsTsvRow<'a> = (&'a str, usize, &'a str, usize, usize, usize, usize, usize, usize, &'a str, &'a str, u32);

pub(crate) struct ChainsTsvDebug(DebugTsvWriter);

impl ChainsTsvDebug {
    const HEADERS: &[&str] = &[
        "read_name", "cluster_id", "row_type", "read_start", "read_end", "read_width",
        "ref_start", "ref_end", "ref_width", "chrom", "strand", "uniqueness",
    ];
    const _CHECK: () = assert!(Self::HEADERS.len() == <ChainsTsvRow<'static> as TsvRow>::NUM_FIELDS);
}

impl DebugOutput for ChainsTsvDebug {
    type Item<'a> = ChainsTsvRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_chains_tsv;
        if path.is_empty() { return None; }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &ChainsTsvRow<'_>) { self.0.append_row(item); }
    fn finish(&self) { self.0.finish(); }
}

enum AlignmentError {
    NoClusters,
    #[allow(dead_code)]
    LowQuality,
}

/// Detected input file format.
enum InputFormat {
    Fastq,
    Bam,
}

/// Detect input format from file extension.
fn detect_input_format(path: &str) -> InputFormat {
    let path_lower = path.to_lowercase();
    if path_lower.ends_with(".bam") {
        InputFormat::Bam
    } else {
        // Default to FASTQ for .fq, .fastq, .fq.gz, .fastq.gz, etc.
        InputFormat::Fastq
    }
}

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

        // Write debug SAM output for seed hits
        if SEEDS_SAM.is_enabled() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                SEEDS_SAM.append(
                    &hit.to_sam_line(read_name, chrom_name, is_reverse, strand_seq, strand_qual),
                );
            }
        }
        // Write debug TSV output for seed hits
        if SEEDS_TSV.is_enabled() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                let strand = if is_reverse { "-" } else { "+" };
                // Convert strand coordinates to forward coordinates
                let (fwd_start, fwd_end) = if is_reverse {
                    (seq_len - hit.read_end(), seq_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                SEEDS_TSV.append(&(
                    read_name, fwd_start, fwd_end, seq_len,
                    chrom_name, hit.ref_pos, hit.ref_end(), strand, hit.match_len,
                ));
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
        index.lookup_batch(
            &self.kmer_batch,
            |read_pos, chrom_id, chrom_pos, kmer_val, hit_count| {
                if hit_count <= max_occ {
                    self.hits.push(SeedHit::new(
                        chrom_id, chrom_pos, read_pos, kmer_val, hit_count, K,
                    ));
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

        // Write debug SAM output for seed hits
        if SEEDS_SAM.is_enabled() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                SEEDS_SAM.append(
                    &hit.to_sam_line(read_name, chrom_name, is_reverse, strand_seq, strand_qual),
                );
            }
        }
        // Write debug TSV output for seed hits
        if SEEDS_TSV.is_enabled() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                let strand = if is_reverse { "-" } else { "+" };
                let (fwd_start, fwd_end) = if is_reverse {
                    (seq_len - hit.read_end(), seq_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                SEEDS_TSV.append(&(
                    read_name, fwd_start, fwd_end, seq_len,
                    chrom_name, hit.ref_pos, hit.ref_end(), strand, hit.match_len,
                ));
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
pub fn align_read<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    writer: &AlignmentWriter,
    read_name: &str,
    seq: &[u8],
    qual: &[u8],
    alignment_params: &AlignParams,
) {
    match align_read_inner(
        index,
        reference,
        writer,
        read_name,
        seq,
        qual,
        alignment_params,
    ) {
        Ok(()) => (),
        Err(AlignmentError::NoClusters) => {
            log::info!(
                "Read {}: no seed clusters found, outputting unmapped",
                read_name
            );
            let record = build_unmapped_record(read_name, seq, qual);
            let _ = writer.write_record(&record);
        }
        Err(AlignmentError::LowQuality) => {
            log::info!(
                "Read {}: all seed clusters filtered as low quality, outputting unmapped",
                read_name
            );
            let record = build_unmapped_record(read_name, seq, qual);
            let _ = writer.write_record(&record);
        }
    }
}

fn align_read_inner<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    writer: &AlignmentWriter,
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
        let aligned_clusters = cluster.align_gaps(
            read_name,
            strand_seq,
            ref_seq,
            min_seed_length,
            &mut aligner,
        );
        new_clusters.extend(aligned_clusters);
    }

    let mut all_clusters = new_clusters;
    all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));

    if all_clusters.is_empty() {
        return Err(AlignmentError::NoClusters);
    }

    // Write debug chains TSV output
    // Each seed and each gap between seeds gets its own row
    if CHAINS_TSV.is_enabled() {
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
                    let read_width = gap_read_end.saturating_sub(gap_read_start);
                    let ref_width = gap_ref_end.saturating_sub(gap_ref_start);
                    CHAINS_TSV.append(&(
                        read_name, i, "gap",
                        gap_read_start, gap_read_end, read_width,
                        gap_ref_start, gap_ref_end, ref_width,
                        chrom_name, strand, 0u32,
                    ));
                }

                // Write seed row
                let read_width = seed.match_len;
                let ref_width = seed.match_len;
                CHAINS_TSV.append(&(
                    read_name, i, "seed",
                    seed.read_pos, seed.read_end(), read_width,
                    seed.ref_pos, seed.ref_end(), ref_width,
                    chrom_name, strand, seed.kmer_uniqueness,
                ));
            }
        }
    } // Write debug chain SAM with SA tags linking seeds
    if seeds::CHAINS_SAM.is_enabled() {
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
        log::debug!(
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
        .map(|set| score_clusters(set.iter(), seq_len, alignment_params))
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
                log::debug!(
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

        log::debug!(
            "{}: After MQ adjustment, primary MQs: {:?}",
            read_name,
            mapqs[0],
        );

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

        // Sort cluster indices by forward-strand read position.
        let mut pos_order: Vec<usize> = (0..set.len()).collect();
        pos_order.sort_by_key(|&j| set[j].fwd_read_range(seq_len).0);

        // Effective ranges in fwd-read coordinates, indexed by position order.
        // Initially the seed ranges; updated as extensions are computed.
        let mut effective_ranges: Vec<(usize, usize)> = pos_order
            .iter()
            .map(|&j| set[j].fwd_read_range(seq_len))
            .collect();

        // Process clusters in quality-descending order.
        // For each, compute the remaining left/right gap from the effective
        // ranges of the neighboring clusters (in positional order).
        let params = AlignParams::default();
        let mut quality_order: Vec<usize> = (0..pos_order.len()).collect(); // indices into pos_order
        quality_order.sort_by_key(|&pi| {
            let j = pos_order[pi];
            OrderedFloat(-set[j].quality(&params).value())
        });

        // Store per-cluster results indexed by original set index j.
        struct ClusterResult {
            alignment: Alignment,
            ref_start_adjustment: usize,
            seq_start: usize,
            seq_end: usize,
        }
        let mut results: Vec<Option<ClusterResult>> = (0..set.len()).map(|_| None).collect();

        for &pi in &quality_order {
            let j = pos_order[pi]; // original index in set

            // Compute left budget: gap from previous cluster's effective end (or read start)
            let left_bound = if pi == 0 {
                0
            } else {
                effective_ranges[pi - 1].1
            };
            let left_budget = effective_ranges[pi].0.saturating_sub(left_bound);

            // Compute right budget: gap to next cluster's effective start (or read end)
            let right_bound = if pi == pos_order.len() - 1 {
                seq_len
            } else {
                effective_ranges[pi + 1].0
            };
            let right_budget = right_bound.saturating_sub(effective_ranges[pi].1);

            let cluster = &set[j];
            let soft_clip = i == 0 && j == 0;
            let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
            let chrom_len = reference.chrom_length(cluster.chrom_id) as usize;
            let ref_seq = reference.get_seq(cluster.chrom_id, 0, chrom_len);

            let (alignment, ref_start_adjustment, seq_start, seq_end) =
                cluster.clone().into_alignment(
                    left_budget,
                    right_budget,
                    soft_clip,
                    strand_seq,
                    ref_seq,
                    &mut aligner,
                );

            // Update effective range: claim the full budget for this cluster.
            // Even if x-drop didn't extend fully, the budget is consumed so
            // neighboring clusters don't re-extend into the same gap.
            effective_ranges[pi] = (
                effective_ranges[pi].0.saturating_sub(left_budget),
                (effective_ranges[pi].1 + right_budget).min(seq_len),
            );

            results[j] = Some(ClusterResult {
                alignment,
                ref_start_adjustment,
                seq_start,
                seq_end,
            });
        }

        // Emit alignments in original set order (j = 0 is primary, j > 0 supplementary).
        for (j, cluster) in set.iter().enumerate() {
            let result = results[j]
                .as_ref()
                .expect("all clusters should have results");

            let mut flags = Flags::empty();
            if cluster.is_reverse {
                flags |= Flags::REVERSE_COMPLEMENTED;
            }
            if i > 0 {
                flags |= Flags::SECONDARY;
            }
            if j > 0 {
                flags |= Flags::SUPPLEMENTARY;
            }

            let (mate_ref_id, mate_pos) = if i == 0 && j > 0 {
                let primary_cluster = &set[0];
                (
                    Some(primary_cluster.chrom_id),
                    Some(primary_cluster.ref_start() + 1),
                )
            } else {
                (None, None)
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

            let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
            let strand_qual = if cluster.is_reverse { &rc_qual } else { qual };
            let chrom_len = reference.chrom_length(cluster.chrom_id) as usize;
            let ref_seq = reference.get_seq(cluster.chrom_id, 0, chrom_len);

            let seq_segment = &strand_seq[result.seq_start..result.seq_end];
            let qual_segment = &strand_qual[result.seq_start..result.seq_end];

            let cigar_str = result.alignment.cigar_string();
            let noodles_cigar: noodles::sam::alignment::record_buf::Cigar =
                result.alignment.cigar.iter().copied().collect();

            // Adjust reference position: the left extension consumes reference bases
            // before the first seed, so subtract the adjustment
            let ref_pos = cluster
                .ref_start()
                .saturating_sub(result.ref_start_adjustment)
                + 1;

            // Validate the alignment against sequences
            let query_start = result.alignment.leading_hard_clip();
            let ref_slice = &ref_seq[ref_pos - 1..]; // ref_pos is 1-based
            if let Err(err) = result
                .alignment
                .validate(ref_slice, strand_seq, query_start)
            {
                log::error!(
                    "Alignment validation failed for read {}: {} | \
                     cluster.read_range={:?}, left_budget={}, right_budget={}, \
                     ref_pos={}, strand_seq.len={}, seq_range={}..{}, cigar={}",
                    read_name,
                    err,
                    cluster.read_range(),
                    0, // budget info not readily available here
                    0,
                    ref_pos,
                    strand_seq.len(),
                    result.seq_start,
                    result.seq_end,
                    cigar_str
                );
                panic!("Alignment validation failed");
            }

            // Build auxiliary data (tags)
            let data: Data = [
                (Tag::try_from(*b"mc").unwrap(), Value::from(mc as i32)),
                (
                    Tag::try_from(*b"SA").unwrap(),
                    Value::from(summary.as_str()),
                ),
            ]
            .into_iter()
            .collect();

            let record = build_record(
                read_name,
                flags,
                cluster.chrom_id,
                ref_pos,
                mapqs[i][j] as u8,
                noodles_cigar,
                mate_ref_id,
                mate_pos,
                seq_segment,
                qual_segment,
                data,
            );
            writer.write_record(&record).expect("write failed");
        }
    }

    metrics::histogram!("analysis_alignment").record(alignment_start.elapsed().as_secs_f64());

    Ok(())
}

type SegmentSet = (RangeSet, Vec<usize>); // (covered read segments, cluster indices)

struct SegmentSetHeap<'a> {
    clusters: &'a [SeedCluster],
    params: AlignParams,
    read_len: usize,
}

impl<'a> Heapable for SegmentSetHeap<'a> {
    type Item = SegmentSet;

    const ORDERING: HeapOrdering = HeapOrdering::Max;

    fn cmp(&self, lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
        let l = score_clusters(
            lhs.1.iter().map(|&i| &self.clusters[i]),
            self.read_len,
            &self.params,
        );
        let r = score_clusters(
            rhs.1.iter().map(|&i| &self.clusters[i]),
            self.read_len,
            &self.params,
        );
        l.partial_cmp(&r).unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Score a segment set: sum of cluster qualities minus gap penalties for
/// all uncovered read regions — leading, internal, and trailing. Each gap
/// is scored as a deletion Op, using the same scoring model as
/// alignment gaps. This ensures that any changes to the gap scoring model
/// (e.g. non-linear penalties) are automatically reflected here.
fn score_clusters<'a>(
    clusters: impl Iterator<Item = &'a SeedCluster>,
    read_len: usize,
    params: &AlignParams,
) -> f64 {
    let mut cluster_score: f64 = 0.0;
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    for cluster in clusters {
        cluster_score += cluster.quality(params).value();
        ranges.push(cluster.fwd_read_range(read_len));
    }
    ranges.sort_unstable();
    let mut gap_penalty: f64 = 0.0;
    // Leading uncovered region: [0, first_start)
    if let Some(&(first_start, _)) = ranges.first() {
        if first_start > 0 {
            gap_penalty += params.quality(Op::new(Kind::Deletion, first_start)).value();
        }
    }
    // Internal gaps between clusters
    for pair in ranges.windows(2) {
        let gap_len = pair[1].0.saturating_sub(pair[0].1);
        if gap_len > 0 {
            gap_penalty += params.quality(Op::new(Kind::Deletion, gap_len)).value();
        }
    }
    // Trailing uncovered region: [last_end, read_len)
    if let Some(&(_, last_end)) = ranges.last() {
        if last_end < read_len {
            gap_penalty += params
                .quality(Op::new(Kind::Deletion, read_len - last_end))
                .value();
        }
    }
    cluster_score + gap_penalty // gap_penalty is already negative
}

fn form_covering_sets(
    clusters: &[SeedCluster],
    read_name: &str,
    read_len: usize,
) -> Vec<Vec<SeedCluster>> {
    let mut order_by_quality: Vec<usize> = (0..clusters.len()).collect();
    let params = AlignParams::default();
    order_by_quality.sort_by_key(|i| OrderedFloat(-clusters[*i].quality(&params).value()));

    let mut segment_set_heap = Heap::new(SegmentSetHeap {
        clusters,
        params,
        read_len,
    });
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

    log::debug!(
        "Read {}: assigned {} clusters to {} segment sets",
        read_name,
        clusters.len(),
        segment_set_heap.len(),
    );

    let segment_sets: Vec<Vec<SeedCluster>> = segment_set_heap
        .drain()
        .enumerate()
        .map(|(set_idx, (_, set))| {
            let score = score_clusters(set.iter().map(|&i| &clusters[i]), read_len, &params);
            log::debug!(
                "  {}: Set {}: score {:.2}, clusters [{}]",
                read_name,
                set_idx,
                score,
                set.iter()
                    .map(|&i| {
                        let c = &clusters[i];
                        let (s, e) = c.fwd_read_range(read_len);
                        format!("{}({}-{},q={:.0})", i, s, e, c.quality(&params).value())
                    })
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            set.iter().map(|&i| clusters[i].clone()).collect()
        })
        .collect();

    segment_sets
}

/// A read to be processed by a worker thread
struct ReadWork {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Process reads from a FASTQ or unaligned BAM file using multiple threads.
///
/// The input format is auto-detected from the file extension:
/// - `.bam` → unaligned BAM
/// - anything else → FASTQ (with optional gzip/bzip2/xz compression)
///
/// Reads are distributed to worker threads via a channel. The InMemoryReference
/// is shared across all threads via Arc (no per-thread cloning needed).
pub fn process_reads_parallel<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    reads: &str,
    sam: Option<&str>,
    num_threads: usize,
    command_line: &str,
    read_group_header: Option<&str>,
    output_format: OutputFormat,
) -> Result<()> {
    use crossbeam::channel::bounded;

    let format = detect_input_format(reads);
    log::info!(
        "Processing reads from {} ({}) using {} threads, output format: {}",
        reads,
        match format {
            InputFormat::Fastq => "FASTQ",
            InputFormat::Bam => "BAM",
        },
        num_threads,
        output_format,
    );

    let now = std::time::Instant::now();
    let mut num_records = 0;

    // Store reference chromosome info for debug SAM headers
    debug::set_reference_info(reference.chromosomes());

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
        AlignmentWriter::builder(output, output_format, reference.to_fasta_repository())
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

        match format {
            InputFormat::Fastq => {
                // Read FASTQ and send to workers
                // Use niffler for transparent decompression (gzip, bzip2, xz)
                let (decompressed_reader, compression) =
                    niffler::from_path(std::path::Path::new(reads))
                        .expect("Failed to open FASTQ file");
                if compression != niffler::Format::No {
                    log::info!("Detected {:?} compression", compression);
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
            }
            InputFormat::Bam => {
                // Read unaligned BAM and send to workers
                let file = std::fs::File::open(reads).expect("Failed to open BAM file");
                let mut reader = noodles::bam::io::Reader::new(file);
                let header = reader.read_header().expect("Failed to read BAM header");

                let mut rc_buf = Vec::new();

                for result in reader.record_bufs(&header) {
                    let record = result.expect("Failed to read BAM record");

                    let raw_seq: Vec<u8> = record.sequence().as_ref().iter().cloned().collect();
                    if raw_seq.is_empty() {
                        continue;
                    }

                    let name = record
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                        .unwrap_or_else(|| format!("unnamed_{}", num_records));

                    let is_reverse = record.flags().is_reverse_complemented();
                    let raw_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

                    let (seq, qual) = if is_reverse {
                        // Undo reverse complement applied by previous aligner:
                        // reverse-complement the sequence and reverse the quality.
                        reverse_complement_into(&raw_seq, &mut rc_buf);
                        let seq = rc_buf.clone();
                        let qual: Vec<u8> = raw_qual
                            .iter()
                            .rev()
                            .map(|&q| q.saturating_add(33))
                            .collect();
                        (seq, qual)
                    } else {
                        // Convert quality from raw Phred (BAM) to Phred+33 (SAM/FASTQ).
                        let qual: Vec<u8> =
                            raw_qual.iter().map(|&q| q.saturating_add(33)).collect();
                        (raw_seq, qual)
                    };

                    // Handle missing quality scores (all 0xFF in BAM → empty after decode)
                    let qual = if qual.is_empty() {
                        vec![b'!'; seq.len()] // Phred 0 + 33 = '!' as placeholder
                    } else {
                        qual
                    };

                    let work = ReadWork { name, seq, qual };
                    sender.send(work).expect("Failed to send work to thread");
                    num_records += 1;
                }
            }
        }

        // Signal completion by dropping sender
        drop(sender);

        // Scoped threads automatically join when scope ends
    })
    .expect("Scoped thread panicked");

    writer.finish()?;

    // Finish all debug files
    DebugFile::<SeedsSamDebug>::finish_all();

    let elapsed = now.elapsed();
    log::info!(
        "Completed processing {} reads from {} in {:.2?}",
        num_records,
        reads,
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
        use crate::align::{Alignment, Kind, Op};
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
                cigar: vec![Op::new(Kind::SequenceMatch, 10)],
            },
            Alignment {
                divergence: DivergenceScore::new(20.0),
                cigar: vec![Op::new(Kind::SequenceMatch, 20)],
            },
            Alignment {
                divergence: DivergenceScore::new(15.0),
                cigar: vec![Op::new(Kind::SequenceMatch, 15)],
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
        use crate::align::{Alignment, Kind, Op};
        use crate::scores::DivergenceScore;

        let seeds = vec![make_hit(0, 100, 0, 20), make_hit(0, 300, 100, 20)];

        let mut cluster = SeedCluster::new(seeds, true, 1).unwrap();

        // Add dummy gap alignment (1 gap for 2 seeds)
        cluster.gap_alignments = vec![Alignment {
            divergence: DivergenceScore::new(10.0),
            cigar: vec![Op::new(Kind::SequenceMatch, 10)],
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
