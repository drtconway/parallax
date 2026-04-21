#[cfg(feature = "conventional")]
use core::panic;
#[cfg(feature = "conventional")]
use std::io::Write;
use std::{sync::Arc, usize};

use noodles::sam::alignment::{
    record::Flags,
    record::data::field::Tag,
    record_buf::{Data, data::field::Value},
};
#[cfg(feature = "conventional")]
use ordered_float::OrderedFloat;

use crate::{
    align::{AlignParams, Aligner, Alignment, Kind, Op},
    config::{self},
    error::Result,
    index::{Index, decode_locus},
    kmers::Kmer,
    reads::{
        builder::{build_record, build_unmapped_record},
        seeds::SeedHit,
    },
    reference::InMemoryReference,
    utils::{
        debug::{self, DebugFile, DebugOutput, DebugTsvWriter, TsvRow},
        hasher::FnvHasher,
        sequence::{complement, reverse_complement_into},
    },
    writer::{AlignmentWriter, OutputFormat},
};

#[cfg(feature = "explanatory")]
use crate::reads::extended::ExtendedSeed;
#[cfg(feature = "conventional")]
use crate::reads::seeds::analyze_gap_fills;
#[cfg(feature = "conventional")]
use crate::reads::seeds::seed_cluster::SeedCluster;
#[cfg(feature = "conventional")]
use crate::scores::compute_mapq_from_diff;
#[cfg(feature = "conventional")]
use crate::utils::GroupByTrait;
#[cfg(feature = "conventional")]
use crate::utils::heap::{Heap, HeapOrdering, Heapable};
#[cfg(feature = "conventional")]
use crate::utils::range_set::RangeSet;

pub mod builder;
#[cfg(feature = "conventional")]
pub mod chains;
pub mod seeds;

// ── Debug file statics ──────────────────────────────────────────────────────

/// Debug SAM file with extended seeds (before clustering).
static SEEDS_SAM: DebugFile<SeedsSamDebug> = DebugFile::new();

/// Debug TSV file with candidate seeds.
static SEEDS_TSV: DebugFile<SeedsTsvDebug> = DebugFile::new();

/// Debug TSV file with seed chains (after chaining, before alignment).
#[cfg(feature = "conventional")]
static CHAINS_TSV: DebugFile<ChainsTsvDebug> = DebugFile::new();

/// Debug SAM file with seed-level WIS groupings.
#[cfg(feature = "conventional")]
static WIS_SAM: DebugFile<WisSamDebug> = DebugFile::new();

// ── Concrete debug types ─────────────────────────────────────────────────────

pub struct SeedsSamDebug(DebugTsvWriter);

impl DebugOutput for SeedsSamDebug {
    type Item<'a> = str;
    fn create() -> Option<Self> {
        let path = &config::get().seeding.debug_seeds_sam;
        if path.is_empty() {
            return None;
        }
        DebugTsvWriter::open(path, debug::sam_header().as_deref())
            .ok()
            .map(Self)
    }
    fn append(&self, item: &str) {
        self.0.append(item);
    }
    fn finish(&self) {
        self.0.finish();
    }
}

type SeedsTsvRow<'a> = (
    &'a str,
    usize,
    usize,
    usize,
    &'a str,
    usize,
    usize,
    &'a str,
    usize,
);

pub(crate) struct SeedsTsvDebug(DebugTsvWriter);

impl SeedsTsvDebug {
    const HEADERS: &[&str] = &[
        "read_name",
        "read_start",
        "read_end",
        "read_len",
        "chrom",
        "ref_start",
        "ref_end",
        "strand",
        "score",
    ];
    const _CHECK: () = assert!(Self::HEADERS.len() == <SeedsTsvRow<'static> as TsvRow>::NUM_FIELDS);
}

impl DebugOutput for SeedsTsvDebug {
    type Item<'a> = SeedsTsvRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_seeds_tsv;
        if path.is_empty() {
            return None;
        }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &SeedsTsvRow<'_>) {
        self.0.append_row(item);
    }
    fn finish(&self) {
        self.0.finish();
    }
}

#[cfg(feature = "conventional")]
type ChainsTsvRow<'a> = (
    &'a str,
    usize,
    &'a str,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    &'a str,
    &'a str,
    u32,
);

#[cfg(feature = "conventional")]
pub(crate) struct ChainsTsvDebug(DebugTsvWriter);

#[cfg(feature = "conventional")]
impl ChainsTsvDebug {
    const HEADERS: &[&str] = &[
        "read_name",
        "cluster_id",
        "row_type",
        "read_start",
        "read_end",
        "read_width",
        "ref_start",
        "ref_end",
        "ref_width",
        "chrom",
        "strand",
        "uniqueness",
    ];
    const _CHECK: () =
        assert!(Self::HEADERS.len() == <ChainsTsvRow<'static> as TsvRow>::NUM_FIELDS);
}

#[cfg(feature = "conventional")]
pub(crate) struct WisSamDebug(DebugTsvWriter);

#[cfg(feature = "conventional")]
impl DebugOutput for WisSamDebug {
    type Item<'a> = str;
    fn create() -> Option<Self> {
        let path = &config::get().seeding.debug_wis_sam;
        if path.is_empty() {
            return None;
        }
        DebugTsvWriter::open(path, debug::sam_header().as_deref())
            .ok()
            .map(Self)
    }
    fn append(&self, item: &str) {
        self.0.append(item);
    }
    fn finish(&self) {
        self.0.finish();
    }
}

#[cfg(feature = "conventional")]
impl DebugOutput for ChainsTsvDebug {
    type Item<'a> = ChainsTsvRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_chains_tsv;
        if path.is_empty() {
            return None;
        }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &ChainsTsvRow<'_>) {
        self.0.append_row(item);
    }
    fn finish(&self) {
        self.0.finish();
    }
}

#[derive(Debug)]
enum AlignmentError {
    #[allow(dead_code)]
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
    /// Scratch space for merging/deduplication
    merge_scratch: Vec<SeedHit>,
    /// Batch buffer for prefetched lookups: (read_pos, kmer_value)
    kmer_batch: Vec<(usize, u64)>,
    /// Deferred mid-frequency seeds: (read_pos, kmer_value, hit_count).
    /// Collected during Phase 1 and selectively rescued into gaps after
    /// merge+extend.
    deferred_seeds: Vec<(usize, u64, u32)>,
}

impl ClusterCollector {
    /// Create a new collector with empty buffers
    fn new() -> Self {
        ClusterCollector {
            hits: Vec::new(),
            merge_scratch: Vec::new(),
            kmer_batch: Vec::new(),
            deferred_seeds: Vec::new(),
        }
    }

    /// Collect seed clusters from a single strand.
    ///
    /// This performs seeding, merging, extension, and DBSCAN clustering, returning
    /// the resulting clusters without building alignments. This separation allows
    /// for cross-strand analysis before alignment construction.

    /// Sort, merge adjacent seeds on the same diagonal, extend exact matches,
    /// and remove duplicates.
    ///
    /// This is the core seed-consolidation pipeline (Phases 2–3c) used after
    /// initial seed collection and again after rescue. It operates in-place on
    /// `self.hits`, using `self.merge_scratch` as temporary storage.
    fn sort_merge_extend<const K: usize>(
        &mut self,
        strand_seq: &[u8],
        reference: &InMemoryReference,
    ) {
        // Sort: SeedHit's Ord gives (chrom_id, diagonal, ref_pos) order
        self.hits.sort_unstable();

        // Merge overlapping/adjacent hits on the same diagonal
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

        // Extend each seed's exact match bidirectionally
        for hit in self.hits.iter_mut() {
            let ref_seq = reference.get_seq(hit.chrom_id, 0, usize::MAX);
            hit.extend_exact(strand_seq, ref_seq);
        }

        // Remove duplicates created by extension
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
    }

    /// Rescue deferred mid-frequency seeds into gaps in the seed layout.
    ///
    /// After the initial merge+extend pass, scans `self.hits` for read-position
    /// gaps larger than `rescue_spacing`. For each gap, selects the deferred seed
    /// with the lowest hit count, re-queries the index, and adds the resulting
    /// hits. Rate-limited to one rescue per `rescue_spacing` bp of gap.
    ///
    /// Returns the number of seeds rescued.
    fn rescue_seeds<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        index: &Index<K, S>,
        reference: &InMemoryReference,
        rescue_spacing: usize,
    ) -> usize {
        if self.deferred_seeds.is_empty() || rescue_spacing == 0 {
            return 0;
        }
        let seq_len = strand_seq.len();

        // Sort deferred seeds by read position for binary search
        self.deferred_seeds.sort_unstable_by_key(|&(pos, _, _)| pos);

        // Build a sorted list of read coverage intervals from current seeds
        // (hits are sorted by chrom_id/diagonal, but we need read_pos order)
        let mut read_positions: Vec<(usize, usize)> = self
            .hits
            .iter()
            .map(|h| (h.read_pos, h.read_end()))
            .collect();
        read_positions.sort_unstable();

        // Merge overlapping read intervals to get coverage
        let mut coverage: Vec<(usize, usize)> = Vec::new();
        for (start, end) in read_positions {
            if let Some(last) = coverage.last_mut() {
                if start <= last.1 {
                    last.1 = last.1.max(end);
                    continue;
                }
            }
            coverage.push((start, end));
        }

        // Find gaps and rescue into them
        let mut rescued = 0usize;
        let mut prev_end = 0usize;
        // Add a sentinel for the end of the read
        let gap_iter = coverage
            .iter()
            .map(|&(s, e)| (s, e))
            .chain(std::iter::once((seq_len, seq_len)));

        for (gap_right, _) in gap_iter {
            let gap_left = prev_end;
            let gap_size = gap_right.saturating_sub(gap_left);

            if gap_size >= rescue_spacing {
                // How many rescues allowed in this gap?
                let max_rescues = gap_size / rescue_spacing;

                // Find deferred seeds that fall in this gap, sorted by hit_count
                let lo = self
                    .deferred_seeds
                    .partition_point(|&(pos, _, _)| pos < gap_left);
                let hi = self
                    .deferred_seeds
                    .partition_point(|&(pos, _, _)| pos < gap_right);

                if lo < hi {
                    // Collect candidates in this gap, pick lowest hit_count first
                    let mut candidates: Vec<_> = self.deferred_seeds[lo..hi].to_vec();
                    candidates.sort_unstable_by_key(|&(_, _, count)| count);

                    let mut rescues_in_gap = 0;
                    let mut last_rescue_pos: Option<usize> = None;

                    for (read_pos, kmer_val, hit_count) in candidates {
                        // Rate-limit: enforce spacing between rescues
                        if let Some(lrp) = last_rescue_pos {
                            if read_pos.abs_diff(lrp) < rescue_spacing {
                                continue;
                            }
                        }
                        if rescues_in_gap >= max_rescues {
                            break;
                        }

                        // Re-query the index for this kmer and decode locations
                        let kmer = Kmer::<K>(kmer_val);
                        index.with(&kmer, |_count, loci| {
                            for &loc in loci {
                                let (chrom_id, chrom_pos) = decode_locus(loc);
                                self.hits.push(SeedHit::new(
                                    chrom_id, chrom_pos, read_pos, kmer_val, hit_count, K,
                                ));
                            }
                        });

                        last_rescue_pos = Some(read_pos);
                        rescues_in_gap += 1;
                        rescued += 1;
                    }
                }
            }

            // Advance prev_end past this coverage interval
            if let Some(last) = coverage.iter().find(|&&(s, _)| s == gap_right) {
                prev_end = last.1;
            }
        }

        if rescued > 0 {
            self.sort_merge_extend::<K>(strand_seq, reference);
        }

        rescued
    }
    #[cfg(feature = "conventional")]
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
            #[cfg(all(feature = "chainer-fenwick", not(feature = "chainer-kruskal")))]
            let chrom_clusters =
                chains::fenwick::collect_chains(&mut seeds, &chrom_name, is_reverse);
            #[cfg(any(feature = "chainer-kruskal", not(feature = "chainer-fenwick")))]
            let chrom_clusters =
                chains::kruskal::collect_chains(&mut seeds, &chrom_name, is_reverse);
            clusters.extend(chrom_clusters);
        }

        clusters
    }

    #[cfg(feature = "conventional")]
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
        self.deferred_seeds.clear();
        let mid_occ = cfg.seeding.mid_seed_occurrences as u32;
        let max_occ = cfg.seeding.max_seed_occurrences as u32;

        // Phase 1: Collect seed hits using forward-only syncmers
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            strand_seq,
            [(); S],
            |pos, kmer| {
                index.with(&kmer, |hit_count, loci| {
                    if mid_occ > 0 && hit_count > mid_occ && hit_count <= max_occ {
                        // Mid-frequency: defer for potential rescue
                        self.deferred_seeds.push((pos, kmer.0, hit_count));
                    } else if hit_count <= max_occ {
                        // Low-frequency (or rescue disabled): collect immediately
                        for &loc in loci {
                            let (chrom_id, chrom_pos) = decode_locus(loc);
                            self.hits
                                .push(SeedHit::new(chrom_id, chrom_pos, pos, kmer.0, hit_count, K));
                        }
                    }
                    // hit_count > max_occ: skip entirely
                });
            },
        );

        let strand_name = if is_reverse { "REV" } else { "FWD" };
        metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase()))
            .record(self.hits.len() as f64);

        // Phases 2–3c: Sort, merge, extend, dedup
        self.sort_merge_extend::<K>(strand_seq, reference);

        // Phase 3d: Rescue deferred mid-frequency seeds into coverage gaps
        let rescued =
            self.rescue_seeds::<K, S>(strand_seq, index, reference, cfg.seeding.rescue_spacing);
        if rescued > 0 {
            log::info!("{read_name} {strand_name}: rescued {rescued} deferred seeds into gaps");
        }

        // Write debug SAM output for seed hits
        if SEEDS_SAM.is_enabled() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                SEEDS_SAM.append(&hit.to_sam_line(
                    read_name,
                    chrom_name,
                    is_reverse,
                    strand_seq,
                    strand_qual,
                ));
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
                    read_name,
                    fwd_start,
                    fwd_end,
                    seq_len,
                    chrom_name,
                    hit.ref_pos,
                    hit.ref_end(),
                    strand,
                    hit.match_len,
                ));
            }
        }

        // Dead debug code — kept for occasional manual use
        #[cfg(any(feature = "conventional", feature = "explanatory"))]
        if false {
            use crate::reads::seeds::{Read, SeedSaver};
            use std::io::Write as _;
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
        self.deferred_seeds.clear();

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
        let mid_occ = cfg.seeding.mid_seed_occurrences as u32;
        index.lookup_batch(&self.kmer_batch, |read_pos, kmer_val, hit_count, loci| {
            if mid_occ > 0 && hit_count > mid_occ && hit_count <= max_occ {
                // Mid-frequency: defer for potential rescue
                self.deferred_seeds.push((read_pos, kmer_val, hit_count));
            } else if hit_count <= max_occ {
                // Low-frequency (or rescue disabled): collect immediately
                for &loc in loci {
                    let (chrom_id, chrom_pos) = decode_locus(loc);
                    self.hits.push(SeedHit::new(
                        chrom_id, chrom_pos, read_pos, kmer_val, hit_count, K,
                    ));
                }
            }
            // hit_count > max_occ: skip entirely
        });

        let strand_name = if is_reverse { "REV" } else { "FWD" };
        metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase()))
            .record(self.hits.len() as f64);

        // Phases 2–3c: Sort, merge, extend, dedup
        self.sort_merge_extend::<K>(strand_seq, reference);

        // Phase 3d: Rescue deferred mid-frequency seeds into coverage gaps
        let rescued =
            self.rescue_seeds::<K, S>(strand_seq, index, reference, cfg.seeding.rescue_spacing);
        if rescued > 0 {
            let strand_name = if is_reverse { "REV" } else { "FWD" };
            log::debug!("{read_name} {strand_name}: rescued {rescued} deferred seeds into gaps");
        }

        // Write debug SAM output for seed hits
        if SEEDS_SAM.is_enabled() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                SEEDS_SAM.append(&hit.to_sam_line(
                    read_name,
                    chrom_name,
                    is_reverse,
                    strand_seq,
                    strand_qual,
                ));
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
                    read_name,
                    fwd_start,
                    fwd_end,
                    seq_len,
                    chrom_name,
                    hit.ref_pos,
                    hit.ref_end(),
                    strand,
                    hit.match_len,
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
    aligner: &mut Aligner,
) {
    #[cfg(feature = "explanatory")]
    {
        align_read_inner(
            index,
            reference,
            writer,
            read_name,
            seq,
            qual,
            alignment_params,
            aligner,
        )
        .expect("alignment failed");
    }

    #[cfg(not(feature = "explanatory"))]
    {
        match align_read_inner(
            index,
            reference,
            writer,
            read_name,
            seq,
            qual,
            alignment_params,
            aligner,
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
}

#[cfg(feature = "explanatory")]
pub mod extended;

#[cfg(feature = "explanatory")]
fn align_read_inner<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    writer: &AlignmentWriter,
    read_name: &str,
    seq: &[u8],
    qual: &[u8],
    alignment_params: &AlignParams,
    aligner: &mut Aligner, // reused across all reads on this thread
) -> std::result::Result<(), AlignmentError> {
    // Note: alignment_params is currently unused in this function, but we include it in the signature
    // because it is used in the conventional alignment pipeline and we want to keep the signatures similar for easier comparison.
    let _ = alignment_params;

    let seq_len = seq.len();

    // Reusable cluster collector
    let mut collector = ClusterCollector::new();

    // Compute reverse complement for reverse strand processing
    let mut rc_seq = Vec::with_capacity(seq_len);
    reverse_complement_into(seq, &mut rc_seq);

    // Reverse quality scores for reverse strand (if available)
    let rc_qual: Vec<u8> = qual.iter().rev().copied().collect();

    // =========================================================================
    // Phase 1: Collect all seeds from both strands
    // =========================================================================

    let mut all_seeds: Vec<ExtendedSeed> = Vec::new();

    collector.gather_seeds_batched::<K, S>(seq, qual, false, index, reference, read_name);
    all_seeds.extend(
        collector
            .hits
            .iter()
            .map(|seed| ExtendedSeed::from_seed_hit(seed, false, seq_len)),
    );
    collector.gather_seeds_batched::<K, S>(&rc_seq, &rc_qual, true, index, reference, read_name);
    all_seeds.extend(
        collector
            .hits
            .iter()
            .map(|seed| ExtendedSeed::from_seed_hit(seed, true, seq_len)),
    );

    // Simplify seeds by merging overlapping ones on the same diagonal
    ExtendedSeed::simplify_seeds(&mut all_seeds);

    let mut groups = ExtendedSeed::form_explanatory_groups(&all_seeds);

    if groups.is_empty() {
        let record = build_unmapped_record(read_name, seq, qual);
        writer.write_record(&record).expect("write failed");
        return Ok(());
    }

    // Assemble segments: each segment is a maximal run of colinear seeds.
    // A None gap (or end of group) terminates the current segment.
    struct Segment {
        first_seed: usize,
        last_seed: usize,
        alignment: Alignment,
    }

    let mut explanations: Vec<Vec<Segment>> = Vec::new();

    for i in 0..groups.len() {
        let group = &mut groups[i];

        if false {
            let total_weight: f64 = group.iter().map(|s| s.weight()).sum();
            let max_weight = group.iter().map(|s| s.weight()).fold(0. / 0., f64::max);
            let max_length: usize = group.iter().map(|s| s.length()).max().unwrap_or(0);
            let min_multiplicity: usize = group.iter().map(|s| s.multiplicity()).min().unwrap_or(0);
            let max_multiplicity: usize = group.iter().map(|s| s.multiplicity()).max().unwrap_or(0);
            println!(
                "Group {}: total weight {:.1}, max weight {:.1}, max length {}, min multiplicity {}, max multiplicity {}, {} seeds",
                i,
                total_weight,
                max_weight,
                max_length,
                min_multiplicity,
                max_multiplicity,
                group.len()
            );
        }

        ExtendedSeed::extend_and_trim(group, seq, reference);

        let gaps = ExtendedSeed::align_gaps(group, seq, reference, aligner);

        let n = group.len();
        let mut segments: Vec<Segment> = Vec::new();
        let mut current_parts: Vec<Alignment> = Vec::new();
        let mut segment_start = 0;
        for j in 0..n {
            if current_parts.is_empty() {
                segment_start = j;
            }
            current_parts.push(group[j].to_alignment());
            match gaps.get(j) {
                Some(Some(aln)) => {
                    current_parts.push(aln.clone());
                }
                None | Some(None) => {
                    segments.push(Segment {
                        first_seed: segment_start,
                        last_seed: j,
                        alignment: Alignment::concat(&std::mem::take(&mut current_parts)),
                    });
                }
            }
        }
        explanations.push(segments);
    }

    if false {
        for (i, segmentss) in explanations.iter().enumerate() {
            let mut query_coverage = 0usize;
            let mut total_score = 0.0f64;
            for segment in segmentss.iter() {
                query_coverage += segment.alignment.query_length();
                total_score += segment.alignment.divergence.0;
            }
            // Treat the uncovered portion of the query as a deletion
            let missing_coverage = seq_len - query_coverage;
            total_score += missing_coverage as f64;
            let coverage_pct = 100.0 * (query_coverage as f64) / (seq_len as f64);
            println!(
                "Group {}: {} segments, total score {}, query coverage {:.1}%",
                i,
                segmentss.len(),
                total_score,
                coverage_pct
            );
        }
    }

    for (i, segments) in explanations.iter().enumerate() {
        let group = &groups[i];

        // Build SA tag summaries for each segment so we can cross-reference.
        // Format per SAM spec: rname,pos,strand,CIGAR,mapQ,NM
        let sa_entries: Vec<String> = segments
            .iter()
            .map(|segment| {
                let first = &group[segment.first_seed];
                let last = &group[segment.last_seed];
                let is_reverse = first.is_reverse();
                let chrom_id = first.ref_chrom_id();
                let chrom_name = reference.chrom_name(chrom_id);

                let ref_pos = if is_reverse {
                    last.ref_start() + 1
                } else {
                    first.ref_start() + 1
                };
                let strand = if is_reverse { "-" } else { "+" };
                let summary_cigar = segment.alignment.summary_cigar(
                    first.read_start(),
                    last.read_end(),
                    seq_len,
                    is_reverse,
                );
                let nm = segment.alignment.mismatch_count();

                format!(
                    "{},{},{},{},255,{}",
                    chrom_name, ref_pos, strand, summary_cigar, nm
                )
            })
            .collect();

        // Pick the best segment (longest query span) as the representative.
        let best_seg_idx = segments
            .iter()
            .enumerate()
            .max_by_key(|(_, seg)| seg.alignment.query_length())
            .map(|(idx, _)| idx)
            .unwrap_or(0);

        for (seg_idx, segment) in segments.iter().enumerate() {
            let first = &group[segment.first_seed];
            let last = &group[segment.last_seed];
            let is_reverse = first.is_reverse();
            let chrom_id = first.ref_chrom_id();

            // SAM POS: leftmost reference position (1-based).
            // Forward: first seed has the leftmost ref position.
            // Reverse: last seed has the leftmost ref position (ref
            // decreases as read advances in a colinear reverse segment).
            let ref_pos = if is_reverse {
                last.ref_start() + 1
            } else {
                first.ref_start() + 1
            };

            // Read range covered by this segment.
            let seg_read_start = first.read_start();
            let seg_read_end = last.read_end();

            // Validate the segment alignment against the reference and query.
            if false {
                let (ref_begin, ref_end) = if is_reverse {
                    (last.ref_start(), first.ref_start() + first.length())
                } else {
                    (first.ref_start(), last.ref_start() + last.length())
                };

                let ref_slice: Vec<u8> = if is_reverse {
                    reference
                        .get_seq(chrom_id, ref_begin, ref_end)
                        .iter()
                        .rev()
                        .map(|&b| complement(b))
                        .collect()
                } else {
                    reference.get_seq(chrom_id, ref_begin, ref_end).to_vec()
                };

                // The alignment was built against seq (forward read),
                // regardless of strand.
                let query_seq = &seq[seg_read_start..seg_read_end];

                if let Err(e) = segment.alignment.validate(&ref_slice, query_seq, 0) {
                    let chrom_name = reference.chrom_name(chrom_id);
                    let strand = if is_reverse { "-" } else { "+" };
                    log::error!(
                        "VALIDATION FAILED: group {} seg {} ({} {}:{}-{} {}): {}",
                        i,
                        seg_idx,
                        read_name,
                        chrom_name,
                        ref_begin,
                        ref_end,
                        strand,
                        e
                    );
                }
            }

            // Flags
            let mut flags = Flags::empty();
            if is_reverse {
                flags |= Flags::REVERSE_COMPLEMENTED;
            }
            if i > 0 {
                flags |= Flags::SECONDARY;
            }
            if seg_idx != best_seg_idx {
                flags |= Flags::SUPPLEMENTARY;
            }

            let is_primary = i == 0 && seg_idx == best_seg_idx;

            // Build CIGAR: primary gets soft clips, secondary/supplementary
            // get hard clips (and a truncated SEQ/QUAL).
            let clip_kind = if is_primary {
                Kind::SoftClip
            } else {
                Kind::HardClip
            };
            let mut cigar = Vec::new();
            if is_reverse {
                // The alignment was built as seq vs rc(ref), but SAM
                // convention is rc_seq vs forward_ref.  We reverse the
                // CIGAR and use rc_seq coordinates for clipping:
                //   rc_start = seq_len - seg_read_end
                //   rc_end   = seq_len - seg_read_start
                let rc_start = seq_len - seg_read_end;
                let rc_end = seq_len - seg_read_start;
                if rc_start > 0 {
                    cigar.push(Op::new(clip_kind, rc_start));
                }
                for &op in segment.alignment.cigar.iter().rev() {
                    cigar.push(op);
                }
                if rc_end < seq_len {
                    cigar.push(Op::new(clip_kind, seq_len - rc_end));
                }
            } else {
                if seg_read_start > 0 {
                    cigar.push(Op::new(clip_kind, seg_read_start));
                }
                cigar.extend_from_slice(&segment.alignment.cigar);
                if seg_read_end < seq_len {
                    cigar.push(Op::new(clip_kind, seq_len - seg_read_end));
                }
            }
            let noodles_cigar: noodles::sam::alignment::record_buf::Cigar =
                cigar.iter().copied().collect();

            // SEQ/QUAL: primary emits the full read; secondary/supplementary
            // emit only the aligned portion.
            let (strand_seq, strand_qual) = if is_reverse {
                (&rc_seq[..], &rc_qual[..])
            } else {
                (seq, qual)
            };
            let (out_seq, out_qual) = if is_primary {
                (strand_seq, strand_qual)
            } else if is_reverse {
                // Non-primary reverse: use rc_seq coordinates.
                let rc_start = seq_len - seg_read_end;
                let rc_end = seq_len - seg_read_start;
                (&rc_seq[rc_start..rc_end], &rc_qual[rc_start..rc_end])
            } else {
                (
                    &strand_seq[seg_read_start..seg_read_end],
                    &strand_qual[seg_read_start..seg_read_end],
                )
            };

            // Build SA tag: list all OTHER segments in this group.
            let sa_value: String = sa_entries
                .iter()
                .enumerate()
                .filter(|&(k, _)| k != seg_idx)
                .map(|(_, entry)| entry.as_str())
                .collect::<Vec<_>>()
                .join(";");

            let data: Data = if segments.len() > 1 {
                vec![(
                    Tag::try_from(*b"SA").unwrap(),
                    Value::from(sa_value.as_str()),
                )]
                .into_iter()
                .collect()
            } else {
                Data::default()
            };

            let record = build_record(
                read_name,
                flags,
                chrom_id,
                ref_pos,
                255, // mapq placeholder
                noodles_cigar,
                None, // mate_ref_id
                None, // mate_pos
                out_seq,
                out_qual,
                data,
            );
            writer.write_record(&record).expect("write failed");
        }

        if i > 2 {
            break;
        }
    }

    Ok(())
}

#[cfg(feature = "conventional")]
fn align_read_inner<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    writer: &AlignmentWriter,
    read_name: &str,
    seq: &[u8],
    qual: &[u8],
    alignment_params: &AlignParams,
    mut aligner: &mut Aligner, // reused across all reads on this thread
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
    let stage_seeding = std::time::Instant::now();
    let fwd_clusters = collector.collect_from_strand(seq, qual, false, index, reference, read_name);
    all_clusters.extend(fwd_clusters);

    // Collect clusters from reverse strand
    let rev_clusters =
        collector.collect_from_strand(&rc_seq, &rc_qual, true, index, reference, read_name);
    all_clusters.extend(rev_clusters);
    metrics::histogram!("stage_seeding").record(stage_seeding.elapsed().as_secs_f64());

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
    // PASS 1.1: Merge colinear deletion clusters
    // =========================================================================
    // Clusters on the same chrom/strand that are colinear and adjacent on the
    // read likely span a deletion too large for the chainer's diagonal band.
    // Merge their seed lists so align_gaps() can bridge the ref gap with a D op.
    let pre_merge = all_clusters.len();
    merge_deletion_clusters(&mut all_clusters, seq_len, K / 2);
    if all_clusters.len() < pre_merge {
        log::info!(
            "Read {}: merged {} clusters into {} (deletion bridging)",
            read_name,
            pre_merge,
            all_clusters.len(),
        );
    }

    // =========================================================================
    // EXPERIMENTAL: Seed-level weighted interval scheduling (WIS)
    // =========================================================================
    // Flatten all seeds, run WIS on fwd-read intervals to find the best
    // non-overlapping seed subset, then group selected seeds into segments
    // by (chrom, strand, colinearity). Results are written to a debug SAM
    // file (one record per seed with XE/XS tags) and then discarded.
    if WIS_SAM.is_enabled() {
        construct_read_explanations(
            &all_clusters,
            seq_len,
            reference,
            read_name,
            seq,
            qual,
            &rc_seq,
            &rc_qual,
        );
    }

    // =========================================================================
    // Use estimated covering sets to select clusters for gap alignment
    // =========================================================================
    // Only gap-align clusters that appear in the top few estimated sets.
    // The covering set algorithm exhaustively assigns every cluster to some
    // set, but we only emit the primary + a handful of secondaries; clusters
    // relegated to low-scoring tail sets don't need expensive gap alignment.
    let estimated_sets_idx = form_covering_sets_estimated(&all_clusters, read_name, seq_len);
    let max_estimated_sets = 3; // primary + up to 2 secondaries for MAPQ
    let mut needed = vec![false; all_clusters.len()];
    for set in estimated_sets_idx.iter().take(max_estimated_sets) {
        for &i in set {
            needed[i] = true;
        }
    }
    let num_needed = needed.iter().filter(|&&b| b).count();
    let num_skipped = all_clusters.len() - num_needed;
    metrics::histogram!("est_clusters_needed").record(num_needed as f64);
    metrics::histogram!("est_clusters_skipped").record(num_skipped as f64);
    log::debug!(
        "Read {}: estimated selection: {} needed, {} skipped out of {} clusters",
        read_name,
        num_needed,
        num_skipped,
        all_clusters.len(),
    );

    // =========================================================================
    // PASS 1.5: Align gaps and split at failed alignments
    // =========================================================================
    // Only gap-align clusters selected by the estimated covering sets.
    // The estimated quality preserves relative ordering well enough (99.97%
    // agreement) that skipping unselected clusters is safe.

    let cfg = config::get();

    let stage_gap_align = std::time::Instant::now();
    let mut new_clusters = Vec::new();
    for (idx, cluster) in all_clusters.into_iter().enumerate() {
        if !needed[idx] {
            // Keep the cluster without gap-aligning it so it can still
            // participate in gap-fill analysis as a potential filler.
            new_clusters.push(cluster);
            continue;
        }
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
    metrics::histogram!("stage_gap_align").record(stage_gap_align.elapsed().as_secs_f64());

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
                        0u32,
                    ));
                }

                // Write seed row
                let read_width = seed.match_len;
                let ref_width = seed.match_len;
                CHAINS_TSV.append(&(
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
                    seed.kmer_uniqueness,
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
    // than bridging them with block aligner.

    let gap_fills = analyze_gap_fills(
        read_name,
        &all_clusters,
        seq_len,
        cfg.seeding.min_gap_for_split,
        2 * K,
        cfg.seeding.gap_fill_tolerance,
        alignment_params,
        &|id| reference.chrom_name(id),
    );

    if !gap_fills.is_empty() {
        log::info!(
            "Read {}: found {} gap fills for potential splitting",
            read_name,
            gap_fills.len(),
        );

        // Group splits by cluster, keeping the best filler per gap.
        // "Best" = highest filler cluster quality score.
        let mut best_by_gap: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        for fill in &gap_fills {
            let key = (fill.cluster_idx, fill.gap_seed_idx);
            let is_better = best_by_gap.get(&key).map_or(true, |&prev| {
                all_clusters[fill.filler_idx]
                    .quality(alignment_params)
                    .value()
                    > all_clusters[prev].quality(alignment_params).value()
            });
            if is_better {
                best_by_gap.insert(key, fill.filler_idx);
            }
        }

        // Collect into per-cluster lists sorted descending by gap_seed_idx
        // so we can split back-to-front without invalidating earlier indices.
        let mut splits_by_cluster: std::collections::HashMap<usize, Vec<(usize, usize)>> =
            std::collections::HashMap::new();
        for ((cluster_idx, gap_seed_idx), filler_idx) in best_by_gap {
            splits_by_cluster
                .entry(cluster_idx)
                .or_default()
                .push((gap_seed_idx, filler_idx));
        }
        for entries in splits_by_cluster.values_mut() {
            entries.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        }

        // Apply splits in descending cluster index order to preserve indices.
        // For each split, tag the head (filler follows) and tail (filler precedes)
        // with the filler's read/ref locus.
        let mut cluster_indices: Vec<_> = splits_by_cluster.keys().copied().collect();
        cluster_indices.sort_unstable_by(|a, b| b.cmp(a));

        for cluster_idx in cluster_indices {
            let entries = &splits_by_cluster[&cluster_idx];

            for &(gap_seed_idx, filler_idx) in entries {
                // Compute filler description before splitting (filler_idx is stable).
                let filler = &all_clusters[filler_idx];
                let (filler_read_start, filler_read_end) = filler.fwd_read_range(seq_len);
                let filler_chrom = reference.chrom_name(filler.chrom_id);
                let filler_strand = if filler.is_reverse { '-' } else { '+' };
                let ref_part = format!(
                    "{filler_chrom}:{}-{}{filler_strand}",
                    filler.ref_start(),
                    filler.ref_end(),
                );

                // Tag for the head piece: filler follows → star after read range
                let head_tag = format!("{filler_read_start}-{filler_read_end}*;{ref_part}",);
                // Tag for the tail piece: filler precedes → star before read range
                let tail_tag = format!("*{filler_read_start}-{filler_read_end};{ref_part}",);

                if let Some((new_cluster, _)) = all_clusters[cluster_idx].split_at_gap(gap_seed_idx)
                {
                    all_clusters[cluster_idx].split_fill_tags.push(head_tag);
                    let tail_idx = all_clusters.len();
                    all_clusters.push(new_cluster);
                    all_clusters[tail_idx].split_fill_tags.push(tail_tag);
                }
            }
        }

        // Re-sort after splitting
        all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));
    }

    // Remove clusters that were kept only for gap-fill analysis but were
    // never gap-aligned. They have multiple seeds but no gap_alignments,
    // so into_alignment() would produce invalid CIGARs.
    all_clusters.retain(|c| c.chain.len() < 2 || !c.gap_alignments.is_empty());

    let stage_covering = std::time::Instant::now();
    let segment_sets = form_covering_sets(&all_clusters, read_name, seq_len);

    let set_scores: Vec<f64> = segment_sets
        .iter()
        .map(|set| score_clusters(set.iter(), seq_len, alignment_params))
        .collect();
    metrics::histogram!("stage_covering_set").record(stage_covering.elapsed().as_secs_f64());

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
            if best_covering_score[k] >= set_scores[0] {
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

    let stage_emit = std::time::Instant::now();
    let emit_negative_primary = cfg.classification.emit_negative_primary;
    for (i, set) in segment_sets.into_iter().enumerate() {
        // Skip segment sets with non-positive scores.
        // These are noise — tiny repeat-derived clusters that can't
        // justify their existence against the uncovered-read penalty.
        if set_scores[i] <= 0.0 {
            if i > 0 {
                continue;
            }
            if !emit_negative_primary {
                return Err(AlignmentError::LowQuality);
            }
        }

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
            let mut tags: Vec<(Tag, Value)> = vec![
                (Tag::try_from(*b"mc").unwrap(), Value::from(mc as i32)),
                (
                    Tag::try_from(*b"SA").unwrap(),
                    Value::from(summary.as_str()),
                ),
            ];

            // If this cluster was split due to a gap-fill event, tag it with
            // the filler's read/ref locus so related segments can be verified.
            let xg_str;
            if !cluster.split_fill_tags.is_empty() {
                xg_str = cluster.split_fill_tags.join(",");
                tags.push((Tag::try_from(*b"XG").unwrap(), Value::from(xg_str.as_str())));
            }

            let data: Data = tags.into_iter().collect();

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

    metrics::histogram!("stage_emit").record(stage_emit.elapsed().as_secs_f64());
    metrics::histogram!("analysis_alignment").record(alignment_start.elapsed().as_secs_f64());

    Ok(())
}

#[cfg(feature = "conventional")]
type SegmentSet = (RangeSet, Vec<usize>, f64); // (covered read segments, cluster indices, cached score)

#[cfg(feature = "conventional")]
struct SegmentSetHeap;

#[cfg(feature = "conventional")]
impl Heapable for SegmentSetHeap {
    type Item = SegmentSet;

    const ORDERING: HeapOrdering = HeapOrdering::Max;

    fn cmp(&self, lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
        lhs.2
            .partial_cmp(&rhs.2)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Score a segment set: sum of cluster qualities minus gap penalties for
/// all uncovered read regions — leading, internal, and trailing. Each gap
/// is scored as a deletion Op, using the same scoring model as
/// alignment gaps. This ensures that any changes to the gap scoring model
/// (e.g. non-linear penalties) are automatically reflected here.
#[cfg(feature = "conventional")]
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

/// Flattened seed with provenance back to the source cluster.
#[cfg(feature = "conventional")]
struct FlatSeed {
    fwd_start: usize,
    fwd_end: usize,
    weight: f64,
    cluster_idx: usize,
    seed_idx: usize,
}

/// Run weighted interval scheduling on a subset of `flat_seeds` (given by
/// `indices` into the flat_seeds array, which must be sorted by fwd_end).
///
/// Returns the indices (into `flat_seeds`) of the selected seeds.
#[cfg(feature = "conventional")]
fn wis_select(flat_seeds: &[FlatSeed], indices: &[usize]) -> Vec<usize> {
    let m = indices.len();
    if m == 0 {
        return Vec::new();
    }

    let mut dp = vec![0.0f64; m];
    let mut best_by_end = vec![0.0f64; m];
    let mut take = vec![false; m];

    for i in 0..m {
        let ui = indices[i];
        let start = flat_seeds[ui].fwd_start;
        let w = flat_seeds[ui].weight;

        let prev_best = if start == 0 {
            0.0
        } else {
            let idx = indices[..i].partition_point(|&j| flat_seeds[j].fwd_end <= start);
            if idx == 0 { 0.0 } else { best_by_end[idx - 1] }
        };

        let take_score = w + prev_best;
        let skip_score = if i > 0 { best_by_end[i - 1] } else { 0.0 };

        if take_score >= skip_score {
            dp[i] = take_score;
            take[i] = true;
        } else {
            dp[i] = skip_score;
            take[i] = false;
        }
        best_by_end[i] = if i > 0 {
            dp[i].max(best_by_end[i - 1])
        } else {
            dp[i]
        };
    }

    // Backtrack
    let mut selected = Vec::new();
    let mut i = m;
    while i > 0 {
        i -= 1;
        if take[i] && dp[i] == best_by_end[i] {
            selected.push(indices[i]);
            let start = flat_seeds[indices[i]].fwd_start;
            while i > 0 && flat_seeds[indices[i - 1]].fwd_end > start {
                i -= 1;
            }
        }
    }
    selected.reverse();
    selected
}

/// Group a set of selected seed indices (sorted by fwd_start) into segments
/// by (chrom_id, is_reverse, reference colinearity).
///
/// Returns per-seed segment IDs and the next available segment number.
#[cfg(feature = "conventional")]
fn group_into_segments(
    selected: &[usize],
    flat_seeds: &[FlatSeed],
    clusters: &[SeedCluster],
    segment_id: &mut [usize],
    start_segment: usize,
) -> usize {
    let mut current = start_segment;
    for (pos, &si) in selected.iter().enumerate() {
        if pos == 0 {
            segment_id[si] = current;
            continue;
        }
        let prev_si = selected[pos - 1];
        let prev = &flat_seeds[prev_si];
        let curr = &flat_seeds[si];
        let prev_cluster = &clusters[prev.cluster_idx];
        let curr_cluster = &clusters[curr.cluster_idx];
        let prev_seed = &prev_cluster.chain[prev.seed_idx];
        let curr_seed = &curr_cluster.chain[curr.seed_idx];

        let same_chrom = prev_cluster.chrom_id == curr_cluster.chrom_id;
        let same_strand = prev_cluster.is_reverse == curr_cluster.is_reverse;
        let colinear = if same_chrom && same_strand {
            // As fwd-read position increases, ref positions increase for
            // fwd-strand seeds and decrease for rev-strand seeds.
            if curr_cluster.is_reverse {
                curr_seed.ref_pos < prev_seed.ref_pos
            } else {
                curr_seed.ref_pos > prev_seed.ref_pos
            }
        } else {
            false
        };

        if !same_chrom || !same_strand || !colinear {
            current += 1;
        }
        segment_id[si] = current;
    }
    current
}

/// Construct read explanations using seed-level weighted interval scheduling.
///
/// Flattens all seeds from all clusters, runs WIS on forward-read intervals
/// to find the best non-overlapping seed set (primary explanation), then
/// repeats on the remaining seeds (secondary explanation). Seeds are grouped
/// into segments by (chrom, strand, colinearity).
///
/// Results are written to the WIS debug SAM file: one record per seed with
/// tags XE (explanation: 0=primary, 1=secondary, -1=unselected),
/// XS (segment id), and XW (weight used in WIS).
#[cfg(feature = "conventional")]
fn construct_read_explanations(
    all_clusters: &[SeedCluster],
    seq_len: usize,
    reference: &InMemoryReference,
    read_name: &str,
    seq: &[u8],
    qual: &[u8],
    rc_seq: &[u8],
    rc_qual: &[u8],
) {
    // Flatten seeds from all clusters with provenance.
    let mut flat_seeds: Vec<FlatSeed> = Vec::new();
    for (ci, cluster) in all_clusters.iter().enumerate() {
        for (si, seed) in cluster.chain.iter().enumerate() {
            let (fwd_start, fwd_end) = seed.fwd_read_range(seq_len, cluster.is_reverse);
            let weight = seed.match_len as f64 / seed.kmer_uniqueness.max(1) as f64;
            flat_seeds.push(FlatSeed {
                fwd_start,
                fwd_end,
                weight,
                cluster_idx: ci,
                seed_idx: si,
            });
        }
    }

    // Sort by fwd_read_end (required for WIS DP).
    flat_seeds.sort_by_key(|s| (s.fwd_end, s.fwd_start));
    let n = flat_seeds.len();

    // Primary WIS: all seeds.
    let all_indices: Vec<usize> = (0..n).collect();
    let primary_selected = wis_select(&flat_seeds, &all_indices);

    // Group primary seeds into segments.
    let mut segment_id = vec![0usize; n];
    let mut explanation_ids = vec![-1i32; n];
    let last_primary_segment = group_into_segments(
        &primary_selected,
        &flat_seeds,
        all_clusters,
        &mut segment_id,
        0,
    );
    for &si in &primary_selected {
        explanation_ids[si] = 0;
    }

    // Secondary WIS: unselected seeds only.
    let primary_set: Vec<bool> = {
        let mut v = vec![false; n];
        for &si in &primary_selected {
            v[si] = true;
        }
        v
    };
    let unselected: Vec<usize> = (0..n).filter(|&i| !primary_set[i]).collect();
    let secondary_selected = wis_select(&flat_seeds, &unselected);

    // Group secondary into segments (numbering continues after primary).
    let mut sec_sorted = secondary_selected.clone();
    sec_sorted.sort_by_key(|&i| flat_seeds[i].fwd_start);
    group_into_segments(
        &sec_sorted,
        &flat_seeds,
        all_clusters,
        &mut segment_id,
        last_primary_segment + 1,
    );
    for &si in &secondary_selected {
        explanation_ids[si] = 1;
    }

    // Write debug SAM: one record per seed.
    for (i, fs) in flat_seeds.iter().enumerate() {
        let cluster = &all_clusters[fs.cluster_idx];
        let seed = &cluster.chain[fs.seed_idx];
        let chrom_name = reference.chrom_name(cluster.chrom_id);
        let strand_seq = if cluster.is_reverse { rc_seq } else { seq };
        let strand_qual = if cluster.is_reverse { rc_qual } else { qual };

        let flag: u16 = if cluster.is_reverse { 0x10 } else { 0 };
        let read_len = strand_seq.len();
        let hclip_start = seed.read_pos;
        let hclip_end = read_len.saturating_sub(seed.read_pos + seed.match_len);
        let cigar = match (hclip_start > 0, hclip_end > 0) {
            (true, true) => format!("{}H{}={}H", hclip_start, seed.match_len, hclip_end),
            (true, false) => format!("{}H{}=", hclip_start, seed.match_len),
            (false, true) => format!("{}={}H", seed.match_len, hclip_end),
            (false, false) => format!("{}=", seed.match_len),
        };
        let mapq = 60 / seed.kmer_uniqueness.max(1) as u8;
        let seq_slice = &strand_seq[seed.read_pos..seed.read_pos + seed.match_len];
        let seq_str = String::from_utf8_lossy(seq_slice);
        let qual_slice = &strand_qual[seed.read_pos..seed.read_pos + seed.match_len];
        let qual_str: String = qual_slice.iter().map(|&q| q as char).collect();

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tXE:i:{}\tXS:i:{}\tXW:f:{:.1}",
            read_name,
            flag,
            chrom_name,
            seed.ref_pos + 1,
            mapq,
            cigar,
            seq_str,
            qual_str,
            explanation_ids[i],
            segment_id[i] as i32,
            fs.weight,
        );
        WIS_SAM.append(&line);
    }

    log::debug!(
        "Read {}: WIS selected {} primary seeds in {} segments, {} total seeds",
        read_name,
        primary_selected.len(),
        last_primary_segment + 1,
        n,
    );
}

#[cfg(feature = "conventional")]
fn form_covering_sets(
    clusters: &[SeedCluster],
    read_name: &str,
    read_len: usize,
) -> Vec<Vec<SeedCluster>> {
    let mut order_by_quality: Vec<usize> = (0..clusters.len()).collect();
    let params = AlignParams::default();
    order_by_quality.sort_by_key(|i| OrderedFloat(-clusters[*i].quality(&params).value()));

    let mut segment_set_heap = Heap::new(SegmentSetHeap);
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
        if let Some((mut ranges, mut set, _)) = wanted_segment_set.take() {
            assert!(!ranges.overlaps(&(read_start, read_end)));
            ranges.add_range(read_start, read_end);
            set.push(i);
            let score = score_clusters(set.iter().map(|&j| &clusters[j]), read_len, &params);
            segment_set_heap.push((ranges, set, score));
        } else {
            let mut ranges = RangeSet::new();
            ranges.add_range(read_start, read_end);
            let set = vec![i];
            let score = score_clusters(set.iter().map(|&j| &clusters[j]), read_len, &params);
            segment_set_heap.push((ranges, set, score));
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
        .map(|(set_idx, (_, set, score))| {
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

/// Score a segment set using the lightweight estimated quality.
/// Same structure as `score_clusters` but uses `estimated_quality()`.
#[cfg(feature = "conventional")]
fn score_clusters_estimated<'a>(
    clusters: impl Iterator<Item = &'a SeedCluster>,
    read_len: usize,
    params: &AlignParams,
) -> f64 {
    let mut cluster_score: f64 = 0.0;
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    for cluster in clusters {
        cluster_score += cluster.estimated_quality(params).value();
        ranges.push(cluster.fwd_read_range(read_len));
    }
    ranges.sort_unstable();
    let mut gap_penalty: f64 = 0.0;
    if let Some(&(first_start, _)) = ranges.first() {
        if first_start > 0 {
            gap_penalty += params.quality(Op::new(Kind::Deletion, first_start)).value();
        }
    }
    for pair in ranges.windows(2) {
        let gap_len = pair[1].0.saturating_sub(pair[0].1);
        if gap_len > 0 {
            gap_penalty += params.quality(Op::new(Kind::Deletion, gap_len)).value();
        }
    }
    if let Some(&(_, last_end)) = ranges.last() {
        if last_end < read_len {
            gap_penalty += params
                .quality(Op::new(Kind::Deletion, read_len - last_end))
                .value();
        }
    }
    cluster_score + gap_penalty
}

/// Form covering sets using lightweight estimated quality (no gap alignments needed).
///
/// This is identical to `form_covering_sets` except it uses `estimated_quality()`
/// instead of `quality()` to rank clusters. Returns indices into the input slice
/// for each set, so the caller can identify which clusters were selected.
#[cfg(feature = "conventional")]
fn form_covering_sets_estimated(
    clusters: &[SeedCluster],
    _read_name: &str,
    read_len: usize,
) -> Vec<Vec<usize>> {
    let mut order_by_quality: Vec<usize> = (0..clusters.len()).collect();
    let params = AlignParams::default();
    order_by_quality
        .sort_by_key(|i| OrderedFloat(-clusters[*i].estimated_quality(&params).value()));

    let mut segment_set_heap = Heap::new(SegmentSetHeap);
    let mut wanted_segment_set: Option<SegmentSet> = None;
    let mut stack: Vec<SegmentSet> = vec![];

    for &i in order_by_quality.iter() {
        let cluster = &clusters[i];
        let (read_start, read_end) = cluster.fwd_read_range(read_len);

        while let Some(segment_set) = segment_set_heap.pop() {
            if segment_set.0.overlaps(&(read_start, read_end)) {
                stack.push(segment_set);
            } else {
                wanted_segment_set = Some(segment_set);
                break;
            }
        }

        if let Some((mut ranges, mut set, _)) = wanted_segment_set.take() {
            assert!(!ranges.overlaps(&(read_start, read_end)));
            ranges.add_range(read_start, read_end);
            set.push(i);
            let score =
                score_clusters_estimated(set.iter().map(|&j| &clusters[j]), read_len, &params);
            segment_set_heap.push((ranges, set, score));
        } else {
            let mut ranges = RangeSet::new();
            ranges.add_range(read_start, read_end);
            let set = vec![i];
            let score =
                score_clusters_estimated(set.iter().map(|&j| &clusters[j]), read_len, &params);
            segment_set_heap.push((ranges, set, score));
        }

        while let Some(segment_set) = stack.pop() {
            segment_set_heap.push(segment_set);
        }
    }

    segment_set_heap.drain().map(|(_, set, _)| set).collect()
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
                let mut aligner = Aligner::new();
                while let Ok(work) = receiver.recv() {
                    align_read(
                        index,
                        reference,
                        &writer,
                        &work.name,
                        &work.seq,
                        &work.qual,
                        &params,
                        &mut aligner,
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

/// Compute the "core" read-end of a cluster by trimming minor tail seeds.
///
/// Cluster boundaries (`read_end`) can be inflated by small noisy seeds beyond
/// the main alignment (e.g., false-positive matches in a deleted region). This
/// function finds the rightmost seed that's part of the cluster's main mass by
/// trimming up to 10% of total seed length from the tail.
#[cfg(feature = "conventional")]
fn core_read_end(cluster: &SeedCluster) -> usize {
    let total: usize = cluster.chain.iter().map(|s| s.match_len).sum();
    if total == 0 {
        return cluster.read_end;
    }
    let trim_budget = total / 10;
    let mut trimmed = 0;
    for seed in cluster.chain.iter().rev() {
        if trimmed + seed.match_len > trim_budget {
            return seed.read_end();
        }
        trimmed += seed.match_len;
    }
    cluster.read_end
}

/// Compute the "core" read-start of a cluster by trimming minor head seeds.
#[cfg(feature = "conventional")]
fn core_read_start(cluster: &SeedCluster) -> usize {
    let total: usize = cluster.chain.iter().map(|s| s.match_len).sum();
    if total == 0 {
        return cluster.read_start;
    }
    let trim_budget = total / 10;
    let mut trimmed = 0;
    for seed in &cluster.chain {
        if trimmed + seed.match_len > trim_budget {
            return seed.read_pos;
        }
        trimmed += seed.match_len;
    }
    cluster.read_start
}

/// Merge clusters on the same chromosome and strand that are colinear and
/// nearly abutting on the read, indicating a deletion that the chainer
/// could not bridge (due to `MAX_DIAGONAL_DIST`).
///
/// For a genuine deletion the read is continuous across the breakpoint, so
/// the two clusters should be almost exactly adjacent on the read. Small
/// overlaps (≤ `del_merge_max_read_overlap`) are permitted for microhomology
/// at breakpoints; small gaps (≤ `del_merge_max_read_gap`) for seed placement
/// granularity.
///
/// Cluster boundaries (`read_start`/`read_end`) can be inflated by noisy
/// tail seeds that mapped into the deleted region, so the adjacency check
/// uses **core boundaries** — the extent of seeds near the cluster's
/// length-weighted median diagonal — rather than outermost seed positions.
///
/// After merging seed lists, `SeedCluster::new()` re-resolves overlaps and
/// filters noisy seeds, so the merged cluster is clean. The resulting gap
/// between the two seed groups will be aligned by `align_gaps()` and typically
/// produces a pure deletion (`D`) CIGAR operation.
#[cfg(feature = "conventional")]
fn merge_deletion_clusters(
    clusters: &mut Vec<SeedCluster>,
    read_len: usize,
    min_seed_length: usize,
) {
    let cfg = config::get();
    let max_read_overlap = cfg.seeding.del_merge_max_read_overlap as i64;
    let max_read_gap = cfg.seeding.del_merge_max_read_gap as i64;
    let max_ref_gap = cfg.seeding.del_merge_max_ref_gap as i64;

    if clusters.len() < 2 {
        return;
    }

    let mut merged_any = true;
    while merged_any {
        merged_any = false;

        // Sort by (chrom_id, is_reverse, read_start) in strand coordinates
        clusters.sort_by_key(|c| (c.chrom_id, c.is_reverse, c.read_start));

        let mut i = 0;
        while i + 1 < clusters.len() {
            let can_merge = {
                let a = &clusters[i];
                let b = &clusters[i + 1];

                if a.chrom_id != b.chrom_id
                    || a.is_reverse != b.is_reverse
                    || a.ref_end() > b.ref_start()
                    || (b.ref_start() as i64 - a.ref_end() as i64) > max_ref_gap
                {
                    false
                } else {
                    // Use core boundaries: the read extent of seeds close to
                    // each cluster's median diagonal. This ignores noisy tail
                    // seeds that mapped into the deleted region.
                    let a_core_end = core_read_end(a);
                    let b_core_start = core_read_start(b);
                    let read_delta = b_core_start as i64 - a_core_end as i64;
                    read_delta >= -max_read_overlap && read_delta <= max_read_gap
                }
            };

            if can_merge {
                let is_reverse = clusters[i].is_reverse;
                let mut combined = clusters[i].chain.clone();
                combined.extend(clusters[i + 1].chain.iter().cloned());

                if let Some(new_cluster) = SeedCluster::new(combined, is_reverse, min_seed_length) {
                    log::info!(
                        "Merged colinear clusters: read [{}-{}]+[{}-{}] -> [{}-{}], \
                         ref [{}-{}]+[{}-{}], gap {}bp",
                        clusters[i].read_start,
                        clusters[i].read_end,
                        clusters[i + 1].read_start,
                        clusters[i + 1].read_end,
                        new_cluster.read_start,
                        new_cluster.read_end,
                        clusters[i].ref_start(),
                        clusters[i].ref_end(),
                        clusters[i + 1].ref_start(),
                        clusters[i + 1].ref_end(),
                        clusters[i + 1].ref_start() as i64 - clusters[i].ref_end() as i64,
                    );
                    clusters[i] = new_cluster;
                    clusters.remove(i + 1);
                    merged_any = true;
                    continue; // retry at same position
                }
            }
            i += 1;
        }
    }

    // Re-sort by fwd_read_range for the rest of the pipeline
    clusters.sort_by_key(|c| c.fwd_read_range(read_len));
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
    // SeedCluster tests (conventional feature only)
    // =========================================================================

    #[cfg(feature = "conventional")]
    mod seed_cluster_tests {
        use super::*;

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

        #[test]
        fn test_filter_misplaced_seeds_removes_bad_region() {
            // Construct a chain where seeds 2 & 3 are misplaced:
            // They introduce a big insertion gap followed by a big deletion gap
            // (or vice versa), indicating contradictory anchor placement.
            //
            // Good seeds are on diagonal ~0 (ref_pos ≈ read_pos).
            // Misplaced seeds shift the diagonal creating simultaneous ins + del.
            //
            // Layout (read_pos, ref_pos) for 5 seeds of length 20:
            //   seed 0: read 0,   ref 0     → gap to seed 1: read_delta 30, ref_delta 30, gap = 0   (ok)
            //   seed 1: read 50,  ref 50    → gap to seed 2: read_delta 30, ref_delta 80, gap = -50 (del)
            //   seed 2: read 100, ref 150   → gap to seed 3: read_delta 80, ref_delta 30, gap = +50 (ins)
            //   seed 3: read 200, ref 200   → gap to seed 4: read_delta 30, ref_delta 30, gap = 0   (ok)
            //   seed 4: read 250, ref 250
            //
            // At the two long gaps: n_del = 50, n_ins = 50 → diff = 2*min(50,50) = 100 > 40.
            // Seeds 2 and 3 (indices 2..4 from long_gap at index 2 to long_gap at index 3)
            // should be removed, but the minimap2 algorithm marks seeds from K[st]..K[en],
            // which means seed at index 2 is removed but seed at index 3 is kept.
            let seeds = vec![
                make_hit(0, 0, 0, 20),
                make_hit(0, 50, 50, 20),
                make_hit(0, 150, 100, 20), // misplaced: shifted +100bp on ref
                make_hit(0, 200, 200, 20),
                make_hit(0, 250, 250, 20),
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();

            // The misplaced seed (ref 150 @ read 100) should have been removed,
            // leaving 4 seeds. The remaining seeds should all be on the good diagonal.
            assert_eq!(cluster.chain.len(), 4);
            assert_eq!(cluster.chain[0].ref_pos, 0);
            assert_eq!(cluster.chain[1].ref_pos, 50);
            assert_eq!(cluster.chain[2].ref_pos, 200);
            assert_eq!(cluster.chain[3].ref_pos, 250);
        }

        #[test]
        fn test_filter_misplaced_seeds_good_chain_unchanged() {
            // A well-behaved chain (all seeds on the same diagonal) should not
            // have any seeds removed.
            let seeds = vec![
                make_hit(0, 100, 0, 20),
                make_hit(0, 150, 50, 20),
                make_hit(0, 200, 100, 20),
                make_hit(0, 250, 150, 20),
                make_hit(0, 300, 200, 20),
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();
            assert_eq!(cluster.chain.len(), 5);
        }

        // ── Jittery seed filter tests ────────────────────────────────────────────

        #[test]
        fn test_filter_jittery_removes_bouncing_seeds() {
            // Simulate a chain with a stable region followed by jittery seeds
            // from different repeat copies. Diagonal = ref_pos - read_pos.
            //
            // Stable region (diagonal ~1000, spread ~2):
            //   seed 0: read 0,    ref 1000, len 100 → diag 1000
            //   seed 1: read 120,  ref 1121, len 80  → diag 1001
            //   seed 2: read 220,  ref 1220, len 60  → diag 1000
            //
            // Jittery region (diagonal bouncing by 40-80bp, short seeds):
            //   seed 3: read 310,  ref 1350, len 25  → diag 1040 (shift +40)
            //   seed 4: read 360,  ref 1340, len 25  → diag  980 (shift -60)
            //   seed 5: read 410,  ref 1460, len 25  → diag 1050 (shift +70)
            //   seed 6: read 460,  ref 1430, len 25  → diag  970 (shift -80)
            //
            // Stable again:
            //   seed 7: read 520,  ref 1520, len 100 → diag 1000
            //   seed 8: read 640,  ref 1641, len 80  → diag 1001
            //
            // In the jittery zone (seeds 3-6), the shifts are 40+60+70+80 = 250
            // over a ref span of ~1430+25 - 1350 = 105 → density = 250/105 ≈ 2.4
            // which far exceeds the default threshold of 0.15.
            let seeds = vec![
                make_hit(0, 1000, 0, 100),
                make_hit(0, 1121, 120, 80),
                make_hit(0, 1220, 220, 60),
                make_hit(0, 1350, 310, 25),
                make_hit(0, 1340, 360, 25),
                make_hit(0, 1460, 410, 25),
                make_hit(0, 1430, 460, 25),
                make_hit(0, 1520, 520, 100),
                make_hit(0, 1641, 640, 80),
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            // The jittery interior seeds (3-6) should be removed.
            // Boundary seeds and stable seeds should remain.
            let remaining_reads: Vec<usize> = chain.iter().map(|s| s.read_pos).collect();
            // Seeds 0,1,2 (stable), 7,8 (stable) should survive.
            // Some boundary seeds from the jittery zone might also survive
            // since we keep the first and last of each window.
            assert!(
                remaining_reads.len() >= 5,
                "should keep at least the 5 stable seeds, got {:?}",
                remaining_reads
            );
            // The core jittery seeds (interior of the window) should be gone.
            // Seeds at read_pos 360 and 410 are always interior to any 4-gap window
            // that covers the jittery region.
            assert!(
                !remaining_reads.contains(&360) || !remaining_reads.contains(&410),
                "at least some jittery seeds should be removed, got {:?}",
                remaining_reads
            );
        }

        #[test]
        fn test_filter_jittery_preserves_stable_chain() {
            // A chain with consistent diagonal (small wobble from real indels)
            // should not lose any seeds.
            //
            // Diagonal ~500, wobble ≤ 5bp.
            let seeds = vec![
                make_hit(0, 500, 0, 50),   // diag 500
                make_hit(0, 553, 50, 40),  // diag 503
                make_hit(0, 600, 98, 30),  // diag 502
                make_hit(0, 635, 130, 35), // diag 505
                make_hit(0, 670, 165, 45), // diag 505
                make_hit(0, 718, 215, 50), // diag 503
                make_hit(0, 770, 268, 40), // diag 502
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            assert_eq!(chain.len(), 7, "stable chain should keep all seeds");
        }

        #[test]
        fn test_filter_jittery_single_shift_preserved() {
            // A single large diagonal shift (real indel) shouldn't trigger
            // removal because one shift doesn't make a high density over the
            // full window.
            //
            // Seeds: stable diagonal 500, then shift to 550 (real 50bp deletion),
            // then stable at 550.
            let seeds = vec![
                make_hit(0, 500, 0, 80),    // diag 500
                make_hit(0, 590, 88, 60),   // diag 502
                make_hit(0, 660, 150, 50),  // diag 510 (slight shift)
                make_hit(0, 770, 210, 40),  // diag 560 (big shift — real indel)
                make_hit(0, 870, 310, 80),  // diag 560
                make_hit(0, 1010, 448, 60), // diag 562
                make_hit(0, 1130, 568, 50), // diag 562
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            assert_eq!(
                chain.len(),
                7,
                "single diagonal shift (real indel) should not trigger removal"
            );
        }

        #[test]
        fn test_filter_jittery_disabled_at_zero() {
            // When threshold is 0.0, filter should be disabled.
            let seeds = vec![
                make_hit(0, 1000, 0, 25),
                make_hit(0, 1100, 30, 25),  // big shift
                make_hit(0, 1000, 60, 25),  // shift back
                make_hit(0, 1100, 90, 25),  // shift again
                make_hit(0, 1000, 120, 25), // shift back
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.0, 4);

            assert_eq!(chain.len(), 5, "filter disabled at threshold 0.0");
        }

        #[test]
        fn test_filter_jittery_too_few_seeds() {
            // Chain shorter than window_size+1 should be untouched.
            let seeds = vec![
                make_hit(0, 1000, 0, 25),
                make_hit(0, 1100, 30, 25),
                make_hit(0, 1000, 60, 25),
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            assert_eq!(chain.len(), 3, "short chain should be untouched");
        }
    } // mod seed_cluster_tests
}
