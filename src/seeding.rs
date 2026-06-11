use crate::reads::seeds::SeedHit;
use parallax::{
    config::SeedingConfig,
    index::{Index, PackedLocus, Strand},
    kmers::Kmer,
    reference::InMemoryReference,
    utils::hasher::FnvHasher,
};
use std::collections::HashMap;

pub struct SeedCollector {
    /// Seed hits collected from k-mer index lookups
    pub hits: Vec<SeedHit>,
    /// Scratch space for merging/deduplication
    merge_scratch: Vec<SeedHit>,
    /// Batch buffer for prefetched lookups: (read_pos, kmer_value)
    kmer_batch: Vec<(usize, u64)>,
    /// Deferred mid-frequency seeds: (read_pos, kmer_value, hit_count).
    /// Collected during Phase 1 and selectively rescued into gaps after
    /// merge+extend.
    deferred_seeds: Vec<(usize, u64, u32)>,
}

impl SeedCollector {
    /// Create a new collector with empty buffers
    pub fn new() -> Self {
        SeedCollector {
            hits: Vec::new(),
            merge_scratch: Vec::new(),
            kmer_batch: Vec::new(),
            deferred_seeds: Vec::new(),
        }
    }

    /// Sort, merge adjacent seeds on the same diagonal, extend exact matches,
    /// and remove duplicates.
    ///
    /// This is the core seed-consolidation pipeline (Phases 2–3c) used after
    /// initial seed collection and again after rescue. It operates in-place on
    /// `self.hits`, using `self.merge_scratch` as temporary storage.
    pub fn sort_merge_extend<const K: usize>(
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
                        hit.read_frequency,
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
    pub fn rescue_seeds<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        is_reverse: bool,
        index: &impl Index<K, S>,
        reference: &InMemoryReference,
        rescue_spacing: usize,
    ) -> usize {
        if self.deferred_seeds.is_empty() || rescue_spacing == 0 {
            return 0;
        }
        let seq_len = strand_seq.len();

        let read_strand = Strand::from_is_reverse(is_reverse);

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
                            for loc in loci {
                                let (chrom_id, chrom_pos, hit_strand) = loc.unpack();
                                let strand = read_strand.combine(&hit_strand);
                                self.hits.push(SeedHit::new(
                                    chrom_id, chrom_pos, read_pos, kmer_val, hit_count, K, strand,
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

    /// Instead of interleaving k-mer generation with index lookups (which causes
    /// serial cache misses in the multi-GB hash tables), this method:
    /// 1. Generates all syncmer k-mers into a batch buffer (sequential writes, cache-friendly)
    /// 2. Looks up all k-mers using `Index::lookup_batch` which prefetches PIPE steps ahead
    ///
    /// The rest of the pipeline (sort, merge, extend, dedup) is identical.
    pub fn gather_seeds_batched<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        is_reverse: bool,
        index: &impl Index<K, S>,
        reference: &InMemoryReference,
        read_name: &str,
        cfg: &SeedingConfig,
    ) {
        self.hits.clear();
        self.kmer_batch.clear();
        self.deferred_seeds.clear();

        let read_strand = Strand::from_is_reverse(is_reverse);

        // Phase 1a: Generate all syncmers into the batch buffer
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            strand_seq,
            [(); S],
            |pos, kmer| {
                self.kmer_batch.push((pos, kmer.0));
            },
        );

        // Phase 1a': Count how many times each kmer value appears in the read.
        let mut read_freq: HashMap<u64, u32> = HashMap::new();
        for &(_, kmer_val) in &self.kmer_batch {
            *read_freq.entry(kmer_val).or_insert(0) += 1;
        }

        // Phase 1b: Batched lookup with prefetching
        let max_occ = cfg.max_seed_occurrences;
        let mid_occ = cfg.mid_seed_occurrences;
        index.lookup_batch(&self.kmer_batch, |read_pos, kmer_val, hit_count, loci| {
            let rf = *read_freq.get(&kmer_val).unwrap_or(&1);
            if mid_occ > 0 && hit_count > mid_occ && hit_count <= max_occ {
                // Mid-frequency: defer for potential rescue
                self.deferred_seeds.push((read_pos, kmer_val, hit_count as u32));
            } else if hit_count <= max_occ {
                // Low-frequency (or rescue disabled): collect immediately
                for loc in loci {
                    let (chrom_id, chrom_pos, hit_strand) = loc.unpack();
                    let strand = read_strand.combine(&hit_strand);
                    self.hits.push(SeedHit::with_read_frequency(
                        chrom_id, chrom_pos, read_pos, kmer_val, hit_count as u32, rf, K, strand,
                    ));
                }
            }
            // hit_count > max_occ: skip entirely
        });

        // Phases 2–3c: Sort, merge, extend, dedup
        self.sort_merge_extend::<K>(strand_seq, reference);

        // Phase 3d: Rescue deferred mid-frequency seeds into coverage gaps
        let rescued = self.rescue_seeds::<K, S>(strand_seq, is_reverse, index, reference, cfg.rescue_spacing);
        if rescued > 0 {
            let strand_name = if is_reverse { "REV" } else { "FWD" };
            log::debug!("{read_name} {strand_name}: rescued {rescued} deferred seeds into gaps");
        }
    }
}
