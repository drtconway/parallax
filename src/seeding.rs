use crate::reads::seeds::SeedHit;
use parallax::{config::SeedingConfig, index::IndexHit, reference::InMemoryReference};
use std::collections::HashMap;

pub struct SeedCollector {
    /// Seed hits collected from k-mer index lookups
    pub hits: Vec<SeedHit>,
    /// Scratch space for merging/deduplication
    merge_scratch: Vec<SeedHit>,
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
            deferred_seeds: Vec::new(),
        }
    }

    /// Sort, merge adjacent seeds on the same diagonal, extend exact matches,
    /// and remove duplicates.
    ///
    /// This is the core seed-consolidation pipeline (Phases 2–3c) used after
    /// initial seed collection and again after rescue. It operates in-place on
    /// `self.hits`, using `self.merge_scratch` as temporary storage.
    pub fn sort_merge_extend(&mut self, strand_seq: &[u8], reference: &InMemoryReference) {
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
                        hit.match_len,
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
    pub fn rescue_seeds(
        &mut self,
        strand_seq: &[u8],
        index: &dyn parallax::index::Index,
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
                        if let Some(hit) = index.lookup_kmer(kmer_val) {
                            let IndexHit { loci, k, .. } = hit;
                            for &locus in loci {
                                let (chrom_id, chrom_pos) = index.unpack_locus(locus);

                                self.hits.push(SeedHit::new(
                                    chrom_id, chrom_pos, read_pos, kmer_val, hit_count, k,
                                ));
                            }
                        }

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
            self.sort_merge_extend(strand_seq, reference);
        }

        rescued
    }

    pub fn gather_seeds_batched(
        &mut self,
        strand_seq: &[u8],
        is_reverse: bool,
        index: &dyn parallax::index::Index,
        reference: &InMemoryReference,
        read_name: &str,
        cfg: &SeedingConfig,
    ) {
        self.hits.clear();
        self.deferred_seeds.clear();

        let max_occ = cfg.max_seed_occurrences;
        let mid_occ = cfg.mid_seed_occurrences;

        // Phase 1: kmerize and look up all hits via find_seeds.
        // read_frequency is filled in as a second pass below, so use new() with rf=1.
        index.find_seeds(strand_seq, &mut |hit| {
            let IndexHit {
                query_pos,
                seed_kmer,
                loci,
                k,
            } = hit;
            let hit_count = loci.len();
            if mid_occ > 0 && hit_count > mid_occ && hit_count <= max_occ {
                self.deferred_seeds
                    .push((query_pos, seed_kmer, hit_count as u32));
            } else if hit_count <= max_occ {
                for &locus in loci {
                    let (chrom_id, chrom_pos) = index.unpack_locus(locus);
                    self.hits.push(SeedHit::new(
                        chrom_id,
                        chrom_pos,
                        query_pos,
                        seed_kmer,
                        hit_count as u32,
                        k,
                    ));
                }
            }
        });

        // Phase 1b: Count how many times each kmer appears among the collected hits
        // and patch read_frequency on each hit.
        let mut read_freq: HashMap<u64, u32> = HashMap::new();
        for hit in &self.hits {
            *read_freq.entry(hit.kmer).or_insert(0) += 1;
        }
        for hit in &mut self.hits {
            hit.read_frequency = *read_freq.get(&hit.kmer).unwrap_or(&1);
        }

        // Phases 2–3c: Sort, merge, extend, dedup
        self.sort_merge_extend(strand_seq, reference);

        // Phase 3d: Rescue deferred mid-frequency seeds into coverage gaps
        let rescued =
            self.rescue_seeds(strand_seq, index, reference, cfg.rescue_spacing);
        if rescued > 0 {
            let strand_name = if is_reverse { "REV" } else { "FWD" };
            log::debug!("{read_name} {strand_name}: rescued {rescued} deferred seeds into gaps");
        }
    }
}
