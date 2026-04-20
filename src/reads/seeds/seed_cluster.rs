
use super::{CHAINS_SAM, ClusterError, FLAG_REVERSE, FLAG_SUPPLEMENTARY, GAP_ALIGNMENTS};
use crate::{align::{AlignParams, Aligner, Alignment, Kind, Op}, config, utils::{GroupsTrait, InterleaveTrait}};
use super::SeedHit;
use crate::scores::QualityScore;

/// A cluster of seed hits with its LIS chain, before alignment building.
/// This intermediate structure allows for cross-strand gap analysis before
/// committing to alignment construction.
#[derive(Clone, Debug)]
pub struct SeedCluster {
    /// Read region covered by this cluster
    pub read_start: usize,
    pub read_end: usize,

    /// The chain of colinear seeds (sorted by read position)
    pub chain: Vec<SeedHit>,

    /// Which strand this cluster came from
    pub is_reverse: bool,

    /// Chromosome this cluster aligns to
    pub chrom_id: usize,

    /// Alignments across gaps between seeds.
    /// Initially empty. After calling `align_gaps()`, contains one entry per gap
    /// (i.e., `chain.len() - 1` entries). Each entry is `Some(alignment)` if the
    /// block aligner alignment succeeded, or `None` if it failed or was skipped.
    ///
    /// This allows gap-splitting decisions to consider actual alignment quality
    /// rather than just seed absence.
    pub gap_alignments: Vec<Alignment>,

    /// Split-fill tags describing the filler(s) that triggered a split.
    /// Each entry encodes the filler's read/ref locus and its position relative
    /// to this segment: `[*]read_start-read_end[*];chrom:ref_start-ref_end[strand]`
    /// where `*` appears before or after the read locus indicating whether the
    /// filler precedes or follows this segment in read space.
    pub split_fill_tags: Vec<String>,
}

impl SeedCluster {
    /// Create a new cluster from a chain of seeds.
    ///
    /// The chain is sorted by read position and overlapping seeds are resolved
    /// by truncating the right-hand seed. Seeds that become too short (less than
    /// `min_seed_length`) are dropped entirely.
    ///
    /// # Arguments
    /// * `chain` - The seeds to form into a cluster
    /// * `is_reverse` - Whether this cluster is on the reverse strand
    /// * `min_seed_length` - Minimum seed length after truncation (typically K/2)
    pub fn new(mut chain: Vec<SeedHit>, is_reverse: bool, min_seed_length: usize) -> Option<Self> {
        if chain.is_empty() {
            return None;
        }

        // Sort by strand read position for alignment building
        // (For forward strand this == forward order; for reverse strand this is RC order)
        chain.sort_by_key(|hit| hit.read_pos);

        // Resolve overlaps by truncating right-hand seeds
        Self::resolve_overlaps(&mut chain, min_seed_length);

        // Remove misplaced anchors that would require simultaneous insertion
        // and deletion during gap alignment (minimap2's mm_filter_bad_seeds heuristic)
        let threshold = config::get().seeding.misplaced_seed_threshold;
        Self::filter_misplaced_seeds(&mut chain, threshold);

        // Remove seeds in jittery regions where the diagonal bounces around,
        // indicating seeds from different copies of a tandem repeat.
        let jitter_cfg = &config::get().seeding;
        Self::filter_jittery_seeds(
            &mut chain,
            jitter_cfg.jitter_density_threshold,
            jitter_cfg.jitter_window,
        );

        if chain.is_empty() {
            return None;
        }

        let read_start = chain.first().map(|h| h.read_pos).unwrap_or(0);
        let read_end = chain.last().map(|h| h.read_end()).unwrap_or(0);
        let chrom_id = chain[0].chrom_id;

        Some(SeedCluster {
            chrom_id,
            read_start,
            read_end,
            chain,
            is_reverse,
            gap_alignments: Vec::new(),
            split_fill_tags: Vec::new(),
        })
    }

    /// Convert the seed cluster to a full alignment with extensions.
    ///
    /// Extension budgets are in *forward-strand read coordinates*:
    /// - `left_budget`: max read bases to extend toward the start of the forward read
    /// - `right_budget`: max read bases to extend toward the end of the forward read
    ///
    /// The method translates these into strand-space operations (a forward-strand
    /// left extension is a strand-space left extension for forward clusters, but a
    /// strand-space right extension for reverse-complement clusters).
    ///
    /// Returns `(Alignment, ref_start_adjustment, seq_start, seq_end)` where:
    /// - `ref_start_adjustment` is the number of reference bases consumed by the left extension.
    ///   The caller should compute the actual reference start as `ref_start() - ref_start_adjustment`.
    /// - `seq_start` and `seq_end` define the range of `strand_seq` that should be in the SEQ field.
    ///   Hard clips are excluded from SEQ; soft clips and aligned bases are included.
    pub fn into_alignment(
        self,
        left_budget: usize,
        right_budget: usize,
        soft_clip: bool,
        strand_seq: &[u8],
        ref_seq: &[u8],
        aligner: &mut Aligner,
    ) -> (Alignment, usize, usize, usize) {
        let mut left_clip: Option<Op> = None;
        let mut left_extension: Option<Vec<Op>> = None;
        let mut right_extension: Option<Vec<Op>> = None;
        let mut right_clip: Option<Op> = None;
        let mut ref_start_adjustment: usize = 0;

        // Track sequence range for SEQ field
        // Start with cluster range, expand if we have real extensions or soft clips
        let mut seq_start = self.read_start;
        let mut seq_end = self.read_end;

        // Translate forward-strand budgets into strand-space budgets.
        // For a reverse-complement cluster, forward-left maps to strand-right and vice versa.
        let (strand_left_budget, strand_right_budget) = if self.is_reverse {
            (right_budget, left_budget)
        } else {
            (left_budget, right_budget)
        };

        // Extend left in strand space (toward lower strand read positions).
        if strand_left_budget > 0 && self.read_start > 0 {
            let available = self.read_start.min(strand_left_budget);
            let read_ext_seq = &strand_seq[self.read_start - available..self.read_start];
            let ref_ext_seq = &ref_seq[..self.ref_start()];
            let ext = aligner
                .extend_left(read_ext_seq, ref_ext_seq)
                .expect("alignment extension failed");
            ref_start_adjustment = ext.reference_consumed();
            let ext_query_len = ext.query_length();
            seq_start = self.read_start - ext_query_len;
            left_extension = Some(ext.cigar);
            // Clip for everything before where the extension reached
            if seq_start > 0 {
                if soft_clip {
                    left_clip = Some(Op::new(Kind::SoftClip, seq_start));
                    seq_start = 0;
                } else {
                    left_clip = Some(Op::new(Kind::HardClip, seq_start));
                }
            }
        }

        // Extend right in strand space (toward higher strand read positions).
        if strand_right_budget > 0 && self.read_end < strand_seq.len() {
            let available = (strand_seq.len() - self.read_end).min(strand_right_budget);
            let read_ext_seq = &strand_seq[self.read_end..self.read_end + available];
            let ref_end = self.chain.last().unwrap().ref_end();
            let ref_ext_seq = &ref_seq[ref_end..];
            let ext = aligner
                .extend_right(read_ext_seq, ref_ext_seq)
                .expect("alignment extension failed");
            let ext_query_len = ext.query_length();
            seq_end = self.read_end + ext_query_len;
            right_extension = Some(ext.cigar);
            // Clip for everything after where the extension reached
            let remaining = strand_seq.len() - seq_end;
            if remaining > 0 {
                if soft_clip {
                    seq_end = strand_seq.len();
                    right_clip = Some(Op::new(Kind::SoftClip, remaining));
                } else {
                    right_clip = Some(Op::new(Kind::HardClip, remaining));
                }
            }
        }

        // Clip for bases before our leftward reach in strand space
        if left_extension.is_none() && left_clip.is_none() && self.read_start > 0 {
            // No left extension was attempted (zero budget on this side)
            if soft_clip {
                left_clip = Some(Op::new(Kind::SoftClip, self.read_start));
                seq_start = 0;
            } else {
                left_clip = Some(Op::new(Kind::HardClip, self.read_start));
            }
        }

        // Clip for bases beyond our rightward reach in strand space
        if right_extension.is_none() && right_clip.is_none() && self.read_end < strand_seq.len() {
            let clip_len = strand_seq.len() - self.read_end;
            if soft_clip {
                seq_end = strand_seq.len();
                right_clip = Some(Op::new(Kind::SoftClip, clip_len));
            } else {
                right_clip = Some(Op::new(Kind::HardClip, clip_len));
            }
        }

        let self_ref_start = self.ref_start();
        let seed_parts = self
            .chain
            .into_iter()
            .map(|hit| vec![Op::new(Kind::SequenceMatch, hit.match_len)]);
        let gap_alignments = self.gap_alignments.into_iter().map(|aln| aln.cigar);
        let interleaved = seed_parts.interleave(gap_alignments);
    
        // Build CIGAR: left_clip + left_extension + seeds/gaps + right_extension + right_clip
        let cigar_ops: Vec<Op> = left_clip
            .into_iter()
            .chain(left_extension.into_iter().flatten())
            .chain(interleaved.flatten())
            .chain(right_extension.into_iter().flatten())
            .chain(right_clip.into_iter())
            .collect();
        let mut alignment = Alignment::from(cigar_ops);
        alignment.normalize();

        // Left-align indels for deterministic placement in tandem repeats.
        // The aligned query is strand_seq[seq_start..seq_end] and the aligned
        // reference starts at ref_start - ref_start_adjustment.
        let ref_start_abs = self_ref_start.saturating_sub(ref_start_adjustment);
        let ref_consumed = alignment.reference_consumed();
        let ref_end_abs = (ref_start_abs + ref_consumed).min(ref_seq.len());
        let query_slice = &strand_seq[seq_start..seq_end];
        let ref_slice = &ref_seq[ref_start_abs..ref_end_abs];
        aligner
            .indel_shifter
            .left_align_indels(&mut alignment, query_slice, ref_slice);

        (alignment, ref_start_adjustment, seq_start, seq_end)
    }

    /// Return a summary for populating the supplementary alignment (SA) tag.
    /// (chrom_id, read_start, read_end, is_reverse, condensed CIGAR string, number of mismatches)
    pub fn summary(&self, read_len: usize) -> (usize, usize, bool, String, usize) {
        let cigar = self.cigar_summary(read_len);
        let mismatch_count = self.mismatch_count();
        (
            self.chrom_id,
            self.ref_start(),
            self.is_reverse,
            cigar,
            mismatch_count,
        )
    }

    /// Reference start position of the seed chain.
    pub fn ref_start(&self) -> usize {
        self.chain.first().map(|h| h.ref_pos).unwrap_or(0)
    }

    /// Reference end position of the seed chain.
    pub fn ref_end(&self) -> usize {
        self.chain.last().map(|h| h.ref_end()).unwrap_or(0)
    }

    pub fn fwd_read_range(&self, read_len: usize) -> (usize, usize) {
        if self.is_reverse {
            (read_len - self.read_end, read_len - self.read_start)
        } else {
            (self.read_start, self.read_end)
        }
    }

    pub fn read_range(&self) -> std::ops::Range<usize> {
        self.read_start..self.read_end
    }

    /// Generate a condensed CIGAR string summarizing the accumulated insertions, deletions, matches and mismatches.
    /// This is used for the SA tag and other summaries where we want a quick representation of the alignment without the full CIGAR.
    pub fn cigar_summary(&self, read_len: usize) -> String {
        let mut matches = 0;
        let mut mismatches = 0;
        let mut indels: i64 = 0;
        for seed in &self.chain {
            matches += seed.match_len;
        }
        for gap in &self.gap_alignments {
            for op in &gap.cigar {
                match op.kind() {
                    Kind::SequenceMatch => matches += op.len(),
                    Kind::SequenceMismatch => mismatches += op.len(),
                    Kind::Insertion => indels += op.len() as i64,
                    Kind::Deletion => indels -= op.len() as i64,
                    _ => {}
                }
            }
        }
        let left_clip = if self.read_start > 0 {
            format!("{}S", self.read_start)
        } else {
            String::new()
        };
        let right_clip = if self.read_end < read_len {
            format!("{}S", read_len - self.read_end)
        } else {
            String::new()
        };
        if indels == 0 {
            format!("{}{}={}X{}", left_clip, matches, mismatches, right_clip)
        } else if indels > 0 {
            format!("{}{}={}X{}I{}", left_clip, matches, mismatches, indels, right_clip)
        } else {
            format!("{}{}={}X{}D{}", left_clip, matches, mismatches, -indels, right_clip)
        }
    }

    /// Calculate the total number of mismatches across all gap alignments in this cluster.
    /// Seed matches are assumed to be perfect and are not counted as mismatches.
    pub fn mismatch_count(&self) -> usize {
        self.gap_alignments
            .iter()
            .map(|align| align.mismatch_count())
            .sum()
    }

    /// Resolve overlapping seeds by truncating the right-hand seed.
    ///
    /// When consecutive seeds overlap (in read coordinates), the second seed
    /// is truncated so it starts where the first seed ends. Seeds that become
    /// shorter than `min_seed_length` are dropped entirely.
    ///
    /// This ensures downstream code never has to deal with overlapping seeds,
    /// simplifying gap calculation and CIGAR generation.
    pub(crate) fn resolve_overlaps(chain: &mut Vec<SeedHit>, min_seed_length: usize) {
        if chain.len() < 2 {
            return;
        }

        let mut write_idx = 1; // First seed always kept
        let mut prev_read_end = chain[0].read_end();
        let mut prev_ref_end = chain[0].ref_end();

        for read_idx in 1..chain.len() {
            // Copy the hit since SeedHit is Copy
            let hit = chain[read_idx];

            // Calculate overlap in both read and reference space
            let read_overlap = prev_read_end.saturating_sub(hit.read_pos);
            let ref_overlap = prev_ref_end.saturating_sub(hit.ref_pos);

            // Use the maximum overlap to ensure consistency
            let overlap = read_overlap.max(ref_overlap);

            if overlap == 0 {
                // No overlap - keep seed as-is
                chain[write_idx] = hit;
                prev_read_end = hit.read_end();
                prev_ref_end = hit.ref_end();
                write_idx += 1;
            } else if overlap < hit.match_len {
                // Partial overlap - truncate the seed
                let new_match_len = hit.match_len - overlap;
                if new_match_len >= min_seed_length {
                    // Keep the truncated seed
                    let mut truncated = hit;
                    truncated.read_pos += overlap;
                    truncated.ref_pos += overlap;
                    truncated.match_len = new_match_len;
                    // Diagonal stays the same (ref_pos - read_pos is unchanged)
                    chain[write_idx] = truncated;
                    prev_read_end = truncated.read_end();
                    prev_ref_end = truncated.ref_end();
                    write_idx += 1;
                }
                // else: truncated seed too short, drop it
            }
            // else: fully overlapped, drop the seed
        }

        chain.truncate(write_idx);
    }

    /// Filter misplaced seeds that would require simultaneous insertion AND
    /// deletion during gap alignment.
    ///
    /// Adapted from minimap2's `mm_filter_bad_seeds()` (align.c). The
    /// algorithm walks consecutive seeds in the chain and identifies
    /// "long‐gap" transitions where |(read_delta − ref_delta)| > 10 bp.
    /// Across sliding windows of these long‐gap positions it accumulates the
    /// total insertion and deletion implied by each transition and computes
    /// `diff = 2 * min(total_ins, total_del)`.  When `diff` exceeds the
    /// configured `threshold` the seeds spanning that window are removed,
    /// leaving the flanking good seeds to be bridged by the block aligner in
    /// `align_gaps()`.
    ///
    /// A threshold of 0 disables the filter.
    pub(crate) fn filter_misplaced_seeds(chain: &mut Vec<SeedHit>, threshold: i64) {
        /// Minimum |read_delta − ref_delta| to count as a "long gap".
        pub(crate) const MIN_GAP: i64 = 10;
        /// Maximum number of long‐gap positions in one sliding window.
        pub(crate) const MAX_WINDOW: usize = 10;

        if chain.len() < 3 || threshold <= 0 {
            return;
        }

        // Step 1 — collect indices where consecutive seeds imply a long gap.
        // `long_gaps[k]` is the index *i* of the second seed in a pair whose
        // `read_delta − ref_delta` exceeds MIN_GAP.
        let mut long_gaps: Vec<usize> = Vec::new();
        for i in 1..chain.len() {
            let read_delta = chain[i].read_pos as i64 - chain[i - 1].read_end() as i64;
            let ref_delta = chain[i].ref_pos as i64 - chain[i - 1].ref_end() as i64;
            if (read_delta - ref_delta).abs() > MIN_GAP {
                long_gaps.push(i);
            }
        }

        if long_gaps.len() < 2 {
            return; // need at least two long gaps for simultaneous ins + del
        }

        // Step 2 — sweep long‐gap positions with a sliding window,
        // accumulating insertion and deletion.  When a window's diff exceeds
        // the threshold, mark the spanned seed indices for removal.
        //
        // This mirrors minimap2's greedy scan: it flushes the current best
        // region when the sweep index passes its end, then continues looking
        // for more non‐overlapping bad regions.
        let n = long_gaps.len();
        let mut remove = vec![false; chain.len()];

        let mut best_diff: i64 = 0;
        let mut best_range: Option<(usize, usize)> = None; // (st, en) in long_gaps indices

        for k in 0..=n {
            // Flush when we've passed the current best region's end (or at
            // the very end of the long‐gap array).
            let flush = k == n || best_range.is_none_or(|(_, en)| k >= en);
            if flush {
                if let Some((st, en)) = best_range {
                    if best_diff > threshold {
                        for idx in long_gaps[st]..long_gaps[en] {
                            remove[idx] = true;
                        }
                    }
                }
                best_diff = 0;
                best_range = None;
                if k == n {
                    break;
                }
            }

            // Accumulate across a forward window of up to MAX_WINDOW
            // long‐gap positions starting at k.
            let i = long_gaps[k];
            let gap = chain[i].read_pos as i64 - chain[i - 1].read_end() as i64
                - (chain[i].ref_pos as i64 - chain[i - 1].ref_end() as i64);

            let mut n_ins = gap.max(0);
            let mut n_del = (-gap).max(0);

            for l in (k + 1)..n.min(k + 1 + MAX_WINDOW) {
                let j = long_gaps[l];
                let g = chain[j].read_pos as i64 - chain[j - 1].read_end() as i64
                    - (chain[j].ref_pos as i64 - chain[j - 1].ref_end() as i64);

                if g > 0 {
                    n_ins += g;
                } else {
                    n_del += -g;
                }

                let diff = n_ins + n_del - (n_ins - n_del).abs(); // = 2 * min(ins, del)
                if diff > threshold && diff > best_diff {
                    best_diff = diff;
                    best_range = Some((k, l));
                }
            }
        }

        // Step 3 — remove flagged seeds.
        if remove.iter().any(|&r| r) {
            let n_removed = remove.iter().filter(|&&r| r).count();
            log::debug!(
                "filter_misplaced_seeds: removing {} of {} seeds (threshold {})",
                n_removed,
                chain.len(),
                threshold,
            );
            let mut idx = 0;
            chain.retain(|_| {
                let keep = !remove[idx];
                idx += 1;
                keep
            });
        }
    }

    /// Filter seeds in regions where the diagonal shifts frequently relative
    /// to the reference span ("jittery" regions).
    ///
    /// Uses a sliding window of `window_size` inter-seed gaps. Within each
    /// window position the *shift density* is computed as:
    ///
    /// ```text
    /// density = sum(|diag[i] - diag[i-1]|) / ref_span
    /// ```
    ///
    /// where `ref_span` is the reference distance from the first seed's
    /// start to the last seed's end in the window. Seed lengths naturally
    /// contribute to `ref_span` without contributing shift, so long seeds
    /// dilute the density while short repeat-zone seeds with large diagonal
    /// jumps amplify it.
    ///
    /// Seeds within any window that exceeds `density_threshold` are marked
    /// for removal. This strips misleading anchors from tandem-repeat
    /// regions and lets the gap-fill DP aligner bridge between the
    /// flanking stable anchors.
    ///
    /// A `density_threshold` of 0.0 disables the filter.
    pub(crate) fn filter_jittery_seeds(chain: &mut Vec<SeedHit>, density_threshold: f64, window_size: usize) {
        let window_size = window_size.max(2); // need at least 2 gaps

        if chain.len() < window_size + 1 || density_threshold <= 0.0 {
            return;
        }

        // Precompute per-gap diagonal shifts.
        let diags: Vec<i64> = chain.iter().map(|s| s.diagonal).collect();
        let shifts: Vec<i64> = diags.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
        // shifts[i] = |diag[i+1] - diag[i]|, length = chain.len()-1

        let n_gaps = shifts.len();
        if n_gaps < window_size {
            return;
        }

        // Phase 1: classify each gap as "high-shift" based on a per-gap
        // criterion.  A gap is high-shift if its diagonal shift exceeds
        // what we'd expect from the seed length:
        //
        //   shift[i] > density_threshold × match_len[i]
        //
        // This prevents stable gaps flanked by long seeds (shift 1-3bp,
        // seed 50-100bp) from being caught in the net.
        let high_shift: Vec<bool> = (0..n_gaps)
            .map(|i| {
                let min_seed_len = chain[i].match_len.min(chain[i + 1].match_len).max(1);
                shifts[i] as f64 > density_threshold * min_seed_len as f64
            })
            .collect();

        // Phase 2: find contiguous runs of high-shift gaps.  For each
        // run of ≥ window_size gaps, verify the aggregate shift-density
        // (total shift / ref span) exceeds the threshold, then remove
        // the interior seeds.
        let mut remove = vec![false; chain.len()];

        let mut g = 0;
        while g < n_gaps {
            if !high_shift[g] {
                g += 1;
                continue;
            }
            let run_start = g;
            while g < n_gaps && high_shift[g] {
                g += 1;
            }
            let run_end = g; // exclusive, covers gaps [run_start..run_end)

            if run_end - run_start < window_size {
                continue;
            }

            // Aggregate density check: total shift / ref span of the run.
            let total_shift: i64 = shifts[run_start..run_end].iter().sum();
            let first_ref = chain[run_start].ref_pos;
            let last_ref_end = chain[run_end].ref_end();
            let ref_span = last_ref_end.saturating_sub(first_ref);
            if ref_span == 0 {
                continue;
            }
            let density = total_shift as f64 / ref_span as f64;
            if density <= density_threshold {
                continue;
            }

            // Remove interior seeds of this run.  Keep seed[run_start]
            // and seed[run_end] as boundary anchors for the DP aligner.
            for i in (run_start + 1)..run_end {
                remove[i] = true;
            }
        }

        // Remove flagged seeds.
        if remove.iter().any(|&r| r) {
            let n_removed = remove.iter().filter(|&&r| r).count();
            log::debug!(
                "filter_jittery_seeds: removing {} of {} seeds (density_threshold {:.2}, window {})",
                n_removed,
                chain.len(),
                density_threshold,
                window_size,
            );
            let mut idx = 0;
            chain.retain(|_| {
                let keep = !remove[idx];
                idx += 1;
                keep
            });
        }
    }

    pub fn diagonal(&self) -> f64 {
        let mut sum = 0i64;
        let mut total_length = 0usize;
        for seed in &self.chain {
            sum += seed.diagonal * seed.match_len as i64;
            total_length += seed.match_len;
        }
        if total_length == 0 {
            0.0
        } else {
            sum as f64 / total_length as f64
        }
    }

    pub fn total_identity(&self) -> usize {
        let total_seed_length: usize = self.chain.iter().map(|h| h.match_len).sum();
        let total_aligned_identities: usize = self
            .gap_alignments
            .iter()
            .map(|align| align.total_identity())
            .sum();
        total_seed_length + total_aligned_identities
    }

    /// Calculate fraction of read covered by this cluster
    pub fn read_coverage(&self, read_len: usize) -> f64 {
        if read_len == 0 {
            return 0.0;
        }
        (self.read_end - self.read_start) as f64 / read_len as f64
    }

    /// Total seed match length in this cluster
    pub fn total_seed_length(&self) -> usize {
        self.chain.iter().map(|h| h.match_len).sum()
    }

    /// Split this cluster at the specified gap, returning the new cluster
    /// containing seeds after the gap.
    ///
    /// `gap_seed_idx` is the index of the seed before the gap (as returned by `gaps()`).
    /// After splitting:
    /// - `self` contains seeds 0..=gap_seed_idx
    /// - The returned cluster contains seeds gap_seed_idx+1..
    ///
    /// Returns `None` if the split would leave either cluster empty.
    pub fn split_at_gap(&mut self, gap_seed_idx: usize) -> Option<(SeedCluster, Alignment)> {
        // Need at least one seed on each side
        if gap_seed_idx + 1 >= self.chain.len() {
            return None;
        }

        // Split the chain
        let tail_chain = self.chain.split_off(gap_seed_idx + 1);

        // Update self's read_end
        self.read_end = self
            .chain
            .last()
            .map(|h| h.read_end())
            .unwrap_or(self.read_start);

        // Split gap_alignments if populated
        // gap_seed_idx corresponds to the gap between seeds gap_seed_idx and gap_seed_idx+1
        // After split: self keeps gaps 0..gap_seed_idx-1, tail gets gaps gap_seed_idx+1..
        // The gap at gap_seed_idx is discarded (it's the split point)
        let tail_gap_alignments = if self.gap_alignments.len() > gap_seed_idx + 1 {
            self.gap_alignments.split_off(gap_seed_idx + 1)
        } else {
            Vec::new()
        };
        // Remove the gap alignment at gap_seed_idx (the split point)
        // After split_off, self has alignments 0..=gap_seed_idx, but needs 0..gap_seed_idx-1
        let dropped_alignment = self.gap_alignments.pop().unwrap();

        // Build the new cluster from the tail
        let tail_read_start = tail_chain.first().map(|h| h.read_pos).unwrap_or(0);
        let tail_read_end = tail_chain.last().map(|h| h.read_end()).unwrap_or(0);

        Some((
            SeedCluster {
                read_start: tail_read_start,
                read_end: tail_read_end,
                chain: tail_chain,
                is_reverse: self.is_reverse,
                chrom_id: self.chrom_id,
                gap_alignments: tail_gap_alignments,
                split_fill_tags: self.split_fill_tags.clone(),
            },
            dropped_alignment,
        ))
    }

    /// Fraction of the chain covered by seeds
    pub fn seed_density(&self) -> f64 {
        let total_length = if let Some(first) = self.chain.first() {
            if let Some(last) = self.chain.last() {
                last.ref_end() - first.ref_pos
            } else {
                0
            }
        } else {
            0
        };
        if total_length == 0 {
            0.0
        } else {
            self.total_seed_length() as f64 / total_length as f64
        }
    }

    /// Write this chain to a debug SAM file with SA tags linking all seeds.
    ///
    /// Each seed in the chain is written as a SAM record. The primary alignment
    /// is the first seed, and all others are marked as supplementary (0x800).
    /// All records include an SA tag listing all other alignments in the chain.
    ///
    /// # Arguments
    /// * `read_name` - Read name for SAM output
    /// * `cluster_id` - Cluster identifier for the cc tag
    /// * `chrom_name` - Chromosome name
    /// * `strand_seq` - Sequence for this strand
    /// * `strand_qual` - Quality scores (optional)
    pub fn write_chain_sam(
        &self,
        read_name: &str,
        cluster_id: usize,
        chrom_name: &str,
        strand_seq: &[u8],
        strand_qual: &[u8],
    ) {
        if self.chain.is_empty() {
            return;
        }

        let read_len = strand_seq.len();
        let strand_char = if self.is_reverse { '-' } else { '+' };

        // Build SA tag entries for all seeds in the chain
        // Format: rname,pos,strand,CIGAR,mapQ,NM;
        let sa_entries: Vec<String> = self
            .chain
            .iter()
            .map(|seed| {
                let cigar = format!("{}M", seed.match_len);
                let mapq = 60 / seed.kmer_uniqueness.max(1) as u8;
                format!(
                    "{},{},{},{},{},0",
                    chrom_name,
                    seed.ref_pos + 1,
                    strand_char,
                    cigar,
                    mapq
                )
            })
            .collect();

        // Write each seed as a SAM record
        for (i, seed) in self.chain.iter().enumerate() {
            // Build SA tag excluding the current seed
            let sa_others: Vec<&str> = sa_entries
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, s)| s.as_str())
                .collect();
            let sa_tag = if sa_others.is_empty() {
                String::new()
            } else {
                format!("\tSA:Z:{};", sa_others.join(";"))
            };

            // Flag: primary (0) for first seed, supplementary (0x800) for rest
            // Plus reverse flag (0x10) if reverse strand
            let mut flag = if i == 0 {
                0u16
            } else {
                FLAG_SUPPLEMENTARY
            };
            if self.is_reverse {
                flag |= FLAG_REVERSE;
            }

            // Build CIGAR with hard clips
            let hclip_start = seed.read_pos;
            let hclip_end = read_len.saturating_sub(seed.read_pos + seed.match_len);
            let cigar = match (hclip_start > 0, hclip_end > 0) {
                (true, true) => format!("{}H{}={}H", hclip_start, seed.match_len, hclip_end),
                (true, false) => format!("{}H{}=", hclip_start, seed.match_len),
                (false, true) => format!("{}={}H", seed.match_len, hclip_end),
                (false, false) => format!("{}=", seed.match_len),
            };

            let mapq = 60 / seed.kmer_uniqueness.max(1) as u8;

            // Extract aligned sequence and quality
            let seq_slice = &strand_seq[seed.read_pos..seed.read_pos + seed.match_len];
            let seq_str = String::from_utf8_lossy(seq_slice);

            let qual_str = {
                let qual_slice = &strand_qual[seed.read_pos..seed.read_pos + seed.match_len];
                qual_slice.iter().map(|&q| q as char).collect::<String>()
            };

            // Write SAM line with SA tag and cluster ID tag
            CHAINS_SAM.append(
                &format!(
                    "{}.{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}{}\tcc:i:{}",
                    read_name,
                    cluster_id,
                    flag,
                    chrom_name,
                    seed.ref_pos + 1,
                    mapq,
                    cigar,
                    seq_str,
                    qual_str,
                    sa_tag,
                    cluster_id
                ),
            );
        }
    }

    pub fn quality(&self, params: &AlignParams) -> QualityScore {
        let mut s = 0.0;
        for (i, seed) in self.chain.iter().enumerate() {
            let op = Op::new(Kind::SequenceMatch, seed.match_len);
            let score = params.quality(op).0;
            log::debug!("Seed {}: length = {}, score = {}", i, seed.match_len, score);
            s += score;
        }
        for (i, aln) in self.gap_alignments.iter().enumerate() {
            let a_len = aln.query_length();
            let a_score = aln.quality(params).0;
            log::debug!(
                "Gap {}: alignment length = {}, score = {}",
                i,
                a_len,
                a_score
            );
            s += a_score;
        }
        QualityScore::from(s)
    }

    /// Lightweight quality estimate that doesn't require gap alignments.
    ///
    /// Uses seed match scores (same as `quality()`) plus a geometric estimate
    /// for each gap between consecutive seeds. The gap estimate assumes the
    /// alignable portion (min of query/ref gap) scores as matches, and the
    /// length difference is scored as an indel.
    ///
    /// This is an optimistic upper bound on the real gap alignment quality,
    /// but preserves relative ordering well enough for covering-set selection.
    pub fn estimated_quality(&self, params: &AlignParams) -> QualityScore {
        let mut s = 0.0;
        for seed in self.chain.iter() {
            let op = Op::new(Kind::SequenceMatch, seed.match_len);
            s += params.quality(op).0;
        }
        // Estimate gap scores from geometry
        for pair in self.chain.windows(2) {
            let query_gap = pair[1].read_pos.saturating_sub(pair[0].read_end());
            let ref_gap = pair[1].ref_pos.saturating_sub(pair[0].ref_end());
            if query_gap == 0 && ref_gap == 0 {
                continue;
            }
            let min_len = query_gap.min(ref_gap);
            let indel_len = query_gap.abs_diff(ref_gap);
            // Optimistic: assume alignable bases are matches
            if min_len > 0 {
                s += params.quality(Op::new(Kind::SequenceMatch, min_len)).0;
            }
            // Indel penalty for the length difference
            if indel_len > 0 {
                s += params.quality(Op::new(Kind::Deletion, indel_len)).0;
            }
        }
        QualityScore::from(s)
    }

    #[allow(dead_code)]
    pub fn cigar_string(&self, read_len: usize) -> String {
        let mut cigar_ops = Vec::new();

        if let Some(first_seed) = self.chain.first() {
            // Handle leading hard clip
            let (fwd_start, _) = first_seed.fwd_read_range(read_len, self.is_reverse);
            if fwd_start > 0 {
                cigar_ops.push(format!("{}H", fwd_start));
            }
        }

        for (i, seed) in self.chain.iter().enumerate() {
            // Add the alignment before this seed if not the first seed
            if i > 0 {
                let gap_aln = &self.gap_alignments[i - 1];
                // Use the gap alignment CIGAR
                cigar_ops.push(gap_aln.cigar_string());
            }

            // Add match operation for this seed
            cigar_ops.push(format!("{}=", seed.match_len));
        }

        if let Some(last_seed) = self.chain.last() {
            // Handle trailing hard clip
            let (_, fwd_end) = last_seed.fwd_read_range(read_len, self.is_reverse);
            if fwd_end < read_len {
                let hclip_len = read_len - fwd_end;
                cigar_ops.push(format!("{}H", hclip_len));
            }
        }

        cigar_ops.join("")
    }

    /// Return an iterator over gaps in the seed coverage that are
    /// at least as big as the given size. Each gap is (start, end) in forward-strand
    /// read coordinates, along with the index of the seed before the gap.
    pub fn gaps(
        &self,
        read_len: usize,
        min_gap_length: usize,
    ) -> impl Iterator<Item = ((usize, usize), usize)> + '_ {
        self.chain
            .windows(2)
            .enumerate()
            .filter_map(move |(i, pair)| {
                let gap_start = pair[0].read_end();
                let gap_end = pair[1].read_pos;
                if gap_end > gap_start && gap_end - gap_start >= min_gap_length {
                    let fwd_read_start = if self.is_reverse {
                        read_len - gap_end
                    } else {
                        gap_start
                    };
                    let fwd_read_end = if self.is_reverse {
                        read_len - gap_start
                    } else {
                        gap_end
                    };
                    Some(((fwd_read_start, fwd_read_end), i))
                } else {
                    None
                }
            })
    }

    /// Align across all gaps between seeds using block aligner.
    ///
    /// Populates `gap_alignments` with one entry per gap (chain.len() - 1 entries).
    /// Each entry is `Some(alignment)` if alignment succeeded, `None` otherwise.
    ///
    /// # Arguments
    /// * `read_seq` - The read sequence (strand-specific, already rev-comped if reverse)
    /// * `ref_seq` - The reference sequence for this chromosome
    /// * `min_seed_length` - Minimum seed length for creating new clusters if splitting is needed
    ///
    /// This should be called before gap analysis to enable alignment-aware
    /// split decisions.
    pub fn align_gaps(
        mut self,
        read_name: &str,
        read_seq: &[u8],
        ref_seq: &[u8],
        min_seed_length: usize,
        aligner: &mut Aligner,
    ) -> Vec<Self> {
        if self.chain.len() < 2 {
            // No gaps to align
            return vec![self];
        }

        let num_gaps = self.chain.len() - 1;

        let mut gap_alignments = Vec::with_capacity(num_gaps);
        for (i, pair) in self.chain.windows(2).enumerate() {
            let seed1 = &pair[0];
            let seed2 = &pair[1];

            // Extract gap regions
            let read_gap_start = seed1.read_end();
            let read_gap_end = seed2.read_pos;
            let ref_gap_start = seed1.ref_end();
            let ref_gap_end = seed2.ref_pos;

            assert!(read_gap_start <= read_gap_end);
            assert!(ref_gap_start <= ref_gap_end);

            // Check for valid gap (one of the sequences must have a gap)
            if read_gap_end == read_gap_start && ref_gap_end == ref_gap_start {
                // No gap or negative gap - shouldn't happen after overlap resolution
                panic!(
                    "Invalid gap between seeds {} and {}: no gap in either read or reference",
                    i,
                    i + 1
                );
            }

            let read_gap = &read_seq[read_gap_start..read_gap_end];
            let ref_gap = &ref_seq[ref_gap_start..ref_gap_end];

            let alignment = aligner.align(read_gap, ref_gap);
            if alignment.is_none() {
                log::debug!(
                    "Gap alignment failed for read {}, after seed {}, read_gap_len={}, ref_gap_len={}",
                    read_name,
                    i,
                    read_gap.len(),
                    ref_gap.len()
                );
                if true || (read_gap.len() < 1000 && ref_gap.len() < 1000) {
                    // Write in FASTA format with descriptive headers
                    GAP_ALIGNMENTS.append(
                        &format!(
                            ">{}:read:{}-{}\n{}\n>{}:ref:{}-{}\n{}",
                            read_name,
                            read_gap_start,
                            read_gap_end,
                            String::from_utf8_lossy(read_gap),
                            read_name,
                            ref_gap_start,
                            ref_gap_end,
                            String::from_utf8_lossy(ref_gap)
                        ),
                    );
                }
            }
            gap_alignments.push(alignment);
        }

        if !gap_alignments.iter().any(|a| a.is_none()) {
            self.gap_alignments = gap_alignments.into_iter().map(|a| a.unwrap()).collect();
            return vec![self];
        }

        let groups: Vec<Vec<Alignment>> = gap_alignments.into_iter().groups().collect();

        let mut clusters = Vec::new();
        let mut seed_idx = 0;

        for group in groups {
            // Each group of n alignments connects n+1 seeds
            let num_seeds = group.len() + 1;
            let segment_seeds: Vec<SeedHit> = self.chain[seed_idx..seed_idx + num_seeds].to_vec();

            if let Some(mut cluster) =
                SeedCluster::new(segment_seeds, self.is_reverse, min_seed_length)
            {
                cluster.gap_alignments = group;
                clusters.push(cluster);
            }

            seed_idx += num_seeds;
        }

        clusters
    }

    /// Get the alignment for a specific gap by index.
    ///
    /// `gap_idx` is the index of the seed before the gap (0-based).
    /// Returns `None` if gaps haven't been aligned yet or if the gap index is invalid.
    pub fn gap_alignment(&self, gap_idx: usize) -> Option<&Alignment> {
        self.gap_alignments.get(gap_idx)
    }

    /// Format a visual diagram of seeds and gaps in this cluster.
    ///
    /// Returns a pair of strings showing the query and reference views:
    /// ```text
    /// QRY: [- 65bp -] <-  8bp -> [- 111bp -] <- 19bp -> [- 44bp -]
    /// REF: [- 65bp -] <- 11bp -> [- 111bp -] <- 17bp -> [- 44bp -]
    /// ```
    ///
    /// Seeds are shown as `[- Nbp -]` and have the same width on both lines
    /// since they represent exact matches. Gaps are shown as `<- Nbp ->`
    /// with padding to keep seeds aligned between the two lines.
    #[allow(dead_code)]
    pub fn format_seed_diagram(&self) -> (String, String) {
        if self.chain.is_empty() {
            return ("QRY:".to_string(), "REF:".to_string());
        }

        let mut qry_parts: Vec<String> = Vec::new();
        let mut ref_parts: Vec<String> = Vec::new();

        for (i, seed) in self.chain.iter().enumerate() {
            // Add the gap before this seed (if not the first seed)
            if i > 0 {
                let prev = &self.chain[i - 1];
                let qry_gap = seed.read_pos.saturating_sub(prev.read_end()) as i64;
                let ref_gap = seed.ref_pos.saturating_sub(prev.ref_end()) as i64;

                // Format gap strings and find the max width for alignment
                let qry_gap_str = format!("{}bp", qry_gap);
                let ref_gap_str = format!("{}bp", ref_gap);
                let max_width = qry_gap_str.len().max(ref_gap_str.len());

                // Pad to align
                qry_parts.push(format!(
                    " <- {:>width$} -> ",
                    qry_gap_str,
                    width = max_width
                ));
                ref_parts.push(format!(
                    " <- {:>width$} -> ",
                    ref_gap_str,
                    width = max_width
                ));
            }

            // Add the seed (same width on both lines)
            let seed_str = format!("[- {}bp -]", seed.match_len);
            qry_parts.push(seed_str.clone());
            ref_parts.push(seed_str);
        }

        let qry_line = format!("QRY: {}", qry_parts.join(""));
        let ref_line = format!("REF: {}", ref_parts.join(""));

        (qry_line, ref_line)
    }

    #[allow(dead_code)]
    pub fn validate(
        &self,
        read_seq: &[u8],
        ref_seq: &[u8],
    ) -> std::result::Result<(), ClusterError> {
        // Validate each seed
        for (_i, seed) in self.chain.iter().enumerate() {
            let seed_read = &read_seq[seed.read_pos..seed.read_end()];
            let seed_ref = &ref_seq[seed.ref_pos..seed.ref_end()];
            if seed_read != seed_ref {
                return Err(ClusterError::SequenceMismatch {
                    read_bases: String::from_utf8_lossy(seed_read).to_string(),
                    ref_bases: String::from_utf8_lossy(seed_ref).to_string(),
                    read_pos: seed.read_pos,
                    ref_pos: seed.ref_pos,
                });
            }
        }
        for (i, aln) in self.gap_alignments.iter().enumerate() {
            let seed1 = &self.chain[i];
            let seed2 = &self.chain[i + 1];
            let read_gap_start = seed1.read_end();
            let read_gap_end = seed2.read_pos;
            let ref_gap_start = seed1.ref_end();
            let ref_gap_end = seed2.ref_pos;

            let read_gap = &read_seq[read_gap_start..read_gap_end];
            let ref_gap = &ref_seq[ref_gap_start..ref_gap_end];

            if let Err(err) = aln.validate(ref_gap, read_gap, 0) {
                return Err(ClusterError::AlignmentMismatch {
                    gap_index: i,
                    read_start: read_gap_start,
                    read_end: read_gap_end,
                    ref_start: ref_gap_start,
                    ref_end: ref_gap_end,
                    error: err,
                });
            }
        }
        Ok(())
    }
}
