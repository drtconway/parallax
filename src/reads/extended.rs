use ordered_float::OrderedFloat;

use crate::{
    align::{Alignment, DpAligner},
    config,
    reads::seeds::SeedHit,
    reference::InMemoryReference,
    utils::sequence::complement,
};

/// Extended seeds with additional metadata for weighted interval scheduling and chaining.
/// NB these seeds are always interpreted as forward strand, with is_reverse flag indicating
/// if they came from the reverse complement.
///
/// If is_reverse is true, the read_start is still the position on the original read, but the
/// ref_start is the position on the reference where the reverse complement of the read matches.
///
/// This allows us to treat all seeds uniformly in the chaining and scheduling steps, while still
/// retaining the necessary information to construct the final alignment correctly.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ExtendedSeed {
    read_start: usize,
    length: usize,
    ref_chrom_id: usize,
    ref_start: usize,
    multiplicity: usize,
    is_reverse: bool,
    weight: OrderedFloat<f64>,
}

impl ExtendedSeed {
    pub fn from_seed_hit(seed: &SeedHit, is_reverse: bool, read_len: usize) -> Self {
        let read_start = if is_reverse {
            read_len - seed.read_pos - seed.match_len
        } else {
            seed.read_pos
        };
        let weight = Self::calculate_weight(seed.match_len as f64, seed.kmer_uniqueness as f64);
        Self {
            read_start,
            length: seed.match_len,
            ref_chrom_id: seed.chrom_id,
            ref_start: seed.ref_pos,
            multiplicity: seed.kmer_uniqueness as usize,
            is_reverse,
            weight: OrderedFloat(weight),
        }
    }

    pub fn weight(&self) -> f64 {
        self.weight.0
    }

    fn calculate_weight(n: f64, m: f64) -> f64 {
        // The weight of a seed is a function of its length and multiplicity.
        // Longer seeds are more informative.
        // Higher multiplicity seeds are less informative, so we want to penalise them,
        // but we take the length into account since very long seeds are probably more
        // unique than the multiplicity suggests.

        const ALPHA: f64 = 0.25;
        const BETA: f64 = 0.5;
        const GAMMA: f64 = 0.25;

        let log_n = n.log10();
        let log_m = m.log10();
        n * (1.0 + ALPHA * log_n) / (1.0 + (BETA * log_m) / (1.0 + GAMMA * log_n))
    }

    pub fn multiplicity(&self) -> usize {
        self.multiplicity
    }

    pub fn read_start(&self) -> usize {
        self.read_start
    }

    pub fn read_end(&self) -> usize {
        self.read_start + self.length
    }

    pub fn length(&self) -> usize {
        self.length
    }

    pub fn ref_chrom_id(&self) -> usize {
        self.ref_chrom_id
    }

    pub fn ref_start(&self) -> usize {
        self.ref_start
    }

    pub fn is_reverse(&self) -> bool {
        self.is_reverse
    }

    /// Convert the seed to an alignment record.
    /// Seeds are exact matches by construction (k-mer matches merged and
    /// extended along exact-match diagonals), so the CIGAR is a single
    /// `=` (SequenceMatch) run covering the full length.
    pub fn to_alignment(&self) -> Alignment {
        Alignment::from_perfect_match(self.length)
    }

    /// Take a vector of extended seeds and merges that are overlapping on the read and reference into single seeds.
    /// Merges seeds that lie on exactly the same diagonal and overlap on the read.
    ///
    /// Diagonal is defined per-strand:
    /// - Forward: `ref_start - read_start` (constant along a gapless match)
    /// - Reverse: `ref_start + read_start` (constant because ref decreases as read advances)
    ///
    /// Seeds are temporarily re-sorted by (chrom, strand, diagonal, read_start)
    /// so that mergeable seeds are adjacent, then merged with a linear scan,
    /// and finally re-sorted into the natural derived ordering.
    pub fn simplify_seeds(seeds: &mut Vec<ExtendedSeed>) {
        if seeds.len() <= 1 {
            return;
        }

        // Sort by (chrom, strand, diagonal, read_start) to make mergeable
        // seeds adjacent.  We encode the diagonal as an isize to handle
        // forward-strand cases where ref_start < read_start.
        seeds.sort_by(|a, b| {
            a.ref_chrom_id
                .cmp(&b.ref_chrom_id)
                .then(a.is_reverse.cmp(&b.is_reverse))
                .then(Self::diagonal(a).cmp(&Self::diagonal(b)))
                .then(a.read_start.cmp(&b.read_start))
        });

        // Linear merge scan.
        let mut write = 0;
        for read in 1..seeds.len() {
            let curr_read_start = seeds[read].read_start;
            let curr_length = seeds[read].length;
            let curr_ref_chrom_id = seeds[read].ref_chrom_id;
            let curr_is_reverse = seeds[read].is_reverse;
            let curr_multiplicity = seeds[read].multiplicity;
            let curr_diagonal = Self::diagonal(&seeds[read]);

            let prev = &seeds[write];
            let prev_diagonal = Self::diagonal(prev);

            let same_diagonal = prev.ref_chrom_id == curr_ref_chrom_id
                && prev.is_reverse == curr_is_reverse
                && prev_diagonal == curr_diagonal;

            let prev_read_end = prev.read_start + prev.length;
            let read_overlaps = curr_read_start < prev_read_end;

            if same_diagonal && read_overlaps {
                // Merge: extend to cover both on the read; ref follows from the diagonal.
                let new_read_end = prev_read_end.max(curr_read_start + curr_length);
                let new_length = new_read_end - seeds[write].read_start;
                let new_multiplicity = prev.multiplicity.min(curr_multiplicity);
                let new_weight = Self::calculate_weight(new_length as f64, new_multiplicity as f64);
                seeds[write].length = new_length;
                seeds[write].multiplicity = new_multiplicity;
                seeds[write].weight = OrderedFloat(new_weight);
                // ref_start: take the min (for reverse strand the earlier
                // read position has the higher ref_start, but we still want
                // the ref interval to cover the union).
                if seeds[read].ref_start < seeds[write].ref_start {
                    seeds[write].ref_start = seeds[read].ref_start;
                }
            } else {
                write += 1;
                if write != read {
                    seeds.swap(write, read);
                }
            }
        }
        seeds.truncate(write + 1);

        // Re-sort into the natural derived ordering.
        seeds.sort();
    }

    /// Compute the diagonal for a seed (same value for seeds on a gapless match).
    fn diagonal(s: &ExtendedSeed) -> isize {
        if s.is_reverse {
            s.ref_start as isize + s.read_start as isize
        } else {
            s.ref_start as isize - s.read_start as isize
        }
    }

    /// The edge penalty is the cost of chaining two seeds together, based on how far apart they are on the read and reference.
    /// The following principles apply:
    /// - Seeds that are close together (but not overlapping) on both the read and reference should have a low penalty, encouraging them to be chained together.
    /// - Seeds that are far apart on both the read and reference should have a high penalty, discouraging them from being chained together.
    /// - Seeds that are close together on the read but somewhat distant on the reference, but on the same strand and in congruent order (representing a deletion) should have a small penalty.
    /// - Seeds on different chromosomes, different strands, or in non-colinear order may mark an SV; these receive a moderate fixed reference-side penalty rather than being rejected.
    /// - Read overlap is never permitted.
    pub fn edge_penalty(&self, other: &ExtendedSeed) -> Option<f64> {
        // Allow a small overlap on the read (k-mer seeding can produce seeds
        // that straddle breakpoints), but reject large overlaps.
        const READ_OVERLAP_TOLERANCE: usize = 5;
        let self_read_end = self.read_start + self.length;
        let read_gap = if other.read_start >= self_read_end {
            (other.read_start - self_read_end) as f64
        } else {
            let overlap = self_read_end - other.read_start;
            if overlap > READ_OVERLAP_TOLERANCE {
                return None;
            }
            (overlap * overlap) as f64 // small overlap → small penalty
            // TODO: instead of an ad-hoc penalty, return the overlap so the
            // DP can deduct trimmed bases from the seed's weight.  Deferred
            // because the overlap is at most READ_OVERLAP_TOLERANCE bases and the
            // weight impact is negligible in practice.
        };

        // Fixed penalty applied when the reference side is discontinuous
        // (different chrom, different strand, or non-colinear).
        let sv_penalty = config::get().seeding.sv_penalty;

        // Maximum reference-vs-read deviation we treat as a simple indel.
        // Beyond this, the gap is more likely a rearrangement (e.g. two seeds
        // happen to be on the same strand but are megabases apart on the
        // reference).  Without this cap the logarithmic penalty would let
        // such pairs chain cheaply — ln(1 + 14M) ≈ 16.5, far less than
        // sv_penalty — causing the DP to prefer a spurious same-strand
        // seed over the correct cross-strand one.
        const MAX_INDEL_DEVIATION: f64 = 100_000.0;

        // Try to compute a colinear reference gap.  If the seeds are on the
        // same chromosome and strand and in the right order, we get a
        // non-negative gap; otherwise we fall back to the SV penalty.
        /// Maximum ref-coordinate overlap we still treat as a small indel
        /// rather than an SV.  Matches the tolerance in `is_colinear`.
        const REF_OVERLAP_TOLERANCE: i64 = 10;

        let ref_penalty =
            if self.ref_chrom_id != other.ref_chrom_id || self.is_reverse != other.is_reverse {
                sv_penalty
            } else {
                let ref_gap = if self.is_reverse {
                    // Reverse strand: reference positions decrease as read advances.
                    let other_ref_end = (other.ref_start + other.length) as i64;
                    self.ref_start as i64 - other_ref_end
                } else {
                    // Forward strand: reference positions increase as read advances.
                    let self_ref_end = (self.ref_start + self.length) as i64;
                    other.ref_start as i64 - self_ref_end
                };

                if ref_gap < -REF_OVERLAP_TOLERANCE {
                    sv_penalty // too large overlap — structural variant
                } else {
                    let deviation = (ref_gap as f64 - read_gap).abs();
                    if deviation > MAX_INDEL_DEVIATION {
                        sv_penalty
                    } else {
                        (1.0 + deviation).ln()
                    }
                }
            };

        Some(read_gap + ref_penalty)
    }

    /// Form explanatory groups by greedy peeling.
    ///
    /// Each group is a complete alternative explanation of the read — a chain of
    /// non-overlapping seeds (on the read) that maximises total weight minus
    /// edge penalties.  Successive groups are extracted by removing consumed
    /// seeds and repeating the DP until the best remaining chain falls below
    /// `MIN_GROUP_WEIGHT` or `MAX_GROUPS` is reached.
    pub fn form_explanatory_groups(seeds: &[ExtendedSeed]) -> Vec<Vec<ExtendedSeed>> {
        const MIN_GROUP_WEIGHT: f64 = 50.0;
        const MAX_GROUPS: usize = 10;

        let mut groups: Vec<Vec<ExtendedSeed>> = Vec::new();
        if seeds.is_empty() {
            return groups;
        }

        // Build read_start–sorted index (seeds may already be in natural
        // order, but we sort explicitly to be safe).
        let mut order: Vec<usize> = (0..seeds.len()).collect();
        order.sort_by_key(|&i| seeds[i].read_start);

        let mut available = vec![true; seeds.len()];

        for g in 0..MAX_GROUPS {
            // Collect indices of seeds not yet consumed, in read_start order.
            let active: Vec<usize> = order.iter().copied().filter(|&i| available[i]).collect();
            if active.is_empty() {
                break;
            }

            let n = active.len();
            let mut dp = vec![0.0f64; n];
            let mut pred = vec![usize::MAX; n];

            for i in 0..n {
                let seed_i = &seeds[active[i]];
                dp[i] = seed_i.weight();

                for j in (0..i).rev() {
                    if let Some(penalty) = seeds[active[j]].edge_penalty(seed_i) {
                        let score = dp[j] + seed_i.weight() - penalty;
                        if score > dp[i] {
                            dp[i] = score;
                            pred[i] = j;
                        }
                    }
                }
            }

            // Find the best endpoint.
            let best = (0..n)
                .max_by(|&a, &b| {
                    dp[a]
                        .partial_cmp(&dp[b])
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap();
            if dp[best] < MIN_GROUP_WEIGHT {
                break;
            }

            log::debug!("Group {g}: best score = {:.2}", dp[best]);

            // Traceback to extract the chain.
            let mut chain = Vec::new();
            let mut cur = best;
            loop {
                let seed_idx = active[cur];
                chain.push(seeds[seed_idx].clone());
                available[seed_idx] = false;
                if pred[cur] == usize::MAX {
                    break;
                }
                cur = pred[cur];
            }
            chain.reverse();
            groups.push(chain);
        }

        for (g, group) in groups.iter().enumerate() {
            let n = group.len();
            for (i, seed) in group.iter().enumerate() {
                let strand = if seed.is_reverse { '-' } else { '+' };
                let edge = if i < n - 1 {
                    if let Some(penalty) = seed.edge_penalty(&group[i + 1]) {
                        format!("{:.2}", -penalty) // show edge score (negative penalty)
                    } else {
                        "NA".to_string()
                    }
                } else {
                    "NA".to_string()
                };
                log::debug!(
                    "{g}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
                    seed.read_start,
                    seed.read_start + seed.length,
                    seed.length,
                    seed.ref_chrom_id,
                    seed.ref_start,
                    strand,
                    seed.weight(),
                    edge
                );
            }
        }

        groups
    }

    /// Extend this seed rightward on the read by up to `limit` bases using
    /// exact matching against the reference.  For a forward-strand seed this
    /// walks rightward on the reference; for reverse-strand it walks leftward.
    fn extend_right(&mut self, limit: usize, read_seq: &[u8], reference: &InMemoryReference) {
        if limit == 0 {
            return;
        }
        let read_start = self.read_start + self.length;
        if self.is_reverse {
            let ref_end = self.ref_start;
            let ref_begin = ref_end.saturating_sub(limit);
            let ref_seq = reference.get_seq(self.ref_chrom_id, ref_begin, ref_end);
            let mut ext = 0;
            for k in 0..limit.min(ref_seq.len()) {
                let ref_idx = ref_seq.len() - 1 - k;
                if complement(read_seq[read_start + k]) == ref_seq[ref_idx] {
                    ext += 1;
                } else {
                    break;
                }
            }
            if ext > 0 {
                self.length += ext;
                self.ref_start -= ext;
            }
        } else {
            let ref_start = self.ref_start + self.length;
            let ref_seq = reference.get_seq(self.ref_chrom_id, ref_start, ref_start + limit);
            let mut ext = 0;
            for k in 0..limit.min(ref_seq.len()) {
                if read_seq[read_start + k] == ref_seq[k] {
                    ext += 1;
                } else {
                    break;
                }
            }
            if ext > 0 {
                self.length += ext;
            }
        }
    }

    /// Extend this seed leftward on the read by up to `limit` bases using
    /// exact matching against the reference.  For a reverse-strand seed this
    /// walks rightward on the reference; for forward-strand it walks leftward.
    fn extend_left(&mut self, limit: usize, read_seq: &[u8], reference: &InMemoryReference) {
        if limit == 0 {
            return;
        }
        if self.is_reverse {
            let ref_start = self.ref_start + self.length;
            let ref_seq = reference.get_seq(self.ref_chrom_id, ref_start, ref_start + limit);
            let mut ext = 0;
            for k in 0..limit.min(ref_seq.len()) {
                let read_pos = self.read_start - 1 - k;
                if complement(read_seq[read_pos]) == ref_seq[k] {
                    ext += 1;
                } else {
                    break;
                }
            }
            if ext > 0 {
                self.read_start -= ext;
                self.length += ext;
            }
        } else {
            let ref_end = self.ref_start;
            let ref_begin = ref_end.saturating_sub(limit);
            let ref_seq = reference.get_seq(self.ref_chrom_id, ref_begin, ref_end);
            let mut ext = 0;
            for k in 0..limit.min(ref_seq.len()) {
                let ref_idx = ref_seq.len() - 1 - k;
                if read_seq[self.read_start - 1 - k] == ref_seq[ref_idx] {
                    ext += 1;
                } else {
                    break;
                }
            }
            if ext > 0 {
                self.read_start -= ext;
                self.length += ext;
                self.ref_start -= ext;
            }
        }
    }

    /// Trim `n` bases from the right end of this seed on the read.
    fn trim_right(&mut self, n: usize) {
        let n = n.min(self.length);
        self.length -= n;
        if self.is_reverse {
            // Right end of read = left end of ref.
            self.ref_start += n;
        }
    }

    /// Trim `n` bases from the left end of this seed on the read.
    fn trim_left(&mut self, n: usize) {
        let n = n.min(self.length);
        self.read_start += n;
        self.length -= n;
        if !self.is_reverse {
            // Left end of read = left end of ref.
            self.ref_start += n;
        }
    }

    /// Extend seeds to fill gaps and trim overlaps within a single group.
    ///
    /// Adjacent seeds in the group (ordered by read position) may have small
    /// gaps or overlaps on the read.  This function resolves both:
    ///
    /// **Gaps**: one seed is extended along its reference diagonal using exact
    /// matching against the read, stopping at the first mismatch or the start
    /// of the next seed.
    ///
    /// **Overlaps**: one seed is trimmed so the two seeds abut without
    /// overlapping on the read.
    ///
    /// Which seed is extended/trimmed is determined by a strand-aware rule that
    /// left-aligns breakpoints on the chromosome:
    ///
    /// - **Both forward**: extend the left seed rightward on read (= rightward
    ///   on ref).  Trim the right seed's left end if needed.
    /// - **Both reverse**: extend the right seed leftward on read (= rightward
    ///   on ref).  Trim the left seed's right end if needed.
    /// - **Opposite strands**: extend the left seed rightward on read (= left-
    ///   to-right on read).  Trim the right seed's left end if needed.
    ///
    /// This ensures that reads from either strand spanning the same SV
    /// breakpoint with microhomology place the breakpoint at the same
    /// chromosomal position.
    pub fn extend_and_trim(
        group: &mut Vec<ExtendedSeed>,
        read_seq: &[u8],
        reference: &InMemoryReference,
    ) {
        if group.len() <= 1 {
            return;
        }

        for i in 0..group.len() - 1 {
            let a_read_end = group[i].read_start + group[i].length;
            let b_read_start = group[i + 1].read_start;

            if a_read_end <= b_read_start {
                // Gap: extend one seed to fill it.
                let gap = b_read_start - a_read_end;
                if gap == 0 {
                    continue;
                }

                match (group[i].is_reverse, group[i + 1].is_reverse) {
                    (false, false) | (true, false) | (false, true) => {
                        // Extend A rightward on read.
                        group[i].extend_right(gap, read_seq, reference);
                    }
                    (true, true) => {
                        // Both reverse: extend B leftward on read (= rightward on ref).
                        group[i + 1].extend_left(gap, read_seq, reference);
                    }
                }
            } else {
                // Overlap: trim one seed to remove it.
                let overlap = a_read_end - b_read_start;

                match (group[i].is_reverse, group[i + 1].is_reverse) {
                    (true, true) => {
                        // Both reverse: trim A's right end on read.
                        group[i].trim_right(overlap);
                    }
                    (false, false) | (true, false) | (false, true) => {
                        // Trim B's left end on read.
                        group[i + 1].trim_left(overlap);
                    }
                }
            }
        }

        // Remove any seeds that got trimmed to zero length.
        group.retain(|s| s.length > 0);

        // ── Resolve reference overlaps ─────────────────────────────────
        //
        // After extending/trimming on the read, adjacent same-chrom
        // same-strand seeds may share reference bases.  This can happen
        // when:
        //   (a) The original seeds already overlapped on ref (e.g. two
        //       k-mer seeds flanking a 1-base insertion share a ref base
        //       due to microhomology).
        //   (b) An extension above pushed a seed past the adjacent seed's
        //       ref boundary.
        //
        // Trimming the overlap creates a small read gap that align_gaps()
        // will fill with the appropriate insertion CIGAR.
        for i in 0..group.len().saturating_sub(1) {
            if group[i].ref_chrom_id != group[i + 1].ref_chrom_id {
                continue;
            }
            if group[i].is_reverse != group[i + 1].is_reverse {
                continue;
            }

            if group[i].is_reverse {
                // Reverse strand: b.ref_end should be ≤ a.ref_start.
                let b_ref_end = group[i + 1].ref_start + group[i + 1].length;
                if b_ref_end > group[i].ref_start {
                    let overlap = b_ref_end - group[i].ref_start;
                    // trim_left on reverse: shrinks read-left / ref-right end
                    group[i + 1].trim_left(overlap);
                }
            } else {
                // Forward strand: a.ref_end should be ≤ b.ref_start.
                let a_ref_end = group[i].ref_start + group[i].length;
                if a_ref_end > group[i + 1].ref_start {
                    let overlap = a_ref_end - group[i + 1].ref_start;
                    // trim_right on forward: shrinks read-right / ref-right end
                    group[i].trim_right(overlap);
                }
            }
        }

        // Remove any seeds that got trimmed to zero length (from ref overlap resolution).
        group.retain(|s| s.length > 0);
    }

    /// Test whether two adjacent seeds (in read order) are colinear — i.e. on
    /// the same chromosome, same strand, in the correct reference order, and
    /// close enough that the gap is a simple indel rather than a rearrangement.
    fn is_colinear(&self, other: &ExtendedSeed) -> bool {
        const MAX_INDEL_DEVIATION: f64 = 100_000.0;
        /// Maximum ref-coordinate overlap (negative gap) we still treat as
        /// colinear.  Small overlaps arise when adjacent k-mer seeds share a
        /// reference base at a microhomology / small-indel boundary.
        const REF_OVERLAP_TOLERANCE: i64 = 10;

        if self.ref_chrom_id != other.ref_chrom_id || self.is_reverse != other.is_reverse {
            return false;
        }

        let self_read_end = self.read_start + self.length;
        let read_gap = if other.read_start >= self_read_end {
            (other.read_start - self_read_end) as f64
        } else {
            return false; // overlapping on read — not a simple gap
        };

        let ref_gap = if self.is_reverse {
            // Reverse strand: reference positions decrease as read advances.
            let other_ref_end = (other.ref_start + other.length) as i64;
            self.ref_start as i64 - other_ref_end
        } else {
            // Forward strand: reference positions increase as read advances.
            let self_ref_end = (self.ref_start + self.length) as i64;
            other.ref_start as i64 - self_ref_end
        };

        if ref_gap < -REF_OVERLAP_TOLERANCE {
            return false; // too large an overlap — not a simple indel
        }

        let deviation = (ref_gap as f64 - read_gap).abs();
        deviation <= MAX_INDEL_DEVIATION
    }

    /// Align the gaps between adjacent seeds in a group.
    ///
    /// Returns a vector of `n - 1` entries (one per gap between consecutive
    /// seeds), where each entry is:
    ///
    /// - `Some(alignment)` for a colinear gap that can be bridged with DP
    ///   alignment (same chrom, same strand, correct order, reasonable size).
    /// - `None` for a structural-variant gap (different chrom, different
    ///   strand, non-colinear, or too large).
    ///
    /// The query for each alignment is the read subsequence in the gap; the
    /// reference is the corresponding genomic interval.  For reverse-strand
    /// seeds the reference is reverse-complemented so the alignment is always
    /// in read-forward orientation.
    pub fn align_gaps(
        group: &[ExtendedSeed],
        read_seq: &[u8],
        reference: &InMemoryReference,
        aligner: &mut DpAligner,
    ) -> Vec<Option<Alignment>> {
        if group.len() <= 1 {
            return Vec::new();
        }

        let mut alignments = Vec::with_capacity(group.len() - 1);

        for i in 0..group.len() - 1 {
            let a = &group[i];
            let b = &group[i + 1];

            if !a.is_colinear(b) {
                alignments.push(None);
                continue;
            }

            // Read subsequence in the gap.
            let a_read_end = a.read_start + a.length;
            let b_read_start = b.read_start;
            let query = &read_seq[a_read_end..b_read_start];

            // Reference subsequence in the gap.
            //
            // After extend_and_trim resolves ref overlaps, ref_begin ≤
            // ref_end should always hold.  We clamp defensively so a
            // residual overlap doesn't panic.
            let ref_slice = if a.is_reverse {
                // Reverse strand: ref positions decrease as read advances.
                // Gap on ref is [b.ref_start + b.length .. a.ref_start).
                let ref_begin = b.ref_start + b.length;
                let ref_end = a.ref_start;
                if ref_begin >= ref_end {
                    Vec::new()
                } else {
                    let fwd = reference.get_seq(a.ref_chrom_id, ref_begin, ref_end);
                    // Reverse-complement so it aligns in read-forward orientation.
                    fwd.iter()
                        .rev()
                        .map(|&base| complement(base))
                        .collect::<Vec<u8>>()
                }
            } else {
                // Forward strand: ref positions increase as read advances.
                let ref_begin = a.ref_start + a.length;
                let ref_end = b.ref_start;
                if ref_begin >= ref_end {
                    Vec::new()
                } else {
                    reference
                        .get_seq(a.ref_chrom_id, ref_begin, ref_end)
                        .to_vec()
                }
            };

            alignments.push(aligner.align(query, &ref_slice));
        }

        alignments
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to construct an ExtendedSeed without needing a SeedHit.
    fn seed(
        read_start: usize,
        length: usize,
        chrom: usize,
        ref_start: usize,
        is_reverse: bool,
    ) -> ExtendedSeed {
        let weight = ExtendedSeed::calculate_weight(length as f64, 1.0);
        ExtendedSeed {
            read_start,
            length,
            ref_chrom_id: chrom,
            ref_start,
            multiplicity: 1,
            is_reverse,
            weight: OrderedFloat(weight)
        }
    }

    // ── simplify_seeds ──────────────────────────────────────────────────

    #[test]
    fn simplify_no_overlap() {
        let mut seeds = vec![seed(0, 10, 0, 100, false), seed(20, 10, 0, 120, false)];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn simplify_merge_two_overlapping() {
        // Overlap on both read [0..10) & [5..15) and ref [100..110) & [105..115)
        let mut seeds = vec![seed(0, 10, 0, 100, false), seed(5, 10, 0, 105, false)];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].read_start, 0);
        assert_eq!(seeds[0].length, 15);
        assert_eq!(seeds[0].ref_start, 100);
    }

    #[test]
    fn simplify_merge_three_overlapping() {
        let mut seeds = vec![
            seed(0, 10, 0, 100, false),
            seed(5, 10, 0, 105, false),
            seed(10, 10, 0, 110, false),
        ];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].read_start, 0);
        assert_eq!(seeds[0].length, 20);
        assert_eq!(seeds[0].ref_start, 100);
    }

    #[test]
    fn simplify_no_merge_different_chrom() {
        let mut seeds = vec![seed(0, 10, 0, 100, false), seed(5, 10, 1, 105, false)];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn simplify_no_merge_different_strand() {
        let mut seeds = vec![seed(0, 10, 0, 100, false), seed(5, 10, 0, 105, true)];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn simplify_no_merge_read_overlap_but_no_ref_overlap() {
        // Overlap on read but ref regions are disjoint.
        let mut seeds = vec![seed(0, 10, 0, 100, false), seed(5, 10, 0, 200, false)];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn simplify_reverse_strand_overlap() {
        // Reverse strand: both seeds map to overlapping ref and read intervals.
        let mut seeds = vec![seed(0, 10, 0, 100, true), seed(5, 10, 0, 95, true)];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].read_start, 0);
        assert_eq!(seeds[0].length, 15);
        assert_eq!(seeds[0].ref_start, 95);
    }

    #[test]
    fn simplify_keeps_min_multiplicity() {
        let w_10_5 = ExtendedSeed::calculate_weight(10.0, 5.0);
        let w_10_2 = ExtendedSeed::calculate_weight(10.0, 2.0);
        let mut seeds = vec![
            ExtendedSeed {
                read_start: 0,
                length: 10,
                ref_chrom_id: 0,
                ref_start: 100,
                multiplicity: 5,
                is_reverse: false,
                weight: OrderedFloat(w_10_5)
            },
            ExtendedSeed {
                read_start: 5,
                length: 10,
                ref_chrom_id: 0,
                ref_start: 105,
                multiplicity: 2,
                is_reverse: false,
                weight: OrderedFloat(w_10_2)
            },
        ];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].multiplicity, 2);
    }

    #[test]
    fn simplify_multi_hit_partial_merge() {
        // Two seeds overlap on the read: [0..10) and [5..15).
        // The second seed has two reference hits: ref 105 (overlaps with
        // the first seed's ref [100..110)) and ref 500 (no ref overlap).
        //
        // Sorted order: (0,10,0,100), (5,10,0,105), (5,10,0,500)
        //   - first pair: read overlap ✓, ref overlap [100..110)∩[105..115) ✓ → merge
        //   - merged (0,15,0,100) vs (5,10,0,500): read overlap ✓, ref overlap? no → keep
        //
        // Expected: 2 seeds after merging.
        let mut seeds = vec![
            seed(0, 10, 0, 100, false),
            seed(5, 10, 0, 105, false),
            seed(5, 10, 0, 500, false),
        ];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 2, "seeds: {seeds:?}");

        // First seed is the merged one covering read [0..15), ref starting at 100.
        assert_eq!(seeds[0].read_start, 0);
        assert_eq!(seeds[0].length, 15);
        assert_eq!(seeds[0].ref_start, 100);

        // Second seed is the unmerged hit at ref 500.
        assert_eq!(seeds[1].read_start, 5);
        assert_eq!(seeds[1].length, 10);
        assert_eq!(seeds[1].ref_start, 500);
    }

    #[test]
    fn simplify_non_adjacent_in_natural_order() {
        // Seeds that should merge are separated in natural sort order by
        // a seed on a different diagonal.
        //
        // Natural order: (0,10,0,100), (0,10,0,200), (5,10,0,105)
        // Diagonals:       100           200           100
        //
        // After diagonal sort: (0,10,0,100), (5,10,0,105), (0,10,0,200)
        //   → first two merge (same diagonal 100, read overlap)
        //   → (0,10,0,200) stays (different diagonal)
        //
        // Expected: 2 seeds.
        let mut seeds = vec![
            seed(0, 10, 0, 100, false),
            seed(0, 10, 0, 200, false),
            seed(5, 10, 0, 105, false),
        ];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 2, "seeds: {seeds:?}");

        // After re-sort into natural order (read_start, length, ...):
        // (0,10,0,200) then (0,15,0,100) — shorter length sorts first.
        assert_eq!(seeds[0].read_start, 0);
        assert_eq!(seeds[0].length, 10);
        assert_eq!(seeds[0].ref_start, 200);
        assert_eq!(seeds[1].read_start, 0);
        assert_eq!(seeds[1].length, 15);
        assert_eq!(seeds[1].ref_start, 100);
    }

    #[test]
    fn simplify_empty_and_single() {
        let mut empty: Vec<ExtendedSeed> = vec![];
        ExtendedSeed::simplify_seeds(&mut empty);
        assert!(empty.is_empty());

        let mut single = vec![seed(0, 10, 0, 100, false)];
        ExtendedSeed::simplify_seeds(&mut single);
        assert_eq!(single.len(), 1);
    }

    // ── edge_penalty ────────────────────────────────────────────────────

    #[test]
    fn penalty_large_read_overlap_returns_none() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(2, 10, 0, 110, false);
        assert!(a.edge_penalty(&b).is_none());
    }

    #[test]
    fn penalty_small_read_overlap_is_tolerated() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(9, 10, 0, 110, false); // 1bp overlap
        assert!(a.edge_penalty(&b).is_some());
    }

    #[test]
    fn penalty_adjacent_colinear_is_small() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 110, false);
        let p = a.edge_penalty(&b).unwrap();
        assert!(p < 1.0, "expected small penalty, got {p}");
    }

    #[test]
    fn penalty_far_apart_is_large() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(1000, 10, 0, 1100, false);
        let p = a.edge_penalty(&b).unwrap();
        assert!(p > 900.0, "expected large penalty, got {p}");
    }

    #[test]
    fn penalty_deletion_is_moderate() {
        // Close on read, far on reference = deletion.
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 50_110, false);
        let p = a.edge_penalty(&b).unwrap();
        // read_gap = 0, ref_gap = 50000, so penalty = 0 + ln(1 + 50000) ≈ 10.8
        assert!(p < 15.0, "expected moderate penalty for deletion, got {p}");
    }

    #[test]
    fn penalty_different_chrom_uses_sv_penalty() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 1, 100, false);
        let p = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_different_strand_uses_sv_penalty() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 100, true);
        let p = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_non_colinear_fwd_uses_sv_penalty() {
        // Forward strand but ref goes backwards.
        let a = seed(0, 10, 0, 200, false);
        let b = seed(10, 10, 0, 100, false);
        let p = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_colinear_reverse_strand() {
        // Reverse strand: ref positions decrease as read advances.
        let a = seed(0, 10, 0, 200, true);
        let b = seed(10, 10, 0, 180, true);
        // ref gap = 200 - (180 + 10) = 10, read gap = 0.
        let p = a.edge_penalty(&b).unwrap();
        assert!(
            p < 5.0,
            "expected small penalty for colinear reverse, got {p}"
        );
    }

    #[test]
    fn penalty_non_colinear_reverse_uses_sv_penalty() {
        // Reverse strand but ref increases — non-colinear.
        let a = seed(0, 10, 0, 100, true);
        let b = seed(10, 10, 0, 200, true);
        let p = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    // ── is_colinear ─────────────────────────────────────────────────────

    #[test]
    fn colinear_forward_abutting() {
        // Forward strand, abutting on both read and ref → colinear.
        let a = seed(0, 100, 0, 1000, false);
        let b = seed(100, 100, 0, 1100, false);
        assert!(a.is_colinear(&b));
    }

    #[test]
    fn colinear_forward_with_gap() {
        // Forward strand, 10bp gap on both read and ref → colinear.
        let a = seed(0, 100, 0, 1000, false);
        let b = seed(110, 100, 0, 1110, false);
        assert!(a.is_colinear(&b));
    }

    #[test]
    fn colinear_forward_1bp_ref_overlap() {
        // Forward strand, 1-base ref overlap: a covers ref [1000,1100),
        // b starts at ref 1099 → ref overlap of 1.
        // Seeds abut on read (read gap = 0).
        // This represents a 1bp insertion at the junction.
        let a = seed(0, 100, 0, 1000, false);
        let b = seed(100, 100, 0, 1099, false);
        assert!(a.is_colinear(&b), "1bp ref overlap should be colinear");
    }

    #[test]
    fn colinear_reverse_abutting() {
        // Reverse strand, abutting: ref gap = 0, read gap = 0.
        let a = seed(0, 100, 0, 500, true); // ref [500, 600)
        let b = seed(100, 100, 0, 400, true); // ref [400, 500)
        assert!(a.is_colinear(&b));
    }

    #[test]
    fn colinear_reverse_with_gap() {
        // Reverse strand, 10bp gap on both read and ref → colinear.
        let a = seed(0, 100, 0, 500, true); // ref [500, 600)
        let b = seed(110, 100, 0, 390, true); // ref [390, 490)
        // ref gap = 500 - 490 = 10, read gap = 10
        assert!(a.is_colinear(&b));
    }

    #[test]
    fn colinear_reverse_1bp_ref_overlap() {
        // Reverse strand, 1-base ref overlap:
        // a: ref [500,600), b: ref [401,501)
        // b.ref_end=501 > a.ref_start=500 → overlap of 1.
        // Seeds abut on read (read gap = 0).
        // This represents a 1bp insertion at the junction.
        let a = seed(0, 100, 0, 500, true);
        let b = seed(100, 100, 0, 401, true);
        assert!(a.is_colinear(&b), "1bp ref overlap should be colinear");
    }

    #[test]
    fn colinear_reverse_3bp_ref_overlap() {
        // Reverse strand, 3-base ref overlap (also seen in real data).
        let a = seed(0, 100, 0, 500, true); // ref [500, 600)
        let b = seed(100, 100, 0, 403, true); // ref [403, 503), overlap = 3
        assert!(a.is_colinear(&b), "3bp ref overlap should be colinear");
    }

    #[test]
    fn colinear_reverse_large_ref_overlap_rejected() {
        // Reverse strand, 20-base ref overlap → too large → not colinear.
        let a = seed(0, 100, 0, 500, true); // ref [500, 600)
        let b = seed(100, 100, 0, 420, true); // ref [420, 520), overlap = 20
        assert!(!a.is_colinear(&b), "large ref overlap should be rejected");
    }

    // ── edge_penalty with small ref overlaps ─────────────────────────

    #[test]
    fn penalty_reverse_1bp_ref_overlap_is_small() {
        // 1-base ref overlap → deviation is 1 → penalty = ln(2) ≈ 0.69.
        let a = seed(0, 100, 0, 500, true);
        let b = seed(100, 100, 0, 401, true);
        let p = a.edge_penalty(&b).unwrap();
        let expected = (2.0f64).ln(); // ln(1 + 1) = ln(2)
        assert!(
            (p - expected).abs() < 0.01,
            "expected ~{expected:.2}, got {p}"
        );
    }

    #[test]
    fn penalty_forward_1bp_ref_overlap_is_small() {
        let a = seed(0, 100, 0, 1000, false);
        let b = seed(100, 100, 0, 1099, false);
        let p = a.edge_penalty(&b).unwrap();
        let expected = (2.0f64).ln();
        assert!(
            (p - expected).abs() < 0.01,
            "expected ~{expected:.2}, got {p}"
        );
    }
}
