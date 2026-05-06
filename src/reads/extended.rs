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
    kmer_uniqueness: u32,
    read_frequency: u32,
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
            kmer_uniqueness: seed.kmer_uniqueness,
            read_frequency: seed.read_frequency,
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

    pub fn kmer_uniqueness(&self) -> u32 {
        self.kmer_uniqueness
    }

    pub fn read_frequency(&self) -> u32 {
        self.read_frequency
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

    pub fn ref_end(&self) -> usize {
        self.ref_start + self.length
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
                .then(a.diagonal().cmp(&b.diagonal()))
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
            let curr_kmer_uniqueness = seeds[read].kmer_uniqueness;
            let curr_read_frequency = seeds[read].read_frequency;
            let curr_diagonal = seeds[read].diagonal();

            let prev = &seeds[write];
            let prev_diagonal = prev.diagonal();

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
                let new_kmer_uniqueness = prev.kmer_uniqueness.min(curr_kmer_uniqueness);
                let new_read_frequency = prev.read_frequency.max(curr_read_frequency);
                let new_weight = Self::calculate_weight(new_length as f64, new_multiplicity as f64);
                seeds[write].length = new_length;
                seeds[write].multiplicity = new_multiplicity;
                seeds[write].kmer_uniqueness = new_kmer_uniqueness;
                seeds[write].read_frequency = new_read_frequency;
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

    /// Identify seeds that are immediate diagonal excursions using the minimap2
    /// neighbour heuristic.
    ///
    /// Seeds must be sorted in increasing read-start order.  `sv_breaks[i]` marks
    /// a break between `seeds[i]` and `seeds[i+1]`; seeds adjacent to a break are
    /// skipped.
    ///
    /// For each interior seed (not adjacent to an SV break, same chrom/strand as
    /// both neighbours), if the diagonal shifts to the immediate predecessor and
    /// successor both exceed `threshold` in magnitude with opposite signs, the seed
    /// is flagged.
    ///
    /// Returns a `Vec<bool>` of length `seeds.len()` where `true` means flagged.
    pub fn find_immediate_diagonal_excursions(
        seeds: &[ExtendedSeed],
        sv_breaks: &[bool],
        threshold: isize,
    ) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];

        if n < 3 {
            return flagged;
        }

        for pos in 1..n - 1 {
            if sv_breaks[pos - 1] || sv_breaks[pos] {
                continue;
            }

            let seed = &seeds[pos];
            let prev = &seeds[pos - 1];
            let next = &seeds[pos + 1];

            if prev.ref_chrom_id != seed.ref_chrom_id || prev.is_reverse != seed.is_reverse {
                continue;
            }
            if next.ref_chrom_id != seed.ref_chrom_id || next.is_reverse != seed.is_reverse {
                continue;
            }

            let diag = seed.diagonal();
            let left_shift = prev.diagonal() - diag;
            let right_shift = next.diagonal() - diag;

            if left_shift.abs() > threshold
                && right_shift.abs() > threshold
                && left_shift.signum() != right_shift.signum()
            {
                flagged[pos] = true;
            }
        }

        flagged
    }

    /// Identify seeds whose weight-to-read-frequency ratio is low and that shift
    /// the diagonal from an adjacent seed.
    ///
    /// Seeds must be sorted in increasing read-start order.  `sv_breaks[i]` marks
    /// a break between `seeds[i]` and `seeds[i+1]`; seeds adjacent to a break are
    /// skipped.
    ///
    /// A seed is flagged if `weight / read_frequency < min_weight_per_frequency` AND
    /// its diagonal deviates from at least one immediate neighbour (on the same
    /// chrom/strand, not separated by an SV break) by more than `threshold`.
    ///
    /// Returns a `Vec<bool>` of length `seeds.len()` where `true` means flagged.
    pub fn find_high_read_frequency_excursions(
        seeds: &[ExtendedSeed],
        sv_breaks: &[bool],
        threshold: isize,
        min_weight_per_frequency: f64,
    ) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];

        if n < 2 {
            return flagged;
        }

        for pos in 0..n {
            // Skip seeds adjacent to an SV break on either side — removing them
            // could destroy a break anchor.
            let has_left_sv = pos == 0 || sv_breaks[pos - 1];
            let has_right_sv = pos == n - 1 || sv_breaks[pos];
            if has_left_sv || has_right_sv {
                continue;
            }

            let seed = &seeds[pos];
            if seed.weight() / seed.read_frequency() as f64 >= min_weight_per_frequency {
                continue;
            }

            let diag = seed.diagonal();

            // Check left neighbour (same chrom/strand guaranteed by no SV break).
            let left_shift = if seeds[pos - 1].ref_chrom_id == seed.ref_chrom_id
                && seeds[pos - 1].is_reverse == seed.is_reverse
            {
                Some((seeds[pos - 1].diagonal() - diag).abs())
            } else {
                None
            };

            // Check right neighbour (same chrom/strand guaranteed by no SV break).
            let right_shift = if seeds[pos + 1].ref_chrom_id == seed.ref_chrom_id
                && seeds[pos + 1].is_reverse == seed.is_reverse
            {
                Some((seeds[pos + 1].diagonal() - diag).abs())
            } else {
                None
            };

            let deviates = left_shift.map_or(false, |s| s > threshold)
                || right_shift.map_or(false, |s| s > threshold);

            if deviates {
                flagged[pos] = true;
            }
        }

        flagged
    }

    /// Prune seeds identified as diagonal excursions, iterating until stable.
    ///
    /// On each pass, calls `detect` to identify seeds to remove, then removes
    /// them and updates `sv_breaks` in sync.  Repeats until no seeds are removed.
    pub fn prune_repetitive_seeds(
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<bool>,
        threshold: isize,
        reference: Option<&InMemoryReference>,
    ) {
        // First pass: high read-frequency seeds that shift the diagonal, iterated to convergence.
        loop {
            let flagged =
                Self::find_high_read_frequency_excursions(seeds, sv_breaks, threshold, 50.0);
            if flagged.iter().all(|&f| !f) {
                break;
            }
            Self::log_removals(&flagged, seeds, reference);
            Self::remove_flagged(seeds, sv_breaks, &flagged);
        }

        // Second pass: immediate-neighbour heuristic, iterated to convergence.
        loop {
            let flagged = Self::find_immediate_diagonal_excursions(seeds, sv_breaks, threshold);
            if flagged.iter().all(|&f| !f) {
                break;
            }
            Self::log_removals(&flagged, seeds, reference);
            Self::remove_flagged(seeds, sv_breaks, &flagged);
        }

        if false {
            // Second pass: windowed median heuristic, catches runs that protected
            // each other from the immediate-neighbour test.
            let flagged =
                Self::find_transient_diagonal_excursions(seeds, sv_breaks, 450, threshold);
            if flagged.iter().any(|&f| f) {
                Self::log_removals(&flagged, seeds, reference);
                Self::remove_flagged(seeds, sv_breaks, &flagged);
            }
        }

        if true {
            // Third pass: terminal trim — trim seeds at segment ends whose diagonal
            // deviates from the segment's global weighted median.  These cannot be
            // caught by the windowed tests because there are no good seeds on the
            // outer side to establish a reference median.
            let flagged = Self::find_terminal_diagonal_excursions(seeds, sv_breaks, threshold);
            if flagged.iter().any(|&f| f) {
                Self::log_removals(&flagged, seeds, reference);
                Self::remove_flagged(seeds, sv_breaks, &flagged);
            }
        }
    }

    /// Remove flagged seeds and update `sv_breaks` in sync.
    /// When seed at pos is removed, its two flanking breaks are merged with OR.
    fn remove_flagged(seeds: &mut Vec<ExtendedSeed>, sv_breaks: &mut Vec<bool>, flagged: &[bool]) {
        let n = seeds.len();
        let mut new_sv_breaks = Vec::with_capacity(sv_breaks.len());
        let mut pending: Option<bool> = None;
        for pos in 0..n {
            if pos < n - 1 {
                let b = sv_breaks[pos];
                if !flagged[pos] {
                    new_sv_breaks.push(pending.take().map_or(b, |p| p || b));
                } else {
                    pending = Some(pending.map_or(b, |p| p || b));
                }
            }
        }
        *sv_breaks = new_sv_breaks;

        let mut i = 0;
        seeds.retain(|_| {
            let result = !flagged[i];
            i += 1;
            result
        });
    }

    fn log_removals(
        flagged: &[bool],
        seeds: &[ExtendedSeed],
        reference: Option<&InMemoryReference>,
    ) {
        if let Some(r) = reference {
            for (pos, seed) in seeds.iter().enumerate() {
                if flagged[pos] {
                    log::debug!(
                        "removing badly placed seed {}-{} {}:{}-{} ({}) mult={} ku={} rf={}",
                        seed.read_start(),
                        seed.read_end(),
                        r.chrom_name(seed.ref_chrom_id()),
                        seed.ref_start(),
                        seed.ref_end(),
                        if seed.is_reverse() { "-" } else { "+" },
                        seed.multiplicity(),
                        seed.kmer_uniqueness(),
                        seed.read_frequency(),
                    );
                }
            }
        }
    }

    /// Identify seeds that are transient diagonal excursions within a colinear segment.
    ///
    /// Seeds must be sorted in increasing read-start order.  `sv_breaks[i]` marks
    /// a permanent step change between `seeds[i]` and `seeds[i+1]`; seeds are only
    /// evaluated within the colinear segment they belong to.
    ///
    /// For each seed, the weighted median diagonal of the `min_window_len` bases
    /// immediately before it and immediately after it (within the same segment) are
    /// computed.  If both windows agree with each other (within `threshold`) but the
    /// seed's own diagonal disagrees with both, the seed is a transient excursion —
    /// it jumped away from the prevailing diagonal and returned — and is flagged.
    ///
    /// Returns a `Vec<bool>` of length `seeds.len()` where `true` means the seed
    /// should be removed.

    /// Identify seeds at the start or end of a colinear segment whose diagonal
    /// deviates from the global weighted median of the segment.
    ///
    /// For each segment, compute the weighted median diagonal over all seeds.
    /// Then trim seeds from each end while they deviate from that median by more
    /// than `threshold`.  This catches terminal clusters of bad seeds that have
    /// no good flanking seeds on the outer side to anchor the windowed test.
    pub fn find_terminal_diagonal_excursions(
        seeds: &[ExtendedSeed],
        sv_breaks: &[bool],
        threshold: isize,
    ) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];

        if n < 3 {
            return flagged;
        }

        let global_weighted_median = |range: std::ops::Range<usize>| -> Option<isize> {
            let mut pairs: Vec<(isize, f64)> = range
                .map(|i| (seeds[i].diagonal(), seeds[i].weight()))
                .collect();
            if pairs.is_empty() {
                return None;
            }
            pairs.sort_by_key(|&(d, _)| d);
            let total: f64 = pairs.iter().map(|(_, w)| w).sum();
            let mut cumulative = 0f64;
            for (d, w) in &pairs {
                cumulative += w;
                if cumulative * 2.0 >= total {
                    return Some(*d);
                }
            }
            pairs.last().map(|&(d, _)| d)
        };

        let mut seg_start = 0;
        loop {
            let seg_end = (seg_start..n - 1)
                .find(|&i| sv_breaks[i])
                .map(|i| i + 1)
                .unwrap_or(n);

            let seg_len = seg_end - seg_start;

            if seg_len >= 3 {
                if let Some(median) = global_weighted_median(seg_start..seg_end) {
                    // Trim from the left end.
                    for pos in seg_start..seg_end {
                        if (seeds[pos].diagonal() - median).abs() > threshold {
                            flagged[pos] = true;
                        } else {
                            break;
                        }
                    }
                    // Trim from the right end.
                    for pos in (seg_start..seg_end).rev() {
                        if flagged[pos] {
                            // Already flagged from the left pass; don't break.
                            continue;
                        }
                        if (seeds[pos].diagonal() - median).abs() > threshold {
                            flagged[pos] = true;
                        } else {
                            break;
                        }
                    }
                }
            }

            if seg_end == n {
                break;
            }
            seg_start = seg_end;
        }

        flagged
    }

    /// Identify seeds whose diagonal is a transient excursion within a colinear segment.
    pub fn find_transient_diagonal_excursions(
        seeds: &[ExtendedSeed],
        sv_breaks: &[bool],
        min_window_len: usize,
        threshold: isize,
    ) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];

        if n < 3 {
            return flagged;
        }

        // Weighted median of diagonals for a slice of seeds, weighted by seed weight.
        let weighted_median = |indices: &[usize]| -> Option<isize> {
            if indices.is_empty() {
                return None;
            }
            let mut pairs: Vec<(isize, f64)> = indices
                .iter()
                .map(|&i| (seeds[i].diagonal(), seeds[i].weight()))
                .collect();
            pairs.sort_by_key(|&(d, _)| d);
            let total: f64 = pairs.iter().map(|(_, w)| w).sum();
            let mut cumulative = 0f64;
            for (d, w) in &pairs {
                cumulative += w;
                if cumulative * 2.0 >= total {
                    return Some(*d);
                }
            }
            pairs.last().map(|&(d, _)| d)
        };

        // Find the segment boundaries (runs of seeds with no sv_break between them).
        // Process each segment independently.
        let mut seg_start = 0;
        loop {
            // Find the end of the current segment.
            let seg_end = (seg_start..n - 1)
                .find(|&i| sv_breaks[i])
                .map(|i| i + 1)
                .unwrap_or(n);

            let seg = seg_start..seg_end;
            let seg_len = seg.len();

            if seg_len >= 3 {
                for pos in seg_start..seg_end {
                    // Build left window: seeds before pos, accumulating weight up to
                    // min_window_len, staying within the segment.
                    let mut left_weight = 0f64;
                    let mut left_indices: Vec<usize> = Vec::new();
                    for j in (seg_start..pos).rev() {
                        left_indices.push(j);
                        left_weight += seeds[j].weight();
                        if left_weight >= min_window_len as f64 {
                            break;
                        }
                    }

                    // Build right window: seeds after pos, accumulating weight up to
                    // min_window_len, staying within the segment.
                    let mut right_weight = 0f64;
                    let mut right_indices: Vec<usize> = Vec::new();
                    for j in pos + 1..seg_end {
                        right_indices.push(j);
                        right_weight += seeds[j].weight();
                        if right_weight >= min_window_len as f64 {
                            break;
                        }
                    }

                    // Need at least one seed on each side to evaluate.
                    let (Some(left_med), Some(right_med)) = (
                        weighted_median(&left_indices),
                        weighted_median(&right_indices),
                    ) else {
                        continue;
                    };

                    let diag = seeds[pos].diagonal();

                    log::info!(
                        "pos = {}, left {}, right {}, window {}",
                        pos,
                        (diag - left_med).abs(),
                        (diag - right_med).abs(),
                        (left_med - right_med).abs()
                    );

                    // Transient: seed deviates from both windows, but the windows agree.
                    if (diag - left_med).abs() > threshold
                        && (diag - right_med).abs() > threshold
                        && (left_med - right_med).abs() <= threshold
                    {
                        flagged[pos] = true;
                    }
                }
            }

            if seg_end == n {
                break;
            }
            seg_start = seg_end;
        }

        flagged
    }

    /// Compute the diagonal for a seed (same value for seeds on a gapless match).
    pub fn diagonal(&self) -> isize {
        if self.is_reverse {
            // read_fwd[read_start + k] <-> ref[ref_start + length - 1 - k]
            // invariant: ref_start + length - 1 + read_start = const
            self.ref_start as isize + self.length as isize - 1 + self.read_start as isize
        } else {
            self.ref_start as isize - self.read_start as isize
        }
    }

    /// The edge penalty is the cost of chaining two seeds together, based on how far apart they are on the read and reference.
    /// The following principles apply:
    /// - Seeds that are close together (but not overlapping) on both the read and reference should have a low penalty, encouraging them to be chained together.
    /// - Seeds that are far apart on both the read and reference should have a high penalty, discouraging them from being chained together.
    /// - Seeds that are close together on the read but somewhat distant on the reference, but on the same strand and in congruent order (representing a deletion) should have a small penalty.
    /// - Seeds on different chromosomes, different strands, or in non-colinear order may mark an SV; these receive a moderate fixed reference-side penalty rather than being rejected.
    /// - Read overlap is never permitted.
    #[inline(never)]
    pub fn edge_penalty(&self, other: &ExtendedSeed) -> Option<(f64, bool)> {
        // Seeds with a large read overlap cannot be chained — the downstream
        // seed would contribute almost no new information.
        const MAX_READ_OVERLAP: usize = 50;

        let self_read_end = self.read_start + self.length;
        let read_overlap = if other.read_start < self_read_end {
            self_read_end - other.read_start
        } else {
            0
        };

        if read_overlap >= other.length || read_overlap > MAX_READ_OVERLAP {
            return None;
        }

        // Signed read gap: positive = bases between seeds, negative = overlap.
        // Used for deviation calculation against the reference gap.
        // The cost of a gap is added below; overlapping bases are free here
        // (their weight deduction is handled separately).
        let read_gap: f64 = if read_overlap == 0 {
            (other.read_start - self_read_end) as f64
        } else {
            -(read_overlap as f64)
        };
        let read_gap_cost = read_gap.max(0.0);

        let cfg = &config::get().seeding;
        let sv_penalty = cfg.sv_penalty;
        let threshold = cfg.gap_linear_threshold;
        let scale = cfg.gap_linear_scale;
        let quad = cfg.read_gap_quad_scale;

        // Quadratic read-gap cost: long unanchored read stretches are weighted
        // more than proportionally, since large insertions are better explained
        // as SV breakpoints than colinear gaps.
        let read_gap_cost = read_gap_cost + quad * read_gap_cost * read_gap_cost;

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
        const REF_OVERLAP_TOLERANCE: i64 = 10;

        let (ref_penalty, is_sv) =
            if self.ref_chrom_id != other.ref_chrom_id || self.is_reverse != other.is_reverse {
                (sv_penalty, true)
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
                    (sv_penalty, true)
                } else {
                    let deviation = (ref_gap as f64 - read_gap).abs();
                    if deviation > MAX_INDEL_DEVIATION {
                        (sv_penalty, true)
                    } else {
                        let log_part = (1.0 + deviation.min(threshold)).ln();
                        let linear_part = scale * (deviation - threshold).max(0.0);
                        (log_part + linear_part, false)
                    }
                }
            };

        Some((read_gap_cost + ref_penalty, is_sv))
    }

    /// Form explanatory groups by greedy peeling.
    ///
    /// Each group is a complete alternative explanation of the read — a chain of
    /// non-overlapping seeds (on the read) that maximises total weight minus
    /// edge penalties.  Successive groups are extracted by removing consumed
    /// seeds and repeating the DP until the best remaining chain falls below
    /// `MIN_GROUP_WEIGHT` or `MAX_GROUPS` is reached.
    #[inline(never)]
    pub fn form_explanatory_groups(seeds: &[ExtendedSeed]) -> Vec<(Vec<ExtendedSeed>, Vec<bool>)> {
        const MIN_GROUP_WEIGHT: f64 = 50.0;
        const MAX_GROUPS: usize = 10;
        // Stop early if the best remaining chain scores below this fraction of
        // group 0's score.  Chains this weak contribute negligibly to MAPQ and
        // are not worth the cost of additional DP iterations.
        const MIN_RELATIVE_SCORE: f64 = 0.05;

        let mut groups: Vec<(Vec<ExtendedSeed>, Vec<bool>)> = Vec::new();
        if seeds.is_empty() {
            return groups;
        }

        // Build read_start–sorted index (seeds may already be in natural
        // order, but we sort explicitly to be safe).
        let mut order: Vec<usize> = (0..seeds.len()).collect();
        order.sort_by_key(|&i| seeds[i].read_start);

        let mut available = vec![true; seeds.len()];
        let mut group0_score = 0.0f64;

        for g in 0..MAX_GROUPS {
            // Collect indices of seeds not yet consumed, in read_start order.
            let active: Vec<usize> = order.iter().copied().filter(|&i| available[i]).collect();
            if active.is_empty() {
                break;
            }

            let n = active.len();
            let mut dp = vec![0.0f64; n];
            let mut pred = vec![usize::MAX; n];
            let mut pred_is_sv = vec![false; n];

            for i in 0..n {
                let seed_i = &seeds[active[i]];
                dp[i] = seed_i.weight();

                for j in (0..i).rev() {
                    if let Some((penalty, is_sv)) = seeds[active[j]].edge_penalty(seed_i) {
                        let score = dp[j] + seed_i.weight() - penalty;
                        if score > dp[i] {
                            dp[i] = score;
                            pred[i] = j;
                            pred_is_sv[i] = is_sv;
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
            let best_score = dp[best];
            if best_score < MIN_GROUP_WEIGHT {
                break;
            }
            if g == 0 {
                group0_score = best_score;
            } else if best_score < MIN_RELATIVE_SCORE * group0_score {
                break;
            }

            log::debug!("Group {g}: best score = {:.2}", best_score);

            // Traceback to extract the chain and SV-break flags.
            // sv_breaks[i] is true if there is an SV break between chain[i] and chain[i+1].
            let mut chain = Vec::new();
            let mut sv_breaks = Vec::new();
            let mut cur = best;
            loop {
                let seed_idx = active[cur];
                chain.push(seeds[seed_idx].clone());
                available[seed_idx] = false;
                if pred[cur] == usize::MAX {
                    break;
                }
                sv_breaks.push(pred_is_sv[cur]);
                cur = pred[cur];
            }
            chain.reverse();
            sv_breaks.reverse();
            groups.push((chain, sv_breaks));
        }

        std::hint::black_box(groups)
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
        sv_breaks: &mut Vec<bool>,
        read_seq: &[u8],
        reference: &InMemoryReference,
    ) {
        if group.len() <= 1 {
            return;
        }

        // Remove zero-length seeds and keep sv_breaks in sync.
        // When seed at pos is removed, merge its two flanking breaks with OR.
        let retain_nonzero = |group: &mut Vec<ExtendedSeed>, sv_breaks: &mut Vec<bool>| {
            let mut new_sv_breaks = Vec::with_capacity(sv_breaks.len());
            let mut pending: Option<bool> = None;
            for pos in 0..group.len() {
                if pos < group.len() - 1 {
                    let b = sv_breaks[pos];
                    if group[pos].length > 0 {
                        new_sv_breaks.push(pending.take().map_or(b, |p| p || b));
                    } else {
                        pending = Some(pending.map_or(b, |p| p || b));
                    }
                }
            }
            *sv_breaks = new_sv_breaks;
            group.retain(|s| s.length > 0);
        };

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

        retain_nonzero(group, sv_breaks);

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
        //
        // A single pass isn't enough: trimming seed i may expose a new
        // overlap between seed i and seed i+2 once seed i+1 is gone, or
        // trimming may leave a cascade.  Repeat until the group is stable.
        loop {
            let mut any_trimmed = false;
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
                        any_trimmed = true;
                    }
                } else {
                    // Forward strand: a.ref_end should be ≤ b.ref_start.
                    let a_ref_end = group[i].ref_start + group[i].length;
                    if a_ref_end > group[i + 1].ref_start {
                        let overlap = a_ref_end - group[i + 1].ref_start;
                        // trim_right on forward: shrinks read-right / ref-right end
                        group[i].trim_right(overlap);
                        any_trimmed = true;
                    }
                }
            }
            if !any_trimmed {
                break;
            }
            retain_nonzero(group, sv_breaks);
        }

        retain_nonzero(group, sv_breaks);

        // Re-evaluate SV breaks: after extension, trimming, and pruning, the
        // colinearity of adjacent seed pairs may have changed.  Recompute every
        // break unconditionally — both clearing breaks that are now simple indels
        // and setting breaks that are now SV-sized gaps (e.g. because an
        // intermediate bridging seed was pruned away).
        for i in 0..sv_breaks.len() {
            match group[i].edge_penalty(&group[i + 1]) {
                Some((_, is_sv)) => sv_breaks[i] = is_sv,
                None => sv_breaks[i] = true,
            }
        }
    }

    /// Align the gaps between adjacent seeds in a group.
    ///
    /// `sv_breaks[i]` must be `true` if the DP placed an SV breakpoint between
    /// `group[i]` and `group[i+1]`.  SV gaps are returned as `None` without
    /// attempting alignment; colinear gaps are aligned with DP.
    ///
    /// Returns a vector of `n - 1` entries (one per gap between consecutive
    /// seeds), where each entry is:
    ///
    /// - `Some(alignment)` for a colinear gap that can be bridged with DP
    ///   alignment (same chrom, same strand, correct order, reasonable size).
    /// - `None` for a structural-variant gap.
    ///
    /// The query for each alignment is the read subsequence in the gap; the
    /// reference is the corresponding genomic interval.  For reverse-strand
    /// seeds the reference is reverse-complemented so the alignment is always
    /// in read-forward orientation.
    pub fn align_gaps(
        group: &[ExtendedSeed],
        sv_breaks: &[bool],
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

            if sv_breaks[i] {
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

pub enum TagValue {
    Str(String),
    Int(i64),
}

pub struct ExtendedSeedDumpItem<'a> {
    reference: &'a InMemoryReference,
    read_id: &'a str,
    read_len: usize,
    seed_num: usize,
    seed: &'a ExtendedSeed,
    seq: &'a str,
    qual: &'a str,
    tags: Vec<(String, TagValue)>,
}

impl<'a> crate::utils::dump::DumpItem for ExtendedSeedDumpItem<'a> {
    type HeaderInfo = InMemoryReference;

    fn write_header(header_info: &Self::HeaderInfo, writer: &mut impl std::io::Write) {
        let res: std::io::Result<()> = (|| {
            writeln!(writer, "@HD\tVN:1.6")?;
            for chrom in header_info.all_chrom_info() {
                writeln!(writer, "@SQ\tSN:{}\tLN:{}", chrom.name, chrom.length)?;
            }
            Ok(())
        })();
        res.expect("writing header failed");
    }

    fn write(&self, writer: &mut impl std::io::Write) {
        let flag: u8 = if self.seed.is_reverse { 0x10 } else { 0x00 };
        let chrom = self.reference.chrom_name(self.seed.ref_chrom_id());
        let mapq = (self.seed.weight().floor() as i32).min(200);
        let read_left = self.seed.read_start();
        let len = self.seed.length();
        let read_right = self.read_len - self.seed.read_end();
        let (left, right) = if self.seed.is_reverse {
            (read_right, read_left)
        } else {
            (read_left, read_right)
        };
        let pos = self.seed.ref_start() + 1;
        let left_clip = if left > 0 {
            format!("{}H", left)
        } else {
            String::new()
        };
        let right_clip = if right > 0 {
            format!("{}H", right)
        } else {
            String::new()
        };
        let seq = if self.seed.is_reverse {
            self.seq
                .chars()
                .map(|c| match c {
                    'A' | 'a' => 'T',
                    'C' | 'c' => 'G',
                    'G' | 'g' => 'C',
                    'T' | 't' => 'A',
                    _ => c,
                })
                .rev()
                .collect()
        } else {
            String::from(self.seq)
        };
        write!(
            writer,
            "{}_{}\t{}\t{}\t{}\t{}\t{}{}={}\t*\t0\t0\t{}\t{}",
            self.read_id,
            self.seed_num,
            flag,
            chrom,
            pos,
            mapq,
            left_clip,
            len,
            right_clip,
            seq,
            self.qual
        )
        .expect("write failed");
        for (tag, value) in self.tags.iter() {
            match value {
                TagValue::Str(s) => {
                    write!(writer, "\t{}:Z:{}", tag, s).expect("write failed");
                }
                TagValue::Int(i) => {
                    write!(writer, "\t{}:i:{}", tag, i).expect("write failed");
                }
            }
        }
        writeln!(writer, "").expect("write failed");
    }
}

impl<'a>
    From<(
        &'a InMemoryReference,
        &'a str,
        usize,
        usize,
        &'a ExtendedSeed,
        &'a str,
        &'a str,
    )> for ExtendedSeedDumpItem<'a>
{
    fn from(
        value: (
            &'a InMemoryReference,
            &'a str,
            usize,
            usize,
            &'a ExtendedSeed,
            &'a str,
            &'a str,
        ),
    ) -> Self {
        ExtendedSeedDumpItem {
            reference: value.0,
            read_id: value.1,
            read_len: value.2,
            seed_num: value.3,
            seed: value.4,
            seq: value.5,
            qual: value.6,
            tags: vec![],
        }
    }
}

impl<'a>
    From<(
        &'a InMemoryReference,
        &'a str,
        usize,
        usize,
        &'a ExtendedSeed,
        &'a str,
        &'a str,
        Vec<(String, TagValue)>,
    )> for ExtendedSeedDumpItem<'a>
{
    fn from(
        value: (
            &'a InMemoryReference,
            &'a str,
            usize,
            usize,
            &'a ExtendedSeed,
            &'a str,
            &'a str,
            Vec<(String, TagValue)>,
        ),
    ) -> Self {
        ExtendedSeedDumpItem {
            reference: value.0,
            read_id: value.1,
            read_len: value.2,
            seed_num: value.3,
            seed: value.4,
            seq: value.5,
            qual: value.6,
            tags: value.7,
        }
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
            kmer_uniqueness: 1,
            read_frequency: 1,
            is_reverse,
            weight: OrderedFloat(weight),
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
                kmer_uniqueness: 5,
                read_frequency: 1,
                is_reverse: false,
                weight: OrderedFloat(w_10_5),
            },
            ExtendedSeed {
                read_start: 5,
                length: 10,
                ref_chrom_id: 0,
                ref_start: 105,
                multiplicity: 2,
                kmer_uniqueness: 2,
                read_frequency: 1,
                is_reverse: false,
                weight: OrderedFloat(w_10_2),
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
        // Overlap >= other.length: other contributes nothing new.
        let a = seed(0, 20, 0, 100, false);
        let b = seed(0, 20, 0, 120, false); // 20-base overlap == other.length
        assert!(a.edge_penalty(&b).is_none());

        // Overlap > MAX_READ_OVERLAP (50): too much redundancy.
        let a2 = seed(0, 100, 0, 100, false);
        let b2 = seed(40, 60, 0, 200, false); // 60-base overlap > 50
        assert!(a2.edge_penalty(&b2).is_none());
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
        let (p, _) = a.edge_penalty(&b).unwrap();
        assert!(p < 1.0, "expected small penalty, got {p}");
    }

    #[test]
    fn penalty_far_apart_is_large() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(1000, 10, 0, 1100, false);
        let (p, _) = a.edge_penalty(&b).unwrap();
        assert!(p > 900.0, "expected large penalty, got {p}");
    }

    #[test]
    fn penalty_small_deletion_is_cheap() {
        // 20bp deletion: deviation = 20, below gap_linear_threshold (50).
        // penalty = ln(1 + 20) ≈ 3.04 — purely logarithmic, no linear term.
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 130, false);
        let (p, _) = a.edge_penalty(&b).unwrap();
        assert!(
            p < 5.0,
            "expected cheap penalty for small deletion, got {p}"
        );
    }

    #[test]
    fn penalty_large_deletion_grows_with_scale() {
        // 50000bp deletion: well above gap_linear_threshold.
        // With default scale=0.08 and threshold=50 the linear term dominates.
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 50_110, false);
        let (p, _) = a.edge_penalty(&b).unwrap();
        let cfg = config::get();
        let expected_min =
            cfg.seeding.gap_linear_scale * (50_000.0 - cfg.seeding.gap_linear_threshold);
        assert!(
            p > expected_min,
            "expected large penalty for big deletion, got {p}"
        );
    }

    #[test]
    fn penalty_different_chrom_uses_sv_penalty() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 1, 100, false);
        let (p, _) = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_different_strand_uses_sv_penalty() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 100, true);
        let (p, _) = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_non_colinear_fwd_uses_sv_penalty() {
        // Forward strand but ref goes backwards.
        let a = seed(0, 10, 0, 200, false);
        let b = seed(10, 10, 0, 100, false);
        let (p, _) = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_colinear_reverse_strand() {
        // Reverse strand: ref positions decrease as read advances.
        let a = seed(0, 10, 0, 200, true);
        let b = seed(10, 10, 0, 180, true);
        // ref gap = 200 - (180 + 10) = 10, read gap = 0.
        let (p, _) = a.edge_penalty(&b).unwrap();
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
        let (p, _) = a.edge_penalty(&b).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    // ── edge_penalty with small ref overlaps ─────────────────────────

    #[test]
    fn penalty_reverse_1bp_ref_overlap_is_small() {
        // 1-base ref overlap → deviation is 1 → penalty = ln(2) ≈ 0.69.
        let a = seed(0, 100, 0, 500, true);
        let b = seed(100, 100, 0, 401, true);
        let (p, _) = a.edge_penalty(&b).unwrap();
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
        let (p, _) = a.edge_penalty(&b).unwrap();
        let expected = (2.0f64).ln();
        assert!(
            (p - expected).abs() < 0.01,
            "expected ~{expected:.2}, got {p}"
        );
    }

    // ── prune_repetitive_seeds ──────────────────────────────────────────

    #[test]
    fn prune_removes_seed_with_opposite_diagonal_shifts() {
        // Seed at read 100 is on diagonal 0 (ref 100 - read 100 = 0).
        // Its neighbors are on diagonals +50 and -50 — opposite shifts of 50,
        // both exceeding the threshold of 10.  The middle seed should be pruned.
        let mut seeds = vec![
            seed(0, 10, 0, 50, false),    // diagonal = 50
            seed(100, 10, 0, 100, false), // diagonal = 0  ← repetitive
            seed(200, 10, 0, 150, false), // diagonal = -50
        ];
        let n = seeds.len();
        ExtendedSeed::prune_repetitive_seeds(&mut seeds, &mut vec![false; n - 1], 10, None);
        assert_eq!(seeds.len(), 2);
        assert_eq!(seeds[0].ref_start, 50);
        assert_eq!(seeds[1].ref_start, 150);
    }

    #[test]
    fn prune_keeps_colinear_seeds() {
        // All three seeds on the same diagonal (20) — no pruning.
        let mut seeds = vec![
            seed(0, 10, 0, 20, false),    // diagonal = 20
            seed(100, 10, 0, 120, false), // diagonal = 20
            seed(200, 10, 0, 220, false), // diagonal = 20
        ];
        let n = seeds.len();
        ExtendedSeed::prune_repetitive_seeds(&mut seeds, &mut vec![false; n - 1], 10, None);
        assert_eq!(seeds.len(), 3, "colinear seeds should not be pruned");
    }

    #[test]
    fn prune_keeps_seed_when_shifts_below_threshold() {
        // Shifts of 5 on both sides — both within the threshold of 10.
        let mut seeds = vec![
            seed(0, 10, 0, 5, false),     // diagonal = 5
            seed(100, 10, 0, 100, false), // diagonal = 0
            seed(200, 10, 0, 195, false), // diagonal = -5
        ];
        let n = seeds.len();
        ExtendedSeed::prune_repetitive_seeds(&mut seeds, &mut vec![false; n - 1], 10, None);
        assert_eq!(seeds.len(), 3, "small shifts should not be pruned");
    }

    #[test]
    fn prune_ignores_neighbors_on_different_chrom_or_strand() {
        // Three segments separated by SV breaks: chrom 1 / chrom 0 / chrom 1.
        // No seed should be pruned — each segment has only one seed.
        let mut seeds = vec![
            seed(0, 10, 1, 50, false),    // chrom 1
            seed(100, 10, 0, 100, false), // chrom 0
            seed(200, 10, 1, 150, false), // chrom 1
        ];
        let mut sv_breaks = vec![true, true]; // SV break on both sides of middle seed
        ExtendedSeed::prune_repetitive_seeds(&mut seeds, &mut sv_breaks, 10, None);
        assert_eq!(seeds.len(), 3);
    }

    #[test]
    fn prune_too_few_seeds() {
        let mut one = vec![seed(0, 10, 0, 100, false)];
        ExtendedSeed::prune_repetitive_seeds(&mut one, &mut vec![], 10, None);
        assert_eq!(one.len(), 1);

        let mut two = vec![seed(0, 10, 0, 100, false), seed(100, 10, 0, 200, false)];
        ExtendedSeed::prune_repetitive_seeds(&mut two, &mut vec![false], 10, None);
        assert_eq!(two.len(), 2);
    }
}
