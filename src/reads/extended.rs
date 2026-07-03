use std::sync::OnceLock;

use noodles::sam::alignment::{
    record::Flags,
    record::cigar::{Op, op::Kind},
    record_buf::{Cigar, Data, RecordBuf, data::field::Value},
    record::data::field::Tag,
};
use ordered_float::OrderedFloat;
use parallax::utils::telemetry::histogram::HistogramRecorder;

use crate::align::{Alignment, DpAligner};
use crate::reads::{builder::build_record, seeds::SeedHit};
use parallax::{
    config::SeedingConfig,
    reference::InMemoryReference,
    utils::{sequence::complement, telemetry::RecorderExt},
};

/// The type of transition between two consecutive seeds in a chain.
///
/// Stored in the `edge_types` vector parallel to the seed vector: `edge_types[i]`
/// describes the edge from `seeds[i]` to `seeds[i+1]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeType {
    /// Seeds are colinear on the same chromosome and strand with a consistent
    /// diagonal — the transition is a simple insertion, deletion, or exact match.
    Continuation,
    /// Seeds cross a structural variant boundary: different chromosome, different
    /// strand, or a large non-colinear reference jump.  The chain segments on
    /// either side of this edge are reported as separate alignment records.
    SvBreak,
    /// Seeds step backward in reference space but remain within a narrow ref
    /// window, indicating the read is traversing extra copies of a tandem repeat
    /// that are absent from the reference.  Runs of consecutive Repeat edges
    /// (flanked by Continuation segments) are candidates for collapsing into a
    /// single insertion event.
    Repeat,
}

impl EdgeType {
    /// Whether this edge represents any kind of discontinuity (not a simple colinear gap).
    pub fn is_break(self) -> bool {
        self != EdgeType::Continuation
    }

    /// Return the more significant of two edge types: Sv > Repeat > Continuation.
    /// Used when merging edges across removed seeds.
    fn max(self, other: EdgeType) -> EdgeType {
        match (self, other) {
            (EdgeType::SvBreak, _) | (_, EdgeType::SvBreak) => EdgeType::SvBreak,
            (EdgeType::Repeat, _) | (_, EdgeType::Repeat) => EdgeType::Repeat,
            _ => EdgeType::Continuation,
        }
    }
}

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
    kmer_frequency: u32,
    read_frequency: u32,
    is_reverse: bool,
    weight: OrderedFloat<f64>,
}

impl ExtendedSeed {
    pub fn from_seed_hit(seed: &SeedHit, read_len: usize, is_reverse: bool) -> Self {
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
            kmer_frequency: seed.kmer_uniqueness,
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

    pub fn kmer_uniqueness(&self) -> u32 {
        self.kmer_frequency
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
        assert!(
            self.length > 0,
            "to_alignment called on zero-length seed (read_start={})",
            self.read_start
        );
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
            let curr_kmer_uniqueness = seeds[read].kmer_frequency;
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
                let new_kmer_uniqueness = prev.kmer_frequency.min(curr_kmer_uniqueness);
                let new_read_frequency = prev.read_frequency.max(curr_read_frequency);
                let new_weight =
                    Self::calculate_weight(new_length as f64, new_kmer_uniqueness as f64);
                seeds[write].length = new_length;
                seeds[write].kmer_frequency = new_kmer_uniqueness;
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

    /// Prune seeds identified as diagonal excursions, iterating until stable.
    pub fn prune_repetitive_seeds(
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        threshold: isize,
        seeding_cfg: &SeedingConfig,
    ) {
        let hf_filter = HighReadFrequencyFilter {
            threshold,
            min_weight_per_frequency: 50.0,
        };
        let excursion_filter = ImmediateDiagonalExcursionFilter { threshold };
        let terminal_filter = TerminalDiagonalExcursionFilter { threshold };
        if let Err(e) = Self::validate_chain(seeds, sv_breaks) {
            log::error!("chain invalid on entry to prune_repetitive_seeds: {e}");
        }
        hf_filter.apply_until_stable(seeds, sv_breaks, seeding_cfg);
        excursion_filter.apply_until_stable(seeds, sv_breaks, seeding_cfg);
        terminal_filter.apply_filter(seeds, sv_breaks, seeding_cfg);
    }

    /// Remove flagged seeds and update `sv_breaks` in sync.
    /// Breaks spanning removed seeds are OR'd together into a single break
    /// between the neighbouring kept seeds.
    fn remove_flagged(
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        flagged: &[bool],
    ) {
        let n = seeds.len();
        let mut write = 0;
        let mut pending_break = EdgeType::Continuation;
        for read in 0..n {
            // Accumulate the most significant edge type seen so far in the removed span.
            if read > 0 {
                pending_break = pending_break.max(sv_breaks[read - 1]);
            }
            if !flagged[read] {
                if write > 0 {
                    // Emit the merged edge between the previous kept seed and this one.
                    sv_breaks[write - 1] = pending_break;
                }
                pending_break = EdgeType::Continuation;
                seeds.swap(write, read);
                write += 1;
            }
        }
        seeds.truncate(write);
        sv_breaks.truncate(write.saturating_sub(1));
    }

    /// Validate the internal consistency of a seed chain.
    ///
    /// For each colinear segment (run of seeds between sv_breaks), checks:
    /// - All seeds share the same chromosome and strand.
    /// - Forward strand: ref_start is strictly increasing between seeds.
    /// - Reverse strand: ref_start is strictly decreasing between seeds.
    /// - No adjacent seeds overlap on the read.
    ///
    /// Returns `Ok(())` if the chain is valid, or `Err(String)` describing the
    /// first violation found.
    pub fn validate_chain(seeds: &[ExtendedSeed], sv_breaks: &[EdgeType]) -> Result<(), String> {
        let n = seeds.len();
        if n == 0 {
            return Ok(());
        }
        if sv_breaks.len() != n - 1 {
            return Err(format!(
                "sv_breaks length {} != seeds.len() - 1 = {}",
                sv_breaks.len(),
                n - 1
            ));
        }
        for i in 0..n - 1 {
            let a = &seeds[i];
            let b = &seeds[i + 1];
            if b.read_start < a.read_start {
                return Err(format!(
                    "seeds[{}] and seeds[{}] are out of read order: read_start {} > {}",
                    i,
                    i + 1,
                    a.read_start,
                    b.read_start,
                ));
            }
            let a_read_end = a.read_start + a.length;
            if b.read_start < a_read_end {
                return Err(format!(
                    "seeds[{}] and seeds[{}] overlap on read: [{},{}) and [{},{})",
                    i,
                    i + 1,
                    a.read_start,
                    a_read_end,
                    b.read_start,
                    b.read_start + b.length,
                ));
            }
            if sv_breaks[i].is_break() {
                continue;
            }
            if a.ref_chrom_id != b.ref_chrom_id {
                return Err(format!(
                    "seeds[{}] and seeds[{}] are on different chroms ({} vs {}) but no sv_break",
                    i,
                    i + 1,
                    a.ref_chrom_id,
                    b.ref_chrom_id
                ));
            }
            if a.is_reverse != b.is_reverse {
                return Err(format!(
                    "seeds[{}] and seeds[{}] are on different strands but no sv_break",
                    i,
                    i + 1
                ));
            }
            if a.is_reverse {
                if b.ref_start >= a.ref_start {
                    return Err(format!(
                        "reverse seeds[{}] and seeds[{}] not anticolinear: ref_start {} >= {}",
                        i,
                        i + 1,
                        b.ref_start,
                        a.ref_start
                    ));
                }
            } else if b.ref_start <= a.ref_start {
                return Err(format!(
                    "forward seeds[{}] and seeds[{}] not colinear: ref_start {} <= {}",
                    i,
                    i + 1,
                    b.ref_start,
                    a.ref_start
                ));
            }
        }
        Ok(())
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
    pub fn edge_penalty(
        &self,
        other: &ExtendedSeed,
        cfg: &parallax::config::SeedingConfig,
    ) -> Option<(f64, EdgeType)> {
        // ── Read-gap check ────────────────────────────────────────────────────
        // Seeds that overlap heavily on the read cannot form a valid colinear
        // chain.  A small tolerance accommodates minor seed-extension
        // inaccuracies; beyond that the pair is rejected.
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
        let read_gap: f64 = if read_overlap == 0 {
            (other.read_start - self_read_end) as f64
        } else {
            -(read_overlap as f64)
        };

        // ── Read-gap cost ─────────────────────────────────────────────────────
        // A positive read gap means bases in the read are unanchored between
        // the two seeds.  The quadratic term makes long unanchored stretches
        // disproportionately expensive: a large insertion is better represented
        // as an SV breakpoint (paying sv_penalty once) than as a colinear gap
        // that leaves many read bases unexplained.
        let read_gap_cost = read_gap.max(0.0);
        let read_gap_cost = read_gap_cost + cfg.read_gap_quad_scale * read_gap_cost * read_gap_cost;

        // ── Reference penalty ─────────────────────────────────────────────────
        // We classify the transition into one of three cases and assign a
        // reference-space penalty accordingly.

        // Case 1 — cross-chromosome or cross-strand: always a hard SV break.
        // Case 2 — same chrom/strand, colinear (ref_gap >= 0): penalty scales
        //   with the deviation between ref gap and read gap (insertion/deletion
        //   size).  A logarithmic base keeps small gaps cheap; a linear term
        //   above `gap_linear_threshold` suppresses implausibly large jumps
        //   (e.g. two seeds coincidentally colinear but millions of bp apart).
        //   If deviation exceeds MAX_GAP_DEVIATION we treat it as an SV
        //   (without the cap, ln(1 + 14M) ≈ 16.5 < sv_penalty, so distant
        //   same-strand pairs would spuriously out-compete genuine SV edges).
        // Case 3 — same chrom/strand, backward ref jump (ref_gap < 0): normally
        //   an SV break, but if both seeds land within a narrow ref window
        //   (`repeat_expansion_max_ref_window`) the backward step is most
        //   likely a tandem repeat traversal — the read contains extra copies
        //   of a short motif all anchored to the same ref region.  We apply
        //   the same log+linear formula as a forward gap of the same size
        //   (Case 2), plus a fixed `repeat_expansion_penalty` additive so that
        //   a backward step always costs more than an equivalent forward gap.
        //   This lets expansion seeds chain through without accumulating a
        //   prohibitive cost, while still preferring colinear chains.
        const MAX_GAP_DEVIATION: f64 = 100_000.0;
        const REF_OVERLAP_TOLERANCE: i64 = 10;

        let (ref_penalty, edge_type) =
            if self.ref_chrom_id != other.ref_chrom_id || self.is_reverse != other.is_reverse {
                // Case 1: cross-chromosome or cross-strand.
                (cfg.sv_penalty, EdgeType::SvBreak)
            } else {
                let ref_gap = if self.is_reverse {
                    let other_ref_end = (other.ref_start + other.length) as i64;
                    self.ref_start as i64 - other_ref_end
                } else {
                    let self_ref_end = (self.ref_start + self.length) as i64;
                    other.ref_start as i64 - self_ref_end
                };

                if ref_gap >= -REF_OVERLAP_TOLERANCE {
                    // Case 2: colinear (or within tolerance).
                    let deviation = (ref_gap as f64 - read_gap).abs();
                    if deviation > MAX_GAP_DEVIATION {
                        (cfg.sv_penalty, EdgeType::SvBreak)
                    } else {
                        let log_part = (1.0 + deviation.min(cfg.gap_linear_threshold)).ln();
                        let linear_part =
                            cfg.gap_linear_scale * (deviation - cfg.gap_linear_threshold).max(0.0);
                        (log_part + linear_part, EdgeType::Continuation)
                    }
                } else {
                    // Case 3: backward ref jump — SV or tandem repeat traversal.
                    let ref_window = cfg.repeat_expansion_max_ref_window;
                    let ref_distance = self.ref_start.abs_diff(other.ref_start);
                    if ref_window > 0 && ref_distance <= ref_window {
                        // Both seeds are within a narrow ref window: treat as a
                        // tandem repeat traversal.  Use the same log+linear
                        // formula as a forward gap of the same size, plus a
                        // fixed additive penalty so backward steps always cost
                        // more than an equivalent forward gap.
                        let d = (-ref_gap) as f64;
                        let log_part = (1.0 + d.min(cfg.gap_linear_threshold)).ln();
                        let linear_part =
                            cfg.gap_linear_scale * (d - cfg.gap_linear_threshold).max(0.0);
                        (
                            cfg.repeat_expansion_penalty + log_part + linear_part,
                            EdgeType::Repeat,
                        )
                    } else {
                        (cfg.sv_penalty, EdgeType::SvBreak)
                    }
                }
            };

        Some((read_gap_cost + ref_penalty, edge_type))
    }

    /// Collinearity-weighted edge penalty.
    ///
    /// Returns `(penalty, edge_type, next_weight_scale)` where `next_weight_scale`
    /// accounts for any read-overlap truncation of `other`.  Returns `None` if
    /// `other` is fully consumed by the overlap with `self`.
    pub fn edge_penalty_v2(
        &self,
        other: &ExtendedSeed,
        cfg: &parallax::config::SeedingConfig,
    ) -> Option<(f64, EdgeType, f64)> {
        const REF_OVERLAP_TOLERANCE: i64 = 10;
        const REF_DEV_THRESHOLD: f64 = 50.0;
        const REF_DEV_COST_HI: f64 = 1.0;
        const REF_DEV_COST_LO: f64 = 0.1;
        const READ_GAP_THRESHOLD: f64 = 15.0;
        const READ_GAP_COST_LO: f64 = 2.0;
        const READ_GAP_COST_HI: f64 = 10.0;

        let self_read_end = self.read_start + self.length;
        let read_overlap = self_read_end.saturating_sub(other.read_start);

        let (next_weight_scale, read_gap_cost, eff_ref_start, eff_ref_end) = if read_overlap > 0 {
            if read_overlap >= other.length {
                return None;
            }
            let scale = (other.length - read_overlap) as f64 / other.length as f64;
            let (ers, ere) = if other.is_reverse {
                (other.ref_start, other.ref_start + other.length - read_overlap)
            } else {
                (other.ref_start + read_overlap, other.ref_start + other.length)
            };
            (scale, 0.0_f64, ers, ere)
        } else {
            let rg = (other.read_start - self_read_end) as f64;
            let cost = READ_GAP_COST_LO * rg
                + (READ_GAP_COST_HI - READ_GAP_COST_LO) * (rg - READ_GAP_THRESHOLD).max(0.0);
            (1.0_f64, cost, other.ref_start, other.ref_start + other.length)
        };

        if self.ref_chrom_id != other.ref_chrom_id || self.is_reverse != other.is_reverse {
            return Some((read_gap_cost + cfg.sv_penalty, EdgeType::SvBreak, next_weight_scale));
        }

        let eff_read_gap = if read_overlap > 0 { 0i64 } else { (other.read_start - self_read_end) as i64 };

        let ref_gap: i64 = if self.is_reverse {
            self.ref_start as i64 - eff_ref_end as i64
        } else {
            eff_ref_start as i64 - (self.ref_start + self.length) as i64
        };

        if ref_gap >= -REF_OVERLAP_TOLERANCE {
            let deviation = (ref_gap - eff_read_gap).unsigned_abs() as f64;
            if deviation > cfg.collinearity_max_gap_deviation {
                return Some((read_gap_cost + cfg.sv_penalty, EdgeType::SvBreak, next_weight_scale));
            }
            let ref_penalty = REF_DEV_COST_HI * deviation
                + (REF_DEV_COST_LO - REF_DEV_COST_HI) * (deviation - REF_DEV_THRESHOLD).max(0.0);
            Some((read_gap_cost + ref_penalty, EdgeType::Continuation, next_weight_scale))
        } else {
            let ref_window = cfg.repeat_expansion_max_ref_window;
            let deviation = (ref_gap - eff_read_gap).unsigned_abs() as f64;
            if ref_window > 0 && deviation <= ref_window as f64 {
                let ref_penalty = REF_DEV_COST_HI * deviation
                    + (REF_DEV_COST_LO - REF_DEV_COST_HI) * (deviation - REF_DEV_THRESHOLD).max(0.0);
                Some((read_gap_cost + cfg.repeat_expansion_penalty + ref_penalty, EdgeType::Repeat, next_weight_scale))
            } else {
                Some((read_gap_cost + cfg.sv_penalty, EdgeType::SvBreak, next_weight_scale))
            }
        }
    }

    /// Compute per-seed collinearity weights.
    ///
    /// For each seed, sums `1 / (1 + d²)` over all seeds on the same chrom/strand
    /// whose diagonal is within `diagonal_cutoff` bp.  Self-contribution is 1.0;
    /// isolated seeds (no neighbours within the window) get weight ≈ 1.0.
    pub fn compute_collinearity_weights(seeds: &[ExtendedSeed], diagonal_cutoff: f64) -> Vec<f64> {
        let n = seeds.len();
        let mut indexed: Vec<(usize, isize)> = (0..n).map(|i| (i, seeds[i].diagonal())).collect();
        indexed.sort_unstable_by_key(|&(i, d)| (seeds[i].ref_chrom_id, seeds[i].is_reverse, d));

        let mut weights = vec![0.0f64; n];
        let cutoff = diagonal_cutoff as isize;
        let mut lo = 0usize;

        for hi in 0..n {
            let (i, d_hi) = indexed[hi];
            // Advance lo to keep window within cutoff.
            while {
                let (j, d_lo) = indexed[lo];
                seeds[j].ref_chrom_id != seeds[i].ref_chrom_id
                    || seeds[j].is_reverse != seeds[i].is_reverse
                    || d_hi - d_lo > cutoff
            } {
                lo += 1;
            }
            for k in lo..=hi {
                let (j, d_j) = indexed[k];
                let d = (d_hi - d_j) as f64;
                let w = 1.0 / (1.0 + d * d);
                weights[i] += w;
                if k != hi {
                    weights[j] += w; // symmetric contribution
                }
            }
        }
        weights
    }

    /// Collinearity-based seed weight: `length * collinearity / sqrt(kmer_frequency)`.
    fn collinearity_seed_weight(&self, collinearity: f64) -> f64 {
        let freq = self.kmer_frequency.max(1) as f64;
        self.length as f64 * collinearity / freq.sqrt()
    }

    /// Return the set of seed indices that have no colinear neighbour within
    /// `max_gap_deviation` bp — seeds that can only ever appear as isolated
    /// single-seed SV-break segments and are safe to prune before the DP.
    pub fn find_isolated_seeds(seeds: &[ExtendedSeed], cfg: &parallax::config::SeedingConfig) -> Vec<bool> {
        let n = seeds.len();
        let dev = cfg.collinearity_max_gap_deviation as i64;
        let tol = REF_OVERLAP_TOLERANCE_STATIC;
        let mut has_neighbour = vec![false; n];

        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by_key(|&i| seeds[i].read_start);

        for rank in 0..n {
            let i = order[rank];
            if has_neighbour[i] {
                continue;
            }
            let s = &seeds[i];
            let s_read_end = s.read_start + s.length;

            for rank2 in rank + 1..n {
                let j = order[rank2];
                let t = &seeds[j];
                let read_gap = t.read_start as i64 - s_read_end as i64;
                if read_gap > dev {
                    break;
                }
                if t.ref_chrom_id != s.ref_chrom_id || t.is_reverse != s.is_reverse {
                    continue;
                }
                let overlap = s_read_end.saturating_sub(t.read_start);
                let (eff_ref_start, eff_ref_end) = if t.is_reverse {
                    (t.ref_start, t.ref_start + t.length - overlap)
                } else {
                    (t.ref_start + overlap, t.ref_start + t.length)
                };
                let ref_gap: i64 = if s.is_reverse {
                    s.ref_start as i64 - eff_ref_end as i64
                } else {
                    eff_ref_start as i64 - (s.ref_start + s.length) as i64
                };
                if ref_gap < -tol {
                    continue;
                }
                let eff_read_gap = read_gap.max(0);
                if (ref_gap - eff_read_gap).abs() <= dev {
                    has_neighbour[i] = true;
                    has_neighbour[j] = true;
                    break;
                }
            }
        }

        // Return true for seeds that are isolated (no neighbour found).
        has_neighbour.iter().map(|&h| !h).collect()
    }
}

const REF_OVERLAP_TOLERANCE_STATIC: i64 = 10;

// ── Compact DP node representation ───────────────────────────────────────────
//
// During `form_explanatory_groups` we operate on a reduced set of nodes.
// Adjacent seeds (in read-start order) that are very close on both the read
// and the reference are merged eagerly before the DP, reducing n and therefore
// the O(n²) work.  A merged node carries the pre-summed weight (constituent
// weights minus intra-merge edge penalties) so the DP score is identical to
// what a full-n DP would compute.

/// A single node in the chaining DP — either one original seed or a group of
/// eagerly merged seeds.
enum DpNode {
    Single(usize),
    Merged(MergedSeed),
}

struct MergedSeed {
    /// Indices into the original `seeds` slice, in read-start order.
    indices: Vec<usize>,
    /// Pre-summed weight: Σ seed.weight() − Σ edge_penalty for intra-merge edges.
    weight: f64,
}

impl DpNode {
    /// The seed that represents the *left* end of this node (for `read_start`,
    /// chrom, strand, and acting as `other` in `edge_penalty`).
    fn left_seed<'a>(&self, seeds: &'a [ExtendedSeed]) -> &'a ExtendedSeed {
        match self {
            DpNode::Single(i) => &seeds[*i],
            DpNode::Merged(m) => &seeds[m.indices[0]],
        }
    }

    /// The seed that represents the *right* end of this node (for ref-end
    /// computation, acting as `self` in `edge_penalty`).
    fn right_seed<'a>(&self, seeds: &'a [ExtendedSeed]) -> &'a ExtendedSeed {
        match self {
            DpNode::Single(i) => &seeds[*i],
            DpNode::Merged(m) => &seeds[*m.indices.last().unwrap()],
        }
    }

    fn weight(&self, seeds: &[ExtendedSeed]) -> f64 {
        match self {
            DpNode::Single(i) => seeds[*i].weight(),
            DpNode::Merged(m) => m.weight,
        }
    }

    /// Expand this node into its constituent seed indices (in read-start order).
    fn indices<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        match self {
            DpNode::Single(i) => std::slice::from_ref(i).iter().copied(),
            DpNode::Merged(m) => m.indices.iter().copied(),
        }
    }
}

fn dp_node_count_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("dp_node_count"))
}

fn dp_merged_node_size_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("dp_merged_node_size"))
}

/// Build the list of `DpNode`s used by the chaining DP.
///
/// Consecutive seeds in `order` (read-start sorted) that satisfy all of:
///   - same chromosome and strand
///   - read gap ≤ `MAX_EAGER_GAP`
///   - ref gap ≤ `MAX_EAGER_GAP`  (unsigned distance between adjacent ref ends/starts)
///
/// are merged into a single `DpNode::Merged`.  All other seeds become
/// `DpNode::Single`.  The edge penalty for each intra-merge pair is subtracted
/// from the node's weight so that the DP score is unchanged relative to
/// operating on the full seed list.
fn build_dp_nodes(
    seeds: &[ExtendedSeed],
    order: &[usize],
    available: &[bool],
    seeding_cfg: &SeedingConfig,
) -> Vec<DpNode> {
    const MAX_EAGER_GAP: usize = 5;

    let active: Vec<usize> = order.iter().copied().filter(|&i| available[i]).collect();
    if active.is_empty() {
        return Vec::new();
    }
    let orig_len = active.len();

    // Sort by (chrom, strand, read_start) so that all merge candidates within a
    // (chrom, strand) group are contiguous.  Seeds that belong to different
    // groups can never be merged, so interleaving them (as the original
    // read-start order does) caused the greedy scan to break prematurely.
    let mut by_group = active.clone();
    by_group.sort_unstable_by_key(|&i| {
        let s = &seeds[i];
        (s.ref_chrom_id, s.is_reverse, s.read_start)
    });

    let mut consumed = vec![false; by_group.len()];
    let mut nodes: Vec<DpNode> = Vec::with_capacity(active.len());

    for i in 0..by_group.len() {
        if consumed[i] {
            continue;
        }

        let head_idx = by_group[i];
        let head = &seeds[head_idx];
        let mut merged_indices = vec![head_idx];
        let mut merged_weight = head.weight();

        // Scan forward within the same (chrom, strand) group merging seeds onto
        // this node.  read_gap from the head is monotone in j, so take_while
        // gives us a bounded window of candidates; find picks the first one
        // whose ref_gap from the current tail is also within threshold.  Seeds
        // that fail ref_gap are skipped (not consumed) so they can start their
        // own nodes later.
        let mut tail_idx = head_idx;
        let mut j = i + 1;

        loop {
            let tail = &seeds[tail_idx];
            let tail_read_end = tail.read_start + tail.length;
            let found = by_group[j..]
                .iter()
                .copied()
                .enumerate()
                .take_while(|&(_, idx)| {
                    let s = &seeds[idx];
                    s.ref_chrom_id == head.ref_chrom_id
                        && s.is_reverse == head.is_reverse
                        && s.read_start < tail_read_end + MAX_EAGER_GAP
                })
                .find(|&(jj, idx)| {
                    !consumed[j + jj]
                        && seeds[idx].read_start >= tail_read_end
                        && if tail.is_reverse {
                            let next_ref_end = seeds[idx].ref_start + seeds[idx].length;
                            tail.ref_start >= next_ref_end
                                && tail.ref_start - next_ref_end <= MAX_EAGER_GAP
                        } else {
                            let tail_ref_end = tail.ref_start + tail.length;
                            seeds[idx].ref_start >= tail_ref_end
                                && seeds[idx].ref_start - tail_ref_end <= MAX_EAGER_GAP
                        }
                });

            match found {
                None => break,
                Some((jj, next_idx)) => {
                    let next = &seeds[next_idx];
                    let penalty = tail
                        .edge_penalty(next, seeding_cfg)
                        .map(|(p, _)| p)
                        .unwrap_or(0.0);
                    merged_weight += next.weight() - penalty;
                    merged_indices.push(next_idx);
                    consumed[j + jj] = true;
                    tail_idx = next_idx;
                    j += jj + 1;
                }
            }
        }

        let node = if merged_indices.len() == 1 {
            dp_merged_node_size_recorder().record(1usize);
            DpNode::Single(merged_indices[0])
        } else {
            dp_merged_node_size_recorder().record(merged_indices.len());
            // Restore read-start order within the merged node (the DP relies on
            // left_seed / right_seed being the true read-order endpoints).
            merged_indices.sort_unstable_by_key(|&k| seeds[k].read_start);
            DpNode::Merged(MergedSeed {
                indices: merged_indices,
                weight: merged_weight,
            })
        };
        nodes.push(node);
    }

    // Restore read-start order across nodes so the DP processes them in the
    // correct sequence (it assumes nodes are ordered by read position).
    nodes.sort_unstable_by_key(|n| n.left_seed(seeds).read_start);

    dp_node_count_recorder().record(orig_len - nodes.len());
    nodes
}

impl ExtendedSeed {
    /// Form explanatory groups by greedy peeling.
    ///
    /// Each group is a complete alternative explanation of the read — a chain of
    /// non-overlapping seeds (on the read) that maximises total weight minus
    /// edge penalties.  Successive groups are extracted by removing consumed
    /// seeds and repeating the DP until the best remaining chain falls below
    /// `MIN_GROUP_WEIGHT` or `MAX_GROUPS` is reached.
    ///
    /// When `seeding_cfg.use_collinearity_weights` is true, collinearity weights
    /// are computed before the DP, isolated seeds are pruned, and `edge_penalty_v2`
    /// is used instead of `edge_penalty`.
    #[inline(never)]
    pub fn form_explanatory_groups(
        seeds: &[ExtendedSeed],
        seeding_cfg: &SeedingConfig,
    ) -> Vec<(Vec<ExtendedSeed>, Vec<EdgeType>)> {
        if seeding_cfg.use_collinearity_weights {
            Self::form_explanatory_groups_v2(seeds, seeding_cfg)
        } else {
            Self::form_explanatory_groups_v1(seeds, seeding_cfg)
        }
    }

    fn form_explanatory_groups_v1(
        seeds: &[ExtendedSeed],
        seeding_cfg: &SeedingConfig,
    ) -> Vec<(Vec<ExtendedSeed>, Vec<EdgeType>)> {
        const MIN_GROUP_WEIGHT: f64 = 50.0;
        const MAX_GROUPS: usize = 10;
        const MIN_RELATIVE_SCORE: f64 = 0.05;

        let mut groups: Vec<(Vec<ExtendedSeed>, Vec<EdgeType>)> = Vec::new();
        if seeds.is_empty() {
            return groups;
        }

        let mut order: Vec<usize> = (0..seeds.len()).collect();
        order.sort_by_key(|&i| seeds[i].read_start);

        let mut available = vec![true; seeds.len()];
        let mut group0_score = 0.0f64;

        for g in 0..MAX_GROUPS {
            let nodes = build_dp_nodes(seeds, &order, &available, seeding_cfg);
            if nodes.is_empty() {
                break;
            }

            let n = nodes.len();
            let mut dp = vec![0.0f64; n];
            let mut pred = vec![usize::MAX; n];
            let mut pred_edge_type = vec![EdgeType::Continuation; n];

            for i in 0..n {
                let node_i = &nodes[i];
                dp[i] = node_i.weight(seeds);
                let left_i = node_i.left_seed(seeds);
                for j in (0..i).rev() {
                    let right_j = nodes[j].right_seed(seeds);
                    if let Some((penalty, edge_type)) = right_j.edge_penalty(left_i, seeding_cfg) {
                        let score = dp[j] + node_i.weight(seeds) - penalty;
                        if score > dp[i] {
                            dp[i] = score;
                            pred[i] = j;
                            pred_edge_type[i] = edge_type;
                        }
                    }
                }
            }

            let best = (0..n)
                .max_by(|&a, &b| dp[a].partial_cmp(&dp[b]).unwrap_or(std::cmp::Ordering::Equal))
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

            let mut chain: Vec<ExtendedSeed> = Vec::new();
            let mut edge_types: Vec<EdgeType> = Vec::new();
            let mut cur = best;
            loop {
                let node = &nodes[cur];
                let node_seeds: Vec<usize> = node.indices().collect();
                for (k, &seed_idx) in node_seeds.iter().enumerate() {
                    if k > 0 {
                        edge_types.push(EdgeType::Continuation);
                    }
                    chain.push(seeds[seed_idx].clone());
                    available[seed_idx] = false;
                }
                if pred[cur] == usize::MAX {
                    break;
                }
                edge_types.push(pred_edge_type[cur]);
                cur = pred[cur];
            }
            chain.reverse();
            edge_types.reverse();
            groups.push((chain, edge_types));
        }

        groups
    }

    fn form_explanatory_groups_v2(
        seeds: &[ExtendedSeed],
        seeding_cfg: &SeedingConfig,
    ) -> Vec<(Vec<ExtendedSeed>, Vec<EdgeType>)> {
        const MIN_GROUP_WEIGHT: f64 = 50.0;
        const MAX_GROUPS: usize = 10;
        const MIN_RELATIVE_SCORE: f64 = 0.05;

        let mut groups: Vec<(Vec<ExtendedSeed>, Vec<EdgeType>)> = Vec::new();
        if seeds.is_empty() {
            return groups;
        }

        // Pre-compute collinearity weights and seed weights once for all seeds.
        let col_weights = Self::compute_collinearity_weights(seeds, seeding_cfg.collinearity_diagonal_cutoff);
        let seed_weights: Vec<f64> = seeds.iter().enumerate()
            .map(|(i, s)| s.collinearity_seed_weight(col_weights[i]))
            .collect();

        // Prune isolated seeds — those with no colinear neighbour — before the DP.
        let isolated = Self::find_isolated_seeds(seeds, seeding_cfg);
        let active_seeds: Vec<&ExtendedSeed> = seeds.iter().enumerate()
            .filter(|&(i, _)| !isolated[i])
            .map(|(_, s)| s)
            .collect();
        let active_weights: Vec<f64> = seeds.iter().enumerate()
            .filter(|&(i, _)| !isolated[i])
            .map(|(i, _)| seed_weights[i])
            .collect();

        if active_seeds.is_empty() {
            return groups;
        }

        let mut order: Vec<usize> = (0..active_seeds.len()).collect();
        order.sort_by_key(|&i| active_seeds[i].read_start);

        let mut available = vec![true; active_seeds.len()];
        let mut group0_score = 0.0f64;

        for g in 0..MAX_GROUPS {
            let active: Vec<usize> = order.iter().copied().filter(|&i| available[i]).collect();
            if active.is_empty() {
                break;
            }

            let n = active.len();
            let mut dp: Vec<f64> = active.iter().map(|&i| active_weights[i]).collect();
            let mut pred = vec![usize::MAX; n];
            let mut pred_edge_type = vec![EdgeType::Continuation; n];

            for i in 0..n {
                let si = active_seeds[active[i]];
                let wi = active_weights[active[i]];
                for j in (0..i).rev() {
                    let sj = active_seeds[active[j]];
                    if let Some((penalty, edge_type, w_scale)) = sj.edge_penalty_v2(si, seeding_cfg) {
                        let score = dp[j] + wi * w_scale - penalty;
                        if score > dp[i] {
                            dp[i] = score;
                            pred[i] = j;
                            pred_edge_type[i] = edge_type;
                        }
                    }
                }
            }

            let best = (0..n)
                .max_by(|&a, &b| dp[a].partial_cmp(&dp[b]).unwrap_or(std::cmp::Ordering::Equal))
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
            log::debug!("Group {g} (v2): best score = {:.2}", best_score);

            let mut chain: Vec<ExtendedSeed> = Vec::new();
            let mut edge_types: Vec<EdgeType> = Vec::new();
            let mut cur = best;
            loop {
                chain.push((*active_seeds[active[cur]]).clone());
                available[active[cur]] = false;
                if pred[cur] == usize::MAX {
                    break;
                }
                edge_types.push(pred_edge_type[cur]);
                cur = pred[cur];
            }
            chain.reverse();
            edge_types.reverse();
            groups.push((chain, edge_types));
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
    pub(crate) fn trim_right(&mut self, n: usize) {
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

    /// Trim `n` bases from the left end of `seeds[j]`, propagating any spill
    /// forward through subsequent seeds to preserve read order.
    ///
    /// If trimming `seeds[j]` by `n` would push its `read_start` past
    /// `seeds[j+1].read_start`, the excess is applied to `seeds[j+1]`, and so
    /// on until the constraint is satisfied or seeds are zeroed out.  Zero-length
    /// seeds must be culled by the caller via `retain_nonzero`.
    fn trim_left_propagate(seeds: &mut Vec<ExtendedSeed>, j: usize, n: usize) {
        seeds[j].trim_left(n);
        let mut k = j;
        while k + 1 < seeds.len() {
            let cur_end = seeds[k].read_start + seeds[k].length;
            let next_start = seeds[k + 1].read_start;
            if cur_end <= next_start {
                break;
            }
            let spill = (cur_end - next_start).min(seeds[k + 1].length);
            seeds[k + 1].trim_left(spill);
            k += 1;
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
        read_name: &str,
        group: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        read_seq: &[u8],
        reference: &InMemoryReference,
        seeding_cfg: &SeedingConfig,
    ) {
        let _ = read_name;

        if group.len() <= 1 {
            return;
        }

        let retain_nonzero = |group: &mut Vec<ExtendedSeed>, sv_breaks: &mut Vec<EdgeType>| {
            let zero_flagged: Vec<bool> = group.iter().map(|s| s.length == 0).collect();
            if zero_flagged.iter().any(|&f| f) {
                Self::remove_flagged(group, sv_breaks, &zero_flagged);
            }
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
            } else if !sv_breaks[i].is_break() {
                // Overlap: trim one seed to remove it.  Only meaningful for
                // colinear pairs — SV-break pairs belong to separate alignment
                // segments, so trimming one based on a read overlap would be
                // incorrect.
                let overlap = a_read_end - b_read_start;

                match (group[i].is_reverse, group[i + 1].is_reverse) {
                    (true, true) => {
                        // Both reverse: trim A's right end on read.
                        // trim_right caps at length internally, so A may become
                        // zero-length; retain_nonzero will remove it.
                        group[i].trim_right(overlap);
                    }
                    (false, false) | (true, false) | (false, true) => {
                        // Trim B's left end on read, propagating any spill forward.
                        Self::trim_left_propagate(group, i + 1, overlap);
                    }
                }
            }
        }

        retain_nonzero(group, sv_breaks);

        Self::recompute_sv_breaks(group, sv_breaks, seeding_cfg);
        Self::resolve_ref_overlaps(group, sv_breaks);
        Self::resolve_read_overlaps(group, sv_breaks);
    }

    /// Resolve reference overlaps between adjacent colinear seeds by trimming.
    ///
    /// Adjacent same-chrom same-strand seeds that are not separated by an SV
    /// break must not share reference bases.  This trims the overlapping end
    /// (right end of the earlier seed for forward strand, left end of the later
    /// seed for reverse strand), then removes any zero-length seeds produced.
    /// Repeated until stable, since one trim can expose a new overlap.
    pub fn resolve_ref_overlaps(seeds: &mut Vec<ExtendedSeed>, sv_breaks: &mut Vec<EdgeType>) {
        loop {
            let mut any_trimmed = false;
            for i in 0..seeds.len().saturating_sub(1) {
                if sv_breaks[i].is_break() {
                    continue;
                }
                if seeds[i].ref_chrom_id != seeds[i + 1].ref_chrom_id {
                    continue;
                }
                if seeds[i].is_reverse != seeds[i + 1].is_reverse {
                    continue;
                }
                if seeds[i].is_reverse {
                    // Reverse: seeds[i].ref_start > seeds[i+1].ref_start.
                    // The right ref end of seeds[i+1] is ref_start_{i+1} + length_{i+1}.
                    // It must not exceed seeds[i].ref_start.
                    // Trim by reducing seeds[i+1]'s right ref end (trim_left on read),
                    // propagating any spill forward to preserve read order.
                    let b_ref_end = seeds[i + 1].ref_start + seeds[i + 1].length;
                    if b_ref_end > seeds[i].ref_start {
                        let overlap = b_ref_end - seeds[i].ref_start;
                        Self::trim_left_propagate(seeds, i + 1, overlap);
                        any_trimmed = true;
                    }
                } else {
                    // Forward: seeds[i].ref_start < seeds[i+1].ref_start.
                    // The right ref end of seeds[i] must not exceed seeds[i+1].ref_start.
                    // Trim by reducing seeds[i]'s right ref end (trim_right).
                    let a_ref_end = seeds[i].ref_start + seeds[i].length;
                    if a_ref_end > seeds[i + 1].ref_start {
                        let overlap = a_ref_end - seeds[i + 1].ref_start;
                        seeds[i].trim_right(overlap);
                        any_trimmed = true;
                    }
                }
            }
            if !any_trimmed {
                break;
            }
            let zero_flagged: Vec<bool> = seeds.iter().map(|s| s.length == 0).collect();
            if zero_flagged.iter().any(|&f| f) {
                Self::remove_flagged(seeds, sv_breaks, &zero_flagged);
            }
        }
    }

    /// Resolve read overlaps between adjacent colinear seeds by trimming.
    ///
    /// After extend-and-trim and ref-overlap resolution, colinear seeds may
    /// still overlap on the read (on different diagonals).  This trims the
    /// overlapping end using the same strand-aware rule as `extend_and_trim`:
    /// for both-reverse pairs trim A's right end; otherwise propagate a
    /// trim_left on B forward.  Zero-length seeds are removed; repeated until
    /// stable.
    pub fn resolve_read_overlaps(seeds: &mut Vec<ExtendedSeed>, sv_breaks: &mut Vec<EdgeType>) {
        loop {
            let mut any_trimmed = false;
            for i in 0..seeds.len().saturating_sub(1) {
                let a_read_end = seeds[i].read_start + seeds[i].length;
                let b_read_start = seeds[i + 1].read_start;
                if a_read_end <= b_read_start {
                    continue;
                }
                let overlap = a_read_end - b_read_start;
                if sv_breaks[i].is_break() {
                    // At an SV break the overlap is microhomology: trim the
                    // right end of seed[i] so seed[i+1]'s locus is unchanged.
                    seeds[i].trim_right(overlap);
                } else if seeds[i].is_reverse && seeds[i + 1].is_reverse {
                    seeds[i].trim_right(overlap);
                } else {
                    Self::trim_left_propagate(seeds, i + 1, overlap);
                }
                any_trimmed = true;
            }
            if !any_trimmed {
                break;
            }
            let zero_flagged: Vec<bool> = seeds.iter().map(|s| s.length == 0).collect();
            if zero_flagged.iter().any(|&f| f) {
                Self::remove_flagged(seeds, sv_breaks, &zero_flagged);
            }
        }
    }

    /// Recompute every SV break from scratch using `edge_penalty`.
    ///
    /// Called after any operation that may change seed positions or remove seeds,
    /// so that breaks that are now simple indels are cleared and new SV-sized gaps
    /// (e.g. exposed by removing a bridging seed) are set.
    pub fn recompute_sv_breaks(
        seeds: &[ExtendedSeed],
        sv_breaks: &mut Vec<EdgeType>,
        seeding_cfg: &SeedingConfig,
    ) {
        for i in 0..sv_breaks.len() {
            match seeds[i].edge_penalty(&seeds[i + 1], seeding_cfg) {
                Some((_, edge_type)) => sv_breaks[i] = edge_type,
                None => sv_breaks[i] = EdgeType::SvBreak,
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
    /// Align the gap between two adjacent seeds.
    ///
    /// Extracts the read subsequence between the end of `a` and the start of
    /// `b`, and the corresponding reference interval, then runs DP alignment.
    /// For reverse-strand seeds the reference is reverse-complemented so the
    /// alignment is always in read-forward orientation.
    ///
    /// Returns `None` if the aligner fails.  If the read or reference interval
    /// is empty the aligner handles it directly (producing a pure deletion or
    /// insertion CIGAR).
    pub fn align_gap(
        a: &ExtendedSeed,
        b: &ExtendedSeed,
        read_seq: &[u8],
        reference: &InMemoryReference,
        aligner: &mut DpAligner,
    ) -> Option<Alignment> {
        let a_end = a.read_start + a.length;
        if a_end > b.read_start {
            log::warn!(
                "align_gap: seeds overlap on read ([{},{}) vs [{},{}))",
                a.read_start,
                a_end,
                b.read_start,
                b.read_start + b.length
            );
            return None;
        }
        let query = &read_seq[a_end..b.read_start];

        let ref_slice = if a.is_reverse {
            // Reverse strand: ref positions decrease as read advances.
            // Gap on ref is [b.ref_start + b.length .. a.ref_start).
            let ref_begin = b.ref_start + b.length;
            let ref_end = a.ref_start;
            if ref_begin >= ref_end {
                vec![]
            } else {
                let fwd = reference.get_seq(a.ref_chrom_id, ref_begin, ref_end);
                fwd.iter().rev().map(|&base| complement(base)).collect()
            }
        } else {
            // Forward strand: ref positions increase as read advances.
            let ref_begin = a.ref_start + a.length;
            let ref_end = b.ref_start;
            if ref_begin >= ref_end {
                vec![]
            } else {
                reference
                    .get_seq(a.ref_chrom_id, ref_begin, ref_end)
                    .to_vec()
            }
        };

        aligner.align(query, &ref_slice)
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
    pub fn align_gaps(
        read_name: &str,
        group: &[ExtendedSeed],
        sv_breaks: &[EdgeType],
        read_seq: &[u8],
        reference: &InMemoryReference,
        aligner: &mut DpAligner,
    ) -> Vec<Option<Alignment>> {
        if group.len() <= 1 {
            return Vec::new();
        }

        let mut alignments = Vec::with_capacity(group.len() - 1);
        for i in 0..group.len() - 1 {
            if sv_breaks[i].is_break() {
                alignments.push(None);
            } else {
                let a = &group[i];
                let b = &group[i + 1];
                if a.read_start + a.length > b.read_start {
                    log::error!(
                        "{read_name}: align_gaps: colinear seeds[{i}] and seeds[{}] overlap on read: [{},{}) vs [{},{})",
                        i + 1,
                        a.read_start,
                        a.read_start + a.length,
                        b.read_start,
                        b.read_start + b.length,
                    );
                }
                alignments.push(Self::align_gap(a, b, read_seq, reference, aligner));
            }
        }
        alignments
    }
}

pub trait SeedFilter {
    /// Identify which seeds should be removed, returning a mask of length `seeds.len()`.
    fn find_seeds_to_remove(&self, seeds: &[ExtendedSeed], sv_breaks: &[EdgeType]) -> Vec<bool>;

    /// Remove flagged seeds, then recompute sv_breaks.
    /// Returns the number of seeds removed.
    fn apply_filter(
        &self,
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        seeding_cfg: &SeedingConfig,
    ) -> usize {
        let flagged = self.find_seeds_to_remove(seeds, sv_breaks);
        let count = flagged.iter().filter(|&&f| f).count();
        if count > 0 {
            ExtendedSeed::remove_flagged(seeds, sv_breaks, &flagged);
            ExtendedSeed::recompute_sv_breaks(seeds, sv_breaks, seeding_cfg);
            ExtendedSeed::resolve_ref_overlaps(seeds, sv_breaks);
            ExtendedSeed::resolve_read_overlaps(seeds, sv_breaks);
        }
        if let Err(e) = ExtendedSeed::validate_chain(seeds, sv_breaks) {
            log::error!("chain invalid after {}: {e}", std::any::type_name::<Self>());
        }
        count
    }

    /// Repeatedly apply the filter until no seeds are removed.
    fn apply_until_stable(
        &self,
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        seeding_cfg: &SeedingConfig,
    ) {
        while self.apply_filter(seeds, sv_breaks, seeding_cfg) > 0 {}
    }
}

/// Removes seeds whose weight-to-read-frequency ratio is low and that shift
/// the diagonal relative to a neighbour.
pub struct HighReadFrequencyFilter {
    pub threshold: isize,
    pub min_weight_per_frequency: f64,
}

impl SeedFilter for HighReadFrequencyFilter {
    fn find_seeds_to_remove(&self, seeds: &[ExtendedSeed], sv_breaks: &[EdgeType]) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];
        if n < 2 {
            return flagged;
        }
        let threshold = self.threshold;
        for pos in 0..n {
            let has_left_sv = pos == 0 || sv_breaks[pos - 1].is_break();
            let has_right_sv = pos == n - 1 || sv_breaks[pos].is_break();
            if has_left_sv || has_right_sv {
                continue;
            }
            let seed = &seeds[pos];
            if seed.weight() / seed.read_frequency() as f64 >= self.min_weight_per_frequency {
                continue;
            }
            let diag = seed.diagonal();
            let left_shift = if seeds[pos - 1].ref_chrom_id == seed.ref_chrom_id
                && seeds[pos - 1].is_reverse == seed.is_reverse
            {
                Some((seeds[pos - 1].diagonal() - diag).abs())
            } else {
                None
            };
            let right_shift = if seeds[pos + 1].ref_chrom_id == seed.ref_chrom_id
                && seeds[pos + 1].is_reverse == seed.is_reverse
            {
                Some((seeds[pos + 1].diagonal() - diag).abs())
            } else {
                None
            };
            if left_shift.map_or(false, |s| s > threshold)
                || right_shift.map_or(false, |s| s > threshold)
            {
                flagged[pos] = true;
            }
        }
        flagged
    }
}

/// Removes seeds that are immediate diagonal excursions (minimap2 neighbour heuristic).
///
/// For each interior seed not adjacent to an SV break, if the diagonal shift to
/// both neighbours exceeds `threshold` with opposite signs, the seed is flagged.
pub struct ImmediateDiagonalExcursionFilter {
    pub threshold: isize,
}

impl SeedFilter for ImmediateDiagonalExcursionFilter {
    fn find_seeds_to_remove(&self, seeds: &[ExtendedSeed], sv_breaks: &[EdgeType]) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];
        if n < 3 {
            return flagged;
        }
        let threshold = self.threshold;
        for pos in 1..n - 1 {
            if sv_breaks[pos - 1].is_break() || sv_breaks[pos].is_break() {
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
}

/// Removes seeds at the ends of colinear segments whose diagonal deviates from
/// the segment's global weighted median.
pub struct TerminalDiagonalExcursionFilter {
    pub threshold: isize,
}

impl SeedFilter for TerminalDiagonalExcursionFilter {
    fn find_seeds_to_remove(&self, seeds: &[ExtendedSeed], sv_breaks: &[EdgeType]) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];
        if n < 3 {
            return flagged;
        }
        let threshold = self.threshold;

        let weighted_median = |range: std::ops::Range<usize>| -> Option<isize> {
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
                .find(|&i| sv_breaks[i].is_break())
                .map(|i| i + 1)
                .unwrap_or(n);
            if seg_end - seg_start >= 3 {
                if let Some(median) = weighted_median(seg_start..seg_end) {
                    for pos in seg_start..seg_end {
                        if (seeds[pos].diagonal() - median).abs() > threshold {
                            flagged[pos] = true;
                        } else {
                            break;
                        }
                    }
                    for pos in (seg_start..seg_end).rev() {
                        if flagged[pos] {
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
}

/// Removes single-seed segments whose seed length is below `min_length`.
pub struct ShortSingleSeedSegmentFilter {
    pub min_length: usize,
}

impl SeedFilter for ShortSingleSeedSegmentFilter {
    fn find_seeds_to_remove(&self, seeds: &[ExtendedSeed], sv_breaks: &[EdgeType]) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];
        if n == 0 {
            return flagged;
        }
        let mut seg_start = 0;
        loop {
            let seg_end = (seg_start..n - 1)
                .find(|&i| sv_breaks[i].is_break())
                .map(|i| i + 1)
                .unwrap_or(n);
            if seg_end - seg_start == 1 && seeds[seg_start].length < self.min_length {
                flagged[seg_start] = true;
            }
            if seg_end == n {
                break;
            }
            seg_start = seg_end;
        }
        flagged
    }
}

/// Removes single-segment excursions where the flanking segments can rejoin.
///
/// A segment is an excursion candidate when:
///   1. It has at most `max_seeds` seeds.
///   2. Its total reference span is at most `max_ref_span` bp.
///   3. It is flanked on both sides by SV-break edges (i.e. it is a non-colinear
///      detour, not an interior segment of a longer colinear run).
///   4. The seeds immediately flanking it on either side can form a valid edge
///      (any `EdgeType`) via `edge_penalty`, meaning the flanks can reconnect
///      without the excursion.
///
/// Equivalent to Python's `remove_excursions` post-DP filter.
pub struct ExcursionSegmentFilter {
    pub max_seeds: usize,
    pub max_ref_span: usize,
}

impl ExcursionSegmentFilter {
    fn find_excursions(
        &self,
        seeds: &[ExtendedSeed],
        sv_breaks: &[EdgeType],
        cfg: &parallax::config::SeedingConfig,
    ) -> Vec<bool> {
        let n = seeds.len();
        let mut flagged = vec![false; n];
        if n < 3 || self.max_seeds == 0 {
            return flagged;
        }

        let mut segments: Vec<(usize, usize)> = Vec::new();
        let mut seg_start = 0;
        for i in 0..n - 1 {
            if sv_breaks[i].is_break() {
                segments.push((seg_start, i));
                seg_start = i + 1;
            }
        }
        segments.push((seg_start, n - 1));

        let num_segs = segments.len();
        if num_segs < 3 {
            return flagged;
        }

        'outer: for seg_idx in 1..num_segs - 1 {
            let (s_start, s_end) = segments[seg_idx];
            let seg_len = s_end - s_start + 1;

            if seg_len > self.max_seeds {
                continue;
            }

            let ref_lo = seeds[s_start].ref_start.min(seeds[s_end].ref_start);
            let ref_hi = (seeds[s_start].ref_start + seeds[s_start].length)
                .max(seeds[s_end].ref_start + seeds[s_end].length);
            if ref_hi - ref_lo > self.max_ref_span {
                continue;
            }

            let left_break = sv_breaks[s_start - 1];
            let right_break = sv_breaks[s_end];
            if !left_break.is_break() || !right_break.is_break() {
                continue;
            }

            // Walk back past already-flagged segments to find true left flank.
            let mut left_seg_idx = seg_idx - 1;
            while left_seg_idx > 0 {
                let (ls, le) = segments[left_seg_idx];
                if !flagged[ls..=le].iter().all(|&f| f) {
                    break;
                }
                left_seg_idx -= 1;
            }
            let (ls, le) = segments[left_seg_idx];
            if flagged[ls..=le].iter().all(|&f| f) {
                continue 'outer;
            }
            let left_seed = &seeds[le];
            let right_seed = &seeds[segments[seg_idx + 1].0];

            // Criterion 4: flanks can form any valid edge.
            if left_seed.edge_penalty(right_seed, cfg).is_none() {
                continue;
            }

            for k in s_start..=s_end {
                flagged[k] = true;
            }
        }
        flagged
    }

    pub fn apply(&self,
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        cfg: &parallax::config::SeedingConfig,
    ) -> usize {
        let flagged = self.find_excursions(seeds, sv_breaks, cfg);
        let count = flagged.iter().filter(|&&f| f).count();
        if count > 0 {
            ExtendedSeed::remove_flagged(seeds, sv_breaks, &flagged);
            ExtendedSeed::recompute_sv_breaks(seeds, sv_breaks, cfg);
            ExtendedSeed::resolve_ref_overlaps(seeds, sv_breaks);
            ExtendedSeed::resolve_read_overlaps(seeds, sv_breaks);
        }
        if let Err(e) = ExtendedSeed::validate_chain(seeds, sv_breaks) {
            log::error!("chain invalid after ExcursionSegmentFilter: {e}");
        }
        count
    }

    pub fn apply_until_stable(
        &self,
        seeds: &mut Vec<ExtendedSeed>,
        sv_breaks: &mut Vec<EdgeType>,
        cfg: &parallax::config::SeedingConfig,
    ) {
        while self.apply(seeds, sv_breaks, cfg) > 0 {}
    }
}

pub enum TagValue {
    Str(String),
    Int(i64),
    Flt(f64)
}

/// Build a `RecordBuf` for a single seed, suitable for writing via a `RecordWriter`.
///
/// The record uses hard clips for the non-seed read portions and a single `=<len>`
/// CIGAR op for the seed itself.  For reverse-strand seeds the sequence is
/// reverse-complemented (matching what IGV expects for FLAG=0x10 records).
pub fn seed_to_record(
    name: &str,
    read_len: usize,
    seed: &ExtendedSeed,
    seq: &[u8],   // forward-strand query bases for seed's [read_start, read_end)
    qual: &[u8],  // Phred+33 quality for the same range
    tags: Vec<(String, TagValue)>,
) -> RecordBuf {
    let len = seed.length();
    let read_left = seed.read_start();
    let read_right = read_len - seed.read_end();
    let (left_clip, right_clip) = if seed.is_reverse {
        (read_right, read_left)
    } else {
        (read_left, read_right)
    };

    let mut cigar_ops: Vec<Op> = Vec::with_capacity(3);
    if left_clip > 0 {
        cigar_ops.push(Op::new(Kind::HardClip, left_clip));
    }
    cigar_ops.push(Op::new(Kind::SequenceMatch, len));
    if right_clip > 0 {
        cigar_ops.push(Op::new(Kind::HardClip, right_clip));
    }
    let cigar: Cigar = cigar_ops.iter().copied().collect();

    let out_seq: Vec<u8> = if seed.is_reverse {
        seq.iter().rev().map(|&b| complement(b)).collect()
    } else {
        seq.to_vec()
    };
    let out_qual: Vec<u8> = if seed.is_reverse {
        qual.iter().rev().copied().collect()
    } else {
        qual.to_vec()
    };

    let mapq = (seed.weight().floor() as u8).min(254);
    let mut flags = Flags::empty();
    if seed.is_reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }

    let mut data_tags: Vec<(Tag, Value)> = Vec::with_capacity(tags.len());
    for (key, value) in tags {
        let bytes = key.as_bytes();
        if bytes.len() == 2 {
            let tag = Tag::from([bytes[0], bytes[1]]);
            let v = match value {
                TagValue::Str(s) => Value::from(s.as_str()),
                TagValue::Int(i) => Value::from(i as i32),
                TagValue::Flt(f) => Value::from(f as f32),
            };
            data_tags.push((tag, v));
        }
    }
    let data: Data = data_tags.into_iter().collect();

    build_record(
        name,
        flags,
        seed.ref_chrom_id(),
        seed.ref_start() + 1,
        mapq,
        cigar,
        None,
        None,
        &out_seq,
        &out_qual,
        data,
    )
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

impl<'a> parallax::utils::dump::DumpItem for ExtendedSeedDumpItem<'a> {
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
            "{}\t{}\t{}\t{}\t{}\t{}{}={}\t*\t0\t0\t{}\t{}\tXN:i:{}",
            self.read_id,
            flag,
            chrom,
            pos,
            mapq,
            left_clip,
            len,
            right_clip,
            seq,
            self.qual,
            self.seed_num
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
                TagValue::Flt(f) => {
                    write!(writer, "\t{}:f:{}", tag, f).expect("write failed");
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
    use parallax::config;

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
            kmer_frequency: 1,
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
                kmer_frequency: 5,
                read_frequency: 1,
                is_reverse: false,
                weight: OrderedFloat(w_10_5),
            },
            ExtendedSeed {
                read_start: 5,
                length: 10,
                ref_chrom_id: 0,
                ref_start: 105,
                kmer_frequency: 2,
                read_frequency: 1,
                is_reverse: false,
                weight: OrderedFloat(w_10_2),
            },
        ];
        ExtendedSeed::simplify_seeds(&mut seeds);
        assert_eq!(seeds.len(), 1);
        assert_eq!(seeds[0].kmer_frequency, 2);
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
        assert!(a.edge_penalty(&b, &config::get().seeding).is_none());

        // Overlap > MAX_READ_OVERLAP (50): too much redundancy.
        let a2 = seed(0, 100, 0, 100, false);
        let b2 = seed(40, 60, 0, 200, false); // 60-base overlap > 50
        assert!(a2.edge_penalty(&b2, &config::get().seeding).is_none());
    }

    #[test]
    fn penalty_small_read_overlap_is_tolerated() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(9, 10, 0, 110, false); // 1bp overlap
        assert!(a.edge_penalty(&b, &config::get().seeding).is_some());
    }

    #[test]
    fn penalty_adjacent_colinear_is_small() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 110, false);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        assert!(p < 1.0, "expected small penalty, got {p}");
    }

    #[test]
    fn penalty_far_apart_is_large() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(1000, 10, 0, 1100, false);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        assert!(p > 900.0, "expected large penalty, got {p}");
    }

    #[test]
    fn penalty_small_deletion_is_cheap() {
        // 20bp deletion: deviation = 20, below gap_linear_threshold (50).
        // penalty = ln(1 + 20) ≈ 3.04 — purely logarithmic, no linear term.
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 130, false);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
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
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
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
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_different_strand_uses_sv_penalty() {
        let a = seed(0, 10, 0, 100, false);
        let b = seed(10, 10, 0, 100, true);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_non_colinear_fwd_uses_sv_penalty() {
        // Forward strand but ref goes backwards by more than repeat_expansion_max_ref_window.
        let a = seed(0, 10, 0, 1000, false);
        let b = seed(10, 10, 0, 100, false);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    #[test]
    fn penalty_colinear_reverse_strand() {
        // Reverse strand: ref positions decrease as read advances.
        let a = seed(0, 10, 0, 200, true);
        let b = seed(10, 10, 0, 180, true);
        // ref gap = 200 - (180 + 10) = 10, read gap = 0.
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        assert!(
            p < 5.0,
            "expected small penalty for colinear reverse, got {p}"
        );
    }

    #[test]
    fn penalty_non_colinear_reverse_uses_sv_penalty() {
        // Reverse strand but ref increases by more than repeat_expansion_max_ref_window — non-colinear.
        let a = seed(0, 10, 0, 100, true);
        let b = seed(10, 10, 0, 1000, true);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
        let sv = config::get().seeding.sv_penalty;
        assert!((p - sv).abs() < 1e-9, "expected sv_penalty ({sv}), got {p}");
    }

    // ── edge_penalty with small ref overlaps ─────────────────────────

    #[test]
    fn penalty_reverse_1bp_ref_overlap_is_small() {
        // 1-base ref overlap → deviation is 1 → penalty = ln(2) ≈ 0.69.
        let a = seed(0, 100, 0, 500, true);
        let b = seed(100, 100, 0, 401, true);
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
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
        let (p, _) = a.edge_penalty(&b, &config::get().seeding).unwrap();
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
        ExtendedSeed::prune_repetitive_seeds(
            &mut seeds,
            &mut vec![EdgeType::Continuation; n - 1],
            10,
            &config::get().seeding,
        );
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
        ExtendedSeed::prune_repetitive_seeds(
            &mut seeds,
            &mut vec![EdgeType::Continuation; n - 1],
            10,
            &config::get().seeding,
        );
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
        ExtendedSeed::prune_repetitive_seeds(
            &mut seeds,
            &mut vec![EdgeType::Continuation; n - 1],
            10,
            &config::get().seeding,
        );
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
        let mut sv_breaks = vec![EdgeType::SvBreak, EdgeType::SvBreak]; // SV break on both sides of middle seed
        ExtendedSeed::prune_repetitive_seeds(
            &mut seeds,
            &mut sv_breaks,
            10,
            &config::get().seeding,
        );
        assert_eq!(seeds.len(), 3);
    }

    // ── ShortSingleSeedSegmentFilter ────────────────────────────────────

    #[test]
    fn filter_removes_short_single_seed_segment() {
        // Three segments: [seed0] --sv-- [seed1] --sv-- [seed2]
        // seed1 is a single-seed segment with length 10, below the threshold of 20.
        let mut seeds = vec![
            seed(0, 30, 0, 100, false),
            seed(100, 10, 1, 200, false),
            seed(200, 30, 0, 300, false),
        ];
        let mut sv_breaks = vec![EdgeType::SvBreak, EdgeType::SvBreak];
        let removed = ShortSingleSeedSegmentFilter { min_length: 20 }.apply_filter(
            &mut seeds,
            &mut sv_breaks,
            &config::get().seeding,
        );
        assert_eq!(removed, 1);
        assert_eq!(seeds.len(), 2);
        assert_eq!(sv_breaks.len(), 1);
    }

    #[test]
    fn filter_keeps_long_single_seed_segment() {
        let mut seeds = vec![
            seed(0, 30, 0, 100, false),
            seed(100, 25, 1, 200, false),
            seed(200, 30, 0, 300, false),
        ];
        let mut sv_breaks = vec![EdgeType::SvBreak, EdgeType::SvBreak];
        let removed = ShortSingleSeedSegmentFilter { min_length: 20 }.apply_filter(
            &mut seeds,
            &mut sv_breaks,
            &config::get().seeding,
        );
        assert_eq!(removed, 0);
        assert_eq!(seeds.len(), 3);
    }

    #[test]
    fn filter_keeps_multi_seed_short_segment() {
        // A segment with two seeds whose individual lengths are short should not be removed.
        let mut seeds = vec![seed(0, 10, 0, 100, false), seed(10, 10, 0, 110, false)];
        let mut sv_breaks = vec![EdgeType::Continuation]; // colinear, one segment
        let removed = ShortSingleSeedSegmentFilter { min_length: 20 }.apply_filter(
            &mut seeds,
            &mut sv_breaks,
            &config::get().seeding,
        );
        assert_eq!(removed, 0);
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn prune_too_few_seeds() {
        let mut one = vec![seed(0, 10, 0, 100, false)];
        ExtendedSeed::prune_repetitive_seeds(&mut one, &mut vec![], 10, &config::get().seeding);
        assert_eq!(one.len(), 1);

        let mut two = vec![seed(0, 10, 0, 100, false), seed(100, 10, 0, 200, false)];
        ExtendedSeed::prune_repetitive_seeds(
            &mut two,
            &mut vec![EdgeType::Continuation],
            10,
            &config::get().seeding,
        );
        assert_eq!(two.len(), 2);
    }
}
