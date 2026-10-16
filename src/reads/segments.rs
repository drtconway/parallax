use crate::reads::compound::{EdgeType, GapComputable, Seed, Weighted};
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use parallax::config::ChainingConfig;
use parallax::utils::{coverage::SetCoverage, union_find::UnionFind};

// ── Traits ────────────────────────────────────────────────────────────────────

/// Scheme for level-2 segment chaining: partitions compound seeds into
/// connected components and scores edges within each component.
///
/// `reachable` determines whether two seeds (lhs before rhs in read coords)
/// can belong to the same chain — used for the union-find partition and
/// isolated-seed pruning.  `edge_cost` computes the actual penalty for an
/// edge in the chaining DP; returns `None` if the pair cannot be connected.
///
/// The read-gap cutoff and all other thresholds are implementation details
/// of the concrete scheme; callers only see this interface.
pub trait SegmentScheme {
    /// Maximum read-space gap (bases) within which two seeds are candidates
    /// for the same component.  Used to bound forward scans; seeds further
    /// apart than this are never reachable from each other.
    fn max_read_gap(&self) -> i64;

    /// Minimum read-space span (from earliest seed start to latest seed end)
    /// that a component must cover to be worth chaining.  Components whose
    /// span is smaller than this are discarded after partitioning.
    fn min_segment_span(&self) -> i64;

    /// Maximum `|ref_gap - read_gap|` for two seeds to be placed in the same
    /// component.  Seeds whose diagonal shift exceeds this are treated as being
    /// on different loci.
    fn max_ref_deviation(&self) -> i64;

    /// Returns `true` if `rhs` is reachable from `lhs` — i.e. they could
    /// belong to the same segment chain.  Called with `lhs.read_pos() ≤
    /// rhs.read_pos()`; callers need not assume symmetry.
    fn reachable<S: Seed>(&self, lhs: &S, rhs: &S, k: usize) -> bool;

    /// Returns the edge penalty for connecting `lhs` → `rhs` in the chaining
    /// DP, or `None` if the pair cannot be connected.
    fn edge_cost<S: Seed + Weighted + GapComputable>(
        &self,
        lhs: &S,
        rhs: &S,
        k: usize,
    ) -> Option<(f64, EdgeType)>;

    /// Returns the node reward for `seed` in the chaining DP.  The flat scheme
    /// returns `seed.weight()` unchanged; a band-biased scheme multiplies by a
    /// Gaussian centred on the band diagonal and returns 0.0 for seeds on the
    /// wrong chrom or strand.
    fn seed_weight<S: Seed + Weighted>(&self, seed: &S, k: usize) -> f64;
}

// ── Partition ─────────────────────────────────────────────────────────────────

/// Partition `seeds` (sorted by `read_pos`) into connected components using
/// `scheme.reachable`.  Returns one `Vec<usize>` of seed indices per
/// component, in read-position order within each component.  Singleton
/// components (isolated seeds) are dropped — they cannot form a chain.
///
/// The forward scan per seed stops as soon as the read gap exceeds
/// `scheme.max_read_gap()`, keeping the overall cost O(n · avg_density)
/// rather than O(n²).
pub fn partition_seeds<S, Scheme>(seeds: &[S], k: usize, scheme: &Scheme) -> Vec<Vec<usize>>
where
    S: Seed,
    Scheme: SegmentScheme,
{
    let n = seeds.len();
    if n == 0 {
        return Vec::new();
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| seeds[i].read_start());

    let uf = UnionFind::new(n);

    for rank in 0..n {
        let i = order[rank];
        let i_read_end = seeds[i].read_end(k) as i64;

        for rank2 in rank + 1..n {
            let j = order[rank2];
            let read_gap = seeds[j].read_start() as i64 - i_read_end;
            if read_gap > scheme.max_read_gap() {
                break;
            }
            if scheme.reachable(&seeds[i], &seeds[j], k) {
                uf.union(i, j);
            }
        }
    }

    // Collect components keyed by root, preserving read-pos order (order is
    // already sorted by read_pos so iterating it gives components in order too).
    let mut components: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for &i in &order {
        let root = uf.find(i);
        components.entry(root).or_default().push(i);
    }

    // Drop components whose read span is below the minimum.  Span is from
    // the read_pos of the first seed to the read_end of the last seed.
    let min_span = scheme.min_segment_span();
    components
        .into_values()
        .filter(|c| {
            let first = seeds[c[0]].read_start() as i64;
            let last = seeds[*c.last().unwrap()].read_end(k) as i64;
            last - first >= min_span
        })
        .collect()
}

// ── Concrete scheme ───────────────────────────────────────────────────────────

/// Configuration for `FullSegmentScheme`.
pub struct SegmentConfig {
    /// Maximum read-space gap for two seeds to be reachable.  Typically 2*k.
    pub max_read_gap: i64,
    /// Maximum |ref_gap - read_gap| for two seeds to be reachable.  Seeds
    /// whose diagonal shift exceeds this are on different loci and should not
    /// be joined into the same component.
    pub max_ref_deviation: i64,
    /// Minimum read span a component must cover to be worth chaining.
    pub min_segment_span: i64,
    /// Cost per base of read-space gap.
    pub read_gap_cost_per_base: f64,
    /// Cost per base of ref-deviation (|ref_gap - read_gap|) on continuation edges.
    pub ref_dev_cost_per_base: f64,
    /// Cost per base of diagonal deviation between consecutive seeds on any edge
    /// where a ref_gap is defined (continuation and repeat).  Applied universally
    /// regardless of band compatibility, so off-band seeds accumulate cost
    /// proportional to how far they shift the diagonal from the previous seed.
    pub diag_dev_cost_per_base: f64,
    /// Fixed penalty for a repeat edge (same chrom/strand, ref gap too negative
    /// to bridge with SW).
    pub repeat_penalty: f64,
}

impl SegmentConfig {
    pub fn default_for_k(k: usize) -> Self {
        Self {
            max_read_gap: 2 * k as i64,
            max_ref_deviation: 200,
            min_segment_span: 2 * k as i64,
            read_gap_cost_per_base: 0.05,
            ref_dev_cost_per_base: 0.01,
            diag_dev_cost_per_base: 0.0,
            repeat_penalty: 0.0,
        }
    }
}

pub struct FullSegmentScheme {
    pub cfg: SegmentConfig,
}

impl FullSegmentScheme {
    pub fn new(cfg: SegmentConfig) -> Self {
        Self { cfg }
    }
}

impl SegmentScheme for FullSegmentScheme {
    fn max_read_gap(&self) -> i64 {
        self.cfg.max_read_gap
    }

    fn min_segment_span(&self) -> i64 {
        self.cfg.min_segment_span
    }

    fn max_ref_deviation(&self) -> i64 {
        self.cfg.max_ref_deviation
    }

    fn reachable<S: Seed>(&self, lhs: &S, rhs: &S, k: usize) -> bool {
        if lhs.chrom_id() != rhs.chrom_id() || lhs.is_reverse() != rhs.is_reverse() {
            return false;
        }
        let lhs_read_end = lhs.read_end(k) as i64;
        let read_gap = rhs.read_start() as i64 - lhs_read_end;
        if read_gap > self.cfg.max_read_gap {
            return false;
        }
        let ref_gap = if lhs.is_reverse() {
            lhs.ref_start() as i64 - rhs.ref_end(k) as i64
        } else {
            rhs.ref_start() as i64 - lhs.ref_end(k) as i64
        };
        let deviation = (ref_gap - read_gap).abs();
        deviation <= self.cfg.max_ref_deviation
    }

    fn edge_cost<S: Seed + Weighted + GapComputable>(
        &self,
        lhs: &S,
        rhs: &S,
        k: usize,
    ) -> Option<(f64, EdgeType)> {
        let (gap, lhs_eff, rhs_eff) = lhs.gap_to(rhs, k)?;
        // Cross-chrom/strand sentinel — not connectable at level 2.
        if gap.ref_gap == i64::MIN {
            return None;
        }
        let read_gap_cost = (gap.read_gap.max(0) as f64) * self.cfg.read_gap_cost_per_base;
        let deviation = (gap.ref_gap - gap.read_gap).unsigned_abs() as f64;
        let diag_dev_cost = deviation * self.cfg.diag_dev_cost_per_base;
        if gap.ref_gap < -(rhs.length(k) as i64) {
            return Some((read_gap_cost + diag_dev_cost + self.cfg.repeat_penalty, EdgeType::Repeat));
        }
        let ref_dev_cost = deviation * self.cfg.ref_dev_cost_per_base;
        let trimmed_cost = (lhs.weight() - lhs_eff.as_ref().unwrap_or(lhs).weight())
            + (rhs.weight() - rhs_eff.as_ref().unwrap_or(rhs).weight());
        Some((read_gap_cost + diag_dev_cost + ref_dev_cost + trimmed_cost, EdgeType::Continuation))
    }

    fn seed_weight<S: Seed + Weighted>(&self, seed: &S, _k: usize) -> f64 {
        seed.weight()
    }
}

// ── Banded Segment Scheme ──────────────────────────────────────────────────────

/// A `SegmentScheme` that anchors chaining to a specific diagonal band.
///
/// Node rewards are full seed weights (`seed.weight()`), unchanged regardless
/// of diagonal position.  Off-diagonal cost is folded into the *edge* penalty
/// instead, keeping node rewards large enough to dominate gap costs.
///
/// Edge costs have two modes:
///
/// * **Compatible** — both seeds share the same chrom/strand as the band and
///   both lie within `cfg.max_ref_deviation` of `central_diagonal`:
pub struct BandedSegmentScheme {
    pub cfg: SegmentConfig,
    pub chrom_id: u32,
    pub is_reverse: bool,
    pub central_diagonal: f64,
    /// Half-width of the diagonal penalty: the deviation at which the per-seed
    /// cost equals half of `max_diag_cost`.  Set to `diagonal_variance.sqrt()`
    /// so the band's own spread defines the transition zone.
    pub diag_half_width: f64,
    /// Maximum per-seed diagonal penalty (approached asymptotically for distant
    /// seeds).  Set to `sv_break_penalty / 2` so a run of off-band seeds is
    /// expensive but not as costly as repeated SV breaks.
    pub max_diag_cost: f64,
    /// Fixed penalty added to any edge that crosses a structural-variant break.
    pub sv_break_penalty: f64,
    /// Fixed penalty added to any edge that crosses a repeat boundary.
    /// Defaults to `sv_break_penalty / sqrt(2)`.
    pub repeat_penalty: f64,
    /// Maximum read-space gap between compound seeds.
    pub max_read_gap: i64,
}

impl BandedSegmentScheme {
    pub fn new(
        cfg: SegmentConfig,
        chrom_id: u32,
        is_reverse: bool,
        central_diagonal: f64,
        diagonal_variance: f64,
        sv_break_penalty: f64,
        max_diag_cost: f64,
        max_read_gap: i64,
    ) -> Self {
        let repeat_penalty = sv_break_penalty / std::f64::consts::SQRT_2;
        let diag_half_width = diagonal_variance.sqrt();
        Self {
            cfg,
            chrom_id,
            is_reverse,
            central_diagonal,
            diag_half_width,
            max_diag_cost,
            sv_break_penalty,
            repeat_penalty,
            max_read_gap,
        }
    }

    /// Returns `true` if `seed` is on the same chrom and strand as this band.
    /// The diagonal deviation is penalized continuously in `edge_cost` rather
    /// than gated here.
    fn is_compatible<S: Seed>(&self, seed: &S) -> bool {
        seed.chrom_id() == self.chrom_id && seed.is_reverse() == self.is_reverse
    }
}

impl SegmentScheme for BandedSegmentScheme {
    fn max_read_gap(&self) -> i64 {
        self.max_read_gap
    }
    fn min_segment_span(&self) -> i64 {
        self.cfg.min_segment_span
    }
    fn max_ref_deviation(&self) -> i64 {
        self.cfg.max_ref_deviation
    }

    fn reachable<S: Seed>(&self, lhs: &S, rhs: &S, k: usize) -> bool {
        let lhs_read_end = lhs.read_end(k) as i64;
        let read_gap = rhs.read_start() as i64 - lhs_read_end;
        read_gap <= self.max_read_gap
    }

    fn edge_cost<S: Seed + Weighted + GapComputable>(
        &self,
        lhs: &S,
        rhs: &S,
        k: usize,
    ) -> Option<(f64, EdgeType)> {
        let read_gap_cost = {
            let lhs_read_end = lhs.read_end(k) as i64;
            (rhs.read_start() as i64 - lhs_read_end).max(0) as f64 * self.cfg.read_gap_cost_per_base
        };

        let (gap, lhs_eff, rhs_eff) = lhs.gap_to(rhs, k)?;

        if gap.ref_gap == i64::MIN {
            // SV break: lhs and rhs are on different chrom or strand from each other.
            Some((read_gap_cost + self.sv_break_penalty, EdgeType::SvBreak))
        } else {
            let deviation = (gap.ref_gap - gap.read_gap).unsigned_abs() as f64;
            let diag_dev_cost = deviation * self.cfg.diag_dev_cost_per_base;
            if gap.ref_gap < -(rhs.length(k) as i64) {
                // Repeat: same chrom/strand but ref displacement too large to bridge with SW.
                return Some((
                    read_gap_cost + diag_dev_cost + self.repeat_penalty,
                    EdgeType::Repeat,
                ));
            }
            let read_gap = rhs.read_start() as i64 - lhs.read_end(k) as i64;
            if read_gap > self.cfg.max_read_gap {
                // Gap too large for SW alignment — treat as an SV break.
                return Some((
                    read_gap_cost + self.sv_break_penalty,
                    EdgeType::SvBreak,
                ));
            }
            // Continuation: normal gap costs, SW-bridgeable.
            let ref_dev_cost = deviation * self.cfg.ref_dev_cost_per_base;
            let trimmed_cost = (lhs.weight() - lhs_eff.as_ref().unwrap_or(lhs).weight())
                + (rhs.weight() - rhs_eff.as_ref().unwrap_or(rhs).weight());
            Some((
                read_gap_cost + diag_dev_cost + ref_dev_cost + trimmed_cost,
                EdgeType::Continuation,
            ))
        }
    }

    fn seed_weight<S: Seed + Weighted>(&self, seed: &S, _k: usize) -> f64 {
        let diag_penalty = if self.is_compatible(seed) {
            let z = seed.diagonal() as f64 - self.central_diagonal;
            let t = z * z / (self.diag_half_width * self.diag_half_width);
            self.max_diag_cost * t.tanh()
        } else {
            self.max_diag_cost
        };
        (seed.weight() - diag_penalty).max(0.0)
    }
}

// ── Segment chain result ──────────────────────────────────────────────────────

/// A single chain extracted from a connected component of compound seeds.
pub struct SegmentChain {
    /// DP score of this chain.
    pub score: f64,
    /// Indices into the original seed slice, in read-position order.
    pub chain: Vec<usize>,
    /// Edge types. The edge_type[i] value denotes the type of edge between chain[i] and chain[i + 1]
    pub edge_type: Vec<EdgeType>
}

// ── Per-component chaining DP ─────────────────────────────────────────────────

/// Run the chaining DP over one component (a slice of seed indices into
/// `seeds`, already sorted by read_pos) and return all chains found by
/// iterative extraction.
///
/// Each iteration runs the DP over the unassigned seeds in the component,
/// extracts the highest-scoring chain, marks those seeds as used, and
/// repeats.  Iteration stops when no remaining seed can start a chain with
/// a positive score (i.e. every remaining seed is either isolated or would
/// lose weight to gap costs).
///
/// Seeds in `component` must be sorted by `read_pos`.  `seeds` is the full
/// seed slice that `component` indexes into.
pub fn chain_component<S, Scheme>(
    seeds: &[S],
    component: &[usize],
    k: usize,
    scheme: &Scheme,
) -> Vec<SegmentChain>
where
    S: Seed + Weighted + GapComputable,
    Scheme: SegmentScheme,
{
    let mut used = vec![false; component.len()];
    let mut result = Vec::new();

    loop {
        // Collect indices of unassigned seeds in this component, in read-pos order.
        let active: Vec<usize> = (0..component.len()).filter(|&r| !used[r]).collect();

        if active.len() < 2 {
            break;
        }

        // DP over active seeds.  dp[r] is the best score ending at active[r].
        let n = active.len();
        let mut dp = vec![0.0f64; n];
        let mut prev = vec![usize::MAX; n];
        let mut prev_edge: Vec<Option<EdgeType>> = (0..n).map(|_| None).collect();

        for rank in 0..n {
            let i = component[active[rank]];
            let w_i = scheme.seed_weight(&seeds[i], k);
            dp[rank] = w_i;

            for r in (0..rank).rev() {
                let j = component[active[r]];

                let edge = scheme.edge_cost(&seeds[j], &seeds[i], k);
                if let Some((cost, edge_type)) = edge {
                    let candidate = dp[r] + w_i - cost;
                    if candidate > dp[rank] {
                        dp[rank] = candidate;
                        prev[rank] = r;
                        prev_edge[rank] = Some(edge_type);
                    }
                }
            }
        }

        // Best endpoint.
        let best_rank = (0..n)
            .max_by(|&a, &b| dp[a].partial_cmp(&dp[b]).unwrap())
            .unwrap();

        if dp[best_rank] <= 0.0 {
            break;
        }

        // Traceback.
        let mut chain_ranks = Vec::new();
        let mut cur = best_rank;
        loop {
            chain_ranks.push(cur);
            let p = prev[cur];
            if p == usize::MAX {
                break;
            }
            cur = p;
        }
        chain_ranks.reverse();

        // Single-seed chains add nothing beyond what level-1 already captured.
        if chain_ranks.len() < 2 {
            break;
        }


        // Collect edge types in forward order: edge_types[i] is the edge from
        // chain_ranks[i] to chain_ranks[i+1], stored in prev_edge[chain_ranks[i+1]].
        let edge_types: Vec<EdgeType> = chain_ranks[1..]
            .iter()
            .map(|&r| prev_edge[r].take().unwrap_or(EdgeType::Continuation))
            .collect();

        // Mark seeds used and record the chain with original seed indices.
        let chain_indices: Vec<usize> = chain_ranks
            .iter()
            .map(|&r| {
                used[active[r]] = true;
                component[active[r]]
            })
            .collect();

        result.push(SegmentChain {
            score: dp[best_rank],
            chain: chain_indices,
            edge_type: edge_types,
        });
    }

    // Emit any remaining unassigned seeds as single-seed segments.  These are
    // either genuine singletons (no reachable neighbour in the component) or
    // seeds left over after all multi-seed chains were extracted.  Both carry
    // real alignment evidence and must be visible to the level-3 pass.
    for r in 0..component.len() {
        if !used[r] {
            let idx = component[r];
            result.push(SegmentChain {
                score: scheme.seed_weight(&seeds[idx], k),
                chain: vec![idx],
                edge_type: vec![],
            });
        }
    }

    result
}

/// Partition `seeds` into components and extract all chains from each
/// component.  Convenience wrapper around `partition_seeds` + `chain_component`.
pub fn find_all_chains<S, Scheme>(seeds: &[S], k: usize, scheme: &Scheme) -> Vec<SegmentChain>
where
    S: Seed + Weighted + GapComputable,
    Scheme: SegmentScheme,
{
    let components = partition_seeds(seeds, k, scheme);
    let mut chains = Vec::new();
    for component in &components {
        chains.extend(chain_component(seeds, component, k, scheme));
    }
    chains
}

/// Extract chains, removing used seeds progressively.
pub fn extract_chains<S, Scheme>(seeds: &[S], k: usize, scheme: &Scheme) -> Vec<SegmentChain>
where
    S: Seed + Weighted + GapComputable,
    Scheme: SegmentScheme,
{
    let n = seeds.len();
    let live_seeds: Vec<usize> = (0..n).collect();
    chain_component(seeds, &live_seeds, k, scheme)
}

// ── Bits-scored chaining ──────────────────────────────────────────────────────

/// log2 of the human genome size: the bits needed to name an arbitrary
/// position, charged when a break changes chromosome.
const GENOME_BITS: f64 = 31.5;

/// Seedless gaps are aligned only from this many of a seed's nearest
/// same-locus predecessors (nearest in read space).
const GAP_ALIGN_PREDECESSORS: usize = 3;

/// A seed's bounds including its verified exact-match extension: the read
/// interval, and the reference coordinate at its read-left (`entry`) and
/// read-right (`exit`) ends.
struct Extended {
    read_start: i64,
    read_end: i64,
    entry: i64,
    exit: i64,
    chrom_id: u32,
    is_reverse: bool,
    diagonal: i64,
}

impl Extended {
    fn new<S: Seed>(seed: &S, k: usize) -> Self {
        let (left, right) = seed.extension();
        let (left, right) = (left as i64, right as i64);
        let (entry, exit) = if seed.is_reverse() {
            (seed.ref_end(k) as i64 + left, seed.ref_start() as i64 - right)
        } else {
            (seed.ref_start() as i64 - left, seed.ref_end(k) as i64 + right)
        };
        Self {
            read_start: seed.read_start() as i64 - left,
            read_end: seed.read_end(k) as i64 + right,
            entry,
            exit,
            chrom_id: seed.chrom_id(),
            is_reverse: seed.is_reverse(),
            diagonal: seed.diagonal(),
        }
    }
}

/// Score an alignment in bits: matched bases earn `match_bits`, mismatches
/// cost `mismatch_bits`, and each indel costs `gap_open_bits + gap_extend_bits * len`.
pub fn alignment_bits(cigar: &[Op], cfg: &ChainingConfig) -> f64 {
    cigar
        .iter()
        .map(|op| {
            let n = op.len() as f64;
            match op.kind() {
                Kind::SequenceMatch | Kind::Match => cfg.match_bits * n,
                Kind::SequenceMismatch => -cfg.mismatch_bits * n,
                Kind::Insertion | Kind::Deletion => -(cfg.gap_open_bits + cfg.gap_extend_bits * n),
                _ => 0.0,
            }
        })
        .sum()
}

/// Find the best explanation of the read as a chain of `seeds`, scored in bits.
///
/// Appending seed `i` after `j` earns `match_bits` per read base that `i`
/// covers beyond the end of `j`, so no read base is credited twice whatever
/// kind of edge joins them; `i` must end after `j`.  The edge costs
/// `-log2 P(move)`, the cheaper of:
///
/// * an indel on the same chrom/strand: `mismatch_bits` for a gap on the same
///   diagonal (extended seeds are maximal exact matches, so a gap holds at
///   least one difference), otherwise `gap_open_bits + gap_extend_bits * |shift|`;
/// * a structural break: `sv_break_bits + log2(1 + |shift|)` on the same
///   chrom/strand, `sv_break_bits + 1 + log2(1 + distance)` across strands,
///   and `sv_break_bits + GENOME_BITS` across chromosomes.
///
/// Seedless gaps between seeds on the same locus are also scored by aligning
/// them: `gap_bits(lhs_read_end, rhs_read_start, ref_lo, ref_hi, chrom_id,
/// is_reverse)` returns the bits of the gap alignment, which is used when it
/// beats treating the gap bases as unexplained.
///
/// Edge types in the result say how the chain is emitted: consecutive seeds
/// on the same locus within `max_bridge` are a `Continuation`, a large
/// backward reference jump is a `Repeat`, and anything else is an `SvBreak`.
pub fn chain_bits<S, F>(
    seeds: &[S],
    k: usize,
    cfg: &ChainingConfig,
    mut gap_bits: F,
) -> Option<SegmentChain>
where
    S: Seed,
    F: FnMut(u32, u32, u32, u32, u32, bool) -> Option<f64>,
{
    let ext: Vec<Extended> = seeds.iter().map(|s| Extended::new(s, k)).collect();
    let mut order: Vec<usize> = (0..seeds.len()).collect();
    order.sort_unstable_by_key(|&i| (ext[i].read_end, ext[i].read_start));

    let n = order.len();
    let mut dp = vec![0.0f64; n];
    let mut prev = vec![usize::MAX; n];

    for r in 0..n {
        let i = &ext[order[r]];
        dp[r] = cfg.match_bits * (i.read_end - i.read_start) as f64;
        // Ranks of the nearest same-locus predecessors across a seedless gap.
        let mut nearest: Vec<usize> = Vec::with_capacity(GAP_ALIGN_PREDECESSORS + 1);

        for q in 0..r {
            let j = &ext[order[q]];
            if j.read_end >= i.read_end {
                continue;
            }
            let new_bases = i.read_end - i.read_start.max(j.read_end);
            let read_gap = i.read_start - j.read_end;
            let cost = if j.chrom_id != i.chrom_id {
                cfg.sv_break_bits + GENOME_BITS
            } else if j.is_reverse != i.is_reverse {
                cfg.sv_break_bits + 1.0 + (1.0 + (i.entry - j.exit).abs() as f64).log2()
            } else {
                let shift = if i.is_reverse {
                    j.diagonal - i.diagonal
                } else {
                    i.diagonal - j.diagonal
                }
                .abs();
                let indel = match (shift, read_gap > 0) {
                    (0, true) => cfg.mismatch_bits,
                    (0, false) => 0.0,
                    _ => cfg.gap_open_bits + cfg.gap_extend_bits * shift as f64,
                };
                let ref_gap = if i.is_reverse { j.exit - i.entry } else { i.entry - j.exit };
                if read_gap >= cfg.gap_align_min as i64
                    && read_gap <= cfg.gap_align_max as i64
                    && ref_gap >= 0
                    && ref_gap <= cfg.gap_align_max as i64
                {
                    nearest.push(q);
                    if nearest.len() > GAP_ALIGN_PREDECESSORS {
                        nearest.remove(0);
                    }
                }
                indel.min(cfg.sv_break_bits + (1.0 + shift as f64).log2())
            };
            let candidate = dp[q] + cfg.match_bits * new_bases as f64 - cost;
            if candidate > dp[r] {
                dp[r] = candidate;
                prev[r] = q;
            }
        }

        for &q in &nearest {
            let j = &ext[order[q]];
            let (ref_lo, ref_hi) = if i.is_reverse { (i.entry, j.exit) } else { (j.exit, i.entry) };
            let Some(bits) = gap_bits(
                j.read_end as u32,
                i.read_start as u32,
                ref_lo as u32,
                ref_hi as u32,
                i.chrom_id,
                i.is_reverse,
            ) else {
                continue;
            };
            let candidate = dp[q] + cfg.match_bits * (i.read_end - i.read_start) as f64 + bits;
            if candidate > dp[r] {
                dp[r] = candidate;
                prev[r] = q;
            }
        }
    }

    let best = (0..n).max_by(|&a, &b| dp[a].total_cmp(&dp[b]))?;
    let mut ranks = Vec::new();
    let mut cur = best;
    loop {
        ranks.push(cur);
        if prev[cur] == usize::MAX {
            break;
        }
        cur = prev[cur];
    }
    ranks.reverse();

    let chain: Vec<usize> = ranks.iter().map(|&r| order[r]).collect();
    let chain = inline_tandem_copies(&chain, &ext, cfg.max_bridge as i64);
    let edge_type = chain
        .windows(2)
        .map(|w| bits_edge_type(&ext[w[0]], &ext[w[1]], cfg.max_bridge as i64))
        .collect();
    Some(SegmentChain {
        score: dp[best],
        chain,
        edge_type,
    })
}

/// Drop chain seeds that lie wholly within reference already passed by the
/// previous kept seed on the same locus, stepping back no more than
/// `max_bridge`: the extra copies of a tandem expansion.  The read bases they
/// explain are then aligned within one record, as an insertion, by the gap
/// aligner.  This shapes the output only; the chain's score is unchanged.
fn inline_tandem_copies(chain: &[usize], ext: &[Extended], max_bridge: i64) -> Vec<usize> {
    let mut kept: Vec<usize> = Vec::with_capacity(chain.len());
    for &s in chain {
        if let Some(&l) = kept.last() {
            let (lhs, rhs) = (&ext[l], &ext[s]);
            if lhs.chrom_id == rhs.chrom_id && lhs.is_reverse == rhs.is_reverse {
                // How far rhs moves the reference frontier on, and how far
                // behind the frontier it starts, in the direction of alignment.
                let (advance, back) = if rhs.is_reverse {
                    (lhs.exit - rhs.exit, rhs.entry - lhs.exit)
                } else {
                    (rhs.exit - lhs.exit, lhs.exit - rhs.entry)
                };
                if advance <= 0 && back <= max_bridge {
                    continue;
                }
            }
        }
        kept.push(s);
    }
    kept
}

/// How consecutive chain seeds are emitted: within one record across a
/// bridgeable gap on the same locus, otherwise as a split.  Uses the extended
/// bounds, as the DP and `compute_chain_bounds` do.
fn bits_edge_type(lhs: &Extended, rhs: &Extended, max_bridge: i64) -> EdgeType {
    if lhs.chrom_id != rhs.chrom_id || lhs.is_reverse != rhs.is_reverse {
        return EdgeType::SvBreak;
    }
    let read_gap = rhs.read_start - lhs.read_end;
    let ref_gap = if lhs.is_reverse { lhs.exit - rhs.entry } else { rhs.entry - lhs.exit };
    if ref_gap < -(rhs.read_end - rhs.read_start) {
        EdgeType::Repeat
    } else if read_gap > max_bridge || (ref_gap - read_gap).abs() > max_bridge {
        EdgeType::SvBreak
    } else {
        EdgeType::Continuation
    }
}

// ── Diagonal support weights ──────────────────────────────────────────────────

/// Compute a diagonal-support weight for each segment chain.
///
/// For each chain i, the support is the sum over all other chains j on the
/// same chrom and strand of `chain_j.score * exp(-0.5 * ((d_i - d_j) / sigma)²)`,
/// where `d` is the diagonal of the first seed of each chain.
///
/// Chains are sorted by `(chrom_id, is_reverse, diagonal)` and the inner scan
/// breaks early when `|d_i - d_j| > cutoff` (default 4σ), keeping the overall
/// cost O(n · avg_neighbours_within_cutoff).
///
/// Returns a `Vec<f64>` parallel to `chains`, with self-contribution included.
pub fn diagonal_support<S>(seeds: &[S], chains: &[SegmentChain], sigma: f64) -> Vec<f64>
where
    S: Seed,
{
    let cutoff = 4.0 * sigma;
    let two_sigma_sq = 2.0 * sigma * sigma;
    let n = chains.len();

    // For each chain, extract the diagonal of its first seed plus chrom/strand.
    let diag: Vec<(u32, bool, i64)> = chains
        .iter()
        .map(|c| {
            let s = &seeds[c.chain[0]];
            (s.chrom_id(), s.is_reverse(), s.diagonal())
        })
        .collect();

    // Index sorted by (chrom_id, is_reverse, diagonal).
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| diag[i]);

    let mut support = vec![0.0f64; n];

    for rank in 0..n {
        let i = order[rank];
        let (chrom_i, rev_i, d_i) = diag[i];
        let score_i = chains[i].score;

        // Scan forward.
        for rank2 in rank + 1..n {
            let j = order[rank2];
            let (chrom_j, rev_j, d_j) = diag[j];
            if chrom_j != chrom_i || rev_j != rev_i {
                break;
            }
            let delta = (d_j - d_i) as f64;
            if delta > cutoff {
                break;
            }
            let w = chains[j].score * (-delta * delta / two_sigma_sq).exp();
            support[i] += w;
            // Symmetric: j also gets i's contribution.
            support[j] += score_i * (-delta * delta / two_sigma_sq).exp();
        }

        // Self-contribution.
        support[i] += score_i;
    }

    support
}

// ── Alignment DP ──────────────────────────────────────────────────────────────

/// The result of the alignment DP: the best non-overlapping subset of segment
/// chains covering the read.
pub struct AlignmentResult {
    /// Total DP score.
    pub score: f64,
    /// Indices into the `chains` slice, in read-position order.
    pub segments: Vec<usize>,
}

/// Select the best non-overlapping subset of `chains` to form an alignment.
///
/// Scores each chain as `chain.score + diag_weight * support[i]` and penalises
/// uncovered gaps between consecutive selected segments at `gap_cost_per_base`
/// per base.  Segments that overlap in read space by more than
/// `overlap_tolerance` bases are treated as incompatible.
///
/// `seeds` is the compound-seed slice that chain indices refer to.
/// `support` must be parallel to `chains` (from `diagonal_support`).
pub fn best_alignment<S>(
    seeds: &[S],
    chains: &[SegmentChain],
    support: &[f64],
    k: usize,
    diag_weight: f64,
    gap_cost_per_base: f64,
    overlap_tolerance: i64,
) -> Option<AlignmentResult>
where
    S: Seed,
{
    let n = chains.len();
    if n == 0 {
        return None;
    }

    // Pre-compute read_start / read_end for each chain from the seed slice.
    let read_start = |i: usize| seeds[chains[i].chain[0]].read_start() as i64;
    let read_end = |i: usize| seeds[*chains[i].chain.last().unwrap()].read_end(k) as i64;

    // Sort by read_end, breaking ties by read_start descending.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| (read_end(i), -(read_start(i))));

    let node_weight = |i: usize| chains[i].score + diag_weight * support[i];

    let inf = f64::NEG_INFINITY;
    let mut dp = vec![inf; n];
    let mut prev = vec![usize::MAX; n];

    for rank in 0..n {
        let i = order[rank];
        let rs_i = read_start(i);

        // Chain i alone.
        dp[rank] = node_weight(i);

        // Try all earlier chains as predecessor.
        for r in (0..rank).rev() {
            let j = order[r];
            let re_j = read_end(j);
            if re_j > rs_i + overlap_tolerance {
                continue; // incompatible overlap
            }
            if dp[r] == inf {
                continue;
            }
            let gap = (rs_i - re_j).max(0) as f64;
            let candidate = dp[r] + node_weight(i) - gap_cost_per_base * gap;
            if candidate > dp[rank] {
                dp[rank] = candidate;
                prev[rank] = r;
            }
        }
    }

    let best_rank = (0..n)
        .filter(|&r| dp[r] != inf)
        .max_by(|&a, &b| dp[a].partial_cmp(&dp[b]).unwrap())?;

    let mut segments = Vec::new();
    let mut cur = best_rank;
    loop {
        segments.push(order[cur]);
        let p = prev[cur];
        if p == usize::MAX {
            break;
        }
        cur = p;
    }
    segments.reverse();

    Some(AlignmentResult {
        score: dp[best_rank],
        segments,
    })
}

// ── Banded assembly ───────────────────────────────────────────────────────────

/// One band's contribution to the final assembled alignment.
pub struct BandAlignment {
    /// Which diagonal band (index into the bands slice).
    pub band_idx: usize,
    /// The alignment result for this band.
    pub result: AlignmentResult,
    /// Read start of the covered span.
    pub read_start: u32,
    /// Read end of the covered span.
    pub read_end: u32,
}

/// Assemble an alignment from per-band optimal chains.
///
/// For each diagonal band, collect only those `chains` whose first seed's
/// diagonal falls within the band, run `best_alignment` over them, and record
/// the covered read span.  Then sort bands by covered read span (descending)
/// and greedily accept bands that don't overlap any already-accepted segment.
/// Stop when the band under consideration covers less than 50% of the total
/// read length.
///
/// Returns the accepted band alignments in the order they were accepted.
pub fn assemble_alignment<S>(
    seeds: &[S],
    chains: &[SegmentChain],
    support: &[f64],
    bands: &[DiagonalBand],
    k: usize,
    read_length: u32,
    diag_weight: f64,
    gap_cost_per_base: f64,
) -> Vec<BandAlignment>
where
    S: Seed,
{
    if bands.is_empty() || chains.is_empty() {
        return Vec::new();
    }

    // For each band, collect chain indices whose first seed's diagonal is in-band.
    let mut band_results: Vec<(usize, AlignmentResult, u32, u32)> = Vec::new();

    for (bi, band) in bands.iter().enumerate() {
        let band_set: std::collections::HashSet<usize> = band.members.iter().copied().collect();

        // Chains whose first seed is a member of this band.
        let band_chain_indices: Vec<usize> = (0..chains.len())
            .filter(|&ci| band_set.contains(&chains[ci].chain[0]))
            .collect();

        if band_chain_indices.is_empty() {
            continue;
        }

        // Build local slices for the band's chains and support.
        let band_chains: Vec<&SegmentChain> =
            band_chain_indices.iter().map(|&ci| &chains[ci]).collect();
        let band_support: Vec<f64> = band_chain_indices.iter().map(|&ci| support[ci]).collect();

        // Build owned SegmentChain vec for best_alignment (needs &[SegmentChain]).
        let owned_chains: Vec<SegmentChain> = band_chains
            .iter()
            .map(|c| SegmentChain {
                score: c.score,
                chain: c.chain.clone(),
                edge_type: c.edge_type.clone(),
            })
            .collect();

        let Some(result) = best_alignment(
            seeds,
            &owned_chains,
            &band_support,
            k,
            diag_weight,
            gap_cost_per_base,
            0, // exact non-overlap for assembly
        ) else {
            continue;
        };

        // Compute covered read span from the alignment's selected chains.
        let rs = result
            .segments
            .iter()
            .map(|&si| seeds[owned_chains[si].chain[0]].read_start())
            .min()
            .unwrap();
        let re = result
            .segments
            .iter()
            .map(|&si| seeds[*owned_chains[si].chain.last().unwrap()].read_end(k) as u32)
            .max()
            .unwrap();

        // Remap segment indices back to global chain indices.
        let remapped = AlignmentResult {
            score: result.score,
            segments: result
                .segments
                .iter()
                .map(|&si| band_chain_indices[si])
                .collect(),
        };

        band_results.push((bi, remapped, rs, re));
    }

    // Sort by covered read span descending.
    band_results.sort_by(|a, b| {
        let span_b = b.3 - b.2;
        let span_a = a.3 - a.2;
        span_b.cmp(&span_a)
    });

    let half_read = read_length / 2;
    let mut accepted: Vec<BandAlignment> = Vec::new();
    let mut occupied: Vec<(u32, u32)> = Vec::new();
    let mut used = vec![false; band_results.len()];

    let intervals_of = |result: &AlignmentResult| -> Vec<(u32, u32)> {
        result
            .segments
            .iter()
            .map(|&ci| {
                let chain = &chains[ci];
                let rs = seeds[chain.chain[0]].read_start();
                let re = seeds[*chain.chain.last().unwrap()].read_end(k) as u32;
                (rs, re)
            })
            .collect()
    };

    let overlaps_any = |ivs: &[(u32, u32)], occupied: &[(u32, u32)]| -> bool {
        ivs.iter()
            .any(|&(s, e)| occupied.iter().any(|&(os, oe)| s < oe && e > os))
    };

    // Outer loop: iterate over bands in descending span order.
    // Halt when the current band covers less than 50% of read length.
    for outer in 0..band_results.len() {
        if used[outer] {
            continue;
        }
        let span = band_results[outer].3 - band_results[outer].2;
        if span < half_read {
            break;
        }

        // Accept the outer band and mark its read intervals as occupied.
        let ivs = intervals_of(&band_results[outer].1);
        occupied.extend_from_slice(&ivs);
        used[outer] = true;
        let (bi, ref result, rs, re) = band_results[outer];
        accepted.push(BandAlignment {
            band_idx: bi,
            result: AlignmentResult {
                score: result.score,
                segments: result.segments.clone(),
            },
            read_start: rs,
            read_end: re,
        });

        // Inner loop: scan remaining bands for non-overlapping gap-fillers.
        for inner in (outer + 1)..band_results.len() {
            if used[inner] {
                continue;
            }
            let ivs = intervals_of(&band_results[inner].1);
            if overlaps_any(&ivs, &occupied) {
                continue;
            }

            occupied.extend_from_slice(&ivs);
            used[inner] = true;
            let (bi, ref result, rs, re) = band_results[inner];
            accepted.push(BandAlignment {
                band_idx: bi,
                result: AlignmentResult {
                    score: result.score,
                    segments: result.segments.clone(),
                },
                read_start: rs,
                read_end: re,
            });
        }
    }

    accepted
}

// ── Diagonal bands ────────────────────────────────────────────────────────────

/// A group of compound seeds that share a common diagonal band.
pub struct DiagonalBand {
    /// Chrom id shared by all seeds in this band.
    pub chrom_id: u32,
    /// Strand shared by all seeds in this band.
    pub is_reverse: bool,
    /// Minimum ref_pos across all seeds in the band.
    pub ref_min: u32,
    /// Maximum ref_end across all seeds in the band.
    pub ref_max: u32,
    /// Weighted-average diagonal (weight = seed weight).
    pub central_diagonal: f64,
    /// Weighted variance of diagonal values across band members.
    pub diagonal_variance: f64,
    /// coverage of the band in read space.
    pub coverage: f64,
    /// Indices into the seed slice.
    pub members: Vec<usize>,
}

/// Partition `seeds` into diagonal bands by sorting on `(chrom_id, is_reverse,
/// diagonal)` and grouping consecutive seeds whose diagonal differs by less
/// than `max_deviation` from the first seed of the current group.
pub fn partition_by_diagonal<S>(
    seeds: &[S],
    k: usize,
    max_deviation: i64,
    read_len: usize,
) -> Vec<DiagonalBand>
where
    S: Seed + Weighted,
{
    let n = seeds.len();
    if n == 0 {
        return Vec::new();
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| {
        (
            seeds[i].chrom_id(),
            seeds[i].is_reverse(),
            seeds[i].diagonal(),
        )
    });

    let mut bands: Vec<DiagonalBand> = Vec::new();
    let mut group_start = 0usize;

    while group_start < n {
        let first = order[group_start];
        let chrom = seeds[first].chrom_id();
        let rev = seeds[first].is_reverse();
        let d0 = seeds[first].diagonal();

        let mut group_end = group_start + 1;
        while group_end < n {
            let j = order[group_end];
            if seeds[j].chrom_id() != chrom || seeds[j].is_reverse() != rev {
                break;
            }
            if (seeds[j].diagonal() - d0).abs() >= max_deviation {
                break;
            }
            group_end += 1;
        }

        let members: Vec<usize> = order[group_start..group_end].to_vec();
        let band = diagonal_band_stats(seeds, &members, k, chrom, rev, read_len);
        bands.push(band);
        group_start = group_end;
    }

    bands
}

/// Compute summary statistics for a group of seeds forming one diagonal band.
fn diagonal_band_stats<S>(
    seeds: &[S],
    members: &[usize],
    k: usize,
    chrom_id: u32,
    is_reverse: bool,
    read_len: usize,
) -> DiagonalBand
where
    S: Seed + Weighted,
{
    let mut ref_min = u32::MAX;
    let mut ref_max = 0u32;
    let mut weight_sum = 0.0f64;
    let mut diag_sum = 0.0f64;
    let mut diag_sum2 = 0.0f64;

    for &i in members {
        let s = &seeds[i];
        let w = s.weight();
        let d = s.diagonal() as f64;
        ref_min = ref_min.min(s.ref_start());
        ref_max = ref_max.max(s.ref_end(k));
        weight_sum += w;
        diag_sum += w * d;
        diag_sum2 += w * d * d;
    }

    let central_diagonal = if weight_sum > 0.0 {
        diag_sum / weight_sum
    } else {
        0.0
    };

    // Weighted variance: E[d^2] - E[d]^2.
    // Floored at 1.0 (1 bp^2) to guard against zero variance when all seeds are
    // perfectly collinear — ensures the normalized diagonal cost stays finite.
    let diagonal_variance = if weight_sum > 0.0 {
        (diag_sum2 / weight_sum - central_diagonal * central_diagonal).max(1.0)
    } else {
        1.0
    };

    let mut coverage = SetCoverage::new();
    for &i in members {
        let s = &seeds[i];
        coverage.insert(s.read_start() as usize..s.read_end(k) as usize);
    }
    let coverage = coverage.coverage() as f64 / (read_len as f64);
    DiagonalBand {
        chrom_id,
        is_reverse,
        ref_min,
        ref_max,
        central_diagonal,
        diagonal_variance,
        coverage,
        members: members.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads::compound::{AtomicSeed, SeedCollection};
    use confique::Config;

    use crate::reads::test_helpers::{load_seeds, rows_to_atomic};

    fn make_seed(read_pos: u32, ref_pos: u32, chrom_id: u32, k: usize) -> AtomicSeed {
        AtomicSeed::new(read_pos, 10000, k, chrom_id, ref_pos, false, 0, 1)
    }

    // Seeds must be spaced > k apart on different diagonals to remain distinct
    // after compound_seeds() — atoms on the same diagonal within k of each
    // other are merged into a single compound.
    #[test]
    fn short_span_component_dropped() {
        let k = 20;
        // Two seeds close together (span = 25+20 = 45 >= min_span=40) → kept.
        // The distant seed (diagonal 100) is isolated; span = k = 20 < 40 → dropped.
        let collection = SeedCollection::new(
            k,
            vec![
                make_seed(0, 0, 0, k),     // diagonal 0
                make_seed(25, 25, 0, k),   // diagonal 0, read gap=5 → not merged
                make_seed(500, 600, 0, k), // diagonal 100, isolated
            ],
        );
        let compounds = collection.compound_seeds();
        let scheme = FullSegmentScheme::new(SegmentConfig::default_for_k(k));
        let components = partition_seeds(&compounds, k, &scheme);
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].len(), 2);
    }

    #[test]
    fn cross_chrom_not_reachable() {
        let k = 20;
        // Two seeds on different chroms: each is an unreachable singleton.
        // Space them > k apart so compound_seeds keeps them separate.
        let collection = SeedCollection::new(
            k,
            vec![
                make_seed(0, 0, 0, k),
                make_seed(25, 25, 1, k), // different chrom
            ],
        );
        let compounds = collection.compound_seeds();
        let scheme = FullSegmentScheme::new(SegmentConfig::default_for_k(k));
        let components = partition_seeds(&compounds, k, &scheme);
        assert!(
            components.is_empty(),
            "cross-chrom seeds should both be singletons"
        );
    }

    #[test]
    fn single_component_yields_one_chain() {
        let k = 20;
        // Four seeds on the same diagonal, evenly spaced within 2k of each other.
        // All should end up in one chain.
        let collection = SeedCollection::new(
            k,
            vec![
                make_seed(0, 0, 0, k),
                make_seed(25, 25, 0, k),
                make_seed(50, 50, 0, k),
                make_seed(75, 75, 0, k),
            ],
        );
        let compounds = collection.compound_seeds();
        let scheme = FullSegmentScheme::new(SegmentConfig::default_for_k(k));
        let chains = find_all_chains(&compounds, k, &scheme);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].chain.len(), 4);
    }

    #[test]
    fn two_interleaved_diagonals_yield_two_chains() {
        let k = 20;
        // Two sets of seeds on distinct diagonals, interleaved in read space.
        // Same-diagonal edges have a small gap (read gap = 1 bp → cost ≈ 0.05)
        // while cross-diagonal edges have deviation = 500 → cost = 5.0, so the
        // DP strongly prefers same-diagonal connections.
        //
        // Diagonal 0:   read=0,21,42   ref=0,21,42    (gap=1 → cost 0.05+0.01=0.06)
        // Diagonal 500: read=10,31,52  ref=510,531,552 (gap=1, same low cost within diag)
        //
        // Cross-diagonal: deviation ≈ 500 → ref_dev cost = 5.0, so cross edges lose weight.
        let collection = SeedCollection::new(
            k,
            vec![
                make_seed(0, 0, 0, k),    // diag 0
                make_seed(10, 510, 0, k), // diag 500
                make_seed(21, 21, 0, k),  // diag 0
                make_seed(31, 531, 0, k), // diag 500
                make_seed(42, 42, 0, k),  // diag 0
                make_seed(52, 552, 0, k), // diag 500
            ],
        );
        let compounds = collection.compound_seeds();
        let scheme = FullSegmentScheme::new(SegmentConfig::default_for_k(k));
        let chains = find_all_chains(&compounds, k, &scheme);
        assert_eq!(chains.len(), 2);
        for c in &chains {
            assert_eq!(c.chain.len(), 3);
        }
    }

    /// Ad-hoc exploration test for level-2 segment chaining.  Loads seeds from
    /// `PARALLAX_TEST_SEEDS`, partitions into connected components, then runs
    /// the iterative chaining DP and prints the resulting segments.
    ///
    /// Run with:
    ///   PARALLAX_TEST_SEEDS=/path/to/seeds.tsv[.gz] \
    ///     cargo test -p parallax segments_adhoc -- --nocapture --ignored
    #[test]
    #[ignore]
    fn segments_adhoc() {
        let path = match std::env::var("PARALLAX_TEST_SEEDS") {
            Ok(p) => p,
            Err(_) => {
                eprintln!("PARALLAX_TEST_SEEDS not set — skipping");
                return;
            }
        };
        if !std::path::Path::new(&path).exists() {
            eprintln!("seed file not found at {path} — skipping");
            return;
        }

        let k = 20;
        let rows = load_seeds(&path);
        let (atoms, chrom_names) = rows_to_atomic(&rows, k);
        let chrom_name = |id: u32| {
            chrom_names
                .get(id as usize)
                .map(|s| s.as_str())
                .unwrap_or("?")
        };
        let read_name = rows[0].read_id.clone();

        let read_length = atoms
            .iter()
            .map(|s| s.read_end(k))
            .max()
            .unwrap_or(0) as usize;

        let collection = SeedCollection::new(k, atoms);
        let compounds = collection.compound_seeds();
        println!("=== Compound seeds: {} ===", compounds.len());
        println!("rank\tread_start\tread_end\tchrom\tstrand\tref_start\tref_end\tatoms\tread_span");
        for (ci, cs) in compounds.iter().enumerate() {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                ci,
                cs.read_start(),
                cs.read_end(k),
                chrom_name(cs.chrom_id()),
                if cs.is_reverse() { "-" } else { "+" },
                cs.ref_start(),
                cs.ref_end(k),
                cs.atoms().len(),
                cs.read_end(k) - cs.read_start()
            );
        }

        let scheme = FullSegmentScheme::new(SegmentConfig::default_for_k(k));

        let components = partition_seeds(&compounds, k, &scheme);
        println!(
            "=== Components (after span filter): {} ===\n",
            components.len()
        );
        println!("comp\trank\tread_start\tread_end\tchrom\tstrand\tref_start\tref_end\tatoms");
        for (ci, comp) in components.iter().enumerate() {
            for (ri, &idx) in comp.iter().enumerate() {
                let cs = &compounds[idx];
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    ci,
                    ri,
                    cs.read_start(),
                    cs.read_end(k),
                    chrom_name(cs.chrom_id()),
                    if cs.is_reverse() { "-" } else { "+" },
                    cs.ref_start(),
                    cs.ref_end(k),
                    cs.atoms().len()
                );
            }
        }

        let mut chains = find_all_chains(&compounds, k, &scheme);
        chains.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

        let sigma = 50.0_f64;
        let support = diagonal_support(&compounds, &chains, sigma);

        // Chain summary table.
        println!("\n=== Segment chains: {} ===", chains.len());
        println!(
            "chain\tscore\tdiag_support\tseeds\tchrom\tstrand\tread_start\tread_end\tread_span\tref_start\tref_end\tref_span"
        );
        for (ci, chain) in chains.iter().enumerate() {
            let first = &compounds[chain.chain[0]];
            let last = &compounds[*chain.chain.last().unwrap()];
            println!(
                "{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                ci,
                chain.score,
                support[ci],
                chain.chain.len(),
                chrom_name(first.chrom_id()),
                if first.is_reverse() { "-" } else { "+" },
                first.read_start(),
                last.read_end(k),
                last.read_end(k) - first.read_start(),
                first.ref_start(),
                last.ref_end(k),
                first.ref_end(k).max(last.ref_end(k)) - first.ref_start().max(last.ref_start())
            );
        }

        // Per-chain seed detail.
        println!(
            "\nchain\trank\tidx\tchrom\tstrand\tread_start\tread_end\tref_start\tref_end\tatoms\tweight\tscore"
        );
        for (ci, chain) in chains.iter().enumerate() {
            for (rank, &idx) in chain.chain.iter().enumerate() {
                let cs = &compounds[idx];
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.2}",
                    ci,
                    rank,
                    idx,
                    chrom_name(cs.chrom_id()),
                    if cs.is_reverse() { "-" } else { "+" },
                    cs.read_start(),
                    cs.read_end(k),
                    cs.ref_start(),
                    cs.ref_end(k),
                    cs.atoms().len(),
                    cs.weight(),
                    chain.score
                );
            }
        }

        // Alignment DP.
        let diag_weight = 0.1;
        let gap_cost_per_base = 0.02;
        let overlap_tolerance = 10;
        if let Some(result) = best_alignment(
            &compounds,
            &chains,
            &support,
            k,
            diag_weight,
            gap_cost_per_base,
            overlap_tolerance,
        ) {
            println!(
                "\n=== Alignment (score={:.3}, {} segments) ===",
                result.score,
                result.segments.len()
            );
            println!(
                "seg\tchain\tscore\tdiag_support\tchrom\tstrand\tread_start\tread_end\tref_start\tref_end"
            );
            for (si, &ci) in result.segments.iter().enumerate() {
                let chain = &chains[ci];
                let first = &compounds[chain.chain[0]];
                let last = &compounds[*chain.chain.last().unwrap()];
                println!(
                    "{}\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}",
                    si,
                    ci,
                    chain.score,
                    support[ci],
                    chrom_name(first.chrom_id()),
                    if first.is_reverse() { "-" } else { "+" },
                    first.read_start(),
                    last.read_end(k),
                    first.ref_start(),
                    last.ref_end(k),
                );
            }
        }

        let used_compounds: std::collections::HashSet<usize> =
            components.iter().flat_map(|c| c.iter().copied()).collect();
        let compounds = compounds
            .into_iter()
            .enumerate()
            .filter(|(idx, _seed)| used_compounds.contains(idx))
            .map(|(_idx, seed)| seed)
            .collect::<Vec<_>>();

        // Diagonal bands.
        let bands = partition_by_diagonal(
            &compounds,
            k,
            scheme.cfg.max_ref_deviation,
            read_length as usize,
        );
        println!("\n=== Diagonal bands: {} ===", bands.len());
        println!("band\tchrom\tstrand\tref_min\tref_max\tdiagonal\tmembers\tweight_sum\tcoverage");
        for (bi, band) in bands.iter().enumerate() {
            let weight_sum: f64 = band.members.iter().map(|&i| compounds[i].weight()).sum();
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{:.3}\t{:.3}",
                bi,
                chrom_name(band.chrom_id),
                if band.is_reverse { "-" } else { "+" },
                band.ref_min,
                band.ref_max,
                band.central_diagonal,
                band.members.len(),
                weight_sum,
                band.coverage
            );
        }

        let mut bands = bands;
        bands.retain(|b| b.coverage >= 0.33);

        // Per-band chaining with Gaussian-weighted scheme.
        println!("\n=== Banded chains (coverage >= 33%) ===");
        let print_sam = true;
        if print_sam {
            println!("@SQ\tSN:chr1\tLN:248956422");
            println!("@SQ\tSN:chr2\tLN:242193529");
            println!("@SQ\tSN:chr3\tLN:198295559");
            println!("@SQ\tSN:chr4\tLN:190214555");
            println!("@SQ\tSN:chr5\tLN:181538259");
            println!("@SQ\tSN:chr6\tLN:170805979");
            println!("@SQ\tSN:chr7\tLN:159345973");
            println!("@SQ\tSN:chr8\tLN:145138636");
            println!("@SQ\tSN:chr9\tLN:138394717");
            println!("@SQ\tSN:chr10\tLN:133797422");
            println!("@SQ\tSN:chr11\tLN:135086622");
            println!("@SQ\tSN:chr12\tLN:133275309");
            println!("@SQ\tSN:chr13\tLN:114364328");
            println!("@SQ\tSN:chr14\tLN:107043718");
            println!("@SQ\tSN:chr15\tLN:101991189");
            println!("@SQ\tSN:chr16\tLN:90338345");
            println!("@SQ\tSN:chr17\tLN:83257441");
            println!("@SQ\tSN:chr18\tLN:80373285");
            println!("@SQ\tSN:chr19\tLN:58617616");
            println!("@SQ\tSN:chr20\tLN:64444167");
            println!("@SQ\tSN:chr21\tLN:46709983");
            println!("@SQ\tSN:chr22\tLN:50818468");
            println!("@SQ\tSN:chrX\tLN:156040895");
            println!("@SQ\tSN:chrY\tLN:57227415");
            println!("@SQ\tSN:chrM\tLN:16569");
        } else {
            println!(
                "band_rank\tband\tchain\tseed\tchrom\tstrand\tread_start\tread_end\tref_start\tref_end\tweight\tedge_penalty\tread_gap\tref_gap"
            );
        }
        for (band_rank, band) in bands.iter().enumerate() {
            let band_label = format!(
                "{}:{}:{}",
                chrom_name(band.chrom_id),
                if band.is_reverse { "-" } else { "+" },
                band.central_diagonal.floor() as i64
            );
            let band_scheme = BandedSegmentScheme::new(
                SegmentConfig {
                    read_gap_cost_per_base: 0.25,
                    ..SegmentConfig::default_for_k(k)
                },
                band.chrom_id,
                band.is_reverse,
                band.central_diagonal,
                band.diagonal_variance,
                15.0,   // sv_break_penalty
                1.0,
                read_length as i64,
            );
            let mut band_chains = find_all_chains(&compounds, k, &band_scheme);
            band_chains.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

            let mut qual = String::new();
            for (ci, chain) in band_chains.iter().enumerate() {
                for (pos, &idx) in chain.chain.iter().enumerate() {
                    let seed = &compounds[idx];
                    let (next_edge, read_gap_str, ref_gap_str) = if pos + 1 < chain.chain.len() {
                        let next = &compounds[chain.chain[pos + 1]];
                        let (rg_str, rfg_str) = match seed.gap_to(next, k) {
                            Some((gap, _, _)) => (
                                format!("{}", gap.read_gap),
                                if gap.ref_gap == i64::MIN {
                                    "NA".to_string()
                                } else {
                                    format!("{}", gap.ref_gap)
                                },
                            ),
                            None => ("NA".to_string(), "NA".to_string()),
                        };
                        let edge_str = match band_scheme.edge_cost(seed, next, k) {
                            Some((c, _)) => format!("{:.3}", c),
                            None => "NA".to_string(),
                        };
                        (edge_str, rg_str, rfg_str)
                    } else {
                        ("NA".to_string(), "NA".to_string(), "NA".to_string())
                    };
                    if print_sam {
                        let flag = if seed.is_reverse() { 0x10 } else { 0x00 } | 0x800;
                        let chrom = chrom_name(seed.chrom_id());
                        let mut cigar_parts = Vec::new();
                        if seed.read_start() > 0 {
                            cigar_parts.push(format!("{}H", seed.read_start()));
                        }
                        cigar_parts.push(format!("{}=", seed.length(k)));
                        let right_clip = read_length - seed.read_end(k) as usize;
                        if right_clip > 0 {
                            cigar_parts.push(format!("{}H", right_clip));
                        }
                        let cigar = cigar_parts.join("");
                        let seq = seed.to_string(k);
                        while qual.len() < seq.len() {
                            qual.push('I');
                        }
                        println!(
                            "{read_name}\t{flag}\t{chrom}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
                            seed.ref_start() + 1,
                            60,
                            cigar,
                            seq,
                            &qual[0..seq.len()]
                        );
                    } else {
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}",
                            band_rank,
                            band_label,
                            ci,
                            pos,
                            chrom_name(seed.chrom_id()),
                            if seed.is_reverse() { "-" } else { "+" },
                            seed.read_start(),
                            seed.read_end(k),
                            seed.ref_start(),
                            seed.ref_end(k),
                            seed.weight(),
                            next_edge,
                            read_gap_str,
                            ref_gap_str,
                        );
                    }
                }
                break;
            }
        }
    }

    #[test]
    fn two_disjoint_components() {
        let k = 20;
        // Two pairs, each pair reachable within 2k, the pairs separated by > 2k.
        // Use distinct diagonals within each pair to avoid merging.
        let collection = SeedCollection::new(
            k,
            vec![
                make_seed(0, 0, 0, k),     // diagonal 0
                make_seed(25, 25, 0, k),   // diagonal 0, gap 5 → same component
                make_seed(200, 210, 0, k), // diagonal 10 — gap from pair 1 end (45) = 155 > 2k=40
                make_seed(225, 235, 0, k), // diagonal 10, gap 5 → same component as above
            ],
        );
        let compounds = collection.compound_seeds();
        let scheme = FullSegmentScheme::new(SegmentConfig::default_for_k(k));
        let components = partition_seeds(&compounds, k, &scheme);
        assert_eq!(components.len(), 2);
        for c in &components {
            assert_eq!(c.len(), 2);
        }
    }

    fn make_rev_seed(read_pos: u32, ref_pos: u32, chrom_id: u32, k: usize) -> AtomicSeed {
        AtomicSeed::new(10000 - read_pos - k as u32, 10000, k, chrom_id, ref_pos, true, 0, 1)
    }

    #[test]
    fn bits_overlapping_seeds_not_double_counted() {
        let k = 20;
        let mut atoms = vec![];
        // Main locus: read 0..60 and 100..160 on one diagonal.
        for r in (0..60).step_by(20).chain((100..160).step_by(20)) {
            atoms.push(make_seed(r, 1000 + r, 0, k));
        }
        // Another chromosome covering read 50..110, overlapping both.  Credited
        // in full it would be worth a detour; credited only for the read bases
        // it adds, it is not.
        for r in (50..110).step_by(20) {
            atoms.push(make_seed(r, 5000 + r, 1, k));
        }
        let collection = SeedCollection::new(k, atoms);
        let compounds = collection.compound_seeds();
        let cfg = ChainingConfig::builder().load().unwrap();
        let chain = chain_bits(&compounds, k, &cfg, |_, _, _, _, _, _| None).unwrap();

        let chroms: Vec<u32> = chain.chain.iter().map(|&i| compounds[i].chrom_id()).collect();
        assert_eq!(chroms, vec![0, 0]);
        assert_eq!(chain.edge_type, vec![EdgeType::Continuation]);
    }

    #[test]
    fn bits_prefers_nearby_inverted_segment() {
        let k = 20;
        let mut atoms = vec![];
        // Main locus, forward: read 0..200, then read 320..520 after a 40bp
        // net deletion.
        for r in (0..200).step_by(20) {
            atoms.push(make_seed(r, 10000 + r, 0, k));
        }
        for r in (320..520).step_by(20) {
            atoms.push(make_seed(r, 10040 + r, 0, k));
        }
        // Read 200..320 matches exactly on the reverse strand ~1kb away, and
        // equally well on another chromosome.
        for r in (200..320).step_by(20) {
            atoms.push(make_rev_seed(r, 9300 - r, 0, k));
            atoms.push(make_rev_seed(r, 50300 - r, 1, k));
        }
        let collection = SeedCollection::new(k, atoms);
        let compounds = collection.compound_seeds();
        let cfg = ChainingConfig::builder().load().unwrap();
        let chain = chain_bits(&compounds, k, &cfg, |_, _, _, _, _, _| None).unwrap();

        let loci: Vec<(u32, bool)> = chain
            .chain
            .iter()
            .map(|&i| (compounds[i].chrom_id(), compounds[i].is_reverse()))
            .collect();
        assert_eq!(loci, vec![(0, false), (0, true), (0, false)]);
        assert_eq!(chain.edge_type, vec![EdgeType::SvBreak, EdgeType::SvBreak]);
    }

    #[test]
    fn bits_tandem_copies_emitted_inline() {
        let k = 20;
        let mut atoms = vec![];
        // Read 0..100 matches ref 1000..1100; read 100..140 and 140..180 are
        // two extra copies of ref 1060..1100; read 190..290 carries on at ref
        // 1110..1210.
        for r in (0..100).step_by(20) {
            atoms.push(make_seed(r, 1000 + r, 0, k));
        }
        for r in (100..180).step_by(20) {
            atoms.push(make_seed(r, 1060 + (r - 100) % 40, 0, k));
        }
        for r in (190..290).step_by(20) {
            atoms.push(make_seed(r, 920 + r, 0, k));
        }
        let collection = SeedCollection::new(k, atoms);
        let compounds = collection.compound_seeds();
        let cfg = ChainingConfig::builder().load().unwrap();
        let chain = chain_bits(&compounds, k, &cfg, |_, _, _, _, _, _| None).unwrap();

        // The copies earn their coverage in the score (more than the 200 bases
        // of flanks alone could), but are emitted as an insertion between the
        // flanking seeds rather than as separate records.
        assert!(chain.score > cfg.match_bits * 200.0);
        let reads: Vec<(u32, u32)> = chain
            .chain
            .iter()
            .map(|&i| (compounds[i].read_start(), compounds[i].read_end(k)))
            .collect();
        assert_eq!(reads, vec![(0, 100), (190, 290)]);
        assert_eq!(chain.edge_type, vec![EdgeType::Continuation]);
    }
}
