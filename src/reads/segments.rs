use super::compound::{GapComputable, Seed, Weighted};
use parallax::utils::union_find::UnionFind;

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
    ) -> Option<f64>;
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
    order.sort_unstable_by_key(|&i| seeds[i].read_pos());

    let uf = UnionFind::new(n);

    for rank in 0..n {
        let i = order[rank];
        let i_read_end = seeds[i].read_end(k) as i64;

        for rank2 in rank + 1..n {
            let j = order[rank2];
            let read_gap = seeds[j].read_pos() as i64 - i_read_end;
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
            let first = seeds[c[0]].read_pos() as i64;
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
    /// Cost per base of ref-deviation (|ref_gap - read_gap|).
    pub ref_dev_cost_per_base: f64,
}

impl SegmentConfig {
    pub fn default_for_k(k: usize) -> Self {
        Self {
            max_read_gap: 2 * k as i64,
            max_ref_deviation: 200,
            min_segment_span: 2 * k as i64,
            read_gap_cost_per_base: 0.05,
            ref_dev_cost_per_base: 0.01,
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
        let read_gap = rhs.read_pos() as i64 - lhs_read_end;
        if read_gap > self.cfg.max_read_gap {
            return false;
        }
        let ref_gap = if lhs.is_reverse() {
            lhs.ref_pos() as i64 - (rhs.ref_end(k) as i64)
        } else {
            rhs.ref_pos() as i64 - (lhs.ref_end(k) as i64)
        };
        let deviation = (ref_gap - read_gap).abs();
        deviation <= self.cfg.max_ref_deviation
    }

    fn edge_cost<S: Seed + Weighted + GapComputable>(
        &self,
        lhs: &S,
        rhs: &S,
        k: usize,
    ) -> Option<f64> {
        let gap = lhs.gap_to(rhs, k)?;
        // Cross-chrom/strand sentinel — not connectable at level 2.
        if gap.ref_gap == i64::MIN {
            return None;
        }
        let read_gap_cost = (gap.read_gap.max(0) as f64) * self.cfg.read_gap_cost_per_base;
        let deviation = (gap.ref_gap - gap.read_gap).unsigned_abs() as f64;
        let ref_dev_cost = deviation * self.cfg.ref_dev_cost_per_base;
        Some(read_gap_cost + ref_dev_cost + gap.weight_trimmed)
    }
}

// ── Segment chain result ──────────────────────────────────────────────────────

/// A single chain extracted from a connected component of compound seeds.
pub struct SegmentChain {
    /// DP score of this chain.
    pub score: f64,
    /// Indices into the original seed slice, in read-position order.
    pub chain: Vec<usize>,
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

        for rank in 0..n {
            let i = component[active[rank]];
            dp[rank] = seeds[i].weight();

            for r in (0..rank).rev() {
                let j = component[active[r]];
                // Early exit: if the read gap is already beyond max_read_gap,
                // no earlier seed can connect either.
                let read_gap = seeds[i].read_pos() as i64 - (seeds[j].read_pos() as i64 + k as i64);
                if read_gap > scheme.max_read_gap() {
                    break;
                }

                if let Some(cost) = scheme.edge_cost(&seeds[j], &seeds[i], k) {
                    let candidate = dp[r] + seeds[i].weight() - cost;
                    if candidate > dp[rank] {
                        dp[rank] = candidate;
                        prev[rank] = r;
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
                score: seeds[idx].weight(),
                chain: vec![idx],
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
pub fn diagonal_support<S>(
    seeds: &[S],
    chains: &[SegmentChain],
    sigma: f64,
) -> Vec<f64>
where
    S: Seed,
{
    let cutoff = 4.0 * sigma;
    let two_sigma_sq = 2.0 * sigma * sigma;
    let n = chains.len();

    // For each chain, extract the diagonal of its first seed plus chrom/strand.
    let diag: Vec<(u32, bool, i64)> = chains.iter().map(|c| {
        let s = &seeds[c.chain[0]];
        (s.chrom_id(), s.is_reverse(), s.diagonal())
    }).collect();

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
    let read_start = |i: usize| seeds[chains[i].chain[0]].read_pos() as i64;
    let read_end   = |i: usize| seeds[*chains[i].chain.last().unwrap()].read_end(k) as i64;

    // Sort by read_end, breaking ties by read_start descending.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| (read_end(i), -(read_start(i))));

    let node_weight = |i: usize| chains[i].score + diag_weight * support[i];

    let inf = f64::NEG_INFINITY;
    let mut dp   = vec![inf; n];
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
        if p == usize::MAX { break; }
        cur = p;
    }
    segments.reverse();

    Some(AlignmentResult { score: dp[best_rank], segments })
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
        let band_set: std::collections::HashSet<usize> =
            band.members.iter().copied().collect();

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
        let band_support: Vec<f64> =
            band_chain_indices.iter().map(|&ci| support[ci]).collect();

        // Build owned SegmentChain vec for best_alignment (needs &[SegmentChain]).
        let owned_chains: Vec<SegmentChain> = band_chains
            .iter()
            .map(|c| SegmentChain { score: c.score, chain: c.chain.clone() })
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
        let rs = result.segments.iter()
            .map(|&si| seeds[owned_chains[si].chain[0]].read_pos())
            .min()
            .unwrap();
        let re = result.segments.iter()
            .map(|&si| seeds[*owned_chains[si].chain.last().unwrap()].read_end(k) as u32)
            .max()
            .unwrap();

        // Remap segment indices back to global chain indices.
        let remapped = AlignmentResult {
            score: result.score,
            segments: result.segments.iter()
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
        result.segments.iter().map(|&ci| {
            let chain = &chains[ci];
            let rs = seeds[chain.chain[0]].read_pos();
            let re = seeds[*chain.chain.last().unwrap()].read_end(k) as u32;
            (rs, re)
        }).collect()
    };

    let overlaps_any = |ivs: &[(u32, u32)], occupied: &[(u32, u32)]| -> bool {
        ivs.iter().any(|&(s, e)| occupied.iter().any(|&(os, oe)| s < oe && e > os))
    };

    // Outer loop: iterate over bands in descending span order.
    // Halt when the current band covers less than 50% of read length.
    for outer in 0..band_results.len() {
        if used[outer] { continue; }
        let span = band_results[outer].3 - band_results[outer].2;
        if span < half_read { break; }

        // Accept the outer band and mark its read intervals as occupied.
        let ivs = intervals_of(&band_results[outer].1);
        occupied.extend_from_slice(&ivs);
        used[outer] = true;
        let (bi, ref result, rs, re) = band_results[outer];
        accepted.push(BandAlignment {
            band_idx: bi,
            result: AlignmentResult { score: result.score, segments: result.segments.clone() },
            read_start: rs,
            read_end: re,
        });

        // Inner loop: scan remaining bands for non-overlapping gap-fillers.
        for inner in (outer + 1)..band_results.len() {
            if used[inner] { continue; }
            let ivs = intervals_of(&band_results[inner].1);
            if overlaps_any(&ivs, &occupied) { continue; }

            occupied.extend_from_slice(&ivs);
            used[inner] = true;
            let (bi, ref result, rs, re) = band_results[inner];
            accepted.push(BandAlignment {
                band_idx: bi,
                result: AlignmentResult { score: result.score, segments: result.segments.clone() },
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
    /// Indices into the seed slice.
    pub members: Vec<usize>,
}

/// Partition `seeds` into diagonal bands by sorting on `(chrom_id, is_reverse,
/// diagonal)` and grouping consecutive seeds whose diagonal differs by less
/// than `max_deviation` from the first seed of the current group.
pub fn partition_by_diagonal<S>(seeds: &[S], k: usize, max_deviation: i64) -> Vec<DiagonalBand>
where
    S: Seed + Weighted,
{
    let n = seeds.len();
    if n == 0 {
        return Vec::new();
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| {
        (seeds[i].chrom_id(), seeds[i].is_reverse(), seeds[i].diagonal())
    });

    let mut bands: Vec<DiagonalBand> = Vec::new();
    let mut group_start = 0usize;

    while group_start < n {
        let first = order[group_start];
        let chrom = seeds[first].chrom_id();
        let rev   = seeds[first].is_reverse();
        let d0    = seeds[first].diagonal();

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
        let band = diagonal_band_stats(seeds, &members, k, chrom, rev);
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
) -> DiagonalBand
where
    S: Seed + Weighted,
{
    let mut ref_min = u32::MAX;
    let mut ref_max = 0u32;
    let mut weight_sum = 0.0f64;
    let mut diag_sum   = 0.0f64;

    for &i in members {
        let s = &seeds[i];
        let w = s.weight();
        ref_min = ref_min.min(s.ref_start());
        ref_max = ref_max.max(s.ref_end(k));
        weight_sum += w;
        diag_sum   += w * s.diagonal() as f64;
    }

    let central_diagonal = if weight_sum > 0.0 { diag_sum / weight_sum } else { 0.0 };

    DiagonalBand {
        chrom_id,
        is_reverse,
        ref_min,
        ref_max,
        central_diagonal,
        members: members.to_vec(),
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::reads::compound::{AtomicSeed, SeedCollection};

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

        let collection = SeedCollection::new(k, atoms);
        let compounds = collection.compound_seeds();
        println!("=== Compound seeds: {} ===", compounds.len());
        println!(
            "rank\tread_start\tread_end\tchrom\tstrand\tref_start\tref_end\tatoms\tread_span"
        );
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
        println!(
            "comp\trank\tread_start\tread_end\tchrom\tstrand\tref_start\tref_end\tatoms"
        );
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
            &compounds, &chains, &support, k,
            diag_weight, gap_cost_per_base, overlap_tolerance,
        ) {
            println!("\n=== Alignment (score={:.3}, {} segments) ===", result.score, result.segments.len());
            println!("seg\tchain\tscore\tdiag_support\tchrom\tstrand\tread_start\tread_end\tref_start\tref_end");
            for (si, &ci) in result.segments.iter().enumerate() {
                let chain = &chains[ci];
                let first = &compounds[chain.chain[0]];
                let last  = &compounds[*chain.chain.last().unwrap()];
                println!(
                    "{}\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}",
                    si, ci,
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

        // Diagonal bands.
        let bands = partition_by_diagonal(&compounds, k, scheme.cfg.max_ref_deviation);
        println!("\n=== Diagonal bands: {} ===", bands.len());
        println!("band\tchrom\tstrand\tref_min\tref_max\tcentral_diagonal\tmembers\tweight_sum");
        for (bi, band) in bands.iter().enumerate() {
            let weight_sum: f64 = band.members.iter().map(|&i| compounds[i].weight()).sum();
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{:.3}",
                bi,
                chrom_name(band.chrom_id),
                if band.is_reverse { "-" } else { "+" },
                band.ref_min,
                band.ref_max,
                band.central_diagonal,
                band.members.len(),
                weight_sum
            );
        }

        // Banded assembly.
        let read_length = rows.iter().map(|r| r.read_pos + k as u32).max().unwrap_or(0);
        let assembly = assemble_alignment(
            &compounds, &chains, &support, &bands,
            k, read_length, diag_weight, gap_cost_per_base,
        );
        println!("\n=== Banded assembly: {} bands accepted ===", assembly.len());
        println!("order\tband\tchrom\tstrand\tread_start\tread_end\tspan\tscore\tsegments");
        for (order, ba) in assembly.iter().enumerate() {
            let band = &bands[ba.band_idx];
            let first_chain = &chains[ba.result.segments[0]];
            let first_seed = &compounds[first_chain.chain[0]];
            let total_span = ba.read_end - ba.read_start;
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
                order,
                ba.band_idx,
                chrom_name(band.chrom_id),
                if band.is_reverse { "-" } else { "+" },
                ba.read_start,
                ba.read_end,
                total_span,
                ba.result.score,
                ba.result.segments.len(),
            );
            println!("  seg\tchain\tchrom\tstrand\tread_start\tread_end\tref_start\tref_end\tscore\tdiag_support");
            for &ci in &ba.result.segments {
                let chain = &chains[ci];
                let fs = &compounds[chain.chain[0]];
                let ls = &compounds[*chain.chain.last().unwrap()];
                println!(
                    "  {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}",
                    ci, ci,
                    chrom_name(fs.chrom_id()),
                    if fs.is_reverse() { "-" } else { "+" },
                    fs.read_start(),
                    ls.read_end(k),
                    fs.ref_start(),
                    ls.ref_end(k),
                    chain.score,
                    support[ci],
                );
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
}
