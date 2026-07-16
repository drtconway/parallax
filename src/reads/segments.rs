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

#[cfg(test)]
mod tests {
    use std::cmp::Reverse;

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

        // Chain summary table.
        println!("\n=== Segment chains: {} ===", chains.len());
        println!(
            "chain\tscore\tseeds\tchrom\tstrand\tread_start\tread_end\tread_span\tref_start\tref_end\tref_span"
        );
        for (ci, chain) in chains.iter().enumerate() {
            let first = &compounds[chain.chain[0]];
            let last = &compounds[*chain.chain.last().unwrap()];
            println!(
                "{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                ci,
                chain.score,
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
