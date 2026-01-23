// =============================================================================
// Minimap2-style RMQ chaining
// =============================================================================

/// A seed anchor for minimap2-style chaining.
/// Must have reference position, query position, and a score contribution.
pub trait ChainAnchor {
    /// Reference position (x coordinate)
    fn ref_pos(&self) -> i64;
    /// Query position (y coordinate)  
    fn query_pos(&self) -> i64;
    /// Diagonal = ref_pos - query_pos
    fn diagonal(&self) -> i64 {
        self.ref_pos() - self.query_pos()
    }
    /// Score contribution of this anchor (e.g., seed length)
    fn weight(&self) -> f64;
}

/// A lightweight proxy anchor used internally for anchor-first chaining.
/// This avoids requiring the Copy trait on user-provided anchor types.
#[derive(Clone, Copy)]
struct IndexProxyAnchor {
    ref_pos: i64,
    query_pos: i64,
    weight: f64,
}

impl ChainAnchor for IndexProxyAnchor {
    fn ref_pos(&self) -> i64 {
        self.ref_pos
    }
    fn query_pos(&self) -> i64 {
        self.query_pos
    }
    fn weight(&self) -> f64 {
        self.weight
    }
}

/// Parameters for minimap2-style chaining.
#[derive(Clone, Debug)]
pub struct ChainParams {
    /// Maximum gap in reference allowed between consecutive anchors
    pub max_gap_ref: i64,
    /// Maximum gap in query allowed between consecutive anchors
    pub max_gap_query: i64,
    /// Linear gap penalty coefficient (α in minimap2)
    pub gap_penalty_linear: f64,
    /// Logarithmic gap penalty coefficient (β in minimap2)
    pub gap_penalty_log: f64,
    /// Linear penalty for diagonal jumps (|ref_gap - query_gap|)
    pub diagonal_penalty_linear: f64,
    /// Bandwidth: max diagonal difference between chained anchors
    pub bandwidth: i64,
}

impl Default for ChainParams {
    fn default() -> Self {
        Self {
            max_gap_ref: 10000,
            max_gap_query: 10000,
            gap_penalty_linear: 0.01,
            gap_penalty_log: 0.5,
            diagonal_penalty_linear: 0.5,
            bandwidth: 500,
        }
    }
}

/// Result of chaining: the best chain and its score.
#[derive(Clone, Debug)]
pub struct ChainResult {
    /// Indices of anchors in the chain (in order)
    pub chain: Vec<usize>,
}

/// Minimap2-style RMQ chainer.
///
/// Uses a segment tree indexed by diagonal for O(log n) predecessor queries.
/// The algorithm:
/// 1. Sort anchors by reference position
/// 2. For each anchor, find the best predecessor within bandwidth
/// 3. Use a segment tree for RMQ on (diagonal -> best score ending at that diagonal)
pub struct RmqChainer {
    /// DP scores for each anchor
    scores: Vec<f64>,
    /// Predecessor index for traceback
    predecessors: Vec<usize>,
    /// Segment tree for RMQ queries
    /// Maps discretized diagonal -> (best_score, anchor_index)
    tree: Vec<(f64, usize)>,
    /// Diagonal offset for indexing
    diag_offset: i64,
    /// Number of diagonal buckets
    n_diags: usize,
}

impl Default for RmqChainer {
    fn default() -> Self {
        Self::new()
    }
}

impl RmqChainer {
    pub fn new() -> Self {
        Self {
            scores: Vec::new(),
            predecessors: Vec::new(),
            tree: Vec::new(),
            diag_offset: 0,
            n_diags: 0,
        }
    }

    /// Chain anchors using minimap2-style DP with RMQ optimization.
    ///
    /// Anchors should implement the ChainAnchor trait.
    /// Returns the best chain found.
    ///
    /// Time complexity: O(n log n) where n is the number of anchors.
    pub fn chain<T: ChainAnchor>(&mut self, anchors: &[T], params: &ChainParams) -> ChainResult {
        let n = anchors.len();
        if n == 0 {
            return ChainResult { chain: Vec::new() };
        }

        // Sort indices by reference position
        let mut sorted_indices: Vec<usize> = (0..n).collect();
        sorted_indices.sort_unstable_by_key(|&i| anchors[i].ref_pos());

        // Find diagonal range for segment tree sizing
        let min_diag = anchors.iter().map(|a| a.diagonal()).min().unwrap();
        let max_diag = anchors.iter().map(|a| a.diagonal()).max().unwrap();
        self.diag_offset = min_diag;
        self.n_diags = (max_diag - min_diag + 1) as usize;

        // Initialize segment tree (size = 2 * next power of 2)
        let tree_size = self.n_diags.next_power_of_two() * 2;
        self.tree.clear();
        self.tree.resize(tree_size, (f64::NEG_INFINITY, usize::MAX));

        // Initialize DP arrays
        self.scores.clear();
        self.scores.resize(n, 0.0);
        self.predecessors.clear();
        self.predecessors.resize(n, usize::MAX);

        let mut best_score = f64::NEG_INFINITY;
        let mut best_end = 0usize;

        // Process anchors in reference position order
        for &i in &sorted_indices {
            let anchor = &anchors[i];
            let ref_pos = anchor.ref_pos();
            let query_pos = anchor.query_pos();
            let diag = anchor.diagonal();
            let weight = anchor.weight();

            // Query range: diagonals within bandwidth
            let diag_lo = diag - params.bandwidth;
            let diag_hi = diag + params.bandwidth;

            // Find best predecessor using RMQ
            let (pred_score, pred_idx) = self.query_range(diag_lo, diag_hi);

            // Compute score if we chain from predecessor
            let mut score = weight;
            if pred_idx != usize::MAX {
                let pred = &anchors[pred_idx];
                let gap_ref = ref_pos - pred.ref_pos();
                let gap_query = query_pos - pred.query_pos();

                // Check gap constraints
                if gap_ref > 0
                    && gap_query > 0
                    && gap_ref <= params.max_gap_ref
                    && gap_query <= params.max_gap_query
                {
                    // Compute gap penalty (minimap2 style)
                    let gap_cost = self.gap_penalty(gap_ref, gap_query, params);
                    let chain_score = pred_score + weight - gap_cost;

                    if chain_score > weight {
                        score = chain_score;
                        self.predecessors[i] = pred_idx;
                    }
                }
            }

            self.scores[i] = score;

            // Update segment tree
            let diag_idx = (diag - self.diag_offset) as usize;
            self.update(diag_idx, score, i);

            // Track best ending position
            if score > best_score {
                best_score = score;
                best_end = i;
            }
        }

        // Traceback to build chain
        let mut chain = Vec::new();
        let mut idx = best_end;
        while idx != usize::MAX {
            chain.push(idx);
            idx = self.predecessors[idx];
        }
        chain.reverse();

        ChainResult { chain }
    }

    /// Compute gap penalty using minimap2's log-linear model.
    fn gap_penalty(&self, gap_ref: i64, gap_query: i64, params: &ChainParams) -> f64 {
        let gap = gap_ref.max(gap_query) as f64;
        let gap_diff = (gap_ref - gap_query).abs() as f64;

        // Linear component + log component for long gaps + diagonal deviation penalty
        params.gap_penalty_linear * gap
            + params.gap_penalty_log * (gap + 1.0).ln()
            + 0.5 * gap_diff.ln().max(0.0)
            + params.diagonal_penalty_linear * gap_diff
    }

    /// Query the segment tree for the maximum score in diagonal range [lo, hi].
    fn query_range(&self, diag_lo: i64, diag_hi: i64) -> (f64, usize) {
        let lo = ((diag_lo - self.diag_offset).max(0) as usize).min(self.n_diags.saturating_sub(1));
        let hi = ((diag_hi - self.diag_offset).max(0) as usize).min(self.n_diags.saturating_sub(1));

        if lo > hi || self.tree.is_empty() {
            return (f64::NEG_INFINITY, usize::MAX);
        }

        let tree_n = self.tree.len() / 2;
        let mut lo = lo + tree_n;
        let mut hi = hi + tree_n + 1;

        let mut best = (f64::NEG_INFINITY, usize::MAX);

        while lo < hi {
            if lo & 1 == 1 {
                if self.tree[lo].0 > best.0 {
                    best = self.tree[lo];
                }
                lo += 1;
            }
            if hi & 1 == 1 {
                hi -= 1;
                if self.tree[hi].0 > best.0 {
                    best = self.tree[hi];
                }
            }
            lo >>= 1;
            hi >>= 1;
        }

        best
    }

    /// Update the segment tree at diagonal index with a new score.
    fn update(&mut self, diag_idx: usize, score: f64, anchor_idx: usize) {
        if self.tree.is_empty() {
            return;
        }
        let tree_n = self.tree.len() / 2;
        let mut i = diag_idx + tree_n;

        if i >= self.tree.len() {
            return;
        }

        // Only update if this score is better
        if score > self.tree[i].0 {
            self.tree[i] = (score, anchor_idx);

            // Propagate up
            while i > 1 {
                i >>= 1;
                let left = self.tree[i * 2];
                let right = self.tree[i * 2 + 1];
                self.tree[i] = if left.0 >= right.0 { left } else { right };
            }
        }
    }

    /// Anchor-first chaining: chain long seeds first, then fill gaps with short seeds.
    ///
    /// This approach addresses the problem where short off-diagonal seeds can
    /// "hijack" the chain away from the correct diagonal. By chaining long seeds
    /// first (which are more trustworthy), we establish the correct diagonal,
    /// then only accept short seeds that are consistent with it.
    ///
    /// # Algorithm
    /// 1. Partition seeds into anchors (length >= threshold) and fillers (shorter)
    /// 2. Chain anchors using tight diagonal tolerance
    /// 3. For each gap between consecutive anchors, add fillers that:
    ///    - Are positioned between the anchors (in ref and query space)
    ///    - Have a diagonal within tolerance of the interpolated anchor diagonal
    ///
    /// # Arguments
    /// * `anchors` - All seed anchors
    /// * `params` - Chaining parameters (bandwidth used for anchor chaining)
    /// * `anchor_threshold` - Minimum seed length to be considered an anchor
    /// * `filler_diagonal_tolerance` - Max diagonal deviation for fillers from expected
    ///
    /// # Returns
    /// Indices of seeds in the final chain (in reference position order)
    pub fn chain_anchor_first<T: ChainAnchor>(
        &mut self,
        seeds: &[T],
        params: &ChainParams,
        anchor_threshold: f64,
        filler_diagonal_tolerance: i64,
    ) -> ChainResult {
        if seeds.is_empty() {
            return ChainResult { chain: Vec::new() };
        }

        // Phase 1: Partition into anchors and fillers
        let anchor_indices: Vec<usize> = seeds
            .iter()
            .enumerate()
            .filter(|(_, s)| s.weight() >= anchor_threshold)
            .map(|(i, _)| i)
            .collect();

        let filler_indices: Vec<usize> = seeds
            .iter()
            .enumerate()
            .filter(|(_, s)| s.weight() < anchor_threshold)
            .map(|(i, _)| i)
            .collect();

        // If no anchors, fall back to regular chaining
        if anchor_indices.is_empty() {
            return self.chain(seeds, params);
        }

        // If only one anchor, use it as the chain and add compatible fillers
        if anchor_indices.len() == 1 {
            let anchor_idx = anchor_indices[0];
            let anchor = &seeds[anchor_idx];
            let anchor_diag = anchor.diagonal();

            // Collect fillers within diagonal tolerance
            let mut compatible_fillers: Vec<usize> = Vec::new();
            for &fi in &filler_indices {
                let filler = &seeds[fi];
                if (filler.diagonal() - anchor_diag).abs() <= filler_diagonal_tolerance {
                    compatible_fillers.push(fi);
                }
            }

            // Sort fillers by reference position
            compatible_fillers.sort_unstable_by_key(|&i| seeds[i].ref_pos());

            // Filter to ensure colinearity: query_pos must be monotonically increasing
            let anchor_ref_pos = anchor.ref_pos();
            let anchor_query_pos = anchor.query_pos();
            let anchor_weight = anchor.weight() as i64;

            // Split fillers into those before and after the anchor
            let (before_anchor, after_anchor): (Vec<_>, Vec<_>) = compatible_fillers
                .into_iter()
                .partition(|&fi| seeds[fi].ref_pos() < anchor_ref_pos);

            // Build chain: fillers before anchor (must have query_pos < anchor's)
            let mut chain: Vec<usize> = Vec::new();
            let mut last_query_end: i64 = i64::MIN;
            for &fi in &before_anchor {
                let filler = &seeds[fi];
                let filler_query_pos = filler.query_pos();
                let filler_query_end = filler_query_pos + filler.weight() as i64;
                // Must be colinear and end before anchor starts
                if filler_query_pos > last_query_end && filler_query_end <= anchor_query_pos {
                    chain.push(fi);
                    last_query_end = filler_query_end;
                }
            }

            // Add the anchor
            chain.push(anchor_idx);

            // Add fillers after anchor (must have query_pos > anchor's end)
            let mut last_query_end = anchor_query_pos + anchor_weight;
            for &fi in &after_anchor {
                let filler = &seeds[fi];
                let filler_query_pos = filler.query_pos();
                if filler_query_pos > last_query_end {
                    chain.push(fi);
                    last_query_end = filler_query_pos + filler.weight() as i64;
                }
            }

            return ChainResult { chain };
        }

        // Phase 2: Chain anchors only, using tighter parameters
        // Create index-proxy anchors that reference back to the original slice
        let anchor_proxies: Vec<IndexProxyAnchor> = anchor_indices
            .iter()
            .map(|&i| {
                let s = &seeds[i];
                IndexProxyAnchor {
                    ref_pos: s.ref_pos(),
                    query_pos: s.query_pos(),
                    weight: s.weight(),
                }
            })
            .collect();

        // Use tighter bandwidth for anchor chaining
        let anchor_params = ChainParams {
            bandwidth: params.bandwidth / 2, // Tighter diagonal tolerance
            diagonal_penalty_linear: params.diagonal_penalty_linear * 2.0, // Penalize diagonal jumps more
            ..params.clone()
        };

        let anchor_chain_result = self.chain(&anchor_proxies, &anchor_params);

        if anchor_chain_result.chain.is_empty() {
            // No valid anchor chain found, fall back to regular chaining
            return self.chain(seeds, params);
        }

        // Map anchor chain indices back to original seed indices
        let anchor_chain: Vec<usize> = anchor_chain_result
            .chain
            .iter()
            .map(|&i| anchor_indices[i])
            .collect();

        // Phase 3: Fill gaps between anchors with consistent fillers
        let mut final_chain: Vec<usize> = Vec::with_capacity(seeds.len());

        for window in anchor_chain.windows(2) {
            let (prev_idx, next_idx) = (window[0], window[1]);
            let prev_anchor = &seeds[prev_idx];
            let next_anchor = &seeds[next_idx];

            // Add the first anchor of this pair
            final_chain.push(prev_idx);

            // Get anchor positions and diagonals for interpolation
            let prev_diag = prev_anchor.diagonal();
            let next_diag = next_anchor.diagonal();

            // Find fillers that are:
            // 1. Between the two anchors in reference space
            // 2. Between the two anchors in query space
            // 3. Within diagonal tolerance of the *interpolated* expected diagonal
            let prev_ref_end = prev_anchor.ref_pos() + prev_anchor.weight() as i64;
            let next_ref_start = next_anchor.ref_pos();
            let prev_query_end = prev_anchor.query_pos() + prev_anchor.weight() as i64;
            let next_query_start = next_anchor.query_pos();

            // Reference span for interpolation
            let ref_span = next_ref_start - prev_ref_end;

            let mut gap_fillers: Vec<usize> = filler_indices
                .iter()
                .copied()
                .filter(|&fi| {
                    let filler = &seeds[fi];
                    let f_ref = filler.ref_pos();
                    let f_query = filler.query_pos();
                    let f_diag = filler.diagonal();

                    // Must be between anchors in both dimensions
                    if f_ref < prev_ref_end
                        || f_ref >= next_ref_start
                        || f_query < prev_query_end
                        || f_query >= next_query_start
                    {
                        return false;
                    }

                    // Interpolate expected diagonal based on filler's ref position
                    // This handles natural diagonal drift between anchors
                    let expected_diag = if ref_span > 0 {
                        let t = (f_ref - prev_ref_end) as f64 / ref_span as f64;
                        prev_diag as f64 + t * (next_diag - prev_diag) as f64
                    } else {
                        (prev_diag + next_diag) as f64 / 2.0
                    };

                    // Check if filler's diagonal is within tolerance of expected
                    (f_diag as f64 - expected_diag).abs() <= filler_diagonal_tolerance as f64
                })
                .collect();

            // Sort fillers by reference position
            gap_fillers.sort_unstable_by_key(|&i| seeds[i].ref_pos());

            // Filter to ensure colinearity: query_pos must be monotonically increasing
            // when sorted by ref_pos. This is a simple greedy approach that keeps
            // fillers that maintain monotonicity with respect to previous filler.
            let mut colinear_fillers: Vec<usize> = Vec::with_capacity(gap_fillers.len());
            let mut last_query_pos = prev_query_end;

            for &fi in &gap_fillers {
                let filler_query_pos = seeds[fi].query_pos();
                if filler_query_pos > last_query_pos {
                    colinear_fillers.push(fi);
                    last_query_pos = filler_query_pos + seeds[fi].weight() as i64;
                }
            }

            // Add colinear gap fillers
            final_chain.extend(colinear_fillers);
        }

        // Add the last anchor
        if let Some(&last_anchor_idx) = anchor_chain.last() {
            final_chain.push(last_anchor_idx);
        }

        // The chain should already be in ref_pos order by construction.
        // We built it by iterating through anchor_chain.windows(2) in order,
        // adding anchors and their gap fillers (sorted by ref_pos) sequentially.
        // 
        // DO NOT re-sort here - that could interleave fillers from different gaps
        // and break query_pos monotonicity (colinearity).

        // Remove duplicates (in case an anchor was added twice)
        // Use a stable approach that preserves order
        let mut seen = std::collections::HashSet::new();
        final_chain.retain(|&idx| seen.insert(idx));

        ChainResult { chain: final_chain }
    }
}

#[cfg(test)]
mod rmq_chainer_tests {
    use super::*;

    #[derive(Clone, Debug)]
    struct TestAnchor {
        ref_pos: i64,
        query_pos: i64,
        weight: f64,
    }

    impl ChainAnchor for TestAnchor {
        fn ref_pos(&self) -> i64 {
            self.ref_pos
        }
        fn query_pos(&self) -> i64 {
            self.query_pos
        }
        fn weight(&self) -> f64 {
            self.weight
        }
    }

    #[test]
    fn test_empty() {
        let mut chainer = RmqChainer::new();
        let anchors: Vec<TestAnchor> = vec![];
        let result = chainer.chain(&anchors, &ChainParams::default());
        assert!(result.chain.is_empty());
    }

    #[test]
    fn test_single_anchor() {
        let mut chainer = RmqChainer::new();
        let anchors = vec![TestAnchor {
            ref_pos: 100,
            query_pos: 50,
            weight: 20.0,
        }];
        let result = chainer.chain(&anchors, &ChainParams::default());
        assert_eq!(result.chain, vec![0]);
    }

    #[test]
    fn test_colinear_chain() {
        let mut chainer = RmqChainer::new();
        // Three anchors on the same diagonal, should chain together
        let anchors = vec![
            TestAnchor {
                ref_pos: 100,
                query_pos: 100,
                weight: 20.0,
            },
            TestAnchor {
                ref_pos: 200,
                query_pos: 200,
                weight: 20.0,
            },
            TestAnchor {
                ref_pos: 300,
                query_pos: 300,
                weight: 20.0,
            },
        ];
        let result = chainer.chain(&anchors, &ChainParams::default());
        assert_eq!(result.chain, vec![0, 1, 2]);
    }

    #[test]
    fn test_non_colinear_filtered() {
        let mut chainer = RmqChainer::new();
        // Anchor 1 has decreasing query pos - should not chain
        let anchors = vec![
            TestAnchor {
                ref_pos: 100,
                query_pos: 200,
                weight: 20.0,
            },
            TestAnchor {
                ref_pos: 200,
                query_pos: 100,
                weight: 20.0,
            }, // Not colinear!
        ];
        let result = chainer.chain(&anchors, &ChainParams::default());
        // Should pick the better single anchor, not chain them
        assert!(result.chain.len() <= 2);
    }

    #[test]
    fn test_respects_bandwidth() {
        let mut chainer = RmqChainer::new();
        let params = ChainParams {
            bandwidth: 10, // Very narrow
            ..Default::default()
        };
        // Anchors with very different diagonals should not chain
        let anchors = vec![
            TestAnchor {
                ref_pos: 100,
                query_pos: 100,
                weight: 20.0,
            }, // diag = 0
            TestAnchor {
                ref_pos: 200,
                query_pos: 150,
                weight: 20.0,
            }, // diag = 50, outside bandwidth
        ];
        let result = chainer.chain(&anchors, &params);
        // Should not chain due to bandwidth constraint
        assert!(result.chain.len() == 1);
    }

    #[test]
    fn test_diagonal_penalty_discourages_short_off_diag_anchor() {
        let anchors = vec![
            // Matches around chr11:116,818,0xx from selected.sam
            // Anchor A: ref 116,818,074; read 162,643; len 297
            TestAnchor {
                ref_pos: 116_818_074,
                query_pos: 162_643,
                weight: 297.0,
            },
            // Anchor B: short, off-diagonal anchor (len 29)
            TestAnchor {
                ref_pos: 116_818_472,
                query_pos: 163_320,
                weight: 29.0,
            },
            // Anchor C: ref 116,818,608; read 163,361; len 640
            TestAnchor {
                ref_pos: 116_818_608,
                query_pos: 163_361,
                weight: 640.0,
            },
        ];

        let mut chainer = RmqChainer::new();
        let params_no_diag = ChainParams {
            max_gap_ref: 10_000,
            max_gap_query: 10_000,
            gap_penalty_linear: 0.0,
            gap_penalty_log: 0.0,
            diagonal_penalty_linear: 0.0,
            bandwidth: 1_000,
        };
        let result_no_diag = chainer.chain(&anchors, &params_no_diag);
        assert_eq!(result_no_diag.chain, vec![0, 1, 2]);

        let mut chainer = RmqChainer::new();
        let params_with_diag = ChainParams {
            max_gap_ref: 10_000,
            max_gap_query: 10_000,
            gap_penalty_linear: 0.0,
            gap_penalty_log: 0.0,
            diagonal_penalty_linear: 0.5,
            bandwidth: 1_000,
        };
        let result_with_diag = chainer.chain(&anchors, &params_with_diag);
        assert_eq!(result_with_diag.chain, vec![0, 2]);
    }

    #[test]
    fn test_diagonal_penalty_prefers_closer_seed_over_29bp_anchor() {
        let anchors = vec![
            // Anchor A: ref 116,818,074; read 162,643; len 297
            TestAnchor {
                ref_pos: 116_818_074,
                query_pos: 162_643,
                weight: 297.0,
            },
            // In-between candidates (from selected.sam)
            TestAnchor {
                ref_pos: 116_818_340,
                query_pos: 162_826,
                weight: 22.0,
            },
            TestAnchor {
                ref_pos: 116_818_405,
                query_pos: 162_870,
                weight: 29.0,
            },
            TestAnchor {
                ref_pos: 116_818_427,
                query_pos: 162_870,
                weight: 41.0,
            },
            TestAnchor {
                ref_pos: 116_818_461,
                query_pos: 162_870,
                weight: 23.0,
            },
            // The problematic 29 bp anchor
            TestAnchor {
                ref_pos: 116_818_472,
                query_pos: 163_320,
                weight: 29.0,
            },
            // Additional candidates at same ref
            TestAnchor {
                ref_pos: 116_818_473,
                query_pos: 163_297,
                weight: 28.0,
            },
            TestAnchor {
                ref_pos: 116_818_473,
                query_pos: 163_264,
                weight: 38.0,
            },
            TestAnchor {
                ref_pos: 116_818_473,
                query_pos: 163_207,
                weight: 28.0,
            },
            TestAnchor {
                ref_pos: 116_818_473,
                query_pos: 163_074,
                weight: 38.0,
            },
            TestAnchor {
                ref_pos: 116_818_473,
                query_pos: 162_929,
                weight: 38.0,
            },
            // Anchor C: ref 116,818,608; read 163,361; len 640
            TestAnchor {
                ref_pos: 116_818_608,
                query_pos: 163_361,
                weight: 640.0,
            },
        ];

        for anchor in &anchors {
            println!(
                "Anchor: ref_pos={}, query_pos={}, diag={}, weight={}",
                anchor.ref_pos,
                anchor.query_pos,
                anchor.diagonal(),
                anchor.weight
            );
        }

        let mut chainer = RmqChainer::new();
        let params_no_diag = ChainParams {
            max_gap_ref: 10_000,
            max_gap_query: 10_000,
            gap_penalty_linear: 0.0,
            gap_penalty_log: 0.0,
            diagonal_penalty_linear: 0.0,
            bandwidth: 1_000,
        };
        let result_no_diag = chainer.chain(&anchors, &params_no_diag);
        assert!(result_no_diag.chain.contains(&5)); // 29bp anchor is index 5

        let mut chainer = RmqChainer::new();
        let params_with_diag = ChainParams {
            max_gap_ref: 10_000,
            max_gap_query: 10_000,
            gap_penalty_linear: 0.0,
            gap_penalty_log: 0.0,
            diagonal_penalty_linear: 0.5,
            bandwidth: 1_000,
        };
        let result_with_diag = chainer.chain(&anchors, &params_with_diag);
        println!("Chained anchors with diagonal penalty:");
        for &idx in &result_with_diag.chain {
            let anchor = &anchors[idx];
            println!(
                "  Anchor: ref_pos={}, query_pos={}, diag={}, weight={}",
                anchor.ref_pos,
                anchor.query_pos,
                anchor.diagonal(),
                anchor.weight
            );
        }
        assert!(!result_with_diag.chain.contains(&5));
        assert!(result_with_diag.chain.contains(&9)); // prefer 38bp anchor at query 163,074
    }

    #[test]
    fn test_anchor_first_chaining_excludes_off_diagonal_short_seeds() {
        // Scenario: Two long anchors on similar diagonals, with short off-diagonal seeds in between
        let anchors = vec![
            // Anchor A: long seed, diagonal = 116655431
            TestAnchor {
                ref_pos: 116_818_074,
                query_pos: 162_643,
                weight: 297.0,
            }, // diag = 116655431
            // Short seeds on very different diagonals (these should be excluded)
            TestAnchor {
                ref_pos: 116_818_400,
                query_pos: 162_700,
                weight: 25.0,
            }, // diag = 116655700 (off by 269)
            TestAnchor {
                ref_pos: 116_818_450,
                query_pos: 162_600,
                weight: 20.0,
            }, // diag = 116655850 (off by 419)
            TestAnchor {
                ref_pos: 116_818_472,
                query_pos: 163_320,
                weight: 29.0,
            }, // diag = 116655152 (off by 279 - way off!)
            // Short seed on the correct diagonal (should be included)
            // Position it so interpolated diagonal matches closely
            TestAnchor {
                ref_pos: 116_818_500,
                query_pos: 163_069,
                weight: 30.0,
            }, // diag = 116655431 (exactly on anchor A's diagonal)
            // Anchor B: long seed, diagonal = 116655431 (same as A for simplicity)
            TestAnchor {
                ref_pos: 116_818_608,
                query_pos: 163_177,
                weight: 640.0,
            }, // diag = 116655431
        ];

        let mut chainer = RmqChainer::new();
        let params = ChainParams {
            max_gap_ref: 10_000,
            max_gap_query: 10_000,
            gap_penalty_linear: 0.01,
            gap_penalty_log: 0.5,
            diagonal_penalty_linear: 0.5,
            bandwidth: 500,
        };

        // Anchor-first with threshold 40 (only 297bp and 640bp are anchors)
        let result = chainer.chain_anchor_first(&anchors, &params, 40.0, 50);

        println!("Anchor-first chain:");
        for &idx in &result.chain {
            let a = &anchors[idx];
            println!(
                "  idx={}, ref={}, query={}, diag={}, weight={}",
                idx,
                a.ref_pos,
                a.query_pos,
                a.diagonal(),
                a.weight
            );
        }

        // Should include both anchors
        assert!(result.chain.contains(&0)); // First anchor
        assert!(result.chain.contains(&5)); // Last anchor

        // Should include the seed on the correct diagonal
        assert!(result.chain.contains(&4)); // Short seed on same diagonal

        // Should NOT include the off-diagonal short seeds
        assert!(!result.chain.contains(&1)); // diag off by 269
        assert!(!result.chain.contains(&2)); // diag off by 419
        assert!(!result.chain.contains(&3)); // diag off by 279
    }

    #[test]
    fn test_anchor_first_falls_back_when_no_anchors() {
        let anchors = vec![
            // All short seeds, no anchors
            TestAnchor {
                ref_pos: 100,
                query_pos: 100,
                weight: 20.0,
            },
            TestAnchor {
                ref_pos: 200,
                query_pos: 200,
                weight: 25.0,
            },
            TestAnchor {
                ref_pos: 300,
                query_pos: 300,
                weight: 30.0,
            },
        ];

        let mut chainer = RmqChainer::new();
        let params = ChainParams::default();

        // With threshold 40, all seeds are fillers
        let result = chainer.chain_anchor_first(&anchors, &params, 40.0, 50);

        // Should fall back to regular chaining and chain all three
        assert_eq!(result.chain.len(), 3);
        assert_eq!(result.chain, vec![0, 1, 2]);
    }
}
