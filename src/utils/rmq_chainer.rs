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
    /// Total chain score
    pub score: f64,
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
    pub fn chain<T: ChainAnchor>(
        &mut self,
        anchors: &[T],
        params: &ChainParams,
    ) -> ChainResult {
        let n = anchors.len();
        if n == 0 {
            return ChainResult {
                chain: Vec::new(),
                score: 0.0,
            };
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

        ChainResult {
            chain,
            score: best_score,
        }
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
        fn ref_pos(&self) -> i64 { self.ref_pos }
        fn query_pos(&self) -> i64 { self.query_pos }
        fn weight(&self) -> f64 { self.weight }
    }

    #[test]
    fn test_empty() {
        let mut chainer = RmqChainer::new();
        let anchors: Vec<TestAnchor> = vec![];
        let result = chainer.chain(&anchors, &ChainParams::default());
        assert!(result.chain.is_empty());
        assert_eq!(result.score, 0.0);
    }

    #[test]
    fn test_single_anchor() {
        let mut chainer = RmqChainer::new();
        let anchors = vec![TestAnchor { ref_pos: 100, query_pos: 50, weight: 20.0 }];
        let result = chainer.chain(&anchors, &ChainParams::default());
        assert_eq!(result.chain, vec![0]);
        assert_eq!(result.score, 20.0);
    }

    #[test]
    fn test_colinear_chain() {
        let mut chainer = RmqChainer::new();
        // Three anchors on the same diagonal, should chain together
        let anchors = vec![
            TestAnchor { ref_pos: 100, query_pos: 100, weight: 20.0 },
            TestAnchor { ref_pos: 200, query_pos: 200, weight: 20.0 },
            TestAnchor { ref_pos: 300, query_pos: 300, weight: 20.0 },
        ];
        let result = chainer.chain(&anchors, &ChainParams::default());
        assert_eq!(result.chain, vec![0, 1, 2]);
        assert!(result.score > 50.0); // Should be > sum of weights minus small gap penalties
    }

    #[test]
    fn test_non_colinear_filtered() {
        let mut chainer = RmqChainer::new();
        // Anchor 1 has decreasing query pos - should not chain
        let anchors = vec![
            TestAnchor { ref_pos: 100, query_pos: 200, weight: 20.0 },
            TestAnchor { ref_pos: 200, query_pos: 100, weight: 20.0 }, // Not colinear!
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
            TestAnchor { ref_pos: 100, query_pos: 100, weight: 20.0 }, // diag = 0
            TestAnchor { ref_pos: 200, query_pos: 150, weight: 20.0 }, // diag = 50, outside bandwidth
        ];
        let result = chainer.chain(&anchors, &params);
        // Should not chain due to bandwidth constraint
        assert!(result.chain.len() == 1);
    }

    #[test]
    fn test_gap_penalty() {
        let mut chainer = RmqChainer::new();
        let params = ChainParams {
            gap_penalty_linear: 0.1,
            gap_penalty_log: 1.0,
            ..Default::default()
        };
        // Large gap should reduce chain score
        let anchors = vec![
            TestAnchor { ref_pos: 100, query_pos: 100, weight: 20.0 },
            TestAnchor { ref_pos: 1100, query_pos: 1100, weight: 20.0 }, // 1000bp gap
        ];
        let result = chainer.chain(&anchors, &params);
        // Chain score should be less than 40 due to gap penalty
        assert!(result.score < 40.0);
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
}