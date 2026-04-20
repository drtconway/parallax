//! Colinear seed chaining via a diagonal-indexed segment tree (O(n log n)).
//!
//! # Algorithm
//!
//! Seeds are processed in ascending `read_pos` order. For each seed `i` we
//! need the predecessor `j` that maximises the score function
//!
//! ```text
//!   f[i] = max over j: f[j] + alpha(i, j) - beta(i, j)
//! ```
//!
//! where alpha is the non-overlapping match bonus and beta is the gap
//! penalty (same scoring as the Kruskal implementation).
//!
//! The banding constraint `|diag_i - diag_j| ≤ MAX_DIAGONAL_DIST` is
//! enforced by mapping each diagonal value to a bucket in a segment tree
//! that supports O(log B) range-max queries and point updates, where B is
//! the number of distinct diagonal buckets.
//!
//! Reference colinearity (`ref_pos_j < ref_pos_i`) is not directly
//! enforced by the segment tree, but because we process seeds in read
//! order and the diagonal band is tight, non-colinear transitions produce
//! a poor (or negative) alpha and will not win the max.
//!
//! After the DP, chains are recovered by traceback from the highest-scored
//! unused endpoint, as in the rmq_dp implementation.

use crate::reads::seeds::{SeedHit, seed_cluster::SeedCluster};

// ── Constants ────────────────────────────────────────────────────────────────

/// Maximum allowed difference between the diagonals of two chained seeds.
const MAX_DIAGONAL_DIST: i64 = 2000;

/// Gap penalty weight (identical to the kruskal implementation).
const W: f64 = 2.0;

/// Minimum chain score to report a chain.
const MIN_CHAIN_SCORE: f64 = 75.0;

/// Maximum number of chains to return per chromosome strand.
const MAX_CHAINS: usize = 100;

// ── Scoring ──────────────────────────────────────────────────────────────────

/// Score the transition from predecessor `j` to seed `i`.
///
/// Returns `None` if the pair is non-colinear or outside the diagonal band.
/// Otherwise returns `f[j] + alpha - beta`.
fn transition_score(
    seeds: &[SeedHit],
    i: usize,
    j: usize,
    fj: f64,
) -> Option<f64> {
    let si = &seeds[i];
    let sj = &seeds[j];

    let q_i = si.read_pos as i64;
    let r_i = si.ref_pos as i64;
    let l_i = si.match_len as i64;
    let q_j = sj.read_pos as i64;
    let r_j = sj.ref_pos as i64;
    let l_j = sj.match_len as i64;

    // Diagonal banding
    let d_i = r_i - q_i;
    let d_j = r_j - q_j;
    if (d_i - d_j).abs() > MAX_DIAGONAL_DIST {
        return None;
    }

    // Colinearity: j must be strictly before i in both coordinates
    if q_j >= q_i || r_j >= r_i {
        return None;
    }

    let gap_q = q_i - q_j;
    let gap_r = r_i - r_j;

    // Match bonus: non-overlapping portion of seed j that fits in the gap
    let alpha = l_i.min(l_j).min(gap_q).min(gap_r) as f64;

    // Gap penalty: affine-like, same formula as kruskal/rmq_dp
    let g = (gap_q - gap_r).unsigned_abs() as usize;
    let diag_dev = (d_i - d_j).abs() as f64;
    let beta = if g == 0 {
        0.0
    } else {
        0.01 * W * g as f64 + 0.5 * (g as f64).log2()
    } + 0.05 * diag_dev;

    Some(fj + alpha - beta)
}

// ── Segment tree (range max) ─────────────────────────────────────────────────
//
// A simple 1-indexed segment tree over B leaves. Each leaf stores the best
// (score, seed_index) pair seen so far for that diagonal bucket. Interior
// nodes store the max.
//
// We store (f64, usize) so we can recover the predecessor index directly.

const NEG_INF: f64 = f64::NEG_INFINITY;

struct SegTree {
    n: usize,            // number of leaves (rounded up to power of 2)
    data: Vec<(f64, usize)>, // (score, seed_index); usize::MAX = no entry
}

impl SegTree {
    fn new(leaves: usize) -> Self {
        // Round up to next power of two for simple 1-indexed layout
        let n = leaves.next_power_of_two();
        Self {
            n,
            data: vec![(NEG_INF, usize::MAX); 2 * n],
        }
    }

    /// Update leaf `pos` if `val` is better than the current value.
    fn update(&mut self, pos: usize, val: (f64, usize)) {
        debug_assert!(pos < self.n);
        let mut idx = self.n + pos;
        if val.0 > self.data[idx].0 {
            self.data[idx] = val;
        }
        idx >>= 1;
        while idx >= 1 {
            let left = self.data[2 * idx];
            let right = self.data[2 * idx + 1];
            let best = if left.0 >= right.0 { left } else { right };
            if best == self.data[idx] {
                break; // no change propagates upward
            }
            self.data[idx] = best;
            idx >>= 1;
        }
    }

    /// Query the maximum (score, seed_index) over leaves `[lo, hi]` (inclusive).
    fn query(&self, lo: usize, hi: usize) -> (f64, usize) {
        let mut lo = lo + self.n;
        let mut hi = hi + self.n + 1; // make exclusive
        let mut best = (NEG_INF, usize::MAX);
        while lo < hi {
            if lo & 1 == 1 {
                if self.data[lo].0 > best.0 {
                    best = self.data[lo];
                }
                lo += 1;
            }
            if hi & 1 == 1 {
                hi -= 1;
                if self.data[hi].0 > best.0 {
                    best = self.data[hi];
                }
            }
            lo >>= 1;
            hi >>= 1;
        }
        best
    }
}

// ── Main DP ──────────────────────────────────────────────────────────────────

/// Run the O(n log n) chaining DP and return per-seed (score, predecessor).
///
/// `seeds` must be sorted by `read_pos` ascending before calling.
///
/// The segment tree is indexed by *rank* of the diagonal value, not by the
/// raw diagonal coordinate.  This compresses the tree to at most `n` leaves
/// regardless of the genomic position of the seeds, avoiding the 256 M-leaf
/// tree that would result from using genomic coordinates directly.
fn chain_dp(seeds: &[SeedHit]) -> (Vec<f64>, Vec<i32>) {
    let n = seeds.len();
    let mut f = vec![0.0f64; n];
    let mut pred = vec![-1i32; n];

    if n == 0 {
        return (f, pred);
    }

    // Coordinate-compress diagonals: collect unique values, sort them, and
    // use each value's rank (index in this sorted list) as the tree bucket.
    // The number of distinct diagonals is at most n, so the tree has at most
    // n leaves instead of up to ~250 M leaves for genomic coordinates.
    let mut sorted_diags: Vec<i64> = seeds.iter().map(|s| s.diagonal).collect();
    sorted_diags.sort_unstable();
    sorted_diags.dedup();
    let num_diags = sorted_diags.len();

    let mut tree = SegTree::new(num_diags);

    for i in 0..n {
        let si = &seeds[i];
        // Base score: just this seed on its own
        f[i] = si.match_len as f64;

        let d_i = si.diagonal;

        // Translate the diagonal band [d_i - D, d_i + D] to a rank range
        // using binary search on the compressed diagonal list.
        let lo_rank = sorted_diags.partition_point(|&d| d < d_i - MAX_DIAGONAL_DIST);
        let hi_rank_excl = sorted_diags.partition_point(|&d| d <= d_i + MAX_DIAGONAL_DIST);

        if lo_rank < hi_rank_excl {
            let (best_score, best_j) = tree.query(lo_rank, hi_rank_excl - 1);
            if best_j != usize::MAX && best_score > NEG_INF {
                if let Some(score) = transition_score(seeds, i, best_j, best_score) {
                    if score > f[i] {
                        f[i] = score;
                        pred[i] = best_j as i32;
                    }
                }
            }
        }

        // Insert seed i at its diagonal rank.
        let rank = sorted_diags.partition_point(|&d| d < d_i);
        tree.update(rank, (f[i], i));
    }

    (f, pred)
}

// ── Chain extraction ─────────────────────────────────────────────────────────

fn extract_chains(
    f: &[f64],
    pred: &[i32],
    seeds: &[SeedHit],
    is_reverse: bool,
) -> Vec<SeedCluster> {
    let n = seeds.len();
    let mut used = vec![false; n];
    let mut chains = Vec::new();

    loop {
        if chains.len() >= MAX_CHAINS {
            break;
        }

        // Find the best unused endpoint.
        let Some((end, _)) = f
            .iter()
            .enumerate()
            .filter(|&(i, _)| !used[i])
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        else {
            break;
        };

        if f[end] < MIN_CHAIN_SCORE {
            break;
        }

        // Traceback.
        let mut chain_indices = Vec::new();
        let mut cur = end as i32;
        while cur >= 0 {
            let idx = cur as usize;
            if used[idx] {
                break;
            }
            chain_indices.push(idx);
            used[idx] = true;
            cur = pred[idx];
        }

        if chain_indices.is_empty() {
            break;
        }

        chain_indices.reverse(); // now in ascending read_pos order

        let chain_seeds: Vec<SeedHit> = chain_indices
            .iter()
            .map(|&idx| seeds[idx].clone())
            .collect();

        if let Some(cluster) = SeedCluster::new(chain_seeds, is_reverse, 8) {
            chains.push(cluster);
        }
    }

    chains
}

// ── Public entry point ───────────────────────────────────────────────────────

pub fn collect_chains(
    seeds: &mut [SeedHit],
    chrom_name: &str,
    is_reverse: bool,
) -> Vec<SeedCluster> {
    // Sort by read position; the DP requires this ordering.
    seeds.sort_unstable_by_key(|s| s.read_pos);

    let (f, pred) = chain_dp(seeds);

    let chains = extract_chains(&f, &pred, seeds, is_reverse);

    log::debug!(
        "Fenwick chaining produced {} chains on {}",
        chains.len(),
        chrom_name
    );

    chains
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(ref_pos: usize, read_pos: usize, match_len: usize) -> SeedHit {
        SeedHit::new(0, ref_pos, read_pos, 0, 1, match_len)
    }

    #[test]
    fn test_seg_tree_single_update_query() {
        let mut t = SegTree::new(16);
        t.update(5, (3.0, 42));
        let (score, idx) = t.query(4, 7);
        assert_eq!(idx, 42);
        assert!((score - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_seg_tree_range_max() {
        let mut t = SegTree::new(16);
        t.update(2, (1.0, 1));
        t.update(7, (5.0, 7));
        t.update(10, (3.0, 10));
        let (score, idx) = t.query(5, 9);
        assert_eq!(idx, 7);
        assert!((score - 5.0).abs() < 1e-9);
        // Outside that range:
        let (score2, _) = t.query(0, 3);
        assert!((score2 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_chain_dp_three_colinear_seeds() {
        // Three perfectly colinear seeds (same diagonal)
        let mut seeds = vec![
            seed(100, 10, 20),
            seed(130, 40, 20),
            seed(160, 70, 20),
        ];
        seeds.sort_unstable_by_key(|s| s.read_pos);
        let (f, pred) = chain_dp(&seeds);
        // The last seed should have the highest score and a valid predecessor chain.
        assert!(f[2] > f[0]);
        assert!(pred[2] >= 0);
    }

    #[test]
    fn test_collect_chains_simple() {
        let mut seeds = vec![
            seed(100, 10, 30),
            seed(130, 40, 30),
            seed(160, 70, 30),
        ];
        let chains = collect_chains(&mut seeds, "chr1", false);
        assert!(!chains.is_empty());
    }

    #[test]
    fn test_transition_score_non_colinear() {
        let seeds = vec![seed(100, 50, 20), seed(50, 80, 20)]; // ref not colinear
        let result = transition_score(&seeds, 1, 0, 20.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_transition_score_diagonal_too_far() {
        // diagonal of seed 0 = 100 - 50 = 50; diagonal of seed 1 = 10000 - 10 = 9990; diff > 2000
        let seeds = vec![seed(100, 50, 20), seed(10000, 10, 20)];
        let result = transition_score(&seeds, 0, 1, 20.0);
        assert!(result.is_none());
    }
}
