#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Selection<T, U> {
    Left(T),
    Right(U),
    Both(T, U),
}

pub mod debug;
pub mod frozen_big_table;
pub mod frozen_table;
pub mod join;
pub mod sequence;
pub mod swiss;
pub mod table;

#[allow(dead_code)]
pub struct GroupByKey<'a, F: Fn(&'a T) -> K, T, K: PartialEq> {
    items: &'a [T],
    key_fn: F,
    begin: usize,
    end: usize,
}

impl<'a, F: Fn(&'a T) -> K, T, K: PartialEq> GroupByKey<'a, F, T, K> {
    #[allow(dead_code)]
    pub fn new(iter: &'a [T], key_fn: F) -> Self {
        GroupByKey {
            items: iter,
            key_fn,
            begin: 0,
            end: 0,
        }
    }
}

impl<'a, F: Fn(&'a T) -> K, T, K: PartialEq> Iterator for GroupByKey<'a, F, T, K> {
    type Item = (K, &'a [T]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.end >= self.items.len() {
            return None;
        }

        self.begin = self.end;
        let key = (self.key_fn)(&self.items[self.begin]);
        self.end += 1;

        while self.end < self.items.len() {
            let next_key = (self.key_fn)(&self.items[self.end]);
            if &next_key != &key {
                break;
            }
            self.end += 1;
        }

        Some((key, &self.items[self.begin..self.end]))
    }
}

#[allow(dead_code)]
pub trait GroupByTrait<F: Fn(&Self::Item) -> K, K: PartialEq> {
    type Item;

    fn group_by(&'_ self, key_fn: F) -> GroupByKey<'_, F, Self::Item, K>;
}

impl<T, F: Fn(&T) -> K, K: PartialEq> GroupByTrait<F, K> for [T] {
    type Item = T;

    fn group_by(&'_ self, key_fn: F) -> GroupByKey<'_, F, T, K> {
        GroupByKey::new(self, key_fn)
    }
}

/// Returns cluster boundaries as indices into `points`.
/// The returned vector always starts with 0 and ends with points.len().
/// A new cluster starts at i when points[i] - points[i-1] > eps.
#[allow(dead_code)]
pub fn dbscan_1d_boundaries<T, F: Fn(&T) -> i64>(
    points: &[T],
    eps: i64,
    key: F,
    cuts: &mut Vec<usize>,
) {
    cuts.clear();
    cuts.push(0);

    let n = points.len();
    if n == 0 {
        return;
    }

    for i in 1..n {
        if key(&points[i]) - key(&points[i - 1]) > eps {
            cuts.push(i);
        }
    }

    cuts.push(n);
}

#[allow(dead_code)]
pub struct Summarizer {
    pub count: usize,
    pub total: f64,
    pub total_sq: f64,
}

impl Summarizer {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Summarizer {
            count: 0,
            total: 0.0,
            total_sq: 0.0,
        }
    }

    #[allow(dead_code)]
    pub fn add(&mut self, value: f64) {
        self.count += 1;
        self.total += value;
        self.total_sq += value * value;
    }

    #[allow(dead_code)]
    pub fn add_multiple(&mut self, value: f64, n: usize) {
        self.count += n;
        self.total += value * (n as f64);
        self.total_sq += value * value * (n as f64);
    }

    #[allow(dead_code)]
    pub fn mean(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.total / (self.count as f64)
        }
    }

    #[allow(dead_code)]
    pub fn stddev(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            let mean = self.mean();
            let variance = (self.total_sq / (self.count as f64)) - (mean * mean);
            variance.sqrt()
        }
    }
}

/// Cluster with adaptive merging based on within-cluster variance.
/// After initial gap-based clustering with `initial_eps`, adjacent clusters
/// are merged if their combined variance stays below `max_var`.
pub fn dbscan_variance_aware<T, F>(
    items: &[T],
    initial_eps: i64,
    max_var: f64,
    key: F,
    cuts: &mut Vec<usize>,
) where
    F: Fn(&T) -> i64,
{
    let n = items.len();
    cuts.clear();
    cuts.push(0);

    if n == 0 {
        return;
    }

    // Initial gap-based clustering with tight epsilon.
    for i in 1..n {
        if key(&items[i]) - key(&items[i - 1]) > initial_eps {
            cuts.push(i);
        }
    }
    cuts.push(n);

    // Merge adjacent clusters if combined variance stays below threshold.
    loop {
        let mut merged = false;
        let mut new_cuts = Vec::with_capacity(cuts.len());
        new_cuts.push(0);

        let mut i = 1;
        while i < cuts.len() {
            let prev_start = *new_cuts.last().unwrap();
            let curr_end = cuts[i];

            if i + 1 < cuts.len() {
                let next_end = cuts[i + 1];
                // Try merging current and next cluster.
                let combined_var = compute_variance(items, prev_start, next_end, &key);
                if combined_var <= max_var {
                    // Skip the intermediate boundary; merge.
                    i += 1;
                    merged = true;
                    continue;
                }
            }

            new_cuts.push(curr_end);
            i += 1;
        }

        *cuts = new_cuts;
        if !merged {
            break;
        }
    }
}

fn compute_variance<T, F>(items: &[T], start: usize, end: usize, key: &F) -> f64
where
    F: Fn(&T) -> i64,
{
    let n = end - start;
    if n < 2 {
        return 0.0;
    }
    let mut sum = 0i64;
    let mut sum_sq = 0i128;
    for i in start..end {
        let v = key(&items[i]);
        sum += v;
        sum_sq += (v as i128) * (v as i128);
    }
    let mean = sum as f64 / n as f64;
    let var = (sum_sq as f64 / n as f64) - mean * mean;
    var.max(0.0)
}

/// Filter a cluster to keep only the longest colinear chain.
/// For forward mapping (colinear=true): ref_pos should increase with read_pos.
/// For reverse mapping (colinear=false): ref_pos should decrease with read_pos.
///
/// `items` should be sorted by read position.
/// `ref_pos_key` extracts the reference position from each item.
///
/// Returns indices into `items` of the longest colinear/anti-colinear chain.
#[allow(dead_code)]
pub fn longest_colinear_chain<T, F>(items: &[T], ref_pos_key: F, forward: bool) -> Vec<usize>
where
    F: Fn(&T) -> i64,
{
    let mut result = Vec::new();
    LongestSubsequence::default().longest_colinear_chain(items, ref_pos_key, forward, &mut result);
    result
}

/// Reusable buffers for computing Longest Increasing/Decreasing Subsequence.
/// Use this when calling LIS/LDS repeatedly to avoid repeated allocations.
#[derive(Default)]
pub struct LongestSubsequence {
    /// tails[i] = index of smallest tail value for LIS of length i+1
    tails: Vec<usize>,
    /// prev[i] = predecessor index for items[i] in the LIS
    prev: Vec<usize>,
}

impl LongestSubsequence {
    /// Clear internal buffers and resize for `n` items.
    fn clear(&mut self, n: usize) {
        self.tails.clear();
        self.prev.clear();
        self.prev.resize(n, usize::MAX);
    }

    /// Compute the Longest Increasing Subsequence (LIS) on the values extracted by `key`.
    /// Results are written to `result` (cleared first).
    /// Items are assumed to already be sorted by a primary key (e.g., read position).
    /// O(n log n) time complexity.
    pub fn longest_increasing_subsequence<T, F>(
        &mut self,
        items: &[T],
        key: F,
        result: &mut Vec<usize>,
    ) where
        F: Fn(&T) -> i64,
    {
        result.clear();
        let n = items.len();
        if n == 0 {
            return;
        }

        self.clear(n);

        for i in 0..n {
            let val = key(&items[i]);
            // Binary search for the position in tails
            let pos = self.tails.partition_point(|&t| key(&items[t]) < val);

            if pos == self.tails.len() {
                self.tails.push(i);
            } else {
                self.tails[pos] = i;
            }

            if pos > 0 {
                self.prev[i] = self.tails[pos - 1];
            }
        }

        // Reconstruct the LIS
        result.reserve(self.tails.len());
        let mut idx = *self.tails.last().unwrap();
        while idx != usize::MAX {
            result.push(idx);
            idx = self.prev[idx];
        }
        result.reverse();
    }

    /// Compute the Longest Decreasing Subsequence (LDS) on the values extracted by `key`.
    /// Results are written to `result` (cleared first).
    /// Items are assumed to already be sorted by a primary key (e.g., read position).
    /// O(n log n) time complexity.
    pub fn longest_decreasing_subsequence<T, F>(
        &mut self,
        items: &[T],
        key: F,
        result: &mut Vec<usize>,
    ) where
        F: Fn(&T) -> i64,
    {
        result.clear();
        let n = items.len();
        if n == 0 {
            return;
        }

        self.clear(n);

        for i in 0..n {
            let val = key(&items[i]);
            // Binary search for position where val would go (decreasing order)
            let pos = self.tails.partition_point(|&t| key(&items[t]) > val);

            if pos == self.tails.len() {
                self.tails.push(i);
            } else {
                self.tails[pos] = i;
            }

            if pos > 0 {
                self.prev[i] = self.tails[pos - 1];
            }
        }

        // Reconstruct the LDS
        result.reserve(self.tails.len());
        let mut idx = *self.tails.last().unwrap();
        while idx != usize::MAX {
            result.push(idx);
            idx = self.prev[idx];
        }
        result.reverse();
    }

    /// Filter a cluster to keep only the longest colinear chain.
    /// For forward mapping (forward=true): ref_pos should increase with read_pos.
    /// For reverse mapping (forward=false): ref_pos should decrease with read_pos.
    ///
    /// `items` should be sorted by read position.
    /// `ref_pos_key` extracts the reference position from each item.
    /// Results are written to `result` (cleared first).
    pub fn longest_colinear_chain<T, F>(
        &mut self,
        items: &[T],
        ref_pos_key: F,
        forward: bool,
        result: &mut Vec<usize>,
    ) where
        F: Fn(&T) -> i64,
    {
        if forward {
            self.longest_increasing_subsequence(items, ref_pos_key, result);
        } else {
            self.longest_decreasing_subsequence(items, ref_pos_key, result);
        }
    }
}

/// Compute the longest common prefix of two byte slices.
///
/// This implementation processes 16 bytes at a time using patterns that
/// LLVM can optimize to SIMD instructions (SSE2/NEON).
#[inline]
pub fn longest_common_prefix(lhs: &[u8], rhs: &[u8]) -> usize {
    const BLOCK: usize = 16;

    let len = lhs.len().min(rhs.len());
    let mut offset = 0;

    // Process 16 bytes at a time - this vectorizes to pcmpeqb + pmovmskb
    while offset + BLOCK <= len {
        // Compare 16 bytes at once, producing a bitmask of mismatches
        let mismatch_mask = compare_block_16(
            &lhs[offset..offset + BLOCK],
            &rhs[offset..offset + BLOCK],
        );

        if mismatch_mask != 0 {
            // Found a mismatch - return position of first differing byte
            return offset + mismatch_mask.trailing_zeros() as usize;
        }

        offset += BLOCK;
    }

    // Handle remaining bytes (< 16)
    while offset < len {
        if lhs[offset] != rhs[offset] {
            return offset;
        }
        offset += 1;
    }

    len
}

/// Compare two 16-byte slices and return a bitmask where bit i is set if bytes differ.
///
/// This function is structured to enable LLVM auto-vectorization.
#[inline(always)]
fn compare_block_16(a: &[u8], b: &[u8]) -> u32 {
    debug_assert!(a.len() >= 16 && b.len() >= 16);

    // Load into fixed-size arrays for better optimization
    let mut block_a = [0u8; 16];
    let mut block_b = [0u8; 16];
    block_a.copy_from_slice(&a[..16]);
    block_b.copy_from_slice(&b[..16]);

    // XOR the blocks - equal bytes become 0, different bytes become non-zero
    // Then check which positions are non-zero
    let mut mismatch_mask = 0u32;
    for i in 0..16 {
        // XOR detects difference, comparison packs into bitmask
        // This pattern compiles to: pxor + pcmpeqb + pmovmskb (inverted)
        mismatch_mask |= ((block_a[i] != block_b[i]) as u32) << i;
    }
    mismatch_mask
}

#[cfg(test)]
mod lcp_tests {
    use super::longest_common_prefix;

    #[test]
    fn test_identical() {
        let a = b"ACGTACGTACGTACGT";
        let b = b"ACGTACGTACGTACGT";
        assert_eq!(longest_common_prefix(a, b), 16);
    }

    #[test]
    fn test_differ_at_start() {
        let a = b"ACGTACGTACGTACGT";
        let b = b"XCGTACGTACGTACGT";
        assert_eq!(longest_common_prefix(a, b), 0);
    }

    #[test]
    fn test_differ_in_middle() {
        let a = b"ACGTACGTACGTACGT";
        let b = b"ACGTACGXACGTACGT";
        assert_eq!(longest_common_prefix(a, b), 7);
    }

    #[test]
    fn test_differ_at_end() {
        let a = b"ACGTACGTACGTACGT";
        let b = b"ACGTACGTACGTACGX";
        assert_eq!(longest_common_prefix(a, b), 15);
    }

    #[test]
    fn test_short_slices() {
        assert_eq!(longest_common_prefix(b"ACG", b"ACG"), 3);
        assert_eq!(longest_common_prefix(b"ACG", b"ACX"), 2);
        assert_eq!(longest_common_prefix(b"ACG", b"XCG"), 0);
    }

    #[test]
    fn test_different_lengths() {
        let a = b"ACGTACGTACGTACGTACGT";
        let b = b"ACGTACGTACGTACGT";
        assert_eq!(longest_common_prefix(a, b), 16);
        assert_eq!(longest_common_prefix(b, a), 16);
    }

    #[test]
    fn test_empty() {
        assert_eq!(longest_common_prefix(b"", b"ACGT"), 0);
        assert_eq!(longest_common_prefix(b"ACGT", b""), 0);
        assert_eq!(longest_common_prefix(b"", b""), 0);
    }

    #[test]
    fn test_longer_than_16() {
        let a = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bytes
        let b = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        assert_eq!(longest_common_prefix(a, b), 32);

        let c = b"ACGTACGTACGTACGTACGTACGTACGTXCGT"; // differs at position 28
        assert_eq!(longest_common_prefix(a, c), 28);
    }

    #[test]
    fn test_boundary_at_16() {
        let a = b"ACGTACGTACGTACGTX";
        let b = b"ACGTACGTACGTACGTY";
        assert_eq!(longest_common_prefix(a, b), 16);
    }
}