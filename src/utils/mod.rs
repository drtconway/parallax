#![allow(dead_code)] // This module contains many utilities that are not used in all configurations

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Selection<T, U> {
    Left(T),
    Right(U),
    Both(T, U),
}

pub mod debug;
pub mod dump;
pub mod frozen_big_table;
pub mod frozen_table;
pub mod hasher;
pub mod heap;
pub mod human;
pub mod join;
pub mod pool;
pub mod progress;
pub mod range_set;
pub mod rope;
pub mod sequence;
pub mod swiss;
pub mod table;
pub mod telemetry;
pub mod union_find;

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

/// Interleave two iterators.
pub struct Interleave<I: Iterator, J: Iterator<Item = I::Item>> {
    iter1: I,
    iter2: J,
    turn: bool,
}

impl<I: Iterator, J: Iterator<Item = I::Item>> Iterator for Interleave<I, J> {
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        self.turn = !self.turn;
        if self.turn {
            self.iter1.next().or_else(|| self.iter2.next())
        } else {
            self.iter2.next().or_else(|| self.iter1.next())
        }
    }
}

pub trait InterleaveTrait: Iterator + Sized {
    fn interleave<J: Iterator<Item = Self::Item>>(self, other: J) -> Interleave<Self, J> {
        Interleave { iter1: self, iter2: other, turn: false }
    }
}

impl<I: Iterator> InterleaveTrait for I {}

#[allow(dead_code)]
pub fn which_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .min_by_key(|&(_, value)| value)
        .map(|(index, _)| index)
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

pub struct PairsIterator<'a, T, I: Iterator<Item = &'a T>> {
    iter: I,
    peeked: Option<&'a T>,
    phantom: std::marker::PhantomData<&'a T>,
}

impl<'a, T, I: Iterator<Item = &'a T>> PairsIterator<'a, T, I> {
    #[allow(dead_code)]
    pub fn new(iter: I) -> Self {
        Self {
            iter,
            peeked: None,
            phantom: std::marker::PhantomData,
        }
    }
}

impl<'a, T, I: Iterator<Item = &'a T>> Iterator for PairsIterator<'a, T, I> {
    type Item = (&'a T, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let first = match self.peeked.take() {
            Some(v) => v,
            None => self.iter.next()?,
        };

        let second = self.iter.next()?;
        self.peeked = Some(second);

        Some((first, second))
    }
}

#[allow(dead_code)]
pub trait PairsTrait<'a, T: 'a, I: Iterator<Item = &'a T>> {
    fn pairs(self) -> PairsIterator<'a, T, I>;
}

impl<'a, T: 'a, I: Iterator<Item = &'a T>> PairsTrait<'a, T, I> for I {
    fn pairs(self) -> PairsIterator<'a, T, I> {
        PairsIterator::new(self)
    }
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

pub struct GroupsIterator<T, I: Iterator<Item = Option<T>>> {
    iter: I,
    phantom: std::marker::PhantomData<T>,
}

impl<T, I: Iterator<Item = Option<T>>> Iterator for GroupsIterator<T, I> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut group = Vec::new();

        while let Some(item) = self.iter.next() {
            match item {
                Some(value) => group.push(value),
                None => {
                    if !group.is_empty() {
                        return Some(group);
                    }
                }
            }
        }

        if !group.is_empty() {
            Some(group)
        } else {
            None
        }
    }
}

pub trait GroupsTrait<T>: Iterator<Item = Option<T>> + Sized {
    fn groups(self) -> GroupsIterator<T, Self> {
        GroupsIterator { iter: self, phantom: std::marker::PhantomData }
    }
}

impl<T, I: Iterator<Item = Option<T>>> GroupsTrait<T> for I {
    fn groups(self) -> GroupsIterator<T, Self> {
        GroupsIterator { iter: self, phantom: std::marker::PhantomData }
    }
}

/// Compute the longest common prefix of two byte slices.
///
/// This implementation processes 16 bytes at a time using patterns that
/// LLVM can optimize to SIMD instructions (SSE2/NEON).
#[allow(dead_code)]
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
#[allow(dead_code)]
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

#[cfg(test)]
mod interleave_tests {
    use super::InterleaveTrait;

    #[test]
    fn test_equal_length() {
        let a = vec![1, 3, 5];
        let b = vec![2, 4, 6];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn test_first_longer() {
        let a = vec![1, 3, 5, 7, 9];
        let b = vec![2, 4];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec![1, 2, 3, 4, 5, 7, 9]);
    }

    #[test]
    fn test_second_longer() {
        let a = vec![1, 3];
        let b = vec![2, 4, 6, 8, 10];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec![1, 2, 3, 4, 6, 8, 10]);
    }

    #[test]
    fn test_first_empty() {
        let a: Vec<i32> = vec![];
        let b = vec![2, 4, 6];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec![2, 4, 6]);
    }

    #[test]
    fn test_second_empty() {
        let a = vec![1, 3, 5];
        let b: Vec<i32> = vec![];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec![1, 3, 5]);
    }

    #[test]
    fn test_both_empty() {
        let a: Vec<i32> = vec![];
        let b: Vec<i32> = vec![];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert!(result.is_empty());
    }

    #[test]
    fn test_single_elements() {
        let a = vec![1];
        let b = vec![2];
        let result: Vec<i32> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec![1, 2]);
    }

    #[test]
    fn test_with_strings() {
        let a = vec!["a", "c", "e"];
        let b = vec!["b", "d", "f"];
        let result: Vec<&str> = a.into_iter().interleave(b.into_iter()).collect();
        assert_eq!(result, vec!["a", "b", "c", "d", "e", "f"]);
    }

    #[test]
    fn test_chained_interleave() {
        let a = vec![1, 4];
        let b = vec![2, 5];
        let c = vec![3, 6];
        // First interleave a and b: [1, 2, 4, 5]
        // Then interleave with c: [1, 3, 2, 6, 4, 5]
        let result: Vec<i32> = a.into_iter()
            .interleave(b.into_iter())
            .interleave(c.into_iter())
            .collect();
        assert_eq!(result, vec![1, 3, 2, 6, 4, 5]);
    }
}
