use std::collections::VecDeque;
use std::ops::{Bound, RangeBounds};

pub struct SetCoverage {
    // Stored as half-open intervals [start, end).
    intervals: VecDeque<(usize, usize)>,
}

impl SetCoverage {
    pub fn new() -> Self {
        SetCoverage {
            intervals: VecDeque::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Insert a range, which may be half-open (`start..end`) or inclusive
    /// (`start..=end`).  Internally all intervals are stored as half-open.
    /// Adjacent intervals (e.g. `0..3` and `3..5`) are merged.
    pub fn insert(&mut self, range: impl RangeBounds<usize>) {
        let start = match range.start_bound() {
            Bound::Included(&s) => s,
            Bound::Excluded(&s) => s + 1,
            Bound::Unbounded => panic!("unbounded start not supported"),
        };
        let end = match range.end_bound() {
            Bound::Excluded(&e) => e,
            Bound::Included(&e) => e + 1,
            Bound::Unbounded => panic!("unbounded end not supported"),
        };
        assert!(start < end, "invalid interval: start ({start}) must be less than end ({end})");

        // Binary search for the first stored interval that touches or overlaps
        // [start, end).  For half-open intervals [s,e), "strictly before our
        // start" means e < start (no shared point, not even touching).
        let mut first = 0;
        let mut count = self.intervals.len();
        while count > 0 {
            let idx = count / 2;
            let mid = first + idx;
            if self.intervals[mid].1 < start {
                first = mid + 1;
                count -= idx + 1;
            } else {
                count = idx;
            }
        }

        // Merge all stored intervals whose start <= end (they touch or overlap).
        if first < self.intervals.len() && self.intervals[first].0 <= end {
            let (mut new_start, mut new_end) = (start, end);
            while first < self.intervals.len() && self.intervals[first].0 <= new_end {
                let (s, e) = self.intervals.remove(first).unwrap();
                new_start = new_start.min(s);
                new_end = new_end.max(e);
            }
            self.intervals.insert(first, (new_start, new_end));
        } else {
            self.intervals.insert(first, (start, end));
        }
    }

    /// Iterate over stored intervals as half-open `(start, end)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = &(usize, usize)> {
        self.intervals.iter()
    }

    pub fn start(&self) -> Option<usize> {
        self.intervals.front().map(|(start, _)| *start)
    }

    pub fn end(&self) -> Option<usize> {
        self.intervals.back().map(|(_, end)| *end)
    }

    pub fn coverage(&self) -> usize {
        self.intervals.iter().map(|(start, end)| end - start).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn intervals(cov: &SetCoverage) -> Vec<(usize, usize)> {
        cov.iter().copied().collect()
    }

    #[test]
    fn single_half_open() {
        let mut cov = SetCoverage::new();
        cov.insert(2..5);
        assert_eq!(intervals(&cov), [(2, 5)]);
        assert_eq!(cov.coverage(), 3);
    }

    #[test]
    fn single_inclusive() {
        let mut cov = SetCoverage::new();
        cov.insert(2..=4);
        assert_eq!(intervals(&cov), [(2, 5)]);
        assert_eq!(cov.coverage(), 3);
    }

    #[test]
    fn non_overlapping_in_order() {
        let mut cov = SetCoverage::new();
        cov.insert(0..3);
        cov.insert(5..8);
        assert_eq!(intervals(&cov), [(0, 3), (5, 8)]);
        assert_eq!(cov.coverage(), 6);
    }

    #[test]
    fn non_overlapping_reverse_order() {
        let mut cov = SetCoverage::new();
        cov.insert(5..8);
        cov.insert(0..3);
        assert_eq!(intervals(&cov), [(0, 3), (5, 8)]);
        assert_eq!(cov.coverage(), 6);
    }

    #[test]
    fn overlapping_intervals_merge() {
        let mut cov = SetCoverage::new();
        cov.insert(0..5);
        cov.insert(3..8);
        assert_eq!(intervals(&cov), [(0, 8)]);
        assert_eq!(cov.coverage(), 8);
    }

    #[test]
    fn adjacent_intervals_merge() {
        // [0,3) and [3,6) share no bases but should merge to [0,6).
        let mut cov = SetCoverage::new();
        cov.insert(0..3);
        cov.insert(3..6);
        assert_eq!(intervals(&cov), [(0, 6)]);
        assert_eq!(cov.coverage(), 6);
    }

    #[test]
    fn contained_interval_is_absorbed() {
        let mut cov = SetCoverage::new();
        cov.insert(0..10);
        cov.insert(3..7);
        assert_eq!(intervals(&cov), [(0, 10)]);
        assert_eq!(cov.coverage(), 10);
    }

    #[test]
    fn chain_merge_collapses_multiple() {
        let mut cov = SetCoverage::new();
        cov.insert(0..3);
        cov.insert(6..9);
        cov.insert(12..15);
        // Now insert something spanning all three.
        cov.insert(1..13);
        assert_eq!(intervals(&cov), [(0, 15)]);
        assert_eq!(cov.coverage(), 15);
    }

    #[test]
    fn start_and_end_accessors() {
        let mut cov = SetCoverage::new();
        assert_eq!(cov.start(), None);
        assert_eq!(cov.end(), None);
        cov.insert(4..7);
        cov.insert(10..14);
        assert_eq!(cov.start(), Some(4));
        assert_eq!(cov.end(), Some(14));
    }

    #[test]
    fn len_counts_disjoint_intervals() {
        let mut cov = SetCoverage::new();
        cov.insert(0..1);
        cov.insert(2..3);
        cov.insert(4..5);
        assert_eq!(cov.len(), 3);
        cov.insert(1..4); // merges all three
        assert_eq!(cov.len(), 1);
    }

    #[test]
    #[should_panic]
    fn empty_range_panics() {
        let mut cov = SetCoverage::new();
        cov.insert(5..5);
    }

    #[test]
    #[should_panic]
    fn inverted_range_panics() {
        let mut cov = SetCoverage::new();
        cov.insert(5..3);
    }
}
