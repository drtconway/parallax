#![allow(dead_code)] // This module is used in some configurations but not others

pub trait Joinable: Sized {

    fn range(&self) -> (usize, usize);

    /// Returns true if self ends strictly before other starts (even with gap tolerance).
    /// This means there's no possible overlap between self and other.
    fn is_before<Other: Joinable>(&self, other: &Other, gap: usize) -> bool {
        let (_, self_end) = self.range();
        let (other_start, _) = other.range();
        // Use < not <= so that boundary cases (self_end + gap == other_start) are checked for overlap
        self_end.saturating_add(gap) < other_start
    }
}

pub fn sorted_join<T: Joinable, U: Joinable, F: Fn(&T, &U) -> bool>(
    left: &[T],
    right: &[U],
    gap: usize,
    pred: F,
) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    let mut j_start = 0;

    for (i, l_item) in left.iter().enumerate() {
        let (l_start, _) = l_item.range();

        // Advance j_start past right items that are completely before l_item
        // (and thus before all future left items, since left is sorted)
        while j_start < right.len() && right[j_start].is_before(l_item, gap) {
            j_start += 1;
        }

        for j in j_start..right.len() {
            let r_item = &right[j];

            // If l_item is completely before r_item, no more matches for this l_item
            // (but don't update j_start - next left item might still match r_item)
            if l_item.is_before(r_item, gap) {
                break;
            }

            // Check for overlap with gap tolerance
            let (_, l_end) = l_item.range();
            let (r_start, r_end) = r_item.range();
            if l_end.saturating_add(gap) >= r_start && r_end.saturating_add(gap) >= l_start {
                if pred(l_item, r_item) {
                    result.push((i, j));
                }
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Simple interval struct for testing
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct Interval {
        start: usize,
        end: usize,
        label: char,
    }

    impl Interval {
        fn new(start: usize, end: usize, label: char) -> Self {
            Self { start, end, label }
        }
    }

    impl Joinable for Interval {
        fn range(&self) -> (usize, usize) {
            (self.start, self.end)
        }
    }

    #[test]
    fn test_empty_inputs() {
        let left: Vec<Interval> = vec![];
        let right: Vec<Interval> = vec![];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert!(result.is_empty());
    }

    #[test]
    fn test_empty_left() {
        let left: Vec<Interval> = vec![];
        let right = vec![Interval::new(0, 10, 'a')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert!(result.is_empty());
    }

    #[test]
    fn test_empty_right() {
        let left = vec![Interval::new(0, 10, 'a')];
        let right: Vec<Interval> = vec![];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert!(result.is_empty());
    }

    #[test]
    fn test_exact_overlap() {
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(10, 20, 'b')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_partial_overlap() {
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(15, 25, 'b')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_adjacent_no_gap() {
        // Adjacent intervals [10,20) and [20,30) DO overlap with gap=0
        // because they touch at boundary: 20 + 0 >= 20 AND 30 + 0 >= 10
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(20, 30, 'b')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_truly_disjoint() {
        // Intervals with a gap between them should NOT overlap with gap=0
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(21, 30, 'b')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert!(result.is_empty());
    }

    #[test]
    fn test_adjacent_with_gap() {
        // Adjacent intervals [10,20) and [20,30) SHOULD overlap with gap=1
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(20, 30, 'b')];
        let result = sorted_join(&left, &right, 1, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_gap_tolerance() {
        // Intervals with a gap of 5 between them: [10,20) and [25,35)
        // For overlap with gap tolerance: l_end + gap >= r_start AND r_end + gap >= l_start
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(25, 35, 'b')];

        // gap=4: 20 + 4 = 24 < 25, so no match
        let result = sorted_join(&left, &right, 4, |_, _| true);
        assert!(result.is_empty());

        // gap=5: 20 + 5 = 25 >= 25, should match
        let result = sorted_join(&left, &right, 5, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);

        // gap=10 should also match
        let result = sorted_join(&left, &right, 10, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_no_overlap_too_far() {
        let left = vec![Interval::new(10, 20, 'a')];
        let right = vec![Interval::new(100, 110, 'b')];
        let result = sorted_join(&left, &right, 10, |_, _| true);
        assert!(result.is_empty());
    }

    #[test]
    fn test_predicate_filters() {
        let left = vec![
            Interval::new(10, 20, 'a'),
            Interval::new(30, 40, 'b'),
        ];
        let right = vec![
            Interval::new(15, 25, 'x'),
            Interval::new(35, 45, 'y'),
        ];

        // Only match if labels satisfy predicate
        let result = sorted_join(&left, &right, 0, |l, r| l.label == 'a' && r.label == 'x');
        assert_eq!(result, vec![(0, 0)]);

        // Match both with always-true predicate
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0), (1, 1)]);

        // Match none with always-false predicate
        let result = sorted_join(&left, &right, 0, |_, _| false);
        assert!(result.is_empty());
    }

    #[test]
    fn test_one_to_many() {
        // One left interval overlaps multiple right intervals
        let left = vec![Interval::new(10, 50, 'a')];
        let right = vec![
            Interval::new(5, 15, 'x'),
            Interval::new(20, 30, 'y'),
            Interval::new(40, 55, 'z'),
        ];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0), (0, 1), (0, 2)]);
    }

    #[test]
    fn test_many_to_one() {
        // Multiple left intervals overlap one right interval
        let left = vec![
            Interval::new(5, 15, 'a'),
            Interval::new(20, 30, 'b'),
            Interval::new(40, 55, 'c'),
        ];
        let right = vec![Interval::new(10, 50, 'x')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0), (1, 0), (2, 0)]);
    }

    #[test]
    fn test_many_to_many() {
        let left = vec![
            Interval::new(0, 10, 'a'),
            Interval::new(20, 30, 'b'),
            Interval::new(40, 50, 'c'),
        ];
        let right = vec![
            Interval::new(5, 25, 'x'),
            Interval::new(25, 45, 'y'),
        ];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        // 'a' [0,10) overlaps 'x' [5,25)
        // 'b' [20,30) overlaps 'x' [5,25) and 'y' [25,45)
        // 'c' [40,50) overlaps 'y' [25,45)
        assert_eq!(result, vec![(0, 0), (1, 0), (1, 1), (2, 1)]);
    }

    #[test]
    fn test_contained_intervals() {
        // One interval completely contains another
        let left = vec![Interval::new(10, 50, 'a')];
        let right = vec![Interval::new(20, 30, 'x')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);

        // Swap sides
        let left = vec![Interval::new(20, 30, 'a')];
        let right = vec![Interval::new(10, 50, 'x')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_sorted_input_required() {
        // Intervals must be sorted by start position for correct results
        let left = vec![
            Interval::new(0, 10, 'a'),
            Interval::new(20, 30, 'b'),
            Interval::new(50, 60, 'c'),
        ];
        let right = vec![
            Interval::new(5, 15, 'x'),
            Interval::new(55, 65, 'y'),
        ];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0), (2, 1)]);
    }

    #[test]
    fn test_is_before() {
        let a = Interval::new(10, 20, 'a');
        let b = Interval::new(30, 40, 'b');

        // a ends at 20, b starts at 30 -> gap of 10
        // is_before uses strict < so boundary case returns false
        assert!(a.is_before(&b, 0));   // 20 + 0 = 20 < 30
        assert!(a.is_before(&b, 9));   // 20 + 9 = 29 < 30
        assert!(!a.is_before(&b, 10)); // 20 + 10 = 30 is NOT < 30
        assert!(!a.is_before(&b, 11)); // 20 + 11 = 31 > 30
    }

    #[test]
    fn test_point_intervals() {
        // Zero-length intervals (points) at the same position
        // [10,10) and [10,10):
        // is_before uses <: 10 + 0 < 10 is false, so we check overlap
        // overlap: 10 + 0 >= 10 AND 10 + 0 >= 10 -> true
        let left = vec![Interval::new(10, 10, 'a')];
        let right = vec![Interval::new(10, 10, 'x')];
        let result = sorted_join(&left, &right, 0, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }

    #[test]
    fn test_large_gap_value() {
        let left = vec![Interval::new(0, 10, 'a')];
        let right = vec![Interval::new(1000, 1010, 'x')];

        // With large enough gap, they should match
        let result = sorted_join(&left, &right, 1000, |_, _| true);
        assert_eq!(result, vec![(0, 0)]);
    }
}