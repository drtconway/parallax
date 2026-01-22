
/// A set of non-overlapping ranges
pub struct RangeSet {
    ranges: Vec<(usize, usize)>,
    sorted_end: usize,
}

impl RangeSet {
    pub fn new() -> Self {
        RangeSet {
            ranges: Vec::new(),
            sorted_end: 0,
        }
    }

    pub fn add_range(&mut self, start: usize, end: usize) {
        assert!(start < end, "Invalid range: start must be less than end");
        assert!(
            !self.contains_overlap(start, end),
            "Range overlaps with existing ranges"
        );
        self.ranges.push((start, end));
        let num_unsorted = self.ranges.len() - self.sorted_end;
        if num_unsorted * num_unsorted > self.ranges.len() {
            self.ranges.sort_unstable();
            self.sorted_end = self.ranges.len();
        }
    }

    #[allow(dead_code)]
    pub fn contains(&self, value: usize) -> bool {
        let sorted_part = &self.ranges[..self.sorted_end];
        let j = sorted_part.partition_point(|(start, end)| *end <= value);
        if j < sorted_part.len() {
            let (start, end) = sorted_part[j];
            if value >= start && value < end {
                return true;
            }
        }

        for (start, end) in &self.ranges[self.sorted_end..] {
            if value >= *start && value < *end {
                return true;
            }
        }
        false
    }

    pub fn contains_overlap(&self, start: usize, end: usize) -> bool {
        let sorted_part = &self.ranges[..self.sorted_end];
        let j = sorted_part.partition_point(|(_, r_end)| *r_end <= start);
        if j < sorted_part.len() {
            let (r_start, r_end) = sorted_part[j];
            if start < r_end && end > r_start {
                return true;
            }
        }

        for &(r_start, r_end) in &self.ranges[self.sorted_end..] {
            if start < r_end && end > r_start {
                return true;
            }
        }
        false
    }
}

    #[cfg(test)]
    mod tests {
        use super::RangeSet;

        #[test]
        fn test_empty_contains_false() {
            let rs = RangeSet::new();
            assert!(!rs.contains(0));
            assert!(!rs.contains(10));
        }

        #[test]
        fn test_add_and_contains_single_range() {
            let mut rs = RangeSet::new();
            rs.add_range(5, 10);

            assert!(!rs.contains(4));
            assert!(rs.contains(5));
            assert!(rs.contains(9));
            assert!(!rs.contains(10));
        }

        #[test]
        fn test_contains_across_sorted_and_unsorted_parts() {
            let mut rs = RangeSet::new();

            // Add two ranges; these stay in the unsorted region initially.
            rs.add_range(20, 25);
            rs.add_range(5, 8);

            assert!(rs.contains(6));
            assert!(rs.contains(22));
            assert!(!rs.contains(8));
            assert!(!rs.contains(19));

            // Add a third range to trigger sort (num_unsorted^2 > len).
            rs.add_range(12, 15);

            // After sorting, contains should still work for all ranges.
            assert!(rs.contains(6));
            assert!(rs.contains(22));
            assert!(rs.contains(12));
            assert!(!rs.contains(9));
            assert!(!rs.contains(25));
        }

        #[test]
        fn test_contains_overlap_detection() {
            let mut rs = RangeSet::new();
            rs.add_range(0, 3);
            rs.add_range(10, 12);

            assert!(rs.contains_overlap(2, 5));
            assert!(rs.contains_overlap(11, 13));
            assert!(rs.contains_overlap(0, 1));
            assert!(rs.contains_overlap(9, 11));
            assert!(!rs.contains_overlap(3, 10));
            assert!(!rs.contains_overlap(12, 15));
        }

        #[test]
        #[should_panic(expected = "start must be less than end")]
        fn test_add_range_invalid_panics() {
            let mut rs = RangeSet::new();
            rs.add_range(5, 5);
        }

        #[test]
        #[should_panic(expected = "overlaps")]
        fn test_add_range_overlap_panics() {
            let mut rs = RangeSet::new();
            rs.add_range(5, 10);
            rs.add_range(9, 12);
        }
    }
