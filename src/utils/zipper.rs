#![allow(dead_code)]
use std::{collections::VecDeque, ops::{Index, IndexMut}};

/// A simple zipper structure that allows local modifications to
/// the interior of a vector by maintaining a front and back vector.
///
/// The zipper allows moving elements from one side to the other,
/// enabling efficient insertions and deletions in the middle of the sequence,
/// while keeping the overall order intact.
#[derive(Debug, Clone)]
pub struct Zipper<T> {
    left: VecDeque<T>,
    right: VecDeque<T>,
}

impl<T> Zipper<T> {
    /// Create a new empty zipper.
    pub fn new() -> Self {
        Self {
            left: VecDeque::new(),
            right: VecDeque::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.left.len() + self.right.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn clear(&mut self) {
        self.left.clear();
        self.right.clear();
    }

    pub fn push_front(&mut self, item: T) {
        self.left.push_front(item);
    }

    pub fn push_back(&mut self, item: T) {
        self.right.push_back(item);
    }

    pub fn pop_front(&mut self) -> Option<T> {
        if self.left.is_empty() {
            self.move_right_to_left();
        }
        self.left.pop_front()
    }

    pub fn pop_back(&mut self) -> Option<T> {
        if self.right.is_empty() {
            self.move_left_to_right();
        }
        self.right.pop_back()
    }

    /// Insert an element at the given index.
    pub fn insert(&mut self, index: usize, item: T) {
        assert!(index <= self.len(), "Index out of bounds");
        while self.left.len() < index {
            self.move_right_to_left();
        }
        while self.left.len() > index {
            self.move_left_to_right();
        }
        self.left.push_back(item);
    }

    /// Remove and return the element at the given index.
    pub fn remove(&mut self, index: usize) -> Option<T> {
        assert!(index < self.len(), "Index out of bounds");
        // Position so that left has exactly `index` elements
        // The element at `index` will be at the front of right
        while self.left.len() < index {
            self.move_right_to_left();
        }
        while self.left.len() > index {
            self.move_left_to_right();
        }
        self.right.pop_front()
    }

    pub fn partition_point<P: FnMut(&T) -> bool>(&self, mut predicate: P) -> usize {
        // Perform binary search over the combined left and right vectors.
        let mut left = 0;
        let mut right = self.len();
        while left < right {
            let mid = (left + right) / 2;
            if predicate(&self[mid]) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        left
    }

    /// Move an element from the right side to the left side.
    fn move_right_to_left(&mut self)  {
        self.right.pop_front().map(|item| {
            self.left.push_back(item);
        });
    }

    /// Move an element from the left side to the right side.
    fn move_left_to_right(&mut self) {
        self.left.pop_back().map(|item| {
            self.right.push_front(item);
        });
    }

    /// Get the current state as a vector.
    pub fn to_vec(self) -> Vec<T> {
        let mut result: Vec<T> = self.left.into_iter().collect();
        result.extend(self.right.into_iter());
        result
    }
}

impl<T> Index<usize> for Zipper<T>
{
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        if index < self.left.len() {
            &self.left[index]
        } else {
            &self.right[index - self.left.len()]
        }
    }
}

impl<T> IndexMut<usize> for Zipper<T>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index < self.left.len() {
            &mut self.left[index]
        } else {
            &mut self.right[index - self.left.len()]
        }
    }
}

impl<T> Default for Zipper<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> From<Vec<T>> for Zipper<T> {
    fn from(data: Vec<T>) -> Self {
        Self {
            left: VecDeque::new(),
            right: VecDeque::from(data),
        }
    }
}

impl<T> FromIterator<T> for Zipper<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let data: Vec<T> = iter.into_iter().collect();
        Self::from(data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_creates_empty() {
        let zipper: Zipper<i32> = Zipper::new();
        assert!(zipper.is_empty());
        assert_eq!(zipper.len(), 0);
    }

    #[test]
    fn test_default_creates_empty() {
        let zipper: Zipper<i32> = Zipper::default();
        assert!(zipper.is_empty());
        assert_eq!(zipper.len(), 0);
    }

    #[test]
    fn test_from_vec() {
        let zipper: Zipper<i32> = Zipper::from(vec![1, 2, 3, 4, 5]);
        assert_eq!(zipper.len(), 5);
        assert!(!zipper.is_empty());
    }

    #[test]
    fn test_from_empty_vec() {
        let zipper: Zipper<i32> = Zipper::from(vec![]);
        assert!(zipper.is_empty());
        assert_eq!(zipper.len(), 0);
    }

    #[test]
    fn test_index_access() {
        let zipper = Zipper::from(vec![10, 20, 30, 40]);
        assert_eq!(zipper[0], 10);
        assert_eq!(zipper[1], 20);
        assert_eq!(zipper[2], 30);
        assert_eq!(zipper[3], 40);
    }

    #[test]
    fn test_index_mut() {
        let mut zipper = Zipper::from(vec![1, 2, 3]);
        zipper[1] = 99;
        assert_eq!(zipper[1], 99);
        assert_eq!(zipper.to_vec(), vec![1, 99, 3]);
    }

    #[test]
    fn test_push_front() {
        let mut zipper = Zipper::from(vec![2, 3, 4]);
        zipper.push_front(1);
        assert_eq!(zipper.len(), 4);
        assert_eq!(zipper[0], 1);
        assert_eq!(zipper.to_vec(), vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_push_back() {
        let mut zipper = Zipper::from(vec![1, 2, 3]);
        zipper.push_back(4);
        assert_eq!(zipper.len(), 4);
        assert_eq!(zipper[3], 4);
        assert_eq!(zipper.to_vec(), vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_pop_front() {
        let mut zipper = Zipper::from(vec![1, 2, 3]);
        assert_eq!(zipper.pop_front(), Some(1));
        assert_eq!(zipper.len(), 2);
        assert_eq!(zipper.to_vec(), vec![2, 3]);
    }

    #[test]
    fn test_pop_front_empty() {
        let mut zipper: Zipper<i32> = Zipper::new();
        assert_eq!(zipper.pop_front(), None);
    }

    #[test]
    fn test_pop_back() {
        let mut zipper = Zipper::from(vec![1, 2, 3]);
        assert_eq!(zipper.pop_back(), Some(3));
        assert_eq!(zipper.len(), 2);
        assert_eq!(zipper.to_vec(), vec![1, 2]);
    }

    #[test]
    fn test_pop_back_empty() {
        let mut zipper: Zipper<i32> = Zipper::new();
        assert_eq!(zipper.pop_back(), None);
    }

    #[test]
    fn test_insert_at_beginning() {
        let mut zipper = Zipper::from(vec![2, 3, 4]);
        zipper.insert(0, 1);
        assert_eq!(zipper.to_vec(), vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_insert_at_end() {
        let mut zipper = Zipper::from(vec![1, 2, 3]);
        zipper.insert(3, 4);
        assert_eq!(zipper.to_vec(), vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_insert_in_middle() {
        let mut zipper = Zipper::from(vec![1, 2, 4, 5]);
        zipper.insert(2, 3);
        assert_eq!(zipper.to_vec(), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_insert_into_empty() {
        let mut zipper: Zipper<i32> = Zipper::new();
        zipper.insert(0, 42);
        assert_eq!(zipper.to_vec(), vec![42]);
    }

    #[test]
    fn test_remove_from_beginning() {
        let mut zipper = Zipper::from(vec![1, 2, 3, 4]);
        assert_eq!(zipper.remove(0), Some(1));
        assert_eq!(zipper.to_vec(), vec![2, 3, 4]);
    }

    #[test]
    fn test_remove_from_end() {
        let mut zipper = Zipper::from(vec![1, 2, 3, 4]);
        assert_eq!(zipper.remove(3), Some(4));
        assert_eq!(zipper.to_vec(), vec![1, 2, 3]);
    }

    #[test]
    fn test_remove_from_middle() {
        let mut zipper = Zipper::from(vec![1, 2, 3, 4, 5]);
        assert_eq!(zipper.remove(2), Some(3));
        assert_eq!(zipper.to_vec(), vec![1, 2, 4, 5]);
    }

    #[test]
    fn test_multiple_inserts_and_removes() {
        let mut zipper = Zipper::from(vec![1, 5]);
        zipper.insert(1, 2);
        zipper.insert(2, 3);
        zipper.insert(3, 4);
        assert_eq!(zipper.clone().to_vec(), vec![1, 2, 3, 4, 5]);

        zipper.remove(2); // remove 3
        zipper.remove(1); // remove 2
        assert_eq!(zipper.to_vec(), vec![1, 4, 5]);
    }

    #[test]
    fn test_clear() {
        let mut zipper = Zipper::from(vec![1, 2, 3]);
        zipper.insert(1, 10); // move some to left
        zipper.clear();
        assert!(zipper.is_empty());
        assert_eq!(zipper.len(), 0);
    }

    #[test]
    fn test_partition_point_all_true() {
        let zipper = Zipper::from(vec![1, 2, 3, 4, 5]);
        let point = zipper.partition_point(|&x| x < 10);
        assert_eq!(point, 5);
    }

    #[test]
    fn test_partition_point_all_false() {
        let zipper = Zipper::from(vec![1, 2, 3, 4, 5]);
        let point = zipper.partition_point(|&x| x < 0);
        assert_eq!(point, 0);
    }

    #[test]
    fn test_partition_point_middle() {
        let zipper = Zipper::from(vec![1, 2, 3, 4, 5]);
        let point = zipper.partition_point(|&x| x < 3);
        assert_eq!(point, 2); // elements 1, 2 satisfy x < 3
    }

    #[test]
    fn test_partition_point_empty() {
        let zipper: Zipper<i32> = Zipper::new();
        let point = zipper.partition_point(|&x| x < 5);
        assert_eq!(point, 0);
    }

    #[test]
    fn test_partition_point_after_modifications() {
        let mut zipper = Zipper::from(vec![10, 20, 40, 50]);
        zipper.insert(2, 30); // [10, 20, 30, 40, 50]
        let point = zipper.partition_point(|&x| x <= 30);
        assert_eq!(point, 3);
    }

    #[test]
    fn test_index_after_insert() {
        let mut zipper = Zipper::from(vec![1, 2, 4, 5]);
        zipper.insert(2, 3);
        // After insert, some elements are in left, some in right
        assert_eq!(zipper[0], 1);
        assert_eq!(zipper[1], 2);
        assert_eq!(zipper[2], 3);
        assert_eq!(zipper[3], 4);
        assert_eq!(zipper[4], 5);
    }

    #[test]
    fn test_index_after_remove() {
        let mut zipper = Zipper::from(vec![1, 2, 3, 4, 5]);
        zipper.remove(2); // remove 3
        assert_eq!(zipper[0], 1);
        assert_eq!(zipper[1], 2);
        assert_eq!(zipper[2], 4);
        assert_eq!(zipper[3], 5);
    }

    #[test]
    fn test_sequential_operations_preserve_order() {
        let mut zipper = Zipper::from(vec![5, 10, 15, 20, 25]);

        // Insert at various positions
        zipper.insert(0, 0);    // [0, 5, 10, 15, 20, 25]
        zipper.insert(3, 12);   // [0, 5, 10, 12, 15, 20, 25]
        zipper.insert(7, 30);   // [0, 5, 10, 12, 15, 20, 25, 30]

        // Remove some
        zipper.remove(1);       // [0, 10, 12, 15, 20, 25, 30]
        zipper.remove(4);       // [0, 10, 12, 15, 25, 30]

        assert_eq!(zipper.to_vec(), vec![0, 10, 12, 15, 25, 30]);
    }

    #[test]
    fn test_alternating_push_pop() {
        let mut zipper = Zipper::from(vec![3]);
        zipper.push_front(2);
        zipper.push_front(1);
        zipper.push_back(4);
        zipper.push_back(5);
        assert_eq!(zipper.clone().to_vec(), vec![1, 2, 3, 4, 5]);

        assert_eq!(zipper.pop_front(), Some(1));
        assert_eq!(zipper.pop_back(), Some(5));
        assert_eq!(zipper.to_vec(), vec![2, 3, 4]);
    }

    #[test]
    fn test_with_strings() {
        let mut zipper = Zipper::from(vec!["b".to_string(), "d".to_string()]);
        zipper.push_front("a".to_string());
        zipper.insert(2, "c".to_string());
        zipper.push_back("e".to_string());

        let result: Vec<String> = zipper.to_vec();
        assert_eq!(result, vec!["a", "b", "c", "d", "e"]);
    }

    #[test]
    fn test_into_syntax() {
        // Test that .into() works
        let data = vec![1, 2, 3];
        let zipper: Zipper<i32> = data.into();
        assert_eq!(zipper.len(), 3);
    }
}