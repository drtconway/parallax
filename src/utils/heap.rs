#![allow(dead_code)]

pub enum HeapOrdering {
    Min,
    Max,
}

pub trait Heapable {
    type Item;

    const ORDERING: HeapOrdering;

    // Compare two items, always returning the natural ordering.
    fn cmp(lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering;

    // Compare two items according to the heap's ordering.
    fn in_order(lhs: &Self::Item, rhs: &Self::Item) -> bool {
        match Self::ORDERING {
            HeapOrdering::Max => Self::cmp(lhs, rhs) == std::cmp::Ordering::Greater,
            HeapOrdering::Min => Self::cmp(lhs, rhs) == std::cmp::Ordering::Less,
        }
    }
}

pub struct Heap<H: Heapable> {
    data: Vec<H::Item>,
}

impl<H: Heapable> Heap<H> {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn push(&mut self, item: H::Item) {
        self.data.push(item);
        self.upheap(self.data.len() - 1);
    }

    pub fn pop(&mut self) -> Option<H::Item> {
        if self.data.is_empty() {
            return None;
        }
        let n = self.data.len() - 1;
        self.data.swap(0, n);
        let item = self.data.pop();
        self.downheap(0);
        item
    }

    fn heapify(&mut self) {
        let len = self.data.len();
        for i in (0..len / 2).rev() {
            self.downheap(i);
        }
    }

    fn upheap(&mut self, mut idx: usize) {

        while idx > 0 {
            let parent = (idx - 1) / 2;
            let should_swap = H::in_order(&self.data[idx], &self.data[parent]);
            if should_swap {
                self.data.swap(idx, parent);
                idx = parent;
            } else {
                break;
            }
        }
    }

    fn downheap(&mut self, mut idx: usize) {

        let len = self.data.len();
        loop {
            let left = 2 * idx + 1;
            let right = 2 * idx + 2;
            let mut best = idx;

            if left < len && H::in_order(&self.data[left], &self.data[best]) {
                best = left;
            }
            if right < len && H::in_order(&self.data[right], &self.data[best]) {
                best = right;
            }

            if best == idx {
                break;
            }
            self.data.swap(idx, best);
            idx = best;
        }
    }
}

impl<H: Heapable> From<Vec<<H as Heapable>::Item>> for Heap<H> {
    fn from(data: Vec<<H as Heapable>::Item>) -> Self {
        let mut heap = Self { data };
        heap.heapify();
        heap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Max heap for integers
    struct MaxHeap;
    impl Heapable for MaxHeap {
        type Item = i32;
        const ORDERING: HeapOrdering = HeapOrdering::Max;
        fn cmp(lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
            lhs.cmp(rhs)
        }
    }

    // Min heap for integers
    struct MinHeap;
    impl Heapable for MinHeap {
        type Item = i32;
        const ORDERING: HeapOrdering = HeapOrdering::Min;
        fn cmp(lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
            lhs.cmp(rhs)
        }
    }

    #[test]
    fn test_new_heap_is_empty() {
        let heap: Heap<MaxHeap> = Heap::new();
        assert_eq!(heap.len(), 0);
    }

    #[test]
    fn test_pop_empty_heap_returns_none() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_single_element() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        heap.push(42);
        assert_eq!(heap.len(), 1);
        assert_eq!(heap.pop(), Some(42));
        assert_eq!(heap.len(), 0);
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_max_heap_ordering() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        heap.push(3);
        heap.push(1);
        heap.push(4);
        heap.push(1);
        heap.push(5);
        heap.push(9);
        heap.push(2);

        assert_eq!(heap.len(), 7);
        assert_eq!(heap.pop(), Some(9));
        assert_eq!(heap.pop(), Some(5));
        assert_eq!(heap.pop(), Some(4));
        assert_eq!(heap.pop(), Some(3));
        assert_eq!(heap.pop(), Some(2));
        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_min_heap_ordering() {
        let mut heap: Heap<MinHeap> = Heap::new();
        heap.push(3);
        heap.push(1);
        heap.push(4);
        heap.push(1);
        heap.push(5);
        heap.push(9);
        heap.push(2);

        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), Some(2));
        assert_eq!(heap.pop(), Some(3));
        assert_eq!(heap.pop(), Some(4));
        assert_eq!(heap.pop(), Some(5));
        assert_eq!(heap.pop(), Some(9));
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_already_sorted_ascending() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        for i in 1..=10 {
            heap.push(i);
        }

        for i in (1..=10).rev() {
            assert_eq!(heap.pop(), Some(i));
        }
    }

    #[test]
    fn test_already_sorted_descending() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        for i in (1..=10).rev() {
            heap.push(i);
        }

        for i in (1..=10).rev() {
            assert_eq!(heap.pop(), Some(i));
        }
    }

    #[test]
    fn test_duplicate_values() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        heap.push(5);
        heap.push(5);
        heap.push(5);
        heap.push(5);

        assert_eq!(heap.len(), 4);
        assert_eq!(heap.pop(), Some(5));
        assert_eq!(heap.pop(), Some(5));
        assert_eq!(heap.pop(), Some(5));
        assert_eq!(heap.pop(), Some(5));
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_negative_numbers() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        heap.push(-5);
        heap.push(-1);
        heap.push(-10);
        heap.push(-3);

        assert_eq!(heap.pop(), Some(-1));
        assert_eq!(heap.pop(), Some(-3));
        assert_eq!(heap.pop(), Some(-5));
        assert_eq!(heap.pop(), Some(-10));
    }

    #[test]
    fn test_interleaved_push_pop() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        heap.push(5);
        heap.push(3);
        assert_eq!(heap.pop(), Some(5));
        heap.push(7);
        heap.push(1);
        assert_eq!(heap.pop(), Some(7));
        assert_eq!(heap.pop(), Some(3));
        heap.push(9);
        assert_eq!(heap.pop(), Some(9));
        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_large_dataset() {
        let values: Vec<i32> = (0..1000).collect();

        // Create heap from vector (in reverse order)
        let reversed: Vec<i32> = values.iter().rev().copied().collect();
        let mut heap: Heap<MinHeap> = Heap::from(reversed);

        // Should come out in ascending order
        for expected in values.iter() {
            assert_eq!(heap.pop(), Some(*expected));
        }
    }

    #[test]
    fn test_with_custom_type() {
        #[derive(Debug, PartialEq, Eq)]
        struct Priority {
            priority: u32,
            value: String,
        }

        struct PriorityHeap;
        impl Heapable for PriorityHeap {
            type Item = Priority;
            const ORDERING: HeapOrdering = HeapOrdering::Max;
            fn cmp(lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
                lhs.priority.cmp(&rhs.priority)
            }
        }

        let mut heap: Heap<PriorityHeap> = Heap::new();
        heap.push(Priority {
            priority: 5,
            value: "medium".to_string(),
        });
        heap.push(Priority {
            priority: 10,
            value: "high".to_string(),
        });
        heap.push(Priority {
            priority: 1,
            value: "low".to_string(),
        });

        assert_eq!(heap.pop().unwrap().value, "high");
        assert_eq!(heap.pop().unwrap().value, "medium");
        assert_eq!(heap.pop().unwrap().value, "low");
    }

    #[test]
    fn test_zero_value() {
        let mut heap: Heap<MaxHeap> = Heap::new();
        heap.push(0);
        heap.push(-1);
        heap.push(1);

        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), Some(0));
        assert_eq!(heap.pop(), Some(-1));
    }
}
