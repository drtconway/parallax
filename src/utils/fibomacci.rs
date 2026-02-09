//! Fibonacci Heap implementation following the API design of the heap module.
//!
//! A Fibonacci heap provides:
//! - O(1) amortized insert
//! - O(1) amortized find-min/max
//! - O(log n) amortized extract-min/max
//! - O(1) amortized decrease-key (with handles)
//! - O(log n) amortized delete (with handles)
//!
//! This makes it ideal for algorithms like Dijkstra's or Prim's where
//! decrease-key operations are frequent.

#![allow(dead_code)]

use super::heap::Heapable;
use std::cell::Cell;
use std::ptr::NonNull;

/// A handle to a node in the Fibonacci heap, allowing O(1) decrease-key operations.
#[derive(Clone, Copy)]
pub struct FibHandle {
    ptr: NonNull<FibNodeInner>,
}

// Safety: FibHandle is only valid while the heap exists and the node hasn't been removed.
// The caller must ensure proper lifetime management.
unsafe impl Send for FibHandle {}

struct FibNodeInner {
    // Intrusive doubly-linked circular list pointers
    left: Cell<NonNull<FibNodeInner>>,
    right: Cell<NonNull<FibNodeInner>>,
    parent: Cell<Option<NonNull<FibNodeInner>>>,
    child: Cell<Option<NonNull<FibNodeInner>>>,
    degree: Cell<usize>,
    marked: Cell<bool>,
}

/// A node in the Fibonacci heap containing the actual item.
#[repr(C)]
struct FibNode<T> {
    inner: FibNodeInner,
    item: T,
}

impl<T> FibNode<T> {
    fn new(item: T) -> Box<Self> {
        let node = Box::new(FibNode {
            inner: FibNodeInner {
                left: Cell::new(NonNull::dangling()),
                right: Cell::new(NonNull::dangling()),
                parent: Cell::new(None),
                child: Cell::new(None),
                degree: Cell::new(0),
                marked: Cell::new(false),
            },
            item,
        });

        // Initialize circular list pointers to point to self
        let ptr = NonNull::from(&node.inner);
        node.inner.left.set(ptr);
        node.inner.right.set(ptr);

        node
    }

    fn inner_ptr(&self) -> NonNull<FibNodeInner> {
        NonNull::from(&self.inner)
    }
}

/// A Fibonacci heap with the same trait-based API as the binary heap.
pub struct FibHeap<H: Heapable> {
    /// Pointer to the minimum (or maximum, depending on ordering) node
    root: Option<NonNull<FibNodeInner>>,
    /// Number of nodes in the heap
    len: usize,
    /// Phantom data for the Heapable type
    _marker: std::marker::PhantomData<H>,
}

impl<H: Heapable> FibHeap<H> {
    /// Create a new empty Fibonacci heap.
    pub fn new() -> Self {
        Self {
            root: None,
            len: 0,
            _marker: std::marker::PhantomData,
        }
    }

    /// Returns the number of items in the heap.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns true if the heap is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Insert an item into the heap. Returns a handle for decrease-key operations.
    pub fn push(&mut self, item: H::Item) -> FibHandle {
        let node = FibNode::<H::Item>::new(item);
        let ptr = Box::into_raw(node);

        // Safety: ptr is valid, we just allocated it
        let inner_ptr = unsafe { NonNull::new_unchecked(&raw mut (*ptr).inner) };

        self.insert_into_root_list(inner_ptr);

        // Update min pointer if needed
        match self.root {
            None => self.root = Some(inner_ptr),
            Some(min_ptr) => {
                if self.is_higher_priority(inner_ptr, min_ptr) {
                    self.root = Some(inner_ptr);
                }
            }
        }

        self.len += 1;

        FibHandle { ptr: inner_ptr }
    }

    /// Peek at the top item without removing it.
    pub fn peek(&self) -> Option<&H::Item> {
        self.root.map(|ptr| {
            // Safety: ptr is valid while heap exists
            unsafe { &self.node_from_inner(ptr).item }
        })
    }

    /// Remove and return the top item (min or max depending on ordering).
    pub fn pop(&mut self) -> Option<H::Item> {
        let min_ptr = self.root?;

        // Safety: min_ptr is valid
        unsafe {
            // Collect all children first (before modifying any pointers)
            let mut children = Vec::new();
            if let Some(child) = (*min_ptr.as_ptr()).child.get() {
                let mut current = child;
                loop {
                    children.push(current);
                    let next = (*current.as_ptr()).right.get();
                    if next == child {
                        break;
                    }
                    current = next;
                }
            }

            // Find a sibling root before removing min (need this for new root)
            let sibling = {
                let right = (*min_ptr.as_ptr()).right.get();
                if right != min_ptr {
                    Some(right)
                } else {
                    None
                }
            };

            // Remove min from root list
            self.remove_from_list(min_ptr);

            // Update self.root BEFORE adding children, so insert_into_root_list works correctly
            self.root = sibling;

            // Add children to root list
            for child in children {
                // Reset child's pointers to be self-referential first
                (*child.as_ptr()).left.set(child);
                (*child.as_ptr()).right.set(child);
                (*child.as_ptr()).parent.set(None);
                self.insert_into_root_list(child);
                // If this is the first node in root list, set it as root
                if self.root.is_none() {
                    self.root = Some(child);
                }
            }

            self.len -= 1;

            if self.len == 0 {
                self.root = None;
            } else if self.root.is_some() {
                self.consolidate();
            }

            // Reconstruct the Box and extract the item
            let node = Box::from_raw(self.node_ptr_from_inner(min_ptr));
            Some(node.item)
        }
    }

    /// Decrease (or increase for max-heap) the key of an item.
    /// The new_item must have higher priority than the current item.
    ///
    /// # Safety
    /// The handle must be valid (node still in heap, not yet popped).
    pub unsafe fn decrease_key(&mut self, handle: FibHandle, new_item: H::Item) {
        // Safety: caller guarantees handle is valid
        unsafe {
            let ptr = handle.ptr;
            let node = self.node_from_inner_mut(ptr);

            // Verify the new item has higher priority
            debug_assert!(
                H::in_order(&new_item, &node.item),
                "decrease_key called with lower priority item"
            );

            node.item = new_item;

            let parent = (*ptr.as_ptr()).parent.get();

            if let Some(parent_ptr) = parent {
                if self.is_higher_priority(ptr, parent_ptr) {
                    self.cut(ptr, parent_ptr);
                    self.cascading_cut(parent_ptr);
                }
            }

            // Update min if needed
            if let Some(min_ptr) = self.root {
                if self.is_higher_priority(ptr, min_ptr) {
                    self.root = Some(ptr);
                }
            }
        }
    }

    /// Delete a node from the heap.
    ///
    /// # Safety
    /// The handle must be valid (node still in heap, not yet popped).
    pub unsafe fn delete(&mut self, handle: FibHandle) -> H::Item {
        // Safety: caller guarantees handle is valid
        unsafe {
            let ptr = handle.ptr;

            // Cut and move to root if not already there
            if let Some(parent_ptr) = (*ptr.as_ptr()).parent.get() {
                self.cut(ptr, parent_ptr);
                self.cascading_cut(parent_ptr);
            }

            // Make this node the minimum
            self.root = Some(ptr);
        }

        // Pop it
        self.pop().expect("node was in heap")
    }

    // =========================================================================
    // Internal helpers
    // =========================================================================

    /// Check if node a has higher priority than node b.
    fn is_higher_priority(&self, a: NonNull<FibNodeInner>, b: NonNull<FibNodeInner>) -> bool {
        unsafe {
            let item_a = &self.node_from_inner(a).item;
            let item_b = &self.node_from_inner(b).item;
            H::in_order(item_a, item_b)
        }
    }

    /// Get a reference to the FibNode from its inner pointer.
    unsafe fn node_from_inner(&self, ptr: NonNull<FibNodeInner>) -> &FibNode<H::Item> {
        // Safety: The inner is at offset 0 in FibNode, caller guarantees ptr is valid
        unsafe { &*(ptr.as_ptr() as *const FibNode<H::Item>) }
    }

    /// Get a mutable reference to the FibNode from its inner pointer.
    unsafe fn node_from_inner_mut(&self, ptr: NonNull<FibNodeInner>) -> &mut FibNode<H::Item> {
        // Safety: caller guarantees ptr is valid and mutable access is allowed
        unsafe { &mut *(ptr.as_ptr() as *mut FibNode<H::Item>) }
    }

    /// Get the raw pointer to FibNode from inner pointer.
    unsafe fn node_ptr_from_inner(&self, ptr: NonNull<FibNodeInner>) -> *mut FibNode<H::Item> {
        ptr.as_ptr() as *mut FibNode<H::Item>
    }

    /// Insert a node into the root list.
    fn insert_into_root_list(&mut self, node: NonNull<FibNodeInner>) {
        unsafe {
            match self.root {
                None => {
                    // Empty root list - node points to itself
                    (*node.as_ptr()).left.set(node);
                    (*node.as_ptr()).right.set(node);
                }
                Some(root) => {
                    // Insert node to the left of root
                    let left = (*root.as_ptr()).left.get();
                    (*node.as_ptr()).left.set(left);
                    (*node.as_ptr()).right.set(root);
                    (*left.as_ptr()).right.set(node);
                    (*root.as_ptr()).left.set(node);
                }
            }
            (*node.as_ptr()).parent.set(None);
        }
    }

    /// Remove a node from its circular list.
    fn remove_from_list(&mut self, node: NonNull<FibNodeInner>) {
        unsafe {
            let left = (*node.as_ptr()).left.get();
            let right = (*node.as_ptr()).right.get();

            if left == node {
                // Node is alone in its list - nothing to do
                return;
            }

            (*left.as_ptr()).right.set(right);
            (*right.as_ptr()).left.set(left);

            // Reset node's pointers to itself
            (*node.as_ptr()).left.set(node);
            (*node.as_ptr()).right.set(node);
        }
    }

    /// Consolidate trees in the root list to have at most one tree of each degree.
    fn consolidate(&mut self) {
        // Maximum possible degree is log_phi(n) ≈ 1.44 * log2(n)
        // For practical purposes, 64 is more than enough
        const MAX_DEGREE: usize = 64;
        let mut degree_table: [Option<NonNull<FibNodeInner>>; MAX_DEGREE] = [None; MAX_DEGREE];

        // Collect all roots first (since we'll be modifying the list)
        let mut roots = Vec::new();
        if let Some(start) = self.root {
            unsafe {
                let mut current = start;
                loop {
                    roots.push(current);
                    current = (*current.as_ptr()).right.get();
                    if current == start {
                        break;
                    }
                }
            }
        }

        // Clear root list
        self.root = None;

        // Process each root
        for mut root in roots {
            unsafe {
                // Reset pointers since we removed from list
                (*root.as_ptr()).left.set(root);
                (*root.as_ptr()).right.set(root);

                let mut degree = (*root.as_ptr()).degree.get();

                while let Some(mut other) = degree_table[degree] {
                    degree_table[degree] = None;

                    // Link the two trees - higher priority becomes parent
                    if self.is_higher_priority(other, root) {
                        std::mem::swap(&mut root, &mut other);
                    }
                    self.link(other, root);

                    degree = (*root.as_ptr()).degree.get();
                }

                degree_table[degree] = Some(root);
            }
        }

        // Rebuild root list and find new min
        self.root = None;
        for slot in degree_table.iter() {
            if let Some(node) = *slot {
                self.insert_into_root_list(node);
                match self.root {
                    None => self.root = Some(node),
                    Some(min) => {
                        if self.is_higher_priority(node, min) {
                            self.root = Some(node);
                        }
                    }
                }
            }
        }
    }

    /// Make y a child of x.
    unsafe fn link(&mut self, y: NonNull<FibNodeInner>, x: NonNull<FibNodeInner>) {
        // Safety: caller guarantees both pointers are valid
        unsafe {
            // Remove y from root list
            self.remove_from_list(y);

            // Make y a child of x
            (*y.as_ptr()).parent.set(Some(x));
            (*y.as_ptr()).marked.set(false);

            match (*x.as_ptr()).child.get() {
                None => {
                    (*x.as_ptr()).child.set(Some(y));
                    (*y.as_ptr()).left.set(y);
                    (*y.as_ptr()).right.set(y);
                }
                Some(child) => {
                    // Insert y into child list
                    let left = (*child.as_ptr()).left.get();
                    (*y.as_ptr()).left.set(left);
                    (*y.as_ptr()).right.set(child);
                    (*left.as_ptr()).right.set(y);
                    (*child.as_ptr()).left.set(y);
                }
            }

            let deg = (*x.as_ptr()).degree.get();
            (*x.as_ptr()).degree.set(deg + 1);
        }
    }

    /// Cut x from its parent y and add to root list.
    unsafe fn cut(&mut self, x: NonNull<FibNodeInner>, y: NonNull<FibNodeInner>) {
        // Safety: caller guarantees both pointers are valid
        unsafe {
            // Remove x from child list of y
            let left = (*x.as_ptr()).left.get();
            let right = (*x.as_ptr()).right.get();

            if left == x {
                // x is the only child
                (*y.as_ptr()).child.set(None);
            } else {
                (*left.as_ptr()).right.set(right);
                (*right.as_ptr()).left.set(left);
                // If y's child pointer was x, update it
                if (*y.as_ptr()).child.get() == Some(x) {
                    (*y.as_ptr()).child.set(Some(right));
                }
            }

            let deg = (*y.as_ptr()).degree.get();
            (*y.as_ptr()).degree.set(deg.saturating_sub(1));

            // Add x to root list
            (*x.as_ptr()).left.set(x);
            (*x.as_ptr()).right.set(x);
            self.insert_into_root_list(x);

            (*x.as_ptr()).parent.set(None);
            (*x.as_ptr()).marked.set(false);
        }
    }

    /// Perform cascading cut up the tree.
    /// Perform cascading cut up the tree.
    unsafe fn cascading_cut(&mut self, y: NonNull<FibNodeInner>) {
        // Safety: caller guarantees y is valid
        unsafe {
            if let Some(parent) = (*y.as_ptr()).parent.get() {
                if !(*y.as_ptr()).marked.get() {
                    (*y.as_ptr()).marked.set(true);
                } else {
                    self.cut(y, parent);
                    self.cascading_cut(parent);
                }
            }
        }
    }
}

impl<H: Heapable> Default for FibHeap<H> {
    fn default() -> Self {
        Self::new()
    }
}

impl<H: Heapable> Drop for FibHeap<H> {
    fn drop(&mut self) {
        // Drain all nodes to properly deallocate
        while self.pop().is_some() {}
    }
}

impl<H: Heapable> From<Vec<H::Item>> for FibHeap<H> {
    fn from(items: Vec<H::Item>) -> Self {
        let mut heap = Self::new();
        for item in items {
            heap.push(item);
        }
        heap
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::heap::HeapOrdering;

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
        let heap: FibHeap<MaxHeap> = FibHeap::new();
        assert_eq!(heap.len(), 0);
        assert!(heap.is_empty());
    }

    #[test]
    fn test_pop_empty_heap_returns_none() {
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_single_element() {
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
        heap.push(42);
        assert_eq!(heap.len(), 1);
        assert_eq!(heap.peek(), Some(&42));
        assert_eq!(heap.pop(), Some(42));
        assert_eq!(heap.len(), 0);
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_max_heap_ordering() {
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
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
        let mut heap: FibHeap<MinHeap> = FibHeap::new();
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
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
        for i in 1..=10 {
            heap.push(i);
        }

        for i in (1..=10).rev() {
            assert_eq!(heap.pop(), Some(i));
        }
    }

    #[test]
    fn test_already_sorted_descending() {
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
        for i in (1..=10).rev() {
            heap.push(i);
        }

        for i in (1..=10).rev() {
            assert_eq!(heap.pop(), Some(i));
        }
    }

    #[test]
    fn test_duplicate_values() {
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
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
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
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
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
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
        let mut heap: FibHeap<MinHeap> = FibHeap::from(reversed);

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

        let mut heap: FibHeap<PriorityHeap> = FibHeap::new();
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
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
        heap.push(0);
        heap.push(-1);
        heap.push(1);

        assert_eq!(heap.pop(), Some(1));
        assert_eq!(heap.pop(), Some(0));
        assert_eq!(heap.pop(), Some(-1));
    }

    #[test]
    fn test_decrease_key() {
        // For a min-heap, decrease_key makes a value smaller (higher priority)
        let mut heap: FibHeap<MinHeap> = FibHeap::new();

        let h1 = heap.push(10);
        let _h2 = heap.push(20);
        let h3 = heap.push(30);

        assert_eq!(heap.peek(), Some(&10));

        // Decrease 30 -> 5 (now highest priority)
        unsafe {
            heap.decrease_key(h3, 5);
        }
        assert_eq!(heap.peek(), Some(&5));
        assert_eq!(heap.pop(), Some(5));

        // Decrease 10 -> 8
        unsafe {
            heap.decrease_key(h1, 8);
        }
        assert_eq!(heap.pop(), Some(8));
        assert_eq!(heap.pop(), Some(20));
    }

    #[test]
    fn test_delete() {
        let mut heap: FibHeap<MinHeap> = FibHeap::new();

        let _h1 = heap.push(10);
        let h2 = heap.push(20);
        let _h3 = heap.push(30);

        assert_eq!(heap.len(), 3);

        // Delete the middle element
        let deleted = unsafe { heap.delete(h2) };
        assert_eq!(deleted, 20);
        assert_eq!(heap.len(), 2);

        assert_eq!(heap.pop(), Some(10));
        assert_eq!(heap.pop(), Some(30));
        assert_eq!(heap.pop(), None);
    }

    #[test]
    fn test_peek() {
        let mut heap: FibHeap<MaxHeap> = FibHeap::new();
        assert_eq!(heap.peek(), None);

        heap.push(5);
        assert_eq!(heap.peek(), Some(&5));

        heap.push(10);
        assert_eq!(heap.peek(), Some(&10));

        heap.push(3);
        assert_eq!(heap.peek(), Some(&10));
    }
}
