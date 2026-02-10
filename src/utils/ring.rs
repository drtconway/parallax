use std::ops::Index;

pub struct Ring<T: Default + Copy, const N: usize> {
    buffer: [T; N],
    head: usize,
    tail: usize,
}

impl<T: Default + Copy, const N: usize> Ring<T, N> {
    pub fn new() -> Self {
        Self {
            buffer: [T::default(); N],
            head: 0,
            tail: 0,
        }
    }

    pub fn push(&mut self, value: T) {
        self.buffer[self.head] = value;
        self.head = (self.head + 1) % N;
        if self.head == self.tail {
            self.tail = (self.tail + 1) % N;
        }
    }

    pub fn pop(&mut self) -> Option<T> {
        if self.head == self.tail {
            None
        } else {
            let value = self.buffer[self.tail];
            self.tail = (self.tail + 1) % N;
            Some(value)
        }
    }

    pub fn len(&self) -> usize {
        (self.head + N - self.tail) % N
    }

    #[allow(dead_code)]
    pub fn is_full(&self) -> bool {
        (self.head + 1) % N == self.tail
    }

    /// Returns the oldest element (front of the queue).
    #[allow(dead_code)]
    pub fn front(&self) -> T {
        self.buffer[self.tail]
    }

    /// Returns the newest element (back of the queue).
    pub fn back(&self) -> T {
        let pos = (self.head + N - 1) % N;
        self.buffer[pos]
    }
}

impl<T: Default + Copy, const N: usize> Index<usize> for Ring<T, N> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.len() {
            panic!("index out of bounds");
        }
        let pos = (self.tail + index) % N;
        &self.buffer[pos]
    }
}

#[cfg(test)]
mod tests {
    use super::Ring;

    #[test]
    fn new_ring_is_empty() {
        let ring: Ring<i32, 4> = Ring::new();
        assert_eq!(ring.len(), 0);
        assert!(!ring.is_full());
    }

    #[test]
    fn push_increases_len() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        assert_eq!(ring.len(), 1);
        ring.push(2);
        assert_eq!(ring.len(), 2);
        ring.push(3);
        assert_eq!(ring.len(), 3);
    }

    #[test]
    fn pop_decreases_len() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        ring.push(2);
        assert_eq!(ring.len(), 2);
        ring.pop();
        assert_eq!(ring.len(), 1);
        ring.pop();
        assert_eq!(ring.len(), 0);
    }

    #[test]
    fn pop_returns_fifo_order() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(10);
        ring.push(20);
        ring.push(30);
        assert_eq!(ring.pop(), Some(10));
        assert_eq!(ring.pop(), Some(20));
        assert_eq!(ring.pop(), Some(30));
    }

    #[test]
    fn pop_empty_returns_none() {
        let mut ring: Ring<i32, 4> = Ring::new();
        assert_eq!(ring.pop(), None);
    }

    #[test]
    fn is_full_when_capacity_reached() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        ring.push(2);
        ring.push(3);
        assert!(ring.is_full());
    }

    #[test]
    fn push_overwrites_oldest_when_full() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        ring.push(2);
        ring.push(3);
        assert!(ring.is_full());
        
        // Push another element, should overwrite oldest (1)
        ring.push(4);
        assert!(ring.is_full());
        
        // First pop should return 2, not 1
        assert_eq!(ring.pop(), Some(2));
        assert_eq!(ring.pop(), Some(3));
        assert_eq!(ring.pop(), Some(4));
        assert_eq!(ring.pop(), None);
    }

    #[test]
    fn index_access() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(10);
        ring.push(20);
        ring.push(30);
        
        assert_eq!(ring[0], 10);
        assert_eq!(ring[1], 20);
        assert_eq!(ring[2], 30);
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn index_out_of_bounds_panics() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        let _ = ring[1]; // Only one element, index 1 is out of bounds
    }

    #[test]
    fn index_after_pop() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(10);
        ring.push(20);
        ring.push(30);
        ring.pop(); // Remove 10
        
        assert_eq!(ring[0], 20);
        assert_eq!(ring[1], 30);
    }

    #[test]
    fn wraparound_behavior() {
        let mut ring: Ring<i32, 4> = Ring::new();
        
        // Fill and empty multiple times to test wraparound
        for round in 0..3 {
            let base = round * 10;
            ring.push(base + 1);
            ring.push(base + 2);
            ring.push(base + 3);
            
            assert_eq!(ring.pop(), Some(base + 1));
            assert_eq!(ring.pop(), Some(base + 2));
            assert_eq!(ring.pop(), Some(base + 3));
            assert_eq!(ring.pop(), None);
        }
    }

    #[test]
    fn interleaved_push_pop() {
        let mut ring: Ring<i32, 4> = Ring::new();
        
        ring.push(1);
        ring.push(2);
        assert_eq!(ring.pop(), Some(1));
        
        ring.push(3);
        ring.push(4);
        assert_eq!(ring.pop(), Some(2));
        assert_eq!(ring.pop(), Some(3));
        
        ring.push(5);
        assert_eq!(ring.pop(), Some(4));
        assert_eq!(ring.pop(), Some(5));
        assert_eq!(ring.pop(), None);
    }

    #[test]
    fn works_with_different_types() {
        // Test with f64
        let mut float_ring: Ring<f64, 3> = Ring::new();
        float_ring.push(1.5);
        float_ring.push(2.5);
        assert_eq!(float_ring.pop(), Some(1.5));
        assert_eq!(float_ring.pop(), Some(2.5));

        // Test with u8
        let mut byte_ring: Ring<u8, 3> = Ring::new();
        byte_ring.push(255);
        byte_ring.push(0);
        assert_eq!(byte_ring.pop(), Some(255));
        assert_eq!(byte_ring.pop(), Some(0));
    }

    #[test]
    fn len_after_overwrite() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        ring.push(2);
        ring.push(3);
        assert_eq!(ring.len(), 3);
        assert!(ring.is_full());
        
        // Overwrite oldest
        ring.push(4);
        assert_eq!(ring.len(), 3);
        assert!(ring.is_full());
        
        ring.push(5);
        assert_eq!(ring.len(), 3);
        assert!(ring.is_full());
    }

    #[test]
    fn front_returns_oldest() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(10);
        assert_eq!(ring.front(), 10);
        ring.push(20);
        assert_eq!(ring.front(), 10);
        ring.push(30);
        assert_eq!(ring.front(), 10);
        
        // After overwrite, front should be second oldest
        ring.push(40);
        assert_eq!(ring.front(), 20);
    }

    #[test]
    fn back_returns_newest() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(10);
        assert_eq!(ring.back(), 10);
        ring.push(20);
        assert_eq!(ring.back(), 20);
        ring.push(30);
        assert_eq!(ring.back(), 30);
        ring.push(40);
        assert_eq!(ring.back(), 40);
    }

    #[test]
    fn front_and_back_after_pop() {
        let mut ring: Ring<i32, 4> = Ring::new();
        ring.push(1);
        ring.push(2);
        ring.push(3);
        
        assert_eq!(ring.front(), 1);
        assert_eq!(ring.back(), 3);
        
        ring.pop();
        assert_eq!(ring.front(), 2);
        assert_eq!(ring.back(), 3);
        
        ring.pop();
        assert_eq!(ring.front(), 3);
        assert_eq!(ring.back(), 3);
    }
}