#![allow(dead_code)]
/// A vector with deletions: O(log(n/64)) indexed access among live elements.
///
/// Maintains a flat vector of elements alongside a bitmap of deleted positions
/// and a prefix-sum array of deleted counts per 64-element block, inspired by
/// rank/select on succinct bitvectors.
///
/// - **Logical index**: position among live (non-deleted) elements (0-based).
/// - **Physical index**: position in the underlying `Vec<T>`.
///
/// Element access by logical index is O(log(n/64)) — binary search over
/// blocks plus a popcount scan within a block. Deletion is O(n/64) due to
/// the prefix-sum update.
///
/// Deleted elements remain in memory until the `Pothole` is dropped.
pub struct Pothole<T> {
    items: Vec<T>,
    /// Bitmap: bit `i` of `deleted[i / 64]` is 1 iff physical position `i`
    /// has been deleted. Bits beyond `items.len()` are always 0.
    deleted: Vec<u64>,
    /// `prefix_deleted[b]` = number of deleted elements in physical positions
    /// `[0, b * 64)`.  Length = `deleted.len() + 1`.
    prefix_deleted: Vec<usize>,
    /// Number of live elements.
    alive: usize,
}

impl<T> Pothole<T> {
    /// Create an empty `Pothole`.
    pub fn new() -> Self {
        Pothole {
            items: Vec::new(),
            deleted: Vec::new(),
            prefix_deleted: vec![0],
            alive: 0,
        }
    }

    /// Create a `Pothole` from an existing `Vec`, with no deletions.
    pub fn from_vec(items: Vec<T>) -> Self {
        let capacity = items.len();
        let num_blocks = (capacity + 63) / 64;
        let deleted = vec![0u64; num_blocks];
        let prefix_deleted = vec![0usize; num_blocks + 1];
        Pothole {
            items,
            deleted,
            prefix_deleted,
            alive: capacity,
        }
    }

    /// Number of live (non-deleted) elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.alive
    }

    /// Whether there are no live elements.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.alive == 0
    }

    /// Number of live elements in physical positions `[0, block * 64)`.
    #[inline]
    fn alive_before_block(&self, block: usize) -> usize {
        let phys = (block * 64).min(self.items.len());
        phys - self.prefix_deleted[block]
    }

    /// Mask of valid bit positions within a block (1 = valid position).
    #[inline]
    fn valid_mask(&self, block: usize) -> u64 {
        let w = {
            let start = block * 64;
            let end = ((block + 1) * 64).min(self.items.len());
            end - start
        };
        if w >= 64 { !0u64 } else { (1u64 << w) - 1 }
    }

    /// Convert a logical index (among live elements) to a physical index.
    ///
    /// Returns `None` if `logical >= self.len()`.
    fn to_physical(&self, logical: usize) -> Option<usize> {
        if logical >= self.alive {
            return None;
        }

        let num_blocks = self.deleted.len();

        // Binary search: find the rightmost block b where
        // alive_before_block(b) <= logical.
        let mut lo = 0usize;
        let mut hi = num_blocks;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.alive_before_block(mid + 1) <= logical {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let block = lo;

        // How many alive elements do we need within this block?
        let remaining = logical - self.alive_before_block(block);

        // Find the `remaining`-th alive (zero) bit in deleted[block],
        // masked to valid positions.
        let alive_bits = !self.deleted[block] & self.valid_mask(block);
        let bit = select_bit(alive_bits, remaining);

        Some(block * 64 + bit)
    }

    /// Access a live element by logical index.
    pub fn get(&self, logical: usize) -> Option<&T> {
        self.to_physical(logical).map(|p| &self.items[p])
    }

    /// Mutably access a live element by logical index.
    pub fn get_mut(&mut self, logical: usize) -> Option<&mut T> {
        self.to_physical(logical).map(|p| &mut self.items[p])
    }

    /// Delete the element at the given logical index.
    ///
    /// Returns `true` if the element was deleted, `false` if the index was
    /// out of range.
    pub fn delete(&mut self, logical: usize) -> bool {
        let Some(physical) = self.to_physical(logical) else {
            return false;
        };
        let block = physical / 64;
        let bit = physical % 64;
        debug_assert!(self.deleted[block] & (1u64 << bit) == 0, "double delete");
        self.deleted[block] |= 1u64 << bit;
        self.alive -= 1;
        // Update prefix sums for all subsequent blocks.
        for i in (block + 1)..self.prefix_deleted.len() {
            self.prefix_deleted[i] += 1;
        }
        true
    }

    /// Delete multiple elements given their logical indices.
    ///
    /// `indices` **must** be sorted in ascending order. All indices are
    /// resolved against the current live set before any modifications,
    /// then the bitmap is updated and prefix sums rebuilt once.
    pub fn delete_many(&mut self, indices: &[usize]) {
        if indices.is_empty() {
            return;
        }
        debug_assert!(
            indices.windows(2).all(|w| w[0] < w[1]),
            "delete_many: indices must be strictly ascending"
        );

        // Resolve all logical indices to physical before modifying anything.
        let physical_indices: Vec<usize> = indices
            .iter()
            .filter_map(|&logical| self.to_physical(logical))
            .collect();

        // Apply all deletions to the bitmap.
        for &phys in &physical_indices {
            debug_assert!(
                self.deleted[phys / 64] & (1u64 << (phys % 64)) == 0,
                "delete_many: element already deleted"
            );
            self.deleted[phys / 64] |= 1u64 << (phys % 64);
        }

        // Rebuild prefix sums once.
        self.prefix_deleted[0] = 0;
        for b in 0..self.deleted.len() {
            self.prefix_deleted[b + 1] =
                self.prefix_deleted[b] + self.deleted[b].count_ones() as usize;
        }
        self.alive = self.items.len() - *self.prefix_deleted.last().unwrap();
    }
}

// ── Helpers ─────────────────────────────────────────────────────────────────

/// Find the position of the `n`-th set bit (0-indexed) in `word`.
///
/// Panics if `word` has fewer than `n + 1` set bits.
#[inline]
fn select_bit(mut word: u64, n: usize) -> usize {
    debug_assert!(
        word.count_ones() as usize > n,
        "select_bit: not enough set bits ({} ones, asked for index {})",
        word.count_ones(),
        n,
    );
    for _ in 0..n {
        word &= word - 1; // clear lowest set bit
    }
    word.trailing_zeros() as usize
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_vec_basic() {
        let p = Pothole::from_vec(vec![10, 20, 30, 40, 50]);
        assert_eq!(p.len(), 5);
        assert_eq!(*p.get(0).unwrap(), 10);
        assert_eq!(*p.get(1).unwrap(), 20);
        assert_eq!(*p.get(4).unwrap(), 50);
    }

    #[test]
    fn test_delete_first() {
        let mut p = Pothole::from_vec(vec![10, 20, 30, 40, 50]);
        assert!(p.delete(0)); // delete 10
        assert_eq!(p.len(), 4);
        assert_eq!(*p.get(0).unwrap(), 20);
        assert_eq!(*p.get(1).unwrap(), 30);
        assert_eq!(*p.get(3).unwrap(), 50);
    }

    #[test]
    fn test_delete_last() {
        let mut p = Pothole::from_vec(vec![10, 20, 30, 40, 50]);
        assert!(p.delete(4)); // delete 50
        assert_eq!(p.len(), 4);
        assert_eq!(*p.get(0).unwrap(), 10);
        assert_eq!(*p.get(3).unwrap(), 40);
        assert!(p.get(4).is_none());
    }

    #[test]
    fn test_delete_middle() {
        let mut p = Pothole::from_vec(vec![10, 20, 30, 40, 50]);
        assert!(p.delete(2)); // delete 30
        assert_eq!(p.len(), 4);
        assert_eq!(*p.get(0).unwrap(), 10);
        assert_eq!(*p.get(1).unwrap(), 20);
        assert_eq!(*p.get(2).unwrap(), 40);
        assert_eq!(*p.get(3).unwrap(), 50);
    }

    #[test]
    fn test_delete_multiple() {
        let mut p = Pothole::from_vec(vec![10, 20, 30, 40, 50]);
        p.delete(1); // delete 20  → [10, 30, 40, 50]
        p.delete(2); // delete 40  → [10, 30, 50]
        assert_eq!(p.len(), 3);
        assert_eq!(*p.get(0).unwrap(), 10);
        assert_eq!(*p.get(1).unwrap(), 30);
        assert_eq!(*p.get(2).unwrap(), 50);
    }

    #[test]
    fn test_delete_all() {
        let mut p = Pothole::from_vec(vec![10, 20, 30]);
        p.delete(0);
        p.delete(0);
        p.delete(0);
        assert!(p.is_empty());
        assert!(p.get(0).is_none());
    }

    #[test]
    fn test_get_mut() {
        let mut p = Pothole::from_vec(vec![10, 20, 30]);
        *p.get_mut(1).unwrap() = 99;
        assert_eq!(*p.get(1).unwrap(), 99);
        p.delete(0);
        assert_eq!(*p.get(0).unwrap(), 99);
    }

    #[test]
    fn test_large_crosses_block_boundary() {
        // 100 elements → 2 blocks (64 + 36)
        let v: Vec<usize> = (0..100).collect();
        let mut p = Pothole::from_vec(v);

        // Delete every other element from the back (keeps logical indices stable).
        // Delete logical indices 99, 97, 95, ..., 1 → removes odd values.
        for i in (1..100).rev().step_by(2) {
            p.delete(i);
        }
        assert_eq!(p.len(), 50);

        // Remaining: 0, 2, 4, ..., 98
        for i in 0..50 {
            assert_eq!(*p.get(i).unwrap(), 2 * i);
        }
    }

    #[test]
    fn test_empty() {
        let p: Pothole<i32> = Pothole::new();
        assert!(p.is_empty());
        assert_eq!(p.len(), 0);
        assert!(p.get(0).is_none());
    }

    #[test]
    fn test_single_element() {
        let mut p = Pothole::from_vec(vec![42]);
        assert_eq!(*p.get(0).unwrap(), 42);
        p.delete(0);
        assert!(p.is_empty());
    }

    #[test]
    fn test_exact_block_boundary() {
        // Exactly 64 elements = 1 full block
        let v: Vec<usize> = (0..64).collect();
        let mut p = Pothole::from_vec(v);
        p.delete(0);
        p.delete(p.len() - 1); // delete last
        assert_eq!(p.len(), 62);
        assert_eq!(*p.get(0).unwrap(), 1);
        assert_eq!(*p.get(61).unwrap(), 62);
    }

    #[test]
    fn test_out_of_range() {
        let mut p = Pothole::from_vec(vec![10, 20, 30]);
        assert!(p.get(3).is_none());
        assert!(!p.delete(3));
        assert!(p.get_mut(3).is_none());
    }

    #[test]
    fn test_select_bit() {
        assert_eq!(select_bit(0b10110, 0), 1);
        assert_eq!(select_bit(0b10110, 1), 2);
        assert_eq!(select_bit(0b10110, 2), 4);
    }

    #[test]
    fn test_delete_many_basic() {
        let mut p = Pothole::from_vec(vec![10, 20, 30, 40, 50]);
        // Delete logical 1 (20) and logical 3 (40) in one call.
        p.delete_many(&[1, 3]);
        assert_eq!(p.len(), 3);
        assert_eq!(*p.get(0).unwrap(), 10);
        assert_eq!(*p.get(1).unwrap(), 30);
        assert_eq!(*p.get(2).unwrap(), 50);
    }

    #[test]
    fn test_delete_many_all() {
        let mut p = Pothole::from_vec(vec![10, 20, 30]);
        p.delete_many(&[0, 1, 2]);
        assert!(p.is_empty());
    }

    #[test]
    fn test_delete_many_empty() {
        let mut p = Pothole::from_vec(vec![10, 20, 30]);
        p.delete_many(&[]);
        assert_eq!(p.len(), 3);
    }

    #[test]
    fn test_delete_many_cross_block() {
        // 100 elements → 2 blocks; delete all odd-indexed ones
        let v: Vec<usize> = (0..100).collect();
        let mut p = Pothole::from_vec(v);
        let odds: Vec<usize> = (1..100).step_by(2).collect(); // 1,3,5,...,99
        p.delete_many(&odds);
        assert_eq!(p.len(), 50);
        for i in 0..50 {
            assert_eq!(*p.get(i).unwrap(), 2 * i);
        }
    }

    #[test]
    fn test_delete_many_after_single_delete() {
        let mut p = Pothole::from_vec(vec![10, 20, 30, 40, 50, 60]);
        p.delete(0); // remove 10 → [20, 30, 40, 50, 60]
        p.delete_many(&[1, 3]); // remove 30 and 50 → [20, 40, 60]
        assert_eq!(p.len(), 3);
        assert_eq!(*p.get(0).unwrap(), 20);
        assert_eq!(*p.get(1).unwrap(), 40);
        assert_eq!(*p.get(2).unwrap(), 60);
    }
}
