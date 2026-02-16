//! Shared Swiss table probing logic for hash tables.
//!
//! This module provides common functionality for Swiss table-style hash tables,
//! avoiding code duplication across `Table`, `FrozenTable`, and `FrozenBigTable`.
//!
//! The implementation uses 16-byte groups and bitmask operations that LLVM can
//! optimize to SIMD instructions (SSE2/NEON) when available.

/// Number of control bytes per group (matches SIMD register width).
const GROUP: usize = 16;

/// Control byte for empty slots.
const EMPTY: u8 = 0xFF;

/// Control byte for deleted slots (high bit set, distinguishes from h2 tags).
const DELETED: u8 = 0x80;

/// A group of 16 control bytes, sized to match SIMD registers.
///
/// Operations on this type are written to enable LLVM auto-vectorization.
#[derive(Clone, Copy)]
#[repr(align(16))]
struct Group([u8; GROUP]);

impl Group {
    /// Load a group from a slice at the given offset.
    /// 
    /// If the group would extend past the end of the slice, it wraps around
    /// (for tables where ctrl has extra bytes) or returns a partial match.
    #[inline]
    fn load(ctrl: &[u8], offset: usize) -> Self {
        let mut group = [EMPTY; GROUP];
        let len = ctrl.len();
        // Unrolled copy - LLVM can vectorize this
        for i in 0..GROUP {
            let idx = (offset + i) % len;
            group[i] = ctrl[idx];
        }
        Group(group)
    }

    /// Load a group from a slice, assuming offset is group-aligned and 
    /// there are at least GROUP bytes available.
    #[inline]
    fn load_aligned(ctrl: &[u8], group_base: usize) -> Self {
        let mut group = [EMPTY; GROUP];
        // This pattern vectorizes well - copy 16 consecutive bytes
        group.copy_from_slice(&ctrl[group_base..group_base + GROUP]);
        Group(group)
    }

    /// Returns a bitmask where bit i is set if ctrl[i] == needle.
    /// 
    /// This is the core SIMD-friendly operation: compare 16 bytes against
    /// a broadcast value and pack results into a bitmask.
    #[inline]
    fn match_byte(&self, needle: u8) -> u32 {
        // Written to enable auto-vectorization:
        // - Fixed iteration count (unrollable)
        // - Simple comparison + shift pattern
        // - No branches inside loop
        let mut mask = 0u32;
        for i in 0..GROUP {
            // The comparison becomes a SIMD compare, the shift+or becomes
            // a movemask instruction on x86
            mask |= ((self.0[i] == needle) as u32) << i;
        }
        mask
    }

    /// Returns a bitmask where bit i is set if ctrl[i] == EMPTY.
    #[inline]
    fn match_empty(&self) -> u32 {
        self.match_byte(EMPTY)
    }

    /// Returns a bitmask where bit i is set if ctrl[i] is EMPTY or DELETED.
    /// (i.e., the slot is available for insertion)
    #[allow(dead_code)]
    #[inline]
    fn match_empty_or_deleted(&self) -> u32 {
        // EMPTY = 0xFF, DELETED = 0x80
        // Both have high bit set, valid h2 values have high bit clear (0x00-0x7F)
        let mut mask = 0u32;
        for i in 0..GROUP {
            mask |= ((self.0[i] & 0x80 != 0) as u32) << i;
        }
        mask
    }
}

/// Iterator over set bits in a bitmask.
#[derive(Clone, Copy)]
struct BitMaskIter {
    mask: u32,
}

impl BitMaskIter {
    #[inline]
    fn new(mask: u32) -> Self {
        BitMaskIter { mask }
    }
}

impl Iterator for BitMaskIter {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.mask == 0 {
            None
        } else {
            let bit = self.mask.trailing_zeros() as usize;
            self.mask &= self.mask - 1; // Clear lowest set bit
            Some(bit)
        }
    }
}

/// Extract h2 (7-bit tag) from a hash value.
#[inline]
pub fn h2(hash: u64) -> u8 {
    ((hash >> 57) as u8) & 0x7F
}

/// Extract h1 (bucket index) from a hash value.
#[inline]
pub fn h1(hash: u64, mask: usize) -> usize {
    (hash as usize) & mask
}

/// Locate a key in a Swiss table for read-only access.
///
/// Returns `Some(idx)` if the key is found, `None` otherwise.
/// This version is for frozen/immutable tables that don't have DELETED slots.
#[inline]
pub fn locate_readonly<K: Eq>(
    ctrl: &[u8],
    keys: &[K],
    key: &K,
    hash: u64,
    bits: usize,
) -> Option<usize> {
    if ctrl.is_empty() || bits == 0 {
        return None;
    }

    let capacity = 1usize << bits;
    let mask = capacity - 1;
    let h2_val = h2(hash);
    let mut probe = 0usize;
    let mut pos = h1(hash, mask);

    loop {
        let group_base = pos & !(GROUP - 1);
        
        // Load 16 control bytes and find matches via SIMD-friendly bitmask ops
        let group = if group_base + GROUP <= ctrl.len() {
            Group::load_aligned(ctrl, group_base)
        } else {
            Group::load(ctrl, group_base)
        };

        // Check for matching h2 tags
        let h2_matches = group.match_byte(h2_val);
        for bit in BitMaskIter::new(h2_matches) {
            let idx = (group_base + bit) & mask;
            if keys[idx] == *key {
                return Some(idx);
            }
        }

        // Check for empty slots - if we hit one, key doesn't exist
        let empty_matches = group.match_empty();
        if empty_matches != 0 {
            return None;
        }

        probe += 1;
        pos = (pos + probe * GROUP) & mask;
        if probe > mask / GROUP + 1 {
            return None;
        }
    }
}

/// Locate a key or insertion point in a Swiss table for read-write access.
///
/// Returns `Some(idx)` where:
/// - If the key exists at `idx`, `ctrl[idx]` will be h2 and `keys[idx] == key`
/// - If the key doesn't exist, `idx` is where it should be inserted
///   (either an EMPTY slot or the first DELETED slot encountered)
///
/// Returns `None` if the table is full (shouldn't happen with proper load factor).
#[inline]
pub fn locate_readwrite<K: Eq>(
    ctrl: &[u8],
    keys: &[K],
    key: &K,
    hash: u64,
    bits: usize,
) -> Option<usize> {
    if ctrl.is_empty() {
        return None;
    }

    let capacity = 1usize << bits;
    let mask = capacity - 1;
    let h2_val = h2(hash);
    let mut probe = 0usize;
    let mut pos = h1(hash, mask);
    let mut first_deleted: Option<usize> = None;

    loop {
        let group_base = pos & !(GROUP - 1);
        
        let group = if group_base + GROUP <= ctrl.len() {
            Group::load_aligned(ctrl, group_base)
        } else {
            Group::load(ctrl, group_base)
        };

        // Check for matching h2 tags - if key exists, return its slot
        let h2_matches = group.match_byte(h2_val);
        for bit in BitMaskIter::new(h2_matches) {
            let idx = (group_base + bit) & mask;
            if keys[idx] == *key {
                return Some(idx);
            }
        }

        // Track first deleted slot for potential insertion
        if first_deleted.is_none() {
            let deleted_matches = group.match_byte(DELETED);
            if deleted_matches != 0 {
                let bit = deleted_matches.trailing_zeros() as usize;
                first_deleted = Some((group_base + bit) & mask);
            }
        }

        // Check for empty slots - key doesn't exist, return insertion point
        let empty_matches = group.match_empty();
        if empty_matches != 0 {
            if let Some(deleted_idx) = first_deleted {
                return Some(deleted_idx);
            }
            let bit = empty_matches.trailing_zeros() as usize;
            return Some((group_base + bit) & mask);
        }

        probe += 1;
        pos = (pos + probe * GROUP) & mask;
        if probe > mask / GROUP + 1 {
            return first_deleted;
        }
    }
}

/// Find an empty slot for insertion (used when building frozen tables).
///
/// This assumes the table has capacity and will always find an EMPTY slot.
#[inline]
pub fn find_empty_slot(ctrl: &[u8], hash: u64, bits: usize) -> usize {
    let capacity = 1usize << bits;
    let mask = capacity - 1;
    let mut pos = h1(hash, mask);

    loop {
        let group_base = pos & !(GROUP - 1);
        
        let group = if group_base + GROUP <= ctrl.len() {
            Group::load_aligned(ctrl, group_base)
        } else {
            Group::load(ctrl, group_base)
        };

        let empty_matches = group.match_empty();
        if empty_matches != 0 {
            let bit = empty_matches.trailing_zeros() as usize;
            return (group_base + bit) & mask;
        }

        pos = (pos + GROUP) & mask;
    }
}

/// Check if a control byte indicates an occupied slot.
#[inline]
pub fn is_occupied(ctrl: u8) -> bool {
    // Occupied slots have h2 values 0x00-0x7F (high bit clear)
    // EMPTY = 0xFF and DELETED = 0x80 both have high bit set
    ctrl & 0x80 == 0
}

/// Compute the initial probe position for a key, returning (group_base, mask).
/// Used by prefetch methods to know which cache line to warm up.
#[inline]
pub fn probe_position(hash: u64, bits: usize) -> (usize, usize) {
    let capacity = 1usize << bits;
    let mask = capacity - 1;
    let pos = h1(hash, mask);
    let group_base = pos & !(GROUP - 1);
    (group_base, mask)
}

/// Constants for use by implementations.
pub const CTRL_EMPTY: u8 = EMPTY;
pub const CTRL_DELETED: u8 = DELETED;
pub const PROBE_GROUP: usize = GROUP;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_match_byte() {
        let ctrl = [0x00, 0x01, 0x42, 0x42, 0xFF, 0x80, 0x42, 0x00,
                    0xFF, 0xFF, 0x42, 0x01, 0x00, 0x80, 0xFF, 0x42];
        let group = Group(ctrl);
        
        // Match 0x42 - should be at positions 2, 3, 6, 10, 15
        let mask = group.match_byte(0x42);
        assert_eq!(mask, 0b1000_0100_0100_1100);
        
        // Iterate over matches
        let positions: Vec<usize> = BitMaskIter::new(mask).collect();
        assert_eq!(positions, vec![2, 3, 6, 10, 15]);
    }

    #[test]
    fn test_group_match_empty() {
        let ctrl = [0x00, 0x01, 0xFF, 0x42, 0xFF, 0x80, 0x42, 0x00,
                    0xFF, 0xFF, 0x42, 0x01, 0x00, 0x80, 0xFF, 0x42];
        let group = Group(ctrl);
        
        let mask = group.match_empty();
        // EMPTY (0xFF) at positions 2, 4, 8, 9, 14
        assert_eq!(mask, 0b0100_0011_0001_0100);
    }

    #[test]
    fn test_group_match_empty_or_deleted() {
        let ctrl = [0x00, 0x01, 0xFF, 0x42, 0xFF, 0x80, 0x42, 0x00,
                    0xFF, 0xFF, 0x42, 0x01, 0x00, 0x80, 0xFF, 0x42];
        let group = Group(ctrl);
        
        let mask = group.match_empty_or_deleted();
        // EMPTY (0xFF) at 2, 4, 8, 9, 14 and DELETED (0x80) at 5, 13
        assert_eq!(mask, 0b0110_0011_0011_0100);
    }

    #[test]
    fn test_bitmask_iter() {
        let mask = 0b1010_0101_0000_1100u32;
        let bits: Vec<usize> = BitMaskIter::new(mask).collect();
        assert_eq!(bits, vec![2, 3, 8, 10, 13, 15]);
    }

    #[test]
    fn test_is_occupied() {
        assert!(is_occupied(0x00));  // Valid h2
        assert!(is_occupied(0x42));  // Valid h2
        assert!(is_occupied(0x7F));  // Valid h2 (max)
        assert!(!is_occupied(0x80)); // DELETED
        assert!(!is_occupied(0xFF)); // EMPTY
    }
}
