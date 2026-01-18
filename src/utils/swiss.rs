//! Shared Swiss table probing logic for hash tables.
//!
//! This module provides common functionality for Swiss table-style hash tables,
//! avoiding code duplication across `Table`, `FrozenTable`, and `FrozenBigTable`.

const GROUP: usize = 16;
const EMPTY: u8 = 0xFF;
const DELETED: u8 = 0x80;

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

    let mask = (1usize << bits) - 1;
    let h2_val = h2(hash);
    let mut probe = 0usize;
    let mut slot = h1(hash, mask);

    loop {
        let group_base = slot & !(GROUP - 1);
        for offset in 0..GROUP {
            let idx = (group_base + offset) & mask;
            let c = ctrl[idx];
            if c == EMPTY {
                return None;
            }
            if c == h2_val && keys[idx] == *key {
                return Some(idx);
            }
        }

        probe += 1;
        slot = (slot + probe * GROUP) & mask;
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

    let mask = (1usize << bits) - 1;
    let h2_val = h2(hash);
    let mut probe = 0usize;
    let mut slot = h1(hash, mask);
    let mut first_deleted: Option<usize> = None;

    loop {
        let group_base = slot & !(GROUP - 1);
        for offset in 0..GROUP {
            let idx = (group_base + offset) & mask;
            let c = ctrl[idx];
            if c == EMPTY {
                return Some(first_deleted.unwrap_or(idx));
            }
            if c == DELETED {
                if first_deleted.is_none() {
                    first_deleted = Some(idx);
                }
                continue;
            }
            if c == h2_val && keys[idx] == *key {
                return Some(idx);
            }
        }

        probe += 1;
        slot = (slot + probe * GROUP) & mask;
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
        let group_start = pos & !(GROUP - 1);
        for i in 0..GROUP {
            let idx = (group_start + i) & mask;
            if ctrl[idx] == EMPTY {
                return idx;
            }
        }
        pos = (pos + GROUP) & mask;
    }
}

/// Check if a control byte indicates an occupied slot.
#[inline]
pub fn is_occupied(ctrl: u8) -> bool {
    ctrl != EMPTY && ctrl != DELETED
}

/// Constants for use by implementations.
pub const CTRL_EMPTY: u8 = EMPTY;
pub const CTRL_DELETED: u8 = DELETED;
pub const PROBE_GROUP: usize = GROUP;
