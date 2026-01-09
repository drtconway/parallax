use bitvec::vec::BitVec;

pub struct Table<K: Default + Clone + Eq + std::hash::Hash, V: Default + Clone> {
    count: usize,
    seed: u64,
    bits: usize,
    used: BitVec,
    deleted: BitVec,
    keys: Vec<K>,
    values: Vec<V>,
}

impl<K: Default + Clone + Eq + std::hash::Hash, V: Default + Clone> Table<K, V> {
    /// Creates a new, empty table.
    pub fn new() -> Self {
        Table {
            count: 0,
            seed: 0xcbf2_9ce4_8422_2325,
            bits: 0,
            used: BitVec::new(),
            deleted: BitVec::new(),
            keys: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Creates a new table with the specified number of bits.
    pub fn new_with_width(bits: usize) -> Self {
        let n = 1usize << bits;
        Table {
            count: 0,
            seed: 0xcbf2_9ce4_8422_2325,
            bits,
            used: BitVec::repeat(false, n),
            deleted: BitVec::repeat(false, n),
            keys: vec![K::default(); n],
            values: vec![V::default(); n],
        }
    }

    /// Returns the number of entries in the table.
    pub fn len(&self) -> usize {
        self.count
    }

    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.count = 0;
        self.used.fill(false);
        self.deleted.fill(false);
    }

    pub fn swap(&mut self, other: &mut Self) {
        std::mem::swap(self, other);
    }

    #[allow(dead_code)]
    pub fn contains_key(&self, key: &K) -> bool {
        let slot = self.locate(key);
        match slot {
            Some(idx) => self.used[idx] && self.keys[idx] == *key,
            None => false,
        }
    }

    pub fn get(&self, key: &K) -> Option<&V> {
        let slot = self.locate(key);
        match slot {
            Some(idx) if self.used[idx] && self.keys[idx] == *key => Some(&self.values[idx]),
            _ => None,
        }
    }

    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if self.used.is_empty() {
            assert_eq!(self.count, 0);
            assert_eq!(self.bits, 0);
            self.bits = 4;
            let n = 1usize << self.bits;
            self.used = BitVec::repeat(false, n);
            self.deleted = BitVec::repeat(false, n);
            self.keys = vec![K::default(); n];
            self.values = vec![V::default(); n];
        }

        // Grow when table is at or above 75% occupancy.
        if self.count * 4 >= self.used.len() * 3 {
            self.rehash();
        }

        let mut slot = self.locate(&key);
        if slot.is_none() {
            self.rehash();
            slot = self.locate(&key);
        }

        let idx = slot.expect("table must have capacity after rehash");
        if self.used[idx] {
            let old_value = std::mem::replace(&mut self.values[idx], value);
            Some(old_value)
        } else {
            self.used.set(idx, true);
            self.deleted.set(idx, false);
            self.keys[idx] = key;
            self.values[idx] = value;
            self.count += 1;
            None
        }
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        let slot = self.locate(key)?;
        if self.used[slot] && self.keys[slot] == *key {
            self.used.set(slot, false);
            self.deleted.set(slot, true);
            self.count -= 1;
            let mut value = V::default();
            std::mem::swap(&mut value, &mut self.values[slot]);
            Some(value)
        } else {
            None
        }
    }

    fn locate(&self, key: &K) -> Option<usize> {
        const GROUP: usize = 16; // match SwissTable group size

        if self.used.is_empty() {
            return None;
        }

        let mask = (1 << self.bits) - 1; // table width is power-of-two
        let mut hash = self.hash(key) as usize;
        let mut probe = 0usize;

        loop {
            // Round down to the start of the current group.
            let group_base = (hash & mask) & !(GROUP - 1);

            // Scan the group for either the key or the first empty slot.
            let mut first_deleted: Option<usize> = None;
            for offset in 0..GROUP {
                let idx = (group_base + offset) & mask;
                if self.used[idx] {
                    if self.keys[idx] == *key {
                        return Some(idx);
                    }
                } else if self.deleted[idx] {
                    if first_deleted.is_none() {
                        first_deleted = Some(idx);
                    }
                } else {
                    // truly empty slot ends the probe
                    return Some(first_deleted.unwrap_or(idx));
                }
            }

            // Advance to next group using SwissTable-style quadratic-ish step.
            probe += 1;
            hash = (hash + probe * GROUP) & mask;

            // Give up if we've wrapped all groups (table is full or corrupted).
            if probe > mask / GROUP + 1 {
                return None;
            }
        }
    }

    fn rehash(&mut self) {
        let mut tbl = Table::new_with_width(self.bits + 1);
        for i in 0..self.used.len() {
            if self.used[i] {
                let mut key = K::default();
                std::mem::swap(&mut key, &mut self.keys[i]);
                let mut value = V::default();
                std::mem::swap(&mut value, &mut self.values[i]);
                tbl.insert(key, value);
            }
        }
        self.swap(&mut tbl);
    }

    /// Hashes a key to produce a u64 hash value.
    fn hash(&self, key: &K) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.seed.hash(&mut hasher);
        key.hash(&mut hasher);
        hasher.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::Table;

    #[test]
    fn insert_remove_round_trip() {
        let mut tbl: Table<u64, u32> = Table::new();
        assert_eq!(tbl.len(), 0);
        assert!(!tbl.contains_key(&1));

        assert_eq!(tbl.insert(1, 10), None);
        assert_eq!(tbl.len(), 1);
        assert!(tbl.contains_key(&1));

        assert_eq!(tbl.remove(&1), Some(10));
        assert_eq!(tbl.len(), 0);
        assert!(!tbl.contains_key(&1));
    }

    #[test]
    fn insert_updates_existing() {
        let mut tbl: Table<u64, u32> = Table::new();
        assert_eq!(tbl.insert(2, 5), None);
        assert_eq!(tbl.insert(2, 7), Some(5));
        assert_eq!(tbl.len(), 1);
        assert_eq!(tbl.remove(&2), Some(7));
    }

    #[test]
    fn rehash_preserves_entries() {
        let mut tbl: Table<u64, u32> = Table::new();
        for i in 0..20u64 {
            assert!(tbl.insert(i, i as u32).is_none());
        }
        assert_eq!(tbl.len(), 20);
        for i in 0..20u64 {
            assert!(tbl.contains_key(&i));
            assert_eq!(tbl.remove(&i), Some(i as u32));
        }
        assert_eq!(tbl.len(), 0);
    }

    #[test]
    fn clear_resets_usage() {
        let mut tbl: Table<u64, u32> = Table::new();
        for i in 0..4u64 {
            tbl.insert(i, i as u32);
        }
        assert_eq!(tbl.len(), 4);
        tbl.clear();
        assert_eq!(tbl.len(), 0);
        for i in 0..4u64 {
            assert!(!tbl.contains_key(&i));
        }
    }
}
