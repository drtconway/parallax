pub struct Table<K: Default + Clone + Eq + std::hash::Hash, V: Default + Clone> {
    count: usize,
    seed: u64,
    bits: usize,
    ctrl: Vec<u8>, // control bytes (SwissTable style)
    keys: Vec<K>,
    values: Vec<V>,
}

impl<K: Default + Clone + Eq + std::hash::Hash, V: Default + Clone> Table<K, V> {
    pub fn new() -> Self {
        Table {
            count: 0,
            seed: 0xcbf2_9ce4_8422_2325,
            bits: 0,
            ctrl: Vec::new(),
            keys: Vec::new(),
            values: Vec::new(),
        }
    }

    pub fn new_with_width(bits: usize) -> Self {
        let n = 1usize << bits;
        Table {
            count: 0,
            seed: 0xcbf2_9ce4_8422_2325,
            bits,
            ctrl: vec![Self::EMPTY; n + Self::GROUP],
            keys: vec![K::default(); n],
            values: vec![V::default(); n],
        }
    }

    pub fn len(&self) -> usize {
        self.count
    }

    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.count = 0;
        self.ctrl.fill(Self::EMPTY);
    }

    pub fn swap(&mut self, other: &mut Self) {
        std::mem::swap(self, other);
    }

    #[allow(dead_code)]
    pub fn contains_key(&self, key: &K) -> bool {
        match self.locate(key) {
            Some(idx) => self.is_occupied(idx) && self.keys[idx] == *key,
            None => false,
        }
    }

    pub fn get(&self, key: &K) -> Option<&V> {
        match self.locate(key) {
            Some(idx) if self.is_occupied(idx) && self.keys[idx] == *key => {
                Some(&self.values[idx])
            }
            _ => None,
        }
    }

    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if self.ctrl.is_empty() {
            assert_eq!(self.count, 0);
            assert_eq!(self.bits, 0);
            self.bits = 4;
            let n = 1usize << self.bits;
            self.ctrl = vec![Self::EMPTY; n + Self::GROUP];
            self.keys = vec![K::default(); n];
            self.values = vec![V::default(); n];
        }

        if self.count * 4 >= self.bucket_len() * 3 {
            self.rehash();
        }

        let mut slot = self.locate(&key);
        if slot.is_none() {
            self.rehash();
            slot = self.locate(&key);
        }

        let idx = slot.expect("table must have capacity after rehash");
        if self.is_occupied(idx) {
            Some(std::mem::replace(&mut self.values[idx], value))
        } else {
            self.ctrl[idx] = self.h2(self.hash(&key));
            self.keys[idx] = key;
            self.values[idx] = value;
            self.count += 1;
            None
        }
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        let slot = self.locate(key)?;
        if !self.is_occupied(slot) || self.keys[slot] != *key {
            return None;
        }

        self.count -= 1;
        self.ctrl[slot] = Self::DELETED;
        let mut value = V::default();
        std::mem::swap(&mut value, &mut self.values[slot]);
        Some(value)
    }

    fn locate(&self, key: &K) -> Option<usize> {
        if self.ctrl.is_empty() {
            return None;
        }

        let mask = self.bucket_len() - 1;
        let hash = self.hash(key);
        let h2 = self.h2(hash);
        let mut probe = 0usize;
        let mut slot = self.h1(hash);
        let mut first_deleted: Option<usize> = None;

        loop {
            let group_base = slot & !(Self::GROUP - 1);
            for offset in 0..Self::GROUP {
                let idx = (group_base + offset) & mask;
                let ctrl = self.ctrl[idx];
                if ctrl == Self::EMPTY {
                    return Some(first_deleted.unwrap_or(idx));
                }
                if ctrl == Self::DELETED {
                    if first_deleted.is_none() {
                        first_deleted = Some(idx);
                    }
                    continue;
                }
                if ctrl == h2 && self.keys[idx] == *key {
                    return Some(idx);
                }
            }

            probe += 1;
            slot = (slot + probe * Self::GROUP) & mask;
            if probe > mask / Self::GROUP + 1 {
                return None;
            }
        }
    }

    fn rehash(&mut self) {
        let mut tbl = Table::new_with_width(self.bits + 1);
        for i in 0..self.bucket_len() {
            if self.is_occupied(i) {
                let mut key = K::default();
                std::mem::swap(&mut key, &mut self.keys[i]);
                let mut value = V::default();
                std::mem::swap(&mut value, &mut self.values[i]);
                tbl.insert(key, value);
            }
        }
        self.swap(&mut tbl);
    }

    fn hash(&self, key: &K) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.seed.hash(&mut hasher);
        key.hash(&mut hasher);
        hasher.finish()
    }

    #[inline]
    fn bucket_len(&self) -> usize {
        1usize << self.bits
    }

    #[inline]
    fn h1(&self, hash: u64) -> usize {
        (hash as usize) & (self.bucket_len() - 1)
    }

    #[inline]
    fn h2(&self, hash: u64) -> u8 {
        ((hash >> 57) as u8) & 0x7F
    }

    #[inline]
    fn is_occupied(&self, idx: usize) -> bool {
        let c = self.ctrl[idx];
        c != Self::EMPTY && c != Self::DELETED
    }

    const GROUP: usize = 16;
    const EMPTY: u8 = 0xFF;
    const DELETED: u8 = 0x80;
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
