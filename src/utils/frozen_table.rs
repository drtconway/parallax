//! A read-only hash table backed by Arrow arrays.
//!
//! This module provides `FrozenTable`, an immutable hash table that uses Apache Arrow
//! arrays for storage. It can be constructed from a mutable `Table` or loaded from
//! Parquet files.

use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, UInt64Array, UInt8Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;

use super::table::Table;

/// A read-only hash table backed by Arrow arrays.
///
/// This is an immutable version of `Table` that stores data in Arrow arrays,
/// enabling efficient serialization to/from Parquet format.
pub struct FrozenTable {
    count: usize,
    seed: u64,
    bits: usize,
    ctrl: UInt8Array,
    keys: UInt64Array,
    values: UInt64Array,
}

impl FrozenTable {
    const GROUP: usize = 16;
    const EMPTY: u8 = 0xFF;
    #[allow(dead_code)]
    const DELETED: u8 = 0x80;

    /// Create a FrozenTable from a mutable Table<u64, u32>.
    ///
    /// The table is "frozen" - no further modifications are possible.
    pub fn from_table(table: &Table<u64, u32>) -> Self {
        let count = table.len();
        if count == 0 {
            return Self::empty();
        }

        // Extract the internal vectors from table via iteration
        // We need to rebuild the hash table structure
        let entries: Vec<(u64, u32)> = table.iter().map(|(&k, &v)| (k, v)).collect();

        // Determine the size needed
        let bits = Self::compute_bits(entries.len());
        let n = 1usize << bits;

        // Build the hash table arrays
        let mut ctrl = vec![Self::EMPTY; n + Self::GROUP];
        let mut keys = vec![0u64; n];
        let mut values = vec![0u64; n];

        let seed = 0xcbf2_9ce4_8422_2325u64;

        for (key, value) in entries {
            let hash = Self::hash_key(seed, key);
            let slot = Self::find_slot(&ctrl, &keys, key, hash, bits);
            ctrl[slot] = Self::h2(hash);
            keys[slot] = key;
            values[slot] = value as u64;
        }

        FrozenTable {
            count,
            seed,
            bits,
            ctrl: UInt8Array::from(ctrl),
            keys: UInt64Array::from(keys),
            values: UInt64Array::from(values),
        }
    }

    /// Create an empty FrozenTable.
    pub fn empty() -> Self {
        FrozenTable {
            count: 0,
            seed: 0xcbf2_9ce4_8422_2325,
            bits: 0,
            ctrl: UInt8Array::from(Vec::<u8>::new()),
            keys: UInt64Array::from(Vec::<u64>::new()),
            values: UInt64Array::from(Vec::<u64>::new()),
        }
    }

    /// Load a FrozenTable from Parquet files in a directory.
    ///
    /// Expects files named `ctrl.parquet`, `keys.parquet`, `values.parquet`,
    /// and `metadata.parquet` in the given directory.
    pub fn load_from_directory<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();

        // Load metadata
        let metadata_path = dir.join("metadata.parquet");
        let (count, seed, bits) = Self::load_metadata(&metadata_path)?;

        if count == 0 {
            return Ok(Self::empty());
        }

        // Load arrays
        let ctrl = Self::load_u8_array(&dir.join("ctrl.parquet"))?;
        let keys = Self::load_u64_array(&dir.join("keys.parquet"))?;
        let values = Self::load_u64_array(&dir.join("values.parquet"))?;

        Ok(FrozenTable {
            count,
            seed,
            bits,
            ctrl,
            keys,
            values,
        })
    }

    /// Save this FrozenTable to Parquet files in a directory.
    ///
    /// Creates files named `ctrl.parquet`, `keys.parquet`, `values.parquet`,
    /// and `metadata.parquet` in the given directory.
    pub fn save_to_directory<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;

        // Save metadata
        self.save_metadata(&dir.join("metadata.parquet"))?;

        // Save arrays
        Self::save_u8_array(&self.ctrl, &dir.join("ctrl.parquet"))?;
        Self::save_u64_array(&self.keys, &dir.join("keys.parquet"))?;
        Self::save_u64_array(&self.values, &dir.join("values.parquet"))?;

        Ok(())
    }

    /// Get a value by key.
    pub fn get(&self, key: u64) -> Option<u32> {
        if self.count == 0 {
            return None;
        }

        let hash = Self::hash_key(self.seed, key);
        let h2 = Self::h2(hash);
        let mask = (1usize << self.bits) - 1;
        let mut probe = 0usize;
        let mut slot = (hash as usize) & mask;

        loop {
            let group_base = slot & !(Self::GROUP - 1);
            for offset in 0..Self::GROUP {
                let idx = (group_base + offset) & mask;
                let ctrl = self.ctrl.value(idx);
                if ctrl == Self::EMPTY {
                    return None;
                }
                if ctrl == h2 && self.keys.value(idx) == key {
                    return Some(self.values.value(idx) as u32);
                }
            }

            probe += 1;
            slot = (slot + probe * Self::GROUP) & mask;
            if probe > mask / Self::GROUP + 1 {
                return None;
            }
        }
    }

    /// Check if a key exists.
    #[allow(dead_code)]
    pub fn contains_key(&self, key: u64) -> bool {
        self.get(key).is_some()
    }

    /// Get the number of entries.
    pub fn len(&self) -> usize {
        self.count
    }

    /// Check if the table is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    // Internal helper functions

    fn compute_bits(count: usize) -> usize {
        // Target ~75% load factor
        let needed = (count * 4 / 3).max(16);
        (usize::BITS - needed.leading_zeros()) as usize
    }

    fn hash_key(seed: u64, key: u64) -> u64 {
        let mut hasher = DefaultHasher::new();
        seed.hash(&mut hasher);
        key.hash(&mut hasher);
        hasher.finish()
    }

    #[inline]
    fn h2(hash: u64) -> u8 {
        ((hash >> 57) as u8) & 0x7F
    }

    fn find_slot(ctrl: &[u8], keys: &[u64], key: u64, hash: u64, bits: usize) -> usize {
        let mask = (1usize << bits) - 1;
        let h2 = Self::h2(hash);
        let mut probe = 0usize;
        let mut slot = (hash as usize) & mask;

        loop {
            let group_base = slot & !(Self::GROUP - 1);
            for offset in 0..Self::GROUP {
                let idx = (group_base + offset) & mask;
                let c = ctrl[idx];
                if c == Self::EMPTY {
                    return idx;
                }
                if c == h2 && keys[idx] == key {
                    return idx;
                }
            }

            probe += 1;
            slot = (slot + probe * Self::GROUP) & mask;
            if probe > mask / Self::GROUP + 1 {
                panic!("FrozenTable: no slot found (table too full)");
            }
        }
    }

    fn load_metadata<P: AsRef<Path>>(path: P) -> std::io::Result<(usize, u64, usize)> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let mut reader = builder
            .build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let batch = reader
            .next()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "empty metadata"))?
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let count = batch
            .column(0)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad count"))?
            .value(0) as usize;

        let seed = batch
            .column(1)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad seed"))?
            .value(0);

        let bits = batch
            .column(2)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "bad bits"))?
            .value(0) as usize;

        Ok((count, seed, bits))
    }

    fn save_metadata<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![
            Field::new("count", DataType::UInt64, false),
            Field::new("seed", DataType::UInt64, false),
            Field::new("bits", DataType::UInt64, false),
        ]));

        let count_array = UInt64Array::from(vec![self.count as u64]);
        let seed_array = UInt64Array::from(vec![self.seed]);
        let bits_array = UInt64Array::from(vec![self.bits as u64]);

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(count_array) as ArrayRef,
                Arc::new(seed_array) as ArrayRef,
                Arc::new(bits_array) as ArrayRef,
            ],
        )
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = ArrowWriter::try_new(file, schema, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn load_u8_array<P: AsRef<Path>>(path: P) -> std::io::Result<UInt8Array> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let mut reader = builder
            .build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut arrays: Vec<UInt8Array> = Vec::new();
        for batch_result in reader.by_ref() {
            let batch =
                batch_result.map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let array = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt8Array>()
                .ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "bad u8 array")
                })?
                .clone();
            arrays.push(array);
        }

        // Concatenate all batches
        if arrays.len() == 1 {
            Ok(arrays.into_iter().next().unwrap())
        } else {
            let refs: Vec<&dyn Array> = arrays.iter().map(|a| a as &dyn Array).collect();
            let concatenated = arrow::compute::concat(&refs)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            Ok(concatenated
                .as_any()
                .downcast_ref::<UInt8Array>()
                .unwrap()
                .clone())
        }
    }

    fn load_u64_array<P: AsRef<Path>>(path: P) -> std::io::Result<UInt64Array> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let mut reader = builder
            .build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut arrays: Vec<UInt64Array> = Vec::new();
        for batch_result in reader.by_ref() {
            let batch =
                batch_result.map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let array = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "bad u64 array")
                })?
                .clone();
            arrays.push(array);
        }

        // Concatenate all batches
        if arrays.len() == 1 {
            Ok(arrays.into_iter().next().unwrap())
        } else {
            let refs: Vec<&dyn Array> = arrays.iter().map(|a| a as &dyn Array).collect();
            let concatenated = arrow::compute::concat(&refs)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            Ok(concatenated
                .as_any()
                .downcast_ref::<UInt64Array>()
                .unwrap()
                .clone())
        }
    }

    fn save_u8_array<P: AsRef<Path>>(array: &UInt8Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new("data", DataType::UInt8, false)]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = ArrowWriter::try_new(file, schema, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn save_u64_array<P: AsRef<Path>>(array: &UInt64Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "data",
            DataType::UInt64,
            false,
        )]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = ArrowWriter::try_new(file, schema, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_from_table_and_get() {
        let mut table: Table<u64, u32> = Table::new();
        for i in 0..100u64 {
            table.insert(i, (i * 10) as u32);
        }

        let frozen = FrozenTable::from_table(&table);
        assert_eq!(frozen.len(), 100);

        for i in 0..100u64 {
            assert_eq!(frozen.get(i), Some((i * 10) as u32));
        }

        assert_eq!(frozen.get(999), None);
    }

    #[test]
    fn test_empty_table() {
        let table: Table<u64, u32> = Table::new();
        let frozen = FrozenTable::from_table(&table);
        assert_eq!(frozen.len(), 0);
        assert!(frozen.is_empty());
        assert_eq!(frozen.get(42), None);
    }

    #[test]
    fn test_save_and_load() {
        let mut table: Table<u64, u32> = Table::new();
        for i in 0..1000u64 {
            table.insert(i, (i * 7) as u32);
        }

        let frozen = FrozenTable::from_table(&table);
        let dir = tempdir().unwrap();

        frozen.save_to_directory(dir.path()).unwrap();
        let loaded = FrozenTable::load_from_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 1000);
        for i in 0..1000u64 {
            assert_eq!(loaded.get(i), Some((i * 7) as u32));
        }
    }

    #[test]
    fn test_save_and_load_empty() {
        let table: Table<u64, u32> = Table::new();
        let frozen = FrozenTable::from_table(&table);
        let dir = tempdir().unwrap();

        frozen.save_to_directory(dir.path()).unwrap();
        let loaded = FrozenTable::load_from_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 0);
        assert!(loaded.is_empty());
    }
}
