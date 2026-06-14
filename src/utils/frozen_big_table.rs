//! Arrow-backed immutable hash table mapping u64 keys to multiple u64 values.
//!
//! This is similar to `FrozenTable` but supports multiple values per key,
//! constructed from a `HashMap<u64, Vec<u64>>` or loaded from Parquet files.

use arrow::array::{Array, UInt8Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::ipc::reader::FileReader as IpcFileReader;
use arrow::ipc::writer::FileWriter as IpcFileWriter;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use super::swiss;

/// An immutable hash table backed by Arrow arrays, mapping keys to multiple values.
/// 
/// Each key maps to a slice of u64 values. The table can be serialized to and
/// deserialized from Parquet files for persistent storage.
pub struct FrozenBigTable {
    /// Number of entries in the table.
    count: usize,
    /// Hash seed for consistent hashing.
    seed: u64,
    /// Number of bits used for indexing (capacity = 1 << bits).
    bits: usize,
    /// Control bytes for probing.
    ctrl: UInt8Array,
    /// Keys array (same length as ctrl).
    keys: UInt64Array,
    /// Offsets into values array (length = capacity + 1).
    /// Slot i's values are at values[offsets[i]..offsets[i+1]].
    offsets: UInt64Array,
    /// All values concatenated.
    values: UInt64Array,
}

impl FrozenBigTable {
    /// Create a FrozenBigTable from a HashMap.
    pub fn from_hashmap(map: HashMap<u64, Vec<u64>>) -> Self {
        if map.is_empty() {
            return Self::empty();
        }

        let count = map.len();
        let bits = Self::compute_bits(count);
        let capacity = 1usize << bits;
        let seed = 0xcbf2_9ce4_8422_2325u64;

        // Initialize arrays
        let mut ctrl_vec = vec![swiss::CTRL_EMPTY; capacity];
        let mut keys_vec = vec![0u64; capacity];
        let mut offsets_vec = vec![0u64; capacity + 1];
        
        // First pass: place keys and compute value counts per slot
        let mut slot_values: Vec<Option<Vec<u64>>> = vec![None; capacity];
        
        for (key, values) in map {
            let hash = Self::hash(key, seed);
            let slot = swiss::find_empty_slot(&ctrl_vec, hash, bits);
            
            ctrl_vec[slot] = swiss::h2(hash);
            keys_vec[slot] = key;
            slot_values[slot] = Some(values);
        }

        // Second pass: compute offsets
        let mut offset = 0u64;
        for i in 0..capacity {
            offsets_vec[i] = offset;
            if let Some(vals) = slot_values[i].as_ref() {
                offset += vals.len() as u64;
            }
        }
        offsets_vec[capacity] = offset;

        // Third pass: collect all values in slot order
        let total_values = offset as usize;
        let mut values_vec = Vec::with_capacity(total_values);
        for slot_val in slot_values {
            if let Some(vals) = slot_val {
                values_vec.extend_from_slice(&vals);
            }
        }

        FrozenBigTable {
            count,
            seed,
            bits,
            ctrl: UInt8Array::from(ctrl_vec),
            keys: UInt64Array::from(keys_vec),
            offsets: UInt64Array::from(offsets_vec),
            values: UInt64Array::from(values_vec),
        }
    }

    /// Create an empty FrozenBigTable.
    pub fn empty() -> Self {
        FrozenBigTable {
            count: 0,
            seed: 0xcbf2_9ce4_8422_2325,
            bits: 0,
            ctrl: UInt8Array::from(Vec::<u8>::new()),
            keys: UInt64Array::from(Vec::<u64>::new()),
            offsets: UInt64Array::from(vec![0u64]),
            values: UInt64Array::from(Vec::<u64>::new()),
        }
    }

    /// Load a FrozenBigTable from Parquet files in a directory.
    ///
    /// Expects files named `ctrl.parquet`, `keys.parquet`, `offsets.parquet`,
    /// `values.parquet`, and `metadata.parquet` in the given directory.
    pub fn load_from_directory<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();

        let (count, seed, bits) = Self::load_metadata(dir.join("metadata.parquet"))?;
        
        if count == 0 {
            return Ok(Self::empty());
        }

        let ctrl = Self::load_u8_array(dir.join("ctrl.parquet"))?;
        let keys = Self::load_u64_array(dir.join("keys.parquet"))?;
        let offsets = Self::load_u64_array(dir.join("offsets.parquet"))?;
        let values = Self::load_u64_array(dir.join("values.parquet"))?;

        Ok(FrozenBigTable {
            count,
            seed,
            bits,
            ctrl,
            keys,
            offsets,
            values,
        })
    }

    /// Save this FrozenBigTable to Parquet files in a directory.
    ///
    /// Creates files named `ctrl.parquet`, `keys.parquet`, `offsets.parquet`,
    /// `values.parquet`, and `metadata.parquet`.
    pub fn save_to_directory<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;

        self.save_metadata(dir.join("metadata.parquet"))?;

        if self.count == 0 {
            return Ok(());
        }

        Self::save_u8_array(&self.ctrl, dir.join("ctrl.parquet"))?;
        Self::save_u64_array(&self.keys, dir.join("keys.parquet"))?;
        Self::save_u64_array(&self.offsets, dir.join("offsets.parquet"))?;
        Self::save_u64_array(&self.values, dir.join("values.parquet"))?;

        Ok(())
    }

    /// Load a FrozenBigTable from Arrow IPC (Feather) files in a directory.
    ///
    /// Expects files named `ctrl.arrow`, `keys.arrow`, `offsets.arrow`,
    /// `values.arrow`, and `metadata.arrow` in the given directory.
    pub fn load_from_feather_directory<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();

        let (count, seed, bits) = Self::load_metadata_feather(dir.join("metadata.arrow"))?;
        
        if count == 0 {
            return Ok(Self::empty());
        }

        let ctrl = Self::load_u8_array_feather(dir.join("ctrl.arrow"))?;
        let keys = Self::load_u64_array_feather(dir.join("keys.arrow"))?;
        let offsets = Self::load_u64_array_feather(dir.join("offsets.arrow"))?;
        let values = Self::load_u64_array_feather(dir.join("values.arrow"))?;

        Ok(FrozenBigTable {
            count,
            seed,
            bits,
            ctrl,
            keys,
            offsets,
            values,
        })
    }

    /// Save this FrozenBigTable to Arrow IPC (Feather) files in a directory.
    ///
    /// Creates files named `ctrl.arrow`, `keys.arrow`, `offsets.arrow`,
    /// `values.arrow`, and `metadata.arrow`.
    ///
    /// Feather format is generally faster to read/write than Parquet but
    /// produces larger files (no compression by default).
    pub fn save_to_feather_directory<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;

        self.save_metadata_feather(dir.join("metadata.arrow"))?;

        if self.count == 0 {
            return Ok(());
        }

        Self::save_u8_array_feather(&self.ctrl, dir.join("ctrl.arrow"))?;
        Self::save_u64_array_feather(&self.keys, dir.join("keys.arrow"))?;
        Self::save_u64_array_feather(&self.offsets, dir.join("offsets.arrow"))?;
        Self::save_u64_array_feather(&self.values, dir.join("values.arrow"))?;

        Ok(())
    }

    /// Look up a key and return its values as a slice.
    pub fn get(&self, key: u64) -> Option<&[u64]> {
        if self.count == 0 {
            return None;
        }

        let hash = Self::hash(key, self.seed);
        let ctrl_slice = self.ctrl.values();
        let keys_slice = self.keys.values();
        let offsets_slice = self.offsets.values();
        let values_slice = self.values.values();

        let slot = swiss::locate_readonly(ctrl_slice, keys_slice, &key, hash, self.bits)?;
        
        let start = offsets_slice[slot] as usize;
        let end = offsets_slice[slot + 1] as usize;
        
        Some(&values_slice[start..end])
    }

    /// Issue prefetch hints for the cache lines that will be touched when
    /// looking up `key`. Call this well ahead of the corresponding `get()`
    /// to allow the memory subsystem to bring the data into L1.
    #[inline]
    pub fn prefetch_key(&self, key: u64) {
        if self.count == 0 || self.bits == 0 {
            return;
        }
        let hash = Self::hash(key, self.seed);
        let (group_base, _mask) = swiss::probe_position(hash, self.bits);
        unsafe {
            let ctrl_ptr = self.ctrl.values().as_ptr().add(group_base) as *const u8;
            let keys_ptr = self.keys.values().as_ptr().add(group_base) as *const u8;
            #[cfg(target_arch = "x86_64")]
            {
                std::arch::x86_64::_mm_prefetch(ctrl_ptr as *const i8, std::arch::x86_64::_MM_HINT_T0);
                std::arch::x86_64::_mm_prefetch(keys_ptr as *const i8, std::arch::x86_64::_MM_HINT_T0);
            }
            #[cfg(target_arch = "aarch64")]
            {
                std::arch::aarch64::_prefetch(ctrl_ptr as *const i8, std::arch::aarch64::_PREFETCH_READ, std::arch::aarch64::_PREFETCH_LOCALITY3);
                std::arch::aarch64::_prefetch(keys_ptr as *const i8, std::arch::aarch64::_PREFETCH_READ, std::arch::aarch64::_PREFETCH_LOCALITY3);
            }
        }
    }

    /// Returns the number of keys in the table.
    pub fn len(&self) -> usize {
        self.count
    }

    /// Iterate over all (key, slot) pairs in hash-table order.
    pub fn iter(&self) -> FrozenBigTableIter<'_> {
        FrozenBigTableIter { table: self, pos: 0 }
    }

    /// Return the values for `slot` as a slice into the Arrow buffer.
    pub fn loci_as_slice(&self, slot: usize) -> &[u64] {
        let start = self.offsets.value(slot) as usize;
        let end = self.offsets.value(slot + 1) as usize;
        &self.values.values()[start..end]
    }

    // --- Private helper methods ---

    fn compute_bits(count: usize) -> usize {
        if count == 0 {
            return 0;
        }
        // Target ~70% load factor
        let needed = (count * 10 / 7).next_power_of_two();
        needed.trailing_zeros() as usize
    }

    fn hash(key: u64, seed: u64) -> u64 {
        // FxHash-style mixing
        let mut h = key.wrapping_mul(0x517cc1b727220a95);
        h ^= seed;
        h.wrapping_mul(0x517cc1b727220a95)
    }

    // --- Parquet I/O helpers ---

    fn load_metadata<P: AsRef<Path>>(path: P) -> std::io::Result<(usize, u64, usize)> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut reader = builder.build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let batch = reader.next()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty metadata"))?
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let count_col = batch.column(0)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid count column"))?;
        
        let seed_col = batch.column(1)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid seed column"))?;
        
        let bits_col = batch.column(2)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid bits column"))?;
        
        Ok((count_col.value(0) as usize, seed_col.value(0), bits_col.value(0) as usize))
    }

    fn save_metadata<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let schema = Schema::new(vec![
            Field::new("count", DataType::UInt64, false),
            Field::new("seed", DataType::UInt64, false),
            Field::new("bits", DataType::UInt64, false),
        ]);

        let count_arr = UInt64Array::from(vec![self.count as u64]);
        let seed_arr = UInt64Array::from(vec![self.seed]);
        let bits_arr = UInt64Array::from(vec![self.bits as u64]);

        let batch = RecordBatch::try_new(
            Arc::new(schema),
            vec![Arc::new(count_arr), Arc::new(seed_arr), Arc::new(bits_arr)],
        ).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        writer.write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn load_u8_array<P: AsRef<Path>>(path: P) -> std::io::Result<UInt8Array> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut reader = builder.build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut all_values = Vec::new();
        for batch_result in reader.by_ref() {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let col = batch.column(0)
                .as_any()
                .downcast_ref::<UInt8Array>()
                .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid u8 array"))?;
            all_values.extend(col.values().iter().copied());
        }
        
        Ok(UInt8Array::from(all_values))
    }

    fn load_u64_array<P: AsRef<Path>>(path: P) -> std::io::Result<UInt64Array> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut reader = builder.build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut all_values = Vec::new();
        for batch_result in reader.by_ref() {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let col = batch.column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid u64 array"))?;
            all_values.extend(col.values().iter().copied());
        }
        
        Ok(UInt64Array::from(all_values))
    }

    fn save_u8_array<P: AsRef<Path>>(array: &UInt8Array, path: P) -> std::io::Result<()> {
        let schema = Schema::new(vec![Field::new("data", DataType::UInt8, false)]);
        let batch = RecordBatch::try_new(Arc::new(schema), vec![Arc::new(array.clone())])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        writer.write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn save_u64_array<P: AsRef<Path>>(array: &UInt64Array, path: P) -> std::io::Result<()> {
        let schema = Schema::new(vec![Field::new("data", DataType::UInt64, false)]);
        let batch = RecordBatch::try_new(Arc::new(schema), vec![Arc::new(array.clone())])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        writer.write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    // =========================================================================
    // Arrow IPC (Feather) helper functions
    // =========================================================================

    fn load_metadata_feather<P: AsRef<Path>>(path: P) -> std::io::Result<(usize, u64, usize)> {
        let file = File::open(path)?;
        let reader = IpcFileReader::try_new(file, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let batch = reader.into_iter().next()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Empty metadata"))?
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let count_col = batch.column(0).as_any().downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid count column"))?;
        let seed_col = batch.column(1).as_any().downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid seed column"))?;
        let bits_col = batch.column(2).as_any().downcast_ref::<UInt64Array>()
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid bits column"))?;
        
        Ok((count_col.value(0) as usize, seed_col.value(0), bits_col.value(0) as usize))
    }

    fn save_metadata_feather<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![
            Field::new("count", DataType::UInt64, false),
            Field::new("seed", DataType::UInt64, false),
            Field::new("bits", DataType::UInt64, false),
        ]));

        let count_arr = UInt64Array::from(vec![self.count as u64]);
        let seed_arr = UInt64Array::from(vec![self.seed]);
        let bits_arr = UInt64Array::from(vec![self.bits as u64]);

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![Arc::new(count_arr), Arc::new(seed_arr), Arc::new(bits_arr)],
        ).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn load_u8_array_feather<P: AsRef<Path>>(path: P) -> std::io::Result<UInt8Array> {
        let file = File::open(path)?;
        let reader = IpcFileReader::try_new(file, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut all_values = Vec::new();
        for batch_result in reader {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let col = batch.column(0).as_any().downcast_ref::<UInt8Array>()
                .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid u8 array"))?;
            all_values.extend(col.values().iter().copied());
        }
        
        Ok(UInt8Array::from(all_values))
    }

    fn load_u64_array_feather<P: AsRef<Path>>(path: P) -> std::io::Result<UInt64Array> {
        let file = File::open(path)?;
        let reader = IpcFileReader::try_new(file, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        let mut all_values = Vec::new();
        for batch_result in reader {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let col = batch.column(0).as_any().downcast_ref::<UInt64Array>()
                .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid u64 array"))?;
            all_values.extend(col.values().iter().copied());
        }
        
        Ok(UInt64Array::from(all_values))
    }

    fn save_u8_array_feather<P: AsRef<Path>>(array: &UInt8Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new("data", DataType::UInt8, false)]));
        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone())])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn save_u64_array_feather<P: AsRef<Path>>(array: &UInt64Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new("data", DataType::UInt64, false)]));
        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone())])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer.finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }
}

pub struct FrozenBigTableIter<'a> {
    table: &'a FrozenBigTable,
    pos: usize,
}

impl<'a> Iterator for FrozenBigTableIter<'a> {
    type Item = (u64, usize); // (key, slot)

    fn next(&mut self) -> Option<Self::Item> {
        let ctrl = self.table.ctrl.values();
        let capacity = ctrl.len();
        while self.pos < capacity {
            let i = self.pos;
            self.pos += 1;
            if swiss::is_occupied(ctrl[i]) {
                return Some((self.table.keys.values()[i], i));
            }
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.table.count))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_from_hashmap_and_get() {
        let mut map = HashMap::new();
        map.insert(100u64, vec![1u64, 2, 3]);
        map.insert(200u64, vec![4u64, 5]);
        map.insert(300u64, vec![6u64]);

        let frozen = FrozenBigTable::from_hashmap(map);

        assert_eq!(frozen.len(), 3);

        assert_eq!(frozen.get(100), Some(&[1u64, 2, 3][..]));
        assert_eq!(frozen.get(200), Some(&[4u64, 5][..]));
        assert_eq!(frozen.get(300), Some(&[6u64][..]));
        assert_eq!(frozen.get(999), None);
    }

    #[test]
    fn test_empty_table() {
        let map: HashMap<u64, Vec<u64>> = HashMap::new();
        let frozen = FrozenBigTable::from_hashmap(map);

        assert_eq!(frozen.len(), 0);
        assert_eq!(frozen.get(100), None);
    }

    #[test]
    fn test_empty_values() {
        let mut map = HashMap::new();
        map.insert(100u64, vec![1u64, 2]);
        map.insert(200u64, Vec::new()); // Empty value list

        let frozen = FrozenBigTable::from_hashmap(map);

        assert_eq!(frozen.len(), 2);
        assert_eq!(frozen.get(100), Some(&[1u64, 2][..]));
        assert_eq!(frozen.get(200), Some(&[][..])); // Empty slice
    }

    #[test]
    fn test_save_and_load() {
        let dir = tempdir().unwrap();

        let mut map = HashMap::new();
        map.insert(100u64, vec![1u64, 2, 3]);
        map.insert(200u64, vec![4u64, 5]);
        map.insert(300u64, vec![6u64, 7, 8, 9]);

        let original = FrozenBigTable::from_hashmap(map);
        original.save_to_directory(dir.path()).unwrap();

        let loaded = FrozenBigTable::load_from_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), original.len());
        assert_eq!(loaded.get(100), Some(&[1u64, 2, 3][..]));
        assert_eq!(loaded.get(200), Some(&[4u64, 5][..]));
        assert_eq!(loaded.get(300), Some(&[6u64, 7, 8, 9][..]));
        assert_eq!(loaded.get(999), None);
    }

    #[test]
    fn test_save_and_load_empty() {
        let dir = tempdir().unwrap();

        let frozen = FrozenBigTable::empty();
        frozen.save_to_directory(dir.path()).unwrap();

        let loaded = FrozenBigTable::load_from_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 0);
        assert_eq!(loaded.get(100), None);
    }

    #[test]
    fn test_many_entries() {
        let mut map = HashMap::new();
        for i in 0..1000u64 {
            map.insert(i, vec![i * 2, i * 2 + 1]);
        }

        let frozen = FrozenBigTable::from_hashmap(map);

        assert_eq!(frozen.len(), 1000);

        for i in 0..1000u64 {
            let expected = [i * 2, i * 2 + 1];
            assert_eq!(frozen.get(i), Some(&expected[..]));
        }
    }

    #[test]
    fn test_save_and_load_feather() {
        let dir = tempdir().unwrap();

        let mut map = HashMap::new();
        map.insert(100u64, vec![1u64, 2, 3]);
        map.insert(200u64, vec![4u64, 5]);
        map.insert(300u64, vec![6u64, 7, 8, 9]);

        let original = FrozenBigTable::from_hashmap(map);
        original.save_to_feather_directory(dir.path()).unwrap();

        let loaded = FrozenBigTable::load_from_feather_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), original.len());
        assert_eq!(loaded.get(100), Some(&[1u64, 2, 3][..]));
        assert_eq!(loaded.get(200), Some(&[4u64, 5][..]));
        assert_eq!(loaded.get(300), Some(&[6u64, 7, 8, 9][..]));
        assert_eq!(loaded.get(999), None);
    }

    #[test]
    fn test_save_and_load_feather_empty() {
        let dir = tempdir().unwrap();

        let frozen = FrozenBigTable::empty();
        frozen.save_to_feather_directory(dir.path()).unwrap();

        let loaded = FrozenBigTable::load_from_feather_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 0);
        assert_eq!(loaded.get(100), None);
    }
}
