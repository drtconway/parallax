//! A read-only hash table backed by Arrow arrays.
//!
//! This module provides `FrozenTable`, an immutable hash table that uses Apache Arrow
//! arrays for storage. It can be constructed from a mutable `Table` or loaded from
//! Parquet or Arrow IPC (Feather) files.

use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, UInt8Array, UInt32Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::ipc::reader::FileReader as IpcFileReader;
use arrow::ipc::writer::FileWriter as IpcFileWriter;
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use super::table::Table;
use super::swiss;

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
    values: UInt32Array,
}

impl FrozenTable {
    #[allow(dead_code)]
    const GROUP: usize = swiss::PROBE_GROUP;
    #[allow(dead_code)]
    const EMPTY: u8 = swiss::CTRL_EMPTY;
    #[allow(dead_code)]
    const DELETED: u8 = swiss::CTRL_DELETED;

    /// Create a FrozenTable from a mutable Table<u64, u32>.
    ///
    /// The table is "frozen" - no further modifications are possible.
    pub fn from_table(table: Table<u64, u32>) -> Self {
        let count = table.len();
        if count == 0 {
            return Self::empty();
        }

        FrozenTable {
            count: table.count,
            seed: table.seed,
            bits: table.bits,
            ctrl: UInt8Array::from(table.ctrl),
            keys: UInt64Array::from(table.keys),
            values: UInt32Array::from(table.values),
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
            values: UInt32Array::from(Vec::<u32>::new()),
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
        let values = Self::load_u32_array(&dir.join("values.parquet"))?;

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
        Self::save_u32_array(&self.values, &dir.join("values.parquet"))?;

        Ok(())
    }

    /// Load a FrozenTable from Arrow IPC (Feather) files in a directory.
    ///
    /// Expects files named `ctrl.arrow`, `keys.arrow`, `values.arrow`,
    /// and `metadata.arrow` in the given directory.
    pub fn load_from_feather_directory<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();

        // Load metadata
        let metadata_path = dir.join("metadata.arrow");
        let (count, seed, bits) = Self::load_metadata_feather(&metadata_path)?;

        if count == 0 {
            return Ok(Self::empty());
        }

        // Load arrays
        let ctrl = Self::load_u8_array_feather(&dir.join("ctrl.arrow"))?;
        let keys = Self::load_u64_array_feather(&dir.join("keys.arrow"))?;
        let values = Self::load_u32_array_feather(&dir.join("values.arrow"))?;

        Ok(FrozenTable {
            count,
            seed,
            bits,
            ctrl,
            keys,
            values,
        })
    }

    /// Save this FrozenTable to Arrow IPC (Feather) files in a directory.
    ///
    /// Creates files named `ctrl.arrow`, `keys.arrow`, `values.arrow`,
    /// and `metadata.arrow` in the given directory.
    ///
    /// Feather format is generally faster to read/write than Parquet but
    /// produces larger files (no compression by default).
    pub fn save_to_feather_directory<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;

        // Save metadata
        self.save_metadata_feather(&dir.join("metadata.arrow"))?;

        // Save arrays
        Self::save_u8_array_feather(&self.ctrl, &dir.join("ctrl.arrow"))?;
        Self::save_u64_array_feather(&self.keys, &dir.join("keys.arrow"))?;
        Self::save_u32_array_feather(&self.values, &dir.join("values.arrow"))?;

        Ok(())
    }

    /// Get a value by key.
    pub fn get(&self, key: u64) -> Option<u32> {
        if self.count == 0 {
            return None;
        }

        let hash = Self::hash_key(self.seed, key);
        let ctrl = self.ctrl.values();
        let keys = self.keys.values();
        swiss::locate_readonly(ctrl, keys, &key, hash, self.bits)
            .map(|idx| self.values.value(idx) as u32)
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
    #[allow(dead_code)]
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
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
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

    fn load_u32_array<P: AsRef<Path>>(path: P) -> std::io::Result<UInt32Array> {
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let mut reader = builder
            .build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut arrays: Vec<UInt32Array> = Vec::new();
        for batch_result in reader.by_ref() {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let array = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "bad u32 array")
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
                .downcast_ref::<UInt32Array>()
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
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
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
        let schema = Arc::new(Schema::new(vec![Field::new(
            "data",
            DataType::UInt8,
            false,
        )]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .close()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn save_u32_array<P: AsRef<Path>>(array: &UInt32Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "data",
            DataType::UInt32,
            false,
        )]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
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
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .close()
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

        let batch = reader
            .into_iter()
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

    fn save_metadata_feather<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
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
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn load_u8_array_feather<P: AsRef<Path>>(path: P) -> std::io::Result<UInt8Array> {
        let file = File::open(path)?;
        let reader = IpcFileReader::try_new(file, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut arrays: Vec<UInt8Array> = Vec::new();
        for batch_result in reader {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
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

    fn load_u32_array_feather<P: AsRef<Path>>(path: P) -> std::io::Result<UInt32Array> {
        let file = File::open(path)?;
        let reader = IpcFileReader::try_new(file, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut arrays: Vec<UInt32Array> = Vec::new();
        for batch_result in reader {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let array = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "bad u32 array")
                })?
                .clone();
            arrays.push(array);
        }

        if arrays.len() == 1 {
            Ok(arrays.into_iter().next().unwrap())
        } else {
            let refs: Vec<&dyn Array> = arrays.iter().map(|a| a as &dyn Array).collect();
            let concatenated = arrow::compute::concat(&refs)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            Ok(concatenated
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .clone())
        }
    }

    fn load_u64_array_feather<P: AsRef<Path>>(path: P) -> std::io::Result<UInt64Array> {
        let file = File::open(path)?;
        let reader = IpcFileReader::try_new(file, None)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let mut arrays: Vec<UInt64Array> = Vec::new();
        for batch_result in reader {
            let batch = batch_result
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
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

    fn save_u8_array_feather<P: AsRef<Path>>(array: &UInt8Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "data",
            DataType::UInt8,
            false,
        )]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn save_u32_array_feather<P: AsRef<Path>>(array: &UInt32Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "data",
            DataType::UInt32,
            false,
        )]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        Ok(())
    }

    fn save_u64_array_feather<P: AsRef<Path>>(array: &UInt64Array, path: P) -> std::io::Result<()> {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "data",
            DataType::UInt64,
            false,
        )]));

        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(array.clone()) as ArrayRef])
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

        let file = File::create(path)?;
        let mut writer = IpcFileWriter::try_new(file, &schema)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writer
            .finish()
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

        let frozen = FrozenTable::from_table(table);
        assert_eq!(frozen.len(), 100);

        for i in 0..100u64 {
            assert_eq!(frozen.get(i), Some((i * 10) as u32));
        }

        assert_eq!(frozen.get(999), None);
    }

    #[test]
    fn test_empty_table() {
        let table: Table<u64, u32> = Table::new();
        let frozen = FrozenTable::from_table(table);
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

        let frozen = FrozenTable::from_table(table);
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
        let frozen = FrozenTable::from_table(table);
        let dir = tempdir().unwrap();

        frozen.save_to_directory(dir.path()).unwrap();
        let loaded = FrozenTable::load_from_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 0);
        assert!(loaded.is_empty());
    }

    #[test]
    fn test_save_and_load_feather() {
        let mut table: Table<u64, u32> = Table::new();
        for i in 0..1000u64 {
            table.insert(i, (i * 7) as u32);
        }

        let frozen = FrozenTable::from_table(table);
        let dir = tempdir().unwrap();

        frozen.save_to_feather_directory(dir.path()).unwrap();
        let loaded = FrozenTable::load_from_feather_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 1000);
        for i in 0..1000u64 {
            assert_eq!(loaded.get(i), Some((i * 7) as u32));
        }
    }

    #[test]
    fn test_save_and_load_feather_empty() {
        let table: Table<u64, u32> = Table::new();
        let frozen = FrozenTable::from_table(table);
        let dir = tempdir().unwrap();

        frozen.save_to_feather_directory(dir.path()).unwrap();
        let loaded = FrozenTable::load_from_feather_directory(dir.path()).unwrap();

        assert_eq!(loaded.len(), 0);
        assert!(loaded.is_empty());
    }
}
