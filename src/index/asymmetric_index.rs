use std::cmp::Reverse;
use std::collections::HashMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::index::PackedLocus;
use crate::index::SyncmerIndex;
use crate::kmers::Kmer;
use crate::reference::ChromInfo;
use crate::reference::InMemoryReference;
use crate::utils::Selection;
use crate::utils::frozen_big_table::FrozenBigTable;
use crate::utils::frozen_table::FrozenTable;
use crate::utils::hasher::FnvHasher;
use crate::utils::hasher::Hasher;
use crate::utils::hasher::Splitmix64Hasher;
use crate::utils::pool::Pool;
use crate::utils::table::Table;

const INDEX_VERSION: u32 = 1;
const INDEX_TYPE: &str = "asymmetric-syncmer";

#[derive(Debug, Serialize, Deserialize)]
pub struct IndexMetadata {
    pub version: u32,
    pub index_type: String,
    pub k: usize,
    pub s: usize,
    pub portable: bool,
    /// Native byte order at build time; only meaningful when portable == false.
    #[serde(default)]
    pub endian: String,
}

impl IndexMetadata {
    fn native_endian() -> &'static str {
        if cfg!(target_endian = "little") {
            "little"
        } else {
            "big"
        }
    }

    fn new_portable(k: usize, s: usize) -> Self {
        Self {
            version: INDEX_VERSION,
            index_type: INDEX_TYPE.to_string(),
            k,
            s,
            portable: true,
            endian: String::new(),
        }
    }

    fn new_native(k: usize, s: usize) -> Self {
        Self {
            version: INDEX_VERSION,
            index_type: INDEX_TYPE.to_string(),
            k,
            s,
            portable: false,
            endian: Self::native_endian().to_string(),
        }
    }

    fn load<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        serde_json::from_str(&content)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    fn save<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let content = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        std::fs::write(path, content)
    }

    fn validate<const K: usize, const S: usize>(&self) -> std::io::Result<()> {
        if self.version != INDEX_VERSION {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!(
                    "Index version mismatch: index has version {}, expected {}",
                    self.version, INDEX_VERSION
                ),
            ));
        }
        if self.index_type != INDEX_TYPE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!(
                    "Index type mismatch: index has type \"{}\", expected \"{}\"",
                    self.index_type, INDEX_TYPE
                ),
            ));
        }
        if self.k != K {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!(
                    "Index k-mer size mismatch: index has k={}, binary expects k={}",
                    self.k, K
                ),
            ));
        }
        if self.s != S {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!(
                    "Index syncmer size mismatch: index has s={}, binary expects s={}",
                    self.s, S
                ),
            ));
        }
        if !self.portable && self.endian != Self::native_endian() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!(
                    "Index endian mismatch: index was built on a {}-endian system, \
                     but this system is {}-endian. Use a portable index instead.",
                    self.endian,
                    Self::native_endian()
                ),
            ));
        }
        Ok(())
    }
}

// repr(transparent) guarantees Locus has the same memory layout as u64.
// This allows lookup_batch to reinterpret the raw &[u64] locus storage as
// &[Locus] without copying. Removing this attribute would make that cast
// undefined behaviour.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Locus(u64);

impl From<u64> for Locus {
    fn from(value: u64) -> Self {
        Locus(value)
    }
}

impl Into<u64> for Locus {
    fn into(self) -> u64 {
        self.0
    }
}

impl PackedLocus for Locus {
    fn pack(chrom: usize, pos: usize, strand: super::Strand) -> Self {
        let chrom: u64 = chrom as u64;
        let pos: u64 = pos as u64;
        let strand = if strand.is_reverse() { 1u64 } else { 0u64 };
        Locus(chrom << 32 | pos << 1 | strand)
    }

    fn unpack(&self) -> (usize, usize, super::Strand) {
        let chrom = (self.0 >> 32) as usize;
        let pos = ((self.0 & 0xFFFFFFFF) >> 1) as usize;
        let strand = if self.0 & 1 == 1 {
            super::Strand::Reverse
        } else {
            super::Strand::Forward
        };
        (chrom, pos, strand)
    }
}

/// Immutable K-mer seed index backed by Arrow arrays.
///
/// This is the production index structure, optimized for fast lookups and
/// supporting save/load to Parquet files for persistence.
pub struct AsymmetricIndex<const K: usize, const S: usize> {
    chrom_info: Vec<ChromInfo>,
    unique_seeds: FrozenTable,
    nonunique_seeds: FrozenBigTable,
    query_buffers: Pool<Vec<(usize, u64)>>
}

impl<const K: usize, const S: usize> AsymmetricIndex<K, S> {
    /// Load an index from a directory, dispatching on the metadata to select
    /// the correct format (portable Feather or native binary).
    ///
    /// Validates k, s, version, index type, and (for non-portable indexes)
    /// endianness before loading any seed data.
    pub fn load<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();
        let meta = IndexMetadata::load(dir.join("metadata.json"))?;
        meta.validate::<K, S>()?;

        if meta.portable {
            Self::load_feather_inner(dir)
        } else {
            Self::load_native_inner(dir)
        }
    }

    fn load_feather_inner(dir: &Path) -> std::io::Result<Self> {
        let now = std::time::Instant::now();
        let chrom_info = Self::load_chrom_info(dir.join("chrom_info.json"))?;
        let unique_seeds = FrozenTable::load_from_feather_directory(dir.join("unique_seeds"))?;
        let nonunique_seeds =
            FrozenBigTable::load_from_feather_directory(dir.join("nonunique_seeds"))?;
        log::info!(
            "Loaded index (portable): {} chromosomes, {} unique seeds, {} nonunique seeds in {:.2}s",
            chrom_info.len(),
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );
        Ok(AsymmetricIndex {
            chrom_info,
            unique_seeds,
            nonunique_seeds,
            query_buffers: Pool::new()
        })
    }

    fn load_native_inner(dir: &Path) -> std::io::Result<Self> {
        let now = std::time::Instant::now();
        let chrom_info = Self::load_chrom_info(dir.join("chrom_info.json"))?;
        let unique_seeds = FrozenTable::load_from_directory(dir.join("unique_seeds"))?;
        let nonunique_seeds = FrozenBigTable::load_from_directory(dir.join("nonunique_seeds"))?;
        log::info!(
            "Loaded index (native): {} chromosomes, {} unique seeds, {} nonunique seeds in {:.2}s",
            chrom_info.len(),
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );
        Ok(AsymmetricIndex {
            chrom_info,
            unique_seeds,
            nonunique_seeds,
            query_buffers: Pool::new()
        })
    }

    /// Save the index in portable (Feather/Arrow IPC) format.
    pub fn save_portable<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;
        IndexMetadata::new_portable(K, S).save(dir.join("metadata.json"))?;
        Self::save_chrom_info(&self.chrom_info, dir.join("chrom_info.json"))?;
        self.unique_seeds
            .save_to_feather_directory(dir.join("unique_seeds"))?;
        self.nonunique_seeds
            .save_to_feather_directory(dir.join("nonunique_seeds"))?;
        log::info!(
            "Saved index (portable): {} chromosomes, {} unique seeds, {} nonunique seeds",
            self.chrom_info.len(),
            self.unique_seeds.len(),
            self.nonunique_seeds.len()
        );
        Ok(())
    }

    /// Save the index in native (endian-specific binary) format.
    pub fn save_native<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;
        IndexMetadata::new_native(K, S).save(dir.join("metadata.json"))?;
        Self::save_chrom_info(&self.chrom_info, dir.join("chrom_info.json"))?;
        self.unique_seeds
            .save_to_directory(dir.join("unique_seeds"))?;
        self.nonunique_seeds
            .save_to_directory(dir.join("nonunique_seeds"))?;
        log::info!(
            "Saved index (native): {} chromosomes, {} unique seeds, {} nonunique seeds",
            self.chrom_info.len(),
            self.unique_seeds.len(),
            self.nonunique_seeds.len()
        );
        Ok(())
    }

    /// Issue prefetch hints for both hash tables for the given k-mer.
    /// Call this ~PIPE iterations before the corresponding `with()` call.
    #[inline]
    pub fn prefetch_kmer(&self, kmer: u64) {
        self.unique_seeds.prefetch_key(kmer);
        // Also prefetch nonunique_seeds — on a miss in unique, we'll
        // need this immediately and the prefetch is essentially free
        // if we end up not needing it (the line just gets evicted).
        self.nonunique_seeds.prefetch_key(kmer);
    }

    pub(crate) fn load_chrom_info<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<ChromInfo>> {
        let content = std::fs::read_to_string(path)?;
        serde_json::from_str(&content)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    pub(crate) fn save_chrom_info<P: AsRef<Path>>(
        chrom_info: &[ChromInfo],
        path: P,
    ) -> std::io::Result<()> {
        let content = serde_json::to_string(chrom_info)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        std::fs::write(path, content)
    }
}

impl<const K: usize, const S: usize> super::Index for AsymmetricIndex<K, S> {
    type LocusType = Locus;

    fn load<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        Self::load(path)
    }

    fn save<P: AsRef<Path>>(&self, path: P, portable: bool) -> std::io::Result<()> {
        if portable {
            self.save_portable(path)
        } else {
            self.save_native(path)
        }
    }

    /// Get all chromosome info.
    fn all_chrom_info(&self) -> &[ChromInfo] {
        &self.chrom_info
    }

    fn chrom_info(&self, chrom_idx: usize) -> &crate::reference::ChromInfo {
        &self.chrom_info[chrom_idx]
    }

    fn find_seeds<F>(&self, seq: &[u8], callback: F)
    where
        F: FnMut(usize, u64, usize, &[Self::LocusType]) {
        let mut kmer_batch = self.query_buffers.acquire();
        kmer_batch.clear();
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            seq,
            [(); S],
            |pos, kmer| {
                kmer_batch.push((pos, kmer.0));
            },
        );
        self.lookup_batch(&kmer_batch, callback);
    }
}

impl<const K: usize, const S: usize> SyncmerIndex<K, S> for AsymmetricIndex<K, S> {
    fn with<F: FnMut(usize, &[Locus])>(&self, kmer: &Kmer<K>, mut f: F) {
        if let Some(loc) = self.unique_seeds.get(kmer.0) {
            let buf = [Locus::from(loc)];
            f(1, &buf);
        } else if let Some(loci) = self.nonunique_seeds.get(kmer.0) {
            let loci =
                unsafe { std::slice::from_raw_parts(loci.as_ptr() as *const Locus, loci.len()) };

            f(loci.len(), loci);
        }
    }

    fn lookup_batch<F>(&self, batch: &[(usize, u64)], mut callback: F)
    where
        F: FnMut(usize, u64, usize, &[Self::LocusType]),
    {
        const PIPE: usize = 16;
        let n = batch.len();

        // Issue initial prefetches for the first PIPE entries
        for i in 0..PIPE.min(n) {
            self.prefetch_kmer(batch[i].1);
        }

        for i in 0..n {
            // Launch prefetch for entry i+PIPE (if it exists)
            if i + PIPE < n {
                self.prefetch_kmer(batch[i + PIPE].1);
            }

            // Harvest: the data for batch[i] was prefetched PIPE iterations ago
            let (read_pos, kmer_val) = batch[i];

            if let Some(loc) = self.unique_seeds.get(kmer_val) {
                let buf = [Locus::from(loc)];
                callback(read_pos, kmer_val, 1, &buf);
            } else if let Some(loci) = self.nonunique_seeds.get(kmer_val) {
                // SAFETY: Locus is repr(transparent) over u64, so &[u64] and
                // &[Locus] have identical size and alignment. We are only
                // changing the type the compiler uses to interpret existing bytes.
                let loci = unsafe {
                    std::slice::from_raw_parts(loci.as_ptr() as *const Locus, loci.len())
                };
                callback(read_pos, kmer_val, loci.len(), loci);
            }
        }
    }
}

/// Mutable builder for constructing a K-mer seed index.
///
/// Use this to build an index from reference sequences, then call `build()`
/// to create an immutable `Index` for fast lookups.
pub struct AsymmetricIndexBuilder<const K: usize, const S: usize> {}

impl<const K: usize, const S: usize> AsymmetricIndexBuilder<K, S> {
    const MAX_SPACING: usize = K - S + 1;
}

impl<const K: usize, const S: usize> super::IndexBuilder<K, S> for AsymmetricIndexBuilder<K, S> {
    type IndexType = AsymmetricIndex<K, S>;

    fn build(reference: &InMemoryReference) -> Self::IndexType {
        // Step 1: Count frequencies of all syncmers in the reference
        //
        let mut frequencies: Table<u64, u32> = Table::new();
        for chrom_idx in 0..reference.num_chroms() {
            let chrom: &str = &reference.chrom_info(chrom_idx).name;

            log::info!(
                "Frequency counting chrom {} \"{}\" for index construction",
                chrom_idx,
                chrom
            );

            let seq = reference.sequence(chrom_idx);
            for (_pos, sel) in Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(seq, [(); S]) {
                match sel {
                    Selection::Left(lhs) => {
                        if let Some(count) = frequencies.get_mut(&lhs.0) {
                            *count += 1;
                        } else {
                            frequencies.insert(lhs.0, 1);
                        }
                    }
                    Selection::Both(lhs, rhs) => {
                        if let Some(count) = frequencies.get_mut(&lhs.0) {
                            *count += 1;
                        } else {
                            frequencies.insert(lhs.0, 1);
                        }
                        if let Some(count) = frequencies.get_mut(&rhs.0) {
                            *count += 1;
                        } else {
                            frequencies.insert(rhs.0, 1);
                        }
                    }
                    Selection::Right(rhs) => {
                        if let Some(count) = frequencies.get_mut(&rhs.0) {
                            *count += 1;
                        } else {
                            frequencies.insert(rhs.0, 1);
                        }
                    }
                }
            }
        }

        // Step 2: Construct the index, 1 chromosome at a time.
        //

        let chrom_info: Vec<ChromInfo> = (0..reference.num_chroms())
            .map(|i| reference.chrom_info(i).clone())
            .collect();
        let mut unique_seeds: Table<u64, u64> = Table::new();
        let mut nonunique_seeds: HashMap<u64, Vec<u64>> = HashMap::new();

        let mut items: Vec<(u64, Locus)> = Vec::new();
        let mut permutation: Vec<u32> = Vec::new();
        let mut keys: Vec<u64> = Vec::new();
        let mut alive: Vec<bool> = Vec::new();

        for chrom_idx in 0..reference.num_chroms() {
            let chrom: &str = &reference.chrom_info(chrom_idx).name;

            items.clear();
            permutation.clear();
            keys.clear();
            alive.clear();

            log::info!(
                "Phase 2 syncmer extraction from chrom {} \"{}\" for index construction",
                chrom_idx,
                chrom
            );

            let seq = reference.sequence(chrom_idx);
            let mut i = 0;
            for (pos, sel) in Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(seq, [(); S]) {
                match sel {
                    Selection::Left(lhs) => {
                        let locus = Locus::pack(chrom_idx, pos, super::Strand::Forward);
                        items.push((lhs.0, locus));
                        permutation.push(i);
                        keys.push(make_key(lhs.0, locus, &frequencies));
                        i += 1;
                    }
                    Selection::Both(lhs, rhs) => {
                        let locus = Locus::pack(chrom_idx, pos, super::Strand::Forward);
                        items.push((lhs.0, locus));
                        permutation.push(i);
                        keys.push(make_key(lhs.0, locus, &frequencies));
                        i += 1;
                        let locus = Locus::pack(chrom_idx, pos, super::Strand::Reverse);
                        items.push((rhs.0, locus));
                        permutation.push(i);
                        keys.push(make_key(rhs.0, locus, &frequencies));
                        i += 1;
                    }
                    Selection::Right(rhs) => {
                        let locus = Locus::pack(chrom_idx, pos, super::Strand::Reverse);
                        items.push((rhs.0, locus));
                        permutation.push(i);
                        keys.push(make_key(rhs.0, locus, &frequencies));
                        i += 1;
                    }
                }
            }

            log::info!(
                "Sorting chrom {} \"{}\" by frequency for index construction",
                chrom_idx,
                chrom
            );

            permutation.sort_unstable_by_key(|&i| Reverse(keys[i as usize]));

            log::info!(
                "Filtering chrom {} \"{}\" for index construction",
                chrom_idx,
                chrom
            );

            let n = items.len();
            alive.resize(n, true);

            for &idx in &permutation {
                if !alive[idx as usize] {
                    continue;
                }
                let (_chrom_id, this_pos, _strand) = items[idx as usize].1.unpack();
                let left = (0..idx as usize).rev().find(|&j| alive[j]);
                if let Some(l) = left {
                    let (_chrom_id, left_pos, _strand) = items[l as usize].1.unpack();
                    if left_pos == this_pos {
                        alive[idx as usize] = false;
                        continue;
                    }
                }
                let right = (idx as usize + 1..n).find(|&j| alive[j]);
                if let Some(r) = right {
                    let (_chrom_id, right_pos, _strand) = items[r as usize].1.unpack();
                    if right_pos == this_pos {
                        alive[idx as usize] = false;
                        continue;
                    }
                }
                let keep = match (left, right) {
                    (None, None) | (None, Some(_)) | (Some(_), None) => true,
                    (Some(l), Some(r)) => {
                        let (_chrom_id, left_pos, _strand) = items[l as usize].1.unpack();
                        let (_chrom_id, right_pos, _strand) = items[r as usize].1.unpack();
                        right_pos - left_pos > Self::MAX_SPACING
                    }
                };
                alive[idx as usize] = keep;
            }

            log::info!(
                "Finalizing chrom {} \"{}\" for index construction",
                chrom_idx,
                chrom
            );

            for (item, keep) in items.iter().zip(alive.iter()) {
                if !*keep {
                    continue;
                }

                let kmer = item.0;
                let locus = item.1;
                let count = frequencies.get(&kmer).copied().unwrap_or(0);
                if count == 1 {
                    unique_seeds.insert(kmer, locus.into());
                } else {
                    nonunique_seeds
                        .entry(kmer)
                        .or_insert_with(Vec::new)
                        .push(locus.into());
                }
            }
        }

        let unique_seeds = FrozenTable::from_table(unique_seeds);
        let nonunique_seeds = FrozenBigTable::from_hashmap(nonunique_seeds);

        AsymmetricIndex {
            chrom_info,
            unique_seeds,
            nonunique_seeds,
            query_buffers: Pool::new()
        }
    }
}

fn hash(x: u64) -> u32 {
    let x = Splitmix64Hasher::hash64(x);
    (x >> 32) as u32 ^ (x & 0xFFFFFFFF) as u32
}

fn make_key(kmer: u64, locus: Locus, frequencies: &Table<u64, u32>) -> u64 {
    let count = frequencies.get(&kmer).copied().unwrap_or(0);
    let h = hash(kmer) ^ hash(locus.into());
    (count as u64) << 32 | h as u64
}
