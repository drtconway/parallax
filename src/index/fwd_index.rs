use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use crossbeam::channel;
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
use crate::utils::pool::Pool;
use crate::utils::table::Table;

use super::BedRegions;

const INDEX_VERSION: u32 = 1;
const INDEX_TYPE: &str = "fwd-syncmer";

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
        if cfg!(target_endian = "little") { "little" } else { "big" }
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
                format!("Index version mismatch: index has version {}, expected {}", self.version, INDEX_VERSION),
            ));
        }
        if self.index_type != INDEX_TYPE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Index type mismatch: index has type \"{}\", expected \"{}\"", self.index_type, INDEX_TYPE),
            ));
        }
        if self.k != K {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Index k-mer size mismatch: index has k={}, binary expects k={}", self.k, K),
            ));
        }
        if self.s != S {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Index syncmer size mismatch: index has s={}, binary expects s={}", self.s, S),
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

// repr(transparent) guarantees FwdLocus has the same memory layout as u64.
// This allows lookup_batch to reinterpret the raw &[u64] locus storage as
// &[FwdLocus] without copying. Removing this attribute would make that cast
// undefined behaviour.
#[repr(transparent)]
pub struct FwdLocus(u64);

impl From<u64> for FwdLocus {
    fn from(value: u64) -> Self {
        FwdLocus(value)
    }
}

impl Into<u64> for FwdLocus {
    fn into(self) -> u64 {
        self.0
    }
}

impl PackedLocus for FwdLocus {
    fn pack(chrom: usize, pos: usize, strand: super::Strand) -> Self {
        debug_assert!(
            strand == super::Strand::Forward,
            "FwdLocus only supports forward strand"
        );
        let chrom: u64 = chrom as u64;
        let pos: u64 = pos as u64;
        FwdLocus(chrom << 32 | pos)
    }

    fn unpack(&self) -> (usize, usize, super::Strand) {
        let chrom = (self.0 >> 32) as usize;
        let pos = (self.0 & 0xFFFFFFFF) as usize;
        (chrom, pos, super::Strand::Forward)
    }
}

/// Immutable K-mer seed index backed by Arrow arrays.
///
/// This is the production index structure, optimized for fast lookups and
/// supporting save/load to Parquet files for persistence.
pub struct FwdIndex<const K: usize, const S: usize> {
    chrom_info: Vec<ChromInfo>,
    unique_seeds: FrozenTable,
    nonunique_seeds: FrozenBigTable,
    query_buffers: Pool<Vec<(usize, u64)>>
}

impl<const K: usize, const S: usize> FwdIndex<K, S> {
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
        Ok(FwdIndex { chrom_info, unique_seeds, nonunique_seeds, query_buffers: Pool::new() })
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
        Ok(FwdIndex { chrom_info, unique_seeds, nonunique_seeds, query_buffers: Pool::new() })
    }

    /// Save the index in portable (Feather/Arrow IPC) format.
    pub fn save_portable<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;
        IndexMetadata::new_portable(K, S).save(dir.join("metadata.json"))?;
        Self::save_chrom_info(&self.chrom_info, dir.join("chrom_info.json"))?;
        self.unique_seeds.save_to_feather_directory(dir.join("unique_seeds"))?;
        self.nonunique_seeds.save_to_feather_directory(dir.join("nonunique_seeds"))?;
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
        self.unique_seeds.save_to_directory(dir.join("unique_seeds"))?;
        self.nonunique_seeds.save_to_directory(dir.join("nonunique_seeds"))?;
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

    /// Get the chromosome name for a given index.
    #[allow(dead_code)]
    pub fn chrom_name(&self, chrom_idx: usize) -> &str {
        &self.chrom_info[chrom_idx].name
    }

    /// Get the full ChromInfo for a given index.
    #[allow(dead_code)]
    pub fn chrom_info(&self, chrom_idx: usize) -> &ChromInfo {
        &self.chrom_info[chrom_idx]
    }

    /// Get the number of unique seeds.
    #[allow(dead_code)]
    pub fn unique_count(&self) -> usize {
        self.unique_seeds.len()
    }

    /// Get the number of nonunique seeds.
    #[allow(dead_code)]
    pub fn nonunique_count(&self) -> usize {
        self.nonunique_seeds.len()
    }

    /// Validate that this index is compatible with the given reference.
    ///
    /// Checks that chromosome names and lengths match between the index
    /// and the reference. Returns an error describing the first mismatch
    /// found, or `Ok(())` if they are compatible.
    pub fn validate_reference(&self, reference: &InMemoryReference) -> std::io::Result<()> {
        let idx_n = self.chrom_info.len();
        let ref_n = reference.num_chroms();

        if idx_n != ref_n {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Index has {} chromosome(s) but reference has {}. \
                     Was the index built with a different reference?",
                    idx_n, ref_n
                ),
            ));
        }

        for i in 0..idx_n {
            let idx_name = &self.chrom_info[i].name;
            let ref_name = reference.chrom_name(i);

            if idx_name != ref_name {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Chromosome {} name mismatch: index has \"{}\" but reference has \"{}\". \
                         Was the index built with a different reference?",
                        i, idx_name, ref_name
                    ),
                ));
            }

            let idx_len = self.chrom_info[i].length;
            let ref_len = reference.chrom_length(i);

            // Skip length check if index predates the length field (stored as 0)
            if idx_len != 0 && idx_len != ref_len {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Chromosome \"{}\" length mismatch: index has {} but reference has {}. \
                         Was the index built with a different reference?",
                        idx_name, idx_len, ref_len
                    ),
                ));
            }
        }

        Ok(())
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

impl<const K: usize, const S: usize> super::Index for FwdIndex<K, S> {
    type LocusType = FwdLocus;

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

impl<const K: usize, const S: usize> super::SyncmerIndex<K, S> for FwdIndex<K, S> {
    fn with<F: FnMut(usize, &[FwdLocus])>(&self, kmer: &Kmer<K>, mut f: F) {
        if let Some(loc) = self.unique_seeds.get(kmer.0) {
            let buf = [FwdLocus::from(loc)];
            f(1, &buf);
        } else if let Some(loci) = self.nonunique_seeds.get(kmer.0) {
            let loci =
                unsafe { std::slice::from_raw_parts(loci.as_ptr() as *const FwdLocus, loci.len()) };

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
                let buf = [FwdLocus::from(loc)];
                callback(read_pos, kmer_val, 1, &buf);
            } else if let Some(loci) = self.nonunique_seeds.get(kmer_val) {
                // SAFETY: FwdLocus is repr(transparent) over u64, so &[u64] and
                // &[FwdLocus] have identical size and alignment. We are only
                // changing the type the compiler uses to interpret existing bytes.
                let loci = unsafe {
                    std::slice::from_raw_parts(loci.as_ptr() as *const FwdLocus, loci.len())
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
pub struct FwdIndexBuilder<const K: usize, const S: usize> {
    pub(crate) chrom_info: Vec<ChromInfo>,
    pub(crate) unique_seeds: Table<u64, u64>,
    pub(crate) nonunique_seeds: HashMap<u64, Vec<u64>>,
}

impl<const K: usize, const S: usize> FwdIndexBuilder<K, S> {
    pub fn new() -> Self {
        FwdIndexBuilder {
            chrom_info: Vec::new(),
            unique_seeds: Table::new(),
            nonunique_seeds: HashMap::new(),
        }
    }

    /// Create a new builder pre-initialized with chromosome info and lengths.
    pub fn new_with_chrom_info(chrom_info: &[(ChromInfo, usize)]) -> Self {
        let chrom_info_vec: Vec<ChromInfo> = chrom_info
            .iter()
            .map(|(info, len)| {
                let mut info = info.clone();
                info.length = *len as u64;
                info
            })
            .collect();

        FwdIndexBuilder {
            chrom_info: chrom_info_vec,
            unique_seeds: Table::new(),
            nonunique_seeds: HashMap::new(),
        }
    }

    /// Add a chromosome sequence by index.
    /// The builder must have been created with `new_with_chrom_info` and
    /// `chrom_idx` must be valid for the chromosome list.
    #[allow(dead_code)]
    pub fn add_chrom<Seq: AsRef<[u8]>>(&mut self, chrom_idx: usize, seq: Seq) {
        self.add_chrom_regions(chrom_idx, seq, None);
    }

    /// Add a chromosome sequence by index with optional region restriction.
    /// If regions are provided (sorted, non-overlapping intervals), only k-mers
    /// starting within those regions will be indexed. The position in the index
    /// is the absolute position in the chromosome (region_start + relative_pos).
    pub fn add_chrom_regions<Seq: AsRef<[u8]>>(
        &mut self,
        chrom_idx: usize,
        seq: Seq,
        regions: Option<&[(usize, usize)]>,
    ) {
        assert!(chrom_idx < self.chrom_info.len(), "chrom_idx out of bounds");

        let seq = seq.as_ref();
        let n = seq.len();
        let mut m = 0usize;

        match regions {
            None => {
                // No regions specified - index the entire chromosome
                for (pos, sel) in Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(seq, [(); S]) {
                    match sel {
                        Selection::Left(kmer) | Selection::Both(kmer, _) => {
                            let loc = FwdLocus::pack(chrom_idx, pos, super::Strand::Forward).into();
                            m += 1;
                            self.insert_kmer(kmer.0, loc);
                        }
                        _ => {}
                    }
                }
                self.chrom_info[chrom_idx].syncmer_count = m as u64;
                let r = (m as f64) / (n as f64);
                log::debug!(
                    "Indexed chrom {} \"{}\" (length {}, {} seeds; {:.4})",
                    chrom_idx,
                    self.chrom_info[chrom_idx].name,
                    n,
                    m,
                    r
                );
            }
            Some(intervals) => {
                // Index only within the specified regions
                let mut region_bases = 0usize;
                for &(start, end) in intervals {
                    let start = start.min(n);
                    let end = end.min(n);
                    if start >= end {
                        continue;
                    }
                    region_bases += end - start;

                    let region_seq = &seq[start..end];
                    for (rel_pos, sel) in
                        Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(region_seq, [(); S])
                    {
                        match sel {
                            Selection::Left(kmer) | Selection::Both(kmer, _) => {
                                // Absolute position = region start + relative position
                                let abs_pos = start + rel_pos;
                                let loc =
                                    FwdLocus::pack(chrom_idx, abs_pos, super::Strand::Forward)
                                        .into();
                                m += 1;
                                self.insert_kmer(kmer.0, loc);
                            }
                            _ => {}
                        }
                    }
                }
                self.chrom_info[chrom_idx].syncmer_count = m as u64;
                let r = if region_bases > 0 {
                    (m as f64) / (region_bases as f64)
                } else {
                    0.0
                };
                log::debug!(
                    "Indexed chrom {} \"{}\" ({} bp in {} regions, {} seeds; {:.4})",
                    chrom_idx,
                    self.chrom_info[chrom_idx].name,
                    region_bases,
                    intervals.len(),
                    m,
                    r
                );
            }
        }
    }

    /// Insert a k-mer into the index, handling unique vs nonunique.
    #[inline]
    pub(crate) fn insert_kmer(&mut self, kmer: u64, loc: u64) {
        if let Some(loc0) = self.unique_seeds.remove(&kmer) {
            self.nonunique_seeds.insert(kmer, vec![loc0, loc]);
        } else if let Some(locs) = self.nonunique_seeds.get_mut(&kmer) {
            locs.push(loc);
        } else {
            self.unique_seeds.insert(kmer, loc);
        }
    }

    pub fn add<Seq: AsRef<[u8]>>(&mut self, mut chrom_info: ChromInfo, seq: Seq) -> usize {
        let chrom_id = self.chrom_info.len();
        let n = seq.as_ref().len();
        chrom_info.length = n as u64;
        self.chrom_info.push(chrom_info);

        let mut m = 0usize;
        for (pos, sel) in Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(seq.as_ref(), [(); S]) {
            match sel {
                Selection::Left(kmer) | Selection::Both(kmer, _) => {
                    let loc = FwdLocus::pack(chrom_id, pos, super::Strand::Forward).into();
                    m += 1;
                    let x = kmer.0;
                    if let Some(loc0) = self.unique_seeds.remove(&x) {
                        self.nonunique_seeds.insert(x, vec![loc0, loc]);
                    } else if let Some(locs) = self.nonunique_seeds.get_mut(&x) {
                        locs.push(loc);
                    } else {
                        self.unique_seeds.insert(x, loc);
                    }
                }
                _ => { /* do nothing for right-only selections */ }
            }
        }

        self.chrom_info[chrom_id].syncmer_count = m as u64;
        let r = (m as f64) / (n as f64);
        log::debug!(
            "Indexed chrom {} \"{}\" (length {}, {} seeds; {:.4})",
            chrom_id,
            self.chrom_info[chrom_id].name,
            n,
            m,
            r
        );
        chrom_id
    }

    /// Build an immutable Index from this builder.
    pub fn build(self) -> FwdIndex<K, S> {
        let now = std::time::Instant::now();

        let unique_seeds = FrozenTable::from_table(self.unique_seeds);
        let nonunique_seeds = FrozenBigTable::from_hashmap(self.nonunique_seeds);

        log::info!(
            "Built frozen index: {} unique seeds, {} nonunique seeds in {:.2}s",
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        FwdIndex {
            chrom_info: self.chrom_info,
            unique_seeds,
            nonunique_seeds,
            query_buffers: Pool::new()
        }
    }

    /// Merge another builder into this one.
    /// Both builders must have been created with the same chromosome list
    /// (via `new_with_chrom_info`), so no locus recoding is needed.
    pub fn merge(&mut self, other: FwdIndexBuilder<K, S>) {
        // Verify chromosome lists match
        debug_assert_eq!(
            self.chrom_info.len(),
            other.chrom_info.len(),
            "Builders must have same chromosome list"
        );

        let now = std::time::Instant::now();

        // Merge unique seeds from other (no recoding needed)
        for (kmer, loc) in other.unique_seeds.iter() {
            if let Some(existing_loc) = self.unique_seeds.remove(kmer) {
                // Was unique in self, now becomes nonunique
                self.nonunique_seeds.insert(*kmer, vec![existing_loc, *loc]);
            } else if let Some(locs) = self.nonunique_seeds.get_mut(kmer) {
                // Already nonunique in self
                locs.push(*loc);
            } else {
                // New unique seed
                self.unique_seeds.insert(*kmer, *loc);
            }
        }

        // Merge nonunique seeds from other (no recoding needed)
        for (kmer, other_locs) in other.nonunique_seeds {
            if let Some(existing_loc) = self.unique_seeds.remove(&kmer) {
                // Was unique in self, merge with other's nonunique
                let mut combined = vec![existing_loc];
                combined.extend(other_locs);
                self.nonunique_seeds.insert(kmer, combined);
            } else if let Some(locs) = self.nonunique_seeds.get_mut(&kmer) {
                // Already nonunique in self
                locs.extend(other_locs);
            } else {
                // New nonunique seed
                self.nonunique_seeds.insert(kmer, other_locs);
            }
        }

        let t = now.elapsed().as_secs_f64();

        if t > 2.0 {
            log::info!(
                "Merged builders, now {} unique and {} nonunique seeds in {:.2}s",
                self.unique_seeds.len(),
                self.nonunique_seeds.len(),
                t
            );
        }
    }

    /// Build an index from an InMemoryReference using multiple threads.
    ///
    /// Each chromosome is indexed in parallel, then the per-chromosome builders
    /// are merged together. All builders share the same chromosome list, so
    /// merging is fast (no locus recoding required).
    ///
    /// If `bed_regions` is provided, only k-mers starting within those regions
    /// will be indexed. Regions are specified as a map from chromosome name to
    /// sorted intervals.
    pub fn build_parallel(
        reference: &InMemoryReference,
        bed_regions: Option<&BedRegions>,
        num_threads: usize,
    ) -> FwdIndex<K, S> {
        let num_chroms = reference.num_chroms();
        if num_chroms == 0 {
            return FwdIndexBuilder::<K, S>::new().build();
        }

        let now = std::time::Instant::now();

        let num_threads = num_threads.max(1);

        // Build chromosome info for pre-initialization
        let chrom_info: Vec<(ChromInfo, usize)> = (0..num_chroms)
            .map(|i| {
                (
                    reference.chrom_info(i).clone(),
                    reference.chrom_length(i) as usize,
                )
            })
            .collect();
        let chrom_info = Arc::new(chrom_info);

        if bed_regions.is_some() {
            log::info!(
                "Building index for {} chromosomes using {} threads (with BED region masking)",
                num_chroms,
                num_threads
            );
        } else {
            log::info!(
                "Building index for {} chromosomes using {} threads",
                num_chroms,
                num_threads
            );
        }

        // Build a map from chrom index to intervals if BED regions are provided
        // Chromosomes not in the BED file get an empty slice (index nothing)
        let chrom_intervals: Option<Arc<Vec<Vec<(usize, usize)>>>> = bed_regions.map(|regions| {
            let intervals: Vec<Vec<(usize, usize)>> = (0..num_chroms)
                .map(|i| {
                    let chrom_name = &reference.chrom_info(i).name;
                    regions.get(chrom_name).cloned().unwrap_or_default()
                })
                .collect();
            Arc::new(intervals)
        });

        // Create work items
        let reference = Arc::new(reference.clone());
        let (work_tx, work_rx) = channel::bounded::<usize>(num_chroms);
        let (result_tx, result_rx) = channel::bounded::<(usize, FwdIndexBuilder<K, S>)>(num_chroms);

        // Spawn worker threads
        let mut handles = Vec::with_capacity(num_threads);
        for _ in 0..num_threads {
            let work_rx = work_rx.clone();
            let result_tx = result_tx.clone();
            let reference = Arc::clone(&reference);
            let chrom_info = Arc::clone(&chrom_info);
            let chrom_intervals = chrom_intervals.clone();

            let handle = std::thread::spawn(move || {
                for chrom_idx in work_rx {
                    let seq = reference.sequence(chrom_idx);

                    // Create builder with shared chromosome info
                    let mut builder = FwdIndexBuilder::<K, S>::new_with_chrom_info(&chrom_info);

                    // Get intervals for this chromosome if BED regions were provided
                    let regions = chrom_intervals
                        .as_ref()
                        .map(|intervals| intervals[chrom_idx].as_slice());

                    builder.add_chrom_regions(chrom_idx, seq, regions);

                    result_tx
                        .send((chrom_idx, builder))
                        .expect("result channel closed");
                }
            });
            handles.push(handle);
        }

        // Drop our copy of result_tx so channel closes when workers finish
        drop(result_tx);

        // Send work items
        for chrom_idx in 0..num_chroms {
            work_tx.send(chrom_idx).expect("work channel closed");
        }
        drop(work_tx);

        // Collect results and merge
        let mut combined: Option<FwdIndexBuilder<K, S>> = None;
        let mut next_chrom = 0;
        let mut buffered_builders: HashMap<usize, FwdIndexBuilder<K, S>> = HashMap::new();
        for (chrom_idx, builder) in result_rx {
            // Buffer out-of-order builders
            if chrom_idx != next_chrom {
                buffered_builders.insert(chrom_idx, builder);
                continue;
            }
            match combined {
                None => {
                    combined = Some(builder);
                    next_chrom += 1;
                }
                Some(ref mut c) => {
                    c.merge(builder);
                    next_chrom += 1;
                }
            }
            while let Some(builder) = buffered_builders.remove(&next_chrom) {
                if let Some(ref mut c) = combined {
                    c.merge(builder);
                    next_chrom += 1;
                }
            }
        }
        assert_eq!(
            next_chrom, num_chroms,
            "Did not receive all chromosome builders"
        );
        assert_eq!(
            buffered_builders.len(),
            0,
            "Buffered builders not empty after merge"
        );
        // Wait for all workers to finish
        for handle in handles {
            handle.join().expect("worker thread panicked");
        }

        let combined =
            combined.unwrap_or_else(|| FwdIndexBuilder::<K, S>::new_with_chrom_info(&chrom_info));

        let t = now.elapsed().as_secs_f64();

        if t > 5.0 {
            log::info!(
                "Index complete: {} unique seeds, {} nonunique seeds in {:.2}s",
                combined.unique_seeds.len(),
                combined.nonunique_seeds.len(),
                t
            );
        }

        combined.build()
    }
}

impl<const K: usize, const S: usize, R: std::io::BufRead> TryFrom<noodles::fasta::io::Reader<R>>
    for FwdIndexBuilder<K, S>
{
    type Error = crate::error::ParallaxError;

    fn try_from(mut reader: noodles::fasta::io::Reader<R>) -> Result<Self, Self::Error> {
        let mut builder = FwdIndexBuilder::new();

        for record in reader.records() {
            let record = record?;
            let name = String::from_utf8(record.name().to_vec())
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let description = record
                .description()
                .map(|d| String::from_utf8(d.to_vec()).unwrap_or_default())
                .unwrap_or_default();
            let chrom_info = ChromInfo::from_header(&name, &description);
            let seq = record.sequence().as_ref().to_owned();
            builder.add(chrom_info, seq);
        }

        if false {
            let mut hist = HashMap::new();
            hist.insert(1usize, builder.unique_seeds.len());
            for locs in builder.nonunique_seeds.values() {
                let count = locs.len();
                *hist.entry(count).or_insert(0usize) += 1;
            }
            let mut hist_vec: Vec<(usize, usize)> = hist.into_iter().collect();
            hist_vec.sort_by_key(|(count, _freq)| *count);
            log::info!("Seed occurrence histogram:");
            for (count, freq) in hist_vec {
                log::info!("  {}\t{}", count, freq);
            }
        }

        Ok(builder)
    }
}
