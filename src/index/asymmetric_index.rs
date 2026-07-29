use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::index::Index;
use crate::index::IndexHit;
use crate::index::LoadableIndex;
use crate::index::PackedLocus;
use crate::index::SyncmerIndex;
use crate::kmers::Kmer;
use crate::reference::ChromInfo;
use crate::reference::InMemoryReference;
use crate::utils::frozen_big_table::FrozenBigTable;
use crate::utils::frozen_table::FrozenTable;
use crate::utils::hasher::FnvHasher;
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

impl Locus {
    pub fn unpack_from_u64(locus: u64) -> (usize, usize) {
        Locus::from(locus).unpack()
    }
}

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
    fn pack(chrom: usize, pos: usize) -> Self {
        let chrom: u64 = chrom as u64;
        let pos: u64 = pos as u64;
        Locus(chrom << 32 | pos)
    }

    fn unpack(&self) -> (usize, usize) {
        let chrom = (self.0 >> 32) as usize;
        let pos = (self.0 & 0xFFFFFFFF) as usize;
        (chrom, pos)
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
    query_buffers: Pool<Vec<(usize, u64)>>,
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
            query_buffers: Pool::new(),
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
            query_buffers: Pool::new(),
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
    fn save(&self, path: &Path, portable: bool) -> std::io::Result<()> {
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

    fn k(&self) -> usize {
        K
    }

    fn lookup_kmer(&self, kmer: u64) -> Option<IndexHit<'_>> {
        if let Some(loci) = self.unique_seeds.get_as_slice(kmer) {
            return Some(IndexHit {
                query_pos: 0,
                seed_kmer: kmer,
                loci,
                k: K,
            });
        }
        if let Some(loci) = self.nonunique_seeds.get(kmer) {
            return Some(IndexHit {
                query_pos: 0,
                seed_kmer: kmer,
                loci,
                k: K,
            });
        }
        None
    }

    fn unpack_locus(&self, locus: u64) -> (usize, usize) {
        Locus::unpack_from_u64(locus)
    }

    /// Unpack a slice of loci
    fn unpack_loci(&self, packed: &[u64], unpacked: &mut Vec<(usize, usize)>) {
        unpacked.clear();
        unpacked.reserve(packed.len());
        for &value in packed {
            let locus = self.unpack_locus(value);
            unpacked.push(locus);
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = IndexHit<'_>> + '_> {
        let unique = self.unique_seeds.iter().map(|(kmer, slot)| IndexHit {
            query_pos: 0,
            seed_kmer: kmer,
            loci: self.unique_seeds.value_as_slice(slot),
            k: K,
        });
        let nonunique = self.nonunique_seeds.iter().map(|(kmer, slot)| IndexHit {
            query_pos: 0,
            seed_kmer: kmer,
            loci: self.nonunique_seeds.loci_as_slice(slot),
            k: K,
        });
        Box::new(unique.chain(nonunique))
    }

    fn find_seeds(&self, seq: &[u8], callback: &mut dyn FnMut(IndexHit<'_>)) {
        let mut query_kmers = self.query_buffers.acquire();
        query_kmers.clear();
        for (pos, fwd, _rev) in
            Kmer::<K>::agnostic_open_syncmer_iter::<S, FnvHasher>(seq.as_ref(), [(); S])
        {
            query_kmers.push((pos, fwd.0));
        }
        self.lookup_batch(&query_kmers, |hit| {
            callback(hit);
        });
    }
}

impl<const K: usize, const S: usize> LoadableIndex for AsymmetricIndex<K, S> {
    fn load(path: &Path) -> std::io::Result<Self> {
        Self::load(path)
    }
}
impl<const K: usize, const S: usize> SyncmerIndex<K, S> for AsymmetricIndex<K, S> {
    fn with<F: FnMut(IndexHit<'_>)>(&self, kmer: &Kmer<K>, mut f: F) {
        let query_pos = 0;
        let seed_kmer = kmer.0;
        if let Some(locus) = self.unique_seeds.get(seed_kmer) {
            let buf = [locus];
            f(IndexHit {
                query_pos,
                seed_kmer,
                loci: &buf,
                k: K,
            });
        } else if let Some(loci) = self.nonunique_seeds.get(seed_kmer) {
            f(IndexHit {
                query_pos,
                seed_kmer,
                loci,
                k: K,
            });
        }
    }

    fn lookup_batch<F>(&self, batch: &[(usize, u64)], mut callback: F)
    where
        F: FnMut(IndexHit<'_>),
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
            let (query_pos, seed_kmer) = batch[i];

            if let Some(locus) = self.unique_seeds.get(seed_kmer) {
                let buf = [locus];
                callback(IndexHit {
                    query_pos,
                    seed_kmer,
                    loci: &buf,
                    k: K,
                });
            } else if let Some(loci) = self.nonunique_seeds.get(seed_kmer) {
                callback(IndexHit {
                    query_pos,
                    seed_kmer,
                    loci,
                    k: K,
                });
            }
        }
    }
}

pub fn load_index(path: &Path) -> std::io::Result<Arc<dyn Index>> {
    let meta = IndexMetadata::load(path.join("metadata.json"))?;
    if meta.index_type != INDEX_TYPE {
        return Err(std::io::Error::other(format!(
            "unexpected index type '{}'",
            meta.index_type
        )));
    }
    if meta.version != INDEX_VERSION {
        return Err(std::io::Error::other(format!(
            "unexpected index version '{}' - rebuild your index",
            meta.version
        )));
    }
    let k = meta.k;
    let s = meta.s;
    match k {
        15 => load_index_inner::<15>(path, s),
        16 => load_index_inner::<16>(path, s),
        17 => load_index_inner::<17>(path, s),
        18 => load_index_inner::<18>(path, s),
        19 => load_index_inner::<19>(path, s),
        20 => load_index_inner::<20>(path, s),
        21 => load_index_inner::<21>(path, s),
        22 => load_index_inner::<22>(path, s),
        23 => load_index_inner::<23>(path, s),
        24 => load_index_inner::<24>(path, s),
        25 => load_index_inner::<25>(path, s),
        _ => Err(std::io::Error::other("unsupported value of K")),
    }
}

fn load_index_inner<const K: usize>(path: &Path, s: usize) -> std::io::Result<Arc<dyn Index>> {
    if s > K {
        return Err(std::io::Error::other(format!(
            "K = {}, so S must be between 10 and {}",
            K, K
        )));
    }
    match s {
        10 => load_index_inner_2::<K, 10>(path),
        11 => load_index_inner_2::<K, 11>(path),
        12 => load_index_inner_2::<K, 12>(path),
        13 => load_index_inner_2::<K, 13>(path),
        14 => load_index_inner_2::<K, 14>(path),
        15 => load_index_inner_2::<K, 15>(path),
        16 => load_index_inner_2::<K, 16>(path),
        17 => load_index_inner_2::<K, 17>(path),
        18 => load_index_inner_2::<K, 18>(path),
        19 => load_index_inner_2::<K, 19>(path),
        20 => load_index_inner_2::<K, 20>(path),
        _ => Err(std::io::Error::other("unsupported value of S")),
    }
}

fn load_index_inner_2<const K: usize, const S: usize>(
    path: &Path,
) -> std::io::Result<Arc<dyn Index>> {
    let index = AsymmetricIndex::<K, S>::load(path)?;
    Ok(Arc::new(index))
}

/// Mutable builder for constructing a K-mer seed index.
///
/// Use this to build an index from reference sequences, then call `build()`
/// to create an immutable `Index` for fast lookups.
pub struct AsymmetricIndexBuilder<'a, const K: usize, const S: usize> {
    reference: &'a InMemoryReference
}

impl<'a, const K: usize, const S: usize> super::IndexBuilder<'a, K, S> for AsymmetricIndexBuilder<'a, K, S> {
    type IndexType = AsymmetricIndex<K, S>;

    fn make(reference: &'a InMemoryReference) -> Self {
        AsymmetricIndexBuilder { reference }
    }

    fn kmers(&'a self, visitor: &mut impl FnMut(Kmer<K>, u32, u32)) {
        for chrom_idx in 0..self.reference.num_chroms() {
            let seq = self.reference.sequence(chrom_idx);
            for (pos, fwd, _rev) in
            Kmer::<K>::agnostic_open_syncmer_iter::<S, FnvHasher>(seq.as_ref(), [(); S]) {
                visitor(fwd, chrom_idx as u32, pos as u32);
            }
        }
    }

    fn build(&'a self) -> Self::IndexType {

        let chrom_info: Vec<ChromInfo> = (0..self.reference.num_chroms())
            .map(|i| self.reference.chrom_info(i).clone())
            .collect();

        let mut unique_seeds: Table<u64, u64> = Table::new();
        let mut nonunique_seeds: HashMap<u64, Vec<u64>> = HashMap::new();

        let mut syncmer_count = 0;

        self.kmers(&mut |fwd, chrom_idx, pos| {
            syncmer_count += 1;
            let locus: u64 = Locus::pack(chrom_idx as usize, pos as usize).into();
            if let Some(other_locus) = unique_seeds.remove(&fwd.0) {
                nonunique_seeds.insert(fwd.0, vec![other_locus, locus]);
            } else if let Some(loci) = nonunique_seeds.get_mut(&fwd.0) {
                loci.push(locus);
            } else {
                unique_seeds.insert(fwd.0, locus);
            }
        });

        let unique_seeds = FrozenTable::from_table(unique_seeds);
        let nonunique_seeds = FrozenBigTable::from_hashmap(nonunique_seeds);

        AsymmetricIndex {
            chrom_info,
            unique_seeds,
            nonunique_seeds,
            query_buffers: Pool::new(),
        }
    }
}

