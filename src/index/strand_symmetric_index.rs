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
use crate::utils::hasher::{FnvHasher, seeded_mix64};
use crate::utils::pool::Pool;
use crate::utils::swiss;

use super::BedRegions;

// Version 2: renamed from "asymmetric-syncmer"; the non-unique seed table's
// hash and probe order changed.
// Version 3: the unique seed table hashes with a fixed function instead of
// `std`'s `DefaultHasher`, which is not stable across Rust releases.
const INDEX_VERSION: u32 = 3;
const INDEX_TYPE: &str = "strand-symmetric-syncmer";

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
/// Strand-symmetric: a reference position is indexed if either its k-mer or
/// that k-mer's reverse complement is an open syncmer, so a read and its
/// reverse complement select mirror-image seed positions.  The forward k-mer
/// is stored; reverse-strand hits are found by querying the read's reverse
/// complement.
///
/// This is the production index structure, optimized for fast lookups and
/// supporting save/load to Parquet files for persistence.
pub struct StrandSymmetricIndex<const K: usize, const S: usize> {
    chrom_info: Vec<ChromInfo>,
    unique_seeds: FrozenTable,
    nonunique_seeds: FrozenBigTable,
    query_buffers: Pool<Vec<(usize, u64)>>,
}

impl<const K: usize, const S: usize> StrandSymmetricIndex<K, S> {
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
        Ok(StrandSymmetricIndex {
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
        Ok(StrandSymmetricIndex {
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
        crate::utils::try_join_all(vec![
            Box::new(|| self.unique_seeds.save_to_feather_directory(dir.join("unique_seeds"))),
            Box::new(|| self.nonunique_seeds.save_to_feather_directory(dir.join("nonunique_seeds"))),
        ])?;
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
        crate::utils::try_join_all(vec![
            Box::new(|| self.unique_seeds.save_to_directory(dir.join("unique_seeds"))),
            Box::new(|| self.nonunique_seeds.save_to_directory(dir.join("nonunique_seeds"))),
        ])?;
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

impl<const K: usize, const S: usize> super::Index for StrandSymmetricIndex<K, S> {
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

impl<const K: usize, const S: usize> LoadableIndex for StrandSymmetricIndex<K, S> {
    fn load(path: &Path) -> std::io::Result<Self> {
        Self::load(path)
    }
}
impl<const K: usize, const S: usize> SyncmerIndex<K, S> for StrandSymmetricIndex<K, S> {
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
    let index = StrandSymmetricIndex::<K, S>::load(path)?;
    Ok(Arc::new(index))
}

/// Builds a `StrandSymmetricIndex` from an in-memory reference.
///
/// The build runs in three parallel phases:
///
/// 1. The reference is cut into pieces, and each piece's syncmers are
///    scattered as `[kmer, locus]` pairs into partitions by the high bits of
///    the k-mer's table hash.  A counting pass first sizes every
///    (partition, piece) run exactly, so the pairs fill one flat array.
/// 2. Each partition is sorted, bringing a k-mer's loci together in genome
///    order: a run of one is a unique seed, longer runs are nonunique.
/// 3. The frozen tables are filled at their final size.  Table slots are
///    chosen by the high bits of the hash too, so each partition's unique
///    seeds land in its own range of slots and the ranges fill concurrently.
pub struct StrandSymmetricIndexBuilder<'a, const K: usize, const S: usize> {
    reference: &'a InMemoryReference
}

/// A stretch of one chromosome whose syncmers start in `start..end`.
#[derive(Clone, Copy)]
struct Piece {
    chrom: usize,
    start: usize,
    end: usize,
    /// End of the sequence the k-mers may extend into (the region's end).
    limit: usize,
}

/// Raw pointer to a slice that threads write into at disjoint indices.
struct SharedSlice<T> {
    ptr: *mut T,
    len: usize,
}

unsafe impl<T: Send> Sync for SharedSlice<T> {}

impl<T> SharedSlice<T> {
    fn new(slice: &mut [T]) -> Self {
        SharedSlice { ptr: slice.as_mut_ptr(), len: slice.len() }
    }

    /// SAFETY: no other thread may access index `i` while the slice is shared.
    unsafe fn write(&self, i: usize, value: T) {
        assert!(i < self.len);
        unsafe { self.ptr.add(i).write(value) }
    }
}

/// Run `f` over `0..n` on `threads` threads, handing out indices in order.
fn parallel_for(n: usize, threads: usize, f: impl Fn(usize) + Sync) {
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads.min(n).max(1) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= n {
                    break;
                }
                f(i);
            });
        }
    });
}

/// Like `parallel_for`, but collects each index's result in order.
fn parallel_map<T: Send>(n: usize, threads: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let results: Vec<std::sync::Mutex<Option<T>>> = (0..n).map(|_| Default::default()).collect();
    parallel_for(n, threads, |i| *results[i].lock().unwrap() = Some(f(i)));
    results
        .into_iter()
        .map(|r| r.into_inner().unwrap().unwrap())
        .collect()
}

/// Seed shared by both frozen tables' hashes.
const TABLE_SEED: u64 = 0xcbf2_9ce4_8422_2325;

/// Partitions are chosen by the top `PARTITION_BITS` bits of the table hash.
const PARTITION_BITS: u32 = 12;
const PARTITIONS: usize = 1 << PARTITION_BITS;

/// Syncmer start positions per piece of reference.
const PIECE_LEN: usize = 1 << 22;

/// Smallest range of unique-table slots filled by one task.
const MIN_RANGE_SLOTS: usize = 1 << 12;

#[inline]
fn partition_of(hash: u64) -> usize {
    (hash >> (64 - PARTITION_BITS)) as usize
}

/// Unique and nonunique seed counts for one sorted partition.
#[derive(Default, Clone, Copy)]
struct PartitionCounts {
    unique: usize,
    nonunique: usize,
    nonunique_loci: usize,
}

/// Call `f(kmer, loci)` for each run of equal k-mers in a sorted partition.
fn for_each_run(pairs: &[[u64; 2]], mut f: impl FnMut(u64, &[[u64; 2]])) {
    let mut i = 0;
    while i < pairs.len() {
        let kmer = pairs[i][0];
        let mut j = i + 1;
        while j < pairs.len() && pairs[j][0] == kmer {
            j += 1;
        }
        f(kmer, &pairs[i..j]);
        i = j;
    }
}

impl<'a, const K: usize, const S: usize> StrandSymmetricIndexBuilder<'a, K, S> {
    /// Cut the reference (or just the BED regions, if given) into pieces.
    fn pieces(&self, bed_regions: Option<&BedRegions>, piece_len: usize) -> Vec<Piece> {
        let mut pieces = Vec::new();
        for chrom in 0..self.reference.num_chroms() {
            let n = self.reference.sequence(chrom).len();
            let whole = [(0, n)];
            let regions: &[(usize, usize)] = match bed_regions {
                None => &whole,
                Some(bed) => bed
                    .get(&self.reference.chrom_info(chrom).name)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]),
            };
            for &(start, end) in regions {
                let (start, end) = (start.min(n), end.min(n));
                for piece_start in (start..end).step_by(piece_len) {
                    pieces.push(Piece {
                        chrom,
                        start: piece_start,
                        end: (piece_start + piece_len).min(end),
                        limit: end,
                    });
                }
            }
        }
        pieces
    }

    /// Visit the syncmers starting in a piece.
    #[inline]
    fn piece_kmers(&self, piece: Piece, mut visitor: impl FnMut(u64, u64)) {
        let seq = self.reference.sequence(piece.chrom);
        let scan_end = (piece.end + K - 1).min(piece.limit);
        let scan = &seq[piece.start..scan_end];
        for (pos, fwd, _rev) in Kmer::<K>::agnostic_open_syncmer_iter::<S, FnvHasher>(scan, [(); S]) {
            let locus: u64 = Locus::pack(piece.chrom, piece.start + pos).into();
            visitor(fwd.0, locus);
        }
    }

    /// Build the index on `threads` threads, indexing only syncmers that lie
    /// entirely within `bed_regions` when it is given.
    pub fn build_parallel(
        &self,
        bed_regions: Option<&BedRegions>,
        threads: usize,
    ) -> StrandSymmetricIndex<K, S> {
        self.build_pieces(bed_regions, threads, PIECE_LEN)
    }

    fn build_pieces(
        &self,
        bed_regions: Option<&BedRegions>,
        threads: usize,
        piece_len: usize,
    ) -> StrandSymmetricIndex<K, S> {
        let threads = threads.max(1);
        let now = std::time::Instant::now();
        let pieces = self.pieces(bed_regions, piece_len);
        log::info!(
            "Building strand-symmetric index: {} pieces on {} threads",
            pieces.len(),
            threads
        );

        // Phase 1a: count each piece's syncmers per partition.
        let piece_counts: Vec<Vec<u32>> = parallel_map(pieces.len(), threads, |p| {
            let mut counts = vec![0u32; PARTITIONS];
            self.piece_kmers(pieces[p], |kmer, _| {
                counts[partition_of(seeded_mix64(kmer, TABLE_SEED))] += 1;
            });
            counts
        });

        // Lay the runs out partition-major, so each partition is contiguous.
        // starts[p * PARTITIONS + q] is where piece p's run for partition q
        // begins, and bounds[q] where partition q begins.
        let mut starts = vec![0usize; pieces.len() * PARTITIONS];
        let mut bounds = Vec::with_capacity(PARTITIONS + 1);
        let mut total = 0usize;
        for q in 0..PARTITIONS {
            bounds.push(total);
            for p in 0..pieces.len() {
                starts[p * PARTITIONS + q] = total;
                total += piece_counts[p][q] as usize;
            }
        }
        bounds.push(total);

        let mut chrom_info: Vec<ChromInfo> = (0..self.reference.num_chroms())
            .map(|i| {
                let mut info = self.reference.chrom_info(i).clone();
                info.length = self.reference.chrom_length(i);
                info.syncmer_count = 0;
                info
            })
            .collect();
        for (piece, counts) in pieces.iter().zip(&piece_counts) {
            chrom_info[piece.chrom].syncmer_count +=
                counts.iter().map(|&c| c as u64).sum::<u64>();
        }
        drop(piece_counts);

        // Phase 1b: scatter the [kmer, locus] pairs into their partitions.
        let mut pairs = vec![[0u64; 2]; total];
        {
            let shared = SharedSlice::new(&mut pairs);
            parallel_for(pieces.len(), threads, |p| {
                let mut cursors = starts[p * PARTITIONS..(p + 1) * PARTITIONS].to_vec();
                self.piece_kmers(pieces[p], |kmer, locus| {
                    let q = partition_of(seeded_mix64(kmer, TABLE_SEED));
                    // SAFETY: piece p's run for partition q is
                    // starts[p * PARTITIONS + q] plus its count, which no
                    // other piece's run overlaps.
                    unsafe { shared.write(cursors[q], [kmer, locus]) };
                    cursors[q] += 1;
                });
            });
        }
        drop(starts);
        log::info!(
            "Collected {} syncmers in {:.2}s",
            total,
            now.elapsed().as_secs_f64()
        );

        // Phase 2: sort each partition and count its unique and nonunique seeds.
        let mut partitions: Vec<&mut [[u64; 2]]> = Vec::with_capacity(PARTITIONS);
        let mut rest = pairs.as_mut_slice();
        for q in 0..PARTITIONS {
            let (part, tail) = rest.split_at_mut(bounds[q + 1] - bounds[q]);
            partitions.push(part);
            rest = tail;
        }
        let partitions: Vec<std::sync::Mutex<&mut [[u64; 2]]>> =
            partitions.into_iter().map(std::sync::Mutex::new).collect();
        let counts: Vec<PartitionCounts> = parallel_map(PARTITIONS, threads, |q| {
            let mut part = partitions[q].lock().unwrap();
            part.sort_unstable();
            let mut counts = PartitionCounts::default();
            for_each_run(&part, |_, run| {
                if run.len() == 1 {
                    counts.unique += 1;
                } else {
                    counts.nonunique += 1;
                    counts.nonunique_loci += run.len();
                }
            });
            counts
        });
        drop(partitions);
        let partitions: Vec<&[[u64; 2]]> =
            (0..PARTITIONS).map(|q| &pairs[bounds[q]..bounds[q + 1]]).collect();
        log::info!("Sorted partitions in {:.2}s", now.elapsed().as_secs_f64());

        // Phase 3: fill the frozen tables, the nonunique one on its own thread
        // while the rest fill the unique one.
        let unique_count: usize = counts.iter().map(|c| c.unique).sum();
        let (unique_seeds, nonunique_seeds) = std::thread::scope(|scope| {
            let nonunique = scope.spawn(|| Self::fill_nonunique(&partitions, &counts));
            let unique = Self::fill_unique(&partitions, unique_count, threads.saturating_sub(1).max(1));
            (unique, nonunique.join().expect("nonunique table thread panicked"))
        });
        drop(partitions);
        drop(pairs);

        log::info!(
            "Built strand-symmetric index: {} unique seeds, {} nonunique seeds in {:.2}s",
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        StrandSymmetricIndex {
            chrom_info,
            unique_seeds,
            nonunique_seeds,
            query_buffers: Pool::new(),
        }
    }

    /// Fill the unique seed table.  Its slots are split into ranges, each
    /// filled by one task from the partitions whose hashes map into it; the
    /// few seeds whose probe sequence leaves their range are placed afterwards.
    fn fill_unique(partitions: &[&[[u64; 2]]], count: usize, threads: usize) -> FrozenTable {
        if count == 0 {
            return FrozenTable::empty();
        }
        let bits = FrozenTable::compute_bits(count);
        let capacity = 1usize << bits;
        let ranges = (capacity / MIN_RANGE_SLOTS).clamp(1, PARTITIONS);
        let range_len = capacity / ranges;
        let partitions_per_range = PARTITIONS / ranges;

        let mut ctrl = vec![swiss::CTRL_EMPTY; capacity + swiss::PROBE_GROUP];
        let mut keys = vec![0u64; capacity];
        let mut values = vec![0u64; capacity];

        let tasks: Vec<_> = ctrl[..capacity]
            .chunks_mut(range_len)
            .zip(keys.chunks_mut(range_len))
            .zip(values.chunks_mut(range_len))
            .map(|((c, k), v)| std::sync::Mutex::new((c, k, v)))
            .collect();
        let spilled: Vec<Vec<(u64, u64)>> = parallel_map(ranges, threads, |r| {
            let mut task = tasks[r].lock().unwrap();
            let (ctrl, keys, values) = &mut *task;
            let lo = r * range_len;
            let mut spilled = Vec::new();
            for part in &partitions[r * partitions_per_range..(r + 1) * partitions_per_range] {
                for_each_run(part, |kmer, run| {
                    if run.len() != 1 {
                        return;
                    }
                    let locus = run[0][1];
                    let hash = seeded_mix64(kmer, TABLE_SEED);
                    match swiss::find_empty_slot_within(ctrl, lo, hash, bits) {
                        Some(slot) => {
                            ctrl[slot - lo] = swiss::h2(hash);
                            keys[slot - lo] = kmer;
                            values[slot - lo] = locus;
                        }
                        None => spilled.push((kmer, locus)),
                    }
                });
            }
            spilled
        });
        drop(tasks);

        log::debug!(
            "{} of {} unique seeds spilled out of their slot range",
            spilled.iter().map(Vec::len).sum::<usize>(),
            count
        );
        for (kmer, locus) in spilled.into_iter().flatten() {
            let hash = seeded_mix64(kmer, TABLE_SEED);
            let slot = swiss::find_empty_slot(&ctrl[..capacity], hash, bits);
            ctrl[slot] = swiss::h2(hash);
            keys[slot] = kmer;
            values[slot] = locus;
        }

        FrozenTable::from_parts(count, TABLE_SEED, bits, ctrl, keys, values)
    }

    /// Fill the nonunique seed table, whose loci are stored contiguously in
    /// slot order.
    fn fill_nonunique(partitions: &[&[[u64; 2]]], counts: &[PartitionCounts]) -> FrozenBigTable {
        let count: usize = counts.iter().map(|c| c.nonunique).sum();
        let total_loci: usize = counts.iter().map(|c| c.nonunique_loci).sum();
        if count == 0 {
            return FrozenBigTable::empty();
        }
        let bits = FrozenBigTable::compute_bits(count);
        let capacity = 1usize << bits;

        // Place the keys, recording each one's slot and its loci count.
        let mut ctrl = vec![swiss::CTRL_EMPTY; capacity];
        let mut keys = vec![0u64; capacity];
        let mut offsets = vec![0u64; capacity + 1];
        let mut slots = Vec::with_capacity(count);
        for part in partitions {
            for_each_run(part, |kmer, run| {
                if run.len() == 1 {
                    return;
                }
                let hash = seeded_mix64(kmer, TABLE_SEED);
                let slot = swiss::find_empty_slot(&ctrl, hash, bits);
                ctrl[slot] = swiss::h2(hash);
                keys[slot] = kmer;
                offsets[slot] = run.len() as u64;
                slots.push(slot);
            });
        }

        // Turn the counts into starting offsets.
        let mut offset = 0u64;
        for o in offsets.iter_mut() {
            let n = *o;
            *o = offset;
            offset += n;
        }
        debug_assert_eq!(offset as usize, total_loci);

        // Copy the loci, visiting the runs in the same order as above.
        let mut values = vec![0u64; total_loci];
        let mut slots = slots.into_iter();
        for part in partitions {
            for_each_run(part, |_, run| {
                if run.len() == 1 {
                    return;
                }
                let start = offsets[slots.next().unwrap()] as usize;
                for (dst, pair) in values[start..start + run.len()].iter_mut().zip(run) {
                    *dst = pair[1];
                }
            });
        }

        FrozenBigTable::from_parts(count, TABLE_SEED, bits, ctrl, keys, offsets, values)
    }
}

impl<'a, const K: usize, const S: usize> super::IndexBuilder<'a, K, S> for StrandSymmetricIndexBuilder<'a, K, S> {
    type IndexType = StrandSymmetricIndex<K, S>;

    fn make(reference: &'a InMemoryReference) -> Self {
        StrandSymmetricIndexBuilder { reference }
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
        let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
        self.build_parallel(None, threads)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::index::IndexBuilder;

    /// The index a straightforward serial build produces: every syncmer lying
    /// within the regions, keyed by k-mer, with loci in genome order.
    fn naive<const K: usize, const S: usize>(
        reference: &InMemoryReference,
        bed: Option<&BedRegions>,
    ) -> HashMap<u64, Vec<u64>> {
        let mut seeds: HashMap<u64, Vec<u64>> = HashMap::new();
        for chrom in 0..reference.num_chroms() {
            let seq = reference.sequence(chrom);
            let whole = vec![(0, seq.len())];
            let regions = match bed {
                None => &whole,
                Some(bed) => match bed.get(&reference.chrom_info(chrom).name) {
                    Some(r) => r,
                    None => continue,
                },
            };
            for &(start, end) in regions {
                let (start, end) = (start.min(seq.len()), end.min(seq.len()));
                for (pos, fwd, _) in Kmer::<K>::agnostic_open_syncmer_iter::<S, FnvHasher>(
                    &seq[start..end],
                    [(); S],
                ) {
                    let locus: u64 = Locus::pack(chrom, start + pos).into();
                    seeds.entry(fwd.0).or_default().push(locus);
                }
            }
        }
        seeds
    }

    fn assert_matches<const K: usize, const S: usize>(
        index: &StrandSymmetricIndex<K, S>,
        expected: &HashMap<u64, Vec<u64>>,
    ) {
        let unique = expected.values().filter(|l| l.len() == 1).count();
        assert_eq!(index.unique_seeds.len(), unique);
        assert_eq!(index.nonunique_seeds.len(), expected.len() - unique);
        for (&kmer, loci) in expected {
            let hit = index.lookup_kmer(kmer).expect("indexed k-mer not found");
            assert_eq!(hit.loci, loci.as_slice(), "loci differ for k-mer {kmer:x}");
        }
        // Every seed the index holds is one the naive build has.
        assert_eq!(index.iter().count(), expected.len());
        // A k-mer absent from the reference misses.
        let absent = (0u64..).find(|k| !expected.contains_key(k)).unwrap();
        assert!(index.lookup_kmer(absent).is_none());
    }

    fn hox() -> InMemoryReference {
        InMemoryReference::load("tests/data/hox.fa", false).unwrap()
    }

    #[test]
    fn parallel_build_matches_naive() {
        let _ = env_logger::builder().is_test(true).try_init();
        let reference = hox();
        let expected = naive::<20, 10>(&reference, None);
        let builder = StrandSymmetricIndexBuilder::<20, 10>::make(&reference);
        // Short pieces put many syncmers near piece boundaries.
        for (threads, piece_len) in [(1, PIECE_LEN), (4, 1000), (7, 37)] {
            let index = builder.build_pieces(None, threads, piece_len);
            assert_matches(&index, &expected);
            let syncmers: u64 = index.chrom_info.iter().map(|c| c.syncmer_count).sum();
            assert_eq!(syncmers as usize, expected.values().map(Vec::len).sum::<usize>());
        }
    }

    #[test]
    fn parallel_build_respects_bed_regions() {
        let reference = hox();
        let name = reference.chrom_info(0).name.clone();
        let mut bed = BedRegions::new();
        bed.insert(name, vec![(0, 5000), (5003, 5030), (20_000, 90_000), (100_000, usize::MAX)]);
        let expected = naive::<20, 10>(&reference, Some(&bed));
        let builder = StrandSymmetricIndexBuilder::<20, 10>::make(&reference);
        for piece_len in [PIECE_LEN, 999] {
            assert_matches(&builder.build_pieces(Some(&bed), 3, piece_len), &expected);
        }
    }

    #[test]
    fn index_round_trips_through_disk() {
        let reference = hox();
        let expected = naive::<20, 10>(&reference, None);
        let index = StrandSymmetricIndexBuilder::<20, 10>::make(&reference).build_parallel(None, 4);
        for portable in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            index.save(dir.path(), portable).unwrap();
            let loaded = StrandSymmetricIndex::<20, 10>::load(dir.path()).unwrap();
            assert_matches(&loaded, &expected);
        }
    }
}
