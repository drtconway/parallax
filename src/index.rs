use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use noodles::bed;

use crate::kmers::Kmer;
use crate::reference::{ChromInfo, InMemoryReference};

// =============================================================================
// BED file handling for region masking
// =============================================================================

/// Regions of interest loaded from a BED file, organized by chromosome name.
/// Each chromosome has a sorted vector of (start, end) intervals (0-based, half-open).
pub type BedRegions = HashMap<String, Vec<(usize, usize)>>;

/// Load BED regions from a file.
/// Returns a map from chromosome name to sorted, non-overlapping intervals.
pub fn load_bed_regions<P: AsRef<Path>>(path: P) -> std::io::Result<BedRegions> {
    let path = path.as_ref();
    log::info!("Loading BED regions from {}", path.display());

    let file = File::open(path)?;
    let mut reader: bed::io::Reader<3, BufReader<File>> =
        bed::io::Reader::new(BufReader::new(file));

    let mut regions: BedRegions = HashMap::new();

    let mut record: bed::Record<3> = bed::Record::default();
    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        let chrom = record.reference_sequence_name().to_string();
        let start: usize = record.feature_start()?.get() - 1; // Convert from 1-based to 0-based
        let end: usize = record.feature_end().unwrap()?.get(); // Already 0-based half-open in BED

        regions.entry(chrom).or_default().push((start, end));
    }

    // Sort and merge overlapping intervals for each chromosome
    let mut total_regions = 0;
    let mut total_bases = 0usize;
    for intervals in regions.values_mut() {
        intervals.sort_by_key(|&(start, _)| start);

        // Merge overlapping intervals
        let mut merged = Vec::with_capacity(intervals.len());
        for &(start, end) in intervals.iter() {
            if let Some(&mut (ref mut last_start, ref mut last_end)) = merged.last_mut() {
                if start <= *last_end {
                    // Overlapping or adjacent, extend
                    *last_end = (*last_end).max(end);
                    continue;
                }
                let _ = last_start; // silence warning
            }
            merged.push((start, end));
        }

        total_regions += merged.len();
        total_bases += merged.iter().map(|(s, e)| e - s).sum::<usize>();
        *intervals = merged;
    }

    log::info!(
        "Loaded {} regions covering {} bp across {} chromosomes",
        total_regions,
        total_bases,
        regions.len()
    );

    Ok(regions)
}

pub trait PackedLocus: From<u64> + Into<u64> + Sized {
    fn pack(chrom: usize, pos: usize) -> Self;

    fn unpack(&self) -> (usize, usize);
}

pub struct IndexHit<'a> {
    pub query_pos: usize,
    pub seed_kmer: u64,
    pub loci: &'a [u64],
    pub k: usize,
}

// =============================================================================
// Index - Immutable, frozen index for fast lookups
// =============================================================================

pub trait Index: Send + Sync {
    /// Save the index to disk. If `portable` is true, writes Arrow IPC (Feather)
    /// format; otherwise writes native endian binary format.
    fn save(&self, path: &Path, portable: bool) -> std::io::Result<()>;

    /// Return metadata for the chromosome at the given index.
    fn chrom_info(&self, chrom_idx: usize) -> &ChromInfo;

    /// Return metadata for all chromosomes in the index.
    fn all_chrom_info(&self) -> &[ChromInfo];

    /// Return the k-mer size for IndexHits.
    fn k(&self) -> usize;

    /// Find all the indexed seed hits for a given sequence, and invoke the
    /// callback with the hits.
    fn find_seeds(&self, seq: &[u8], callback: &mut dyn FnMut(IndexHit<'_>));

    /// Look up a single k-mer by its raw u64 value.
    ///
    /// The kmer value is implicitly treated as a forward-strand query. Callers
    /// that are processing a reverse-complement sequence must combine the query
    /// strand with the hit strand from `unpack_locus` themselves.
    fn lookup_kmer(&self, kmer: u64) -> Option<IndexHit<'_>>;

    /// Unpack a locus returned by the index, returning (chrom_id, position)
    /// always with respect to the forward strand.
    fn unpack_locus(&self, locus: u64) -> (usize, usize);

    /// Unpack a slice of loci
    fn unpack_loci(&self, packed: &[u64], unpacked: &mut Vec<(usize, usize)>) {
        unpacked.clear();
        unpacked.reserve(packed.len());
        for &value in packed {
            let locus = self.unpack_locus(value);
            unpacked.push(locus);
        }
    }

    /// Iterate over all the seeds in the index in unspecified order.
    fn iter(&self) -> Box<dyn Iterator<Item = IndexHit<'_>> + '_>;

    /// Validate that the index is compatible with the given reference.
    fn validate_reference(&self, reference: &InMemoryReference) -> std::io::Result<()> {
        let chrom_info = self.all_chrom_info();
        let idx_n = chrom_info.len();
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
            let idx_name = &chrom_info[i].name;
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
            let idx_len = chrom_info[i].length;
            let ref_len = reference.chrom_length(i);
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
}

pub trait LoadableIndex: Index + Sized {
    /// Load the index from disk. Reads metadata.json to validate and dispatch format.
    fn load(path: &Path) -> std::io::Result<Self>;
}

pub trait SyncmerIndex<const K: usize, const S: usize>: Index {
    /// Look up a single kmer. The callback receives (hit_count, loci).
    fn with<F: FnMut(IndexHit<'_>)>(&self, kmer: &Kmer<K>, f: F);

    /// Look up a batch of kmers, calling the callback for each hit.
    /// The callback receives (read_pos, kmer_val, hit_count, loci).
    fn lookup_batch<F>(&self, batch: &[(usize, u64)], callback: F)
    where
        F: FnMut(IndexHit<'_>);
}

pub trait IndexBuilder<'a, const K: usize, const S: usize> {
    type IndexType: Index;

    fn make(reference: &'a InMemoryReference) -> Self;

    /// Visit the kmers that will be indexed.
    fn kmers(&'a self, visitor: &mut impl FnMut(Kmer<K>, u32, u32));

    /// Build the index from the reference sequence.
    fn build(&'a self) -> Self::IndexType;
}

/// Read `metadata.json` from an index directory and return the `index_type`
/// field if present. Returns `Ok(None)` if the file doesn't exist or has no
/// `index_type` field. Returns an error only on I/O or JSON parse failure.
pub fn probe_index_kind(path: &Path) -> std::io::Result<Option<String>> {
    let metadata_path = path.join("metadata.json");
    if !metadata_path.exists() {
        return Ok(None);
    }
    let content = std::fs::read_to_string(&metadata_path)?;
    let value: serde_json::Value = serde_json::from_str(&content)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    Ok(value
        .get("index_type")
        .and_then(|v| v.as_str())
        .map(str::to_owned))
}

pub fn load_index(path: &Path) -> std::io::Result<Arc<dyn Index>> {
    if let Some(kind) = probe_index_kind(path)? {
        if kind == "fwd-syncmer" {
            let index = fwd_index::load_index(path)?;
            return Ok(index);
        }
        if kind == "asymmetric-syncmer" {
            let index = asymmetric_index::load_index(path)?;
            return Ok(index)
        }
    }
    Err(std::io::Error::other("couldn't determine index type"))
}

pub mod asymmetric_index;
pub mod fwd_index;
