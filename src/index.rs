use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

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

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Combine two strand orientations.
    #[inline]
    pub fn combine(&self, other: &Strand) -> Strand {
        match (self, other) {
            (Strand::Forward, Strand::Forward) | (Strand::Reverse, Strand::Reverse) => {
                Strand::Forward
            }
            (Strand::Forward, Strand::Reverse) | (Strand::Reverse, Strand::Forward) => {
                Strand::Reverse
            }
        }
    }

    #[inline]
    pub fn as_char(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }

    #[inline]
    pub fn is_reverse(self) -> bool {
        self == Strand::Reverse
    }

    #[inline]
    pub fn from_is_reverse(is_reverse: bool) -> Strand {
        if is_reverse {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
}

pub trait PackedLocus: From<u64> + Into<u64> + Sized {
    fn pack(chrom: usize, pos: usize, strand: Strand) -> Self;

    fn unpack(&self) -> (usize, usize, Strand);
}

// =============================================================================
// Index - Immutable, frozen index for fast lookups
// =============================================================================

pub trait Index: Sized + Send + Sync {
    type LocusType: PackedLocus;

    /// Load the index from disk. Reads metadata.json to validate and dispatch format.
    fn load<P: AsRef<Path>>(path: P) -> std::io::Result<Self>;

    /// Save the index to disk. If `portable` is true, writes Arrow IPC (Feather)
    /// format; otherwise writes native endian binary format.
    fn save<P: AsRef<Path>>(&self, path: P, portable: bool) -> std::io::Result<()>;

    /// Return metadata for the chromosome at the given index.
    fn chrom_info(&self, chrom_idx: usize) -> &ChromInfo;

    /// Return metadata for all chromosomes in the index.
    fn all_chrom_info(&self) -> &[ChromInfo];

    /// Find all the indexed seed hits for a given sequence, and invoke the
    /// callback with the hits.
    fn find_seeds<F>(&self, seq: &[u8], callback: F)
    where
        F: FnMut(usize, u64, usize, &[Self::LocusType]);
}

pub trait SyncmerIndex<const K: usize, const S: usize>: Index {
    /// Look up a single kmer. The callback receives (hit_count, loci).
    fn with<F: FnMut(usize, &[Self::LocusType])>(&self, kmer: &Kmer<K>, f: F);

    /// Look up a batch of kmers, calling the callback for each hit.
    /// The callback receives (read_pos, kmer_val, hit_count, loci).
    fn lookup_batch<F>(&self, batch: &[(usize, u64)], callback: F)
    where
        F: FnMut(usize, u64, usize, &[Self::LocusType]);
}

pub trait IndexBuilder<const K: usize, const S: usize> {
    type IndexType: Index;

    /// Build the index from the reference sequence.
    fn build(reference: &InMemoryReference) -> Self::IndexType;
}

/// Read `metadata.json` from an index directory and return the `index_type`
/// field if present. Returns `Ok(None)` if the file doesn't exist or has no
/// `index_type` field. Returns an error only on I/O or JSON parse failure.
pub fn probe_index_kind<P: AsRef<Path>>(path: P) -> std::io::Result<Option<String>> {
    let metadata_path = path.as_ref().join("metadata.json");
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

pub mod asymmetric_index;
pub mod fwd_index;
