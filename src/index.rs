use std::collections::HashMap;
use std::io::BufReader;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use crossbeam::channel;
use noodles::bed;

use crate::{
    kmers::Kmer,
    reference::{ChromInfo, InMemoryReference},
    utils::{
        Selection, frozen_big_table::FrozenBigTable, frozen_table::FrozenTable, hasher::{FnvHasher}, table::Table
    },
};

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
    let mut reader: bed::io::Reader<3, BufReader<File>> = bed::io::Reader::new(BufReader::new(file));
    
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

// =============================================================================
// Index - Immutable, frozen index for fast lookups
// =============================================================================

/// Immutable K-mer seed index backed by Arrow arrays.
///
/// This is the production index structure, optimized for fast lookups and
/// supporting save/load to Parquet files for persistence.
pub struct Index<const K: usize, const S: usize> {
    chrom_info: Vec<ChromInfo>,
    cumulative_lengths: Vec<u32>,
    unique_seeds: FrozenTable,
    nonunique_seeds: FrozenBigTable,
}

impl<const K: usize, const S: usize> Index<K, S> {
    /// Load an index from a directory containing Parquet files.
    pub fn load<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();

        let now = std::time::Instant::now();

        // Load chromosome metadata
        let chrom_info = Self::load_chrom_info(dir.join("chrom_info.json"))?;
        let cumulative_lengths = Self::load_cumulative_lengths(dir.join("cumulative_lengths.bin"))?;

        // Load seed tables
        let unique_seeds = FrozenTable::load_from_directory(dir.join("unique_seeds"))?;
        let nonunique_seeds = FrozenBigTable::load_from_directory(dir.join("nonunique_seeds"))?;

        log::info!(
            "Loaded index: {} chromosomes, {} unique seeds, {} nonunique seeds in {:.2?}s",
            chrom_info.len(),
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        Ok(Index {
            chrom_info,
            cumulative_lengths,
            unique_seeds,
            nonunique_seeds,
        })
    }

    /// Save the index to a directory as Parquet files.
    pub fn save<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;

        // Save chromosome metadata
        Self::save_chrom_info(&self.chrom_info, dir.join("chrom_info.json"))?;
        Self::save_cumulative_lengths(&self.cumulative_lengths, dir.join("cumulative_lengths.bin"))?;

        // Save seed tables
        self.unique_seeds.save_to_directory(dir.join("unique_seeds"))?;
        self.nonunique_seeds.save_to_directory(dir.join("nonunique_seeds"))?;

        log::info!(
            "Saved index: {} chromosomes, {} unique seeds, {} nonunique seeds",
            self.chrom_info.len(),
            self.unique_seeds.len(),
            self.nonunique_seeds.len()
        );

        Ok(())
    }

    /// Load an index from a directory containing Arrow IPC (Feather) files.
    pub fn load_feather<P: AsRef<Path>>(dir: P) -> std::io::Result<Self> {
        let dir = dir.as_ref();

        let now = std::time::Instant::now();

        // Load chromosome metadata (same format as Parquet)
        let chrom_info = Self::load_chrom_info(dir.join("chrom_info.json"))?;
        let cumulative_lengths = Self::load_cumulative_lengths(dir.join("cumulative_lengths.bin"))?;

        // Load seed tables from Feather format
        let unique_seeds = FrozenTable::load_from_feather_directory(dir.join("unique_seeds"))?;
        let nonunique_seeds = FrozenBigTable::load_from_feather_directory(dir.join("nonunique_seeds"))?;

        log::info!(
            "Loaded index (feather): {} chromosomes, {} unique seeds, {} nonunique seeds in {:.2?}s",
            chrom_info.len(),
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        Ok(Index {
            chrom_info,
            cumulative_lengths,
            unique_seeds,
            nonunique_seeds,
        })
    }

    /// Save the index to a directory as Arrow IPC (Feather) files.
    ///
    /// Feather format is generally faster to read/write than Parquet but
    /// produces larger files (no compression by default).
    pub fn save_feather<P: AsRef<Path>>(&self, dir: P) -> std::io::Result<()> {
        let dir = dir.as_ref();
        std::fs::create_dir_all(dir)?;

        // Save chromosome metadata (same format as Parquet)
        Self::save_chrom_info(&self.chrom_info, dir.join("chrom_info.json"))?;
        Self::save_cumulative_lengths(&self.cumulative_lengths, dir.join("cumulative_lengths.bin"))?;

        // Save seed tables in Feather format
        self.unique_seeds.save_to_feather_directory(dir.join("unique_seeds"))?;
        self.nonunique_seeds.save_to_feather_directory(dir.join("nonunique_seeds"))?;

        log::info!(
            "Saved index (feather): {} chromosomes, {} unique seeds, {} nonunique seeds",
            self.chrom_info.len(),
            self.unique_seeds.len(),
            self.nonunique_seeds.len()
        );

        Ok(())
    }

    /// Look up a k-mer and call `f` for each matching locus.
    pub fn with<F: FnMut(usize, usize)>(&self, kmer: &Kmer<K>, mut f: F) {
        if let Some(loc) = self.unique_seeds.get(kmer.0) {
            let (chrom_idx, pos) = self.decode_locus(loc);
            f(chrom_idx, pos);
        } else if let Some(locs) = self.nonunique_seeds.get(kmer.0) {
            for &loc in locs {
                let (chrom_idx, pos) = self.decode_locus(loc);
                f(chrom_idx, pos);
            }
        }
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

    /// Get all chromosome info.
    #[allow(dead_code)]
    pub fn all_chrom_info(&self) -> &[ChromInfo] {
        &self.chrom_info
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

    fn decode_locus(&self, abs_pos: u32) -> (usize, usize) {
        let chrom_idx = match self.cumulative_lengths.binary_search(&abs_pos) {
            Ok(idx) => idx,
            Err(idx) => idx - 1,
        };
        let pos = abs_pos - self.cumulative_lengths[chrom_idx];
        (chrom_idx, pos as usize)
    }

    fn load_chrom_info<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<ChromInfo>> {
        let content = std::fs::read_to_string(path)?;
        serde_json::from_str(&content)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    fn save_chrom_info<P: AsRef<Path>>(chrom_info: &[ChromInfo], path: P) -> std::io::Result<()> {
        let content = serde_json::to_string(chrom_info)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        std::fs::write(path, content)
    }

    fn load_cumulative_lengths<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<u32>> {
        let bytes = std::fs::read(path)?;
        let lengths: Vec<u32> = bytes
            .chunks_exact(4)
            .map(|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
            .collect();
        Ok(lengths)
    }

    fn save_cumulative_lengths<P: AsRef<Path>>(lengths: &[u32], path: P) -> std::io::Result<()> {
        let bytes: Vec<u8> = lengths.iter().flat_map(|&n| n.to_le_bytes()).collect();
        std::fs::write(path, bytes)
    }
}

// =============================================================================
// IndexBuilder - Mutable builder for constructing an Index
// =============================================================================

/// Mutable builder for constructing a K-mer seed index.
///
/// Use this to build an index from reference sequences, then call `build()`
/// to create an immutable `Index` for fast lookups.
pub struct IndexBuilder<const K: usize, const S: usize> {
    chrom_info: Vec<ChromInfo>,
    cumulative_lengths: Vec<u32>,
    unique_seeds: Table<u64, u32>,
    nonunique_seeds: HashMap<u64, Vec<u32>>,
}

impl<const K: usize, const S: usize> IndexBuilder<K, S> {
    pub fn new() -> Self {
        IndexBuilder {
            chrom_info: Vec::new(),
            cumulative_lengths: vec![0],
            unique_seeds: Table::new(),
            nonunique_seeds: HashMap::new(),
        }
    }

    /// Create a new builder pre-initialized with chromosome info and lengths.
    /// All locus encoding will use these pre-computed cumulative lengths.
    pub fn new_with_chrom_info(chrom_info: &[(ChromInfo, usize)]) -> Self {
        let chrom_info_vec: Vec<ChromInfo> = chrom_info.iter().map(|(info, _)| info.clone()).collect();
        let mut cumulative_lengths = Vec::with_capacity(chrom_info_vec.len() + 1);
        cumulative_lengths.push(0u32);

        let mut cumulative = 0u32;
        for (_, len) in chrom_info {
            cumulative += *len as u32;
            cumulative_lengths.push(cumulative);
        }

        IndexBuilder {
            chrom_info: chrom_info_vec,
            cumulative_lengths,
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
    pub fn add_chrom_regions<Seq: AsRef<[u8]>>(&mut self, chrom_idx: usize, seq: Seq, regions: Option<&[(usize, usize)]>) {
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
                            let loc = self.encode_locus(chrom_idx, pos);
                            m += 1;
                            self.insert_kmer(kmer.0, loc);
                        }
                        _ => {}
                    }
                }
                let r = (m as f64) / (n as f64);
                log::info!(
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
                    for (rel_pos, sel) in Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(region_seq, [(); S]) {
                        match sel {
                            Selection::Left(kmer) | Selection::Both(kmer, _) => {
                                // Absolute position = region start + relative position
                                let abs_pos = start + rel_pos;
                                let loc = self.encode_locus(chrom_idx, abs_pos);
                                m += 1;
                                self.insert_kmer(kmer.0, loc);
                            }
                            _ => {}
                        }
                    }
                }
                let r = if region_bases > 0 { (m as f64) / (region_bases as f64) } else { 0.0 };
                log::info!(
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
    fn insert_kmer(&mut self, kmer: u64, loc: u32) {
        if let Some(loc0) = self.unique_seeds.remove(&kmer) {
            self.nonunique_seeds.insert(kmer, vec![loc0, loc]);
        } else if let Some(locs) = self.nonunique_seeds.get_mut(&kmer) {
            locs.push(loc);
        } else {
            self.unique_seeds.insert(kmer, loc);
        }
    }

    pub fn add<Seq: AsRef<[u8]>>(&mut self, chrom_info: ChromInfo, seq: Seq) -> usize {
        let idx = self.chrom_info.len();
        self.chrom_info.push(chrom_info);

        let n = seq.as_ref().len();
        let l = n as u32 + self.cumulative_lengths.last().copied().unwrap_or(0);
        self.cumulative_lengths.push(l);

        let mut m = 0usize;
        for (pos, sel) in Kmer::<K>::open_syncmer_iter::<S, FnvHasher>(seq.as_ref(), [(); S]) {
            match sel {
                Selection::Left(kmer) | Selection::Both(kmer, _) => {
                    let loc = self.encode_locus(idx, pos);
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

        let r = (m as f64) / (n as f64);
        log::info!(
            "Indexed chrom {} \"{}\" (length {}, {} seeds; {})",
            idx,
            self.chrom_info[idx].name,
            n,
            m,
            r
        );
        idx
    }

    /// Build an immutable Index from this builder.
    pub fn build(self) -> Index<K, S> {
        let now = std::time::Instant::now();

        let unique_seeds = FrozenTable::from_table(self.unique_seeds);
        let nonunique_seeds = FrozenBigTable::from_hashmap(self.nonunique_seeds);

        log::info!(
            "Built frozen index: {} unique seeds, {} nonunique seeds in {:.2}s",
            unique_seeds.len(),
            nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        Index {
            chrom_info: self.chrom_info,
            cumulative_lengths: self.cumulative_lengths,
            unique_seeds,
            nonunique_seeds,
        }
    }

    fn encode_locus(&self, chrom_idx: usize, pos: usize) -> u32 {
        self.cumulative_lengths[chrom_idx] + pos as u32
    }

    /// Merge another builder into this one.
    /// Both builders must have been created with the same chromosome list
    /// (via `new_with_chrom_info`), so no locus recoding is needed.
    pub fn merge(&mut self, other: IndexBuilder<K, S>) {
        // Verify chromosome lists match
        debug_assert_eq!(
            self.chrom_info.len(),
            other.chrom_info.len(),
            "Builders must have same chromosome list"
        );
        debug_assert_eq!(
            self.cumulative_lengths, other.cumulative_lengths,
            "Builders must have same cumulative lengths"
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

        log::info!(
            "Merged builders, now {} unique and {} nonunique seeds in {:.2}s",
            self.unique_seeds.len(),
            self.nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );
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
    pub fn build_parallel(reference: &InMemoryReference, bed_regions: Option<&BedRegions>, num_threads: usize) -> Index<K, S> {
        let num_chroms = reference.num_chroms();
        if num_chroms == 0 {
            return IndexBuilder::<K, S>::new().build();
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
        let (result_tx, result_rx) = channel::bounded::<(usize, IndexBuilder<K, S>)>(num_chroms);

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
                    let mut builder = IndexBuilder::<K, S>::new_with_chrom_info(&chrom_info);
                    
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
        let mut combined: Option<IndexBuilder<K, S>> = None;
        let mut next_chrom = 0;
        let mut buffered_builders: HashMap<usize, IndexBuilder<K, S>> = HashMap::new();
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

        let combined = combined.unwrap_or_else(|| IndexBuilder::<K, S>::new_with_chrom_info(&chrom_info));

        log::info!(
            "Index complete: {} unique seeds, {} nonunique seeds in {:.2}s",
            combined.unique_seeds.len(),
            combined.nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        combined.build()
    }
}

impl<const K: usize, const S: usize, R: std::io::BufRead> TryFrom<noodles::fasta::io::Reader<R>>
    for IndexBuilder<K, S>
{
    type Error = crate::error::ParallaxError;

    fn try_from(mut reader: noodles::fasta::io::Reader<R>) -> Result<Self, Self::Error> {
        let mut builder = IndexBuilder::new();

        for record in reader.records() {
            let record = record?;
            let name = String::from_utf8(record.name().to_vec())
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let description = record
                .description()
                .map(|d| {
                    String::from_utf8(d.to_vec())
                        .unwrap_or_default()
                })
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
