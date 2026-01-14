use std::collections::HashMap;
use std::sync::Arc;

use crossbeam::channel;

use crate::{
    kmers::Kmer,
    reference::InMemoryReference,
    utils::{Selection, table::Table},
};

/// K-mer seed index for fast sequence lookup.
///
/// This struct contains only the k-mer positions, not the actual sequences.
/// Use `Reference` for accessing the underlying sequence data.
pub struct Index<const K: usize, const S: usize> {
    chroms: Vec<String>,
    cumulative_lengths: Vec<u32>,
    unique_seeds: Table<u64, u32>,
    nonunique_seeds: HashMap<u64, Vec<u32>>,
}

pub struct LocusIter<'a, 'b, const K: usize, const S: usize> {
    index: &'a Index<K, S>,
    inner: LocusInner<'b>,
}

enum LocusInner<'b> {
    One(Option<u32>),
    Many(std::slice::Iter<'b, u32>),
}

impl<'a, 'b, const K: usize, const S: usize> Iterator for LocusIter<'a, 'b, K, S> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.inner {
            LocusInner::One(slot) => slot.take().map(|abs| self.index.decode_locus(abs)),
            LocusInner::Many(iter) => iter.next().map(|&abs| self.index.decode_locus(abs)),
        }
    }
}

impl<const K: usize, const S: usize> Index<K, S> {
    pub fn new() -> Self {
        Index {
            chroms: Vec::new(),
            cumulative_lengths: vec![0],
            unique_seeds: Table::new(),
            nonunique_seeds: HashMap::new(),
        }
    }

    /// Create a new index pre-initialized with chromosome names and lengths.
    /// All locus encoding will use these pre-computed cumulative lengths.
    pub fn new_with_chroms(chrom_info: &[(String, usize)]) -> Self {
        let chroms: Vec<String> = chrom_info.iter().map(|(name, _)| name.clone()).collect();
        let mut cumulative_lengths = Vec::with_capacity(chroms.len() + 1);
        cumulative_lengths.push(0u32);

        let mut cumulative = 0u32;
        for (_, len) in chrom_info {
            cumulative += *len as u32;
            cumulative_lengths.push(cumulative);
        }

        Index {
            chroms,
            cumulative_lengths,
            unique_seeds: Table::new(),
            nonunique_seeds: HashMap::new(),
        }
    }

    /// Add a chromosome sequence by index.
    /// The index must have been created with `new_with_chroms` and
    /// `chrom_idx` must be valid for the chromosome list.
    pub fn add_chrom<Seq: AsRef<[u8]>>(&mut self, chrom_idx: usize, seq: Seq) {
        assert!(chrom_idx < self.chroms.len(), "chrom_idx out of bounds");

        let seq = seq.as_ref();
        let n = seq.len();
        let mut m = 0usize;

        for (pos, sel) in Kmer::<K>::open_syncmer_iter(seq, [(); S]) {
            match sel {
                Selection::Left(kmer) | Selection::Both(kmer, _) => {
                    let loc = self.encode_locus(chrom_idx, pos);
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
            "Indexed chrom {} \"{}\" (length {}, {} seeds; {:.4})",
            chrom_idx,
            self.chroms[chrom_idx],
            n,
            m,
            r
        );
    }

    pub fn add<Seq: AsRef<[u8]>>(&mut self, chrom: String, seq: Seq) -> usize {
        let idx = self.chroms.len();
        self.chroms.push(chrom);

        let n = seq.as_ref().len();
        let l = n as u32 + self.cumulative_lengths.last().copied().unwrap_or(0);
        self.cumulative_lengths.push(l);

        let mut m = 0usize;
        for (pos, sel) in Kmer::<K>::open_syncmer_iter(seq.as_ref(), [(); S]) {
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
            self.chroms[idx],
            n,
            m,
            r
        );
        idx
    }

    #[allow(dead_code)]
    pub fn get(&self, kmer: &Kmer<K>) -> Option<LocusIter<'_, '_, K, S>> {
        let x = kmer.0;
        if let Some(locs) = self.nonunique_seeds.get(&x) {
            Some(LocusIter {
                index: self,
                inner: LocusInner::Many(locs.iter()),
            })
        } else if let Some(loc) = self.unique_seeds.get(&x) {
            Some(LocusIter {
                index: self,
                inner: LocusInner::One(Some(*loc)),
            })
        } else {
            None
        }
    }

    pub fn with<F: FnMut(usize, usize)>(&self, kmer: &Kmer<K>, mut f: F) {
        if let Some(loc) = self.unique_seeds.get(&kmer.0) {
            let (chrom_idx, pos) = self.decode_locus(*loc);
            f(chrom_idx, pos);
        } else if let Some(locs) = self.nonunique_seeds.get(&kmer.0) {
            for &loc in locs {
                let (chrom_idx, pos) = self.decode_locus(loc);
                f(chrom_idx, pos);
            }
        }
    }

    fn encode_locus(&self, chrom_idx: usize, pos: usize) -> u32 {
        self.cumulative_lengths[chrom_idx] + pos as u32
    }

    fn decode_locus(&self, abs_pos: u32) -> (usize, usize) {
        let chrom_idx = match self.cumulative_lengths.binary_search(&abs_pos) {
            Ok(idx) => idx,
            Err(idx) => idx - 1,
        };
        let pos = abs_pos - self.cumulative_lengths[chrom_idx];
        (chrom_idx, pos as usize)
    }

    /// Get the chromosome name
    pub fn chrom_name(&self, chrom_idx: usize) -> &str {
        &self.chroms[chrom_idx]
    }

    /// Merge another index into this one.
    /// Both indexes must have been created with the same chromosome list
    /// (via `new_with_chroms`), so no locus recoding is needed.
    pub fn merge(&mut self, other: Index<K, S>) {
        // Verify chromosome lists match
        debug_assert_eq!(
            self.chroms.len(),
            other.chroms.len(),
            "Indexes must have same chromosome list"
        );
        debug_assert_eq!(
            self.cumulative_lengths, other.cumulative_lengths,
            "Indexes must have same cumulative lengths"
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
            "Merged indexes, now {} unique and {} nonunique seeds in {:.2}s",
            self.unique_seeds.len(),
            self.nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );
    }

    /// Build an index from an InMemoryReference using multiple threads.
    ///
    /// Each chromosome is indexed in parallel, then the per-chromosome indexes
    /// are merged together. All indexes share the same chromosome list, so
    /// merging is fast (no locus recoding required).
    pub fn build_parallel(reference: &InMemoryReference, num_threads: usize) -> Self {
        let num_chroms = reference.num_chroms();
        if num_chroms == 0 {
            return Index::new();
        }

        let now = std::time::Instant::now();

        let num_threads = num_threads.max(1);

        // Build chromosome info for pre-initialization
        let chrom_info: Vec<(String, usize)> = (0..num_chroms)
            .map(|i| {
                (
                    reference.chrom_name(i).to_string(),
                    reference.chrom_length(i) as usize,
                )
            })
            .collect();
        let chrom_info = Arc::new(chrom_info);
        log::info!(
            "Building index for {} chromosomes using {} threads",
            num_chroms,
            num_threads
        );

        // Create work items
        let reference = Arc::new(reference.clone());
        let (work_tx, work_rx) = channel::bounded::<usize>(num_chroms);
        let (result_tx, result_rx) = channel::bounded::<(usize, Index<K, S>)>(num_chroms);

        // Spawn worker threads
        let mut handles = Vec::with_capacity(num_threads);
        for _ in 0..num_threads {
            let work_rx = work_rx.clone();
            let result_tx = result_tx.clone();
            let reference = Arc::clone(&reference);
            let chrom_info = Arc::clone(&chrom_info);

            let handle = std::thread::spawn(move || {
                for chrom_idx in work_rx {
                    let seq = reference.sequence(chrom_idx);

                    // Create index with shared chromosome info
                    let mut index = Index::<K, S>::new_with_chroms(&chrom_info);
                    index.add_chrom(chrom_idx, seq);

                    result_tx
                        .send((chrom_idx, index))
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
        let mut combined: Option<Index<K, S>> = None;
        let mut next_chrom = 0;
        let mut buffered_indexes: HashMap<usize, Index<K, S>> = HashMap::new();
        for (chrom_idx, index) in result_rx {
            // Buffer out-of-order indexes
            if chrom_idx != next_chrom {
                buffered_indexes.insert(chrom_idx, index);
                continue;
            }
            match combined {
                None => {
                    combined = Some(index);
                    next_chrom += 1;
                }
                Some(ref mut c) => {
                    c.merge(index);
                    next_chrom += 1;
                }
            }
            while let Some(index) = buffered_indexes.remove(&next_chrom) {
                if let Some(ref mut c) = combined {
                    c.merge(index);
                    next_chrom += 1;
                }
            }
        }
        assert_eq!(
            next_chrom, num_chroms,
            "Did not receive all chromosome indexes"
        );
        assert_eq!(
            buffered_indexes.len(),
            0,
            "Buffered indexes not empty after merge"
        );
        // Wait for all workers to finish
        for handle in handles {
            handle.join().expect("worker thread panicked");
        }

        let combined = combined.unwrap_or_else(|| Index::<K, S>::new_with_chroms(&chrom_info));

        log::info!(
            "Index complete: {} unique seeds, {} nonunique seeds in {:.2}s",
            combined.unique_seeds.len(),
            combined.nonunique_seeds.len(),
            now.elapsed().as_secs_f64()
        );

        combined
    }
}

impl<const K: usize, const S: usize, R: std::io::BufRead> TryFrom<noodles::fasta::io::Reader<R>>
    for Index<K, S>
{
    type Error = crate::error::ParallaxError;

    fn try_from(mut reader: noodles::fasta::io::Reader<R>) -> Result<Self, Self::Error> {
        let mut index = Index::new();

        for record in reader.records() {
            let record = record?;
            let header = String::from_utf8(record.name().to_vec())
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let seq = record.sequence().as_ref().to_owned();
            index.add(header, seq);
        }

        if false {
            let mut hist = HashMap::new();
            hist.insert(1usize, index.unique_seeds.len());
            for locs in index.nonunique_seeds.values() {
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

        Ok(index)
    }
}
