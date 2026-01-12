use std::collections::HashMap;

use crate::{
    kmers::Kmer,
    utils::{Selection, table::Table},
};

pub struct Index<const K: usize, const S: usize> {
    chroms: Vec<String>,
    sequences: Vec<Vec<u8>>,
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
            sequences: Vec::new(),
            cumulative_lengths: vec![0],
            unique_seeds: Table::new(),
            nonunique_seeds: HashMap::new(),
        }
    }

    pub fn add<Seq: AsRef<[u8]>>(&mut self, chrom: String, seq: Seq) -> usize {
        let idx = self.chroms.len();
        self.chroms.push(chrom);
        self.sequences.push(seq.as_ref().to_vec());

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

    /// Get a slice of the reference sequence for a given chromosome
    pub fn get_seq(&self, chrom_idx: usize, start: usize, end: usize) -> &[u8] {
        let seq = &self.sequences[chrom_idx];
        let end = end.min(seq.len());
        &seq[start..end]
    }

    /// Get the chromosome name
    pub fn chrom_name(&self, chrom_idx: usize) -> &str {
        &self.chroms[chrom_idx]
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
