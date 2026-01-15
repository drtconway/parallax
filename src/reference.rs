//! Reference sequence access using indexed FASTA files.
//!
//! This module provides random access to reference sequences using
//! a FASTA file with an associated .fai index.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::core::Region;
use noodles::fasta;

use crate::error::Result;

/// Reference sequence reader using indexed FASTA.
///
/// Provides random access to reference sequences without
/// loading them entirely into memory.
pub struct Reference {
    reader: fasta::io::IndexedReader<BufReader<File>>,
    /// Chromosome names in order (matching Index order)
    chrom_names: Vec<String>,
    /// Chromosome lengths in order
    chrom_lengths: Vec<u64>,
}

impl Reference {
    /// Open a reference FASTA file with its .fai index.
    ///
    /// The .fai file must exist at `{fasta_path}.fai`.
    #[allow(dead_code)]
    pub fn open<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let fasta_path = fasta_path.as_ref();
        let fai_path = fasta_path.with_extension("fa.fai");

        // Try .fa.fai first, then .fasta.fai, then just .fai appended
        let fai_path = if fai_path.exists() {
            fai_path
        } else {
            let fai_path = fasta_path.with_extension("fasta.fai");
            if fai_path.exists() {
                fai_path
            } else {
                // Just append .fai to the original path
                let mut p = fasta_path.as_os_str().to_owned();
                p.push(".fai");
                std::path::PathBuf::from(p)
            }
        };

        log::info!("Loading FASTA index from {}", fai_path.display());
        let index = fasta::fai::fs::read(&fai_path)?;

        // Extract chromosome names and lengths from the index
        let chrom_names: Vec<String> = index
            .as_ref()
            .iter()
            .map(|record| String::from_utf8_lossy(record.name()).into_owned())
            .collect();

        let chrom_lengths: Vec<u64> = index
            .as_ref()
            .iter()
            .map(|record| record.length())
            .collect();

        log::info!("Loaded {} chromosome(s) from index", chrom_names.len());

        let file = File::open(fasta_path)?;
        let reader = fasta::io::IndexedReader::new(BufReader::new(file), index);

        Ok(Self {
            reader,
            chrom_names,
            chrom_lengths,
        })
    }

    /// Get the chromosome name by index.
    #[allow(dead_code)]
    pub fn chrom_name(&self, chrom_idx: usize) -> &str {
        &self.chrom_names[chrom_idx]
    }

    /// Get the number of chromosomes.
    #[allow(dead_code)]
    pub fn num_chroms(&self) -> usize {
        self.chrom_names.len()
    }

    /// Get the chromosome length by index.
    #[allow(dead_code)]
    pub fn chrom_length(&self, chrom_idx: usize) -> u64 {
        self.chrom_lengths[chrom_idx]
    }

    /// Iterate over all chromosomes, yielding (name, length) pairs.
    /// Useful for generating SAM headers.
    #[allow(dead_code)]
    pub fn chromosomes(&self) -> impl Iterator<Item = (&str, u64)> {
        self.chrom_names
            .iter()
            .zip(self.chrom_lengths.iter())
            .map(|(n, &l)| (n.as_str(), l))
    }

    /// Fetch a sequence region into the provided buffer.
    ///
    /// Uses 0-based, half-open coordinates [start, end).
    #[allow(dead_code)]
    pub fn get_seq_into(
        &mut self,
        chrom_idx: usize,
        start: usize,
        end: usize,
        buf: &mut Vec<u8>,
    ) -> Result<()> {
        if start >= end {
            buf.clear();
            return Ok(());
        }

        let chrom_name = &self.chrom_names[chrom_idx];

        // noodles uses 1-based, closed interval [start, end]
        let region: Region = format!("{}:{}-{}", chrom_name, start + 1, end)
            .parse()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("{}", e)))?;

        let record = self.reader.query(&region)?;

        buf.clear();
        buf.extend_from_slice(record.sequence().as_ref());

        Ok(())
    }
}

/// Reference sequence reader that loads entire sequences into memory.
///
/// Unlike `Reference`, this loads all chromosome sequences into memory
/// at construction time and converts lowercase bases to uppercase.
/// This avoids repeated file I/O and case conversion during alignment.
#[derive(Clone)]
pub struct InMemoryReference {
    /// Chromosome names in order
    chrom_names: Vec<String>,
    /// Chromosome sequences (uppercase)
    sequences: Vec<Vec<u8>>,
}

impl InMemoryReference {
    /// Load a reference FASTA file entirely into memory.
    ///
    /// All sequences are converted to uppercase during loading.
    pub fn load<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let fasta_path = fasta_path.as_ref();

        log::info!(
            "Loading reference from {} into memory",
            fasta_path.display()
        );

        let file = File::open(fasta_path)?;
        let mut reader = fasta::io::Reader::new(BufReader::new(file));

        let mut chrom_names = Vec::new();
        let mut sequences = Vec::new();

        for result in reader.records() {
            let record = result?;
            let name = String::from_utf8_lossy(record.name()).into_owned();

            // Convert sequence to uppercase
            let seq: Vec<u8> = record
                .sequence()
                .as_ref()
                .iter()
                .map(|&b| b.to_ascii_uppercase())
                .collect();

            log::debug!("Loaded chromosome {} ({} bp)", name, seq.len());
            chrom_names.push(name);
            sequences.push(seq);
        }

        log::info!(
            "Loaded {} chromosome(s), total {} bp",
            chrom_names.len(),
            sequences.iter().map(|s| s.len()).sum::<usize>()
        );

        Ok(Self {
            chrom_names,
            sequences,
        })
    }

    /// Get the chromosome name by index.
    pub fn chrom_name(&self, chrom_idx: usize) -> &str {
        &self.chrom_names[chrom_idx]
    }

    /// Get the number of chromosomes.
    pub fn num_chroms(&self) -> usize {
        self.chrom_names.len()
    }

    /// Get the chromosome length by index.
    pub fn chrom_length(&self, chrom_idx: usize) -> u64 {
        self.sequences[chrom_idx].len() as u64
    }

    /// Iterate over all chromosomes, yielding (name, length) pairs.
    /// Useful for generating SAM headers.
    pub fn chromosomes(&self) -> impl Iterator<Item = (&str, u64)> {
        self.chrom_names
            .iter()
            .zip(self.sequences.iter())
            .map(|(n, s)| (n.as_str(), s.len() as u64))
    }

    /// Get the full sequence for a chromosome.
    ///
    /// Returns a slice to the in-memory sequence (no copying).
    #[inline]
    pub fn sequence(&self, chrom_idx: usize) -> &[u8] {
        &self.sequences[chrom_idx]
    }

    /// Get a subsequence for a chromosome.
    ///
    /// Uses 0-based, half-open coordinates [start, end).
    /// Returns a slice to the in-memory sequence (no copying).
    #[inline]
    pub fn get_seq(&self, chrom_idx: usize, start: usize, end: usize) -> &[u8] {
        let seq = &self.sequences[chrom_idx];
        let end = end.min(seq.len());
        let start = start.min(end);
        &seq[start..end]
    }
}
