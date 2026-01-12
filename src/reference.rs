//! Reference sequence access using indexed FASTA files.
//!
//! This module provides random access to reference sequences using
//! a FASTA file with an associated .fai index.

use std::io::BufReader;
use std::fs::File;
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
    /// Buffer for fetched sequences (reused to avoid allocations)
    seq_buf: Vec<u8>,
}

impl Reference {
    /// Open a reference FASTA file with its .fai index.
    /// 
    /// The .fai file must exist at `{fasta_path}.fai`.
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
        
        // Extract chromosome names from the index
        let chrom_names: Vec<String> = index
            .as_ref()
            .iter()
            .map(|record| String::from_utf8_lossy(record.name()).into_owned())
            .collect();
        
        log::info!("Loaded {} chromosome(s) from index", chrom_names.len());
        
        let file = File::open(fasta_path)?;
        let reader = fasta::io::IndexedReader::new(BufReader::new(file), index);
        
        Ok(Self {
            reader,
            chrom_names,
            seq_buf: Vec::new(),
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

    /// Fetch a sequence region into the internal buffer and return a reference.
    /// 
    /// Uses 0-based, half-open coordinates [start, end).
    /// The returned slice is valid until the next call to `get_seq`.
    pub fn get_seq(&mut self, chrom_idx: usize, start: usize, end: usize) -> Result<&[u8]> {
        if start >= end {
            self.seq_buf.clear();
            return Ok(&self.seq_buf);
        }

        let chrom_name = &self.chrom_names[chrom_idx];
        
        // noodles uses 1-based, closed interval [start, end]
        let region: Region = format!("{}:{}-{}", chrom_name, start + 1, end)
            .parse()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("{}", e)))?;
        
        let record = self.reader.query(&region)?;
        
        self.seq_buf.clear();
        self.seq_buf.extend_from_slice(record.sequence().as_ref());
        
        Ok(&self.seq_buf)
    }

    /// Fetch a sequence region into the provided buffer.
    /// 
    /// Uses 0-based, half-open coordinates [start, end).
    pub fn get_seq_into(&mut self, chrom_idx: usize, start: usize, end: usize, buf: &mut Vec<u8>) -> Result<()> {
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
