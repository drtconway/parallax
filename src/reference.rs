//! Reference sequence access using indexed FASTA files.
//!
//! This module provides random access to reference sequences using
//! a FASTA file with an associated .fai index.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use niffler;
use noodles::core::Region;
use noodles::fasta;

use crate::error::Result;

/// Chromosome metadata parsed from FASTA headers.
///
/// GRCh38 FASTA headers contain key:value metadata, e.g.:
/// `>chr1_KI270762v1_alt  AC:KI270762.1  rg:chr1:2448811-2791270  rl:alt-scaffold`
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ChromInfo {
    /// Short name (e.g., "chr1", "chr1_KI270762v1_alt")
    pub name: String,
    /// Sequence length in bases (0 if unknown, e.g. from an older index)
    #[serde(default)]
    pub length: u64,
    /// Number of indexed syncmers for this sequence (0 if unknown, e.g. from an older index)
    #[serde(default)]
    pub syncmer_count: u64,
    /// Region localization type
    pub localization: Localization,
    /// Reference group - the primary chromosome this contig is associated with
    /// None for primary chromosomes and unplaced contigs
    pub reference_group: Option<String>,

    pub metadata: Vec<(String, String)>,
}

/// Chromosome localization type from the rl: metadata field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Localization {
    /// Primary chromosome (rl:Chromosome or rl:Mitochondrion)
    Chromosome,
    /// Unlocalized scaffold - associated with a chromosome but position unknown
    Unlocalized,
    /// ALT scaffold - alternative haplotype at a specific location
    AltScaffold,
    /// Unplaced - not associated with any chromosome
    Unplaced,
    /// Decoy sequences (e.g., chrUn_*_decoy)
    Decoy,
    /// Other sequences (e.g., HLA, EBV, etc.)
    Other,
    /// Unknown/missing localization
    Unknown,
}

impl ChromInfo {
    /// Parse chromosome info from a FASTA header.
    ///
    /// The header format is: `name  key:value  key:value  ...`
    pub fn from_header(name: &str, description: &str) -> Self {
        let mut metadata = Vec::new();
        for field in description.split_whitespace() {
            if let Some((key, value)) = field.split_once(':') {
                metadata.push((key.to_string(), value.to_string()));
            }
        }

        let mut localization = Localization::Unknown;
        let mut reference_group = None;

        if let Some(value) = metadata
            .iter()
            .find_map(|(k, v)| if k == "rl" { Some(v) } else { None })
        {
            localization = match value.as_str() {
                "Chromosome" => Localization::Chromosome,
                "Mitochondrion" => Localization::Chromosome, // chrM is considered primary
                "unlocalized" => Localization::Unlocalized,
                "alt-scaffold" => Localization::AltScaffold,
                "unplaced" => Localization::Unplaced,
                _ => Localization::Unknown,
            };
        }

        if let Some(value) = metadata
            .iter()
            .find_map(|(k, v)| if k == "rg" { Some(v) } else { None })
        {
            // rg can be "chr1" or "chr1:start-end"
            // Extract just the chromosome name
            let chrom = value.split(':').next().unwrap_or(value);
            reference_group = Some(chrom.to_string());
        }

        // Infer localization from name if not set
        if localization == Localization::Unknown {
            if name.starts_with("chrUn") {
                localization = Localization::Unplaced;
            } else if name.contains("_alt") {
                localization = Localization::AltScaffold;
            } else if name.contains("_random") {
                localization = Localization::Unlocalized;
            } else if name.contains("_decoy") {
                localization = Localization::Decoy;
            } else if name.starts_with("HLA-") {
                // HLA contigs (e.g., HLA-A*01:01:01:01) are separate sequences
                localization = Localization::Other;
            } else if name.starts_with("chr") && !name.contains('_') {
                // Only chr1, chr2, ..., chrX, chrY, chrM without underscores are primary
                localization = Localization::Chromosome;
            } else {
                // Anything else is unknown/other
                localization = Localization::Other;
            }
        }

        ChromInfo {
            name: name.to_string(),
            length: 0,
            syncmer_count: 0,
            localization,
            reference_group,
            metadata,
        }
    }

    /// Returns true if this is a primary chromosome.
    pub fn is_primary(&self) -> bool {
        self.localization == Localization::Chromosome
    }

    /// Returns true if this is an ALT scaffold.
    #[allow(dead_code)]
    pub fn is_alt(&self) -> bool {
        self.localization == Localization::AltScaffold
    }

    /// Returns the primary chromosome this contig is associated with.
    /// For primary chromosomes, returns their own name.
    #[allow(dead_code)]
    pub fn primary_chrom(&self) -> &str {
        self.reference_group.as_deref().unwrap_or(&self.name)
    }
}

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
    /// Chromosome metadata in order
    chrom_info: Vec<ChromInfo>,
    /// Chromosome sequences (uppercase)
    sequences: Vec<Vec<u8>>,
}

impl InMemoryReference {
    /// Load a reference FASTA file entirely into memory.
    ///
    /// All sequences are converted to uppercase during loading.
    /// If `primary_only` is true, only primary chromosomes are loaded
    /// (ALT scaffolds, unlocalized, and unplaced contigs are skipped).
    pub fn load<P: AsRef<Path>>(fasta_path: P, primary_only: bool) -> Result<Self> {
        let fasta_path = fasta_path.as_ref();

        log::info!(
            "Loading reference from {} into memory{}",
            fasta_path.display(),
            if primary_only {
                " (primary contigs only)"
            } else {
                ""
            }
        );

        let (reader, fmt) = niffler::from_path(fasta_path)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        if fmt != niffler::Format::No {
            log::info!("Detected {:?} compression on reference input", fmt);
        }
        let mut reader = fasta::io::Reader::new(BufReader::new(reader));

        let mut chrom_info = Vec::new();
        let mut sequences = Vec::new();

        for result in reader.records() {
            let record = result?;
            let name = String::from_utf8_lossy(record.name()).into_owned();
            let description = record
                .description()
                .map(|d| String::from_utf8_lossy(d).into_owned());

            let mut info = ChromInfo::from_header(&name, description.as_deref().unwrap_or(""));

            // Skip non-primary contigs if primary_only is set
            if primary_only && !info.is_primary() {
                log::debug!(
                    "Skipping non-primary contig {} ({:?})",
                    name,
                    info.localization
                );
                continue;
            }

            // Convert sequence to uppercase
            let seq: Vec<u8> = record
                .sequence()
                .as_ref()
                .iter()
                .map(|&b| b.to_ascii_uppercase())
                .collect();

            info.length = seq.len() as u64;

            log::debug!(
                "Loaded chromosome {} ({} bp, {:?})",
                name,
                seq.len(),
                info.localization
            );
            chrom_info.push(info);
            sequences.push(seq);
        }

        log::info!(
            "Loaded {} chromosome(s), total {} bp",
            chrom_info.len(),
            sequences.iter().map(|s| s.len()).sum::<usize>()
        );

        Ok(Self {
            chrom_info,
            sequences,
        })
    }

    /// Get the chromosome name by index.
    pub fn chrom_name(&self, chrom_idx: usize) -> &str {
        &self.chrom_info[chrom_idx].name
    }

    /// Get the chromosome info by index.
    pub fn chrom_info(&self, chrom_idx: usize) -> &ChromInfo {
        &self.chrom_info[chrom_idx]
    }

    /// Get all chromosome info.
    #[allow(dead_code)]
    pub fn all_chrom_info(&self) -> &[ChromInfo] {
        &self.chrom_info
    }

    /// Get the number of chromosomes.
    pub fn num_chroms(&self) -> usize {
        self.chrom_info.len()
    }

    /// Get the chromosome length by index.
    pub fn chrom_length(&self, chrom_idx: usize) -> u64 {
        self.sequences[chrom_idx].len() as u64
    }

    /// Iterate over all chromosomes, yielding (name, length) pairs.
    /// Useful for generating SAM headers.
    pub fn chromosomes(&self) -> impl Iterator<Item = (&str, u64)> {
        self.chrom_info
            .iter()
            .zip(self.sequences.iter())
            .map(|(info, s)| (info.name.as_str(), s.len() as u64))
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

    /// Build a noodles FASTA repository backed by the in-memory sequences.
    ///
    /// The returned `Repository` is used by the CRAM writer to resolve
    /// reference sequences during encoding.
    pub fn to_fasta_repository(&self) -> fasta::Repository {
        use std::collections::HashMap;
        let map: HashMap<Vec<u8>, Vec<u8>> = self
            .chrom_info
            .iter()
            .zip(self.sequences.iter())
            .map(|(info, seq)| (info.name.as_bytes().to_vec(), seq.clone()))
            .collect();
        fasta::Repository::new(InMemoryAdapter(map))
    }
}

/// Adapter implementing `fasta::repository::Adapter` for O(1) lookups
/// from an in-memory map of chromosome name → sequence.
struct InMemoryAdapter(std::collections::HashMap<Vec<u8>, Vec<u8>>);

impl fasta::repository::Adapter for InMemoryAdapter {
    fn get(&mut self, name: &[u8]) -> Option<std::io::Result<fasta::Record>> {
        let seq = self.0.get(name)?;
        let record = fasta::Record::new(
            fasta::record::Definition::new(std::str::from_utf8(name).unwrap_or("?"), None),
            fasta::record::Sequence::from(seq.clone()),
        );
        Some(Ok(record))
    }
}
