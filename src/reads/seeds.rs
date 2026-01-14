use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::{Mutex, OnceLock},
};

use crate::{error::Result, writer::AlignmentWriter};

/// Global debug SAM file for writing seed hits
static DEBUG_SAM_FILE: OnceLock<Mutex<BufWriter<File>>> = OnceLock::new();

/// Initialize the debug SAM file. Call once at startup if debugging is needed.
/// Returns Ok(()) if successful or if already initialized.
#[allow(dead_code)]
pub fn init_debug_sam(path: &str) -> std::io::Result<()> {
    DEBUG_SAM_FILE.get_or_init(|| {
        let file = File::create(path).expect("Failed to create debug SAM file");
        Mutex::new(BufWriter::new(file))
    });
    Ok(())
}

/// Write a line to the debug SAM file if it's been initialized.
/// Silently does nothing if debug SAM was not initialized.
#[allow(dead_code)]
pub fn write_debug_sam(line: &str) {
    if let Some(mutex) = DEBUG_SAM_FILE.get() {
        if let Ok(mut writer) = mutex.lock() {
            let _ = writeln!(writer, "{}", line);
        }
    }
}

/// Flush the debug SAM file if it's been initialized.
#[allow(dead_code)]
pub fn flush_debug_sam() {
    if let Some(mutex) = DEBUG_SAM_FILE.get() {
        if let Ok(mut writer) = mutex.lock() {
            let _ = writer.flush();
        }
    }
}

/// A seed hit representing a k-mer match between read and reference
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeedHit {
    /// Chromosome/contig index in the reference
    pub chrom_id: usize,
    /// Diagonal: ref_pos - read_pos (constant for colinear matches)
    pub diagonal: i64,
    /// Position in the reference sequence
    pub ref_pos: usize,
    /// Position in the read sequence
    pub read_pos: usize,
    /// Initial kmer
    pub kmer: u64,
    /// Length of the match (initially k, may be extended)
    pub match_len: usize,
}

impl SeedHit {
    /// Create a new seed hit
    pub fn new(
        chrom_id: usize,
        ref_pos: usize,
        read_pos: usize,
        kmer: u64,
        match_len: usize,
    ) -> Self {
        Self {
            chrom_id,
            diagonal: ref_pos as i64 - read_pos as i64,
            ref_pos,
            read_pos,
            kmer,
            match_len,
        }
    }

    /// End position in the read
    pub fn read_end(&self) -> usize {
        self.read_pos + self.match_len
    }

    /// End position in the reference
    pub fn ref_end(&self) -> usize {
        self.ref_pos + self.match_len
    }

    /// Attempt to the seed hit if the new k-mer extends the current match
    /// or return a new seed hit if the new k-mer does not overlap.
    pub fn extend(
        &mut self,
        chrom_id: usize,
        chrom_pos: usize,
        read_pos: usize,
        kmer: u64,
        k: usize,
    ) -> Option<SeedHit> {
        if chrom_id == self.chrom_id
            && chrom_pos >= self.ref_pos
            && read_pos >= self.read_pos
            && chrom_pos - self.ref_pos == read_pos - self.read_pos
            && chrom_pos - self.ref_pos < self.match_len + k
        {
            // Overlaps or extends current match
            let new_end = (chrom_pos - self.ref_pos) + k;
            if new_end > self.match_len {
                self.match_len = new_end;
            }
            None
        } else {
            // Does not overlap - return new seed hit
            Some(SeedHit::new(chrom_id, chrom_pos, read_pos, kmer, k))
        }
    }

    #[allow(dead_code)]
    pub fn write_sam<W: std::io::Write>(
        &self,
        writer: &mut AlignmentWriter<W>,
        read_id: &str,
        is_reverse: bool,
    ) -> Result<()> {
        let flag = if is_reverse {
            super::FLAG_REVERSE
        } else {
            0u16
        };
        writer.write_alignment(
            read_id,
            flag,
            &format!("chr{}", self.chrom_id),
            self.ref_pos + 1, // 1-based
            255,              // placeholder MAPQ
            &format!("{}M", self.match_len),
            "*",
            0,
            0,
            "*",
            "*",
            &[],
        )?;
        Ok(())
    }

    /// Format as SAM line string for debug output
    #[allow(dead_code)]
    pub fn to_sam_line(&self, read_id: &str, is_reverse: bool) -> String {
        let flag = if is_reverse {
            super::FLAG_REVERSE
        } else {
            0u16
        };
        format!(
            "{}\t{}\tchr{}\t{}\t255\t{}M\t*\t0\t0\t*\t*",
            read_id,
            flag,
            self.chrom_id,
            self.ref_pos + 1,
            self.match_len
        )
    }

    /// Extend the seed match bidirectionally as far as sequences match exactly.
    /// This is the minimap2-style seed extension that maximizes anchor length.
    /// 
    /// Returns the number of bases extended (backward + forward).
    pub fn extend_exact(&mut self, read_seq: &[u8], ref_seq: &[u8]) -> usize {
        let mut extended = 0;

        // Extend backward
        let mut back_ext = 0;
        while self.read_pos > back_ext
            && self.ref_pos > back_ext
            && read_seq[self.read_pos - back_ext - 1] == ref_seq[self.ref_pos - back_ext - 1]
        {
            back_ext += 1;
        }
        if back_ext > 0 {
            self.read_pos -= back_ext;
            self.ref_pos -= back_ext;
            self.match_len += back_ext;
            self.diagonal = self.ref_pos as i64 - self.read_pos as i64;
            extended += back_ext;
        }

        // Extend forward
        let mut fwd_ext = 0;
        while self.ref_pos + self.match_len + fwd_ext < ref_seq.len()
            && self.read_pos + self.match_len + fwd_ext < read_seq.len()
            && read_seq[self.read_pos + self.match_len + fwd_ext]
                == ref_seq[self.ref_pos + self.match_len + fwd_ext]
        {
            fwd_ext += 1;
        }
        if fwd_ext > 0 {
            self.match_len += fwd_ext;
            extended += fwd_ext;
        }

        extended
    }
}
