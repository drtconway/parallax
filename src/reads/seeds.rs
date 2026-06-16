use serde::{Deserialize, Serialize};

/// A seed hit representing a k-mer match between read and reference
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct SeedHit {
    /// Chromosome/contig index in the reference
    pub chrom_id: usize,
    /// Diagonal: ref_pos - read_pos (constant for colinear matches)
    pub diagonal: i64,
    /// Position in the reference sequence (always 5' end on forward strand)
    pub ref_pos: usize,
    /// Position in the read sequence (in strand-local coordinates)
    pub read_pos: usize,
    /// Initial kmer
    pub kmer: u64,
    /// Uniqueness of the most unique k-mer incorporated in this seed (lower is more unique)
    pub kmer_uniqueness: u32,
    /// Number of times this k-mer appears in the read (1 = unique within read)
    pub read_frequency: u32,
    /// Length of the match (initially k, may be extended)
    pub match_len: usize,
}

impl SeedHit {
    /// Create a new seed hit with read_frequency defaulting to 1.
    pub fn new(
        chrom_id: usize,
        ref_pos: usize,
        read_pos: usize,
        kmer: u64,
        kmer_uniqueness: u32,
        match_len: usize,
    ) -> Self {
        Self {
            chrom_id,
            diagonal: ref_pos as i64 - read_pos as i64,
            ref_pos,
            read_pos,
            kmer,
            kmer_uniqueness,
            read_frequency: 1,
            match_len,
        }
    }

    /// Create a new seed hit with an explicit read_frequency.
    pub fn with_read_frequency(
        chrom_id: usize,
        ref_pos: usize,
        read_pos: usize,
        kmer: u64,
        kmer_uniqueness: u32,
        read_frequency: u32,
        match_len: usize,
    ) -> Self {
        Self {
            chrom_id,
            diagonal: ref_pos as i64 - read_pos as i64,
            ref_pos,
            read_pos,
            kmer,
            kmer_uniqueness,
            read_frequency,
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
}

impl SeedHit {
    /// Attempt to extend the seed hit if the new k-mer overlaps the current match
    /// or return a new seed hit if the new k-mer does not overlap.
    ///
    /// Two k-mers can be merged only if they overlap - i.e., the new k-mer starts
    /// within the current match region. This ensures no unverified gaps exist.
    pub fn extend(
        &mut self,
        chrom_id: usize,
        chrom_pos: usize,
        read_pos: usize,
        kmer: u64,
        kmer_uniqueness: u32,
        read_frequency: u32,
        k: usize,
    ) -> Option<SeedHit> {
        if chrom_id == self.chrom_id
            && chrom_pos >= self.ref_pos
            && read_pos >= self.read_pos
            && chrom_pos - self.ref_pos == read_pos - self.read_pos
            // Only merge if the new k-mer starts within our current match region
            // (i.e., it overlaps). This prevents creating gaps with unverified bases.
            && chrom_pos - self.ref_pos <= self.match_len
        {
            // Overlaps or extends current match
            let new_end = (chrom_pos - self.ref_pos) + k;
            if new_end > self.match_len {
                self.match_len = new_end;
            }
            if kmer_uniqueness < self.kmer_uniqueness {
                self.kmer_uniqueness = kmer_uniqueness;
            }
            if read_frequency > self.read_frequency {
                self.read_frequency = read_frequency;
            }
            None
        } else {
            // Does not overlap - return new seed hit
            Some(SeedHit::with_read_frequency(
                chrom_id,
                chrom_pos,
                read_pos,
                kmer,
                kmer_uniqueness,
                read_frequency,
                k,
            ))
        }
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

#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Read {
    pub name: String,
    pub seq: String,
    pub qual: String,
}

impl Read {
    #[allow(dead_code)]
    pub fn new(name: &str, seq: &[u8], qual: &[u8]) -> Self {
        Self {
            name: name.to_string(),
            seq: String::from_utf8_lossy(seq).to_string(),
            qual: String::from_utf8_lossy(qual).to_string(),
        }
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedSaver {
    pub read: Read,
    pub is_reverse: bool,
    pub seeds: Vec<SeedHit>,
}
