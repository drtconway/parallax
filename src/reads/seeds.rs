use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::{Mutex, OnceLock},
};

use crate::{error::Result, writer::AlignmentWriter};

/// Global debug SAM file for writing seed hits
static DEBUG_SAM_FILE: OnceLock<Mutex<BufWriter<File>>> = OnceLock::new();

/// Initialize the debug SAM file with a proper SAM header.
/// Call once at startup if debugging is needed.
/// 
/// # Arguments
/// * `path` - Path to the output SAM file
/// * `chromosomes` - Iterator of (name, length) pairs for @SQ headers
pub fn init_debug_sam<'a>(
    path: &str,
    chromosomes: impl Iterator<Item = (&'a str, u64)>,
) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    
    // Write SAM header
    writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")?;
    for (name, length) in chromosomes {
        writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, length)?;
    }
    writeln!(writer, "@PG\tID:parallax\tPN:parallax\tVN:0.1.0\tCL:debug seeds")?;
    
    DEBUG_SAM_FILE.get_or_init(|| Mutex::new(writer));
    Ok(())
}

/// Write a line to the debug SAM file if it's been initialized.
/// Silently does nothing if debug SAM was not initialized.
pub fn write_debug_sam(line: &str) {
    if let Some(mutex) = DEBUG_SAM_FILE.get() {
        if let Ok(mut writer) = mutex.lock() {
            let _ = writeln!(writer, "{}", line);
        }
    }
}

/// Flush the debug SAM file if it's been initialized.
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

    /// Validate that this seed is an exact match between read and reference.
    /// Returns true if all bases match exactly, false otherwise.
    /// If there are mismatches, logs them for debugging.
    #[allow(dead_code)]
    pub fn validate_exact_match(&self, read_seq: &[u8], ref_seq: &[u8], context: &str) -> bool {
        let read_end = self.read_pos + self.match_len;
        let ref_end = self.ref_pos + self.match_len;

        // Bounds check
        if read_end > read_seq.len() {
            log::warn!(
                "[{}] Seed read range [{}, {}) exceeds read length {}",
                context,
                self.read_pos,
                read_end,
                read_seq.len()
            );
            return false;
        }
        if ref_end > ref_seq.len() {
            log::warn!(
                "[{}] Seed ref range [{}, {}) exceeds ref length {}",
                context,
                self.ref_pos,
                ref_end,
                ref_seq.len()
            );
            return false;
        }

        let read_slice = &read_seq[self.read_pos..read_end];
        let ref_slice = &ref_seq[self.ref_pos..ref_end];

        let mut mismatches = Vec::new();
        for (i, (&r, &q)) in ref_slice.iter().zip(read_slice.iter()).enumerate() {
            if r != q {
                mismatches.push((i, r as char, q as char));
            }
        }

        if !mismatches.is_empty() {
            log::warn!(
                "[{}] Seed at ref_pos={} read_pos={} match_len={} has {} mismatches:",
                context,
                self.ref_pos,
                self.read_pos,
                self.match_len,
                mismatches.len()
            );
            for (i, ref_base, read_base) in mismatches.iter().take(10) {
                log::warn!("  Position {}: ref={} read={}", i, ref_base, read_base);
            }
            return false;
        }

        true
    }

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

    /// Format as SAM line string for debug output with proper hard clips.
    ///
    /// # Arguments
    /// * `read_id` - Read name
    /// * `chrom_name` - Chromosome name (not ID)
    /// * `is_reverse` - Whether this seed is on the reverse strand
    /// * `strand_seq` - The sequence for this strand (already rev-comped if reverse)
    /// * `strand_qual` - Quality scores for this strand (already reversed if reverse), or None
    pub fn to_sam_line(
        &self,
        read_id: &str,
        chrom_name: &str,
        is_reverse: bool,
        strand_seq: &[u8],
        strand_qual: Option<&[u8]>,
    ) -> String {
        let read_len = strand_seq.len();
        let flag = if is_reverse {
            super::FLAG_REVERSE
        } else {
            0u16
        };

        // Build CIGAR with hard clips
        // For forward strand: read_pos bases before, match_len bases aligned, rest after
        // For reverse strand: coordinates are already in forward read space after revcomp
        let hclip_start = self.read_pos;
        let hclip_end = read_len.saturating_sub(self.read_pos + self.match_len);

        let cigar = match (hclip_start > 0, hclip_end > 0) {
            (true, true) => format!("{}H{}={}H", hclip_start, self.match_len, hclip_end),
            (true, false) => format!("{}H{}=", hclip_start, self.match_len),
            (false, true) => format!("{}={}H", self.match_len, hclip_end),
            (false, false) => format!("{}=", self.match_len),
        };

        // Extract the aligned portion of the sequence
        let seq_slice = &strand_seq[self.read_pos..self.read_pos + self.match_len];
        let seq_str = String::from_utf8_lossy(seq_slice);

        // Extract the aligned portion of the quality scores, or use * if not available
        let qual_str = match strand_qual {
            Some(qual) => {
                let qual_slice = &qual[self.read_pos..self.read_pos + self.match_len];
                // Convert Phred+33 quality scores to ASCII
                qual_slice.iter().map(|&q| q as char).collect::<String>()
            }
            None => "*".to_string(),
        };

        format!(
            "{}\t{}\t{}\t{}\t255\t{}\t*\t0\t0\t{}\t{}",
            read_id,
            flag,
            chrom_name,
            self.ref_pos + 1, // 1-based
            cigar,
            seq_str,
            qual_str
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

/// A cluster of seed hits with its LIS chain, before alignment building.
/// This intermediate structure allows for cross-strand gap analysis before
/// committing to alignment construction.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SeedCluster {
    /// Read region covered by this cluster
    pub read_start: usize,
    pub read_end: usize,

    /// The LIS chain of colinear seeds (sorted by read position)
    pub chain: Vec<SeedHit>,

    /// Which strand this cluster came from
    pub is_reverse: bool,

    /// Chromosome this cluster aligns to
    pub chrom_id: usize,
}

impl SeedCluster {
    /// Create a new cluster from a chain of seeds
    pub fn new(mut chain: Vec<SeedHit>, is_reverse: bool) -> Option<Self> {
        if chain.is_empty() {
            return None;
        }

        // Sort by read position for alignment building
        chain.sort_by_key(|hit| hit.read_pos);

        let read_start = chain.first().map(|h| h.read_pos).unwrap_or(0);
        let read_end = chain.last().map(|h| h.read_end()).unwrap_or(0);
        let chrom_id = chain[0].chrom_id;

        Some(SeedCluster {
            chrom_id,
            read_start,
            read_end,
            chain,
            is_reverse,
        })
    }

    pub fn fwd_read_range(&self, read_len: usize) -> (usize, usize) {
        if self.is_reverse {
            (read_len - self.read_end, read_len - self.read_start)
        } else {
            (self.read_start, self.read_end)
        }
    }

    /// Calculate fraction of read covered by this cluster
    pub fn read_coverage(&self, read_len: usize) -> f64 {
        if read_len == 0 {
            return 0.0;
        }
        (self.read_end - self.read_start) as f64 / read_len as f64
    }

    /// Total seed match length in this cluster
    pub fn total_seed_length(&self) -> usize {
        self.chain.iter().map(|h| h.match_len).sum()
    }

    /// Split this cluster at the specified gap, returning the new cluster
    /// containing seeds after the gap.
    ///
    /// `gap_seed_idx` is the index of the seed before the gap (as returned by `gaps()`).
    /// After splitting:
    /// - `self` contains seeds 0..=gap_seed_idx
    /// - The returned cluster contains seeds gap_seed_idx+1..
    ///
    /// Returns `None` if the split would leave either cluster empty.
    pub fn split_at_gap(&mut self, gap_seed_idx: usize) -> Option<SeedCluster> {
        // Need at least one seed on each side
        if gap_seed_idx + 1 >= self.chain.len() {
            return None;
        }

        // Split the chain
        let tail_chain = self.chain.split_off(gap_seed_idx + 1);

        // Update self's read_end
        self.read_end = self.chain.last().map(|h| h.read_end()).unwrap_or(self.read_start);

        // Build the new cluster from the tail
        let tail_read_start = tail_chain.first().map(|h| h.read_pos).unwrap_or(0);
        let tail_read_end = tail_chain.last().map(|h| h.read_end()).unwrap_or(0);

        Some(SeedCluster {
            read_start: tail_read_start,
            read_end: tail_read_end,
            chain: tail_chain,
            is_reverse: self.is_reverse,
            chrom_id: self.chrom_id,
        })
    }

    /// Fraction of the chain covered by seeds
    pub fn seed_density(&self) -> f64 {
        let total_length = if let Some(first) = self.chain.first() {
            if let Some(last) = self.chain.last() {
                last.ref_end() - first.ref_pos
            } else {
                0
            }
        } else {
            0
        };
        if total_length == 0 {
            0.0
        } else {
            self.total_seed_length() as f64 / total_length as f64
        }
    }

    /// Return an iterator over gaps in the seed coverage that are
    /// at least as big as the given size. Each gap is (start, end) in forward-strand
    /// read coordinates, along with the index of the seed before the gap.
    pub fn gaps(
        &self,
        read_len: usize,
        min_gap_length: usize,
    ) -> impl Iterator<Item = ((usize, usize), usize)> + '_ {
        self.chain
            .windows(2)
            .enumerate()
            .filter_map(move |(i, pair)| {
                let gap_start = pair[0].read_end();
                let gap_end = pair[1].read_pos;
                if gap_end > gap_start && gap_end - gap_start >= min_gap_length {
                    let fwd_read_start = if self.is_reverse {
                        read_len - gap_end
                    } else {
                        gap_start
                    };
                    let fwd_read_end = if self.is_reverse {
                        read_len - gap_start
                    } else {
                        gap_end
                    };
                    Some(((fwd_read_start, fwd_read_end), i))
                } else {
                    None
                }
            })
    }
}

/// A gap in one cluster that is filled by another cluster
#[derive(Debug, Clone)]
pub struct GapFill {
    /// Index of the cluster containing the gap
    pub cluster_idx: usize,
    /// Gap range in forward-strand read coordinates
    #[allow(dead_code)]
    pub gap: (usize, usize),
    /// Index of the seed before this gap in the cluster's chain
    pub gap_seed_idx: usize,
    /// Index of the cluster that fills this gap
    #[allow(dead_code)]
    pub filler_idx: usize,
    /// Fraction of the gap covered by the filler
    #[allow(dead_code)]
    pub coverage: f64,
}

/// Pre-sorted cluster index for efficient gap-fill queries
pub struct ClusterIndex {
    /// Clusters sorted by read_start: (read_start, read_end, original_idx)
    by_start: Vec<(usize, usize, usize)>,
}

impl ClusterIndex {
    /// Build an index from clusters, converting all to forward-strand coordinates
    pub fn new(clusters: &[SeedCluster], read_len: usize) -> Self {
        let mut by_start: Vec<_> = clusters
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let (start, end) = c.fwd_read_range(read_len);
                (start, end, i)
            })
            .collect();

        by_start.sort_unstable_by_key(|&(start, _, _)| start);

        ClusterIndex { by_start }
    }

    /// Find clusters that overlap [gap_start - tolerance, gap_end + tolerance]
    /// Returns iterator of (original_idx, overlap_fraction)
    pub fn find_overlapping(
        &self,
        gap_start: usize,
        gap_end: usize,
        tolerance: usize,
        exclude_idx: usize,
    ) -> impl Iterator<Item = (usize, f64)> + '_ {
        let query_start = gap_start.saturating_sub(tolerance);
        let query_end = gap_end.saturating_add(tolerance);
        let gap_len = gap_end - gap_start;

        // Binary search: first cluster where read_end > query_start
        // (clusters ending before our query can't overlap)
        let first = self
            .by_start
            .partition_point(|&(_, end, _)| end <= query_start);

        self.by_start[first..]
            .iter()
            .take_while(move |&&(start, _, _)| start < query_end)
            .filter(move |&&(_, _, idx)| idx != exclude_idx)
            .map(move |&(start, end, idx)| {
                let overlap_start = gap_start.max(start.saturating_sub(tolerance));
                let overlap_end = gap_end.min(end.saturating_add(tolerance));
                let overlap = overlap_end.saturating_sub(overlap_start);
                let fraction = overlap as f64 / gap_len as f64;
                (idx, fraction)
            })
    }
}

/// Analyze all gaps across all clusters to find chimeric breakpoints
///
/// Returns a list of gaps that are filled by other clusters, indicating
/// potential chimeric read breakpoints.
pub fn analyze_gap_fills(
    clusters: &[SeedCluster],
    read_len: usize,
    min_gap: usize,
    tolerance: usize,
    min_coverage: f64,
) -> Vec<GapFill> {
    if clusters.len() < 2 {
        return Vec::new();
    }

    let index = ClusterIndex::new(clusters, read_len);
    let mut fills = Vec::new();

    for (cluster_idx, cluster) in clusters.iter().enumerate() {
        for ((fwd_gap_start, fwd_gap_end), gap_seed_idx) in cluster.gaps(read_len, min_gap) {
            for (filler_idx, coverage) in
                index.find_overlapping(fwd_gap_start, fwd_gap_end, tolerance, cluster_idx)
            {
                if coverage >= min_coverage {
                    fills.push(GapFill {
                        cluster_idx,
                        gap: (fwd_gap_start, fwd_gap_end),
                        gap_seed_idx,
                        filler_idx,
                        coverage,
                    });
                }
            }
        }
    }

    fills
}
