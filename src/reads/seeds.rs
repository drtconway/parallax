use ordered_float::OrderedFloat;

use crate::align::Alignment;
use crate::utils::join::{Joinable, sorted_join};
use crate::{error::Result, writer::AlignmentWriter};

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
    /// Uniqueness of the most unique k-mer incorporated in this seed (lower is more unique)
    pub kmer_uniqueness: u32,
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

    pub fn fwd_read_range(&self, read_len: usize, is_reverse: bool) -> (usize, usize) {
        if is_reverse {
            (read_len - self.read_end(), read_len - self.read_pos)
        } else {
            (self.read_pos, self.read_end())
        }
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
        kmer_uniqueness: u32,
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
            None
        } else {
            // Does not overlap - return new seed hit
            Some(SeedHit::new(
                chrom_id,
                chrom_pos,
                read_pos,
                kmer,
                kmer_uniqueness,
                k,
            ))
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

        let u = self.kmer_uniqueness;
        let mapq = 60 / u as u8;

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
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            read_id,
            flag,
            chrom_name,
            self.ref_pos + 1, // 1-based
            mapq,
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
#[derive(Clone, Debug)]
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

    /// Alignments across gaps between seeds.
    /// Initially empty. After calling `align_gaps()`, contains one entry per gap
    /// (i.e., `chain.len() - 1` entries). Each entry is `Some(alignment)` if the
    /// WFA alignment succeeded, or `None` if it failed or was skipped.
    ///
    /// This allows gap-splitting decisions to consider actual alignment quality
    /// rather than just seed absence.
    pub gap_alignments: Vec<Option<Alignment>>,
}

impl SeedCluster {
    /// Create a new cluster from a chain of seeds.
    ///
    /// The chain is sorted by read position and overlapping seeds are resolved
    /// by truncating the right-hand seed. Seeds that become too short (less than
    /// `min_seed_length`) are dropped entirely.
    ///
    /// # Arguments
    /// * `chain` - The seeds to form into a cluster
    /// * `is_reverse` - Whether this cluster is on the reverse strand
    /// * `min_seed_length` - Minimum seed length after truncation (typically K/2)
    pub fn new(mut chain: Vec<SeedHit>, is_reverse: bool, min_seed_length: usize) -> Option<Self> {
        if chain.is_empty() {
            return None;
        }

        // Sort by strand read position for alignment building
        // (For forward strand this == forward order; for reverse strand this is RC order)
        chain.sort_by_key(|hit| hit.read_pos);

        // Resolve overlaps by truncating right-hand seeds
        Self::resolve_overlaps(&mut chain, min_seed_length);

        if chain.is_empty() {
            return None;
        }

        let read_start = chain.first().map(|h| h.read_pos).unwrap_or(0);
        let read_end = chain.last().map(|h| h.read_end()).unwrap_or(0);
        let chrom_id = chain[0].chrom_id;

        Some(SeedCluster {
            chrom_id,
            read_start,
            read_end,
            chain,
            is_reverse,
            gap_alignments: Vec::new(),
        })
    }

    pub fn ref_start(&self) -> usize {
        self.chain.iter().map(|h| h.ref_pos).min().unwrap_or(0)
    }

    pub fn ref_end(&self) -> usize {
        self.chain.iter().map(|h| h.ref_end()).max().unwrap_or(0)
    }

    pub fn fwd_read_range(&self, read_len: usize) -> (usize, usize) {
        if self.is_reverse {
            (read_len - self.read_end, read_len - self.read_start)
        } else {
            (self.read_start, self.read_end)
        }
    }

    /// Resolve overlapping seeds by truncating the right-hand seed.
    ///
    /// When consecutive seeds overlap (in read coordinates), the second seed
    /// is truncated so it starts where the first seed ends. Seeds that become
    /// shorter than `min_seed_length` are dropped entirely.
    ///
    /// This ensures downstream code never has to deal with overlapping seeds,
    /// simplifying gap calculation and CIGAR generation.
    fn resolve_overlaps(chain: &mut Vec<SeedHit>, min_seed_length: usize) {
        if chain.len() < 2 {
            return;
        }

        let mut write_idx = 1; // First seed always kept
        let mut prev_read_end = chain[0].read_end();
        let mut prev_ref_end = chain[0].ref_end();

        for read_idx in 1..chain.len() {
            // Copy the hit since SeedHit is Copy
            let hit = chain[read_idx];

            // Calculate overlap in both read and reference space
            let read_overlap = prev_read_end.saturating_sub(hit.read_pos);
            let ref_overlap = prev_ref_end.saturating_sub(hit.ref_pos);

            // Use the maximum overlap to ensure consistency
            let overlap = read_overlap.max(ref_overlap);

            if overlap == 0 {
                // No overlap - keep seed as-is
                chain[write_idx] = hit;
                prev_read_end = hit.read_end();
                prev_ref_end = hit.ref_end();
                write_idx += 1;
            } else if overlap < hit.match_len {
                // Partial overlap - truncate the seed
                let new_match_len = hit.match_len - overlap;
                if new_match_len >= min_seed_length {
                    // Keep the truncated seed
                    let mut truncated = hit;
                    truncated.read_pos += overlap;
                    truncated.ref_pos += overlap;
                    truncated.match_len = new_match_len;
                    // Diagonal stays the same (ref_pos - read_pos is unchanged)
                    chain[write_idx] = truncated;
                    prev_read_end = truncated.read_end();
                    prev_ref_end = truncated.ref_end();
                    write_idx += 1;
                }
                // else: truncated seed too short, drop it
            }
            // else: fully overlapped, drop the seed
        }

        chain.truncate(write_idx);
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
        self.read_end = self
            .chain
            .last()
            .map(|h| h.read_end())
            .unwrap_or(self.read_start);

        // Split gap_alignments if populated
        // gap_seed_idx corresponds to the gap between seeds gap_seed_idx and gap_seed_idx+1
        // After split: self keeps gaps 0..gap_seed_idx-1, tail gets gaps gap_seed_idx+1..
        // The gap at gap_seed_idx is discarded (it's the split point)
        let tail_gap_alignments = if self.gap_alignments.len() > gap_seed_idx + 1 {
            self.gap_alignments.split_off(gap_seed_idx + 1)
        } else {
            Vec::new()
        };
        // Remove the gap alignment at gap_seed_idx (the split point)
        // After split_off, self has alignments 0..=gap_seed_idx, but needs 0..gap_seed_idx-1
        self.gap_alignments.pop();

        // Build the new cluster from the tail
        let tail_read_start = tail_chain.first().map(|h| h.read_pos).unwrap_or(0);
        let tail_read_end = tail_chain.last().map(|h| h.read_end()).unwrap_or(0);

        Some(SeedCluster {
            read_start: tail_read_start,
            read_end: tail_read_end,
            chain: tail_chain,
            is_reverse: self.is_reverse,
            chrom_id: self.chrom_id,
            gap_alignments: tail_gap_alignments,
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

    /// Compute chain score using minimap2-style gap penalties.
    ///
    /// Score = Σ(anchor_length) - Σ(gap_penalty)
    /// Gap penalty = α * gap_len + β * log2(gap_len)
    ///
    /// This rewards long seed matches while penalizing gaps between seeds.
    /// A single long seed scores well (no gaps), while many small scattered
    /// seeds with large gaps score poorly.
    ///
    /// # Arguments
    /// * `gap_penalty_linear` - α coefficient (typically 0.01)
    /// * `gap_penalty_log` - β coefficient (typically 0.5)
    pub fn chain_score(&self, gap_penalty_linear: f64, gap_penalty_log: f64) -> f64 {
        if self.chain.is_empty() {
            return 0.0;
        }

        // Sum of all anchor lengths
        let anchor_score: f64 = self.chain.iter().map(|h| h.match_len as f64).sum();

        // Sum of gap penalties between consecutive seeds
        let gap_penalty: f64 = self
            .chain
            .windows(2)
            .map(|pair| {
                // Gap in read space
                let read_gap = if pair[1].read_pos > pair[0].read_end() {
                    pair[1].read_pos - pair[0].read_end()
                } else {
                    0
                };
                // Gap in reference space
                let ref_gap = if pair[1].ref_pos > pair[0].ref_end() {
                    pair[1].ref_pos - pair[0].ref_end()
                } else {
                    0
                };
                // Use the larger gap (handles indels)
                let gap = read_gap.max(ref_gap);

                if gap <= 1 {
                    0.0
                } else {
                    gap_penalty_linear * gap as f64 + gap_penalty_log * (gap as f64).log2()
                }
            })
            .sum();

        anchor_score - gap_penalty
    }

    pub fn score(&self) -> f64 {
        let mut s = 0.0;
        for seed in &self.chain {
            let l = seed.match_len as f64;
            s += l * (1.0 + l.log2());
        }
        for aln in &self.gap_alignments {
            if let Some(a) = aln {
                s += a.score();
            }
        }
        s
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

    /// Align across all gaps between seeds using WFA.
    ///
    /// Populates `gap_alignments` with one entry per gap (chain.len() - 1 entries).
    /// Each entry is `Some(alignment)` if alignment succeeded, `None` otherwise.
    ///
    /// # Arguments
    /// * `read_seq` - The read sequence (strand-specific, already rev-comped if reverse)
    /// * `ref_seq` - The reference sequence for this chromosome
    ///
    /// This should be called before gap analysis to enable alignment-aware
    /// split decisions.
    pub fn align_gaps(&mut self, read_seq: &[u8], ref_seq: &[u8]) {
        use crate::align::align;

        let num_gaps = self.chain.len().saturating_sub(1);
        self.gap_alignments = Vec::with_capacity(num_gaps);

        for pair in self.chain.windows(2) {
            let seed1 = &pair[0];
            let seed2 = &pair[1];

            // Extract gap regions
            let read_gap_start = seed1.read_end();
            let read_gap_end = seed2.read_pos;
            let ref_gap_start = seed1.ref_end();
            let ref_gap_end = seed2.ref_pos;

            assert!(read_gap_start <= read_gap_end);
            assert!(ref_gap_start <= ref_gap_end);

            // Check for valid gap (one of the sequences must have a gap)
            if read_gap_end == read_gap_start && ref_gap_end == ref_gap_start {
                // No gap or negative gap - shouldn't happen after overlap resolution
                self.gap_alignments.push(None);
                continue;
            }

            let read_gap = &read_seq[read_gap_start..read_gap_end];
            let ref_gap = &ref_seq[ref_gap_start..ref_gap_end];

            // Align the gap regions
            let alignment = align(read_gap, ref_gap);
            self.gap_alignments.push(alignment);
        }
    }

    /// Get the alignment for a specific gap by index.
    ///
    /// `gap_idx` is the index of the seed before the gap (0-based).
    /// Returns `None` if gaps haven't been aligned yet or if the gap index is invalid.
    pub fn gap_alignment(&self, gap_idx: usize) -> Option<&Alignment> {
        self.gap_alignments.get(gap_idx).and_then(|a| a.as_ref())
    }

    /// Check if gap alignments have been computed.
    pub fn has_gap_alignments(&self) -> bool {
        !self.gap_alignments.is_empty()
    }

    /// Split this cluster at all gaps where alignment failed (returned None).
    ///
    /// Must be called after `align_gaps()`. Modifies `self` in place to contain
    /// only the seeds before the first failed gap, and returns any new clusters
    /// created from splits.
    ///
    /// Returns an empty Vec if no splits are needed (common case).
    ///
    /// The returned clusters will preserve their gap alignments for gaps that succeeded.
    ///
    /// # Arguments
    /// * `_min_seed_length` - Previously used for validation, now unused but kept for API compatibility
    pub fn split_at_failed_alignments(&mut self, _min_seed_length: usize) -> Vec<SeedCluster> {
        if !self.has_gap_alignments() || self.chain.len() < 2 {
            // No gap alignments computed or only one seed - nothing to split
            return Vec::new();
        }

        // Find indices of gaps where alignment failed (None)
        let failed_gaps: Vec<usize> = self
            .gap_alignments
            .iter()
            .enumerate()
            .filter_map(|(i, alignment)| if alignment.is_none() { Some(i) } else { None })
            .collect();

        if failed_gaps.is_empty() {
            // All alignments succeeded - nothing to split
            return Vec::new();
        }

        // Split the chain at each failed gap, working backwards to keep indices valid
        let mut result = Vec::with_capacity(failed_gaps.len());

        // Process from last failed gap to first
        for &gap_idx in failed_gaps.iter().rev() {
            if let Some(tail) = self.split_at_gap(gap_idx) {
                // split_at_gap already correctly partitions gap_alignments:
                // - self keeps gaps 0..gap_idx (excluding the failed gap)
                // - tail gets gaps gap_idx+1.. (excluding the failed gap)
                if !tail.chain.is_empty() {
                    result.push(tail);
                }
            }
        }

        // Reverse to get clusters in read-position order
        result.reverse();
        result
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

#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct ClusterGap {
    gap_start: usize,
    gap_end: usize,
    cluster_idx: usize,
    gap_seed_idx: usize,
    gap_score: f64,
}

impl ClusterGap {
    fn new(gap_start: usize, gap_end: usize, cluster_idx: usize, gap_seed_idx: usize, gap_score: f64) -> Self {
        Self {
            gap_start,
            gap_end,
            cluster_idx,
            gap_seed_idx,
            gap_score,
        }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.gap_end - self.gap_start
    }
}

impl Joinable for ClusterGap {
    fn range(&self) -> (usize, usize) {
        (self.gap_start, self.gap_end)
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct ClusterAlignment {
    cluster_start: usize,
    cluster_end: usize,
    cluster_idx: usize,
    cluster_score: f64,
}

impl ClusterAlignment {
    fn new(cluster_start: usize, cluster_end: usize, cluster_idx: usize, cluster_score: f64) -> Self {
        Self {
            cluster_start,
            cluster_end,
            cluster_idx,
            cluster_score,
        }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.cluster_end - self.cluster_start
    }
}

impl Joinable for ClusterAlignment {
    fn range(&self) -> (usize, usize) {
        (self.cluster_start, self.cluster_end)
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
    min_fill: usize,
    tolerance: usize,
    min_coverage: f64,
) -> Vec<GapFill> {
    if clusters.len() < 2 {
        return Vec::new();
    }

    let start = std::time::Instant::now();

    let mut gaps: Vec<ClusterGap> = Vec::new();
    let mut alignments: Vec<ClusterAlignment> = Vec::new();
    for (cluster_idx, cluster) in clusters.iter().enumerate() {
        // Add alignment region for this cluster
        let (fwd_start, fwd_end) = cluster.fwd_read_range(read_len);
        if fwd_end > fwd_start && fwd_end - fwd_start >= min_fill {
            log::info!(
                "Cluster {} ({}:{}-{}) alignment region: ({},{}) - {}bp",
                cluster_idx,
                cluster.chrom_id,
                cluster.ref_start(),
                cluster.ref_end(),
                fwd_start,
                fwd_end,
                fwd_end - fwd_start
            );
            alignments.push(ClusterAlignment::new(fwd_start, fwd_end, cluster_idx, cluster.score()));
        }

        // Find gaps in this cluster
        // Seeds are stored in strand order, so for reverse strand the forward coords
        // are in reverse order. We need to compute gaps in forward coordinate space.
        if cluster.chain.len() < 2 {
            continue;
        }
        let n = cluster.chain.len();
        for i in 0..n - 1 {
            let seed1 = &cluster.chain[i];
            let seed1_range = seed1.fwd_read_range(read_len, cluster.is_reverse);
            let seed2 = &cluster.chain[i + 1];
            let seed2_range = seed2.fwd_read_range(read_len, cluster.is_reverse);
            
            // Compute gap in forward coordinates
            // For forward strand: gap is from seed1.end to seed2.start
            // For reverse strand: seeds are in reverse fwd order, so gap is from seed2.end to seed1.start
            let (gap_start, gap_end) = if cluster.is_reverse {
                (seed2_range.1, seed1_range.0)
            } else {
                (seed1_range.1, seed2_range.0)
            };
            
            if cluster_idx == 0 {
                log::info!(
                    "Cluster {} ({}:{}-{}) gap between seeds {} and {}: ({},{}) - {}bp (reverse={})",
                    cluster_idx,
                    cluster.chrom_id,
                    cluster.ref_start(),
                    cluster.ref_end(),
                    i,
                    i + 1,
                    gap_start,
                    gap_end,
                    gap_end as i64 - gap_start as i64,
                    cluster.is_reverse
                );
            }
            
            if gap_end > gap_start && gap_end - gap_start >= min_gap {
                log::info!(
                    "Cluster {} ({}:{}-{}) large gap between seeds {} and {}: ({},{}) - {}bp",
                    cluster_idx,
                    cluster.chrom_id,
                    cluster.ref_start(),
                    cluster.ref_end(),
                    i,
                    i + 1,
                    gap_start,
                    gap_end,
                    gap_end - gap_start
                );
                let aln_score = cluster.gap_alignment(i).map(|a| a.score()).unwrap_or(0.0);
                gaps.push(ClusterGap::new(gap_start, gap_end, cluster_idx, i, aln_score));
            }
        }
    }
    gaps.sort_by_key(|gap| {
        (gap.gap_start, gap.gap_end, gap.cluster_idx)
    });
    alignments.sort_by_key(|aln| {  
        (aln.cluster_start, aln.cluster_end, aln.cluster_idx)
    });

    // Use sorted_join to find gap-alignment pairs from different clusters
    let mut pairs = sorted_join(&gaps, &alignments, tolerance, |gap, aln| {
        // Must be from different clusters
        gap.cluster_idx != aln.cluster_idx
    });

    pairs.sort_by_key(|(gap_idx, aln_idx)| {
        let gap = &gaps[*gap_idx];
        let aln = &alignments[*aln_idx];
        OrderedFloat(aln.cluster_score - gap.gap_score)
    });

    // Convert pairs to GapFill entries with coverage calculation
    let fills: Vec<GapFill> = pairs
        .into_iter()
        .filter_map(|(gap_idx, aln_idx)| {
            let gap = &gaps[gap_idx];
            let aln = &alignments[aln_idx];
            let gap_len = gap.gap_end - gap.gap_start;
            let aln_len = aln.cluster_end - aln.cluster_start;
            let ratio = aln_len as f64 / gap_len as f64;
            let diff = gap_len as isize - aln_len as isize;
            log::info!(
                "Gap fill: cluster {} gap ({},{} - {}bp, score {:.2}) filled by cluster {} aln ({},{} - {}bp, score {:.2}) | ratio {:.2}, diff {}",
                gap.cluster_idx,
                gap.gap_start,
                gap.gap_end,
                gap_len,
                gap.gap_score,
                aln.cluster_idx,
                aln.cluster_start,
                aln.cluster_end,
                aln_len,
                aln.cluster_score,
                ratio,
                diff
            );

            // Calculate coverage: fraction of gap covered by the alignment
            let gap_len = gap.gap_end - gap.gap_start;
            let overlap_start = gap
                .gap_start
                .max(aln.cluster_start.saturating_sub(tolerance));
            let overlap_end = gap.gap_end.min(aln.cluster_end.saturating_add(tolerance));
            let overlap = overlap_end.saturating_sub(overlap_start);
            let coverage = overlap as f64 / gap_len as f64;

            if coverage >= min_coverage {
                Some(GapFill {
                    cluster_idx: gap.cluster_idx,
                    gap: (gap.gap_start, gap.gap_end),
                    gap_seed_idx: gap.gap_seed_idx,
                    filler_idx: aln.cluster_idx,
                    coverage,
                })
            } else {
                None
            }
        })
        .collect();

    metrics::histogram!("analysis.gap_fills").record(start.elapsed().as_micros() as f64);

    fills
}
