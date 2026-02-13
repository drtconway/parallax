use core::panic;
use std::collections::HashMap;

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::align::block::BlockAligner;
use crate::align::{AlignParams, Alignment, CigarOp};
use crate::scores::QualityScore;
use crate::utils::debug::{self, DebugFile};
use crate::utils::join::{Joinable, sorted_join};
use crate::utils::{GroupsTrait, InterleaveTrait};
use crate::{error::Result, writer::AlignmentWriter};

#[derive(Debug)]
pub enum ClusterError {
    SequenceMismatch {
        read_bases: String,
        ref_bases: String,
        read_pos: usize,
        ref_pos: usize,
    },
    AlignmentMismatch {
        gap_index: usize,
        read_start: usize,
        read_end: usize,
        ref_start: usize,
        ref_end: usize,
        error: String,
    },
}

impl std::fmt::Display for ClusterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClusterError::SequenceMismatch {
                read_bases,
                ref_bases,
                read_pos,
                ref_pos,
            } => write!(
                f,
                "Seed at read_pos={} ref_pos={} does not match reference: read_bases='{}' ref_bases='{}'",
                read_pos, ref_pos, read_bases, ref_bases
            ),
            ClusterError::AlignmentMismatch {
                gap_index,
                read_start,
                read_end,
                ref_start,
                ref_end,
                error,
            } => write!(
                f,
                "Alignment mismatch at gap {}: read[{}-{}] ref[{}-{}]: {}",
                gap_index, read_start, read_end, ref_start, ref_end, error
            ),
        }
    }
}

impl std::error::Error for ClusterError {}

/// A seed hit representing a k-mer match between read and reference
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
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
}

impl SeedHit {
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
            "",
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
    /// * `strand_qual` - Quality scores for this strand (already reversed if reverse)
    pub fn to_sam_line(
        &self,
        read_id: &str,
        chrom_name: &str,
        is_reverse: bool,
        strand_seq: &[u8],
        strand_qual: &[u8],
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
        let qual_str = {
            let qual_slice = &strand_qual[self.read_pos..self.read_pos + self.match_len];
            // Convert Phred+33 quality scores to ASCII
            qual_slice.iter().map(|&q| q as char).collect::<String>()
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

    /// The chain of colinear seeds (sorted by read position)
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
    pub gap_alignments: Vec<Alignment>,
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

    pub fn into_alignment(
        self,
        fwd_extend_left: bool,
        fwd_extend_right: bool,
        soft_clip: bool,
        strand_seq: &[u8],
        ref_seq: &[u8],
        aligner: &mut BlockAligner,
    ) -> Alignment {
        let mut left_extension: Option<Vec<CigarOp>> = None;
        let mut right_extension: Option<Vec<CigarOp>> = None;

        // Extend left with respect to the forward strand, if requested.
        // This means:
        // - For forward clusters, extend left from the first seed's start.
        // - For reverse clusters, extend left from the last seed's end (which is the leftmost point on the forward strand).
        // In both cases, the extension is performed in the forward direction of the strand sequence, which corresponds to leftward extension on the forward strand.
        // The alignment is performed using the x-drop extension algorithm.
        if fwd_extend_left {
            if self.is_reverse {
                // Extend right with respect to the strand sequence.
                if self.read_end < strand_seq.len() {
                    let read_ext_seq = &strand_seq[self.read_end..];
                    let ref_ext_seq = &ref_seq[self.chain.last().unwrap().ref_end()..];
                    let ext = aligner
                        .extend_right(read_ext_seq, ref_ext_seq)
                        .expect("alignment extension failed");
                    right_extension = Some(ext.cigar);
                }
            } else {
                // Extend left with respect to the strand sequence.
                if self.read_start > 0 {
                    let read_ext_seq = &strand_seq[..self.read_start];
                    let ref_ext_seq = &ref_seq[..self.ref_start()];
                    let ext = aligner
                        .extend_left(read_ext_seq, ref_ext_seq)
                        .expect("alignment extension failed");
                    left_extension = Some(ext.cigar);
                }
            }
        }

        // Extend right with respect to the forward strand, if requested.
        // This means:
        // - For forward clusters, extend right from the last seed's end.
        // - For reverse clusters, extend right from the first seed's start (which is the rightmost point on the forward strand).
        // In both cases, the extension is performed in the forward direction of the strand sequence, which corresponds to rightward extension on the forward strand.
        // The alignment is performed using the x-drop extension algorithm.
        if fwd_extend_right {
            if self.is_reverse {
                // Extend left with respect to the strand sequence.
                if self.read_start > 0 {
                    let read_ext_seq = &strand_seq[..self.read_start];
                    let ref_ext_seq = &ref_seq[..self.ref_start()];
                    let ext = aligner
                        .extend_left(read_ext_seq, ref_ext_seq)
                        .expect("alignment extension failed");
                    left_extension = Some(ext.cigar);
                }
            } else {
                // Extend right with respect to the strand sequence.
                if self.read_end < strand_seq.len() {
                    let read_ext_seq = &strand_seq[self.read_end..];
                    let ref_ext_seq = &ref_seq[self.chain.last().unwrap().ref_end()..];
                    let ext = aligner
                        .extend_right(read_ext_seq, ref_ext_seq)
                        .expect("alignment extension failed");
                    right_extension = Some(ext.cigar);
                }
            }
        }

        if left_extension.is_none() && self.read_start > 0 {
            left_extension = if soft_clip {
                Some(vec![CigarOp::SoftClip(self.read_start as u32)])
            } else {
                Some(vec![CigarOp::HardClip(self.read_start as u32)])
            };
        }

        if right_extension.is_none() && self.read_end < strand_seq.len() {
            let clip_len = strand_seq.len() - self.read_end;
            right_extension = if soft_clip {
                Some(vec![CigarOp::SoftClip(clip_len as u32)])
            } else {
                Some(vec![CigarOp::HardClip(clip_len as u32)])
            };
        }

        let seed_parts = self
            .chain
            .into_iter()
            .map(|hit| vec![CigarOp::Match(hit.match_len as u32)]);
        let gap_alignments = self.gap_alignments.into_iter().map(|aln| aln.cigar);
        let interleaved = seed_parts.interleave(gap_alignments);
        let cigar_ops: Vec<CigarOp> = left_extension
            .into_iter()
            .chain(interleaved)
            .chain(right_extension.into_iter())
            .flatten()
            .collect();
        Alignment::from(cigar_ops)
    }

    /// Return a summary for populating the supplementary alignment (SA) tag.
    /// (chrom_id, read_start, read_end, is_reverse, condensed CIGAR string, number of mismatches)
    pub fn summary(&self, read_len: usize) -> (usize, usize, bool, String, usize) {
        let cigar = self.cigar_summary(read_len);
        let mismatch_count = self.mismatch_count();
        (
            self.chrom_id,
            self.ref_start(),
            self.is_reverse,
            cigar,
            mismatch_count,
        )
    }

    /// Reference start position of the seed chain.
    pub fn ref_start(&self) -> usize {
        self.chain.first().map(|h| h.ref_pos).unwrap_or(0)
    }

    /// Reference end position of the seed chain.
    pub fn ref_end(&self) -> usize {
        self.chain.last().map(|h| h.ref_end()).unwrap_or(0)
    }

    pub fn fwd_read_range(&self, read_len: usize) -> (usize, usize) {
        if self.is_reverse {
            (read_len - self.read_end, read_len - self.read_start)
        } else {
            (self.read_start, self.read_end)
        }
    }

    pub fn read_range(&self) -> std::ops::Range<usize> {
        self.read_start..self.read_end
    }

    /// Generate a condensed CIGAR string summarizing the accumulated insertions, deletions, matches and mismatches.
    /// This is used for the SA tag and other summaries where we want a quick representation of the alignment without the full CIGAR.
    pub fn cigar_summary(&self, read_len: usize) -> String {
        let mut matches = 0;
        let mut mismatches = 0;
        let mut indels: i64 = 0;
        for seed in &self.chain {
            matches += seed.match_len;
        }
        for gap in &self.gap_alignments {
            for op in &gap.cigar {
                match op {
                    CigarOp::Match(n) => matches += *n as usize,
                    CigarOp::Mismatch(n) => mismatches += *n as usize,
                    CigarOp::Ins(n) => indels += *n as i64,
                    CigarOp::Del(n) => indels -= *n as i64,
                    _ => {}
                }
            }
        }
        let left_clip = if self.read_start > 0 {
            format!("{}S", self.read_start)
        } else {
            String::new()
        };
        let right_clip = if self.read_end < read_len {
            format!("{}S", read_len - self.read_end)
        } else {
            String::new()
        };
        if indels == 0 {
            format!("{}{}={}X{}", left_clip, matches, mismatches, right_clip)
        } else if indels > 0 {
            format!("{}{}={}X{}I{}", left_clip, matches, mismatches, indels, right_clip)
        } else {
            format!("{}{}={}X{}D{}", left_clip, matches, mismatches, -indels, right_clip)
        }
    }

    /// Calculate the total number of mismatches across all gap alignments in this cluster.
    /// Seed matches are assumed to be perfect and are not counted as mismatches.
    pub fn mismatch_count(&self) -> usize {
        self.gap_alignments
            .iter()
            .map(|align| align.mismatch_count())
            .sum()
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

    pub fn diagonal(&self) -> f64 {
        let mut sum = 0i64;
        let mut total_length = 0usize;
        for seed in &self.chain {
            sum += seed.diagonal * seed.match_len as i64;
            total_length += seed.match_len;
        }
        if total_length == 0 {
            0.0
        } else {
            sum as f64 / total_length as f64
        }
    }

    pub fn total_identity(&self) -> usize {
        let total_seed_length: usize = self.chain.iter().map(|h| h.match_len).sum();
        let total_aligned_identities: usize = self
            .gap_alignments
            .iter()
            .map(|align| align.total_identity())
            .sum();
        total_seed_length + total_aligned_identities
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
    pub fn split_at_gap(&mut self, gap_seed_idx: usize) -> Option<(SeedCluster, Alignment)> {
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
        let dropped_alignment = self.gap_alignments.pop().unwrap();

        // Build the new cluster from the tail
        let tail_read_start = tail_chain.first().map(|h| h.read_pos).unwrap_or(0);
        let tail_read_end = tail_chain.last().map(|h| h.read_end()).unwrap_or(0);

        Some((
            SeedCluster {
                read_start: tail_read_start,
                read_end: tail_read_end,
                chain: tail_chain,
                is_reverse: self.is_reverse,
                chrom_id: self.chrom_id,
                gap_alignments: tail_gap_alignments,
            },
            dropped_alignment,
        ))
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

    /// Write this chain to a debug SAM file with SA tags linking all seeds.
    ///
    /// Each seed in the chain is written as a SAM record. The primary alignment
    /// is the first seed, and all others are marked as supplementary (0x800).
    /// All records include an SA tag listing all other alignments in the chain.
    ///
    /// # Arguments
    /// * `read_name` - Read name for SAM output
    /// * `cluster_id` - Cluster identifier for the cc tag
    /// * `chrom_name` - Chromosome name
    /// * `strand_seq` - Sequence for this strand
    /// * `strand_qual` - Quality scores (optional)
    pub fn write_chain_sam(
        &self,
        read_name: &str,
        cluster_id: usize,
        chrom_name: &str,
        strand_seq: &[u8],
        strand_qual: &[u8],
    ) {
        use crate::utils::debug::{self, DebugFile};

        if self.chain.is_empty() {
            return;
        }

        let read_len = strand_seq.len();
        let strand_char = if self.is_reverse { '-' } else { '+' };

        // Build SA tag entries for all seeds in the chain
        // Format: rname,pos,strand,CIGAR,mapQ,NM;
        let sa_entries: Vec<String> = self
            .chain
            .iter()
            .map(|seed| {
                let cigar = format!("{}M", seed.match_len);
                let mapq = 60 / seed.kmer_uniqueness.max(1) as u8;
                format!(
                    "{},{},{},{},{},0",
                    chrom_name,
                    seed.ref_pos + 1,
                    strand_char,
                    cigar,
                    mapq
                )
            })
            .collect();

        // Write each seed as a SAM record
        for (i, seed) in self.chain.iter().enumerate() {
            // Build SA tag excluding the current seed
            let sa_others: Vec<&str> = sa_entries
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, s)| s.as_str())
                .collect();
            let sa_tag = if sa_others.is_empty() {
                String::new()
            } else {
                format!("\tSA:Z:{};", sa_others.join(";"))
            };

            // Flag: primary (0) for first seed, supplementary (0x800) for rest
            // Plus reverse flag (0x10) if reverse strand
            let mut flag = if i == 0 {
                0u16
            } else {
                super::FLAG_SUPPLEMENTARY
            };
            if self.is_reverse {
                flag |= super::FLAG_REVERSE;
            }

            // Build CIGAR with hard clips
            let hclip_start = seed.read_pos;
            let hclip_end = read_len.saturating_sub(seed.read_pos + seed.match_len);
            let cigar = match (hclip_start > 0, hclip_end > 0) {
                (true, true) => format!("{}H{}={}H", hclip_start, seed.match_len, hclip_end),
                (true, false) => format!("{}H{}=", hclip_start, seed.match_len),
                (false, true) => format!("{}={}H", seed.match_len, hclip_end),
                (false, false) => format!("{}=", seed.match_len),
            };

            let mapq = 60 / seed.kmer_uniqueness.max(1) as u8;

            // Extract aligned sequence and quality
            let seq_slice = &strand_seq[seed.read_pos..seed.read_pos + seed.match_len];
            let seq_str = String::from_utf8_lossy(seq_slice);

            let qual_str = {
                let qual_slice = &strand_qual[seed.read_pos..seed.read_pos + seed.match_len];
                qual_slice.iter().map(|&q| q as char).collect::<String>()
            };

            // Write SAM line with SA tag and cluster ID tag
            debug::write(
                DebugFile::ChainsSam,
                &format!(
                    "{}.{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}{}\tcc:i:{}",
                    read_name,
                    cluster_id,
                    flag,
                    chrom_name,
                    seed.ref_pos + 1,
                    mapq,
                    cigar,
                    seq_str,
                    qual_str,
                    sa_tag,
                    cluster_id
                ),
            );
        }
    }

    pub fn quality(&self, params: &AlignParams) -> QualityScore {
        let mut s = 0.0;
        for (i, seed) in self.chain.iter().enumerate() {
            let op = CigarOp::Match(seed.match_len as u32);
            let score = op.quality(params).0;
            log::debug!("Seed {}: length = {}, score = {}", i, seed.match_len, score);
            s += score;
        }
        for (i, aln) in self.gap_alignments.iter().enumerate() {
            let a_len = aln.query_length();
            let a_score = aln.quality(params).0;
            log::debug!(
                "Gap {}: alignment length = {}, score = {}",
                i,
                a_len,
                a_score
            );
            s += a_score;
        }
        QualityScore::from(s)
    }

    #[allow(dead_code)]
    pub fn cigar_string(&self, read_len: usize) -> String {
        let mut cigar_ops = Vec::new();

        if let Some(first_seed) = self.chain.first() {
            // Handle leading hard clip
            let (fwd_start, _) = first_seed.fwd_read_range(read_len, self.is_reverse);
            if fwd_start > 0 {
                cigar_ops.push(format!("{}H", fwd_start));
            }
        }

        for (i, seed) in self.chain.iter().enumerate() {
            // Add the alignment before this seed if not the first seed
            if i > 0 {
                let gap_aln = &self.gap_alignments[i - 1];
                // Use the gap alignment CIGAR
                cigar_ops.push(gap_aln.cigar_string());
            }

            // Add match operation for this seed
            cigar_ops.push(format!("{}=", seed.match_len));
        }

        if let Some(last_seed) = self.chain.last() {
            // Handle trailing hard clip
            let (_, fwd_end) = last_seed.fwd_read_range(read_len, self.is_reverse);
            if fwd_end < read_len {
                let hclip_len = read_len - fwd_end;
                cigar_ops.push(format!("{}H", hclip_len));
            }
        }

        cigar_ops.join("")
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
    /// * `min_seed_length` - Minimum seed length for creating new clusters if splitting is needed
    ///
    /// This should be called before gap analysis to enable alignment-aware
    /// split decisions.
    pub fn align_gaps(
        mut self,
        read_name: &str,
        read_seq: &[u8],
        ref_seq: &[u8],
        min_seed_length: usize,
    ) -> Vec<Self> {
        use crate::align::align;

        if self.chain.len() < 2 {
            // No gaps to align
            return vec![self];
        }

        let num_gaps = self.chain.len() - 1;

        let mut gap_alignments = Vec::with_capacity(num_gaps);
        for (i, pair) in self.chain.windows(2).enumerate() {
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
                panic!(
                    "Invalid gap between seeds {} and {}: no gap in either read or reference",
                    i,
                    i + 1
                );
            }

            let read_gap = &read_seq[read_gap_start..read_gap_end];
            let ref_gap = &ref_seq[ref_gap_start..ref_gap_end];

            // chr11:116,817,620-116,819,109
            //if ref_gap_start > 116_817_620 && ref_gap_end < 116_819_109  && read_gap.len() > 100 && ref_gap.len() > 100 {
            //    log::info!(
            //        "Aligning:\n{}\n{}",
            //        String::from_utf8_lossy(read_gap),
            //        String::from_utf8_lossy(ref_gap)
            //    );
            //}

            // Align the gap regions
            let alignment = align(read_gap, ref_gap);
            if alignment.is_none() {
                log::debug!(
                    "Gap alignment failed for read {}, after seed {}, read_gap_len={}, ref_gap_len={}",
                    read_name,
                    i,
                    read_gap.len(),
                    ref_gap.len()
                );
                if true || (read_gap.len() < 1000 && ref_gap.len() < 1000) {
                    // Write in FASTA format with descriptive headers
                    debug::write(
                        DebugFile::GapAlignments,
                        &format!(
                            ">{}:read:{}-{}\n{}\n>{}:ref:{}-{}\n{}",
                            read_name,
                            read_gap_start,
                            read_gap_end,
                            String::from_utf8_lossy(read_gap),
                            read_name,
                            ref_gap_start,
                            ref_gap_end,
                            String::from_utf8_lossy(ref_gap)
                        ),
                    );
                }
            }
            gap_alignments.push(alignment);
        }

        if !gap_alignments.iter().any(|a| a.is_none()) {
            self.gap_alignments = gap_alignments.into_iter().map(|a| a.unwrap()).collect();
            return vec![self];
        }

        let groups: Vec<Vec<Alignment>> = gap_alignments.into_iter().groups().collect();

        let mut clusters = Vec::new();
        let mut seed_idx = 0;

        for group in groups {
            // Each group of n alignments connects n+1 seeds
            let num_seeds = group.len() + 1;
            let segment_seeds: Vec<SeedHit> = self.chain[seed_idx..seed_idx + num_seeds].to_vec();

            if let Some(mut cluster) =
                SeedCluster::new(segment_seeds, self.is_reverse, min_seed_length)
            {
                cluster.gap_alignments = group;
                clusters.push(cluster);
            }

            seed_idx += num_seeds;
        }

        clusters
    }

    /// Get the alignment for a specific gap by index.
    ///
    /// `gap_idx` is the index of the seed before the gap (0-based).
    /// Returns `None` if gaps haven't been aligned yet or if the gap index is invalid.
    pub fn gap_alignment(&self, gap_idx: usize) -> Option<&Alignment> {
        self.gap_alignments.get(gap_idx)
    }

    /// Format a visual diagram of seeds and gaps in this cluster.
    ///
    /// Returns a pair of strings showing the query and reference views:
    /// ```text
    /// QRY: [- 65bp -] <-  8bp -> [- 111bp -] <- 19bp -> [- 44bp -]
    /// REF: [- 65bp -] <- 11bp -> [- 111bp -] <- 17bp -> [- 44bp -]
    /// ```
    ///
    /// Seeds are shown as `[- Nbp -]` and have the same width on both lines
    /// since they represent exact matches. Gaps are shown as `<- Nbp ->`
    /// with padding to keep seeds aligned between the two lines.
    pub fn format_seed_diagram(&self) -> (String, String) {
        if self.chain.is_empty() {
            return ("QRY:".to_string(), "REF:".to_string());
        }

        let mut qry_parts: Vec<String> = Vec::new();
        let mut ref_parts: Vec<String> = Vec::new();

        for (i, seed) in self.chain.iter().enumerate() {
            // Add the gap before this seed (if not the first seed)
            if i > 0 {
                let prev = &self.chain[i - 1];
                let qry_gap = seed.read_pos.saturating_sub(prev.read_end()) as i64;
                let ref_gap = seed.ref_pos.saturating_sub(prev.ref_end()) as i64;

                // Format gap strings and find the max width for alignment
                let qry_gap_str = format!("{}bp", qry_gap);
                let ref_gap_str = format!("{}bp", ref_gap);
                let max_width = qry_gap_str.len().max(ref_gap_str.len());

                // Pad to align
                qry_parts.push(format!(
                    " <- {:>width$} -> ",
                    qry_gap_str,
                    width = max_width
                ));
                ref_parts.push(format!(
                    " <- {:>width$} -> ",
                    ref_gap_str,
                    width = max_width
                ));
            }

            // Add the seed (same width on both lines)
            let seed_str = format!("[- {}bp -]", seed.match_len);
            qry_parts.push(seed_str.clone());
            ref_parts.push(seed_str);
        }

        let qry_line = format!("QRY: {}", qry_parts.join(""));
        let ref_line = format!("REF: {}", ref_parts.join(""));

        (qry_line, ref_line)
    }

    #[allow(dead_code)]
    pub fn validate(
        &self,
        read_seq: &[u8],
        ref_seq: &[u8],
    ) -> std::result::Result<(), ClusterError> {
        // Validate each seed
        for (_i, seed) in self.chain.iter().enumerate() {
            let seed_read = &read_seq[seed.read_pos..seed.read_end()];
            let seed_ref = &ref_seq[seed.ref_pos..seed.ref_end()];
            if seed_read != seed_ref {
                return Err(ClusterError::SequenceMismatch {
                    read_bases: String::from_utf8_lossy(seed_read).to_string(),
                    ref_bases: String::from_utf8_lossy(seed_ref).to_string(),
                    read_pos: seed.read_pos,
                    ref_pos: seed.ref_pos,
                });
            }
        }
        for (i, aln) in self.gap_alignments.iter().enumerate() {
            let seed1 = &self.chain[i];
            let seed2 = &self.chain[i + 1];
            let read_gap_start = seed1.read_end();
            let read_gap_end = seed2.read_pos;
            let ref_gap_start = seed1.ref_end();
            let ref_gap_end = seed2.ref_pos;

            let read_gap = &read_seq[read_gap_start..read_gap_end];
            let ref_gap = &ref_seq[ref_gap_start..ref_gap_end];

            if let Err(err) = aln.validate(ref_gap, read_gap, 0) {
                return Err(ClusterError::AlignmentMismatch {
                    gap_index: i,
                    read_start: read_gap_start,
                    read_end: read_gap_end,
                    ref_start: ref_gap_start,
                    ref_end: ref_gap_end,
                    error: err,
                });
            }
        }
        Ok(())
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
    gap_score: QualityScore,
}

impl ClusterGap {
    fn new(
        gap_start: usize,
        gap_end: usize,
        cluster_idx: usize,
        gap_seed_idx: usize,
        gap_score: QualityScore,
    ) -> Self {
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
    cluster_score: QualityScore,
    identity_pct: f64,
}

impl ClusterAlignment {
    fn new(
        cluster_start: usize,
        cluster_end: usize,
        cluster_idx: usize,
        cluster_score: QualityScore,
        identity_pct: f64,
    ) -> Self {
        Self {
            cluster_start,
            cluster_end,
            cluster_idx,
            cluster_score,
            identity_pct,
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
    read_name: &str,
    clusters: &[SeedCluster],
    read_len: usize,
    min_gap: usize,
    min_fill: usize,
    tolerance: usize,
    params: &AlignParams,
) -> Vec<GapFill> {
    if clusters.len() < 2 {
        return Vec::new();
    }
    let start = std::time::Instant::now();

    let dump = debug::is_enabled(DebugFile::GapFills);

    let mut gaps: Vec<ClusterGap> = Vec::new();
    let mut alignments: Vec<ClusterAlignment> = Vec::new();
    for (cluster_idx, cluster) in clusters.iter().enumerate() {
        // Add alignment region for this cluster
        let (fwd_start, fwd_end) = cluster.fwd_read_range(read_len);
        if fwd_end > fwd_start && fwd_end - fwd_start >= min_fill {
            let score = cluster.quality(params);
            let identity = cluster.total_identity();
            let total_read_length = fwd_end - fwd_start;
            let identity_pct = if total_read_length > 0 {
                100.0 * identity as f64 / total_read_length as f64
            } else {
                0.0
            };
            if identity_pct > 90.0 {
                log::debug!(
                    "Cluster {} alignment: read ({},{}) ref ({}-{}) length {} score {:.2} identity {:.2}%",
                    cluster_idx,
                    fwd_start,
                    fwd_end,
                    cluster.ref_start(),
                    cluster.ref_end(),
                    total_read_length,
                    score,
                    identity_pct
                );
                alignments.push(ClusterAlignment::new(
                    fwd_start,
                    fwd_end,
                    cluster_idx,
                    score,
                    identity_pct,
                ));
                if dump {
                    let chrom_name = format!("chr{}", cluster.chrom_id + 1); // it's a lie, but close enough for debugging
                    let strand = if cluster.is_reverse { '-' } else { '+' };
                    let _fill_len = total_read_length;
                    debug::write(
                        DebugFile::GapFills,
                        &format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            read_name,
                            read_len,
                            fwd_start,
                            fwd_end,
                            total_read_length,
                            -(cluster_idx as i64 + 1),
                            score,
                            chrom_name,
                            cluster.ref_start(),
                            cluster.ref_end(),
                            strand,
                        ),
                    );
                }
            }
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

            if gap_end > gap_start && gap_end - gap_start >= min_gap {
                let aln = cluster.gap_alignment(i).unwrap();
                let aln_score = cluster
                    .gap_alignment(i)
                    .map(|a| a.quality(params).0)
                    .unwrap_or(0.0);
                let ref_start = seed1.ref_end();
                let ref_end = seed2.ref_pos;
                let identity_pct = if aln.query_length() > 0 {
                    100.0 * aln.total_identity() as f64 / aln.query_length() as f64
                } else {
                    0.0
                };
                log::debug!(
                    "Cluster {} gap between seeds {} and {}: read ({},{}) ref ({}-{}) length {} score {:.2} identity {:.2}% | CIGAR: {}",
                    cluster_idx,
                    i,
                    i + 1,
                    gap_start,
                    gap_end,
                    ref_start,
                    ref_end,
                    gap_end - gap_start,
                    aln_score,
                    identity_pct,
                    aln.cigar_string()
                );
                gaps.push(ClusterGap::new(
                    gap_start,
                    gap_end,
                    cluster_idx,
                    i,
                    QualityScore::new(aln_score),
                ));
                if dump {
                    let chrom_name = format!("chr{}", cluster.chrom_id + 1); // it's a lie, but close enough for debugging
                    let strand = if cluster.is_reverse { '-' } else { '+' };
                    let fill_len = gap_end - gap_start;
                    let (ref_start, ref_end) = if cluster.is_reverse {
                        (seed2.ref_end(), seed1.ref_pos)
                    } else {
                        (seed1.ref_end(), seed2.ref_pos)
                    };
                    debug::write(
                        DebugFile::GapFills,
                        &format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            read_name,
                            read_len,
                            gap_start,
                            gap_end,
                            fill_len,
                            cluster_idx + 1,
                            aln_score,
                            chrom_name,
                            ref_start,
                            ref_end,
                            strand,
                        ),
                    );
                }
            }
        }
    }
    gaps.sort_by_key(|gap| (gap.gap_start, gap.gap_end, gap.cluster_idx));
    alignments.sort_by_key(|aln| (aln.cluster_start, aln.cluster_end, aln.cluster_idx));

    // Use sorted_join to find gap-alignment pairs from different clusters
    let mut pairs = sorted_join(&gaps, &alignments, tolerance, |gap, aln| {
        let gap_len = (gap.gap_end - gap.gap_start) as isize;
        let aln_len = (aln.cluster_end - aln.cluster_start) as isize;

        if gap.cluster_idx == aln.cluster_idx {
            return false;
        }
        if !gap.gap_score.is_worse_than(QualityScore::ZERO) {
            return false;
        }
        if !aln.cluster_score.is_better_than(QualityScore::ZERO) {
            return false;
        }

        let overlap_start = gap.gap_start.max(aln.cluster_start);
        let overlap_end = gap.gap_end.min(aln.cluster_end);
        if overlap_end <= overlap_start {
            return false;
        }
        let overlap = overlap_end - overlap_start;
        let gap_overlap_ratio = overlap as f64 / gap_len as f64;
        let fill_overlap_ratio = aln_len as f64 / overlap as f64;
        let ratio = gap_overlap_ratio * fill_overlap_ratio;
        if ratio >= 0.5 && ratio <= 1.1 {
            log::debug!(
                "Considering gap fill: gap cluster {} gap ({},{}) length {} score {:.2} | aln cluster {} aln ({},{}) length {} score {:.2} identity {:.2}% | gap ratio {:.2} fill ratio {:.2} | ratio {:.2}",
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
                aln.identity_pct,
                gap_overlap_ratio,
                fill_overlap_ratio,
                ratio
            );
        }
        // Must be from different clusters
        gap.cluster_idx != aln.cluster_idx
            && gap.gap_score.is_worse_than(QualityScore::ZERO)
            && aln.cluster_score.is_better_than(QualityScore::ZERO)
            && ratio >= 0.5
            && ratio <= 1.1
    });

    pairs.sort_by_key(|(gap_idx, aln_idx)| {
        let gap = &gaps[*gap_idx];
        let aln = &alignments[*aln_idx];
        OrderedFloat(gap.gap_score.0 - aln.cluster_score.0)
    });

    let mut selected = HashMap::new();
    for (gap_idx, aln_idx) in &pairs {
        if selected.contains_key(gap_idx) {
            continue;
        }
        selected.insert(*gap_idx, *aln_idx);
        let gap = &gaps[*gap_idx];
        let aln = &alignments[*aln_idx];
        log::debug!(
            "Selected gap fill pair: gap cluster {} gap ({},{}) length {} score {:.2} | aln cluster {} aln ({},{}) length {} score {:.2} | score diff {:.2}",
            gap.cluster_idx,
            gap.gap_start,
            gap.gap_end,
            gap.len(),
            gap.gap_score,
            aln.cluster_idx,
            aln.cluster_start,
            aln.cluster_end,
            aln.len(),
            aln.cluster_score,
            aln.cluster_score.0 - gap.gap_score.0
        );
    }
    let pairs: Vec<(usize, usize)> = selected.into_iter().collect();

    // Convert pairs to GapFill entries with coverage calculation
    let fills: Vec<GapFill> = pairs
        .into_iter()
        .map(|(gap_idx, aln_idx)| {
            let gap = &gaps[gap_idx];
            let aln = &alignments[aln_idx];
            let gap_len = gap.gap_end - gap.gap_start;
            let aln_len = aln.cluster_end - aln.cluster_start;
            let ratio = aln_len as f64 / gap_len as f64;
            let diff = gap_len as isize - aln_len as isize;

            let qual = gap.gap_score.0 - aln.cluster_score.0;

            log::debug!(
                "Gap fill: cluster {} gap ({},{} - {}bp, score {:.2}) filled by cluster {} aln ({},{} - {}bp, score {:.2}) | ratio {:.2}, diff {}, qual {:.2}",
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
                diff,
                qual
            );

            // Calculate coverage: fraction of gap covered by the alignment
            let gap_len = gap.gap_end - gap.gap_start;
            let overlap_start = gap
                .gap_start
                .max(aln.cluster_start.saturating_sub(tolerance));
            let overlap_end = gap.gap_end.min(aln.cluster_end.saturating_add(tolerance));
            let overlap = overlap_end.saturating_sub(overlap_start);
            let coverage = overlap as f64 / gap_len as f64;

                GapFill {
                    cluster_idx: gap.cluster_idx,
                    gap: (gap.gap_start, gap.gap_end),
                    gap_seed_idx: gap.gap_seed_idx,
                    filler_idx: aln.cluster_idx,
                    coverage,
                }
        })
        .collect();

    metrics::histogram!("analysis_gap_fills").record(start.elapsed().as_micros() as f64);

    fills
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Read {
    pub name: String,
    pub seq: String,
    pub qual: String,
}

impl Read {
    pub fn new(name: &str, seq: &[u8], qual: &[u8]) -> Self {
        Self {
            name: name.to_string(),
            seq: String::from_utf8_lossy(seq).to_string(),
            qual: String::from_utf8_lossy(qual).to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedSaver {
    pub read: Read,
    pub is_reverse: bool,
    pub seeds: Vec<SeedHit>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_seed_diagram() {
        // Create a chain with 3 seeds and gaps between them
        let chain = vec![
            SeedHit::new(0, 100, 0, 0, 1, 65),   // 65bp seed at ref 100, read 0
            SeedHit::new(0, 176, 73, 0, 1, 111), // 111bp seed at ref 176, read 73 (gap: qry=8, ref=11)
            SeedHit::new(0, 304, 203, 0, 1, 44), // 44bp seed at ref 304, read 203 (gap: qry=19, ref=17)
        ];

        let cluster = SeedCluster::new(chain, false, 10).unwrap();
        let (qry, ref_line) = cluster.format_seed_diagram();

        println!("{}", qry);
        println!("{}", ref_line);

        // Check that seeds are aligned (same width on both lines)
        assert!(qry.contains("[- 65bp -]"));
        assert!(ref_line.contains("[- 65bp -]"));
        assert!(qry.contains("[- 111bp -]"));
        assert!(ref_line.contains("[- 111bp -]"));
        assert!(qry.contains("[- 44bp -]"));
        assert!(ref_line.contains("[- 44bp -]"));

        // Check that gaps show different values
        assert!(qry.contains("8bp"));
        assert!(ref_line.contains("11bp"));
        assert!(qry.contains("19bp"));
        assert!(ref_line.contains("17bp"));

        // Both lines should have the same length (seeds aligned)
        assert_eq!(qry.len(), ref_line.len());
    }

    #[test]
    fn test_format_seed_diagram_single_seed() {
        let chain = vec![SeedHit::new(0, 100, 0, 0, 1, 50)];
        let cluster = SeedCluster::new(chain, false, 10).unwrap();
        let (qry, ref_line) = cluster.format_seed_diagram();

        assert_eq!(qry, "QRY: [- 50bp -]");
        assert_eq!(ref_line, "REF: [- 50bp -]");
    }

    #[test]
    fn test_format_seed_diagram_empty() {
        // An empty chain returns None from SeedCluster::new
        // but we can test the method directly by calling format_seed_diagram
        // on a cluster we construct manually
        let cluster = SeedCluster {
            read_start: 0,
            read_end: 0,
            chain: vec![],
            is_reverse: false,
            chrom_id: 0,
            gap_alignments: vec![],
        };
        let (qry, ref_line) = cluster.format_seed_diagram();

        assert_eq!(qry, "QRY:");
        assert_eq!(ref_line, "REF:");
    }
}
