use std::sync::Arc;

use crate::align::{
    Alignment, CigarOp, ContextAwareParams, ContextAwareScore, align, context_aware_score,
};
use crate::config;
use crate::error::{ParallaxError, Result};
use crate::index::Index;
use crate::kmers::Kmer;
use crate::reads::seeds::{SeedCluster, SeedHit, analyze_gap_fills};
use crate::reference::{ChromInfo, InMemoryReference};
use crate::utils::debug::{self, DebugFile};
use crate::utils::sequence::reverse_complement_into;
use crate::utils::{GroupByTrait, LongestSubsequence, dbscan_variance_aware};
use crate::writer::AlignmentWriter;

pub mod seeds;

/// SAM flags
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Compute the query (read) length consumed by a CIGAR string.
///
/// The query length is the sum of lengths for operations that consume query bases:
/// M, I, S, =, X (but not D, N, H, P).
///
/// Returns the query length, or 0 if the CIGAR is invalid or "*".
fn cigar_query_length(cigar: &str) -> usize {
    if cigar == "*" {
        return 0;
    }

    let mut len = 0usize;
    let mut num = 0usize;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as usize - '0' as usize);
        } else {
            // Operations that consume query: M, I, S, =, X
            // Operations that don't consume query: D, N, H, P
            match c {
                'M' | 'I' | 'S' | '=' | 'X' => len += num,
                'D' | 'N' | 'H' | 'P' => {} // Don't consume query
                _ => {}                     // Unknown op, ignore
            }
            num = 0;
        }
    }

    len
}

/// A candidate alignment with all necessary metadata for SAM output
#[derive(Clone)]
struct CandidateAlignment {
    chrom_id: usize,
    ref_start: usize,
    ref_end: usize,
    read_start: usize,
    read_end: usize,
    is_reverse: bool,
    alignment: Alignment,
    /// Context-aware score accounting for homopolymers, STRs, and sublinear gap extension
    context_score: ContextAwareScore,
}

impl CandidateAlignment {
    /// Calculate the fraction of the read covered by this alignment
    fn read_coverage(&self, read_len: usize) -> f64 {
        (self.read_end - self.read_start) as f64 / read_len as f64
    }

    /// Calculate edit distance (NM tag): mismatches + insertions + deletions
    fn edit_distance(&self) -> u32 {
        let mut nm = 0u32;
        for op in &self.alignment.cigar {
            match op {
                CigarOp::Mismatch(n) | CigarOp::Ins(n) | CigarOp::Del(n) => nm += n,
                CigarOp::Match(_) | CigarOp::SoftClip(_) => {}
            }
        }
        nm
    }

    /// Calculate alignment identity (matches / aligned length)
    /// Uses the pre-computed context_score for efficiency
    fn identity(&self) -> f64 {
        self.context_score.identity
    }

    /// Get the aligned length (excluding soft clips)
    fn aligned_length(&self) -> u32 {
        self.context_score.matches + self.context_score.mismatches + self.context_score.gap_bases
    }

    /// Calculate score per aligned base
    fn score_per_base(&self) -> f64 {
        let aligned = self.aligned_length();
        if aligned == 0 {
            f64::INFINITY
        } else {
            self.context_score.score as f64 / aligned as f64
        }
    }

    /// Calculate a minimap2-style alignment score for ranking
    /// Higher is better. Combines matches with penalties for errors.
    /// Uses: matches * match_bonus - mismatches * mismatch_penalty - gap_penalty
    fn ranking_score(&self) -> i64 {
        let matches = self.context_score.matches as i64;
        let mismatches = self.context_score.mismatches as i64;
        let gap_bases = self.context_score.gap_bases as i64;

        // Scoring: +2 per match, -4 per mismatch, -2 per gap base
        // This gives higher scores to longer, more accurate alignments
        matches * 2 - mismatches * 4 - gap_bases * 2
    }

    /// Calculate an information-theoretic alignment score based on CIGAR.
    ///
    /// This scoring function rewards consecutive match runs using N*log2(N+1),
    /// which reflects the statistical significance of contiguous matches vs
    /// scattered ones. Uses affine gap penalties.
    ///
    /// Parameters:
    /// - mismatch_penalty: penalty per mismatch base (e.g., 4.0)
    /// - gap_open: penalty for starting a gap (e.g., 6.0)
    /// - gap_extend: penalty per gap base (e.g., 1.0)
    ///
    /// Returns a score where higher is better.
    fn information_score(&self, mismatch_penalty: f64, gap_open: f64, gap_extend: f64) -> f64 {
        use crate::align::CigarOp;

        let mut score = 0.0;

        for op in &self.alignment.cigar {
            match op {
                CigarOp::Match(n) => {
                    // N * log2(N + 1) for match runs
                    // Using N+1 to handle N=1 gracefully (1 * log2(2) = 1)
                    let n = *n as f64;
                    score += n * (n + 1.0).log2();
                }
                CigarOp::Mismatch(n) => {
                    // Each mismatch is a discrete event that breaks a match run
                    score -= (*n as f64) * mismatch_penalty;
                }
                CigarOp::Ins(n) | CigarOp::Del(n) => {
                    // Affine gap penalty: open + extend * length
                    score -= gap_open + (*n as f64) * gap_extend;
                }
                CigarOp::SoftClip(_) => {
                    // Soft clips don't contribute to alignment score
                }
            }
        }

        score
    }

    /// Get read coordinates in forward orientation (for overlap detection).
    /// Strand coordinates are stored internally, but for comparing overlaps
    /// between alignments on different strands, we need forward coordinates.
    fn forward_read_coords(&self, seq_len: usize) -> (usize, usize) {
        if self.is_reverse {
            (seq_len - self.read_end, seq_len - self.read_start)
        } else {
            (self.read_start, self.read_end)
        }
    }
}

/// Classification of an alignment
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum AlignmentClass {
    Primary,
    Secondary,
    Supplementary,
    SecondarySupplementary, // Both 0x100 and 0x800
    LowQuality,
}

/// Classified alignment ready for SAM output
struct ClassifiedAlignment {
    candidate: CandidateAlignment,
    class: AlignmentClass,
    mapq: u8,
}

impl ClassifiedAlignment {
    fn sam_flag(&self) -> u16 {
        let mut flag = 0u16;
        if self.candidate.is_reverse {
            flag |= FLAG_REVERSE;
        }
        match self.class {
            AlignmentClass::Secondary => flag |= FLAG_SECONDARY,
            AlignmentClass::Supplementary => flag |= FLAG_SUPPLEMENTARY,
            AlignmentClass::SecondarySupplementary => flag |= FLAG_SECONDARY | FLAG_SUPPLEMENTARY,
            AlignmentClass::Primary | AlignmentClass::LowQuality => {}
        }
        flag
    }
}

/// Build a CandidateAlignment from a SeedCluster using pre-computed gap alignments.
///
/// This function uses the gap alignments already stored in the cluster (from
/// `align_gaps()`) rather than re-running WFA. If gap alignments haven't been
/// computed, it falls back to aligning gaps on demand.
///
/// Unlike `build_alignment_from_chain`, this function does NOT split at large gaps
/// since any necessary splitting should have already been done via
/// `split_at_failed_alignments()`.
///
/// # Arguments
/// * `read_id` - Read identifier for logging
/// * `cluster` - The seed cluster with pre-computed gap alignments
/// * `seq` - Read sequence (already reverse-complemented for reverse strand)
/// * `seq_len` - Length of the original read
/// * `reference` - Reference genome
///
/// # Returns
/// A CandidateAlignment if successful, None if the alignment couldn't be built.
fn build_alignment_from_cluster(
    read_id: &str,
    cluster: &SeedCluster,
    seq: &[u8],
    seq_len: usize,
    reference: &InMemoryReference,
) -> Option<CandidateAlignment> {
    let chain = &cluster.chain;
    if chain.is_empty() {
        return None;
    }

    let cfg = config::get();

    // Require either multiple seeds, or a single seed that's long enough
    if chain.len() == 1 && chain[0].match_len < cfg.seeding.min_single_seed_length {
        return None;
    }

    let chrom_id = cluster.chrom_id;
    let is_reverse = cluster.is_reverse;
    let mut full_cigar: Vec<CigarOp> = Vec::new();
    let mut total_score = 0i32;

    // Compute alignment span from first/last seeds
    let first = chain.first().unwrap();
    let last = chain.last().unwrap();

    let ref_start = first.ref_pos;
    let mut ref_end = last.ref_end();
    let read_start = first.read_pos;
    let mut read_end = last.read_end();

    // Add soft-clip for unaligned prefix
    if read_start > 0 {
        full_cigar.push(CigarOp::SoftClip(read_start as u32));
    }

    // Process seeds and gaps
    for (j, hit) in chain.iter().enumerate() {
        // Align gap before this seed (if not first seed)
        if j > 0 {
            let prev = &chain[j - 1];
            let read_gap_start = prev.read_end();
            let read_gap_end = hit.read_pos;
            let ref_gap_start = prev.ref_end();
            let ref_gap_end = hit.ref_pos;

            let read_gap_len = read_gap_end.saturating_sub(read_gap_start);
            let ref_gap_len = ref_gap_end.saturating_sub(ref_gap_start);

            if read_gap_len > 0 || ref_gap_len > 0 {
                // Try to use pre-computed gap alignment
                let gap_idx = j - 1;
                if let Some(gap_aln) = cluster.gap_alignment(gap_idx) {
                    // Use the pre-computed alignment
                    total_score += gap_aln.score;
                    full_cigar.extend(gap_aln.cigar.iter().copied());
                } else if read_gap_len > 0 && ref_gap_len > 0 {
                    // No pre-computed alignment - align on demand (fallback)
                    let ref_slice = reference.get_seq(chrom_id, ref_gap_start, ref_gap_end);
                    let read_slice = &seq[read_gap_start..read_gap_end];

                    if let Some(aln) = align(read_slice, ref_slice) {
                        total_score += aln.score;
                        full_cigar.extend(aln.cigar);
                    } else {
                        // Alignment failed - emit as I+D
                        full_cigar.push(CigarOp::Ins(read_gap_len as u32));
                        full_cigar.push(CigarOp::Del(ref_gap_len as u32));
                    }
                } else if read_gap_len > 0 {
                    full_cigar.push(CigarOp::Ins(read_gap_len as u32));
                } else if ref_gap_len > 0 {
                    full_cigar.push(CigarOp::Del(ref_gap_len as u32));
                }
            }
        }

        // Add the seed match
        full_cigar.push(CigarOp::Match(hit.match_len as u32));

        // Update endpoints
        read_end = hit.read_end();
        ref_end = hit.ref_end();
    }

    // Add soft-clip for unaligned suffix
    if read_end < seq_len {
        full_cigar.push(CigarOp::SoftClip((seq_len - read_end) as u32));
    }

    // Check we produced a valid CIGAR
    if full_cigar.is_empty()
        || !full_cigar
            .iter()
            .any(|op| matches!(op, CigarOp::Match(_) | CigarOp::Mismatch(_)))
    {
        return None;
    }

    let mut alignment = Alignment {
        score: total_score,
        cigar: full_cigar,
    };
    alignment.normalize();

    // Get the aligned portions for context-aware scoring
    let query_for_scoring = &seq[read_start..read_end];
    let ref_for_scoring = reference.get_seq(chrom_id, ref_start, ref_end);

    // Validate the alignment CIGAR against actual sequences
    if log::log_enabled!(log::Level::Debug) {
        if let Err(e) = alignment.validate(ref_for_scoring, seq, 0) {
            log::debug!(
                "Read {}: alignment validation issue at {}:{}-{}: {}",
                read_id,
                reference.chrom_name(chrom_id),
                ref_start,
                ref_end,
                e
            );
        }
    }

    // Compute context-aware score
    let params = ContextAwareParams::default();
    let context_score = context_aware_score(&alignment, ref_for_scoring, query_for_scoring, &params);

    Some(CandidateAlignment {
        chrom_id,
        ref_start,
        ref_end,
        read_start,
        read_end,
        is_reverse,
        alignment,
        context_score,
    })
}

/// Build alignment from a chain of seed matches, filling gaps with WFA.
///
/// The chain should be sorted by read position. Both sequences (read and reference)
/// are assumed to be in the same orientation - for reverse strand alignments,
/// the caller should pass the reverse-complemented read sequence.
///
/// When a gap exceeds `max_gap_length`, the alignment is split into multiple
/// separate alignments rather than emitting a dubious I+D operation.
///
/// # Arguments
/// * `read_id` - Read identifier for logging
/// * `chain` - Sorted chain of seed hits
/// * `seq` - Read sequence (already reverse-complemented for reverse strand)
/// * `seq_len` - Length of the original read
/// * `reference` - Reference genome
/// * `is_reverse` - Whether this is a reverse strand alignment (for marking in result)
///
/// # Returns
/// A vector of candidate alignments. Usually one, but may be multiple if the chain
/// was split at large gaps.
fn build_alignment_from_chain(
    read_id: &str,
    chain: &[SeedHit],
    seq: &[u8],
    seq_len: usize,
    reference: &InMemoryReference,
    is_reverse: bool,
) -> Vec<CandidateAlignment> {
    let cfg = config::get();

    if chain.is_empty() {
        return Vec::new();
    }

    // Require either multiple seeds, or a single seed that's long enough
    if chain.len() == 1 && chain[0].match_len < cfg.seeding.min_single_seed_length {
        return Vec::new();
    }

    // Split chain into segments at large gaps, then process each segment
    let segments = split_chain_at_large_gaps(chain, cfg.seeding.max_gap_length);

    let mut all_alignments = Vec::new();

    for (seg_idx, segment) in segments.iter().enumerate() {
        if segment.is_empty() {
            continue;
        }

        // For split chains, require each segment to have sufficient seed coverage
        if segments.len() > 1 {
            let segment_seed_len: usize = segment.iter().map(|h| h.match_len).sum();
            if segment_seed_len < cfg.seeding.min_single_seed_length {
                log::debug!(
                    "Read {}: skipping small segment {} with only {} bp of seeds",
                    read_id,
                    seg_idx,
                    segment_seed_len
                );
                continue;
            }
        }

        if let Some(aln) = build_alignment_from_segment(
            read_id,
            segment,
            seq,
            seq_len,
            reference,
            is_reverse,
            seg_idx == 0,                     // is_first_segment
            seg_idx == segments.len() - 1,    // is_last_segment
            cfg,
        ) {
            all_alignments.push(aln);
        }
    }

    all_alignments
}

/// Split a chain of seeds into segments at gaps that exceed max_gap_length.
///
/// Returns a vector of seed slices. Each slice represents a contiguous segment
/// that should become a separate alignment.
fn split_chain_at_large_gaps(chain: &[SeedHit], max_gap_length: usize) -> Vec<&[SeedHit]> {
    if chain.is_empty() {
        return Vec::new();
    }

    let mut segments = Vec::new();
    let mut segment_start = 0;

    for i in 1..chain.len() {
        let prev = &chain[i - 1];
        let curr = &chain[i];

        // Calculate gap in both read and reference space
        let read_gap = if curr.read_pos > prev.read_end() {
            curr.read_pos - prev.read_end()
        } else {
            0
        };
        let ref_gap = if curr.ref_pos > prev.ref_end() {
            curr.ref_pos - prev.ref_end()
        } else {
            0
        };

        let max_gap = read_gap.max(ref_gap);

        if max_gap > max_gap_length {
            // Split here: emit the segment up to (and including) prev
            if segment_start < i {
                segments.push(&chain[segment_start..i]);
            }
            segment_start = i;
        }
    }

    // Don't forget the final segment
    if segment_start < chain.len() {
        segments.push(&chain[segment_start..]);
    }

    segments
}

/// Build a single alignment from a segment of the chain.
///
/// This is the inner function that processes seeds without splitting.
/// Seeds are assumed to be non-overlapping (resolved during cluster creation).
/// Soft-clipping is adjusted based on whether this is the first/last segment.
fn build_alignment_from_segment(
    read_id: &str,
    segment: &[SeedHit],
    seq: &[u8],
    seq_len: usize,
    reference: &InMemoryReference,
    is_reverse: bool,
    _is_first_segment: bool,
    _is_last_segment: bool,
    cfg: &config::ParallaxConfig,
) -> Option<CandidateAlignment> {
    if segment.is_empty() {
        return None;
    }

    let chrom_id = segment[0].chrom_id;
    let mut full_cigar: Vec<CigarOp> = Vec::new();
    let mut total_score = 0i32;

    // Compute alignment span from first/last seeds (already sorted by read_pos)
    let first = segment.first().unwrap();
    let last = segment.last().unwrap();

    let ref_start = first.ref_pos;
    let mut ref_end = last.ref_end();
    let read_start = first.read_pos;
    let mut read_end = last.read_end();

    // Add soft-clip for unaligned prefix
    if read_start > 0 {
        full_cigar.push(CigarOp::SoftClip(read_start as u32));
    }

    // Process seeds - they are guaranteed non-overlapping after resolve_overlaps()
    for (j, hit) in segment.iter().enumerate() {
        // Align gap before this seed (if not first seed)
        if j > 0 {
            let prev = &segment[j - 1];
            let read_gap_start = prev.read_end();
            let read_gap_end = hit.read_pos;
            let ref_gap_start = prev.ref_end();
            let ref_gap_end = hit.ref_pos;

            // Gap lengths are guaranteed non-negative since seeds don't overlap
            let read_gap_len = read_gap_end - read_gap_start;
            let ref_gap_len = ref_gap_end - ref_gap_start;

            if read_gap_len > 0 && ref_gap_len > 0 {
                // Both have gaps - need to align
                // (Large gaps should have been handled by split_chain_at_large_gaps,
                // but we still check here as a safety measure)
                let max_gap = read_gap_len.max(ref_gap_len);

                if max_gap > cfg.seeding.max_gap_length {
                    log::warn!(
                        "Read {}: unexpected large gap within segment (read: {}, ref: {})",
                        read_id,
                        read_gap_len,
                        ref_gap_len
                    );
                    // Emit as I+D and continue (shouldn't happen normally)
                    full_cigar.push(CigarOp::Ins(read_gap_len as u32));
                    full_cigar.push(CigarOp::Del(ref_gap_len as u32));
                } else {
                    // Get reference and read slices
                    let ref_slice = reference.get_seq(chrom_id, ref_gap_start, ref_gap_end);
                    let read_slice = &seq[read_gap_start..read_gap_end];

                    if log::log_enabled!(log::Level::Debug) {
                        if read_slice.len() >= 150 || ref_slice.len() >= 150 {
                            log::debug!(
                                "Aligning read {} gap of size {} to ref gap of size {}: read pos {}-{}, ref pos {}-{}",
                                read_id,
                                read_slice.len(),
                                ref_slice.len(),
                                read_gap_start,
                                read_gap_end,
                                ref_gap_start,
                                ref_gap_end,
                            );
                        }
                    }
                    if let Some(aln) = align(read_slice, ref_slice) {
                        total_score += aln.score;
                        full_cigar.extend(aln.cigar);
                    } else {
                        // Alignment failed, emit as insertions/deletions
                        full_cigar.push(CigarOp::Ins(read_gap_len as u32));
                        full_cigar.push(CigarOp::Del(ref_gap_len as u32));
                    }
                }
            } else if read_gap_len > 0 {
                // Only read has gap - pure insertion
                full_cigar.push(CigarOp::Ins(read_gap_len as u32));
            } else if ref_gap_len > 0 {
                // Only reference has gap - pure deletion
                full_cigar.push(CigarOp::Del(ref_gap_len as u32));
            }
            // else: both zero - adjacent seeds, no gap to process
        }

        // Add the seed match
        full_cigar.push(CigarOp::Match(hit.match_len as u32));

        // Update endpoints (for the final seed)
        read_end = hit.read_end();
        ref_end = hit.ref_end();
    }

    // Add soft-clip for unaligned suffix
    if read_end < seq_len {
        full_cigar.push(CigarOp::SoftClip((seq_len - read_end) as u32));
    }

    // Check we actually produced a valid CIGAR
    if full_cigar.is_empty() || !full_cigar.iter().any(|op| matches!(op, CigarOp::Match(_) | CigarOp::Mismatch(_))) {
        return None;
    }

    let mut alignment = Alignment {
        score: total_score,
        cigar: full_cigar,
    };
    alignment.normalize();

    // Get the aligned portions for context-aware scoring
    let query_for_scoring = &seq[read_start..read_end];
    let ref_for_scoring = reference.get_seq(chrom_id, ref_start, ref_end);

    // Validate the alignment CIGAR against actual sequences
    if log::log_enabled!(log::Level::Debug) {
        if let Err(e) = alignment.validate(ref_for_scoring, seq, 0) {
            log::debug!(
                "Read {}: alignment validation issue at {}:{}-{}: {}",
                read_id,
                reference.chrom_name(chrom_id),
                ref_start,
                ref_end,
                e
            );
        }
    }

    // Compute context-aware score
    let params = ContextAwareParams::default();
    let context_score =
        context_aware_score(&alignment, ref_for_scoring, query_for_scoring, &params);

    Some(CandidateAlignment {
        chrom_id,
        ref_start,
        ref_end,
        read_start,
        read_end,
        is_reverse,
        alignment,
        context_score,
    })
}

/// Calculate minimap2-style MAPQ based on score ratio
///
/// MAPQ ≈ 40 * (1 - s2/s1) * min(1, aligned_len/100) * log2(s1)
/// Where s1 is best score and s2 is second-best score for overlapping region
fn compute_mapq(
    best_score: i64,
    second_best_score: Option<i64>,
    aligned_len: u32,
    identity: f64,
) -> u8 {
    if best_score <= 0 {
        return 0;
    }

    // Score ratio component: how much better is this than alternatives?
    let ratio = match second_best_score {
        Some(s2) if s2 > 0 => 1.0 - (s2 as f64 / best_score as f64),
        Some(_) => 1.0, // second best is non-positive, we're unique
        None => 1.0,    // no alternative, we're unique
    };

    // Length component: longer alignments get higher confidence
    let len_factor = (aligned_len as f64 / 100.0).min(1.0);

    // Score magnitude component: higher scores get higher confidence
    let score_factor = (best_score as f64).log2().max(1.0) / 10.0;

    // Identity component: better identity = higher confidence
    let identity_factor = identity;

    // Combine: base of 40, scaled by all factors
    let mapq = 40.0 * ratio * len_factor * score_factor * identity_factor;

    (mapq.round() as u8).min(60)
}

/// Classify candidate alignments into primary, secondary, supplementary, and low quality.
///
/// Classification rules per SAM spec:
/// 1. Group alignments by overlapping read regions (clusters)
/// 2. Primary: Best alignment overall (best score from best cluster)
/// 3. Secondary (0x100): Alternative mappings in the same cluster as primary
/// 4. Supplementary (0x800): Best alignment from a different cluster (chimeric)
/// 5. Secondary+Supplementary (0x100|0x800): Alternative mappings in a supplementary cluster
/// 6. Low Quality: Alignments below score/coverage/identity thresholds
///
/// ALT contig handling: When computing MAPQ, alignments to related chromosomes
/// (e.g., chr1 and chr1_KI270762v1_alt) are not treated as competing alignments
/// since they represent the same genomic location.
fn classify_alignments(
    candidates: Vec<CandidateAlignment>,
    read_len: usize,
    _chrom_info: &[ChromInfo],
) -> Vec<ClassifiedAlignment> {
    if candidates.is_empty() {
        return Vec::new();
    }

    let cfg = config::get();

    // Helper to check if an alignment passes quality thresholds
    let passes_quality = |c: &CandidateAlignment| -> bool {
        let coverage = c.read_coverage(read_len);
        let aligned_len = c.aligned_length();
        let passes_coverage = coverage >= cfg.filtering.min_read_coverage
            || aligned_len >= cfg.filtering.min_aligned_length;
        passes_coverage
            && c.identity() >= cfg.filtering.min_identity
            && c.score_per_base() <= cfg.filtering.max_score_per_base
    };

    // Filter to only quality alignments for set building
    let quality_indices: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| passes_quality(c))
        .map(|(i, _)| i)
        .collect();

    if quality_indices.is_empty() {
        // No quality alignments - mark all as low quality
        return candidates
            .into_iter()
            .map(|c| ClassifiedAlignment {
                mapq: 0,
                class: AlignmentClass::LowQuality,
                candidate: c,
            })
            .collect();
    }

    let mut candidates = candidates;
    candidates.sort_by(|a, b| a.read_start.cmp(&b.read_start));

    // Build read-covering sets using greedy algorithm
    // The DP approach doesn't work well when we allow some overlap between intervals
    let alignment_sets = build_covering_sets(
        &candidates,
        &quality_indices,
        cfg.classification.overlap_threshold,
        read_len,
        cfg.classification.set_gap_open,
        cfg.classification.set_gap_extend,
        cfg.classification.use_information_score,
        cfg.classification.info_mismatch_penalty,
        cfg.classification.info_gap_open,
        cfg.classification.info_gap_extend,
    );

    // Log alignment details with forward coordinates for debugging
    if log::log_enabled!(log::Level::Debug) {
        // Debug logging
        log::debug!(
            "Built {} read-covering sets from {} quality alignments:",
            alignment_sets.len(),
            quality_indices.len()
        );
        for (set_idx, set) in alignment_sets.iter().enumerate() {
            let mut indices = set.alignment_indices.clone();
            indices.sort_unstable();
            log::debug!(
                "  Set {}: score={:.1} coverage={:.1}% alignments={:?}",
                set_idx,
                set.total_score,
                set.read_coverage * 100.0,
                indices
            );
        }

        for (i, c) in candidates.iter().enumerate() {
            let in_sets: Vec<usize> = alignment_sets
                .iter()
                .enumerate()
                .filter(|(_, s)| s.alignment_indices.contains(&i))
                .map(|(idx, _)| idx)
                .collect();
            let (fwd_start, fwd_end) = c.forward_read_coords(read_len);
            log::debug!(
                "  Alignment {}: read [{}, {}] fwd [{}, {}] (len {}) ref {}:[{}, {}] (len {}) strand={} score={:.1} in_sets={:?} (M={} X={} gaps={} id={:.1}%)",
                i,
                c.read_start,
                c.read_end,
                fwd_start,
                fwd_end,
                c.read_end - c.read_start,
                c.chrom_id,
                c.ref_start,
                c.ref_end,
                c.ref_end - c.ref_start,
                if c.is_reverse { "-" } else { "+" },
                c.information_score(
                    cfg.classification.info_mismatch_penalty,
                    cfg.classification.info_gap_open,
                    cfg.classification.info_gap_extend
                ),
                in_sets,
                c.context_score.matches,
                c.context_score.mismatches,
                c.context_score.gap_bases,
                c.identity() * 100.0
            );
        }
    }

    // The best set determines primary/supplementary
    let best_set = &alignment_sets[0];
    let primary_idx = best_set.alignment_indices[0]; // Best alignment in best set

    // Build set membership for classification
    let mut alignment_to_best_set: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    for (set_idx, set) in alignment_sets.iter().enumerate() {
        for &aln_idx in &set.alignment_indices {
            alignment_to_best_set.entry(aln_idx).or_insert(set_idx);
        }
    }

    // Collect scores for MAPQ calculation
    let mut set_scores: Vec<f64> = alignment_sets.iter().map(|s| s.total_score).collect();
    set_scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let second_best_set_score = set_scores.get(1).copied();

    let mut classified = Vec::with_capacity(candidates.len());

    for (i, candidate) in candidates.into_iter().enumerate() {
        if !passes_quality(&candidate) {
            classified.push(ClassifiedAlignment {
                mapq: 0,
                class: AlignmentClass::LowQuality,
                candidate,
            });
            continue;
        }

        let score = candidate.ranking_score();
        let aligned_len = candidate.aligned_length();
        let identity = candidate.identity();
        let in_best_set = best_set.alignment_indices.contains(&i);
        let is_primary = i == primary_idx;

        if is_primary {
            // Primary alignment - best in the best set
            // MAPQ should reflect confidence that this set is correct vs alternatives
            // With information-theoretic scores, the DIFFERENCE in scores is meaningful:
            // score_diff represents log-likelihood ratio, so:
            // - error_prob ≈ 1 / (1 + 2^score_diff)
            // - MAPQ = -10 * log10(error_prob) ≈ score_diff * 10 * log10(2) for large diffs
            // For small diffs, use the exact formula
            let mapq = match second_best_set_score {
                Some(s2) if best_set.total_score > 0.0 => {
                    let score_diff = best_set.total_score - s2;
                    // Convert score difference to error probability
                    // P(error) = 1 / (1 + 2^(diff/20))
                    let error_prob = 1.0 / (1.0 + (2.0_f64).powf(score_diff / 20.0));
                    // MAPQ = -10 * log10(error_prob), capped at 60
                    let raw_mapq = -10.0 * error_prob.log10();
                    (raw_mapq.round() as u8).min(60)
                }
                _ => 60, // No second-best set, unique mapping
            };
            classified.push(ClassifiedAlignment {
                mapq,
                class: AlignmentClass::Primary,
                candidate,
            });
        } else if in_best_set {
            // Other alignments in best set -> Supplementary (chimeric pieces)
            // These share the same MAPQ confidence as primary since they're part of the same solution
            let mapq = match second_best_set_score {
                Some(s2) if best_set.total_score > 0.0 => {
                    let score_diff = best_set.total_score - s2;
                    let error_prob = 1.0 / (1.0 + (2.0_f64).powf(score_diff / 20.0));
                    let raw_mapq = -10.0 * error_prob.log10();
                    (raw_mapq.round() as u8).min(60)
                }
                _ => 60,
            };
            classified.push(ClassifiedAlignment {
                mapq,
                class: AlignmentClass::Supplementary,
                candidate,
            });
        } else {
            // Not in best set - check if it's the best alignment for its read region
            let my_set_idx = alignment_to_best_set.get(&i).copied().unwrap_or(usize::MAX);
            let is_best_in_my_set = alignment_sets
                .get(my_set_idx)
                .map(|s| s.alignment_indices.first() == Some(&i))
                .unwrap_or(false);

            if is_best_in_my_set {
                // Best in an alternative set -> Secondary (alternative mapping)
                let my_set_score = alignment_sets.get(my_set_idx).map(|s| s.total_score as i64);
                let mapq = compute_mapq(score, my_set_score, aligned_len, identity);
                classified.push(ClassifiedAlignment {
                    mapq,
                    class: AlignmentClass::Secondary,
                    candidate,
                });
            } else {
                // Not best in any set -> Secondary+Supplementary
                let my_set_score = alignment_sets.get(my_set_idx).map(|s| s.total_score as i64);
                let mapq = compute_mapq(score, my_set_score, aligned_len, identity);
                classified.push(ClassifiedAlignment {
                    mapq,
                    class: AlignmentClass::SecondarySupplementary,
                    candidate,
                });
            }
        }
    }

    // Sort so primary comes first, then supplementary, then secondary, then secondary+supplementary
    classified.sort_by_key(|c| match c.class {
        AlignmentClass::Primary => 0,
        AlignmentClass::Supplementary => 1,
        AlignmentClass::Secondary => 2,
        AlignmentClass::SecondarySupplementary => 3,
        AlignmentClass::LowQuality => 4,
    });

    classified
}

/// A set of non-overlapping alignments that together cover the read
#[derive(Debug)]
struct AlignmentSet {
    /// Indices into the candidates array, sorted by score (best first)
    alignment_indices: Vec<usize>,
    /// Combined score for the set (may be f64 for information-theoretic scoring)
    total_score: f64,
    /// Fraction of read covered by this set
    read_coverage: f64,
}

/// Build read-covering sets using a greedy algorithm
///
/// For each starting alignment, greedily add non-overlapping alignments
/// to maximize coverage. Score the resulting set and keep track of
/// unique sets.
///
/// When `use_info_score` is true, uses information-theoretic scoring (N*log2(N) for match runs).
/// Otherwise uses the traditional linear scoring.
#[allow(dead_code)]
fn build_covering_sets(
    candidates: &[CandidateAlignment],
    quality_indices: &[usize],
    overlap_threshold: f64,
    seq_len: usize,
    gap_open: i64,
    gap_extend: i64,
    use_info_score: bool,
    info_mismatch: f64,
    info_gap_open: f64,
    info_gap_extend: f64,
) -> Vec<AlignmentSet> {
    if quality_indices.is_empty() {
        return Vec::new();
    }

    // Scoring function - either information-theoretic or traditional
    let score_alignment = |c: &CandidateAlignment| -> f64 {
        if use_info_score {
            c.information_score(info_mismatch, info_gap_open, info_gap_extend)
        } else {
            c.ranking_score() as f64
        }
    };

    // Helper to check if two alignments can coexist in the same set
    // They must not significantly overlap in read coordinates (using forward coords)
    let can_coexist = |i: usize, j: usize| -> bool {
        let ci = &candidates[i];
        let cj = &candidates[j];

        let (ci_start, ci_end) = ci.forward_read_coords(seq_len);
        let (cj_start, cj_end) = cj.forward_read_coords(seq_len);

        let overlap_start = ci_start.max(cj_start);
        let overlap_end = ci_end.min(cj_end);

        if overlap_start >= overlap_end {
            return true; // No overlap
        }

        let overlap_len = (overlap_end - overlap_start) as f64;
        let len_i = (ci_end - ci_start) as f64;
        let len_j = (cj_end - cj_start) as f64;

        // They can coexist if overlap is small for both
        overlap_len / len_i <= overlap_threshold && overlap_len / len_j <= overlap_threshold
    };

    // Sort quality indices by score (descending)
    let mut sorted_indices = quality_indices.to_vec();
    sorted_indices.sort_by(|&a, &b| {
        score_alignment(&candidates[b])
            .partial_cmp(&score_alignment(&candidates[a]))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Build sets using greedy algorithm starting from each alignment
    let mut all_sets: Vec<AlignmentSet> = Vec::new();
    let mut seen_set_signatures: std::collections::HashSet<Vec<usize>> =
        std::collections::HashSet::new();

    for &start_idx in &sorted_indices {
        // Build a set starting from this alignment
        let mut set_indices = vec![start_idx];
        let mut raw_score = score_alignment(&candidates[start_idx]);

        // Greedily add compatible alignments in score order
        for &candidate_idx in &sorted_indices {
            if set_indices.contains(&candidate_idx) {
                continue;
            }

            // Check if this candidate can coexist with all current set members
            let compatible = set_indices
                .iter()
                .all(|&existing| can_coexist(existing, candidate_idx));

            if compatible {
                set_indices.push(candidate_idx);
                raw_score += score_alignment(&candidates[candidate_idx]);
            }
        }

        // Sort set indices for consistent signature
        set_indices.sort();

        // Skip if we've already seen this exact set
        if seen_set_signatures.contains(&set_indices) {
            continue;
        }
        seen_set_signatures.insert(set_indices.clone());

        // Calculate read coverage and gap penalty for this set
        let read_coverage = calculate_set_coverage(candidates, &set_indices, seq_len);
        let (num_breaks, uncovered_bases) = calculate_set_gaps(candidates, &set_indices, seq_len);

        // Apply affine gap penalty: raw_score - gap_open * breaks - gap_extend * uncovered
        let gap_penalty =
            (gap_open * num_breaks as i64 + gap_extend * uncovered_bases as i64) as f64;
        let total_score = raw_score - gap_penalty;

        // Sort by score (best first) for the final set
        set_indices.sort_by(|&a, &b| {
            score_alignment(&candidates[b])
                .partial_cmp(&score_alignment(&candidates[a]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        all_sets.push(AlignmentSet {
            alignment_indices: set_indices,
            total_score,
            read_coverage,
        });
    }

    // Sort sets by total score (descending), with coverage as tiebreaker
    all_sets.sort_by(|a, b| {
        b.total_score
            .partial_cmp(&a.total_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.read_coverage
                    .partial_cmp(&a.read_coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    all_sets
}

/// Calculate the fraction of read covered by a set of alignments
fn calculate_set_coverage(
    candidates: &[CandidateAlignment],
    indices: &[usize],
    seq_len: usize,
) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }

    // Merge overlapping intervals to get total covered bases
    // Use forward coordinates for consistent overlap calculation
    let mut intervals: Vec<(usize, usize)> = indices
        .iter()
        .map(|&i| candidates[i].forward_read_coords(seq_len))
        .collect();
    intervals.sort_by_key(|&(start, _)| start);

    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        } else {
            merged.push((start, end));
        }
    }

    let covered: usize = merged.iter().map(|(s, e)| e - s).sum();

    // Coverage is fraction of read length covered
    if seq_len == 0 {
        0.0
    } else {
        covered as f64 / seq_len as f64
    }
}

/// Calculate the affine gap penalty for a set of alignments.
///
/// Returns (num_breaks, uncovered_bases) where:
/// - num_breaks: number of gaps between alignment intervals (alignment_count - 1 for non-overlapping)
/// - uncovered_bases: total bases in the read not covered by any alignment
///
/// Uses forward coordinates for consistent calculation across strands.
fn calculate_set_gaps(
    candidates: &[CandidateAlignment],
    indices: &[usize],
    seq_len: usize,
) -> (usize, usize) {
    if indices.is_empty() {
        return (0, seq_len);
    }

    // Merge overlapping intervals to count breaks and uncovered bases
    // Use forward coordinates for consistent overlap calculation
    let mut intervals: Vec<(usize, usize)> = indices
        .iter()
        .map(|&i| candidates[i].forward_read_coords(seq_len))
        .collect();
    intervals.sort_by_key(|&(start, _)| start);

    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        } else {
            merged.push((start, end));
        }
    }

    // Number of breaks = number of merged intervals - 1 (gaps between them)
    let num_breaks = if merged.len() > 1 {
        merged.len() - 1
    } else {
        0
    };

    // Uncovered bases = seq_len - covered bases
    let covered: usize = merged.iter().map(|(s, e)| e - s).sum();
    let uncovered = seq_len.saturating_sub(covered);

    (num_breaks, uncovered)
}

/// Build read-covering sets using weighted interval scheduling DP
///
/// This is an optimal algorithm for finding the maximum-weight set of
/// non-overlapping intervals. It uses dynamic programming on intervals
/// sorted by end position.
///
/// Time complexity: O(n log n) for the optimal set, O(kn log n) for k sets
///
/// Returns sets in order of decreasing total score.
#[allow(dead_code)]
fn build_covering_sets_dp(
    candidates: &[CandidateAlignment],
    quality_indices: &[usize],
    overlap_threshold: f64,
    seq_len: usize,
) -> Vec<AlignmentSet> {
    if quality_indices.is_empty() {
        return Vec::new();
    }

    // Helper to check if two alignments can coexist (same logic as greedy)
    // Uses forward coordinates so overlaps are computed consistently regardless of strand
    let can_coexist = |i: usize, j: usize| -> bool {
        let ci = &candidates[i];
        let cj = &candidates[j];

        let (ci_start, ci_end) = ci.forward_read_coords(seq_len);
        let (cj_start, cj_end) = cj.forward_read_coords(seq_len);

        let overlap_start = ci_start.max(cj_start);
        let overlap_end = ci_end.min(cj_end);

        if overlap_start >= overlap_end {
            return true; // No overlap
        }

        let overlap_len = (overlap_end - overlap_start) as f64;
        let len_i = (ci_end - ci_start) as f64;
        let len_j = (cj_end - cj_start) as f64;

        overlap_len / len_i <= overlap_threshold && overlap_len / len_j <= overlap_threshold
    };

    // Sort indices by read_end position in forward coordinates (required for interval scheduling DP)
    let mut sorted_by_end: Vec<usize> = quality_indices.to_vec();
    sorted_by_end.sort_by_key(|&i| candidates[i].forward_read_coords(seq_len).1);

    let n = sorted_by_end.len();

    // For each alignment i, find the rightmost alignment j where end_j <= start_i
    // (i.e., j is compatible with i and ends before i starts)
    // We use binary search for O(log n) per query
    let find_last_compatible = |idx: usize| -> Option<usize> {
        let (start_i, _) = candidates[sorted_by_end[idx]].forward_read_coords(seq_len);
        // Binary search for largest j where candidates[sorted_by_end[j]].read_end <= start_i
        // and can_coexist is satisfied
        let mut best: Option<usize> = None;
        for j in (0..idx).rev() {
            let cand_j = sorted_by_end[j];
            let (_, end_j) = candidates[cand_j].forward_read_coords(seq_len);
            if end_j <= start_i && can_coexist(sorted_by_end[idx], cand_j) {
                best = Some(j);
                break;
            }
            // Check overlap threshold compatibility even if there's some overlap
            if can_coexist(sorted_by_end[idx], cand_j) {
                // This one is compatible, but might not be the best
                // Keep looking for a truly non-overlapping one
                if best.is_none() {
                    best = Some(j);
                }
            }
        }
        best
    };

    // DP arrays:
    // dp[i] = maximum score achievable using alignments from 0..=i
    // choice[i] = true if we include alignment i in the optimal solution
    let mut dp: Vec<i64> = vec![0; n];
    let mut choice: Vec<bool> = vec![false; n];
    let mut prev: Vec<Option<usize>> = vec![None; n]; // For backtracking

    for i in 0..n {
        let score_i = candidates[sorted_by_end[i]].ranking_score();
        let prev_best = if i > 0 { dp[i - 1] } else { 0 };

        // Option 1: Don't include alignment i
        let exclude_score = prev_best;

        // Option 2: Include alignment i
        let include_score = if let Some(j) = find_last_compatible(i) {
            score_i + dp[j]
        } else {
            score_i
        };

        if include_score >= exclude_score {
            dp[i] = include_score;
            choice[i] = true;
            prev[i] = find_last_compatible(i);
        } else {
            dp[i] = exclude_score;
            choice[i] = false;
            prev[i] = if i > 0 { Some(i - 1) } else { None };
        }
    }

    // Backtrack to find the optimal set
    let mut all_sets: Vec<AlignmentSet> = Vec::new();
    let mut used: Vec<bool> = vec![false; n];

    // Extract multiple sets by repeatedly finding optimal over unused alignments
    loop {
        // Find optimal set among unused alignments
        let mut best_set_indices: Vec<usize> = Vec::new();
        let mut i = n;

        // Find the last unused alignment
        while i > 0 {
            i -= 1;
            if !used[i] {
                break;
            }
        }
        if i == 0 && used[0] {
            break; // All alignments used
        }

        // Recompute DP only over unused alignments
        let unused_indices: Vec<usize> = (0..n).filter(|&j| !used[j]).collect();
        if unused_indices.is_empty() {
            break;
        }

        // Simple DP over unused alignments
        let m = unused_indices.len();
        let mut dp2: Vec<i64> = vec![0; m];
        let mut choice2: Vec<bool> = vec![false; m];

        for ii in 0..m {
            let orig_i = unused_indices[ii];
            let score_i = candidates[sorted_by_end[orig_i]].ranking_score();
            let prev_best = if ii > 0 { dp2[ii - 1] } else { 0 };

            // Find last compatible among unused
            let mut last_compat: Option<usize> = None;
            for jj in (0..ii).rev() {
                let orig_j = unused_indices[jj];
                if can_coexist(sorted_by_end[orig_i], sorted_by_end[orig_j]) {
                    let (_, end_j) = candidates[sorted_by_end[orig_j]].forward_read_coords(seq_len);
                    let (start_i, _) =
                        candidates[sorted_by_end[orig_i]].forward_read_coords(seq_len);
                    if end_j <= start_i {
                        last_compat = Some(jj);
                        break;
                    }
                    if last_compat.is_none() {
                        last_compat = Some(jj);
                    }
                }
            }

            let include_score = if let Some(jj) = last_compat {
                score_i + dp2[jj]
            } else {
                score_i
            };

            if include_score >= prev_best {
                dp2[ii] = include_score;
                choice2[ii] = true;
            } else {
                dp2[ii] = prev_best;
                choice2[ii] = false;
            }
        }

        // Backtrack to extract the set
        let mut ii = m;
        while ii > 0 {
            ii -= 1;
            if choice2[ii] {
                let orig_i = unused_indices[ii];
                best_set_indices.push(sorted_by_end[orig_i]);
                used[orig_i] = true;

                // Jump to last compatible
                let (start_i, _) = candidates[sorted_by_end[orig_i]].forward_read_coords(seq_len);
                while ii > 0 {
                    ii -= 1;
                    let orig_j = unused_indices[ii];
                    let (_, end_j) = candidates[sorted_by_end[orig_j]].forward_read_coords(seq_len);
                    if end_j <= start_i && can_coexist(sorted_by_end[orig_i], sorted_by_end[orig_j])
                    {
                        break;
                    }
                }
            }
        }

        if best_set_indices.is_empty() {
            break;
        }

        // Calculate set metrics
        let total_score: f64 = best_set_indices
            .iter()
            .map(|&i| candidates[i].ranking_score() as f64)
            .sum();
        let read_coverage = calculate_set_coverage(candidates, &best_set_indices, seq_len);

        // Sort by score (best first) for the final set
        best_set_indices.sort_by(|&a, &b| {
            candidates[b]
                .ranking_score()
                .cmp(&candidates[a].ranking_score())
        });

        all_sets.push(AlignmentSet {
            alignment_indices: best_set_indices,
            total_score,
            read_coverage,
        });

        // Continue to find more sets
    }

    // Sort sets by total score (descending)
    all_sets.sort_by(|a, b| {
        b.total_score
            .partial_cmp(&a.total_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.read_coverage
                    .partial_cmp(&a.read_coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    all_sets
}

/// Collector for seed clusters with reusable buffers.
///
/// This struct holds all the intermediate buffers needed for seeding,
/// merging, extension, and clustering. Reusing these buffers across
/// multiple calls avoids repeated allocation.
struct ClusterCollector {
    /// Seed hits collected from k-mer index lookups
    hits: Vec<SeedHit>,
    /// Temporary buffer for index lookups
    hit_vec: Vec<(usize, usize)>,
    /// Scratch space for merging/deduplication
    merge_scratch: Vec<SeedHit>,
    /// DBSCAN cluster boundaries
    cuts: Vec<usize>,
    /// LIS computation helper
    longest_subsequence: LongestSubsequence,
    /// Indices of seeds in LIS chain
    chain_indices: Vec<usize>,
}

impl ClusterCollector {
    /// Create a new collector with empty buffers
    fn new() -> Self {
        ClusterCollector {
            hits: Vec::new(),
            hit_vec: Vec::new(),
            merge_scratch: Vec::new(),
            cuts: Vec::new(),
            longest_subsequence: LongestSubsequence::default(),
            chain_indices: Vec::new(),
        }
    }

    /// Collect seed clusters from a single strand.
    ///
    /// This performs seeding, merging, extension, and DBSCAN clustering, returning
    /// the resulting clusters without building alignments. This separation allows
    /// for cross-strand analysis before alignment construction.
    fn collect_from_strand<const K: usize, const S: usize>(
        &mut self,
        strand_seq: &[u8],
        strand_qual: Option<&[u8]>,
        is_reverse: bool,
        index: &Index<K, S>,
        reference: &InMemoryReference,
        read_name: &str,
    ) -> Vec<SeedCluster> {
        let cfg = config::get();
        let seq_len = strand_seq.len();
        let max_var = (seq_len as f64 * cfg.seeding.variance_coef).powi(2);

        self.hits.clear();

        // Phase 1: Collect seed hits using forward-only syncmers
        Kmer::<K>::kmerize_open_syncmers_fwd(strand_seq, [(); S], |pos, kmer| {
            self.hit_vec.clear();
            index.with(&kmer, |chrom_id, chrom_pos| {
                self.hit_vec.push((chrom_id, chrom_pos));
            });
            let kmer_uniqueness = self.hit_vec.len() as u32;
            // Use seeds up to occurrence threshold
            if self.hit_vec.len() <= cfg.seeding.max_seed_occurrences {
                for &(chrom_id, chrom_pos) in self.hit_vec.iter() {
                    self.hits
                        .push(SeedHit::new(chrom_id, chrom_pos, pos, kmer.0, kmer_uniqueness, K));
                }
            }
        });

        let strand_name = if is_reverse { "REV" } else { "FWD" };
        metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase()))
            .record(self.hits.len() as f64);

        // Phase 2: Sort hits - SeedHit's Ord gives us (chrom_id, diagonal, ref_pos) order
        self.hits.sort_unstable();

        // Phase 3: Merge overlapping/adjacent hits on same diagonal
        self.merge_scratch.clear();
        for hit in self.hits.drain(..) {
            if let Some(last) = self.merge_scratch.last_mut() {
                if last
                    .extend(hit.chrom_id, hit.ref_pos, hit.read_pos, hit.kmer, hit.kmer_uniqueness, K)
                    .is_none()
                {
                    continue; // Successfully merged
                }
            }
            self.merge_scratch.push(hit);
        }
        std::mem::swap(&mut self.hits, &mut self.merge_scratch);

        // Phase 3b: Extend each seed's exact match bidirectionally
        // This is the minimap2-style extension that maximizes anchor length
        for hit in self.hits.iter_mut() {
            let ref_seq = reference.get_seq(hit.chrom_id, 0, usize::MAX);
            hit.extend_exact(strand_seq, ref_seq);
        }

        // Phase 3c: Remove duplicates created by extension
        self.merge_scratch.clear();
        for hit in self.hits.drain(..) {
            if let Some(last) = self.merge_scratch.last() {
                if hit.chrom_id == last.chrom_id
                    && hit.diagonal == last.diagonal
                    && hit.ref_pos == last.ref_pos
                    && hit.match_len == last.match_len
                {
                    continue; // Duplicate, skip
                }
            }
            self.merge_scratch.push(hit);
        }
        std::mem::swap(&mut self.hits, &mut self.merge_scratch);

        // Write debug SAM output for seed hits (if debug file is initialized)
        if debug::is_enabled(DebugFile::Seeds) {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                debug::write(DebugFile::Seeds, &hit.to_sam_line(
                    read_name,
                    chrom_name,
                    is_reverse,
                    strand_seq,
                    strand_qual,
                ));
            }
        }
        // Write debug TSV output for seed hits (if debug file is initialized)
        if debug::is_enabled(DebugFile::Alignments) {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                let strand = if is_reverse { "-" } else { "+" };
                // Convert strand coordinates to forward coordinates
                let (fwd_start, fwd_end) = if is_reverse {
                    (seq_len - hit.read_end(), seq_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                debug::write(DebugFile::Alignments, &format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_name,
                    fwd_start,
                    fwd_end,
                    seq_len,
                    chrom_name,
                    hit.ref_pos,
                    hit.ref_end(),
                    strand,
                    hit.match_len,
                ));
            }
        }

        // Phase 4 & 5: Cluster hits by diagonal using DBSCAN, then build LIS chains.
        // Important: We must partition by chromosome first, since hits from different
        // chromosomes should never be clustered together. Hits are already sorted by
        // (chrom_id, diagonal, ref_pos), so we find chromosome boundaries and process
        // each partition separately.

        let mut clusters = Vec::new();

        // Process each chromosome partition
        for (_chrom_id, partition) in self.hits.group_by(|seed| seed.chrom_id) {
            if partition.is_empty() {
                continue;
            }

            // Run DBSCAN on this chromosome's hits
            self.cuts.clear();
            dbscan_variance_aware(
                partition,
                cfg.seeding.min_seed_cluster_distance,
                max_var,
                |hit| hit.diagonal,
                &mut self.cuts,
            );

            // Build LIS chains for each cluster within this partition
            for i in 1..self.cuts.len() {
                let begin = self.cuts[i - 1];
                let end = self.cuts[i];

                // Sort cluster hits by ref_pos for LIS.
                // DBSCAN groups by diagonal, but within a cluster, different diagonals
                // may interleave in ref_pos order. Sorting by ref_pos ensures we can
                // find the true longest colinear chain.
                let mut cluster_hits: Vec<SeedHit> = partition[begin..end].to_vec();
                cluster_hits.sort_unstable_by_key(|h| h.ref_pos);

                if cluster_hits.len() == 1 {
                    let seed = &cluster_hits[0];
                    if seed.match_len < cfg.seeding.min_single_seed_length {
                        continue; // Skip tiny single-seed clusters
                    }
                }

                if log::log_enabled!(log::Level::Debug) {
                    log::debug!(
                        "  Cluster {}:{}-{} {} ({} seeds) before LIS:",
                        reference.chrom_name(cluster_hits[0].chrom_id),
                        begin,
                        end,
                        (if is_reverse { "-" } else { "+" }),
                        cluster_hits.len()
                    );
                    for seed in cluster_hits.iter() {
                        log::debug!(
                            "    Seed: chrom {}, ref {}-{}, read {}-{}, len {}",
                            reference.chrom_name(seed.chrom_id),
                            seed.ref_pos,
                            seed.ref_end(),
                            seed.read_pos,
                            seed.read_end(),
                            seed.match_len,
                        );
                    }
                }

                // Use LIS on read_pos to ensure the chain is colinear.
                // Since we've sorted by ref_pos, the LIS on read_pos will select
                // seeds that are colinear in both reference and read space.
                self.longest_subsequence.longest_colinear_chain(
                    &cluster_hits,
                    |hit| hit.read_pos as i64,
                    true,
                    &mut self.chain_indices,
                );

                let mut chain: Vec<SeedHit> = self
                    .chain_indices
                    .iter()
                    .map(|&i| cluster_hits[i])
                    .collect();
                chain.sort_by_key(|h| h.fwd_read_range(seq_len, is_reverse));

                metrics::histogram!(format!("{}_chain_length", strand_name.to_lowercase()))
                    .record(chain.len() as f64);

                // Minimum seed length after overlap truncation is K/2
                let min_seed_length = K / 2;

                if let Some(cluster) = SeedCluster::new(chain, is_reverse, min_seed_length) {
                    // Compute chain score with gap penalties
                    let score = cluster.chain_score(
                        cfg.seeding.gap_penalty_linear,
                        cfg.seeding.gap_penalty_log,
                    );

                    // Filter by minimum chain score
                    if score < cfg.seeding.min_chain_score {
                        log::debug!(
                            "  Skipping cluster with low chain score: {:.1} < {:.1} (seeds: {}, total_len: {})",
                            score,
                            cfg.seeding.min_chain_score,
                            cluster.chain.len(),
                            cluster.total_seed_length()
                        );
                        continue;
                    }

                    clusters.push(cluster);
                }
            }
        }

        clusters
    }
}

/// Build alignments from collected seed clusters.
///
/// This converts SeedClusters into CandidateAlignments by running WFA on gaps.
/// Clusters may be split into multiple alignments if they contain large gaps.
#[allow(dead_code)]
fn build_alignments_from_clusters(
    clusters: &[SeedCluster],
    read_name: &str,
    fwd_seq: &[u8],
    rc_seq: &[u8],
    seq_len: usize,
    reference: &InMemoryReference,
) -> Vec<CandidateAlignment> {
    let mut candidates = Vec::new();

    for cluster in clusters {
        let strand_seq = if cluster.is_reverse { rc_seq } else { fwd_seq };

        let alignments = build_alignment_from_chain(
            read_name,
            &cluster.chain,
            strand_seq,
            seq_len,
            reference,
            cluster.is_reverse,
        );
        candidates.extend(alignments);
    }

    candidates
}

/// Align a single read and output SAM record(s) using the provided writer.
///
/// This is the core alignment function that can be called from FASTQ or uBAM readers.
///
/// # Arguments
/// * `index` - The k-mer index for seed lookup
/// * `reference` - The reference genome
/// * `writer` - The AlignmentWriter to output SAM records
/// * `read_name` - Name of the read
/// * `seq` - Read sequence (forward strand)
/// * `qual` - Quality scores (same orientation as seq), or None if unavailable
pub fn align_read<const K: usize, const S: usize, W: std::io::Write>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    writer: &AlignmentWriter<W>,
    read_name: &str,
    seq: &[u8],
    qual: Option<&[u8]>,
) {
    let seq_len = seq.len();

    // Reusable cluster collector
    let mut collector = ClusterCollector::new();

    // Compute reverse complement for reverse strand processing
    let mut rc_seq = Vec::with_capacity(seq_len);
    reverse_complement_into(seq, &mut rc_seq);

    // Reverse quality scores for reverse strand (if available)
    let rc_qual: Option<Vec<u8>> = qual.map(|q| q.iter().rev().copied().collect());

    // =========================================================================
    // PASS 1: Collect all seed clusters from both strands
    // =========================================================================

    let mut all_clusters: Vec<SeedCluster> = Vec::new();

    // Collect clusters from forward strand
    let fwd_clusters = collector.collect_from_strand(seq, qual, false, index, reference, read_name);
    all_clusters.extend(fwd_clusters);

    // Collect clusters from reverse strand
    let rev_clusters = collector.collect_from_strand(
        &rc_seq,
        rc_qual.as_deref(),
        true,
        index,
        reference,
        read_name,
    );
    all_clusters.extend(rev_clusters);

    all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));

    log::info!(
        "Read {}: collected {} seed clusters from both strands (coverage {:.2}%)",
        read_name,
        all_clusters.len(),
        all_clusters.iter().map(|c| c.read_coverage(seq_len)).sum::<f64>() * 100.0,
    );

    // =========================================================================
    // PASS 1.5: Align gaps and split at failed alignments
    // =========================================================================
    // For each cluster, run WFA on all gaps. If any gap alignment fails (None),
    // split the cluster at that point. This must happen before gap-fill analysis
    // so we only consider clusters with valid internal alignments.

    let cfg = config::get();
    let min_seed_length = K / 2;

    let mut new_clusters_from_splits = Vec::new();
    for cluster in all_clusters.iter_mut() {
        let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
        let chrom_len = reference.chrom_length(cluster.chrom_id) as usize;
        let ref_seq = reference.get_seq(cluster.chrom_id, 0, chrom_len);

        // Align all gaps in the cluster
        cluster.align_gaps(strand_seq, ref_seq);

        // Split at any gaps where alignment failed
        let split_off = cluster.split_at_failed_alignments(min_seed_length);
        if !split_off.is_empty() {
            log::debug!(
                "Read {}: split cluster into {} additional clusters due to failed gap alignments",
                read_name,
                split_off.len(),
            );
            new_clusters_from_splits.extend(split_off);
        }
    }

    // Add the split-off clusters and re-sort
    if !new_clusters_from_splits.is_empty() {
        all_clusters.extend(new_clusters_from_splits);
        all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));
        log::info!(
            "Read {}: after alignment-based splitting, have {} clusters",
            read_name,
            all_clusters.len(),
        );
    }

    // Write debug chains TSV output (if debug file is initialized)
    if debug::is_enabled(DebugFile::Chains) {
        let cfg = config::get();
        for (i, cluster) in all_clusters.iter().enumerate() {
            let (read_start, read_end) = cluster.fwd_read_range(seq_len);
            let strand = if cluster.is_reverse { "-" } else { "+" };
            let chrom_name = reference.chrom_name(cluster.chrom_id);
            // Get reference range from the chain
            let ref_start = cluster.chain.first().map(|h| h.ref_pos).unwrap_or(0);
            let ref_end = cluster.chain.last().map(|h| h.ref_end()).unwrap_or(0);
            let chain_score = cluster.chain_score(
                cfg.seeding.gap_penalty_linear,
                cfg.seeding.gap_penalty_log,
            );
            debug::write(DebugFile::Chains, &format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.2}",
                read_name,
                i,
                read_start,
                read_end,
                seq_len,
                chrom_name,
                ref_start,
                ref_end,
                strand,
                cluster.chain.len(),
                cluster.total_seed_length(),
                cluster.read_coverage(seq_len),
                cluster.seed_density(),
                chain_score,
            ));
        }
    }

    if log::log_enabled!(log::Level::Debug) {
        let mut all_gaps: Vec<(usize, usize, usize, usize, usize)> = Vec::new();
        for (i, cluster) in all_clusters.iter().enumerate() {
            let (read_start, read_end) = cluster.fwd_read_range(seq_len);
            log::debug!(
                "  Cluster {}: {}-{} {} seeds on {} strand (chrom: {},seed length: {},coverage {:.2}%, density {:.2})",
                i + 1,
                read_start,
                read_end,
                cluster.chain.len(),
                if cluster.is_reverse { "REV" } else { "FWD" },
                reference.chrom_name(cluster.chrom_id),
                cluster.total_seed_length(),
                cluster.read_coverage(seq_len) * 100.0,
                cluster.seed_density(),
            );

            for ((gap_start, gap_end), gap_index) in cluster.gaps(seq_len, 2 * K) {
                all_gaps.push((gap_start, gap_end, i, gap_index, cluster.chrom_id));
            }
        }

        all_gaps.sort_unstable();

        for (gap_start, gap_end, cluster_index, gap_index, chrom_id) in all_gaps.iter() {
            log::debug!(
                "  Gap: read {}-{} (length {}) in cluster {}/{} (chrom {})",
                gap_start,
                gap_end,
                gap_end - gap_start,
                cluster_index,
                gap_index,
                reference.chrom_name(*chrom_id),
            );
        }
    }

    // =========================================================================
    // PASS 1.6: Split clusters at gaps filled by other clusters
    // =========================================================================
    // Identify gaps where another cluster provides coverage, indicating a
    // potential chimeric breakpoint. Split the cluster at such gaps rather
    // than bridging them with WFA.

    let gap_fills = analyze_gap_fills(
        &all_clusters,
        seq_len,
        cfg.seeding.min_gap_for_split,
        2*K,
        cfg.seeding.gap_fill_tolerance,
        cfg.seeding.min_gap_fill_coverage,
    );

    if !gap_fills.is_empty() {
        log::debug!(
            "Read {}: found {} gap fills for potential splitting",
            read_name,
            gap_fills.len(),
        );

        // Group splits by cluster and sort by gap index descending
        // so we can apply splits from back to front without invalidating indices
        let mut splits_by_cluster: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for fill in &gap_fills {
            splits_by_cluster
                .entry(fill.cluster_idx)
                .or_default()
                .push(fill.gap_seed_idx);
        }

        // Sort each cluster's splits in descending order
        for indices in splits_by_cluster.values_mut() {
            indices.sort_unstable_by(|a, b| b.cmp(a));
            indices.dedup(); // Remove duplicate split points
        }

        // Apply splits (in descending cluster index order to preserve indices)
        let mut cluster_indices: Vec<_> = splits_by_cluster.keys().copied().collect();
        cluster_indices.sort_unstable_by(|a, b| b.cmp(a));

        for cluster_idx in cluster_indices {
            let split_indices = &splits_by_cluster[&cluster_idx];
            for &gap_seed_idx in split_indices {
                if let Some(new_cluster) = all_clusters[cluster_idx].split_at_gap(gap_seed_idx) {
                    log::debug!(
                        "Read {}: split cluster {} at gap {}, new cluster has {} seeds, original now has {} seeds",
                        read_name,
                        cluster_idx,
                        gap_seed_idx,
                        new_cluster.chain.len(),
                        all_clusters[cluster_idx].chain.len(),
                    );
                    all_clusters.push(new_cluster);
                }
            }
        }

        // Re-sort after splitting
        all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));

        log::debug!(
            "Read {}: after gap-fill splitting, have {} clusters",
            read_name,
            all_clusters.len(),
        );
    }

    // =========================================================================
    // PASS 2: Build alignments from clusters
    // =========================================================================
    // Use the new build_alignment_from_cluster which uses pre-computed gap
    // alignments from PASS 1.5. Clusters split by gap-fills won't have gap
    // alignments, so those will be computed on demand.

    let mut candidates = Vec::with_capacity(all_clusters.len());
    for cluster in &all_clusters {
        let strand_seq = if cluster.is_reverse { &rc_seq } else { seq };
        if let Some(candidate) = build_alignment_from_cluster(
            read_name,
            cluster,
            strand_seq,
            seq_len,
            reference,
        ) {
            candidates.push(candidate);
        }
    }

    log::debug!(
        "Read {}: built {} candidate alignments",
        read_name,
        candidates.len(),
    );

    // Classify all candidate alignments
    let classified = classify_alignments(candidates, seq_len, reference.all_chrom_info());

    log::debug!(
        "Read {}: classified into {} alignments",
        read_name,
        classified.len()
    );

    // Check if we have any usable (non-LowQuality) alignments
    let has_usable_alignments = classified
        .iter()
        .any(|aln| aln.class != AlignmentClass::LowQuality);

    if !has_usable_alignments {
        // Output unmapped read (either no candidates or all filtered as low quality)
        let _ = writer.write_alignment(
            read_name,
            FLAG_UNMAPPED,
            "*",
            0,
            0,
            "*",
            "*",
            0,
            0,
            std::str::from_utf8(seq).unwrap_or("*"),
            "*",
            &[],
        );
    } else {
        // Output classified alignments
        for aln in &classified {
            let class_str = match aln.class {
                AlignmentClass::Primary => "primary",
                AlignmentClass::Secondary => "secondary",
                AlignmentClass::Supplementary => "supplementary",
                AlignmentClass::SecondarySupplementary => "secondary+supplementary",
                AlignmentClass::LowQuality => "lowqual",
            };
            let strand = if aln.candidate.is_reverse {
                "REV"
            } else {
                "FWD"
            };

            log::debug!(
                "Read {}: {} {} to {}:{}-{} (read {}..{}), mapq={}, raw_score={}, ctx_score={}, identity={:.1}%, homo_gaps={}, CIGAR={}",
                read_name,
                class_str,
                strand,
                reference.chrom_name(aln.candidate.chrom_id),
                aln.candidate.ref_start,
                aln.candidate.ref_end,
                aln.candidate.read_start,
                aln.candidate.read_end,
                aln.mapq,
                aln.candidate.alignment.score,
                aln.candidate.context_score.score,
                aln.candidate.context_score.identity * 100.0,
                aln.candidate.context_score.homopolymer_gap_bases,
                aln.candidate.alignment.cigar_string(),
            );

            // Output SAM record (skip low quality unless there's no primary)
            if aln.class != AlignmentClass::LowQuality {
                let chrom_name = reference.chrom_name(aln.candidate.chrom_id);
                let pos = aln.candidate.ref_start + 1; // SAM is 1-based
                // Use hard clips for secondary/supplementary to reduce output size
                let cigar = if aln.class == AlignmentClass::Primary {
                    aln.candidate.alignment.cigar_string()
                } else {
                    aln.candidate.alignment.cigar_string().replace('S', "H")
                };

                // For primary alignments (soft clips): full sequence
                // For secondary/supplementary (hard clips): just the aligned portion
                let (seq_str, qual_str, expected_query_len) = if aln.class
                    == AlignmentClass::Primary
                {
                    // Validate CIGAR: query length from CIGAR must match sequence length
                    let cigar_query_len = aln.candidate.alignment.query_length() as usize;
                    if cigar_query_len != seq_len {
                        log::error!(
                            "CIGAR/SEQ mismatch for {}: CIGAR query_len={}, seq_len={}, CIGAR={}",
                            read_name,
                            cigar_query_len,
                            seq_len,
                            cigar
                        );
                        // Skip this alignment to avoid producing invalid SAM
                        continue;
                    }

                    let seq_out = if aln.candidate.is_reverse {
                        String::from_utf8_lossy(&rc_seq).into_owned()
                    } else {
                        String::from_utf8_lossy(seq).into_owned()
                    };

                    let qual_out = qual
                        .and_then(|q| std::str::from_utf8(q).ok())
                        .map(|s| {
                            if aln.candidate.is_reverse {
                                s.chars().rev().collect::<String>()
                            } else {
                                s.to_string()
                            }
                        })
                        .unwrap_or_else(|| "*".to_string());

                    (seq_out, qual_out, seq_len)
                } else {
                    // Hard clips: output only the aligned portion
                    // The aligned portion is defined by the CIGAR, not by read_start/read_end
                    // We need to compute the actual start/end from the soft clips in the CIGAR
                    let cigar_ops = &aln.candidate.alignment.cigar;

                    // Count leading soft clip
                    let leading_clip = match cigar_ops.first() {
                        Some(CigarOp::SoftClip(n)) => *n as usize,
                        _ => 0,
                    };

                    // Count trailing soft clip
                    let trailing_clip = match cigar_ops.last() {
                        Some(CigarOp::SoftClip(n)) if cigar_ops.len() > 1 => *n as usize,
                        _ => 0,
                    };

                    // The aligned portion excludes soft clips
                    let aligned_start = leading_clip;
                    let aligned_end = seq_len - trailing_clip;
                    let aligned_len = aln.candidate.alignment.query_consumed() as usize;

                    // Verify consistency
                    if aligned_end - aligned_start != aligned_len {
                        log::warn!(
                            "Read {}: aligned region mismatch: clip calc gives {}-{}={}, query_consumed={}",
                            read_name,
                            aligned_start,
                            aligned_end,
                            aligned_end - aligned_start,
                            aligned_len
                        );
                    }

                    let seq_out = if aln.candidate.is_reverse {
                        // For reverse strand, use the pre-computed rc_seq
                        String::from_utf8_lossy(&rc_seq[aligned_start..aligned_end]).into_owned()
                    } else {
                        String::from_utf8_lossy(&seq[aligned_start..aligned_end]).into_owned()
                    };

                    let qual_out = qual
                        .and_then(|q| std::str::from_utf8(q).ok())
                        .map(|s| {
                            let chars: Vec<char> = s.chars().collect();
                            if aln.candidate.is_reverse {
                                // Reverse the aligned portion
                                chars[aligned_start..aligned_end]
                                    .iter()
                                    .rev()
                                    .collect::<String>()
                            } else {
                                chars[aligned_start..aligned_end].iter().collect::<String>()
                            }
                        })
                        .unwrap_or_else(|| "*".to_string());

                    (seq_out, qual_out, aligned_len)
                };

                // Validate that SEQ length matches expected
                if seq_str != "*" && seq_str.len() != expected_query_len {
                    log::error!(
                        "SEQ length mismatch for {}: seq_len={}, expected={}",
                        read_name,
                        seq_str.len(),
                        expected_query_len
                    );
                    continue;
                }

                // Build optional tags
                let nm = aln.candidate.edit_distance();
                let as_score = -aln.candidate.alignment.score; // Negate since lower is better internally

                // Build SA tag for supplementary alignments (points back to primary)
                let sa_tag = if aln.class == AlignmentClass::Supplementary {
                    // Find the primary alignment to reference in SA tag
                    if let Some(primary) = classified
                        .iter()
                        .find(|a| a.class == AlignmentClass::Primary)
                    {
                        let p_chrom = reference.chrom_name(primary.candidate.chrom_id);
                        let p_pos = primary.candidate.ref_start + 1;
                        let p_strand = if primary.candidate.is_reverse {
                            '-'
                        } else {
                            '+'
                        };
                        let p_cigar = primary.candidate.alignment.cigar_string();
                        let p_mapq = primary.mapq;
                        let p_nm = primary.candidate.edit_distance();
                        format!(
                            "\tSA:Z:{},{},{},{},{},{}",
                            p_chrom, p_pos, p_strand, p_cigar, p_mapq, p_nm
                        )
                    } else {
                        String::new()
                    }
                } else if aln.class == AlignmentClass::Primary {
                    // For primary, list all supplementary alignments
                    let supps: Vec<String> = classified
                        .iter()
                        .filter(|a| a.class == AlignmentClass::Supplementary)
                        .map(|s| {
                            let s_chrom = reference.chrom_name(s.candidate.chrom_id);
                            let s_pos = s.candidate.ref_start + 1;
                            let s_strand = if s.candidate.is_reverse { '-' } else { '+' };
                            let s_cigar = s.candidate.alignment.cigar_string();
                            let s_mapq = s.mapq;
                            let s_nm = s.candidate.edit_distance();
                            format!(
                                "{},{},{},{},{},{}",
                                s_chrom, s_pos, s_strand, s_cigar, s_mapq, s_nm
                            )
                        })
                        .collect();
                    if supps.is_empty() {
                        String::new()
                    } else {
                        format!("\tSA:Z:{}", supps.join(";"))
                    }
                } else {
                    String::new()
                };

                // Build tags vector
                let mut tags: Vec<(String, String)> = vec![
                    ("NM".to_string(), format!("i:{}", nm)),
                    ("AS".to_string(), format!("i:{}", as_score)),
                ];
                if !sa_tag.is_empty() {
                    // sa_tag starts with \tSA:Z:, extract the value
                    if let Some(value) = sa_tag.strip_prefix("\tSA:Z:") {
                        tags.push(("SA".to_string(), format!("Z:{}", value)));
                    }
                }

                // Final validation: CIGAR query length must match SEQ length
                let cigar_len = cigar_query_length(&cigar);
                let seq_len_actual = if seq_str == "*" { 0 } else { seq_str.len() };
                if cigar_len != seq_len_actual {
                    log::error!(
                        "CIGAR/SEQ length mismatch for read '{}': CIGAR implies {} bases, SEQ has {} bases. \
                         CIGAR={}, class={:?}, strand={}, chrom={}, pos={}",
                        read_name,
                        cigar_len,
                        seq_len_actual,
                        cigar,
                        aln.class,
                        if aln.candidate.is_reverse {
                            "REV"
                        } else {
                            "FWD"
                        },
                        chrom_name,
                        pos
                    );
                    // Skip writing this invalid alignment
                    continue;
                }

                let _ = writer.write_alignment(
                    read_name,
                    aln.sam_flag(),
                    chrom_name,
                    pos - 1, // write_alignment adds 1, so subtract here
                    aln.mapq,
                    &cigar,
                    "*",
                    0,
                    0,
                    &seq_str,
                    &qual_str,
                    &tags,
                );
            }
        }
    }
}

/// Write SAM header using the provided writer
pub fn write_sam_header<W: std::io::Write>(
    writer: &AlignmentWriter<W>,
    reference: &InMemoryReference,
    input_file: &str,
) -> std::io::Result<()> {
    // @SQ - Sequence dictionary (one per reference sequence)
    for (name, length) in reference.chromosomes() {
        writer.write_contig_header(name, length as usize)?;
    }

    // @PG - Program record
    writer.write_command_header(&format!("parallax align <reference> {}", input_file))?;

    Ok(())
}

/// Process reads from a FASTQ file (handles gzip, bzip2, xz compression transparently)
#[allow(dead_code)]
pub fn process_reads_from_fastq<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    fastq: &str,
) -> Result<()> {
    log::info!("Processing reads from {}", fastq);

    let stdout = std::io::stdout();
    let writer = AlignmentWriter::new(stdout.lock());

    write_sam_header(&writer, reference, fastq)?;

    let (decompressed_reader, format) = niffler::from_path(std::path::Path::new(fastq))
        .map_err(|e| ParallaxError::Other(Box::new(e)))?;
    if format != niffler::Format::No {
        log::info!("Detected {:?} compression", format);
    }
    let reader = std::io::BufReader::new(decompressed_reader);
    let mut reader = noodles::fastq::io::Reader::new(reader);

    for record in reader.records() {
        let record = record?;
        let read_name = std::str::from_utf8(record.name()).unwrap_or("?");
        let seq: &[u8] = record.sequence().as_ref();
        let qual: &[u8] = record.quality_scores().as_ref();

        align_read(index, reference, &writer, read_name, seq, Some(qual));
    }

    writer.flush()?;
    Ok(())
}

/// A read to be processed by a worker thread
struct ReadWork {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Process reads from a FASTQ file using multiple threads.
///
/// Reads are distributed to worker threads via a channel. The InMemoryReference
/// is shared across all threads via Arc (no per-thread cloning needed).
pub fn process_reads_parallel<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    fastq: &str,
    sam: Option<&str>,
    num_threads: usize,
) -> Result<()> {
    use crossbeam::channel::bounded;

    log::info!(
        "Processing reads from {} using {} threads",
        fastq,
        num_threads
    );

    let now = std::time::Instant::now();
    let mut num_records = 0;

    // Initialize debug files from config
    let cfg = config::get();
    debug::init(&cfg, reference)?;

    // Create writer - either to file or stdout
    let output: Box<dyn std::io::Write + Send> = match sam {
        Some(path) => {
            log::info!("Writing output to {}", path);
            Box::new(std::fs::File::create(path)?)
        }
        None => Box::new(std::io::stdout()),
    };
    let writer = Arc::new(AlignmentWriter::new(output));

    write_sam_header(&writer, reference, fastq)?;

    // Create a bounded channel for backpressure
    let (sender, receiver) = bounded::<ReadWork>(num_threads * 100);

    // Use crossbeam's scoped threads to safely borrow index, reference, and writer
    crossbeam::scope(|scope| {
        // Spawn worker threads
        for _ in 0..num_threads {
            let receiver = receiver.clone();
            let writer = writer.clone();
            scope.spawn(move |_| {
                while let Ok(work) = receiver.recv() {
                    align_read(
                        index,
                        reference,
                        &writer,
                        &work.name,
                        &work.seq,
                        Some(&work.qual),
                    );
                }
            });
        }

        // Read FASTQ and send to workers (in main thread)
        // Use niffler for transparent decompression (gzip, bzip2, xz)
        let (decompressed_reader, format) =
            niffler::from_path(std::path::Path::new(fastq)).expect("Failed to open FASTQ file");
        if format != niffler::Format::No {
            log::info!("Detected {:?} compression", format);
        }
        let reader = std::io::BufReader::new(decompressed_reader);
        let mut reader = noodles::fastq::io::Reader::new(reader);

        for record in reader.records() {
            let record = record.expect("Failed to read FASTQ record");
            let seq: &[u8] = record.sequence().as_ref();
            let qual: &[u8] = record.quality_scores().as_ref();
            let work = ReadWork {
                name: String::from_utf8_lossy(record.name()).into_owned(),
                seq: seq.to_vec(),
                qual: qual.to_vec(),
            };
            sender.send(work).expect("Failed to send work to thread");
            num_records += 1;
        }

        // Signal completion by dropping sender
        drop(sender);

        // Scoped threads automatically join when scope ends
    })
    .expect("Scoped thread panicked");

    writer.flush()?;

    // Flush all debug files
    debug::flush_all();

    let elapsed = now.elapsed();
    log::info!(
        "Completed processing reads {} from {} in {:.2?}",
        num_records,
        fastq,
        elapsed
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create a SeedHit with a dummy kmer value
    fn make_hit(chrom_id: usize, ref_pos: usize, read_pos: usize, match_len: usize) -> SeedHit {
        SeedHit::new(chrom_id, ref_pos, read_pos, 0, 1, match_len)
    }

    #[test]
    fn test_seed_hit_new() {
        let hit = SeedHit::new(1, 100, 50, 12345, 1, 20);
        assert_eq!(hit.chrom_id, 1);
        assert_eq!(hit.ref_pos, 100);
        assert_eq!(hit.read_pos, 50);
        assert_eq!(hit.kmer, 12345);
        assert_eq!(hit.kmer_uniqueness, 1);
        assert_eq!(hit.match_len, 20);
        // diagonal = ref_pos - read_pos = 100 - 50 = 50
        assert_eq!(hit.diagonal, 50);
    }

    #[test]
    fn test_seed_hit_read_end() {
        let hit = make_hit(0, 100, 50, 20);
        assert_eq!(hit.read_end(), 70); // 50 + 20
    }

    #[test]
    fn test_seed_hit_ref_end() {
        let hit = make_hit(0, 100, 50, 20);
        assert_eq!(hit.ref_end(), 120); // 100 + 20
    }

    #[test]
    fn test_seed_hit_diagonal_calculation() {
        // Forward diagonal (ref ahead of read)
        let hit1 = make_hit(0, 1000, 100, 20);
        assert_eq!(hit1.diagonal, 900);

        // Negative diagonal (read ahead of ref)
        let hit2 = make_hit(0, 100, 1000, 20);
        assert_eq!(hit2.diagonal, -900);

        // Zero diagonal (same position)
        let hit3 = make_hit(0, 500, 500, 20);
        assert_eq!(hit3.diagonal, 0);
    }

    #[test]
    fn test_extend_overlapping_same_diagonal() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);

        // Second seed at read pos 10, ref pos 110, length 20
        // This is on the same diagonal (110-10 = 100-0 = 100)
        // And overlaps: ref 110 < ref_end 120, gap is 10 < k=20
        let k = 20;
        let result = hit.extend(0, 110, 10, 0, 1, k);

        assert!(
            result.is_none(),
            "Should extend in place, not return new hit"
        );
        // New end should be: (110 - 100) + 20 = 30
        assert_eq!(hit.match_len, 30);
        assert_eq!(hit.ref_end(), 130);
        assert_eq!(hit.read_end(), 30);
    }

    #[test]
    fn test_extend_adjacent_same_diagonal() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);

        // Second seed starts exactly where first ends
        // read pos 20, ref pos 120, still same diagonal
        let k = 20;
        let result = hit.extend(0, 120, 20, 0, 1, k);

        // Gap is exactly 20, which equals k, so should NOT extend
        // because condition is: chrom_pos - self.ref_pos < self.match_len + k
        // 120 - 100 = 20 < 20 + 20 = 40, so it SHOULD extend
        assert!(result.is_none(), "Should extend in place");
        assert_eq!(hit.match_len, 40);
    }

    #[test]
    fn test_extend_gap_too_large() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Second seed with large gap (beyond match_len + k)
        // ref pos 200, read pos 100 (same diagonal = 100)
        // Gap check: 200 - 100 = 100 >= 20 + 20 = 40
        let k = 20;
        let result = hit.extend(0, 200, 100, 999, 1, k);

        assert!(result.is_some(), "Should return new hit due to large gap");
        assert_eq!(hit.match_len, original_len, "Original should be unchanged");

        let new_hit = result.unwrap();
        assert_eq!(new_hit.ref_pos, 200);
        assert_eq!(new_hit.read_pos, 100);
        assert_eq!(new_hit.kmer, 999);
    }

    #[test]
    fn test_extend_different_chromosome() {
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Same positions but different chromosome
        let k = 20;
        let result = hit.extend(1, 110, 10, 0, 1, k);

        assert!(
            result.is_some(),
            "Different chromosome should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
        assert_eq!(result.unwrap().chrom_id, 1);
    }

    #[test]
    fn test_extend_different_diagonal() {
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Different diagonal: ref_pos - read_pos = 111 - 10 = 101 != 100
        let k = 20;
        let result = hit.extend(0, 111, 10, 0, 1, k);

        assert!(result.is_some(), "Different diagonal should create new hit");
        assert_eq!(hit.match_len, original_len);

        let new_hit = result.unwrap();
        assert_eq!(new_hit.diagonal, 101);
    }

    #[test]
    fn test_extend_backwards_ref_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;

        // New ref_pos before current ref_pos
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, 1, k);

        assert!(
            result.is_some(),
            "Backwards ref position should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_backwards_read_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;

        // New read_pos before current read_pos (even if same diagonal)
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, 1, k);

        assert!(
            result.is_some(),
            "Backwards read position should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_fully_contained() {
        // Seed covering positions 0-20 in read, 100-120 in ref
        let mut hit = make_hit(0, 100, 0, 20);

        // New seed at pos 5-25 overlaps significantly
        // ref 105, read 5, same diagonal (100)
        let k = 20;
        let result = hit.extend(0, 105, 5, 0, 1, k);

        assert!(result.is_none(), "Overlapping hit should extend");
        // New end: (105 - 100) + 20 = 25
        assert_eq!(hit.match_len, 25);
    }

    #[test]
    fn test_extend_no_length_change_if_contained() {
        // Seed covering positions 0-30 in read
        let mut hit = make_hit(0, 100, 0, 30);

        // New seed fully contained within existing match
        // ref 110, read 10, k=20 means it ends at read 30, ref 130
        // That's exactly where the original ends, so no extension needed
        let k = 20;
        let result = hit.extend(0, 110, 10, 0, 1, k);

        assert!(result.is_none(), "Contained hit should not create new hit");
        // (110 - 100) + 20 = 30, which equals original, so no change
        assert_eq!(hit.match_len, 30);
    }

    #[test]
    fn test_extend_sequence_of_hits() {
        let k = 20;
        let mut hit = make_hit(0, 100, 0, k);

        // Simulate a sequence of overlapping syncmers ~6 bases apart
        // All on the same diagonal
        for i in 1..10 {
            let read_pos = i * 6;
            let ref_pos = 100 + i * 6;
            let result = hit.extend(0, ref_pos, read_pos, 0, 1, k);
            assert!(result.is_none(), "Hit {} should extend in place", i);
        }

        // Final length should cover from 0 to (9*6 + 20) = 74
        assert_eq!(hit.match_len, 9 * 6 + k);
        assert_eq!(hit.read_end(), 74);
        assert_eq!(hit.ref_end(), 174);
    }

    // =========================================================================
    // SeedCluster tests
    // =========================================================================

    #[test]
    fn test_seed_cluster_new_sorts_by_read_pos() {
        // Create seeds in reverse read_pos order
        let seeds = vec![
            make_hit(0, 300, 200, 20), // read_pos = 200
            make_hit(0, 100, 0, 20),   // read_pos = 0
            make_hit(0, 200, 100, 20), // read_pos = 100
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Should be sorted by read_pos
        assert_eq!(cluster.chain[0].read_pos, 0);
        assert_eq!(cluster.chain[1].read_pos, 100);
        assert_eq!(cluster.chain[2].read_pos, 200);

        assert_eq!(cluster.read_start, 0);
        assert_eq!(cluster.read_end, 220); // 200 + 20
    }

    #[test]
    fn test_seed_cluster_empty_returns_none() {
        let seeds: Vec<SeedHit> = vec![];
        assert!(SeedCluster::new(seeds, false, 1).is_none());
    }

    #[test]
    fn test_seed_cluster_fwd_read_range_forward_strand() {
        let seeds = vec![make_hit(0, 100, 50, 20)];
        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        let (start, end) = cluster.fwd_read_range(1000);
        assert_eq!(start, 50);
        assert_eq!(end, 70);
    }

    #[test]
    fn test_seed_cluster_fwd_read_range_reverse_strand() {
        // For reverse strand, coordinates need to be flipped
        let seeds = vec![make_hit(0, 100, 50, 20)];
        let cluster = SeedCluster::new(seeds, true, 1).unwrap();

        // read_start=50, read_end=70, read_len=1000
        // fwd_start = 1000 - 70 = 930, fwd_end = 1000 - 50 = 950
        let (start, end) = cluster.fwd_read_range(1000);
        assert_eq!(start, 930);
        assert_eq!(end, 950);
    }

    #[test]
    fn test_seed_cluster_split_at_gap() {
        // Create a chain with a gap between seeds 1 and 2
        let seeds = vec![
            make_hit(0, 100, 0, 20),   // 0-20
            make_hit(0, 200, 50, 20),  // 50-70
            make_hit(0, 400, 200, 20), // 200-220 (gap here)
            make_hit(0, 500, 250, 20), // 250-270
        ];

        let mut cluster = SeedCluster::new(seeds, false, 1).unwrap();
        assert_eq!(cluster.chain.len(), 4);

        // Split at gap between index 1 and 2
        let tail = cluster.split_at_gap(1).unwrap();

        // Original cluster should have seeds 0, 1
        assert_eq!(cluster.chain.len(), 2);
        assert_eq!(cluster.read_start, 0);
        assert_eq!(cluster.read_end, 70);

        // Tail cluster should have seeds 2, 3
        assert_eq!(tail.chain.len(), 2);
        assert_eq!(tail.read_start, 200);
        assert_eq!(tail.read_end, 270);
        assert_eq!(tail.is_reverse, cluster.is_reverse);
    }

    #[test]
    fn test_seed_cluster_split_preserves_strand() {
        let seeds = vec![make_hit(0, 100, 0, 20), make_hit(0, 300, 100, 20)];

        let mut cluster = SeedCluster::new(seeds, true, 1).unwrap();
        let tail = cluster.split_at_gap(0).unwrap();

        assert!(cluster.is_reverse);
        assert!(tail.is_reverse);
    }

    #[test]
    fn test_seed_cluster_gaps_forward_strand() {
        let seeds = vec![
            make_hit(0, 100, 0, 20),   // 0-20
            make_hit(0, 200, 50, 20),  // 50-70, gap of 30
            make_hit(0, 400, 200, 20), // 200-220, gap of 130
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();
        let gaps: Vec<_> = cluster.gaps(1000, 50).collect();

        // Only the gap of 130 should be returned (min_gap=50)
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].0, (70, 200)); // gap from 70 to 200
        assert_eq!(gaps[0].1, 1); // seed index 1 before this gap
    }

    #[test]
    fn test_seed_cluster_gaps_reverse_strand() {
        // For reverse strand, the chain is in RC coordinates
        // but gaps() should return forward-strand coordinates
        let seeds = vec![
            make_hit(0, 100, 0, 20),   // RC pos 0-20
            make_hit(0, 200, 50, 20),  // RC pos 50-70, gap of 30
            make_hit(0, 400, 200, 20), // RC pos 200-220, gap of 130
        ];

        let cluster = SeedCluster::new(seeds, true, 1).unwrap();
        let read_len = 1000;
        let gaps: Vec<_> = cluster.gaps(read_len, 50).collect();

        // Gap in RC coords: 70-200
        // In forward coords: (1000-200, 1000-70) = (800, 930)
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].0, (800, 930));
    }

    // =========================================================================
    // Chain colinearity tests
    // =========================================================================

    #[test]
    fn test_colinear_chain_both_dimensions_increasing() {
        // A proper colinear chain should have both ref_pos and read_pos increasing
        let seeds = vec![
            make_hit(0, 100, 10, 20),
            make_hit(0, 200, 50, 20),
            make_hit(0, 300, 100, 20),
            make_hit(0, 400, 150, 20),
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Verify both dimensions are strictly increasing
        for i in 1..cluster.chain.len() {
            assert!(
                cluster.chain[i].read_pos > cluster.chain[i - 1].read_end() - 1,
                "read_pos should be increasing: {} vs {}",
                cluster.chain[i].read_pos,
                cluster.chain[i - 1].read_end()
            );
            assert!(
                cluster.chain[i].ref_pos >= cluster.chain[i - 1].ref_end(),
                "ref_pos should be increasing: {} vs {}",
                cluster.chain[i].ref_pos,
                cluster.chain[i - 1].ref_end()
            );
        }
    }

    #[test]
    fn test_chain_ref_pos_monotonic_after_read_sort() {
        // This tests the invariant that should hold after SeedCluster::new
        // Even if seeds come in arbitrary order, after sorting by read_pos,
        // ref_pos should also be increasing for a proper colinear chain
        let seeds = vec![
            make_hit(0, 400, 150, 20), // Will be last after sort
            make_hit(0, 100, 10, 20),  // Will be first after sort
            make_hit(0, 300, 100, 20), // Will be third after sort
            make_hit(0, 200, 50, 20),  // Will be second after sort
        ];

        let cluster = SeedCluster::new(seeds, false, 1).unwrap();

        // Verify ref_pos is monotonically increasing
        for i in 1..cluster.chain.len() {
            assert!(
                cluster.chain[i].ref_pos >= cluster.chain[i - 1].ref_pos,
                "ref_pos not monotonic at {}: {} < {}",
                i,
                cluster.chain[i].ref_pos,
                cluster.chain[i - 1].ref_pos
            );
        }
    }

    // =========================================================================
    // LIS chain building tests
    // =========================================================================

    #[test]
    fn test_lis_on_ref_pos_produces_colinear_chain() {
        use crate::utils::LongestSubsequence;

        // Simulate DBSCAN cluster sorted by (diagonal, ref_pos)
        // All on same diagonal (100), so just sorted by ref_pos
        let hits = vec![
            make_hit(0, 100, 0, 20),   // diagonal = 100
            make_hit(0, 200, 100, 20), // diagonal = 100
            make_hit(0, 300, 200, 20), // diagonal = 100
        ];

        let mut lis = LongestSubsequence::default();
        let mut indices = Vec::new();

        // Old approach: LIS on ref_pos
        lis.longest_colinear_chain(&hits, |h| h.ref_pos as i64, true, &mut indices);

        let chain: Vec<_> = indices.iter().map(|&i| hits[i]).collect();

        // Should select all seeds since they're already colinear
        assert_eq!(chain.len(), 3);

        // Verify colinearity after sorting by read_pos
        let mut sorted = chain.clone();
        sorted.sort_by_key(|h| h.read_pos);

        for i in 1..sorted.len() {
            assert!(
                sorted[i].ref_pos >= sorted[i - 1].ref_pos,
                "ref_pos not monotonic after read_pos sort"
            );
        }
    }

    #[test]
    fn test_lis_on_read_pos_produces_colinear_chain() {
        use crate::utils::LongestSubsequence;

        // Simulate cluster sorted by ref_pos
        let hits = vec![
            make_hit(0, 100, 0, 20),
            make_hit(0, 200, 100, 20),
            make_hit(0, 300, 200, 20),
        ];

        let mut lis = LongestSubsequence::default();
        let mut indices = Vec::new();

        // New approach: sort by ref_pos, LIS on read_pos
        lis.longest_colinear_chain(&hits, |h| h.read_pos as i64, true, &mut indices);

        let chain: Vec<_> = indices.iter().map(|&i| hits[i]).collect();

        assert_eq!(chain.len(), 3);

        // Verify colinearity after sorting by read_pos
        let mut sorted = chain.clone();
        sorted.sort_by_key(|h| h.read_pos);

        for i in 1..sorted.len() {
            assert!(
                sorted[i].ref_pos >= sorted[i - 1].ref_pos,
                "ref_pos not monotonic after read_pos sort"
            );
        }
    }

    #[test]
    fn test_lis_filters_non_colinear_seeds() {
        use crate::utils::LongestSubsequence;

        // Seeds sorted by ref_pos, but one has backwards read_pos (non-colinear)
        let hits = vec![
            make_hit(0, 100, 50, 20),  // ref 100, read 50
            make_hit(0, 200, 30, 20),  // ref 200, read 30 - BACKWARDS!
            make_hit(0, 300, 150, 20), // ref 300, read 150
        ];

        let mut lis = LongestSubsequence::default();
        let mut indices = Vec::new();

        // LIS on read_pos should find an increasing sequence
        lis.longest_colinear_chain(&hits, |h| h.read_pos as i64, true, &mut indices);

        let chain: Vec<_> = indices.iter().map(|&i| hits[i]).collect();

        // Should select 2 seeds (either [0,2] or [1,2] - both are valid LIS of length 2)
        assert_eq!(chain.len(), 2);

        // Key invariant: after LIS, read_pos should be strictly increasing
        for i in 1..chain.len() {
            assert!(
                chain[i].read_pos > chain[i - 1].read_pos,
                "read_pos not strictly increasing: {} <= {}",
                chain[i].read_pos,
                chain[i - 1].read_pos
            );
        }

        // And since input was sorted by ref_pos, ref_pos should also be increasing
        for i in 1..chain.len() {
            assert!(
                chain[i].ref_pos > chain[i - 1].ref_pos,
                "ref_pos not strictly increasing: {} <= {}",
                chain[i].ref_pos,
                chain[i - 1].ref_pos
            );
        }
    }

    #[test]
    fn test_lis_on_ref_pos_may_break_colinearity() {
        use crate::utils::LongestSubsequence;

        // Seeds sorted by diagonal - this is what DBSCAN produces
        // Different diagonals can have different ref_pos orderings relative to read_pos
        let hits = vec![
            make_hit(0, 100, 50, 20),  // diagonal = 50
            make_hit(0, 150, 120, 20), // diagonal = 30 - different!
            make_hit(0, 200, 80, 20),  // diagonal = 120 - different!
        ];

        let mut lis = LongestSubsequence::default();
        let mut indices = Vec::new();

        // Old approach: LIS on ref_pos
        lis.longest_colinear_chain(&hits, |h| h.ref_pos as i64, true, &mut indices);

        let chain: Vec<_> = indices.iter().map(|&i| hits[i]).collect();

        // Sort by read_pos (what SeedCluster::new does)
        let mut sorted = chain.clone();
        sorted.sort_by_key(|h| h.read_pos);

        // Check if ref_pos is still monotonic - it might not be!
        let mut is_monotonic = true;
        for i in 1..sorted.len() {
            if sorted[i].ref_pos < sorted[i - 1].ref_pos {
                is_monotonic = false;
                break;
            }
        }

        // This test documents that the old approach CAN break colinearity
        // when seeds have different diagonals
        // (The actual assertion depends on what we want to enforce)
        println!("Old approach monotonic: {}", is_monotonic);
        println!(
            "Chain after read_pos sort: {:?}",
            sorted
                .iter()
                .map(|h| (h.ref_pos, h.read_pos))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_compare_lis_approaches_same_diagonal() {
        use crate::utils::LongestSubsequence;

        // Seeds all on the SAME diagonal (what DBSCAN should produce)
        // diagonal = 100 for all
        let hits_by_diag = vec![
            make_hit(0, 100, 0, 20),   // diagonal = 100
            make_hit(0, 200, 100, 20), // diagonal = 100
            make_hit(0, 300, 200, 20), // diagonal = 100
            make_hit(0, 400, 300, 20), // diagonal = 100
        ];

        let mut lis = LongestSubsequence::default();
        let mut indices_old = Vec::new();
        let mut indices_new = Vec::new();

        // Old approach: LIS on ref_pos (input sorted by diagonal, ref_pos)
        lis.longest_colinear_chain(&hits_by_diag, |h| h.ref_pos as i64, true, &mut indices_old);

        // New approach: sort by ref_pos, LIS on read_pos
        let mut hits_by_ref = hits_by_diag.clone();
        hits_by_ref.sort_by_key(|h| h.ref_pos);
        lis.longest_colinear_chain(&hits_by_ref, |h| h.read_pos as i64, true, &mut indices_new);

        println!("Same diagonal - Old approach indices: {:?}", indices_old);
        println!("Same diagonal - New approach indices: {:?}", indices_new);

        // Both should select all seeds
        assert_eq!(indices_old.len(), 4);
        assert_eq!(indices_new.len(), 4);
    }

    #[test]
    fn test_compare_lis_approaches_varying_diagonals() {
        use crate::utils::LongestSubsequence;

        // Seeds with VARYING diagonals (within DBSCAN tolerance)
        // This simulates real data where diagonals aren't exactly equal
        let hits_by_diag = vec![
            make_hit(0, 100, 0, 20),   // diagonal = 100
            make_hit(0, 205, 100, 20), // diagonal = 105 (slight variation)
            make_hit(0, 295, 200, 20), // diagonal = 95 (slight variation)
            make_hit(0, 400, 300, 20), // diagonal = 100
        ];

        // Sort by (diagonal, ref_pos) as DBSCAN would produce
        let mut hits_dbscan_order = hits_by_diag.clone();
        hits_dbscan_order.sort_by_key(|h| (h.diagonal, h.ref_pos));

        println!(
            "DBSCAN order: {:?}",
            hits_dbscan_order
                .iter()
                .map(|h| (h.ref_pos, h.read_pos, h.diagonal))
                .collect::<Vec<_>>()
        );

        let mut lis = LongestSubsequence::default();
        let mut indices_old = Vec::new();
        let mut indices_new = Vec::new();

        // Old approach: LIS on ref_pos with DBSCAN order
        lis.longest_colinear_chain(
            &hits_dbscan_order,
            |h| h.ref_pos as i64,
            true,
            &mut indices_old,
        );
        let chain_old: Vec<_> = indices_old.iter().map(|&i| hits_dbscan_order[i]).collect();

        // New approach: sort by ref_pos, LIS on read_pos
        let mut hits_by_ref = hits_by_diag.clone();
        hits_by_ref.sort_by_key(|h| h.ref_pos);
        lis.longest_colinear_chain(&hits_by_ref, |h| h.read_pos as i64, true, &mut indices_new);
        let chain_new: Vec<_> = indices_new.iter().map(|&i| hits_by_ref[i]).collect();

        println!(
            "Old approach chain: {:?}",
            chain_old
                .iter()
                .map(|h| (h.ref_pos, h.read_pos))
                .collect::<Vec<_>>()
        );
        println!(
            "New approach chain: {:?}",
            chain_new
                .iter()
                .map(|h| (h.ref_pos, h.read_pos))
                .collect::<Vec<_>>()
        );

        // Check colinearity after read_pos sort for old approach
        let mut sorted_old = chain_old.clone();
        sorted_old.sort_by_key(|h| h.read_pos);

        println!(
            "Old approach after read_pos sort: {:?}",
            sorted_old
                .iter()
                .map(|h| (h.ref_pos, h.read_pos))
                .collect::<Vec<_>>()
        );

        let mut old_colinear = true;
        for i in 1..sorted_old.len() {
            if sorted_old[i].ref_pos < sorted_old[i - 1].ref_pos {
                old_colinear = false;
                println!(
                    "Old approach breaks colinearity at index {}: ref {} < {}",
                    i,
                    sorted_old[i].ref_pos,
                    sorted_old[i - 1].ref_pos
                );
            }
        }

        // New approach should always be colinear (by construction)
        let mut new_colinear = true;
        for i in 1..chain_new.len() {
            if chain_new[i].ref_pos < chain_new[i - 1].ref_pos {
                new_colinear = false;
            }
        }

        println!("Old approach colinear: {}", old_colinear);
        println!("New approach colinear: {}", new_colinear);
    }

    // ==================== Tests for build_covering_sets ====================

    /// Helper to create a CandidateAlignment for testing
    fn make_candidate(
        chrom_id: usize,
        ref_start: usize,
        ref_end: usize,
        read_start: usize,
        read_end: usize,
        is_reverse: bool,
        matches: u32,
        mismatches: u32,
        gap_bases: u32,
    ) -> CandidateAlignment {
        use crate::align::{Alignment, ContextAwareScore};

        let aligned_len = matches + mismatches + gap_bases;
        let identity = if aligned_len > 0 {
            matches as f64 / aligned_len as f64
        } else {
            0.0
        };

        CandidateAlignment {
            chrom_id,
            ref_start,
            ref_end,
            read_start,
            read_end,
            is_reverse,
            alignment: Alignment {
                score: 0,
                cigar: vec![],
            },
            context_score: ContextAwareScore {
                score: 0,
                matches,
                mismatches,
                gap_bases,
                homopolymer_gap_bases: 0,
                identity,
            },
        }
    }

    #[test]
    fn test_build_covering_sets_single_alignment() {
        // Single alignment covering most of the read
        let candidates = vec![make_candidate(0, 1000, 2000, 0, 1000, false, 950, 25, 25)];
        let quality_indices = vec![0];
        let seq_len = 1000;

        let sets = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );

        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].alignment_indices, vec![0]);
        assert_eq!(sets[0].total_score, candidates[0].ranking_score() as f64);
    }

    #[test]
    fn test_build_covering_sets_non_overlapping_chimeric() {
        // Two non-overlapping alignments (chimeric read scenario)
        // These should be combinable into a single set
        let candidates = vec![
            make_candidate(0, 1000, 2000, 0, 1000, false, 950, 25, 25), // read [0, 1000]
            make_candidate(0, 3000, 4000, 1100, 2100, false, 950, 25, 25), // read [1100, 2100], no overlap
        ];
        let quality_indices = vec![0, 1];
        let seq_len = 2100;

        let sets = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );

        // Should create sets that combine both alignments
        assert!(!sets.is_empty());

        // The best set should contain both alignments since they don't overlap
        let best_set = &sets[0];
        assert!(best_set.alignment_indices.contains(&0));
        assert!(best_set.alignment_indices.contains(&1));

        // Combined score should be sum of individual scores
        let expected_score = (candidates[0].ranking_score() + candidates[1].ranking_score()) as f64;
        assert_eq!(best_set.total_score, expected_score);
    }

    #[test]
    fn test_build_covering_sets_overlapping_alternatives() {
        // Two overlapping alignments (alternative mappings)
        // These should NOT be in the same set
        let candidates = vec![
            make_candidate(0, 1000, 2000, 0, 1000, false, 950, 25, 25), // read [0, 1000]
            make_candidate(1, 5000, 6000, 100, 1100, false, 900, 50, 50), // read [100, 1100], overlaps significantly
        ];
        let quality_indices = vec![0, 1];
        let seq_len = 1100;

        let sets = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );

        // Should create separate sets for each since they overlap
        assert!(sets.len() >= 2);

        // No set should contain both alignments
        for set in &sets {
            let has_both = set.alignment_indices.contains(&0) && set.alignment_indices.contains(&1);
            assert!(
                !has_both,
                "Overlapping alignments should not be in the same set"
            );
        }
    }

    #[test]
    fn test_build_covering_sets_chimeric_vs_single_long() {
        // This test mirrors the real case we're debugging:
        // - Two high-identity alignments covering different parts of the read (chr16:2.1M case)
        // - One lower-identity alignment covering more of the read (chr16:18M case)
        //
        // Current behavior: The algorithm combines alignments purely by read position overlap
        // and score, without considering genomic location coherence or identity weighting.

        // Alignment A: chr16:2.1M first half, read [226, 4561], ~98.8% identity
        // matches=4300, mismatches=20, gaps=34 -> score = 4300*2 - 20*4 - 34*2 = 8452
        let aln_a = make_candidate(0, 2109126, 2113465, 226, 4561, true, 4300, 20, 34);

        // Alignment B: chr16:2.1M second half, read [4686, 8142], ~98.7% identity
        // matches=3427, mismatches=22, gaps=23 -> score = 3427*2 - 22*4 - 23*2 = 6720
        let aln_b = make_candidate(0, 2113632, 2117097, 4686, 8142, true, 3427, 22, 23);

        // Alignment C: chr16:18M, read [226, 6299], ~94.7% identity (longer but lower quality)
        // matches=5864, mismatches=153, gaps=175 -> score = 5864*2 - 153*4 - 175*2 = 10766
        let aln_c = make_candidate(0, 18385686, 18391822, 226, 6299, true, 5864, 153, 175);

        let candidates = vec![aln_a.clone(), aln_b.clone(), aln_c.clone()];
        let quality_indices = vec![0, 1, 2];

        println!(
            "Alignment A score: {} (id={:.1}%)",
            aln_a.ranking_score(),
            aln_a.identity() * 100.0
        );
        println!(
            "Alignment B score: {} (id={:.1}%)",
            aln_b.ranking_score(),
            aln_b.identity() * 100.0
        );
        println!(
            "Alignment C score: {} (id={:.1}%)",
            aln_c.ranking_score(),
            aln_c.identity() * 100.0
        );
        println!(
            "A + B combined: {}",
            aln_a.ranking_score() + aln_b.ranking_score()
        );

        let seq_len = 8142; // max read_end in test data
        let sets = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );

        println!("Built {} sets:", sets.len());
        for (i, set) in sets.iter().enumerate() {
            println!(
                "  Set {}: score={:.1} indices={:?}",
                i, set.total_score, set.alignment_indices
            );
        }

        // Check overlap between C and A (should NOT coexist)
        // C: [226, 6299], A: [226, 4561]
        // Overlap: 4335bp, C_len: 6073, A_len: 4335
        // C_ratio: 0.714 > 0.5, A_ratio: 1.0 > 0.5 -> Cannot coexist

        // Check overlap between C and B (CAN coexist!)
        // C: [226, 6299], B: [4686, 8142]
        // Overlap: 1613bp, C_len: 6073, B_len: 3456
        // C_ratio: 0.266, B_ratio: 0.467 -> Both <= 0.5, CAN coexist

        // So the algorithm correctly finds that C+B has higher score than A+B
        // This is mathematically correct but biologically suboptimal
        // because C is at a different genomic location with lower identity.

        // Verify A + B combined score > C alone (this still holds)
        let combined_score = aln_a.ranking_score() + aln_b.ranking_score();
        let single_score = aln_c.ranking_score();
        assert!(
            combined_score > single_score,
            "Combined A+B ({}) should beat C alone ({})",
            combined_score,
            single_score
        );

        // Current behavior: best set is C+B because they can coexist and have higher combined score
        let best_set = &sets[0];

        // The greedy algorithm starts from highest-scoring alignment (C) and adds compatible ones
        // C cannot coexist with A (too much overlap), but CAN coexist with B
        // So best set is {C, B} with score 10766 + 6720 = 17486

        // A+B set has score 8452 + 6720 = 15172, which is less than C+B
        assert!(
            best_set.alignment_indices.contains(&2),
            "Best set should contain C (alignment 2)"
        );
        assert!(
            best_set.alignment_indices.contains(&1),
            "Best set should contain B (alignment 1)"
        );

        // There should also be an A+B set
        let ab_set = sets
            .iter()
            .find(|s| s.alignment_indices.contains(&0) && s.alignment_indices.contains(&1));
        assert!(ab_set.is_some(), "Should have an A+B set");

        // TODO: Future enhancement - weight scores by identity to prefer high-identity alignments
        // and/or add genomic coherence scoring to prefer alignments at nearby genomic locations
    }

    #[test]
    fn test_build_covering_sets_respects_overlap_threshold() {
        // Two alignments with exactly 50% overlap
        // With threshold 0.5, they should be in separate sets
        // With threshold 0.6, they should be combinable

        // read [0, 1000] and read [500, 1500] have 500bp overlap
        // overlap/len = 500/1000 = 0.5 for both
        let candidates = vec![
            make_candidate(0, 1000, 2000, 0, 1000, false, 950, 25, 25),
            make_candidate(0, 3000, 4000, 500, 1500, false, 950, 25, 25),
        ];
        let quality_indices = vec![0, 1];
        let seq_len = 1500;

        // With threshold 0.5: overlap ratio = 0.5, which is NOT > 0.5, so they CAN coexist
        let sets_low_threshold = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );
        let best_set_low = &sets_low_threshold[0];
        let both_in_set_low = best_set_low.alignment_indices.contains(&0)
            && best_set_low.alignment_indices.contains(&1);
        assert!(
            both_in_set_low,
            "With threshold 0.5, 50% overlap should allow coexistence"
        );

        // With threshold 0.4: overlap ratio = 0.5 > 0.4, so they should NOT coexist
        let sets_high_threshold = build_covering_sets(
            &candidates,
            &quality_indices,
            0.4,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );
        for set in &sets_high_threshold {
            let both_in_set =
                set.alignment_indices.contains(&0) && set.alignment_indices.contains(&1);
            assert!(
                !both_in_set,
                "With threshold 0.4, 50% overlap should prevent coexistence"
            );
        }
    }

    #[test]
    fn test_build_covering_sets_ordering() {
        // Sets should be ordered by total score (highest first)
        let candidates = vec![
            make_candidate(0, 1000, 2000, 0, 1000, false, 900, 50, 50), // lower score
            make_candidate(0, 3000, 4000, 1100, 2100, false, 950, 25, 25), // higher score
        ];
        let quality_indices = vec![0, 1];
        let seq_len = 2100;

        let sets = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );

        // Sets should be in descending score order
        for i in 1..sets.len() {
            assert!(
                sets[i - 1].total_score >= sets[i].total_score,
                "Sets should be ordered by descending score"
            );
        }
    }

    #[test]
    fn test_build_covering_sets_deduplication() {
        // The greedy algorithm starting from different alignments might produce
        // the same set. These should be deduplicated.
        let candidates = vec![
            make_candidate(0, 1000, 2000, 0, 1000, false, 950, 25, 25),
            make_candidate(0, 3000, 4000, 1100, 2100, false, 900, 50, 50),
        ];
        let quality_indices = vec![0, 1];
        let seq_len = 2100;

        let sets = build_covering_sets(
            &candidates,
            &quality_indices,
            0.5,
            seq_len,
            0,
            0,
            false,
            4.0,
            6.0,
            1.0,
        );

        // Check for duplicate sets
        let mut seen_signatures: std::collections::HashSet<Vec<usize>> =
            std::collections::HashSet::new();
        for set in &sets {
            let mut signature = set.alignment_indices.clone();
            signature.sort();
            assert!(
                !seen_signatures.contains(&signature),
                "Found duplicate set: {:?}",
                signature
            );
            seen_signatures.insert(signature);
        }
    }
}
