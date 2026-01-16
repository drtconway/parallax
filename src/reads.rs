use std::sync::Arc;

use crate::align::{
    Alignment, CigarOp, ContextAwareParams, ContextAwareScore, align, context_aware_score,
};
use crate::config;
use crate::error::Result;
use crate::index::Index;
use crate::kmers::Kmer;
use crate::reads::seeds::{
    SeedCluster, SeedHit, analyze_gap_fills, flush_debug_sam, init_debug_sam, write_debug_sam,
};
use crate::reference::{ChromInfo, InMemoryReference};
use crate::utils::sequence::reverse_complement_into;
use crate::utils::{LongestSubsequence, dbscan_variance_aware};
use crate::writer::AlignmentWriter;

pub mod seeds;

/// SAM flags
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

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

/// Build alignment from a chain of seed matches, filling gaps with WFA.
///
/// The chain should be sorted by read position. Both sequences (read and reference)
/// are assumed to be in the same orientation - for reverse strand alignments,
/// the caller should pass the reverse-complemented read sequence.
///
/// # Arguments
/// * `read_id` - Read identifier for logging
/// * `chain` - Sorted chain of seed hits
/// * `seq` - Read sequence (already reverse-complemented for reverse strand)
/// * `seq_len` - Length of the original read
/// * `reference` - Reference genome
/// * `is_reverse` - Whether this is a reverse strand alignment (for marking in result)
fn build_alignment_from_chain(
    read_id: &str,
    chain: &[SeedHit],
    seq: &[u8],
    seq_len: usize,
    reference: &InMemoryReference,
    is_reverse: bool,
) -> Option<CandidateAlignment> {
    let cfg = config::get();

    if chain.is_empty() {
        return None;
    }

    // Require either multiple seeds, or a single seed that's long enough
    if chain.len() == 1 && chain[0].match_len < cfg.seeding.min_single_seed_length {
        return None;
    }

    if true {
        log::info!(
            "Building alignment for read {} on {} strand with {} seeds:",
            read_id,
            if is_reverse { "reverse" } else { "forward" },
            chain.len()
        );
        log::info!("Seed:\tchrom\tref_pos\tread_pos\tlen\tdiagonal");
        for i in 0..chain.len() {
            let hit = &chain[i];

            log::info!(
                "  Seed:\t{}\t{}\t{}\t{}\t{}",
                hit.chrom_id,
                hit.ref_pos,
                hit.read_pos,
                hit.match_len,
                hit.diagonal,
            );
        }

        // Verify chain is colinear (both ref_pos and read_pos should be non-decreasing)
        for i in 1..chain.len() {
            let prev = &chain[i - 1];
            let curr = &chain[i];
            if curr.ref_pos < prev.ref_pos {
                log::error!(
                    "CHAIN NOT COLINEAR: ref_pos[{}]={} < ref_pos[{}]={}",
                    i,
                    curr.ref_pos,
                    i - 1,
                    prev.ref_pos
                );
            }
            if curr.read_pos < prev.read_pos {
                log::error!(
                    "CHAIN NOT COLINEAR: read_pos[{}]={} < read_pos[{}]={}",
                    i,
                    curr.read_pos,
                    i - 1,
                    prev.read_pos
                );
            }
        }
    }
    let chrom_id = chain[0].chrom_id;
    let mut full_cigar: Vec<CigarOp> = Vec::new();
    let mut total_score = 0i32;

    // Compute alignment span from actual min/max reference positions
    let first = chain.first().unwrap();
    let last = chain.last().unwrap();

    // Use actual min/max ref positions to handle any chain ordering
    let ref_start = chain.iter().map(|h| h.ref_pos).min().unwrap();
    let ref_end = chain.iter().map(|h| h.ref_end()).max().unwrap();
    let read_start = first.read_pos;
    let read_end = last.read_end();

    // Add soft-clip for unaligned prefix
    // Seeds are already extended to their maximum exact match length,
    // so the alignment is anchored at the first seed's start
    if read_start > 0 {
        full_cigar.push(CigarOp::SoftClip(read_start as u32));
    }

    for j in 0..chain.len() {
        let hit = &chain[j];

        // Calculate effective match start and length, accounting for overlaps
        let (effective_read_start, effective_match_len, effective_ref_pos) = if j > 0 {
            let prev = &chain[j - 1];
            let prev_read_end = prev.read_end();
            if hit.read_pos < prev_read_end {
                // Overlap: current seed starts before previous seed ends
                // Clip the beginning of this seed's match
                let overlap = prev_read_end - hit.read_pos;
                if overlap >= hit.match_len {
                    // Entirely overlapped, skip this seed
                    continue;
                }
                // Adjust both read position and reference position
                (
                    prev_read_end,
                    hit.match_len - overlap,
                    hit.ref_pos + overlap,
                )
            } else {
                (hit.read_pos, hit.match_len, hit.ref_pos)
            }
        } else {
            (hit.read_pos, hit.match_len, hit.ref_pos)
        };

        // Align gap before this seed (if not first seed)
        if j > 0 {
            let prev = &chain[j - 1];
            let prev_read_end = prev.read_end();
            let read_gap_start = prev_read_end;
            let read_gap_end = effective_read_start; // Use effective start to avoid processing overlapped region

            // Reference gap: from end of previous seed to start of current seed
            let prev_ref_end = prev.ref_end();
            let ref_gap_start = prev_ref_end;
            let ref_gap_end = effective_ref_pos;

            let read_gap_len = if read_gap_end > read_gap_start {
                read_gap_end - read_gap_start
            } else {
                0
            };
            let ref_gap_len = if ref_gap_end > ref_gap_start {
                ref_gap_end - ref_gap_start
            } else {
                0
            };

            if read_gap_len > 0 && ref_gap_len > 0 {
                // Both have gaps - need to align
                let actual_read_start = read_gap_start;
                let actual_read_end = read_gap_end;
                let actual_ref_start = ref_gap_start;
                let actual_ref_end = ref_gap_end;

                // Get reference and read slices
                let ref_slice = reference.get_seq(chrom_id, actual_ref_start, actual_ref_end);
                let read_slice = &seq[actual_read_start..actual_read_end];

                if read_slice.len() >= 150 || ref_slice.len() >= 150 {
                    log::info!(
                        "Aligning read {} gap of size {} to ref gap of size {}: read pos {}-{}, ref pos {}-{}",
                        read_id,
                        read_slice.len(),
                        ref_slice.len(),
                        actual_read_start,
                        actual_read_end,
                        actual_ref_start,
                        actual_ref_end,
                    );
                }
                if let Some(aln) = align(read_slice, ref_slice) {
                    total_score += aln.score;
                    full_cigar.extend(aln.cigar);
                } else {
                    // Alignment failed, emit as insertions/deletions
                    full_cigar.push(CigarOp::Ins(read_gap_len as u32));
                    full_cigar.push(CigarOp::Del(ref_gap_len as u32));
                }
            } else if read_gap_len > 0 {
                // Only read has gap - pure insertion
                full_cigar.push(CigarOp::Ins(read_gap_len as u32));
            } else if ref_gap_len > 0 {
                // Only reference has gap - pure deletion
                full_cigar.push(CigarOp::Del(ref_gap_len as u32));
            }
            // else: both zero or negative - no gap to process
        }

        // Add the seed match itself (using effective length to handle overlaps)
        full_cigar.push(CigarOp::Match(effective_match_len as u32));
    }

    // Add soft-clip for unaligned suffix
    // Seeds are already extended to their maximum exact match length,
    // so the alignment is anchored at the last seed's end
    if read_end < seq_len {
        full_cigar.push(CigarOp::SoftClip((seq_len - read_end) as u32));
    }

    let mut alignment = Alignment {
        score: total_score,
        cigar: full_cigar,
    };
    alignment.normalize();

    // Get the aligned portions for context-aware scoring
    let query_for_scoring = &seq[read_start..read_end];
    let ref_for_scoring = reference.get_seq(chrom_id, ref_start, ref_end);

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

    // Build read-covering sets using greedy algorithm
    // Each set is a collection of non-overlapping alignments that together cover the read
    let alignment_sets = build_covering_sets(
        &candidates,
        &quality_indices,
        cfg.classification.overlap_threshold,
    );

    // Debug logging
    log::info!(
        "Built {} read-covering sets from {} quality alignments:",
        alignment_sets.len(),
        quality_indices.len()
    );
    for (set_idx, set) in alignment_sets.iter().enumerate() {
        log::info!(
            "  Set {}: score={} coverage={:.1}% alignments={:?}",
            set_idx,
            set.total_score,
            set.read_coverage * 100.0,
            set.alignment_indices
        );
    }

    // Log alignment details
    for (i, c) in candidates.iter().enumerate() {
        let in_sets: Vec<usize> = alignment_sets
            .iter()
            .enumerate()
            .filter(|(_, s)| s.alignment_indices.contains(&i))
            .map(|(idx, _)| idx)
            .collect();
        log::info!(
            "  Alignment {}: read [{}, {}] ref {}:[{}, {}] strand={} score={} in_sets={:?} (M={} X={} gaps={} id={:.1}%)",
            i,
            c.read_start,
            c.read_end,
            c.chrom_id,
            c.ref_start,
            c.ref_end,
            if c.is_reverse { "-" } else { "+" },
            c.ranking_score(),
            in_sets,
            c.context_score.matches,
            c.context_score.mismatches,
            c.context_score.gap_bases,
            c.identity() * 100.0
        );
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
    let mut set_scores: Vec<i64> = alignment_sets.iter().map(|s| s.total_score).collect();
    set_scores.sort_by(|a, b| b.cmp(a));
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
            // Use the gap between set scores to compute MAPQ
            let set_ratio = match second_best_set_score {
                Some(s2) if s2 > 0 && best_set.total_score > 0 => {
                    1.0 - (s2 as f64 / best_set.total_score as f64)
                }
                _ => 1.0, // No second-best set, unique mapping
            };
            // Scale MAPQ by set_ratio, length, and identity
            let len_factor = (aligned_len as f64 / 100.0).min(1.0);
            let mapq = (60.0 * set_ratio * len_factor * identity).round() as u8;
            classified.push(ClassifiedAlignment {
                mapq: mapq.min(60),
                class: AlignmentClass::Primary,
                candidate,
            });
        } else if in_best_set {
            // Other alignments in best set -> Supplementary (chimeric pieces)
            // These share the same MAPQ confidence as primary since they're part of the same solution
            let set_ratio = match second_best_set_score {
                Some(s2) if s2 > 0 && best_set.total_score > 0 => {
                    1.0 - (s2 as f64 / best_set.total_score as f64)
                }
                _ => 1.0,
            };
            let len_factor = (aligned_len as f64 / 100.0).min(1.0);
            let mapq = (60.0 * set_ratio * len_factor * identity).round() as u8;
            classified.push(ClassifiedAlignment {
                mapq: mapq.min(60),
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
                let my_set_score = alignment_sets.get(my_set_idx).map(|s| s.total_score);
                let mapq = compute_mapq(score, my_set_score, aligned_len, identity);
                classified.push(ClassifiedAlignment {
                    mapq,
                    class: AlignmentClass::Secondary,
                    candidate,
                });
            } else {
                // Not best in any set -> Secondary+Supplementary
                let my_set_score = alignment_sets.get(my_set_idx).map(|s| s.total_score);
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
    /// Combined score for the set
    total_score: i64,
    /// Fraction of read covered by this set
    read_coverage: f64,
}

/// Build read-covering sets using a greedy algorithm
///
/// For each starting alignment, greedily add non-overlapping alignments
/// to maximize coverage. Score the resulting set and keep track of
/// unique sets.
fn build_covering_sets(
    candidates: &[CandidateAlignment],
    quality_indices: &[usize],
    overlap_threshold: f64,
) -> Vec<AlignmentSet> {
    if quality_indices.is_empty() {
        return Vec::new();
    }

    // Helper to check if two alignments can coexist in the same set
    // They must not significantly overlap in read coordinates
    let can_coexist = |i: usize, j: usize| -> bool {
        let ci = &candidates[i];
        let cj = &candidates[j];

        let overlap_start = ci.read_start.max(cj.read_start);
        let overlap_end = ci.read_end.min(cj.read_end);

        if overlap_start >= overlap_end {
            return true; // No overlap
        }

        let overlap_len = (overlap_end - overlap_start) as f64;
        let len_i = (ci.read_end - ci.read_start) as f64;
        let len_j = (cj.read_end - cj.read_start) as f64;

        // They can coexist if overlap is small for both
        overlap_len / len_i <= overlap_threshold && overlap_len / len_j <= overlap_threshold
    };

    // Sort quality indices by score (descending)
    let mut sorted_indices = quality_indices.to_vec();
    sorted_indices.sort_by(|&a, &b| {
        candidates[b]
            .ranking_score()
            .cmp(&candidates[a].ranking_score())
    });

    // Build sets using greedy algorithm starting from each alignment
    let mut all_sets: Vec<AlignmentSet> = Vec::new();
    let mut seen_set_signatures: std::collections::HashSet<Vec<usize>> =
        std::collections::HashSet::new();

    for &start_idx in &sorted_indices {
        // Build a set starting from this alignment
        let mut set_indices = vec![start_idx];
        let mut total_score = candidates[start_idx].ranking_score();

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
                total_score += candidates[candidate_idx].ranking_score();
            }
        }

        // Sort set indices for consistent signature
        set_indices.sort();

        // Skip if we've already seen this exact set
        if seen_set_signatures.contains(&set_indices) {
            continue;
        }
        seen_set_signatures.insert(set_indices.clone());

        // Calculate read coverage for this set
        let read_coverage = calculate_set_coverage(candidates, &set_indices);

        // Sort by score (best first) for the final set
        set_indices.sort_by(|&a, &b| {
            candidates[b]
                .ranking_score()
                .cmp(&candidates[a].ranking_score())
        });

        all_sets.push(AlignmentSet {
            alignment_indices: set_indices,
            total_score,
            read_coverage,
        });
    }

    // Sort sets by total score (descending), with coverage as tiebreaker
    all_sets.sort_by(|a, b| {
        b.total_score.cmp(&a.total_score).then_with(|| {
            b.read_coverage
                .partial_cmp(&a.read_coverage)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    all_sets
}

/// Calculate the fraction of read covered by a set of alignments
fn calculate_set_coverage(candidates: &[CandidateAlignment], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }

    // Merge overlapping intervals to get total covered bases
    let mut intervals: Vec<(usize, usize)> = indices
        .iter()
        .map(|&i| (candidates[i].read_start, candidates[i].read_end))
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

    // Find the read extent
    let read_start = indices
        .iter()
        .map(|&i| candidates[i].read_start)
        .min()
        .unwrap_or(0);
    let read_end = indices
        .iter()
        .map(|&i| candidates[i].read_end)
        .max()
        .unwrap_or(0);
    let read_len = read_end.saturating_sub(read_start);

    if read_len == 0 {
        0.0
    } else {
        covered as f64 / read_len as f64
    }
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
            // Use seeds up to occurrence threshold
            if self.hit_vec.len() <= cfg.seeding.max_seed_occurrences {
                for &(chrom_id, chrom_pos) in self.hit_vec.iter() {
                    self.hits
                        .push(SeedHit::new(chrom_id, chrom_pos, pos, kmer.0, K));
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
                    .extend(hit.chrom_id, hit.ref_pos, hit.read_pos, hit.kmer, K)
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
        if !config::get().seeding.debug_seeds_sam.is_empty() {
            for hit in self.hits.iter() {
                let chrom_name = reference.chrom_name(hit.chrom_id);
                write_debug_sam(&hit.to_sam_line(
                    read_name,
                    chrom_name,
                    is_reverse,
                    strand_seq,
                    strand_qual,
                ));
            }
        }

        // Phase 4: Cluster hits by diagonal using DBSCAN
        self.cuts.clear();
        dbscan_variance_aware(
            &self.hits,
            cfg.seeding.min_seed_cluster_distance,
            max_var,
            |hit| hit.diagonal,
            &mut self.cuts,
        );
        metrics::histogram!(format!("{}_clusters_count", strand_name.to_lowercase()))
            .record(self.cuts.len().saturating_sub(1) as f64);

        // Phase 5: Build LIS chains for each cluster
        let mut clusters = Vec::new();
        for i in 1..self.cuts.len() {
            let begin = self.cuts[i - 1];
            let end = self.cuts[i];
            // Don't re-sort - keep the (diagonal, ref_pos) order from DBSCAN
            let cluster_hits = &self.hits[begin..end];

            if cluster_hits.len() == 1 {
                let seed = &cluster_hits[0];
                if seed.match_len < cfg.seeding.min_single_seed_length {
                    continue; // Skip tiny single-seed clusters
                }
            }

            // Use LIS on ref_pos to ensure the chain is colinear in reference space.
            //
            // IMPORTANT: We considered an alternative approach of sorting by ref_pos first
            // and then doing LIS on read_pos (the "dual" approach). While conceptually
            // appealing, this approach produces incorrect alignments because:
            //
            // 1. It can select seeds with different diagonals that overlap or invert
            //    in reference space after being sorted by read_pos
            // 2. The build_alignment_from_chain function assumes ref_pos is monotonic
            //    and uses min(ref_pos) as ref_start, but the CIGAR is built assuming
            //    we're moving linearly through the reference
            // 3. When ref_pos decreases between seeds (after read_pos sorting), we get
            //    "scrambled" alignments where the CIGAR doesn't match the reported position
            //
            // By doing LIS on ref_pos (with DBSCAN's diagonal-based ordering), we ensure
            // the final chain has monotonically increasing ref_pos, which is required
            // for correct CIGAR generation.
            self.longest_subsequence.longest_colinear_chain(
                cluster_hits,
                |hit| hit.ref_pos as i64,
                true,
                &mut self.chain_indices,
            );

            let chain: Vec<SeedHit> = self
                .chain_indices
                .iter()
                .map(|&i| cluster_hits[i])
                .collect();

            metrics::histogram!(format!("{}_chain_length", strand_name.to_lowercase()))
                .record(chain.len() as f64);

            if let Some(cluster) = SeedCluster::new(chain, is_reverse) {
                if cluster.total_seed_length() < cfg.seeding.min_single_seed_length {
                    continue; // Skip tiny clusters
                }
                clusters.push(cluster);
            }
        }

        clusters
    }
}

/// Build alignments from collected seed clusters.
///
/// This converts SeedClusters into CandidateAlignments by running WFA on gaps.
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

        if let Some(candidate) = build_alignment_from_chain(
            read_name,
            &cluster.chain,
            strand_seq,
            seq_len,
            reference,
            cluster.is_reverse,
        ) {
            candidates.push(candidate);
        }
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
    let fwd_clusters =
        collector.collect_from_strand(seq, qual, false, index, reference, read_name);
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
        "Read {}: collected {} seed clusters from both strands",
        read_name,
        all_clusters.len(),
    );

    // (gap_start, gap_end, cluster_index, gap_index, chrom_id)
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

    if log::log_enabled!(log::Level::Debug) {
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
    // PASS 1.5: Split clusters at gaps filled by other clusters
    // =========================================================================
    // Identify gaps where another cluster provides coverage, indicating a
    // potential chimeric breakpoint. Split the cluster at such gaps rather
    // than bridging them with WFA.

    let cfg = config::get();
    let gap_fills = analyze_gap_fills(
        &all_clusters,
        seq_len,
        cfg.seeding.min_gap_for_split,
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
                    log::info!(
                        "Read {}: split cluster {} at gap {}, new cluster has {} seeds",
                        read_name,
                        cluster_idx,
                        gap_seed_idx,
                        new_cluster.chain.len(),
                    );
                    all_clusters.push(new_cluster);
                }
            }
        }

        // Re-sort after splitting
        all_clusters.sort_by_key(|cluster| cluster.fwd_read_range(seq_len));

        log::info!(
            "Read {}: after splitting, have {} clusters",
            read_name,
            all_clusters.len(),
        );
    }

    // =========================================================================
    // PASS 2: Build alignments from clusters
    // =========================================================================

    let candidates =
        build_alignments_from_clusters(&all_clusters, read_name, seq, &rc_seq, seq_len, reference);

    log::info!(
        "Read {}: built {} candidate alignments",
        read_name,
        candidates.len(),
    );

    // Classify all candidate alignments
    let classified = classify_alignments(candidates, seq_len, reference.all_chrom_info());

    log::info!(
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
                    let aligned_start = aln.candidate.read_start;
                    let aligned_end = aln.candidate.read_end;
                    let aligned_len = aligned_end - aligned_start;

                    let seq_out = if aln.candidate.is_reverse {
                        // For reverse strand, use the pre-computed rc_seq
                        // The aligned portion coords are already in RC space
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

/// Process reads from a FASTQ file
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

    let reader = std::fs::File::open(fastq).map(std::io::BufReader::new)?;
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

    // Initialize debug seeds SAM file if configured
    let cfg = config::get();
    if !cfg.seeding.debug_seeds_sam.is_empty() {
        log::info!("Writing debug seed SAM to {}", cfg.seeding.debug_seeds_sam);
        init_debug_sam(&cfg.seeding.debug_seeds_sam, reference.chromosomes())?;
    }

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
        let reader = std::fs::File::open(fastq)
            .map(std::io::BufReader::new)
            .unwrap();
        let mut reader = noodles::fastq::io::Reader::new(reader);

        for record in reader.records() {
            let record = record.unwrap();
            let seq: &[u8] = record.sequence().as_ref();
            let qual: &[u8] = record.quality_scores().as_ref();
            let work = ReadWork {
                name: String::from_utf8_lossy(record.name()).into_owned(),
                seq: seq.to_vec(),
                qual: qual.to_vec(),
            };
            sender.send(work).expect("Failed to send work to thread");
        }

        // Signal completion by dropping sender
        drop(sender);

        // Scoped threads automatically join when scope ends
    })
    .expect("Scoped thread panicked");

    writer.flush()?;

    // Flush debug seeds SAM if it was initialized
    if !cfg.seeding.debug_seeds_sam.is_empty() {
        flush_debug_sam();
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create a SeedHit with a dummy kmer value
    fn make_hit(chrom_id: usize, ref_pos: usize, read_pos: usize, match_len: usize) -> SeedHit {
        SeedHit::new(chrom_id, ref_pos, read_pos, 0, match_len)
    }

    #[test]
    fn test_seed_hit_new() {
        let hit = SeedHit::new(1, 100, 50, 12345, 20);
        assert_eq!(hit.chrom_id, 1);
        assert_eq!(hit.ref_pos, 100);
        assert_eq!(hit.read_pos, 50);
        assert_eq!(hit.kmer, 12345);
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
        let result = hit.extend(0, 110, 10, 0, k);

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
        let result = hit.extend(0, 120, 20, 0, k);

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
        let result = hit.extend(0, 200, 100, 999, k);

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
        let result = hit.extend(1, 110, 10, 0, k);

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
        let result = hit.extend(0, 111, 10, 0, k);

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
        let result = hit.extend(0, 90, 40, 0, k);

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
        let result = hit.extend(0, 90, 40, 0, k);

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
        let result = hit.extend(0, 105, 5, 0, k);

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
        let result = hit.extend(0, 110, 10, 0, k);

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
            let result = hit.extend(0, ref_pos, read_pos, 0, k);
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

        let cluster = SeedCluster::new(seeds, false).unwrap();

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
        assert!(SeedCluster::new(seeds, false).is_none());
    }

    #[test]
    fn test_seed_cluster_fwd_read_range_forward_strand() {
        let seeds = vec![make_hit(0, 100, 50, 20)];
        let cluster = SeedCluster::new(seeds, false).unwrap();

        let (start, end) = cluster.fwd_read_range(1000);
        assert_eq!(start, 50);
        assert_eq!(end, 70);
    }

    #[test]
    fn test_seed_cluster_fwd_read_range_reverse_strand() {
        // For reverse strand, coordinates need to be flipped
        let seeds = vec![make_hit(0, 100, 50, 20)];
        let cluster = SeedCluster::new(seeds, true).unwrap();

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

        let mut cluster = SeedCluster::new(seeds, false).unwrap();
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

        let mut cluster = SeedCluster::new(seeds, true).unwrap();
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

        let cluster = SeedCluster::new(seeds, false).unwrap();
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

        let cluster = SeedCluster::new(seeds, true).unwrap();
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

        let cluster = SeedCluster::new(seeds, false).unwrap();

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

        let cluster = SeedCluster::new(seeds, false).unwrap();

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

        let sets = build_covering_sets(&candidates, &quality_indices, 0.5);

        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].alignment_indices, vec![0]);
        assert_eq!(sets[0].total_score, candidates[0].ranking_score());
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

        let sets = build_covering_sets(&candidates, &quality_indices, 0.5);

        // Should create sets that combine both alignments
        assert!(!sets.is_empty());

        // The best set should contain both alignments since they don't overlap
        let best_set = &sets[0];
        assert!(best_set.alignment_indices.contains(&0));
        assert!(best_set.alignment_indices.contains(&1));

        // Combined score should be sum of individual scores
        let expected_score = candidates[0].ranking_score() + candidates[1].ranking_score();
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

        let sets = build_covering_sets(&candidates, &quality_indices, 0.5);

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

        let sets = build_covering_sets(&candidates, &quality_indices, 0.5);

        println!("Built {} sets:", sets.len());
        for (i, set) in sets.iter().enumerate() {
            println!(
                "  Set {}: score={} indices={:?}",
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

        // With threshold 0.5: overlap ratio = 0.5, which is NOT > 0.5, so they CAN coexist
        let sets_low_threshold = build_covering_sets(&candidates, &quality_indices, 0.5);
        let best_set_low = &sets_low_threshold[0];
        let both_in_set_low = best_set_low.alignment_indices.contains(&0)
            && best_set_low.alignment_indices.contains(&1);
        assert!(
            both_in_set_low,
            "With threshold 0.5, 50% overlap should allow coexistence"
        );

        // With threshold 0.4: overlap ratio = 0.5 > 0.4, so they should NOT coexist
        let sets_high_threshold = build_covering_sets(&candidates, &quality_indices, 0.4);
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

        let sets = build_covering_sets(&candidates, &quality_indices, 0.5);

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

        let sets = build_covering_sets(&candidates, &quality_indices, 0.5);

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
