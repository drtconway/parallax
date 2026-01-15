use std::sync::Arc;

use crate::align::{
    Alignment, CigarOp, ContextAwareParams, ContextAwareScore, align, context_aware_score,
};
use crate::error::Result;
use crate::index::Index;
use crate::kmers::Kmer;
use crate::reads::seeds::{SeedHit, flush_debug_sam, init_debug_sam, write_debug_sam};
use crate::reference::InMemoryReference;
use crate::utils::sequence::reverse_complement_into;
use crate::utils::{LongestSubsequence, dbscan_variance_aware};
use crate::writer::AlignmentWriter;

pub mod seeds;

/// SAM flags
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Minimum alignment identity (matches / aligned_length) for a valid alignment
const MIN_ALIGNMENT_IDENTITY: f64 = 0.5;

/// Maximum context-aware score per aligned base (higher = more errors)
const MAX_SCORE_PER_BASE: f64 = 0.3;

/// Minimum fraction of read covered for a valid alignment
const MIN_READ_COVERAGE: f64 = 0.1;

/// Minimum aligned length (bp) - alignments meeting this threshold bypass coverage check
/// This handles chimeric reads where a small portion aligns elsewhere
const MIN_ALIGNED_LENGTH: u32 = 50;

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
    // Minimum length for a single seed to be accepted as a valid chain
    const MIN_SINGLE_SEED_LENGTH: usize = 50;

    if chain.is_empty() {
        return None;
    }

    // Require either multiple seeds, or a single seed that's long enough
    if chain.len() == 1 && chain[0].match_len < MIN_SINGLE_SEED_LENGTH {
        return None;
    }

    if true {
        log::info!(
            "Building alignment for read {} on {} strand with {} seeds:",
            read_id,
            if is_reverse { "reverse" } else { "forward" },
            chain.len()
        );
        log::info!("Seed:\tchrom\tpos\tread\tlen");
        for i in 0..chain.len() {
            let hit = &chain[i];

            log::info!(
                "  Seed:\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}",
                hit.chrom_id,
                hit.ref_pos as i64,
                hit.read_pos as i64,
                hit.match_len,
                (hit.ref_pos as i64) as f64 / 20.0,
                (hit.read_pos as i64) as f64 / 20.0,
                hit.match_len as f64 / 20.0
            );
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

/// Cluster alignments by overlapping read regions using union-find.
/// Returns a vector of cluster indices, one per alignment.
fn cluster_by_read_region(candidates: &[CandidateAlignment], overlap_threshold: f64) -> Vec<usize> {
    let n = candidates.len();
    if n == 0 {
        return Vec::new();
    }

    // Union-find parent array
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]]; // Path compression
            i = parent[i];
        }
        i
    }

    fn union(parent: &mut [usize], i: usize, j: usize) {
        let pi = find(parent, i);
        let pj = find(parent, j);
        if pi != pj {
            parent[pi] = pj;
        }
    }

    // Check each pair for significant overlap
    for i in 0..n {
        for j in (i + 1)..n {
            let ci = &candidates[i];
            let cj = &candidates[j];

            let overlap_start = ci.read_start.max(cj.read_start);
            let overlap_end = ci.read_end.min(cj.read_end);

            if overlap_start < overlap_end {
                let overlap_len = (overlap_end - overlap_start) as f64;
                let len_i = (ci.read_end - ci.read_start) as f64;
                let len_j = (cj.read_end - cj.read_start) as f64;

                // BOTH alignments must have >threshold overlap to be in the same group
                // This prevents transitive chaining (A overlaps B, B overlaps C -> A,B,C together)
                if overlap_len / len_i > overlap_threshold
                    && overlap_len / len_j > overlap_threshold
                {
                    union(&mut parent, i, j);
                }
            }
        }
    }

    // Flatten parent array
    for i in 0..n {
        find(&mut parent, i);
    }

    // Renumber clusters to be contiguous
    let mut cluster_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut next_cluster = 0;
    let mut clusters = Vec::with_capacity(n);

    for i in 0..n {
        let root = find(&mut parent, i);
        let cluster = *cluster_map.entry(root).or_insert_with(|| {
            let c = next_cluster;
            next_cluster += 1;
            c
        });
        clusters.push(cluster);
    }

    clusters
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
fn classify_alignments(
    mut candidates: Vec<CandidateAlignment>,
    read_len: usize,
) -> Vec<ClassifiedAlignment> {
    if candidates.is_empty() {
        return Vec::new();
    }

    // Helper to check if an alignment passes quality thresholds
    let passes_quality = |c: &CandidateAlignment| -> bool {
        let coverage = c.read_coverage(read_len);
        let aligned_len = c.aligned_length();
        let passes_coverage = coverage >= MIN_READ_COVERAGE || aligned_len >= MIN_ALIGNED_LENGTH;
        passes_coverage
            && c.identity() >= MIN_ALIGNMENT_IDENTITY
            && c.score_per_base() <= MAX_SCORE_PER_BASE
    };

    // Sort by ranking_score (higher is better) - this naturally prefers longer, accurate alignments
    candidates.sort_by(|a, b| b.ranking_score().cmp(&a.ranking_score()));

    // Cluster alignments by overlapping read regions
    let clusters = cluster_by_read_region(&candidates, 0.5);

    // Find the best alignment in each cluster (first one after sorting by score)
    let mut cluster_best: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    for (i, &cluster) in clusters.iter().enumerate() {
        if passes_quality(&candidates[i]) {
            cluster_best.entry(cluster).or_insert(i);
        }
    }

    // The primary is the best alignment overall (index 0 if it passes quality)
    let primary_idx = candidates
        .iter()
        .enumerate()
        .find(|(_, c)| passes_quality(c))
        .map(|(i, _)| i);

    let primary_cluster = primary_idx.map(|i| clusters[i]);

    // Build map of cluster -> best score for MAPQ calculation
    let mut cluster_best_score: std::collections::HashMap<usize, i64> =
        std::collections::HashMap::new();
    for (i, &cluster) in clusters.iter().enumerate() {
        if passes_quality(&candidates[i]) {
            cluster_best_score
                .entry(cluster)
                .or_insert_with(|| candidates[i].ranking_score());
        }
    }

    // Collect second-best scores per cluster for MAPQ
    let mut cluster_scores: std::collections::HashMap<usize, Vec<i64>> =
        std::collections::HashMap::new();
    for (i, &cluster) in clusters.iter().enumerate() {
        if passes_quality(&candidates[i]) {
            cluster_scores
                .entry(cluster)
                .or_default()
                .push(candidates[i].ranking_score());
        }
    }

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
        let cluster = clusters[i];

        // Get second-best score in this cluster for MAPQ
        let scores = cluster_scores
            .get(&cluster)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let second_best = if scores.len() > 1 {
            Some(scores[1])
        } else {
            None
        };

        if Some(i) == primary_idx {
            // Primary alignment
            let mapq = compute_mapq(score, second_best, aligned_len, identity);
            classified.push(ClassifiedAlignment {
                mapq,
                class: AlignmentClass::Primary,
                candidate,
            });
        } else if Some(cluster) == primary_cluster {
            // Same cluster as primary -> Secondary
            let best_score = cluster_best_score.get(&cluster).copied();
            let mapq = compute_mapq(score, best_score, aligned_len, identity);
            classified.push(ClassifiedAlignment {
                mapq,
                class: AlignmentClass::Secondary,
                candidate,
            });
        } else if cluster_best.get(&cluster) == Some(&i) {
            // Best in a non-primary cluster -> Supplementary
            let mapq = compute_mapq(score, second_best, aligned_len, identity);
            classified.push(ClassifiedAlignment {
                mapq,
                class: AlignmentClass::Supplementary,
                candidate,
            });
        } else {
            // Non-best in a non-primary cluster -> Secondary+Supplementary
            let best_score = cluster_best_score.get(&cluster).copied();
            let mapq = compute_mapq(score, best_score, aligned_len, identity);
            classified.push(ClassifiedAlignment {
                mapq,
                class: AlignmentClass::SecondarySupplementary,
                candidate,
            });
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

    // Maximum occurrences for a seed to be used (filters highly repetitive k-mers)
    const MAX_SEED_OCCURRENCES: usize = 50;

    // Reusable buffers for seeding
    let mut hits: Vec<SeedHit> = Vec::new();
    let mut hit_vec: Vec<(usize, usize)> = Vec::new();
    let mut merge_scratch: Vec<SeedHit> = Vec::new();
    let mut cuts = Vec::new();
    let mut candidates: Vec<CandidateAlignment> = Vec::new();
    let mut longest_subsequence = LongestSubsequence::default();
    let mut chain_indices = Vec::new();
    let mut chain = Vec::new();
    let max_var = (seq_len as f64 * 0.01).powi(2);

    // Helper closure to process one strand (seeding, merging, clustering, alignment)
    let mut process_strand =
        |strand_seq: &[u8], is_reverse: bool, candidates: &mut Vec<CandidateAlignment>| {
            hits.clear();

            // Phase 1: Collect seed hits using forward-only syncmers
            Kmer::<K>::kmerize_open_syncmers_fwd(strand_seq, [(); S], |pos, kmer| {
                hit_vec.clear();
                index.with(&kmer, |chrom_id, chrom_pos| {
                    hit_vec.push((chrom_id, chrom_pos));
                });
                // Use seeds up to occurrence threshold
                if hit_vec.len() <= MAX_SEED_OCCURRENCES {
                    for &(chrom_id, chrom_pos) in &hit_vec {
                        hits.push(SeedHit::new(chrom_id, chrom_pos, pos, kmer.0, K));
                    }
                }
            });

            let strand_name = if is_reverse { "REV" } else { "FWD" };
            metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase()))
                .record(hits.len() as f64);

            // Phase 2: Sort hits - SeedHit's Ord gives us (chrom_id, diagonal, ref_pos) order
            hits.sort_unstable();

            // Phase 3: Merge overlapping/adjacent hits on same diagonal
            merge_scratch.clear();
            for hit in hits.drain(..) {
                if let Some(last) = merge_scratch.last_mut() {
                    if last
                        .extend(hit.chrom_id, hit.ref_pos, hit.read_pos, hit.kmer, K)
                        .is_none()
                    {
                        continue; // Successfully merged
                    }
                }
                merge_scratch.push(hit);
            }
            std::mem::swap(&mut hits, &mut merge_scratch);

            // Phase 3b: Extend each seed's exact match bidirectionally
            // This is the minimap2-style extension that maximizes anchor length
            for hit in &mut hits {
                let ref_seq = reference.get_seq(hit.chrom_id, 0, usize::MAX);
                hit.extend_exact(strand_seq, ref_seq);
            }

            // Phase 3c: Remove duplicates created by extension
            // When gaps between seeds were due to filtered repetitive k-mers (not mismatches),
            // both adjacent seeds extend to the same flanking mismatches, producing identical
            // (chrom_id, diagonal, ref_pos, match_len). Since extension preserves sort order
            // (a later seed can only reach an earlier seed's start if they converge to the same
            // position), we just deduplicate adjacent entries.
            merge_scratch.clear();
            for hit in hits.drain(..) {
                if let Some(last) = merge_scratch.last() {
                    // All fields except kmer should match for duplicates
                    if hit.chrom_id == last.chrom_id
                        && hit.diagonal == last.diagonal
                        && hit.ref_pos == last.ref_pos
                        && hit.match_len == last.match_len
                    {
                        continue; // Duplicate, skip
                    }
                }
                merge_scratch.push(hit);
            }
            std::mem::swap(&mut hits, &mut merge_scratch);

            if false {
                // Write debug SAM output for seed hits (if debug file is initialized)
                for hit in &hits {
                    write_debug_sam(&hit.to_sam_line(read_name, is_reverse));
                }
            }

            // Phase 4: Cluster hits by diagonal, then build chains and alignments
            cuts.clear();
            dbscan_variance_aware(&hits, 100, max_var, |hit| hit.diagonal, &mut cuts);
            metrics::histogram!(format!("{}_clusters_count", strand_name.to_lowercase()))
                .record(cuts.len().saturating_sub(1) as f64);

            for i in 1..cuts.len() {
                let begin = cuts[i - 1];
                let end = cuts[i];
                let cluster = &hits[begin..end];

                // Use LIS (increasing ref positions) - same for both strands now!
                longest_subsequence.longest_colinear_chain(
                    cluster,
                    |hit| hit.ref_pos as i64,
                    true,
                    &mut chain_indices,
                );
                chain.clear();
                chain.extend(chain_indices.iter().map(|&i| cluster[i]));
                // Sort by read position to ensure proper order for gap alignment
                chain.sort_by_key(|hit| hit.read_pos);
                metrics::histogram!(format!("{}_chain_length", strand_name.to_lowercase()))
                    .record(chain.len() as f64);

                if let Some(candidate) = build_alignment_from_chain(
                    read_name, &chain, strand_seq, seq_len, reference, is_reverse,
                ) {
                    candidates.push(candidate);
                }
            }
        };

    // Process forward strand
    process_strand(seq, false, &mut candidates);

    // Compute reverse complement and process reverse strand
    let mut rc_seq = Vec::with_capacity(seq_len);
    reverse_complement_into(seq, &mut rc_seq);
    process_strand(&rc_seq, true, &mut candidates);

    log::info!(
        "Read {}: found {} candidate alignments",
        read_name,
        candidates.len(),
    );

    // Classify all candidate alignments
    let classified = classify_alignments(candidates, seq_len);

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

    if false {
        init_debug_sam("seeds.sam")?;
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

    if false {
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
}
