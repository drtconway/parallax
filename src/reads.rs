use std::sync::Arc;

use crate::align::{
    Alignment, CigarOp, ContextAwareParams, ContextAwareScore, align, context_aware_score,
};
use crate::error::Result;
use crate::index::Index;
use crate::kmers::Kmer;
use crate::reference::InMemoryReference;
use crate::utils::{dbscan_variance_aware, longest_colinear_chain};
use crate::writer::AlignmentWriter;

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

/// A seed hit representing a k-mer match between read and reference
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SeedHit {
    /// Chromosome/contig index in the reference
    chrom_id: usize,
    /// Diagonal: ref_pos - read_pos (constant for colinear matches)
    diagonal: i64,
    /// Position in the reference sequence
    ref_pos: usize,
    /// Position in the read sequence
    read_pos: usize,
    /// Initial kmer
    kmer: u64,
    /// Length of the match (initially k, may be extended)
    match_len: usize,
}

impl SeedHit {
    /// Create a new seed hit
    fn new(chrom_id: usize, ref_pos: usize, read_pos: usize, kmer: u64, match_len: usize) -> Self {
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
    fn read_end(&self) -> usize {
        self.read_pos + self.match_len
    }

    /// End position in the reference
    fn ref_end(&self) -> usize {
        self.ref_pos + self.match_len
    }

    /// Attempt to the seed hit if the new k-mer extends the current match
    /// or return a new seed hit if the new k-mer does not overlap.
    fn extend(
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

    /// Calculate mapping quality (rough approximation)
    fn mapq(&self, read_len: usize, is_unique: bool) -> u8 {
        let coverage = self.read_coverage(read_len);
        let identity = self.identity();

        // Base quality from identity and coverage
        let base_q = (identity * coverage * 60.0) as u8;

        // Reduce quality if not unique
        if is_unique {
            base_q.min(60)
        } else {
            base_q.min(30)
        }
    }

    /// Check if this alignment overlaps another on the read
    fn read_overlaps(&self, other: &CandidateAlignment) -> bool {
        self.read_start < other.read_end && other.read_start < self.read_end
    }

    /// Calculate read overlap fraction with another alignment
    fn read_overlap_fraction(&self, other: &CandidateAlignment, read_len: usize) -> f64 {
        let overlap_start = self.read_start.max(other.read_start);
        let overlap_end = self.read_end.min(other.read_end);
        if overlap_start >= overlap_end {
            0.0
        } else {
            (overlap_end - overlap_start) as f64 / read_len as f64
        }
    }
}

/// Classification of an alignment
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum AlignmentClass {
    Primary,
    Secondary,
    Supplementary,
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
            _ => {}
        }
        flag
    }
}

/// Complement a single nucleotide
#[inline]
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Reverse complement a sequence into the provided buffer
fn reverse_complement_into(seq: &[u8], buf: &mut Vec<u8>) {
    buf.clear();
    buf.reserve(seq.len());
    for &base in seq.iter().rev() {
        buf.push(complement(base));
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
    if chain.len() < 2 {
        return None;
    }

    if false {
        log::info!(
            "Building alignment for read {} on {} strand with {} seeds:",
            read_id,
            if is_reverse { "reverse" } else { "forward" },
            chain.len()
        );
        log::info!("Seed:\tchrom\tpos\tread\tlen");
        for i in 1..chain.len() {
            let prev = &chain[i - 1];
            let hit = &chain[i];

            log::info!(
                "  Seed:\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}",
                hit.chrom_id,
                hit.ref_pos as i64 - prev.ref_pos as i64,
                hit.read_pos as i64 - prev.read_pos as i64,
                hit.match_len,
                (hit.ref_pos as i64 - prev.ref_pos as i64) as f64 / 20.0,
                (hit.read_pos as i64 - prev.read_pos as i64) as f64 / 20.0,
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

/// Classify candidate alignments into primary, secondary, supplementary, and low quality.
///
/// Classification rules:
/// 1. Primary: The best alignment by score that covers a reasonable fraction of the read
/// 2. Supplementary: Other high-quality alignments that don't overlap the primary on the read
///    (indicating a chimeric read)
/// 3. Secondary: Alternative alignments that overlap with primary on the read
///    (indicating multi-mapping)
/// 4. Low Quality: Alignments below score/coverage thresholds
fn classify_alignments(
    mut candidates: Vec<CandidateAlignment>,
    read_len: usize,
) -> Vec<ClassifiedAlignment> {
    if candidates.is_empty() {
        return Vec::new();
    }

    // Sort by context-aware score (lower is better) then by coverage
    candidates.sort_by(|a, b| {
        a.context_score
            .score
            .cmp(&b.context_score.score)
            .then_with(|| {
                // Higher coverage is better
                let cov_a = a.read_coverage(read_len);
                let cov_b = b.read_coverage(read_len);
                cov_b
                    .partial_cmp(&cov_a)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    let mut classified = Vec::with_capacity(candidates.len());

    // First pass: identify primary alignment (best score that meets quality criteria)
    let primary_idx = candidates
        .iter()
        .enumerate()
        .find(|(_, c)| {
            c.read_coverage(read_len) >= MIN_READ_COVERAGE
                && c.identity() >= MIN_ALIGNMENT_IDENTITY
                && c.score_per_base() <= MAX_SCORE_PER_BASE
        })
        .map(|(i, _)| i);

    let primary = primary_idx.map(|i| candidates[i].clone());

    for (i, candidate) in candidates.into_iter().enumerate() {
        let coverage = candidate.read_coverage(read_len);
        let identity = candidate.identity();
        let score_per_base = candidate.score_per_base();

        // Check if this is low quality using multiple criteria
        let is_low_quality = coverage < MIN_READ_COVERAGE
            || identity < MIN_ALIGNMENT_IDENTITY
            || score_per_base > MAX_SCORE_PER_BASE;

        if is_low_quality {
            classified.push(ClassifiedAlignment {
                mapq: 0,
                class: AlignmentClass::LowQuality,
                candidate,
            });
            continue;
        }

        // Primary alignment
        if Some(i) == primary_idx {
            let is_unique = classified
                .iter()
                .filter(|c| c.class != AlignmentClass::LowQuality)
                .count()
                == 0;
            classified.push(ClassifiedAlignment {
                mapq: candidate.mapq(read_len, is_unique),
                class: AlignmentClass::Primary,
                candidate,
            });
            continue;
        }

        // Compare with primary to determine secondary vs supplementary
        if let Some(ref prim) = primary {
            let overlap = candidate.read_overlap_fraction(prim, read_len);

            if overlap < 0.1 {
                // Low overlap with primary - this is a supplementary (chimeric) alignment
                classified.push(ClassifiedAlignment {
                    mapq: candidate.mapq(read_len, false),
                    class: AlignmentClass::Supplementary,
                    candidate,
                });
            } else {
                // Overlaps with primary - this is a secondary (multi-mapping) alignment
                classified.push(ClassifiedAlignment {
                    mapq: candidate.mapq(read_len, false),
                    class: AlignmentClass::Secondary,
                    candidate,
                });
            }
        } else {
            // No primary, so this is secondary
            classified.push(ClassifiedAlignment {
                mapq: candidate.mapq(read_len, false),
                class: AlignmentClass::Secondary,
                candidate,
            });
        }
    }

    // Sort so primary comes first, then supplementary, then secondary
    classified.sort_by_key(|c| match c.class {
        AlignmentClass::Primary => 0,
        AlignmentClass::Supplementary => 1,
        AlignmentClass::Secondary => 2,
        AlignmentClass::LowQuality => 3,
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
    let max_var = (seq_len as f64 * 0.01).powi(2);

    // Helper closure to process one strand (seeding, merging, clustering, alignment)
    let mut process_strand = |strand_seq: &[u8], is_reverse: bool, candidates: &mut Vec<CandidateAlignment>| {
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
        metrics::histogram!(format!("{}_hits_count", strand_name.to_lowercase())).record(hits.len() as f64);

        // Phase 2: Sort hits - SeedHit's Ord gives us (chrom_id, diagonal, ref_pos) order
        hits.sort_unstable();

        // Phase 3: Merge overlapping/adjacent hits on same diagonal
        merge_scratch.clear();
        for hit in hits.drain(..) {
            if let Some(last) = merge_scratch.last_mut() {
                if last.extend(hit.chrom_id, hit.ref_pos, hit.read_pos, hit.kmer, K).is_none() {
                    continue; // Successfully merged
                }
            }
            merge_scratch.push(hit);
        }
        std::mem::swap(&mut hits, &mut merge_scratch);

        log::info!("{}:\tchrom\tdiag\tref\tread\tlen\tkmer", strand_name);
        for hit in &hits {
            log::info!("{}:\t{}\t{}\t{}\t{}\t{}\t{}", strand_name, hit.chrom_id, hit.diagonal, hit.ref_pos, hit.read_pos, hit.match_len, Kmer::<K>(hit.kmer).to_string());
        }

        // Phase 4: Cluster hits by diagonal, then build chains and alignments
        cuts.clear();
        dbscan_variance_aware(&hits, 100, max_var, |hit| hit.diagonal, &mut cuts);
        metrics::histogram!(format!("{}_clusters_count", strand_name.to_lowercase())).record(cuts.len().saturating_sub(1) as f64);
        
        for i in 1..cuts.len() {
            let begin = cuts[i - 1];
            let end = cuts[i];
            let cluster = &hits[begin..end];

            // Use LIS (increasing ref positions) - same for both strands now!
            let chain_indices = longest_colinear_chain(cluster, |hit| hit.ref_pos as i64, true);
            let mut chain: Vec<_> = chain_indices.iter().map(|&i| cluster[i]).collect();
            // Sort by read position to ensure proper order for gap alignment
            chain.sort_by_key(|hit| hit.read_pos);
            metrics::histogram!(format!("{}_chain_length", strand_name.to_lowercase())).record(chain.len() as f64);

            if let Some(candidate) = build_alignment_from_chain(
                read_name,
                &chain,
                strand_seq,
                seq_len,
                reference,
                is_reverse,
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
    writer.write_command_header(&format!("parallax index <reference> {}", input_file))?;

    Ok(())
}

/// Process reads from a FASTQ file
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
        
        assert!(result.is_none(), "Should extend in place, not return new hit");
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
        
        assert!(result.is_some(), "Different chromosome should create new hit");
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
        
        assert!(result.is_some(), "Backwards ref position should create new hit");
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_backwards_read_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;
        
        // New read_pos before current read_pos (even if same diagonal)
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, k);
        
        assert!(result.is_some(), "Backwards read position should create new hit");
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
