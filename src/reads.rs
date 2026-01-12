use crate::align::{Alignment, CigarOp, align};
use crate::error::Result;
use crate::index::Index;
use crate::kmers::Kmer;
use crate::reference::Reference;
use crate::utils::{Selection, dbscan_variance_aware, longest_colinear_chain};

/// SAM flags
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Minimum alignment score threshold (alignments below this are considered low quality)
const MIN_ALIGNMENT_SCORE: i32 = 500;

/// Minimum fraction of read covered for a valid alignment
const MIN_READ_COVERAGE: f64 = 0.1;

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
}

impl CandidateAlignment {
    /// Calculate the fraction of the read covered by this alignment
    fn read_coverage(&self, read_len: usize) -> f64 {
        (self.read_end - self.read_start) as f64 / read_len as f64
    }

    /// Calculate alignment identity (matches / aligned length)
    fn identity(&self) -> f64 {
        let mut matches = 0u64;
        let mut aligned = 0u64;
        for op in &self.alignment.cigar {
            match op {
                CigarOp::Match(n) => {
                    matches += *n as u64;
                    aligned += *n as u64;
                }
                CigarOp::Mismatch(n) => {
                    aligned += *n as u64;
                }
                CigarOp::Ins(n) | CigarOp::Del(n) => {
                    aligned += *n as u64;
                }
                CigarOp::SoftClip(_) => {}
            }
        }
        if aligned == 0 {
            0.0
        } else {
            matches as f64 / aligned as f64
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
/// Build alignment from a chain of seed matches, filling gaps with WFA.
///
/// For reverse-strand alignments, the read slices are reverse-complemented before
/// aligning to the forward reference.
fn build_alignment_from_chain(
    chain: &[(usize, i64, usize, usize, usize)],
    seq: &[u8],
    seq_len: usize,
    reference: &mut Reference,
    is_reverse: bool,
    rc_buf: &mut Vec<u8>,
    ref_buf: &mut Vec<u8>,
) -> Option<CandidateAlignment> {
    if chain.len() < 2 {
        return None;
    }

    let chrom_id = chain[0].0;
    let mut full_cigar: Vec<CigarOp> = Vec::new();
    let mut total_score = 0i32;

    // Compute alignment span from actual min/max reference positions
    let first = chain.first().unwrap();
    let last = chain.last().unwrap();

    // Use actual min/max ref positions to handle any chain ordering
    let ref_start = chain.iter().map(|h| h.2).min().unwrap();
    let ref_end = chain.iter().map(|h| h.2 + h.4).max().unwrap();
    let read_start = first.3;
    let read_end = last.3 + last.4;

    // Add soft-clip for unaligned prefix
    if read_start > 0 {
        full_cigar.push(CigarOp::SoftClip(read_start as u32));
    }

    for j in 0..chain.len() {
        let (_cid, _d, ref_pos, read_pos, match_len) = chain[j];

        // Align gap before this seed (if not first seed)
        if j > 0 {
            let prev = chain[j - 1];
            let prev_read_end = prev.3 + prev.4;
            let read_gap_start = prev_read_end;
            let read_gap_end = read_pos;

            // Reference gap depends on strand
            let (ref_gap_start, ref_gap_end) = if is_reverse {
                // Reverse: previous ref_pos is higher, current is lower
                // Gap is from (ref_pos + match_len) to prev.2
                (ref_pos + match_len, prev.2)
            } else {
                let prev_ref_end = prev.2 + prev.4;
                (prev_ref_end, ref_pos)
            };

            if ref_gap_end > ref_gap_start || read_gap_end > read_gap_start {
                // Handle overlapping chain elements - clamp gaps to avoid negative ranges
                let actual_read_start = read_gap_start.min(read_gap_end);
                let actual_read_end = read_gap_start.max(read_gap_end);
                let actual_ref_start = ref_gap_start.min(ref_gap_end);
                let actual_ref_end = ref_gap_start.max(ref_gap_end);

                // Fetch reference sequence into buffer
                if reference
                    .get_seq_into(chrom_id, actual_ref_start, actual_ref_end, ref_buf)
                    .is_err()
                {
                    continue;
                }
                let ref_slice = ref_buf.as_slice();
                let read_slice = &seq[actual_read_start..actual_read_end];

                // For reverse strand, reverse-complement the read slice
                let query_slice: &[u8] = if is_reverse {
                    reverse_complement_into(read_slice, rc_buf);
                    rc_buf.as_slice()
                } else {
                    read_slice
                };

                if let Some(aln) = align(query_slice, ref_slice) {
                    total_score += aln.score;
                    // For reverse strand, we need to reverse the CIGAR operations
                    // since we're building CIGAR in read order but aligned in rev-comp
                    if is_reverse {
                        for op in aln.cigar.into_iter().rev() {
                            full_cigar.push(op);
                        }
                    } else {
                        full_cigar.extend(aln.cigar);
                    }
                } else {
                    // Alignment failed, emit as insertions/deletions
                    if !read_slice.is_empty() {
                        full_cigar.push(CigarOp::Ins(read_slice.len() as u32));
                    }
                    if !ref_slice.is_empty() {
                        full_cigar.push(CigarOp::Del(ref_slice.len() as u32));
                    }
                }
            }
        }

        // Add the seed match itself
        full_cigar.push(CigarOp::Match(match_len as u32));
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

    Some(CandidateAlignment {
        chrom_id,
        ref_start,
        ref_end,
        read_start,
        read_end,
        is_reverse,
        alignment,
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

    // Sort by score (lower is better in our scoring system) then by coverage
    candidates.sort_by(|a, b| {
        a.alignment.score.cmp(&b.alignment.score).then_with(|| {
            // Higher coverage is better
            let cov_a = a.read_coverage(read_len);
            let cov_b = b.read_coverage(read_len);
            cov_b
                .partial_cmp(&cov_a)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    let mut classified = Vec::with_capacity(candidates.len());

    // First pass: identify primary alignment
    let primary_idx = candidates
        .iter()
        .enumerate()
        .find(|(_, c)| c.read_coverage(read_len) >= MIN_READ_COVERAGE)
        .map(|(i, _)| i);

    let primary = primary_idx.map(|i| candidates[i].clone());

    for (i, candidate) in candidates.into_iter().enumerate() {
        let coverage = candidate.read_coverage(read_len);

        // Check if this is low quality
        if candidate.alignment.score > MIN_ALIGNMENT_SCORE || coverage < MIN_READ_COVERAGE {
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

pub fn process_reads<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &mut Reference,
    fastq: &str,
) -> Result<()> {
    log::info!("Processing reads from {}", fastq);

    let reader = std::fs::File::open(fastq).map(std::io::BufReader::new)?;
    let mut reader = noodles::fastq::io::Reader::new(reader);

    for record in reader.records() {
        let record = record?;
        let seq: &[u8] = record.sequence().as_ref();
        let seq_len = seq.len();

        // Hit tuple: (chrom_id, d, chrom_pos, read_pos, match_len)
        let mut fwd_hits: Vec<(usize, i64, usize, usize, usize)> = Vec::new();
        let mut rev_hits: Vec<(usize, i64, usize, usize, usize)> = Vec::new();
        let mut hit_vec: Vec<(usize, usize)> = Vec::new();

        // Helper to merge or push a hit, extending if it overlaps the last one
        fn merge_or_push(
            hits: &mut Vec<(usize, i64, usize, usize, usize)>,
            chrom_id: usize,
            d: i64,
            chrom_pos: usize,
            read_pos: usize,
            k: usize,
        ) {
            if let Some(last) = hits.last_mut() {
                // Same chrom + diagonal, and overlaps/adjacent in read coords?
                if last.0 == chrom_id && last.1 == d && read_pos < last.3 + last.4 {
                    // Extend match: new end is read_pos + k
                    let new_end = read_pos + k;
                    let old_end = last.3 + last.4;
                    if new_end > old_end {
                        last.4 = new_end - last.3;
                    }
                    return;
                }
            }
            hits.push((chrom_id, d, chrom_pos, read_pos, k));
        }

        for (pos, selection) in Kmer::<K>::open_syncmer_iter(seq, [(); S]) {
            let fwd: Option<Kmer<K>> = match &selection {
                Selection::Left(kmer) => Some(*kmer),
                Selection::Both(kmer, _) => Some(*kmer),
                _ => None,
            };
            if let Some(kmer) = fwd {
                hit_vec.clear();
                index.with(&kmer, |chrom_id, chrom_pos| {
                    hit_vec.push((chrom_id, chrom_pos));
                });
                if hit_vec.len() == 1 {
                    let (chrom_id, chrom_pos) = hit_vec[0];
                    let d = chrom_pos as i64 - pos as i64;
                    merge_or_push(&mut fwd_hits, chrom_id, d, chrom_pos, pos, K);
                }
            }

            let rev: Option<Kmer<K>> = match &selection {
                Selection::Right(kmer) => Some(*kmer),
                Selection::Both(_, kmer) => Some(*kmer),
                _ => None,
            };
            if let Some(kmer) = rev {
                hit_vec.clear();
                index.with(&kmer, |chrom_id, chrom_pos| {
                    hit_vec.push((chrom_id, chrom_pos));
                });
                if hit_vec.len() == 1 {
                    let (chrom_id, chrom_pos) = hit_vec[0];
                    let d = chrom_pos as i64 - pos as i64;
                    merge_or_push(&mut rev_hits, chrom_id, d, chrom_pos, pos, K);
                }
            }
        }

        fwd_hits.sort_unstable();
        rev_hits.sort_unstable();

        let max_var = (seq_len as f64 * 0.01).powi(2);
        let mut cuts = Vec::new();
        let mut rc_buf = Vec::new(); // Buffer for reverse-complement
        let mut ref_buf = Vec::new(); // Buffer for reference sequence
        let mut candidates: Vec<CandidateAlignment> = Vec::new();

        // Process forward strand hits
        dbscan_variance_aware(&fwd_hits, 100, max_var, |hit| hit.1, &mut cuts);
        for i in 1..cuts.len() {
            let begin = cuts[i - 1];
            let end = cuts[i];
            let cluster = &fwd_hits[begin..end];

            let chain_indices = longest_colinear_chain(cluster, |hit| hit.2 as i64, true);
            let mut chain: Vec<_> = chain_indices.iter().map(|&i| cluster[i]).collect();
            // Sort by read position to ensure proper order for gap alignment
            chain.sort_by_key(|hit| hit.3);

            if let Some(candidate) = build_alignment_from_chain(
                &chain,
                seq,
                seq_len,
                reference,
                false,
                &mut rc_buf,
                &mut ref_buf,
            ) {
                candidates.push(candidate);
            }
        }

        // Process reverse strand hits
        cuts.clear();
        dbscan_variance_aware(&rev_hits, 100, max_var, |hit| hit.1, &mut cuts);
        for i in 1..cuts.len() {
            let begin = cuts[i - 1];
            let end = cuts[i];
            let cluster = &rev_hits[begin..end];

            // For reverse strand, we use LDS (decreasing ref positions as read position increases)
            let chain_indices = longest_colinear_chain(cluster, |hit| hit.2 as i64, false);
            let mut chain: Vec<_> = chain_indices.iter().map(|&i| cluster[i]).collect();
            // Sort by read position to ensure proper order for gap alignment
            chain.sort_by_key(|hit| hit.3);

            if let Some(candidate) = build_alignment_from_chain(
                &chain,
                seq,
                seq_len,
                reference,
                true,
                &mut rc_buf,
                &mut ref_buf,
            ) {
                candidates.push(candidate);
            }
        }

        // Classify all candidate alignments
        let classified = classify_alignments(candidates, seq_len);
        let read_name = std::str::from_utf8(record.name()).unwrap_or("?");

        if classified.is_empty() {
            // Output unmapped read
            println!(
                "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t*",
                read_name,
                FLAG_UNMAPPED,
                std::str::from_utf8(seq).unwrap_or("*"),
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
                    "Read {}: {} {} to {}:{}-{} (read {}..{}), mapq={}, score={}, CIGAR={}",
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
                    aln.candidate.alignment.cigar_string(),
                );

                // Output SAM record (skip low quality unless there's no primary)
                if aln.class != AlignmentClass::LowQuality {
                    let chrom_name = reference.chrom_name(aln.candidate.chrom_id);
                    let pos = aln.candidate.ref_start + 1; // SAM is 1-based
                    let cigar = aln.candidate.alignment.cigar_string();

                    // For reverse strand, we should output reverse-complemented sequence
                    // but for secondary/supplementary, SAM spec says use * or same as primary
                    let seq_str = if aln.class == AlignmentClass::Primary {
                        if aln.candidate.is_reverse {
                            rc_buf.clear();
                            reverse_complement_into(seq, &mut rc_buf);
                            String::from_utf8_lossy(&rc_buf).into_owned()
                        } else {
                            String::from_utf8_lossy(seq).into_owned()
                        }
                    } else {
                        "*".to_string()
                    };

                    let qual = if aln.class == AlignmentClass::Primary {
                        std::str::from_utf8(record.quality_scores().as_ref())
                            .map(|s| {
                                if aln.candidate.is_reverse {
                                    s.chars().rev().collect::<String>()
                                } else {
                                    s.to_string()
                                }
                            })
                            .unwrap_or_else(|_| "*".to_string())
                    } else {
                        "*".to_string()
                    };

                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
                        read_name,
                        aln.sam_flag(),
                        chrom_name,
                        pos,
                        aln.mapq,
                        cigar,
                        seq_str,
                        qual,
                    );
                }
            }
        }
    }

    Ok(())
}
