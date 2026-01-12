use crate::align::{align, Alignment, CigarOp};
use crate::error::Result;
use crate::index::Index;
use crate::kmers::Kmer;
use crate::utils::{Selection, dbscan_variance_aware, longest_colinear_chain};

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
/// Returns (chrom_id, ref_start, ref_end, read_start, read_end, is_reverse, alignment)
/// 
/// For reverse-strand alignments, the read slices are reverse-complemented before
/// aligning to the forward reference.
fn build_alignment_from_chain<const K: usize, const S: usize>(
    chain: &[(usize, i64, usize, usize, usize)],
    seq: &[u8],
    seq_len: usize,
    index: &Index<K, S>,
    is_reverse: bool,
    rc_buf: &mut Vec<u8>,
) -> Option<(usize, usize, usize, usize, usize, bool, Alignment)> {
    if chain.len() < 2 {
        return None;
    }

    let chrom_id = chain[0].0;
    let mut full_cigar: Vec<CigarOp> = Vec::new();
    let mut total_score = 0i32;

    // Compute alignment span
    let first = chain.first().unwrap();
    let last = chain.last().unwrap();
    
    // For reverse strand, ref positions decrease as read positions increase,
    // so the ref_start is actually from the last chain element
    let (ref_start, ref_end) = if is_reverse {
        (last.2, first.2 + first.4)
    } else {
        (first.2, last.2 + last.4)
    };
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
                let ref_slice = index.get_seq(chrom_id, ref_gap_start, ref_gap_end);
                let read_slice = &seq[read_gap_start..read_gap_end];

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

    Some((chrom_id, ref_start, ref_end, read_start, read_end, is_reverse, alignment))
}

pub fn process_reads<const K: usize, const S: usize>(
    index: &Index<K, S>,
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

        // Process forward strand hits
        dbscan_variance_aware(&fwd_hits, 100, max_var, |hit| hit.1, &mut cuts);
        for i in 1..cuts.len() {
            let begin = cuts[i - 1];
            let end = cuts[i];
            let cluster = &fwd_hits[begin..end];

            let chain_indices = longest_colinear_chain(cluster, |hit| hit.2 as i64, true);
            let chain: Vec<_> = chain_indices.iter().map(|&i| cluster[i]).collect();

            if let Some((chrom_id, ref_start, ref_end, read_start, read_end, _is_rev, alignment)) =
                build_alignment_from_chain(&chain, seq, seq_len, index, false, &mut rc_buf)
            {
                log::info!(
                    "Read {}: FWD align to {}:{}-{} (read {}..{}), score={}, CIGAR={}",
                    std::str::from_utf8(record.name()).unwrap_or("?"),
                    index.chrom_name(chrom_id),
                    ref_start,
                    ref_end,
                    read_start,
                    read_end,
                    alignment.score,
                    alignment.cigar_string(),
                );
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
            let chain: Vec<_> = chain_indices.iter().map(|&i| cluster[i]).collect();

            if let Some((chrom_id, ref_start, ref_end, read_start, read_end, _is_rev, alignment)) =
                build_alignment_from_chain(&chain, seq, seq_len, index, true, &mut rc_buf)
            {
                log::info!(
                    "Read {}: REV align to {}:{}-{} (read {}..{}), score={}, CIGAR={}",
                    std::str::from_utf8(record.name()).unwrap_or("?"),
                    index.chrom_name(chrom_id),
                    ref_start,
                    ref_end,
                    read_start,
                    read_end,
                    alignment.score,
                    alignment.cigar_string(),
                );
            }
        }
    }

    Ok(())
}
