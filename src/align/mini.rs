use crate::align::{
    AlignParams, Alignment, WfAligner, kmer_anchors::select_kmer_anchors,
    lcs_anchors::select_lcs_anchors, wfa::WfaFailure,
};

#[derive(Debug)]
pub enum MiniAlignError {
    NoSeeds,
    WfaFailure {
        error: WfaFailure,
        query_start: usize,
        query_end: usize,
        ref_start: usize,
        ref_end: usize,
    },
    PartialAlignment,
}

impl std::fmt::Display for MiniAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MiniAlignError::NoSeeds => write!(f, "No usable seeds found for alignment."),
            MiniAlignError::WfaFailure {
                error,
                query_start,
                query_end,
                ref_start,
                ref_end,
            } => write!(
                f,
                "WFA alignment failed: {:?} (query: {}-{}, ref: {}-{})",
                error, query_start, query_end, ref_start, ref_end
            ),
            MiniAlignError::PartialAlignment => write!(f, "Alignment did not consume all bases."),
        }
    }
}

impl std::error::Error for MiniAlignError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            MiniAlignError::WfaFailure { error, .. } => Some(error),
            _ => None,
        }
    }
}

const USE_KMER_SEEDS: bool = false;

pub fn align(
    read_seq: &[u8],
    ref_seq: &[u8],
    min_length: usize,
) -> Result<Alignment, MiniAlignError> {
    let anchors = if USE_KMER_SEEDS {
        select_kmer_anchors(read_seq, ref_seq, min_length)
    } else {
        select_lcs_anchors(read_seq, ref_seq, min_length)
    };

    if anchors.is_empty() {
        return Err(MiniAlignError::NoSeeds);
    }

    let mut alignments: Vec<Alignment> = Vec::new();
    let n = anchors.len();

    for i in 0..n {
        let seed = &anchors[i];

        // Align up to this seed
        let (q_start, r_start) = if i == 0 {
            (0, 0)
        } else {
            let prev = &anchors[i - 1];
            (prev.query_pos + prev.length, prev.ref_pos + prev.length)
        };
        let q_end = seed.query_pos;
        let r_end = seed.ref_pos;
        if q_end > q_start || r_end > r_start {
            let query = &read_seq[q_start..q_end];
            let reference = &ref_seq[r_start..r_end];
            let aln = WfAligner::new(AlignParams::default())
                .align(query, reference)
                .map_err(|error| MiniAlignError::WfaFailure {
                    error,
                    query_start: q_start,
                    query_end: q_end,
                    ref_start: r_start,
                    ref_end: r_end,
                })?;
            alignments.push(aln);
        }

        // Add the seed match as a perfect alignment
        let seed_query = &read_seq[seed.query_pos..(seed.query_pos + seed.length)];
        let seed_ref = &ref_seq[seed.ref_pos..(seed.ref_pos + seed.length)];
        assert_eq!(seed_query, seed_ref);
        let seed_aln = Alignment::from_perfect_match(seed_query.len());
        alignments.push(seed_aln);

        if i == n - 1 {
            // For the last seed, align to end of read/ref
            let q_start = seed.query_pos + seed.length;
            let q_end = read_seq.len();
            let r_start = seed.ref_pos + seed.length;
            let r_end = ref_seq.len();
            if q_end > q_start || r_end > r_start {
                let query = &read_seq[q_start..q_end];
                let reference = &ref_seq[r_start..r_end];
                let aln = WfAligner::new(AlignParams::default())
                    .align(query, reference)
                    .map_err(|error| MiniAlignError::WfaFailure {
                        error,
                        query_start: q_start,
                        query_end: q_end,
                        ref_start: r_start,
                        ref_end: r_end,
                    })?;
                alignments.push(aln);
            }
        }
    }

    assert!(!alignments.is_empty());

    let mut final_alignment = Alignment::concat(&alignments);
    final_alignment.normalize();

    let query_consumed = final_alignment.query_consumed();
    let ref_consumed = final_alignment.reference_consumed();
    if query_consumed != read_seq.len() || ref_consumed != ref_seq.len() {
        log::error!(
            "MiniAligner produced partial alignment: query_consumed={}, ref_consumed={}, query_len={}, ref_len={} with seeds={}",
            query_consumed,
            ref_consumed,
            read_seq.len(),
            ref_seq.len(),
            anchors.len(),
        );
        return Err(MiniAlignError::PartialAlignment);
    }

    log::debug!(
        "MiniAligner produced alignment: score={}, query_len={}, ref_len={}, seeds={}",
        final_alignment.score,
        read_seq.len(),
        ref_seq.len(),
        anchors.len(),
    );

    Ok(final_alignment)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn align_test_1() {
        let query = b"ATTTTTATATATACTTATATATTTATATATATTTTTATATATACTCATATATTTATATATATTTTATATATACTTATTTATATATATATATTTTTATATATATTTAATTTTTACATATATTTATATTTTTATATATTTATATATTTATATATTTTTATATTTTATATATATGTTTATATATTTATATATTATATATATTTATATATATTTATATATTTATATATTATATATTTATATATATTTATATATTTATATATTATATATATTTATATATTTATATATTTATATATTACATATATTTATATATATTTATATATTTATATATGTTTATATATTTATATATTATATATATTTATATATATTTATATTATATATATACTTATATATTTATATATATTTTTATATATACTTATATATTTATATATATTTTTATATATACTTATATATTTATATATATTTTTATATATACTTATATATATTTTTTATATATTTATATATTTTTATATATATTTAATTTTTAT";
        let reference = b"TTTTATATATACTTATATATTTATATATATTTTTATATATACTCATATATTTATATATATTTTATATATACTTATTTTATATATATATATTTTTATATATATTTAATTTTTAC";
        let result = align(query, reference, 11);
        let alignment = result.expect("alignment failed");
        alignment
            .validate(reference, query, 0)
            .expect("Alignment validation failed.");
        let (r, a, q) = alignment.blast_style(reference, query);
        println!("REF: {}", r);
        println!("ALN: {}", a);
        println!("QRY: {}", q);
        if USE_KMER_SEEDS {
               assert_eq!(alignment.score, 838);
        } else {
               assert_eq!(alignment.score, 842);       
        }
    }

    #[test]
    fn align_test_2() {
        let query = b"ACCTGACTGTCAGAAGGAAAACTAACAAACAGAAAGGAATAGCATCAACATCAACAAAAAGGACATCCCACACCAAAACCCCATCTGTGGGTCACCATCATCGAAGACCAAAGGTAGATAAAACCACAAAGATGGGGAGAAACCAGAGCAGAAAAGCTGAAAATTCTAAAAACCACAGCATCCCTTCTCCTGCAAAGGATTACAGCTCCTTGCCAGCAACGGAATAAAGCTGGACGGAGAATGACTTTGATGAGCTGACAGAAGTAGGCTTCAGAAGGTTGGTAATAACAAACTTCTCTGAGTTAAGGGAGGTATTCGAACCCATCACAAGGAAGCTAAAAGCCTTGCAAACAGATTAGATGAATGGCTAACTAGAATAAGCAGTGTGGAGACCTTAAATGACCTGATGGAGCTGAAAACCATGGAACAAGAACTATGTGACACATGCACAAGCTTCAGTAGTTGATTCAATAAACTGAAAGAAAGGGTATCAGTGATTGAAGATCAAATT";
        let reference = b"TCCTGACTGTTAAAAGGAAAACTAGCAAACAGAAAGGACATCCACACCAAAACCCCATCTGTATGTCACCATCATCAAAGACCAAAGGTAGATAAAACCACAAAGATGGGGAAAAAACAGAACAGAAAAAACTGAAAATTCTAAAAATCAGAGCTCCTCTCCTCCTCCAAAGGAACACAGCTCCTCACCAGCAATGGAACAAAGCTGGACAGAGAATGACTTTGACGAGTTGAGAGAAGAAGGCTTCAGACAATCAAACTTCTCTGAGCTAAAGGAGGAAGTTTGAACCCATGGCAAAGAAGTTAAAAACCTTGAAAAAAGATTAGATGAATGGCTAACTAGAATAAGCAATGCCAAGAAGTCCTTAAAGGACCTAATGGAGCTGAAAACCACAGAACGAGAACTACATGACGAATGCACAAGCTTCAGTAGCCGATTCGATCAACTGGAAGAAAGGGTGTCAGTGATTGAAGATCAAATGAATGA";
        let result = align(query, reference, 15);
        let alignment = result.expect("alignment failed");
        alignment
            .validate(reference, query, 0)
            .expect("Alignment validation failed.");
        let (r, a, q) = alignment.blast_style(reference, query);
        println!("REF: {}", r);
        println!("ALN: {}", a);
        println!("QRY: {}", q);
        assert_eq!(alignment.score, 388);
    }
}
