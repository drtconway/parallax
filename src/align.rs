use crate::config;

pub mod block;
pub mod kmer_anchors;
pub mod lcs_anchors;
pub mod mini;
pub mod wfa;

pub use wfa::WfAligner;

/// Alignment scoring parameters
#[derive(Clone, Copy, Debug)]
pub struct AlignParams {
    pub match_score: i32,
    /// Mismatch penalty (positive value)
    pub mismatch: i32,
    /// Gap open penalty (positive value)
    pub gap_open: i32,
    /// Gap extend penalty (positive value)
    pub gap_extend: i32,
}

impl Default for AlignParams {
    fn default() -> Self {
        let cfg = config::get();
        Self {
            match_score: cfg.alignment.match_score,
            mismatch: cfg.alignment.mismatch,
            gap_open: cfg.alignment.gap_open,
            gap_extend: cfg.alignment.gap_extend,
        }
    }
}

impl AlignParams {
    /// Score an individual alignment operation
    pub fn score_op(&self, op: &CigarOp) -> i32 {
        match op {
            CigarOp::Match(n) => self.match_score * (*n as i32),
            CigarOp::Mismatch(n) => -self.mismatch * (*n as i32),
            CigarOp::Ins(n) | CigarOp::Del(n) => {
                -self.gap_open - self.gap_extend * (*n as i32)
            }
            CigarOp::SoftClip(_) => 0,
        }
    }
}

#[derive(Debug)]
pub enum AlignmentError {
    BlockError(block::BlockAlignerError),
    WFError(wfa::WfaFailure),
    MiniError(mini::MiniAlignError),
}

impl std::fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentError::BlockError(e) => write!(f, "Block aligner error: {}", e),
            AlignmentError::WFError(e) => write!(f, "WFA alignment error: {}", e),
            AlignmentError::MiniError(e) => write!(f, "Mini alignment error: {}", e),
        }
    }
}

impl std::error::Error for AlignmentError {}

impl From<wfa::WfaFailure> for AlignmentError {
    fn from(err: wfa::WfaFailure) -> Self {
        AlignmentError::WFError(err)
    }
}

impl From<block::BlockAlignerError> for AlignmentError {
    fn from(err: block::BlockAlignerError) -> Self {
        AlignmentError::BlockError(err)
    }
}

impl From<mini::MiniAlignError> for AlignmentError {
    fn from(err: mini::MiniAlignError) -> Self {
        AlignmentError::MiniError(err)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Anchor {
    pub query_pos: usize,
    pub ref_pos: usize,
    pub length: usize,
}

impl Anchor {
    pub fn new(query_pos: usize, ref_pos: usize, length: usize) -> Self {
        Self {
            query_pos,
            ref_pos,
            length,
        }
    }

    pub fn diagonal(&self) -> isize {
        self.ref_pos as isize - self.query_pos as isize
    }

    #[allow(dead_code)]
    fn order_by_length(a: &Anchor, b: &Anchor) -> std::cmp::Ordering {
        let res = b.length.cmp(&a.length);
        if res == std::cmp::Ordering::Equal {
            a.query_pos.cmp(&b.query_pos)
        } else {
            res
        }
    }

    fn order_by_query_pos(a: &Anchor, b: &Anchor) -> std::cmp::Ordering {
        a.query_pos
            .cmp(&b.query_pos)
            .then(b.ref_pos.cmp(&a.ref_pos))
    }
}

/// CIGAR operation
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),    // '='
    Mismatch(u32), // 'X'
    Ins(u32),      // 'I' - insertion in query
    Del(u32),      // 'D' - deletion in query
    SoftClip(u32), // 'S' - soft clipped bases
}

impl CigarOp {
    /// Format as SAM-style CIGAR string
    pub fn to_string(&self) -> String {
        match self {
            CigarOp::Match(n) => format!("{}=", n),
            CigarOp::Mismatch(n) => format!("{}X", n),
            CigarOp::Ins(n) => format!("{}I", n),
            CigarOp::Del(n) => format!("{}D", n),
            CigarOp::SoftClip(n) => format!("{}S", n),
        }
    }


    pub fn score(&self, params: &AlignParams) -> f64 {
        params.score_op(self) as f64
    }

    #[allow(dead_code)]
    pub fn make(cig: &str) -> Option<Vec<CigarOp>> {
        let mut cigar = Vec::new();
        let mut count = 0;

        for c in cig.chars() {
            if let Some(d) = c.to_digit(10) {
                count = count * 10 + d;
            } else {
                let op = match c {
                    '=' => CigarOp::Match(count),
                    'X' => CigarOp::Mismatch(count),
                    'I' => CigarOp::Ins(count),
                    'D' => CigarOp::Del(count),
                    'S' => CigarOp::SoftClip(count),
                    _ => return None,
                };
                cigar.push(op);
                count = 0;
            }
        }

        Some(cigar)
    }
}

/// Alignment result
#[derive(Clone, Debug)]
pub struct Alignment {
    /// Alignment score (edit distance style - lower is better)
    pub score: i32,
    /// CIGAR operations
    pub cigar: Vec<CigarOp>,
}

impl Alignment {
    /// Create a new perfect match alignment
    pub fn from_perfect_match(length: usize) -> Self {
        Self {
            score: 0,
            cigar: vec![CigarOp::Match(length as u32)],
        }
    }

    /// Format CIGAR as string
    pub fn cigar_string(&self) -> String {
        self.cigar.iter().map(|op| op.to_string()).collect()
    }

    /// Format CIGAR in the basic format merging = and X into M
    /// (e.g., 10M1I5M2D3M)
    #[allow(dead_code)]
    pub fn basic_cigar_string(&self) -> String {
        let mut merged: Vec<CigarOp> = Vec::new();
        for op in &self.cigar {
            match op {
                CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                    if let Some(last) = merged.last_mut() {
                        match last {
                            CigarOp::Match(m) | CigarOp::Mismatch(m) => *m += n,
                            _ => merged.push(CigarOp::Match(*n)),
                        }
                    } else {
                        merged.push(CigarOp::Match(*n));
                    }
                }
                _ => merged.push(*op),
            }
        }
        merged
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => format!("{}M", n),
                CigarOp::Mismatch(n) => format!("{}M", n),
                CigarOp::Ins(n) => format!("{}I", n),
                CigarOp::Del(n) => format!("{}D", n),
                CigarOp::SoftClip(n) => format!("{}S", n),
            })
            .collect()
    }

    /// Compute the query (read) length consumed by this CIGAR.
    /// This is the sum of M, I, S, =, X operations.
    /// For valid SAM, this must equal the length of the SEQ field.
    pub fn query_length(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as usize,
                CigarOp::Mismatch(n) => *n as usize,
                CigarOp::Ins(n) => *n as usize,
                CigarOp::SoftClip(n) => *n as usize,
                CigarOp::Del(_) => 0,
            })
            .sum()
    }

    /// Compute the reference span consumed by this CIGAR.
    /// This is the sum of M, D, N, =, X operations.
    #[allow(dead_code)]
    pub fn reference_span(&self) -> u64 {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as u64,
                CigarOp::Mismatch(n) => *n as u64,
                CigarOp::Del(n) => *n as u64,
                CigarOp::Ins(_) => 0,
                CigarOp::SoftClip(_) => 0,
            })
            .sum()
    }

    /// Compute the query (read) bases consumed by alignment operations (excluding soft clips).
    /// This is the sum of M, I, =, X operations only.
    #[allow(dead_code)]
    pub fn query_consumed(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as usize,
                CigarOp::Mismatch(n) => *n as usize,
                CigarOp::Ins(n) => *n as usize,
                CigarOp::Del(_) => 0,
                CigarOp::SoftClip(_) => 0,
            })
            .sum()
    }

    /// Compute the reference bases consumed by alignment operations.
    /// This is the sum of M, D, =, X operations.
    #[allow(dead_code)]
    pub fn reference_consumed(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as usize,
                CigarOp::Mismatch(n) => *n as usize,
                CigarOp::Del(n) => *n as usize,
                CigarOp::Ins(_) => 0,
                CigarOp::SoftClip(_) => 0,
            })
            .sum()
    }

    pub fn total_identity(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as usize,
                _ => 0,
            })
            .sum()
    }

    /// Merge adjacent operations of same type
    pub fn normalize(&mut self) {
        if self.cigar.is_empty() {
            return;
        }
        let mut merged = Vec::with_capacity(self.cigar.len());
        for op in self.cigar.drain(..) {
            if let Some(last) = merged.last_mut() {
                match (last, op) {
                    (CigarOp::Match(n), CigarOp::Match(m)) => *n += m,
                    (CigarOp::Mismatch(n), CigarOp::Mismatch(m)) => *n += m,
                    (CigarOp::Ins(n), CigarOp::Ins(m)) => *n += m,
                    (CigarOp::Del(n), CigarOp::Del(m)) => *n += m,
                    (CigarOp::SoftClip(n), CigarOp::SoftClip(m)) => *n += m,
                    _ => merged.push(op),
                }
            } else {
                merged.push(op);
            }
        }
        self.cigar = merged;
    }

    /// Format the alignment in BLAST style with three lines:
    /// 1. Reference sequence with '-' for insertions (query has extra bases)
    /// 2. Match line: '|' for matches, ' ' for mismatches/gaps
    /// 3. Query sequence with '-' for deletions (reference has extra bases)
    ///
    /// Soft-clipped bases are not shown in the output.
    ///
    /// Returns (ref_line, match_line, query_line)
    #[allow(dead_code)]
    pub fn blast_style(&self, reference: &[u8], query: &[u8]) -> (String, String, String) {
        let mut ref_line = String::new();
        let mut match_line = String::new();
        let mut query_line = String::new();

        let mut ref_pos = 0usize;
        let mut query_pos = 0usize;

        for op in &self.cigar {
            match op {
                CigarOp::Match(n) => {
                    for _ in 0..*n {
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        ref_line.push(r as char);
                        query_line.push(q as char);
                        match_line.push('|');
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                CigarOp::Mismatch(n) => {
                    for _ in 0..*n {
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        ref_line.push(r as char);
                        query_line.push(q as char);
                        match_line.push(' ');
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                CigarOp::Ins(n) => {
                    // Insertion in query: query has bases, reference has gap
                    for _ in 0..*n {
                        let q = query[query_pos];
                        ref_line.push('-');
                        query_line.push(q as char);
                        match_line.push(' ');
                        query_pos += 1;
                    }
                }
                CigarOp::Del(n) => {
                    // Deletion in query: reference has bases, query has gap
                    for _ in 0..*n {
                        let r = reference[ref_pos];
                        ref_line.push(r as char);
                        query_line.push('-');
                        match_line.push(' ');
                        ref_pos += 1;
                    }
                }
                CigarOp::SoftClip(n) => {
                    // Soft clips consume query but aren't shown in alignment
                    query_pos += *n as usize;
                }
            }
        }

        (ref_line, match_line, query_line)
    }

    /// Validate that the CIGAR correctly describes the alignment between query and reference.
    ///
    /// Returns Ok(()) if valid, or Err with a description of the first mismatch found.
    ///
    /// This checks:
    /// - Match operations ('=') actually have matching bases
    /// - Mismatch operations ('X') actually have different bases
    /// - Position tracking is correct
    ///
    /// `query_start` is the 0-based position in the query where the alignment starts
    /// (after any leading soft clips).
    pub fn validate(
        &self,
        reference: &[u8],
        query: &[u8],
        query_start: usize,
    ) -> Result<(), String> {
        let mut ref_pos = 0usize;
        let mut query_pos = query_start;

        // First check the cigar makes sense internally
        let n = self.cigar.len();
        for i in 0..n {
            if let CigarOp::Ins(0)
            | CigarOp::Del(0)
            | CigarOp::Match(0)
            | CigarOp::Mismatch(0)
            | CigarOp::SoftClip(0) = self.cigar[i]
            {
                return Err(format!(
                    "CIGAR op {} has zero length: {:?}",
                    i, self.cigar[i]
                ));
            }
            if i > 0 {
                if std::mem::discriminant(&self.cigar[i])
                    == std::mem::discriminant(&self.cigar[i - 1])
                {
                    return Err(format!(
                        "CIGAR ops {} and {} are adjacent and of same type: {:?}, {:?}",
                        i - 1,
                        i,
                        self.cigar[i - 1],
                        self.cigar[i]
                    ));
                }
            }
            if i > 0 && i < n - 1 {
                if let CigarOp::SoftClip(_) = self.cigar[i] {
                    return Err(format!(
                        "CIGAR op {} is a soft clip in the middle of the alignment: {:?}",
                        i, self.cigar[i]
                    ));
                }
            }
        }

        // Check each CIGAR operation against the sequences
        for (op_idx, op) in self.cigar.iter().enumerate() {
            match op {
                CigarOp::Match(n) => {
                    for i in 0..*n as usize {
                        if ref_pos >= reference.len() {
                            return Err(format!(
                                "CIGAR op {} ({}=): ref_pos {} exceeds reference length {} at offset {}",
                                op_idx,
                                n,
                                ref_pos,
                                reference.len(),
                                i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}=): query_pos {} exceeds query length {} at offset {}",
                                op_idx,
                                n,
                                query_pos,
                                query.len(),
                                i
                            ));
                        }
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        if r != q {
                            return Err(format!(
                                "CIGAR op {} ({}=): expected match at ref_pos {} query_pos {}, but ref='{}' query='{}' (offset {})",
                                op_idx, n, ref_pos, query_pos, r as char, q as char, i
                            ));
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                CigarOp::Mismatch(n) => {
                    for i in 0..*n as usize {
                        if ref_pos >= reference.len() {
                            return Err(format!(
                                "CIGAR op {} ({}X): ref_pos {} exceeds reference length {} at offset {}",
                                op_idx,
                                n,
                                ref_pos,
                                reference.len(),
                                i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}X): query_pos {} exceeds query length {} at offset {}",
                                op_idx,
                                n,
                                query_pos,
                                query.len(),
                                i
                            ));
                        }
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        if r == q {
                            return Err(format!(
                                "CIGAR op {} ({}X): expected mismatch at ref_pos {} query_pos {}, but both are '{}' (offset {})",
                                op_idx, n, ref_pos, query_pos, r as char, i
                            ));
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                CigarOp::Ins(n) => {
                    // Insertion in query: query has bases, reference doesn't
                    let new_pos = query_pos + *n as usize;
                    if new_pos > query.len() {
                        return Err(format!(
                            "CIGAR op {} ({}I): query_pos {} + {} exceeds query length {}",
                            op_idx,
                            n,
                            query_pos,
                            n,
                            query.len()
                        ));
                    }
                    query_pos = new_pos;
                }
                CigarOp::Del(n) => {
                    // Deletion in query: reference has bases, query doesn't
                    let new_pos = ref_pos + *n as usize;
                    if new_pos > reference.len() {
                        return Err(format!(
                            "CIGAR op {} ({}D): ref_pos {} + {} exceeds reference length {}",
                            op_idx,
                            n,
                            ref_pos,
                            n,
                            reference.len()
                        ));
                    }
                    ref_pos = new_pos;
                }
                CigarOp::SoftClip(n) => {
                    // Soft clips only consume query
                    let new_pos = query_pos + *n as usize;
                    if new_pos > query.len() {
                        return Err(format!(
                            "CIGAR op {} ({}S): query_pos {} + {} exceeds query length {}",
                            op_idx,
                            n,
                            query_pos,
                            n,
                            query.len()
                        ));
                    }
                    query_pos = new_pos;
                }
            }
        }

        Ok(())
    }

    /// Compute an information based score.
    pub fn score(&self, params: &AlignParams) -> f64 {
        let mut score = 0.0;
        for op in &self.cigar {
            score += op.score(params) as f64;
        }
        score
    }

    pub fn concat(alignments: &[Alignment]) -> Alignment {
        let mut total_score = 0;
        let mut combined_cigar = Vec::new();

        for aln in alignments {
            total_score += aln.score;
            combined_cigar.extend_from_slice(&aln.cigar);
        }

        Alignment {
            score: total_score,
            cigar: combined_cigar,
        }
    }
}

/// Convenience function for quick alignment with default parameters
pub fn align(query: &[u8], reference: &[u8]) -> Option<Alignment> {
    match align_inner(query, reference) {
        Ok(aln) => Some(aln),
        Err(_) => None,
    }
}

const USE_WFA: bool = false;

/// Convenience function for quick alignment with default parameters
pub fn align_inner(
    query: &[u8],
    reference: &[u8],
) -> std::result::Result<Alignment, AlignmentError> {
    metrics::histogram!("align_ref_len").record(reference.len() as f64);
    metrics::histogram!("align_query_len").record(query.len() as f64);

    // First try the WFA aligner
    let start = std::time::Instant::now();
    let result = if USE_WFA {
        WfAligner::new(AlignParams::default())
            .align(query, reference)
            .map_err(|error| AlignmentError::WFError(error))
    } else {
        block::align(query, reference).map_err(|error| AlignmentError::BlockError(error))
    };
    let elapsed = start.elapsed();
    metrics::histogram!("wf_align_time_us").record(elapsed.as_micros() as f64);
    let mut alignment = match result {
        Ok(aln) => Ok(aln),
        Err(error) => {
            metrics::histogram!("align_fail_ref").record(reference.len() as f64);
            metrics::histogram!("align_fail_query").record(query.len() as f64);
            metrics::histogram!("align_fail_time").record(elapsed.as_secs_f64());

            log::debug!(
                "WFA alignment failed (score too high?): {}. Falling back to MiniAlign.",
                error
            );

            let start = std::time::Instant::now();
            let result = mini::align(query, reference, 15);
            let elapsed = start.elapsed();
            metrics::histogram!("mini_align_time_us").record(elapsed.as_micros() as f64);
            match result {
                Ok(aln) => Ok(aln),
                Err(error) => {
                    metrics::histogram!("mini_fail_ref").record(reference.len() as f64);
                    metrics::histogram!("mini_fail_query").record(query.len() as f64);
                    metrics::histogram!("mini_fail_time").record(elapsed.as_secs_f64());
                    Err(AlignmentError::from(error))
                }
            }
        }
    };

    // Ensure the CIGAR is normalized (adjacent same-type ops merged)
    if let Ok(ref mut aln) = alignment {
        aln.normalize();
    }

    alignment
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical() {
        let result = align(b"ACGT", b"ACGT").unwrap();
        assert_eq!(result.score, 0);
        assert_eq!(result.cigar_string(), "4=");
    }

    #[test]
    fn test_single_mismatch() {
        let expected = if USE_WFA { 4 } else { 2 };
        let result = align(b"ACGT", b"ACTT").unwrap();
        assert_eq!(result.score, expected); // mismatch penalty
        assert_eq!(result.cigar_string(), "2=1X1=");
    }

    #[test]
    fn test_single_insertion() {
        let expected = if USE_WFA { 4 } else { 7 };
        let result = align(b"ACGT", b"ACT").unwrap();
        // query has extra G
        assert_eq!(result.score, expected);
        assert!(result.cigar_string().contains('I'));
    }

    #[test]
    fn test_single_deletion() {
        let expected = if USE_WFA { 4 } else { 6 };
        let result = align(b"ACT", b"ACGT").unwrap();
        // query missing G
        assert_eq!(result.score, expected);
        assert!(result.cigar_string().contains('D'));
    }

    #[test]
    fn test_empty() {
        let result = align(b"", b"").unwrap();
        assert_eq!(result.score, 0);
        assert!(result.cigar.is_empty());
    }

    #[test]
    fn test_query_empty() {
        let result = align(b"", b"ACGT").unwrap();
        assert_eq!(result.cigar_string(), "4D");
    }

    #[test]
    fn test_reference_empty() {
        let result = align(b"ACGT", b"").unwrap();
        assert_eq!(result.cigar_string(), "4I");
    }

    #[test]
    fn test_longer_sequences() {
        let query = b"ACGTACGTACGT";
        let reference = b"ACGTACGTACGT";
        let result = align(query, reference).unwrap();
        assert_eq!(result.score, 0);
        assert_eq!(result.cigar_string(), "12=");
    }

    #[test]
    fn test_with_gaps() {
        let query = b"ACGTACGT";
        let reference = b"ACGTTTTACGT";
        let result = align(query, reference).unwrap();
        assert!(result.score > 0);
        // Should have a deletion in the middle
        println!("CIGAR: {}", result.cigar_string());
    }

    #[test]
    fn test_from_data_1() {
        let reference = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGAATAGTTAGTGCTTCTTGTTGAATTGATCCCTTTACCATTATGTAAATGGCTTTCTTTGTCTCTTTTGATCTTTTTGGTTTAAAATCTGTTTTATCAGAGACTATGACAGCAATCCCTGCTTTTTTTGCTTTTCATTTGCTTGGTAGATCTTCCTCTGTCCCTTTATTTTGAGCCCATGTGTGTGTCTGCACATGAGATGGGTCTCCTGAATACAGCACGTTGATGGGTCTCGAATCTTTGTCCAGTTTGTCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTGTTTTGCTCATTAGTTGATGCCGTTTCTTCCTAGCATCAACGGTCTTTACAATTTGGCCTGTTTTTGCAGTGGCTGGTACCAGTTGTTCCTTTCCATGTTTAATGCTTCCTTCAGGAGCTCTTGTAAGGCAGGCCTGGTGGTGACAAAATCTCTCAGGATTTGCTTGTCTGCAAAGGATTTTGTTTCTCCTTCACTTATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGTAAATTATTTTCTTTGAGAATGTTGAATATTGGCCCCCACTCTCATCTGGCTTGTAGGGTTTCTGCCGAGAGATCTGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTATCCCAGCCTTTCTCTCTGGCTGACCTTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTGGGGGTTGCTCTTCTCAAGGAGTGTCTTTGTAGTGTTCTCTGTATTTCTTGAATTTGAATGTTGGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGAATAATATCCTGAAGAGTGTTTTCCAGCTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGATCTTCTCACATAGTTCCATACTTCTTGGAGGCTTTGTTTGTTTCATTTTACTCTTTTTCTCTAAACTTCTCTTCTTGCTTCATTTCATTAATTTGATCTTCAGTCACTGAAACCCTTTCTTCCATTGATCGAATCAGCTACTGAAGCTTCTGTGTGTGTCACGTAGTTCTCGTGTCATGGTTTTCAGCTCCTTCAGGTCATTTAAGGTTTTCTCTACACTGGTTACTCTAGTTAGCCTTTTGTCTAATCTTTTTTCAAGGTTTTTAGCTTCCTTGCGGTGGGTTTGAATATCCTCCTTTAGCTCAGAGAAATTTGTTTTTACCGACCTTCTGAAGCCTAATTCTGTCAACTCGTCAAAGTCATTCTCCATCCAGCCTTGTTCTGTTGCTGGCGAGGAGCTGTGATCCTTTGGAGGAGAAGAGGCACTCTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTTGTTTTTATCTACCTTTGGTCTTTGATGATGGTGACCTACAGATGGGGTTATTGGTGTGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTAACAGTCAGGTCCCTCAGCTGCAGGTCTGTTGGAGTTTGCTAGAGGTCCACTCCAGACACTGTTTGCCTGGGTATCACCTTTGGAGGCTGCAGAACAGCAAATATTGCAGAACAGCAAATATTGCTGCCTGATCCTTCCTCTGGAAGCATCGTCCCAGAGGGGCATACGGCAGCATGAGATGTCAGTCAGCCCCCACTGGGAGGTGTCTCCCTGTTAGGCTACACGGGGGTCAGGGACCCACTTGAGGAGGCAGTCTGTCCGTTCTCAGAGCTCAAACACCGTGCTGGGAGAACCACTGCTCTCTTCAGAGCAGTGCAGACAAGGACATTTAAGTCTGCAGAAGTTTCTGCTGCCTTTTGTTCAGCTATGCCCTGCCACCAGAGGTGGAGTCTATAGAGGCAGCAAGCCTTGTGGTGCTGTGGTGGGCTCTGCCAAGTTCGAGCTTCCTGGCTGCTTTGTTTACCTACTTAAGCCTCAGCAATGGTGGACGCCCCTTCCCCAGCCAGGCTGC";
        let query = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGGGTAGTTAGCTCCCTTTACCATTATGTAATGGCCTTCTTTGTCCCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACAACACCTGCTTTTTTTGTTTTCCATTTGCTTGGTAGGTCTTCCTCCATCCCTTTATTTTGAGCCTATGTGTGTGTCTGCACATGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTCATCTAACTTGCCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATTCTGTCATTATGATGTTAGCTGGTTATTTTTCCCGTTAGTTGATGCAGTTTCTTCCTAGCATCGATGGTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGTGCTTCCTTCAGGAGCTCTTGTAAGGCAGGGCTGGTGGTGACAAAATCTCTCAGCATTTGCTGGTCTATAAAGGATTTTATTTCTCCTTATGAAGCTTTGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTAAGAATGTTGAATATTGGTGCCCACTCTCTTCTGACTTGTAGAGTTTCTGTTGAGAGATCCACTGTTAGTCTGATGGGCTTCCCTTTGTGGCTAACTCGACCTTTCTCTCTGGGTGCCATTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTTGGGGTTGCTCCTCTCGAGGAGCATCTTGGTAGTGTTCTCTGTATTTCCTGAGTTTGAATGTTTGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGGACAATATCCTGAAGAGTGTTTTCGAACTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGGTGTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGTTCATTCTTTTTACTCTTTTTTCTCTAAACTTCTCACTTCATTAATTTGATCTTCAATCACTGATACCCTTTCTTTCAGTTTATTGAATCAACTACTGAAGCTTGTGCATGTGTCACATAGTTCTTGTTCCATGGTTTTCAGCTCCATCAGGTCATTTAAGGTCTCCACACTGCTTATTCTAGTTAGCCATTCATCTAATCTGTTTGCAAGGCTTTTAGCTTCCTTGTGATGGGTTCGAATACCTCCCTTAACTCAGAGAAGTTTGTTATTACCAACCTTCTGAAGCCTACTTCTGTCAGCTCATCAAAGTCATTCTCCGTCCAGCTTTATTCCGTTGCTGGCAAGGAGCTGTAATCCTTTGCAGGAGAAGGGATGCTGTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTGGTTTTATCTACCTTTGGTCTTCGATGATGGTGACCCACAGATGGGGTTTTGGTGTGGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTGACAGTCAGGTCCCTCAGCTGCAGATCTGTTGGAGTTTGCTGGAGGTCCACTCCAGACTCTGTTTACCTGTGTATCACCAGCAGAGGCTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAACAGCAAATATTACTGCCTGATCCTTACTCTGGAAGCTTCATTTCAGAGGGGCACCCAGCTCTATGAGGTGTCATTCGGCCCCTACTGGGAGATGTCTCCCAGTTAGGCTACACAGGGGTCAGGGACACACTTGAGGAGGCAGTCTGTCCATTCTCAGAGCTCAAACTCCATGCTAGGAGAACCACTGCTCTCTTCAGAGCTGTCAGATAGGGACATTTAAGTCTGCAGAAGTTTCTGCTGCCTTTTGTTCAGCTATGCCCTGCCCCCAGAGGTGGAGTCTACAGAGTCAGGCAGGCCTCCTTGAGCTGTGGTGGGCTCCACCCAGTTCGAGCTTCCCAGCCGCTTTGTTTACCTACTCAAGCTTCAGCAATGGCGGACGCCCCTTCCCCAGCCAGGCTGC";
        let alignment = align(query, reference).unwrap();

        println!(
            "Reference len: {}, Query len: {}",
            reference.len(),
            query.len()
        );
        println!("Actual CIGAR:   {}", alignment.cigar_string());
        if let Err(e) = alignment.validate(reference, query, 0) {
            println!("Validation error: {}", e);
        }
        alignment
            .validate(reference, query, 0)
            .expect("Alignment validation failed");
    }

    #[test]
    fn test_from_data_1_short() {
        // Truncated version of test_from_data_1 for debugging
        // Cut at position ~1550 in reference (just before error at ref_pos 1648)
        let reference = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGAATAGTTAGTGCTTCTTGTTGAATTGATCCCTTTACCATTATGTAAATGGCTTTCTTTGTCTCTTTTGATCTTTTTGGTTTAAAATCTGTTTTATCAGAGACTATGACAGCAATCCCTGCTTTTTTTGCTTTTCATTTGCTTGGTAGATCTTCCTCTGTCCCTTTATTTTGAGCCCATGTGTGTGTCTGCACATGAGATGGGTCTCCTGAATACAGCACGTTGATGGGTCTCGAATCTTTGTCCAGTTTGTCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTGTTTTGCTCATTAGTTGATGCCGTTTCTTCCTAGCATCAACGGTCTTTACAATTTGGCCTGTTTTTGCAGTGGCTGGTACCAGTTGTTCCTTTCCATGTTTAATGCTTCCTTCAGGAGCTCTTGTAAGGCAGGCCTGGTGGTGACAAAATCTCTCAGGATTTGCTTGTCTGCAAAGGATTTTGTTTCTCCTTCACTTATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGTAAATTATTTTCTTTGAGAATGTTGAATATTGGCCCCCACTCTCATCTGGCTTGTAGGGTTTCTGCCGAGAGATCTGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTATCCCAGCCTTTCTCTCTGGCTGACCTTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTGGGGGTTGCTCTTCTCAAGGAGTGTCTTTGTAGTGTTCTCTGTATTTCTTGAATTTGAATGTTGGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGAATAATATCCTGAAGAGTGTTTTCCAGCTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGATCTTCTCACATAGTTCCATACTTCTTGGAGGCTTTGTTTGTTTCATTTTACTCTTTTTCTCTAAACTTCTCTTCTTGCTTCATTTCATTAATTTGATCTTCAGTCACTGAAACCCTTTCTTCCATTGATCGAATCAGCTACTGAAGCTTCTGTGTGTGTCACGTAGTTCTCGTGTCATGGTTTTCAGCTCCTTCAGGTCATTTAAGGTTTTCTCTACACTGGTTACTCTAGTTAGCCTTTTGTCTAATCTTTTTTCAAGGTTTTTAGCTTCCTTGCGGTGGGTTTGAATATCCTCCTTTAGCTCAGAGAAATTTGTTTTTACCGACCTTCTGAAGCCTAATTCTGTCAACTCGTCAAAGTCATTCTCCATCCAGCCTTGTTCTGTTGCTGGCGAGGAGCTGTGATCCTTTGGAGGAGAAGAGGCACTCTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTTGTTTTTATCTACCTTTGGTCTTTGATGATGGTGACCTACAGATGGGGTTATTGGTGTGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTAACAGTCAGGTCCCTCAGCTGCAGGTCTGTTGGAGTTTGCTAGAGGTCCACTCCAGACACTGTTTGCCTGGGTATCACCTTTGGAGGCTGCAGAACAGCAAATATTGCAGAACAGCAAATATTGC";
        let query = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGGGTAGTTAGCTCCCTTTACCATTATGTAATGGCCTTCTTTGTCCCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACAACACCTGCTTTTTTTGTTTTCCATTTGCTTGGTAGGTCTTCCTCCATCCCTTTATTTTGAGCCTATGTGTGTGTCTGCACATGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTCATCTAACTTGCCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATTCTGTCATTATGATGTTAGCTGGTTATTTTTCCCGTTAGTTGATGCAGTTTCTTCCTAGCATCGATGGTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGTGCTTCCTTCAGGAGCTCTTGTAAGGCAGGGCTGGTGGTGACAAAATCTCTCAGCATTTGCTGGTCTATAAAGGATTTTATTTCTCCTTATGAAGCTTTGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTAAGAATGTTGAATATTGGTGCCCACTCTCTTCTGACTTGTAGAGTTTCTGTTGAGAGATCCACTGTTAGTCTGATGGGCTTCCCTTTGTGGCTAACTCGACCTTTCTCTCTGGGTGCCATTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTTGGGGTTGCTCCTCTCGAGGAGCATCTTGGTAGTGTTCTCTGTATTTCCTGAGTTTGAATGTTTGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGGACAATATCCTGAAGAGTGTTTTCGAACTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGGTGTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGTTCATTCTTTTTACTCTTTTTTCTCTAAACTTCTCACTTCATTAATTTGATCTTCAATCACTGATACCCTTTCTTTCAGTTTATTGAATCAACTACTGAAGCTTGTGCATGTGTCACATAGTTCTTGTTCCATGGTTTTCAGCTCCATCAGGTCATTTAAGGTCTCCACACTGCTTATTCTAGTTAGCCATTCATCTAATCTGTTTGCAAGGCTTTTAGCTTCCTTGTGATGGGTTCGAATACCTCCCTTAACTCAGAGAAGTTTGTTATTACCAACCTTCTGAAGCCTACTTCTGTCAGCTCATCAAAGTCATTCTCCGTCCAGCTTTATTCCGTTGCTGGCAAGGAGCTGTAATCCTTTGCAGGAGAAGGGATGCTGTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTGGTTTTATCTACCTTTGGTCTTCGATGATGGTGACCCACAGATGGGGTTTTGGTGTGGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTGACAGTCAGGTCCCTCAGCTGCAGATCTGTTGGAGTTTGCTGGAGGTCCACTCCAGACTCTGTTTACCTGTGTATCACCAGCAGAGGCTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAACAGCAAATATTAC";
        let alignment = align(query, reference).unwrap();

        println!(
            "Reference len: {}, Query len: {}",
            reference.len(),
            query.len()
        );
        println!("Actual CIGAR:   {}", alignment.cigar_string());
        if let Err(e) = alignment.validate(reference, query, 0) {
            println!("Validation error: {}", e);
        }
        alignment
            .validate(reference, query, 0)
            .expect("Alignment validation failed (short)");
    }

    #[test]
    fn test_gap_before_seed_35() {
        // This is the gap that causes the error:
        // query[1610..1682] (len 72), ref[1648..1686] (len 38)
        // From the full test_from_data_1_short sequences

        let full_query = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGGGTAGTTAGCTCCCTTTACCATTATGTAATGGCCTTCTTTGTCCCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACAACACCTGCTTTTTTTGTTTTCCATTTGCTTGGTAGGTCTTCCTCCATCCCTTTATTTTGAGCCTATGTGTGTGTCTGCACATGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTCATCTAACTTGCCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATTCTGTCATTATGATGTTAGCTGGTTATTTTTCCCGTTAGTTGATGCAGTTTCTTCCTAGCATCGATGGTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGTGCTTCCTTCAGGAGCTCTTGTAAGGCAGGGCTGGTGGTGACAAAATCTCTCAGCATTTGCTGGTCTATAAAGGATTTTATTTCTCCTTATGAAGCTTTGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTAAGAATGTTGAATATTGGTGCCCACTCTCTTCTGACTTGTAGAGTTTCTGTTGAGAGATCCACTGTTAGTCTGATGGGCTTCCCTTTGTGGCTAACTCGACCTTTCTCTCTGGGTGCCATTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTTGGGGTTGCTCCTCTCGAGGAGCATCTTGGTAGTGTTCTCTGTATTTCCTGAGTTTGAATGTTTGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGGACAATATCCTGAAGAGTGTTTTCGAACTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGGTGTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGTTCATTCTTTTTACTCTTTTTTCTCTAAACTTCTCACTTCATTAATTTGATCTTCAATCACTGATACCCTTTCTTTCAGTTTATTGAATCAACTACTGAAGCTTGTGCATGTGTCACATAGTTCTTGTTCCATGGTTTTCAGCTCCATCAGGTCATTTAAGGTCTCCACACTGCTTATTCTAGTTAGCCATTCATCTAATCTGTTTGCAAGGCTTTTAGCTTCCTTGTGATGGGTTCGAATACCTCCCTTAACTCAGAGAAGTTTGTTATTACCAACCTTCTGAAGCCTACTTCTGTCAGCTCATCAAAGTCATTCTCCGTCCAGCTTTATTCCGTTGCTGGCAAGGAGCTGTAATCCTTTGCAGGAGAAGGGATGCTGTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTGGTTTTATCTACCTTTGGTCTTCGATGATGGTGACCCACAGATGGGGTTTTGGTGTGGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTGACAGTCAGGTCCCTCAGCTGCAGATCTGTTGGAGTTTGCTGGAGGTCCACTCCAGACTCTGTTTACCTGTGTATCACCAGCAGAGGCTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAACAGCAAATATTAC";
        let full_ref = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGAATAGTTAGTGCTTCTTGTTGAATTGATCCCTTTACCATTATGTAAATGGCTTTCTTTGTCTCTTTTGATCTTTTTGGTTTAAAATCTGTTTTATCAGAGACTATGACAGCAATCCCTGCTTTTTTTGCTTTTCATTTGCTTGGTAGATCTTCCTCTGTCCCTTTATTTTGAGCCCATGTGTGTGTCTGCACATGAGATGGGTCTCCTGAATACAGCACGTTGATGGGTCTCGAATCTTTGTCCAGTTTGTCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTGTTTTGCTCATTAGTTGATGCCGTTTCTTCCTAGCATCAACGGTCTTTACAATTTGGCCTGTTTTTGCAGTGGCTGGTACCAGTTGTTCCTTTCCATGTTTAATGCTTCCTTCAGGAGCTCTTGTAAGGCAGGCCTGGTGGTGACAAAATCTCTCAGGATTTGCTTGTCTGCAAAGGATTTTGTTTCTCCTTCACTTATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGTAAATTATTTTCTTTGAGAATGTTGAATATTGGCCCCCACTCTCATCTGGCTTGTAGGGTTTCTGCCGAGAGATCTGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTATCCCAGCCTTTCTCTCTGGCTGACCTTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTGGGGGTTGCTCTTCTCAAGGAGTGTCTTTGTAGTGTTCTCTGTATTTCTTGAATTTGAATGTTGGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGAATAATATCCTGAAGAGTGTTTTCCAGCTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGATCTTCTCACATAGTTCCATACTTCTTGGAGGCTTTGTTTGTTTCATTTTACTCTTTTTCTCTAAACTTCTCTTCTTGCTTCATTTCATTAATTTGATCTTCAGTCACTGAAACCCTTTCTTCCATTGATCGAATCAGCTACTGAAGCTTCTGTGTGTGTCACGTAGTTCTCGTGTCATGGTTTTCAGCTCCTTCAGGTCATTTAAGGTTTTCTCTACACTGGTTACTCTAGTTAGCCTTTTGTCTAATCTTTTTTCAAGGTTTTTAGCTTCCTTGCGGTGGGTTTGAATATCCTCCTTTAGCTCAGAGAAATTTGTTTTTACCGACCTTCTGAAGCCTAATTCTGTCAACTCGTCAAAGTCATTCTCCATCCAGCCTTGTTCTGTTGCTGGCGAGGAGCTGTGATCCTTTGGAGGAGAAGAGGCACTCTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTTGTTTTTATCTACCTTTGGTCTTTGATGATGGTGACCTACAGATGGGGTTATTGGTGTGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTAACAGTCAGGTCCCTCAGCTGCAGGTCTGTTGGAGTTTGCTAGAGGTCCACTCCAGACACTGTTTGCCTGGGTATCACCTTTGGAGGCTGCAGAACAGCAAATATTGCAGAACAGCAAATATTGC";

        // Extract the gap region
        let gap_query = &full_query[1610..1682];
        let gap_ref = &full_ref[1648..1686];

        println!(
            "Gap query (len {}): {}",
            gap_query.len(),
            String::from_utf8_lossy(gap_query)
        );
        println!(
            "Gap ref (len {}): {}",
            gap_ref.len(),
            String::from_utf8_lossy(gap_ref)
        );

        // Align with WFA directly
        let aln = WfAligner::new(AlignParams::default())
            .align(gap_query, gap_ref)
            .unwrap();
        println!("WFA CIGAR: {}", aln.cigar_string());

        // Validate
        if let Err(e) = aln.validate(gap_ref, gap_query, 0) {
            println!("Validation error: {}", e);
        }
        aln.validate(gap_ref, gap_query, 0)
            .expect("Gap alignment validation failed");
    }

    #[test]
    #[ignore]
    fn test_from_data_2() {
        let reference = b"TTAGGCAGTGCCCCAGTGGGGACGCTGTGTGGGGTCTCCAATCCCACATTTCCCTTCTGCACTGCTCTAGCAGAGGTTCTCCATGAGGGCTCTTACCCTGCAGCAAACTTCTGCCTGGGCATTCAGGCATTTCTGTACAACCTCTGAAATCTAGGTGGAAGTTCCCAAACCTCAATTCTTGACTTCTGTGCACCCACAGGCTCAACACCACATAGGAGCTGCCAAGGCTTGGGGCTTGCACCCTCTGAAGCCACAGCCTGAGCTGTACTTTGGCTCCTTGTAGCCATGGCTAGAGTGGCTGGGACACAGGACACCAAATCCCTAGGCTGCACACAGCAGGTGGGCCCTAGGCCCTGCCCACAAAACAATTTTTTCCTCCTAGGACTCTGGGCCTGTGATGGGAAGGGCTGCTGTGAAGAAGACCTCTGACATGCCCTGGAGACATTTTCCCTATTGTCTTGGCAATTAACATTTAGCTCCTCATTAATCATGCAAATTTTTGCAGCCAGCTTGAATTTCTCCTCAGAAAATGGGTTTTTCTTTCCTATTACATTGTTAGGCTGCAAAATTTCCAAACTTTTATGCTCTGTTTCCCTTTTAAACTGAATGCTTTTTAACAGCACTCAAGTCACCTCTTGAATGCTTTGCTGCTTAGAAATTTCTTCTGCCAGATGCCCTAAATCATCTCCCTCAAGTTCAAAGTTCCACAGATCTTTAGGGTGGGAGCAAAATGTCACCAGTCTGTTTGCTAAAACATAGCAAGAGTCACCATTTCTCCCC";
        let query = b"CTAGGCAAGGGGATTCTCTGTGGGGGCTCACTCCCCATATTTCCCTTCCACATGGCCAGTAGAGGTTCTCCATGAGGGCTCTGCCCCTGCAGCAAACTTCTGCTTGGACATGCAGGCATTTCCATATGTCCTCTGAAATCTAGGCGGAGGTTCCCAAACCTCAATTCTTGACTCTGTGCACTGACAGGCTCAGCATCACATGAAAATCACAAGGCTTGGGGAGGGCTTATCCTTCAAGCAATGACCTCTTTGAGCCAGCTGGAGCTGAAGCAGCGGGAGTGAGGCACCATGTCCTGAGGCGGCACAGAGCAGGGCAGCCCTGGGCTCAGCCCAGGAAACCATTTTTCCCTACTAGGTGTCTGTGCCTGTGATGGGAGGGGTGGCCATGAAGACCTCTACATGGCCTGGAGACATTTTCTCCATTGCCTTAATGATTAACATTTGGCACATTTTTCAGATACATGTGGCTGTAAATATGTATCAATATCAGTGATATTTGCAGCTGGCTTGAATTTCTCTTCAAACAATGGGTTTTTCTTTAGTATTGCATCATCAGTACTGTAATTTTAAGTTTTTTTGCTTGTGCTTCCTCTTTTCATGCTTTCTTGAGAAATTTCTTCTGCCAGATACCCTAAATCATATCTGTCTCAAGTTCAACGTTCCACAGTCTCCAGGGCAGGGCAAGTCACCAGCACCAGTCTCTTCTGCCAAAGCATGCAAGAGTCACCTTTGCTCCAG";
        let result = align_inner(query, reference);
        if let Err(e) = &result {
            println!("Validation error: {}", e);
        }
        let alignment = result.unwrap();
        alignment
            .validate(reference, query, 0)
            .expect("Alignment validation failed");
    }

    #[test]
    fn test_cigar_scoring() {
        let cig = "2D1=1X1=1X1=1X2=2X2=1X1=1X2=1X2=1X1=2X2=12X1=1X1=3X1=2X1=1X1=2X5=1X1=2X1=1X2=1X3=21D3=1X6=2X1=1X1=2X1=3X2=1X1=2X3=1X1=3X1=1D2=1X2=1X3=2X2=3X3=1X18D176=1X82=1D";
        let cigar = CigarOp::make(cig).expect("bad cigar string");
        let params = AlignParams::default();
        let mut total = 0.0;
        for op in &cigar {
            total += op.score(&params) as f64;
            println!("{:?} -> {:.2} (running total: {:.2})", op, op.score(&params) as f64, total);
        }
        let alignment = Alignment { score: 0, cigar };
        assert_eq!(alignment.score(&params), 290.0)
    }
}
