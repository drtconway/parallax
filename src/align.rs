pub use noodles::sam::alignment::record::cigar::op::Kind;
pub use noodles::sam::alignment::record::cigar::Op;

use crate::config;
use crate::scores::{DivergenceScore, QualityScore};

pub mod block;

#[cfg(feature = "attic")]
pub mod anchor;
#[cfg(feature = "attic")]
pub mod kmer_anchors;
#[cfg(feature = "attic")]
pub mod lcs_anchors;
#[cfg(feature = "attic")]
pub mod mini;
#[cfg(feature = "attic")]
pub mod wfa;

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
    /// Divergence is calculated as the sum of mismatch and gap penalties (lower is better).
    #[allow(dead_code)]
    pub fn divergence(&self, op: Op) -> DivergenceScore {
        let n = op.len() as f64;
        (match op.kind() {
            Kind::SequenceMatch => 0.0,
            Kind::SequenceMismatch => self.mismatch as f64 * n,
            Kind::Insertion | Kind::Deletion => self.gap_open as f64 + self.gap_extend as f64 * n,
            Kind::SoftClip | Kind::HardClip => 0.0,
            _ => 0.0,
        })
        .into()
    }

    /// Score an individual alignment operation
    pub fn quality(&self, op: Op) -> QualityScore {
        let n = op.len() as i32;
        (match op.kind() {
            Kind::SequenceMatch => self.match_score * n,
            Kind::SequenceMismatch => -self.mismatch * n,
            Kind::Insertion | Kind::Deletion => -self.gap_open - self.gap_extend * n,
            Kind::SoftClip | Kind::HardClip => 0,
            _ => 0,
        } as f64)
            .into()
    }
}

/// High-level aligner abstraction that wraps the underlying alignment engine.
///
/// Provides a stable API for global alignment, left extension, and right extension,
/// allowing the underlying aligner implementation to be swapped without changing callers.
/// Configuration is read from the global config infrastructure.
pub struct Aligner {
    inner: block::BlockAligner,
}

impl Aligner {
    /// Create a new Aligner using the global config infrastructure.
    pub fn new() -> Self {
        let cfg = config::get();
        Self {
            inner: block::BlockAligner::new(&cfg.block_aligner),
        }
    }

    /// Create an Aligner with explicit configuration.
    #[allow(dead_code)]
    pub fn with_config(config: &crate::config::BlockAlignerConfig) -> Self {
        Self {
            inner: block::BlockAligner::new(config),
        }
    }

    /// Create an Aligner with default configuration (no global config required).
    ///
    /// Useful for tests where the global config may not be initialized.
    pub fn with_defaults() -> Self {
        Self {
            inner: block::BlockAligner::with_defaults(),
        }
    }

    /// Align two sequences using global alignment with fallback to mini-aligner.
    ///
    /// Returns `None` if all alignment strategies fail.
    pub fn align(&mut self, query: &[u8], reference: &[u8]) -> Option<Alignment> {
        match self.align_inner(query, reference) {
            Ok(aln) => Some(aln),
            Err(_) => None,
        }
    }

    /// Align two sequences, returning a detailed error on failure.
    pub fn align_inner(
        &mut self,
        query: &[u8],
        reference: &[u8],
    ) -> std::result::Result<Alignment, AlignmentError> {
        metrics::histogram!("align_ref_len").record(reference.len() as f64);
        metrics::histogram!("align_query_len").record(query.len() as f64);

        let start = std::time::Instant::now();
        let result = self
            .inner
            .align(query, reference)
            .map_err(AlignmentError::BlockError);
        let elapsed = start.elapsed();
        metrics::histogram!("align_time_us").record(elapsed.as_micros() as f64);

        let mut alignment = result?;

        alignment.normalize();

        Ok(alignment)
    }

    /// Extend alignment rightward (forward) with X-drop early termination.
    pub fn extend_right(
        &mut self,
        query: &[u8],
        reference: &[u8],
    ) -> std::result::Result<Alignment, AlignmentError> {
        self.inner
            .extend_right(query, reference)
            .map_err(AlignmentError::BlockError)
    }

    /// Extend alignment leftward (backward) with X-drop early termination.
    pub fn extend_left(
        &mut self,
        query: &[u8],
        reference: &[u8],
    ) -> std::result::Result<Alignment, AlignmentError> {
        self.inner
            .extend_left(query, reference)
            .map_err(AlignmentError::BlockError)
    }
}

impl Default for Aligner {
    fn default() -> Self {
        Self::with_defaults()
    }
}

#[derive(Debug)]
pub enum AlignmentError {
    BlockError(block::BlockAlignerError),
    #[cfg(feature = "attic")]
    WFError(wfa::WfaFailure),
    #[cfg(feature = "attic")]
    MiniError(mini::MiniAlignError),
}

impl std::fmt::Display for AlignmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentError::BlockError(e) => write!(f, "Block aligner error: {}", e),
            #[cfg(feature = "attic")]
            AlignmentError::WFError(e) => write!(f, "WFA alignment error: {}", e),
            #[cfg(feature = "attic")]
            AlignmentError::MiniError(e) => write!(f, "Mini alignment error: {}", e),
        }
    }
}

impl std::error::Error for AlignmentError {}

#[cfg(feature = "attic")]
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

#[cfg(feature = "attic")]
impl From<mini::MiniAlignError> for AlignmentError {
    fn from(err: mini::MiniAlignError) -> Self {
        AlignmentError::MiniError(err)
    }
}

/// Character code for a CIGAR operation kind.
fn op_char(kind: Kind) -> char {
    match kind {
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Match => 'M',
        Kind::Skip => 'N',
        Kind::Pad => 'P',
    }
}

/// Alignment result
#[derive(Clone, Debug)]
pub struct Alignment {
    /// Edit distance (lower is better)
    #[allow(dead_code)]
    pub divergence: DivergenceScore,
    /// CIGAR operations
    pub cigar: Vec<Op>,
}

impl Alignment {
    /// Create a new perfect match alignment
    #[allow(dead_code)]
    pub fn from_perfect_match(length: usize) -> Self {
        Self {
            divergence: DivergenceScore::ZERO,
            cigar: vec![Op::new(Kind::SequenceMatch, length)],
        }
    }

    /// Format CIGAR as string
    pub fn cigar_string(&self) -> String {
        self.cigar
            .iter()
            .map(|op| format!("{}{}", op.len(), op_char(op.kind())))
            .collect()
    }

    /// Format CIGAR in the basic format merging = and X into M
    /// (e.g., 10M1I5M2D3M)
    #[allow(dead_code)]
    pub fn basic_cigar_string(&self) -> String {
        let mut merged: Vec<Op> = Vec::new();
        for &op in &self.cigar {
            let n = op.len();
            match op.kind() {
                Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if let Some(last) = merged.last_mut() {
                        match last.kind() {
                            Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match => {
                                *last = Op::new(Kind::Match, last.len() + n);
                            }
                            _ => merged.push(Op::new(Kind::Match, n)),
                        }
                    } else {
                        merged.push(Op::new(Kind::Match, n));
                    }
                }
                _ => merged.push(op),
            }
        }
        merged
            .iter()
            .map(|op| format!("{}{}", op.len(), op_char(op.kind())))
            .collect()
    }

    /// Compute the query (read) length consumed by this CIGAR.
    /// This is the sum of M, I, S, =, X operations.
    /// For valid SAM, this must equal the length of the SEQ field.
    pub fn query_length(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| op.kind().consumes_read())
            .map(|op| op.len())
            .sum()
    }

    /// Compute the reference span consumed by this CIGAR.
    /// This is the sum of M, D, N, =, X operations.
    #[allow(dead_code)]
    pub fn reference_span(&self) -> u64 {
        self.cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len() as u64)
            .sum()
    }

    /// Compute the query (read) bases consumed by alignment operations (excluding soft clips).
    /// This is the sum of M, I, =, X operations only.
    #[allow(dead_code)]
    pub fn query_consumed(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match | Kind::Insertion
                )
            })
            .map(|op| op.len())
            .sum()
    }

    /// Compute the reference bases consumed by alignment operations.
    /// This is the sum of M, D, =, X operations.
    pub fn reference_consumed(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum()
    }

    /// Return the size of the leading hard clip, if any.
    pub fn leading_hard_clip(&self) -> usize {
        match self.cigar.first() {
            Some(op) if op.kind() == Kind::HardClip => op.len(),
            _ => 0,
        }
    }

    pub fn total_identity(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| op.kind() == Kind::SequenceMatch)
            .map(|op| op.len())
            .sum()
    }

    /// Merge adjacent operations of same type
    pub fn normalize(&mut self) {
        if self.cigar.is_empty() {
            return;
        }
        let mut merged: Vec<Op> = Vec::with_capacity(self.cigar.len());
        for op in self.cigar.drain(..) {
            if let Some(last) = merged.last_mut() {
                if last.kind() == op.kind() {
                    *last = Op::new(op.kind(), last.len() + op.len());
                } else {
                    merged.push(op);
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

        for &op in &self.cigar {
            let n = op.len();
            match op.kind() {
                Kind::SequenceMatch => {
                    for _ in 0..n {
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        ref_line.push(r as char);
                        query_line.push(q as char);
                        match_line.push('|');
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Kind::SequenceMismatch => {
                    for _ in 0..n {
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        ref_line.push(r as char);
                        query_line.push(q as char);
                        match_line.push(' ');
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Kind::Insertion => {
                    for _ in 0..n {
                        let q = query[query_pos];
                        ref_line.push('-');
                        query_line.push(q as char);
                        match_line.push(' ');
                        query_pos += 1;
                    }
                }
                Kind::Deletion => {
                    for _ in 0..n {
                        let r = reference[ref_pos];
                        ref_line.push(r as char);
                        query_line.push('-');
                        match_line.push(' ');
                        ref_pos += 1;
                    }
                }
                Kind::SoftClip => {
                    query_pos += n;
                }
                Kind::HardClip => {}
                _ => {}
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
        let len = self.cigar.len();
        for i in 0..len {
            let op = self.cigar[i];
            if op.len() == 0 && op.kind() != Kind::HardClip {
                return Err(format!(
                    "CIGAR op {} has zero length: {:?}",
                    i, op
                ));
            }
            if i > 0 && self.cigar[i].kind() == self.cigar[i - 1].kind() {
                return Err(format!(
                    "CIGAR ops {} and {} are adjacent and of same type: {:?}, {:?}",
                    i - 1, i, self.cigar[i - 1], self.cigar[i]
                ));
            }
            if i > 0 && i < len - 1 {
                if op.kind() == Kind::SoftClip {
                    return Err(format!(
                        "CIGAR op {} is a soft clip in the middle of the alignment: {:?}",
                        i, op
                    ));
                }
                if op.kind() == Kind::HardClip {
                    return Err(format!(
                        "CIGAR op {} is a hard clip in the middle of the alignment: {:?}",
                        i, op
                    ));
                }
            }
        }

        // Check each CIGAR operation against the sequences
        for (op_idx, &op) in self.cigar.iter().enumerate() {
            let n = op.len();
            let ch = op_char(op.kind());
            match op.kind() {
                Kind::SequenceMatch => {
                    for i in 0..n {
                        if ref_pos >= reference.len() {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): ref_pos {} exceeds reference length {} at offset {}",
                                op_idx, n, ref_pos, reference.len(), i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): query_pos {} exceeds query length {} at offset {}",
                                op_idx, n, query_pos, query.len(), i
                            ));
                        }
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        if r != q {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): expected match at ref_pos {} query_pos {}, but ref='{}' query='{}' (offset {})",
                                op_idx, n, ref_pos, query_pos, r as char, q as char, i
                            ));
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Kind::SequenceMismatch => {
                    for i in 0..n {
                        if ref_pos >= reference.len() {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): ref_pos {} exceeds reference length {} at offset {}",
                                op_idx, n, ref_pos, reference.len(), i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): query_pos {} exceeds query length {} at offset {}",
                                op_idx, n, query_pos, query.len(), i
                            ));
                        }
                        let r = reference[ref_pos];
                        let q = query[query_pos];
                        if r == q {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): expected mismatch at ref_pos {} query_pos {}, but both are '{}' (offset {})",
                                op_idx, n, ref_pos, query_pos, r as char, i
                            ));
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Kind::Insertion => {
                    let new_pos = query_pos + n;
                    if new_pos > query.len() {
                        return Err(format!(
                            "CIGAR op {} ({}{ch}): query_pos {} + {} exceeds query length {}",
                            op_idx, n, query_pos, n, query.len()
                        ));
                    }
                    query_pos = new_pos;
                }
                Kind::Deletion => {
                    let new_pos = ref_pos + n;
                    if new_pos > reference.len() {
                        return Err(format!(
                            "CIGAR op {} ({}{ch}): ref_pos {} + {} exceeds reference length {}",
                            op_idx, n, ref_pos, n, reference.len()
                        ));
                    }
                    ref_pos = new_pos;
                }
                Kind::SoftClip => {
                    let new_pos = query_pos + n;
                    if new_pos > query.len() {
                        return Err(format!(
                            "CIGAR op {} ({}{ch}): query_pos {} + {} exceeds query length {}",
                            op_idx, n, query_pos, n, query.len()
                        ));
                    }
                    query_pos = new_pos;
                }
                Kind::HardClip => {}
                _ => {}
            }
        }

        Ok(())
    }

    /// Compute a quality score from the CIGAR (higher is better).
    pub fn quality(&self, params: &AlignParams) -> QualityScore {
        let mut score = 0.0;
        for &op in &self.cigar {
            score += params.quality(op).0;
        }
        QualityScore::new(score)
    }

    #[allow(dead_code)]
    pub fn concat(alignments: &[Alignment]) -> Alignment {
        let mut total_divergence = DivergenceScore::ZERO;
        let mut combined_cigar = Vec::new();

        for aln in alignments {
            total_divergence = DivergenceScore::new(total_divergence.0 + aln.divergence.0);
            combined_cigar.extend_from_slice(&aln.cigar);
        }

        Alignment {
            divergence: total_divergence,
            cigar: combined_cigar,
        }
    }

    pub fn mismatch_count(&self) -> usize {
        self.cigar
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion
                )
            })
            .map(|op| op.len())
            .sum()
    }
}

impl From<Vec<Op>> for Alignment {
    fn from(cigar: Vec<Op>) -> Self {
        let mut divergence = 0usize;
        for &op in &cigar {
            match op.kind() {
                Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion => {
                    divergence += op.len()
                }
                _ => {}
            }
        }
        Self {
            divergence: DivergenceScore::new(divergence as f64),
            cigar,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parse a SAM CIGAR string into a vector of Ops.
    pub fn parse_cigar(cig: &str) -> Option<Vec<Op>> {
        let mut cigar = Vec::new();
        let mut count = 0usize;

        for c in cig.chars() {
            if let Some(d) = c.to_digit(10) {
                count = count * 10 + d as usize;
            } else {
                let kind = match c {
                    '=' => Kind::SequenceMatch,
                    'X' => Kind::SequenceMismatch,
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    'S' => Kind::SoftClip,
                    'H' => Kind::HardClip,
                    'M' => Kind::Match,
                    _ => return None,
                };
                cigar.push(Op::new(kind, count));
                count = 0;
            }
        }

        Some(cigar)
    }
    #[test]
    fn test_identical() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"ACGT", b"ACGT").unwrap();
        assert_eq!(result.divergence.0, 0.0);
        assert_eq!(result.cigar_string(), "4=");
    }

    #[test]
    fn test_single_mismatch() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"ACGT", b"ACTT").unwrap();
        assert_eq!(result.divergence.0, 2.0); // mismatch penalty
        assert_eq!(result.cigar_string(), "2=1X1=");
    }

    #[test]
    fn test_single_insertion() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"ACGT", b"ACT").unwrap();
        // query has extra G
        assert_eq!(result.divergence.0, 7.0);
        assert!(result.cigar_string().contains('I'));
    }

    #[test]
    fn test_single_deletion() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"ACT", b"ACGT").unwrap();
        // query missing G
        assert_eq!(result.divergence.0, 6.0);
        assert!(result.cigar_string().contains('D'));
    }

    #[test]
    fn test_empty() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"", b"").unwrap();
        assert_eq!(result.divergence.0, 0.0);
        assert!(result.cigar.is_empty());
    }

    #[test]
    fn test_query_empty() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"", b"ACGT").unwrap();
        assert_eq!(result.cigar_string(), "4D");
    }

    #[test]
    fn test_reference_empty() {
        let mut aligner = Aligner::with_defaults();
        let result = aligner.align(b"ACGT", b"").unwrap();
        assert_eq!(result.cigar_string(), "4I");
    }

    #[test]
    fn test_longer_sequences() {
        let mut aligner = Aligner::with_defaults();
        let query = b"ACGTACGTACGT";
        let reference = b"ACGTACGTACGT";
        let result = aligner.align(query, reference).unwrap();
        assert_eq!(result.divergence.0, 0.0);
        assert_eq!(result.cigar_string(), "12=");
    }

    #[test]
    fn test_with_gaps() {
        let mut aligner = Aligner::with_defaults();
        let query = b"ACGTACGT";
        let reference = b"ACGTTTTACGT";
        let result = aligner.align(query, reference).unwrap();
        assert!(result.divergence.0 > 0.0);
        // Should have a deletion in the middle
        println!("CIGAR: {}", result.cigar_string());
    }

    #[test]
    fn test_from_data_1() {
        let mut aligner = Aligner::with_defaults();
        let reference = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGAATAGTTAGTGCTTCTTGTTGAATTGATCCCTTTACCATTATGTAAATGGCTTTCTTTGTCTCTTTTGATCTTTTTGGTTTAAAATCTGTTTTATCAGAGACTATGACAGCAATCCCTGCTTTTTTTGCTTTTCATTTGCTTGGTAGATCTTCCTCTGTCCCTTTATTTTGAGCCCATGTGTGTGTCTGCACATGAGATGGGTCTCCTGAATACAGCACGTTGATGGGTCTCGAATCTTTGTCCAGTTTGTCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTGTTTTGCTCATTAGTTGATGCCGTTTCTTCCTAGCATCAACGGTCTTTACAATTTGGCCTGTTTTTGCAGTGGCTGGTACCAGTTGTTCCTTTCCATGTTTAATGCTTCCTTCAGGAGCTCTTGTAAGGCAGGCCTGGTGGTGACAAAATCTCTCAGGATTTGCTTGTCTGCAAAGGATTTTGTTTCTCCTTCACTTATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGTAAATTATTTTCTTTGAGAATGTTGAATATTGGCCCCCACTCTCATCTGGCTTGTAGGGTTTCTGCCGAGAGATCTGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTATCCCAGCCTTTCTCTCTGGCTGACCTTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTGGGGGTTGCTCTTCTCAAGGAGTGTCTTTGTAGTGTTCTCTGTATTTCTTGAATTTGAATGTTGGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGAATAATATCCTGAAGAGTGTTTTCCAGCTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGATCTTCTCACATAGTTCCATACTTCTTGGAGGCTTTGTTTGTTTCATTTTACTCTTTTTCTCTAAACTTCTCTTCTTGCTTCATTTCATTAATTTGATCTTCAGTCACTGAAACCCTTTCTTCCATTGATCGAATCAGCTACTGAAGCTTCTGTGTGTGTCACGTAGTTCTCGTGTCATGGTTTTCAGCTCCTTCAGGTCATTTAAGGTTTTCTCTACACTGGTTACTCTAGTTAGCCTTTTGTCTAATCTTTTTTCAAGGTTTTTAGCTTCCTTGCGGTGGGTTTGAATATCCTCCTTTAGCTCAGAGAAATTTGTTTTTACCGACCTTCTGAAGCCTAATTCTGTCAACTCGTCAAAGTCATTCTCCATCCAGCCTTGTTCTGTTGCTGGCGAGGAGCTGTGATCCTTTGGAGGAGAAGAGGCACTCTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTTGTTTTTATCTACCTTTGGTCTTTGATGATGGTGACCTACAGATGGGGTTATTGGTGTGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTAACAGTCAGGTCCCTCAGCTGCAGGTCTGTTGGAGTTTGCTAGAGGTCCACTCCAGACACTGTTTGCCTGGGTATCACCTTTGGAGGCTGCAGAACAGCAAATATTGCAGAACAGCAAATATTGCTGCCTGATCCTTCCTCTGGAAGCATCGTCCCAGAGGGGCATACGGCAGCATGAGATGTCAGTCAGCCCCCACTGGGAGGTGTCTCCCTGTTAGGCTACACGGGGGTCAGGGACCCACTTGAGGAGGCAGTCTGTCCGTTCTCAGAGCTCAAACACCGTGCTGGGAGAACCACTGCTCTCTTCAGAGCAGTGCAGACAAGGACATTTAAGTCTGCAGAAGTTTCTGCTGCCTTTTGTTCAGCTATGCCCTGCCACCAGAGGTGGAGTCTATAGAGGCAGCAAGCCTTGTGGTGCTGTGGTGGGCTCTGCCAAGTTCGAGCTTCCTGGCTGCTTTGTTTACCTACTTAAGCCTCAGCAATGGTGGACGCCCCTTCCCCAGCCAGGCTGC";
        let query = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGGGTAGTTAGCTCCCTTTACCATTATGTAATGGCCTTCTTTGTCCCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACAACACCTGCTTTTTTTGTTTTCCATTTGCTTGGTAGGTCTTCCTCCATCCCTTTATTTTGAGCCTATGTGTGTGTCTGCACATGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTCATCTAACTTGCCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATTCTGTCATTATGATGTTAGCTGGTTATTTTTCCCGTTAGTTGATGCAGTTTCTTCCTAGCATCGATGGTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGTGCTTCCTTCAGGAGCTCTTGTAAGGCAGGGCTGGTGGTGACAAAATCTCTCAGCATTTGCTGGTCTATAAAGGATTTTATTTCTCCTTATGAAGCTTTGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTAAGAATGTTGAATATTGGTGCCCACTCTCTTCTGACTTGTAGAGTTTCTGTTGAGAGATCCACTGTTAGTCTGATGGGCTTCCCTTTGTGGCTAACTCGACCTTTCTCTCTGGGTGCCATTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTTGGGGTTGCTCCTCTCGAGGAGCATCTTGGTAGTGTTCTCTGTATTTCCTGAGTTTGAATGTTTGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGGACAATATCCTGAAGAGTGTTTTCGAACTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGGTGTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGTTCATTCTTTTTACTCTTTTTTCTCTAAACTTCTCACTTCATTAATTTGATCTTCAATCACTGATACCCTTTCTTTCAGTTTATTGAATCAACTACTGAAGCTTGTGCATGTGTCACATAGTTCTTGTTCCATGGTTTTCAGCTCCATCAGGTCATTTAAGGTCTCCACACTGCTTATTCTAGTTAGCCATTCATCTAATCTGTTTGCAAGGCTTTTAGCTTCCTTGTGATGGGTTCGAATACCTCCCTTAACTCAGAGAAGTTTGTTATTACCAACCTTCTGAAGCCTACTTCTGTCAGCTCATCAAAGTCATTCTCCGTCCAGCTTTATTCCGTTGCTGGCAAGGAGCTGTAATCCTTTGCAGGAGAAGGGATGCTGTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTGGTTTTATCTACCTTTGGTCTTCGATGATGGTGACCCACAGATGGGGTTTTGGTGTGGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTGACAGTCAGGTCCCTCAGCTGCAGATCTGTTGGAGTTTGCTGGAGGTCCACTCCAGACTCTGTTTACCTGTGTATCACCAGCAGAGGCTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAACAGCAAATATTACTGCCTGATCCTTACTCTGGAAGCTTCATTTCAGAGGGGCACCCAGCTCTATGAGGTGTCATTCGGCCCCTACTGGGAGATGTCTCCCAGTTAGGCTACACAGGGGTCAGGGACACACTTGAGGAGGCAGTCTGTCCATTCTCAGAGCTCAAACTCCATGCTAGGAGAACCACTGCTCTCTTCAGAGCTGTCAGATAGGGACATTTAAGTCTGCAGAAGTTTCTGCTGCCTTTTGTTCAGCTATGCCCTGCCCCCAGAGGTGGAGTCTACAGAGTCAGGCAGGCCTCCTTGAGCTGTGGTGGGCTCCACCCAGTTCGAGCTTCCCAGCCGCTTTGTTTACCTACTCAAGCTTCAGCAATGGCGGACGCCCCTTCCCCAGCCAGGCTGC";
        let alignment = aligner.align(query, reference).unwrap();

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
        let mut aligner = Aligner::with_defaults();
        // Truncated version of test_from_data_1 for debugging
        // Cut at position ~1550 in reference (just before error at ref_pos 1648)
        let reference = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGAATAGTTAGTGCTTCTTGTTGAATTGATCCCTTTACCATTATGTAAATGGCTTTCTTTGTCTCTTTTGATCTTTTTGGTTTAAAATCTGTTTTATCAGAGACTATGACAGCAATCCCTGCTTTTTTTGCTTTTCATTTGCTTGGTAGATCTTCCTCTGTCCCTTTATTTTGAGCCCATGTGTGTGTCTGCACATGAGATGGGTCTCCTGAATACAGCACGTTGATGGGTCTCGAATCTTTGTCCAGTTTGTCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTGTTTTGCTCATTAGTTGATGCCGTTTCTTCCTAGCATCAACGGTCTTTACAATTTGGCCTGTTTTTGCAGTGGCTGGTACCAGTTGTTCCTTTCCATGTTTAATGCTTCCTTCAGGAGCTCTTGTAAGGCAGGCCTGGTGGTGACAAAATCTCTCAGGATTTGCTTGTCTGCAAAGGATTTTGTTTCTCCTTCACTTATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGTAAATTATTTTCTTTGAGAATGTTGAATATTGGCCCCCACTCTCATCTGGCTTGTAGGGTTTCTGCCGAGAGATCTGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTATCCCAGCCTTTCTCTCTGGCTGACCTTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTGGGGGTTGCTCTTCTCAAGGAGTGTCTTTGTAGTGTTCTCTGTATTTCTTGAATTTGAATGTTGGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGAATAATATCCTGAAGAGTGTTTTCCAGCTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGATCTTCTCACATAGTTCCATACTTCTTGGAGGCTTTGTTTGTTTCATTTTACTCTTTTTCTCTAAACTTCTCTTCTTGCTTCATTTCATTAATTTGATCTTCAGTCACTGAAACCCTTTCTTCCATTGATCGAATCAGCTACTGAAGCTTCTGTGTGTGTCACGTAGTTCTCGTGTCATGGTTTTCAGCTCCTTCAGGTCATTTAAGGTTTTCTCTACACTGGTTACTCTAGTTAGCCTTTTGTCTAATCTTTTTTCAAGGTTTTTAGCTTCCTTGCGGTGGGTTTGAATATCCTCCTTTAGCTCAGAGAAATTTGTTTTTACCGACCTTCTGAAGCCTAATTCTGTCAACTCGTCAAAGTCATTCTCCATCCAGCCTTGTTCTGTTGCTGGCGAGGAGCTGTGATCCTTTGGAGGAGAAGAGGCACTCTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTTGTTTTTATCTACCTTTGGTCTTTGATGATGGTGACCTACAGATGGGGTTATTGGTGTGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTAACAGTCAGGTCCCTCAGCTGCAGGTCTGTTGGAGTTTGCTAGAGGTCCACTCCAGACACTGTTTGCCTGGGTATCACCTTTGGAGGCTGCAGAACAGCAAATATTGCAGAACAGCAAATATTGC";
        let query = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGGGTAGTTAGCTCCCTTTACCATTATGTAATGGCCTTCTTTGTCCCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACAACACCTGCTTTTTTTGTTTTCCATTTGCTTGGTAGGTCTTCCTCCATCCCTTTATTTTGAGCCTATGTGTGTGTCTGCACATGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTCATCTAACTTGCCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATTCTGTCATTATGATGTTAGCTGGTTATTTTTCCCGTTAGTTGATGCAGTTTCTTCCTAGCATCGATGGTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGTGCTTCCTTCAGGAGCTCTTGTAAGGCAGGGCTGGTGGTGACAAAATCTCTCAGCATTTGCTGGTCTATAAAGGATTTTATTTCTCCTTATGAAGCTTTGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTAAGAATGTTGAATATTGGTGCCCACTCTCTTCTGACTTGTAGAGTTTCTGTTGAGAGATCCACTGTTAGTCTGATGGGCTTCCCTTTGTGGCTAACTCGACCTTTCTCTCTGGGTGCCATTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTTGGGGTTGCTCCTCTCGAGGAGCATCTTGGTAGTGTTCTCTGTATTTCCTGAGTTTGAATGTTTGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGGACAATATCCTGAAGAGTGTTTTCGAACTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGGTGTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGTTCATTCTTTTTACTCTTTTTTCTCTAAACTTCTCACTTCATTAATTTGATCTTCAATCACTGATACCCTTTCTTTCAGTTTATTGAATCAACTACTGAAGCTTGTGCATGTGTCACATAGTTCTTGTTCCATGGTTTTCAGCTCCATCAGGTCATTTAAGGTCTCCACACTGCTTATTCTAGTTAGCCATTCATCTAATCTGTTTGCAAGGCTTTTAGCTTCCTTGTGATGGGTTCGAATACCTCCCTTAACTCAGAGAAGTTTGTTATTACCAACCTTCTGAAGCCTACTTCTGTCAGCTCATCAAAGTCATTCTCCGTCCAGCTTTATTCCGTTGCTGGCAAGGAGCTGTAATCCTTTGCAGGAGAAGGGATGCTGTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTGGTTTTATCTACCTTTGGTCTTCGATGATGGTGACCCACAGATGGGGTTTTGGTGTGGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTGACAGTCAGGTCCCTCAGCTGCAGATCTGTTGGAGTTTGCTGGAGGTCCACTCCAGACTCTGTTTACCTGTGTATCACCAGCAGAGGCTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAACAGCAAATATTAC";
        let alignment = aligner.align(query, reference).unwrap();

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
    fn test_cigar_scoring() {
        let cig = "2D1=1X1=1X1=1X2=2X2=1X1=1X2=1X2=1X1=2X2=12X1=1X1=3X1=2X1=1X1=2X5=1X1=2X1=1X2=1X3=21D3=1X6=2X1=1X1=2X1=3X2=1X1=2X3=1X1=3X1=1D2=1X2=1X3=2X2=3X3=1X18D176=1X82=1D";
        let cigar = parse_cigar(cig).expect("bad cigar string");
        let params = AlignParams::default();
        let mut total = 0.0;
        for &op in &cigar {
            let q = params.quality(op).0;
            total += q;
            println!(
                "{:?} -> {:.2} (running total: {:.2})",
                op, q, total
            );
        }
        let alignment = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        assert_eq!(alignment.quality(&params).0, 290.0)
    }
}
