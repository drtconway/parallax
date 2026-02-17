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
    pub fn divergence(&self, op: &CigarOp) -> DivergenceScore {
        (match op {
            CigarOp::Match(_) => 0.0,
            CigarOp::Mismatch(n) => self.mismatch as f64 * (*n as f64),
            CigarOp::Ins(n) | CigarOp::Del(n) => {
                self.gap_open as f64 + self.gap_extend as f64 * (*n as f64)
            }
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0.0,
        })
        .into()
    }

    /// Score an individual alignment operation
    pub fn quality(&self, op: &CigarOp) -> QualityScore {
        (match op {
            CigarOp::Match(n) => self.match_score * (*n as i32),
            CigarOp::Mismatch(n) => -self.mismatch * (*n as i32),
            CigarOp::Ins(n) | CigarOp::Del(n) => -self.gap_open - self.gap_extend * (*n as i32),
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0,
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

/// CIGAR operation
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),    // '='
    Mismatch(u32), // 'X'
    Ins(u32),      // 'I' - insertion in query
    Del(u32),      // 'D' - deletion in query
    SoftClip(u32), // 'S' - soft clipped bases
    HardClip(u32), // 'H' - hard clipped bases (not consumed by alignment)
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
            CigarOp::HardClip(n) => format!("{}H", n),
        }
    }

    #[allow(dead_code)]
    pub fn divergence(&self, params: &AlignParams) -> DivergenceScore {
        params.divergence(self)
    }

    pub fn quality(&self, params: &AlignParams) -> QualityScore {
        params.quality(self)
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
    /// Edit distance (lower is better)
    #[allow(dead_code)]
    pub divergence: DivergenceScore,
    /// CIGAR operations
    pub cigar: Vec<CigarOp>,
}

impl Alignment {
    /// Create a new perfect match alignment
    #[allow(dead_code)]
    pub fn from_perfect_match(length: usize) -> Self {
        Self {
            divergence: DivergenceScore::ZERO,
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
        merged.iter().map(|op| op.to_string()).collect()
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
                CigarOp::Del(_) | CigarOp::HardClip(_) => 0,
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
                CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0,
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
                CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0,
            })
            .sum()
    }

    /// Compute the reference bases consumed by alignment operations.
    /// This is the sum of M, D, =, X operations.
    pub fn reference_consumed(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as usize,
                CigarOp::Mismatch(n) => *n as usize,
                CigarOp::Del(n) => *n as usize,
                CigarOp::Ins(_) => 0,
                CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0,
            })
            .sum()
    }

    /// Return the size of the leading hard clip, if any.
    pub fn leading_hard_clip(&self) -> usize {
        match self.cigar.first() {
            Some(CigarOp::HardClip(n)) => *n as usize,
            _ => 0,
        }
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
                CigarOp::HardClip(_) => {
                    // Hard clips aren't shown in alignment and don't consume query
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
                if let CigarOp::HardClip(_) = self.cigar[i] {
                    return Err(format!(
                        "CIGAR op {} is a hard clip in the middle of the alignment: {:?}",
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
                CigarOp::HardClip(_) => {
                    // Hard clips don't consume query or reference, so no position change
                }
            }
        }

        Ok(())
    }

    /// Compute a quality score from the CIGAR (higher is better).
    pub fn quality(&self, params: &AlignParams) -> QualityScore {
        let mut score = 0.0;
        for op in &self.cigar {
            score += op.quality(params).0;
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
            .map(|op| match op {
                CigarOp::Mismatch(n) | CigarOp::Ins(n) | CigarOp::Del(n) => *n as usize,
                _ => 0,
            })
            .sum()
    }
}

impl From<Vec<CigarOp>> for Alignment {
    fn from(cigar: Vec<CigarOp>) -> Self {
        // Compute divergence from CIGAR: count mismatches + indels
        let mut divergence = 0u32;
        for op in &cigar {
            match op {
                CigarOp::Mismatch(n) | CigarOp::Ins(n) | CigarOp::Del(n) => divergence += n,
                CigarOp::Match(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {}
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
        let cigar = CigarOp::make(cig).expect("bad cigar string");
        let params = AlignParams::default();
        let mut total = 0.0;
        for op in &cigar {
            total += op.quality(&params).0;
            println!(
                "{:?} -> {:.2} (running total: {:.2})",
                op,
                op.quality(&params).0,
                total
            );
        }
        let alignment = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        assert_eq!(alignment.quality(&params).0, 290.0)
    }
}
