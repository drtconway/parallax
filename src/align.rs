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
            mismatch: cfg.alignment.mismatch,
            gap_open: cfg.alignment.gap_open,
            gap_extend: cfg.alignment.gap_extend,
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
    pub fn score(&self) -> f64 {
        let mut score = 0.0;
        for op in &self.cigar {
            match op {
                CigarOp::Match(m) => {
                    let m = *m as f64;
                    score += m * (1.0 + m.log2());
                }
                CigarOp::Mismatch(m) => {
                    let m = *m as f64;
                    score -= m * (1.0 + m.log2());
                }
                CigarOp::Ins(i) => {
                    let i = *i as f64;
                    score -= 4.0 + i * (1.0 + i.log2());
                }
                CigarOp::Del(d) => {
                    let d = *d as f64;
                    score -= 4.0 + d * (1.0 + d.log2());
                }
                CigarOp::SoftClip(_) => {
                    // Soft clips do not contribute to score
                }
            }
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

/// Parameters for context-aware alignment scoring.
///
/// This scoring model accounts for:
/// 1. Sublinear gap extension costs (long gaps penalized less per base)
/// 2. Reduced penalties for indels in homopolymer/STR regions (systematic sequencing errors)
#[derive(Clone, Copy, Debug)]
pub struct ContextAwareParams {
    /// Mismatch penalty
    pub mismatch: i32,
    /// Gap open penalty (first base of gap)
    pub gap_open: i32,
    /// Gap extend penalty for short gaps (linear portion)
    pub gap_extend: i32,
    /// Gap length threshold where sublinear scaling kicks in
    pub sublinear_threshold: u32,
    /// Sublinear coefficient: penalty = gap_open + gap_extend * threshold + sublinear_coef * log2(len - threshold + 1)
    pub sublinear_coef: f64,
    /// Minimum homopolymer length to trigger reduced penalty
    pub homopolymer_min_len: usize,
    /// Penalty multiplier for gaps in homopolymer context (0.0 - 1.0)
    pub homopolymer_discount: f64,
    /// Minimum repeat unit count for STR discount (e.g., 3 means ATATAT or CAGCAGCAG)
    pub str_min_repeats: usize,
    /// Penalty multiplier for gaps in STR context (0.0 - 1.0)
    pub str_discount: f64,
}

impl Default for ContextAwareParams {
    fn default() -> Self {
        let cfg = config::get();
        Self {
            mismatch: cfg.alignment.mismatch,
            gap_open: cfg.alignment.gap_open,
            gap_extend: cfg.alignment.gap_extend,
            sublinear_threshold: cfg.alignment.sublinear_threshold,
            sublinear_coef: cfg.alignment.sublinear_coef,
            homopolymer_min_len: cfg.alignment.homopolymer_min_len,
            homopolymer_discount: cfg.alignment.homopolymer_discount,
            str_min_repeats: cfg.alignment.str_min_repeats,
            str_discount: cfg.alignment.str_discount,
        }
    }
}

/// Result of context-aware scoring with detailed breakdown
#[derive(Clone, Debug)]
pub struct ContextAwareScore {
    /// Total adjusted score (lower is better)
    pub score: i32,
    /// Number of matches
    pub matches: u32,
    /// Number of mismatches
    pub mismatches: u32,
    /// Total gap bases (insertions + deletions)
    pub gap_bases: u32,
    /// Gap bases in homopolymer context
    pub homopolymer_gap_bases: u32,
    /// Gap bases in STR context
    pub identity: f64,
}

/// Detect if position is in a homopolymer run in the reference.
/// Returns the length of the homopolymer if >= min_len, otherwise 0.
fn detect_homopolymer(seq: &[u8], pos: usize, min_len: usize) -> usize {
    if pos >= seq.len() {
        return 0;
    }

    let base = seq[pos];

    // Look backward for start of run
    let mut start = pos;
    while start > 0 && seq[start - 1] == base {
        start -= 1;
    }

    // Look forward for end of run
    let mut end = pos;
    while end + 1 < seq.len() && seq[end + 1] == base {
        end += 1;
    }

    let len = end - start + 1;
    if len >= min_len { len } else { 0 }
}

/// Detect if position is in a short tandem repeat (STR) region.
/// Checks for di- and tri-nucleotide repeats.
/// Returns (unit_size, repeat_count) if >= min_repeats, otherwise (0, 0).
fn detect_str(seq: &[u8], pos: usize, min_repeats: usize) -> (usize, usize) {
    if pos >= seq.len() {
        return (0, 0);
    }

    // Try dinucleotide repeats
    if let Some(count) = count_repeats(seq, pos, 2, min_repeats) {
        return (2, count);
    }

    // Try trinucleotide repeats
    if let Some(count) = count_repeats(seq, pos, 3, min_repeats) {
        return (3, count);
    }

    (0, 0)
}

/// Count repeat units of given size surrounding position
fn count_repeats(seq: &[u8], pos: usize, unit_size: usize, min_repeats: usize) -> Option<usize> {
    if pos + unit_size > seq.len() {
        return None;
    }

    // Get the repeat unit at this position
    let unit = &seq[pos..pos + unit_size];

    // Don't count homopolymers as STRs (handled separately)
    if unit.iter().all(|&b| b == unit[0]) {
        return None;
    }

    // Find start of repeat region
    let mut start = pos;
    while start >= unit_size {
        let prev_start = start - unit_size;
        if &seq[prev_start..start] == unit {
            start = prev_start;
        } else {
            break;
        }
    }

    // Count repeats forward from start
    let mut count = 0;
    let mut p = start;
    while p + unit_size <= seq.len() && &seq[p..p + unit_size] == unit {
        count += 1;
        p += unit_size;
    }

    if count >= min_repeats {
        Some(count)
    } else {
        None
    }
}

/// Calculate gap penalty with sublinear extension for long gaps.
///
/// For gaps <= threshold: gap_open + gap_extend * len
/// For gaps > threshold:  gap_open + gap_extend * threshold + sublinear_coef * log2(len - threshold + 1)
fn sublinear_gap_penalty(len: u32, params: &ContextAwareParams) -> f64 {
    if len == 0 {
        return 0.0;
    }

    let len = len as f64;
    let threshold = params.sublinear_threshold as f64;

    if len <= threshold {
        params.gap_open as f64 + params.gap_extend as f64 * len
    } else {
        let linear_part = params.gap_open as f64 + params.gap_extend as f64 * threshold;
        let extra = len - threshold;
        linear_part + params.sublinear_coef * (extra + 1.0).log2()
    }
}

/// Re-score an alignment with context-aware penalties.
///
/// This function walks through the CIGAR and reference/query sequences,
/// applying:
/// - Sublinear gap extension costs for long gaps
/// - Reduced penalties for gaps in homopolymer regions
/// - Reduced penalties for gaps in STR regions
///
/// # Arguments
/// * `alignment` - The alignment to re-score
/// * `reference` - Reference sequence (the portion covered by the alignment)
/// * `query` - Query sequence (the portion covered by the alignment, excluding soft clips)
/// * `params` - Scoring parameters
///
/// # Returns
/// Detailed scoring breakdown including adjusted score
pub fn context_aware_score(
    alignment: &Alignment,
    reference: &[u8],
    query: &[u8],
    params: &ContextAwareParams,
) -> ContextAwareScore {
    let mut score = 0.0f64;
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    let mut gap_bases = 0u32;
    let mut homopolymer_gap_bases = 0u32;
    let mut str_gap_bases = 0u32;

    let mut ref_pos = 0usize;
    let mut query_pos = 0usize;

    for op in &alignment.cigar {
        match op {
            CigarOp::Match(n) => {
                matches += n;
                ref_pos += *n as usize;
                query_pos += *n as usize;
            }
            CigarOp::Mismatch(n) => {
                mismatches += n;
                score += (*n as f64) * (params.mismatch as f64);
                ref_pos += *n as usize;
                query_pos += *n as usize;
            }
            CigarOp::Ins(n) => {
                // Insertion: query has extra bases, check query context
                gap_bases += n;

                // Check context at insertion point in query
                let homo_len = detect_homopolymer(query, query_pos, params.homopolymer_min_len);
                let (str_unit, str_count) = detect_str(query, query_pos, params.str_min_repeats);

                let base_penalty = sublinear_gap_penalty(*n, params);

                let discount = if homo_len > 0 {
                    homopolymer_gap_bases += n;
                    params.homopolymer_discount
                } else if str_unit > 0 && str_count >= params.str_min_repeats {
                    str_gap_bases += n;
                    params.str_discount
                } else {
                    1.0
                };

                score += base_penalty * discount;
                query_pos += *n as usize;
            }
            CigarOp::Del(n) => {
                // Deletion: reference has extra bases, check reference context
                gap_bases += n;

                // Check context at deletion point in reference
                let homo_len = detect_homopolymer(reference, ref_pos, params.homopolymer_min_len);
                let (str_unit, str_count) = detect_str(reference, ref_pos, params.str_min_repeats);

                let base_penalty = sublinear_gap_penalty(*n, params);

                let discount = if homo_len > 0 {
                    homopolymer_gap_bases += n;
                    params.homopolymer_discount
                } else if str_unit > 0 && str_count >= params.str_min_repeats {
                    str_gap_bases += n;
                    params.str_discount
                } else {
                    1.0
                };

                score += base_penalty * discount;
                ref_pos += *n as usize;
            }
            CigarOp::SoftClip(n) => {
                // Soft clips don't contribute to alignment score
                query_pos += *n as usize;
            }
        }
    }

    let aligned_length = matches + mismatches + gap_bases;
    let identity = if aligned_length > 0 {
        matches as f64 / aligned_length as f64
    } else {
        0.0
    };

    ContextAwareScore {
        score: score.round() as i32,
        matches,
        mismatches,
        gap_bases,
        homopolymer_gap_bases,
        identity,
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
    fn test_detect_homopolymer() {
        // AAAAA is a 5bp homopolymer
        let seq = b"ACGAAAAAACGT";
        assert_eq!(detect_homopolymer(seq, 3, 4), 6); // pos 3 is in AAAAAA
        assert_eq!(detect_homopolymer(seq, 5, 4), 6); // pos 5 is also in it
        assert_eq!(detect_homopolymer(seq, 0, 4), 0); // A at pos 0 is alone
        assert_eq!(detect_homopolymer(seq, 1, 4), 0); // C at pos 1 is alone
    }

    #[test]
    fn test_detect_str() {
        // ATATAT is 3 repeats of AT
        let seq = b"ACGATATATCGT";
        let (unit, count) = detect_str(seq, 3, 3);
        assert_eq!(unit, 2);
        assert_eq!(count, 3);

        // CAGCAGCAG is 3 repeats of CAG
        let seq2 = b"ACGCAGCAGCAGTTT";
        let (unit2, count2) = detect_str(seq2, 3, 3);
        assert_eq!(unit2, 3);
        assert_eq!(count2, 3);

        // Not enough repeats
        let seq3 = b"ACGATATTTT";
        let (unit3, count3) = detect_str(seq3, 3, 3);
        assert_eq!(unit3, 0);
        assert_eq!(count3, 0);
    }

    #[test]
    fn test_sublinear_gap_penalty() {
        let params = ContextAwareParams::default();

        // Short gap: linear
        let short_penalty = sublinear_gap_penalty(5, &params);
        assert_eq!(short_penalty, 6.0 + 2.0 * 5.0); // gap_open + gap_extend * len

        // At threshold: still linear
        let at_thresh = sublinear_gap_penalty(10, &params);
        assert_eq!(at_thresh, 6.0 + 2.0 * 10.0);

        // Long gap: sublinear
        let long_penalty = sublinear_gap_penalty(20, &params);
        let expected = 6.0 + 2.0 * 10.0 + 4.0 * (11.0f64).log2();
        assert!((long_penalty - expected).abs() < 0.001);

        // Very long gap shouldn't be proportionally more expensive
        let very_long = sublinear_gap_penalty(100, &params);
        // 100bp gap should cost much less than 10x a 10bp gap
        assert!(very_long < at_thresh * 4.0);
    }

    #[test]
    fn test_context_aware_score_basic() {
        // Simple alignment with no gaps
        let alignment = Alignment {
            score: 0,
            cigar: vec![CigarOp::Match(10)],
        };
        let reference = b"ACGTACGTAC";
        let query = b"ACGTACGTAC";
        let params = ContextAwareParams::default();

        let result = context_aware_score(&alignment, reference, query, &params);
        assert_eq!(result.score, 0);
        assert_eq!(result.matches, 10);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.gap_bases, 0);
        assert!((result.identity - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_context_aware_score_homopolymer_discount() {
        // Deletion in homopolymer region should get discount
        let alignment = Alignment {
            score: 10, // original WFA score
            cigar: vec![
                CigarOp::Match(3),
                CigarOp::Del(2), // deletion in AAAAA region
                CigarOp::Match(5),
            ],
        };
        // Reference has AAAAA, query is missing 2 A's
        let reference = b"ACGAAAAACGT";
        let query = b"ACGAAACGT"; // this is what the alignment covers
        let params = ContextAwareParams::default();

        let result = context_aware_score(&alignment, reference, query, &params);
        assert_eq!(result.homopolymer_gap_bases, 2);

        // Compare with non-homopolymer deletion
        let alignment2 = Alignment {
            score: 10,
            cigar: vec![CigarOp::Match(3), CigarOp::Del(2), CigarOp::Match(5)],
        };
        let reference2 = b"ACGXYZYZXCGT"; // no homopolymer
        let result2 = context_aware_score(&alignment2, reference2, query, &params);
        assert_eq!(result2.homopolymer_gap_bases, 0);

        // Homopolymer gap should have lower score
        assert!(result.score < result2.score);
    }

    #[test]
    fn test_context_aware_score_long_gap_sublinear() {
        let params = ContextAwareParams::default();

        // Short gap
        let short_gap = Alignment {
            score: 0,
            cigar: vec![CigarOp::Match(10), CigarOp::Del(5), CigarOp::Match(10)],
        };
        let short_ref = b"ACGTACGTACXXXXXACGTACGTAC";
        let query = b"ACGTACGTACACGTACGTAC";

        let short_result = context_aware_score(&short_gap, short_ref, query, &params);

        // Long gap (50bp)
        let long_gap = Alignment {
            score: 0,
            cigar: vec![CigarOp::Match(10), CigarOp::Del(50), CigarOp::Match(10)],
        };
        let long_ref =
            b"ACGTACGTACXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXACGTACGTAC";

        let long_result = context_aware_score(&long_gap, long_ref, query, &params);

        // Long gap should cost less than 10x short gap (due to sublinear)
        assert!(long_result.score < short_result.score * 5);
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
}
