//! Wavefront Alignment Algorithm (WFA) for sequence alignment.
//!
//! Implements affine-gap WFA as described by Santiago Marco-Sola et al.
//! Time complexity: O(ns) where n is sequence length and s is the alignment score.
//! This is very fast for similar sequences where s << n.

use std::cmp::{max, min};

use crate::config;

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
    pub fn query_length(&self) -> u64 {
        self.cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) => *n as u64,
                CigarOp::Mismatch(n) => *n as u64,
                CigarOp::Ins(n) => *n as u64,
                CigarOp::SoftClip(n) => *n as u64,
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
    #[allow(dead_code)]
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

        for (op_idx, op) in self.cigar.iter().enumerate() {
            match op {
                CigarOp::Match(n) => {
                    for i in 0..*n as usize {
                        if ref_pos >= reference.len() {
                            return Err(format!(
                                "CIGAR op {} ({}=): ref_pos {} exceeds reference length {} at offset {}",
                                op_idx, n, ref_pos, reference.len(), i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}=): query_pos {} exceeds query length {} at offset {}",
                                op_idx, n, query_pos, query.len(), i
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
                                op_idx, n, ref_pos, reference.len(), i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}X): query_pos {} exceeds query length {} at offset {}",
                                op_idx, n, query_pos, query.len(), i
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
                            op_idx, n, query_pos, n, query.len()
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
                            op_idx, n, ref_pos, n, reference.len()
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
                            op_idx, n, query_pos, n, query.len()
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
                },
                CigarOp::Mismatch(m) => {
                    let m = *m as f64;
                    score -= m * (1.0 + m.log2());
                },
                CigarOp::Ins(i) => {
                    let i = *i as f64;
                    score -= 4.0 + i * (1.0 + i.log2());
                },
                CigarOp::Del(d) => {
                    let d = *d as f64;
                    score -= 4.0 + d * (1.0 + d.log2());
                },
                CigarOp::SoftClip(_) => {
                    // Soft clips do not contribute to score
                },
            }
        }
        score
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
    if len >= min_len {
        len
    } else {
        0
    }
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

/// Wavefront for a single score, storing furthest-reaching point on each diagonal.
/// Diagonal k = row - col, so row = offset[k] and col = offset[k] - k.
#[derive(Clone)]
struct Wavefront {
    /// Furthest row reached on each diagonal k in [lo, hi]
    offsets: Vec<i32>,
    /// Minimum diagonal index
    lo: i32,
    /// Maximum diagonal index
    hi: i32,
}

impl Wavefront {
    fn new() -> Self {
        Self {
            offsets: Vec::new(),
            lo: 0,
            hi: -1, // empty
        }
    }

    fn with_diagonals(lo: i32, hi: i32) -> Self {
        let size = (hi - lo + 1) as usize;
        Self {
            offsets: vec![i32::MIN / 2; size], // -inf sentinel
            lo,
            hi,
        }
    }

    fn is_empty(&self) -> bool {
        self.hi < self.lo
    }

    fn get(&self, k: i32) -> i32 {
        if k < self.lo || k > self.hi {
            i32::MIN / 2
        } else {
            self.offsets[(k - self.lo) as usize]
        }
    }

    fn set(&mut self, k: i32, val: i32) {
        if k >= self.lo && k <= self.hi {
            self.offsets[(k - self.lo) as usize] = val;
        }
    }
}

/// Backtrace state for CIGAR reconstruction
#[derive(Clone, Copy, Debug)]
enum TraceOp {
    Mismatch,
    InsOpen,
    InsExt,
    DelOpen,
    DelExt,
}

/// Affine-gap WFA aligner
pub struct WfAligner {
    params: AlignParams,
    /// Maximum score to explore before giving up (prevents runaway on very different sequences)
    max_score: i32,
}

impl WfAligner {
    pub fn new(params: AlignParams) -> Self {
        Self {
            params,
            max_score: 10000,
        }
    }

    /// Align query to reference, returning the alignment.
    /// Uses affine-gap WFA with traceback.
    pub fn align(&self, query: &[u8], reference: &[u8]) -> Option<Alignment> {
        let n = query.len() as i32;
        let m = reference.len() as i32;

        if n == 0 && m == 0 {
            return Some(Alignment {
                score: 0,
                cigar: vec![],
            });
        }
        if n == 0 {
            return Some(Alignment {
                score: self.params.gap_open + self.params.gap_extend * m,
                cigar: vec![CigarOp::Del(m as u32)],
            });
        }
        if m == 0 {
            return Some(Alignment {
                score: self.params.gap_open + self.params.gap_extend * n,
                cigar: vec![CigarOp::Ins(n as u32)],
            });
        }

        let x = self.params.mismatch;
        let o = self.params.gap_open;
        let e = self.params.gap_extend;

        // Wavefronts for M (match/mismatch), I (insertion), D (deletion)
        // Indexed by score
        let mut wf_m: Vec<Wavefront> = Vec::new();
        let mut wf_i: Vec<Wavefront> = Vec::new();
        let mut wf_d: Vec<Wavefront> = Vec::new();

        // Traceback storage: for each (score, diagonal, component) -> operation
        let mut trace: Vec<Vec<(i32, TraceOp)>> = Vec::new(); // score -> [(k, op), ...]

        // Initialize score 0: start at origin on diagonal 0
        wf_m.push(Wavefront::with_diagonals(0, 0));
        wf_m[0].set(0, 0);
        wf_i.push(Wavefront::new());
        wf_d.push(Wavefront::new());
        trace.push(vec![]);

        // Extend matches at score 0
        self.extend(query, reference, &mut wf_m[0]);

        // Target: reach (n, m) which is diagonal k = n - m, offset = n
        let target_k = n - m;

        // Check if already done
        if wf_m[0].get(target_k) >= n {
            return Some(self.traceback(query, reference, 0, &wf_m, &wf_i, &wf_d, &trace));
        }

        for s in 1..=self.max_score {
            let s_idx = s as usize;

            // Determine wavefront bounds
            let mut lo = i32::MAX;
            let mut hi = i32::MIN;

            // From M[s-x], I[s-e], D[s-e], M[s-o-e], etc.
            if s >= x && !wf_m[s_idx - x as usize].is_empty() {
                lo = min(lo, wf_m[s_idx - x as usize].lo);
                hi = max(hi, wf_m[s_idx - x as usize].hi);
            }
            if s >= e && !wf_i[(s - e) as usize].is_empty() {
                lo = min(lo, wf_i[(s - e) as usize].lo - 1);
                hi = max(hi, wf_i[(s - e) as usize].hi - 1);
            }
            if s >= e && !wf_d[(s - e) as usize].is_empty() {
                lo = min(lo, wf_d[(s - e) as usize].lo + 1);
                hi = max(hi, wf_d[(s - e) as usize].hi + 1);
            }
            if s >= o + e && !wf_m[(s - o - e) as usize].is_empty() {
                lo = min(lo, wf_m[(s - o - e) as usize].lo - 1);
                hi = max(hi, wf_m[(s - o - e) as usize].hi + 1);
            }

            if lo > hi {
                // No valid wavefront at this score
                wf_m.push(Wavefront::new());
                wf_i.push(Wavefront::new());
                wf_d.push(Wavefront::new());
                trace.push(vec![]);
                continue;
            }

            // Clamp to valid diagonal range
            lo = max(lo, -m);
            hi = min(hi, n);

            let mut new_m = Wavefront::with_diagonals(lo, hi);
            let mut new_i = Wavefront::with_diagonals(lo, hi);
            let mut new_d = Wavefront::with_diagonals(lo, hi);
            let mut trace_s = Vec::new();

            for k in lo..=hi {
                // Insertion: I[s][k] = max(M[s-o-e][k-1] + 1, I[s-e][k-1] + 1)
                let i_from_m = if s >= o + e {
                    wf_m[(s - o - e) as usize].get(k - 1) + 1
                } else {
                    i32::MIN / 2
                };
                let i_from_i = if s >= e {
                    wf_i[(s - e) as usize].get(k - 1) + 1
                } else {
                    i32::MIN / 2
                };
                let i_val = max(i_from_m, i_from_i);
                new_i.set(k, i_val);

                // Deletion: D[s][k] = max(M[s-o-e][k+1], D[s-e][k+1])
                let d_from_m = if s >= o + e {
                    wf_m[(s - o - e) as usize].get(k + 1)
                } else {
                    i32::MIN / 2
                };
                let d_from_d = if s >= e {
                    wf_d[(s - e) as usize].get(k + 1)
                } else {
                    i32::MIN / 2
                };
                let d_val = max(d_from_m, d_from_d);
                new_d.set(k, d_val);

                // Match/Mismatch: M[s][k] = max(M[s-x][k] + 1, I[s][k], D[s][k])
                let m_from_m = if s >= x {
                    wf_m[(s - x) as usize].get(k) + 1
                } else {
                    i32::MIN / 2
                };
                let m_val = max(max(m_from_m, i_val), d_val);
                new_m.set(k, m_val);

                // Record traceback
                if m_val > i32::MIN / 2 {
                    let op = if m_val == i_val {
                        if i_val == i_from_i {
                            TraceOp::InsExt
                        } else {
                            TraceOp::InsOpen
                        }
                    } else if m_val == d_val {
                        if d_val == d_from_d {
                            TraceOp::DelExt
                        } else {
                            TraceOp::DelOpen
                        }
                    } else {
                        TraceOp::Mismatch
                    };
                    trace_s.push((k, op));
                }
            }

            // Extend matches
            self.extend(query, reference, &mut new_m);

            wf_m.push(new_m);
            wf_i.push(new_i);
            wf_d.push(new_d);
            trace.push(trace_s);

            // Check if we reached the target
            if wf_m[s_idx].get(target_k) >= n {
                return Some(self.traceback(query, reference, s, &wf_m, &wf_i, &wf_d, &trace));
            }
        }

        None // Exceeded max score
    }

    /// Greedy match extension on all diagonals
    fn extend(&self, query: &[u8], reference: &[u8], wf: &mut Wavefront) {
        use crate::utils::longest_common_prefix;

        let n = query.len() as i32;
        let m = reference.len() as i32;

        for k in wf.lo..=wf.hi {
            let mut row = wf.get(k);
            let col = row - k;

            // Check bounds before calling LCP
            if row >= 0 && row < n && col >= 0 && col < m {
                let query_slice = &query[row as usize..];
                let ref_slice = &reference[col as usize..];
                let lcp = longest_common_prefix(query_slice, ref_slice);
                row += lcp as i32;
            }

            wf.set(k, row);
        }
    }

    /// Traceback to produce CIGAR using stored trace information
    fn traceback(
        &self,
        query: &[u8],
        reference: &[u8],
        final_score: i32,
        wf_m: &[Wavefront],
        wf_i: &[Wavefront],
        wf_d: &[Wavefront],
        trace: &[Vec<(i32, TraceOp)>],
    ) -> Alignment {
        let n = query.len() as i32;
        let m = reference.len() as i32;

        let x = self.params.mismatch;
        let o = self.params.gap_open;
        let e = self.params.gap_extend;

        let mut cigar = Vec::new();
        let mut row = n;
        let mut col = m;
        let mut s = final_score;

        // State: which wavefront are we currently in?
        #[derive(Clone, Copy, PartialEq)]
        enum State {
            M,
            I,
            D,
        }
        let mut state = State::M;

        while row > 0 || col > 0 {
            let k = row - col;

            if s == 0 {
                // All remaining must be matches (score 0 means we extended from origin)
                if row > 0 && col > 0 {
                    // Verify they actually match
                    let match_count = row.min(col);
                    cigar.push(CigarOp::Match(match_count as u32));
                }
                break;
            }

            if s < 0 {
                break;
            }

            let s_idx = s as usize;

            match state {
                State::M => {
                    let current_row = wf_m[s_idx].get(k);

                    // First, check how many matches were extended at this score
                    let row_before_extend =
                        self.row_before_extend_full(s, k, wf_m, wf_i, wf_d, x, o, e);
                    let matches = (current_row - row_before_extend).max(0).min(row.min(col));

                    if matches > 0 {
                        cigar.push(CigarOp::Match(matches as u32));
                        row -= matches;
                        col -= matches;
                    }

                    if row <= 0 && col <= 0 {
                        break;
                    }

                    // Now find what operation brought us to row_before_extend
                    // Look up the trace for this score/k
                    if let Some(op) = self.find_trace_op(s, k, trace) {
                        match op {
                            TraceOp::Mismatch => {
                                if s >= x && row > 0 && col > 0 {
                                    cigar.push(CigarOp::Mismatch(1));
                                    row -= 1;
                                    col -= 1;
                                    s -= x;
                                    // Stay in M state
                                } else {
                                    break;
                                }
                            }
                            TraceOp::InsOpen => {
                                // This M cell came from I cell at same score
                                // The I cell was opened from M[s-o-e][k-1]
                                if row > 0 {
                                    cigar.push(CigarOp::Ins(1));
                                    row -= 1;
                                    s -= o + e;
                                    // Go back to M state (gap open comes from M)
                                } else {
                                    break;
                                }
                            }
                            TraceOp::InsExt => {
                                // This M cell came from I cell at same score
                                // The I cell was extended from I[s-e][k-1]
                                if row > 0 {
                                    cigar.push(CigarOp::Ins(1));
                                    row -= 1;
                                    s -= e;
                                    state = State::I; // Continue in I state
                                } else {
                                    break;
                                }
                            }
                            TraceOp::DelOpen => {
                                if col > 0 {
                                    cigar.push(CigarOp::Del(1));
                                    col -= 1;
                                    s -= o + e;
                                    // Go back to M state
                                } else {
                                    break;
                                }
                            }
                            TraceOp::DelExt => {
                                if col > 0 {
                                    cigar.push(CigarOp::Del(1));
                                    col -= 1;
                                    s -= e;
                                    state = State::D; // Continue in D state
                                } else {
                                    break;
                                }
                            }
                        }
                    } else {
                        // No trace found - try to infer from wavefronts
                        if s >= x && row > 0 && col > 0 {
                            let prev_row = wf_m[(s - x) as usize].get(k);
                            if prev_row + 1 == row_before_extend {
                                cigar.push(CigarOp::Mismatch(1));
                                row -= 1;
                                col -= 1;
                                s -= x;
                                continue;
                            }
                        }
                        // Try gap transitions
                        if s >= o + e && row > 0 {
                            cigar.push(CigarOp::Ins(1));
                            row -= 1;
                            s -= o + e;
                            continue;
                        }
                        if s >= o + e && col > 0 {
                            cigar.push(CigarOp::Del(1));
                            col -= 1;
                            s -= o + e;
                            continue;
                        }
                        break;
                    }
                }
                State::I => {
                    // We're tracing back through insertion states
                    // I[s][k] came from either I[s-e][k-1] (extend) or M[s-o-e][k-1] (open)
                    let i_from_i = if s >= e {
                        wf_i[(s - e) as usize].get(k - 1)
                    } else {
                        i32::MIN / 2
                    };
                    let i_from_m = if s >= o + e {
                        wf_m[(s - o - e) as usize].get(k - 1) + 1
                    } else {
                        i32::MIN / 2
                    };

                    if i_from_i >= i_from_m && s >= e {
                        // Extended from I
                        if row > 0 {
                            cigar.push(CigarOp::Ins(1));
                            row -= 1;
                            s -= e;
                            // Stay in I state
                        } else {
                            break;
                        }
                    } else if s >= o + e {
                        // Opened from M
                        if row > 0 {
                            cigar.push(CigarOp::Ins(1));
                            row -= 1;
                            s -= o + e;
                            state = State::M;
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                State::D => {
                    // D[s][k] came from either D[s-e][k+1] (extend) or M[s-o-e][k+1] (open)
                    let d_from_d = if s >= e {
                        wf_d[(s - e) as usize].get(k + 1)
                    } else {
                        i32::MIN / 2
                    };
                    let d_from_m = if s >= o + e {
                        wf_m[(s - o - e) as usize].get(k + 1)
                    } else {
                        i32::MIN / 2
                    };

                    if d_from_d >= d_from_m && s >= e {
                        // Extended from D
                        if col > 0 {
                            cigar.push(CigarOp::Del(1));
                            col -= 1;
                            s -= e;
                            // Stay in D state
                        } else {
                            break;
                        }
                    } else if s >= o + e {
                        // Opened from M
                        if col > 0 {
                            cigar.push(CigarOp::Del(1));
                            col -= 1;
                            s -= o + e;
                            state = State::M;
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
            }
        }

        cigar.reverse();
        let mut alignment = Alignment {
            score: final_score,
            cigar,
        };
        alignment.normalize();
        alignment
    }

    /// Find the trace operation for a given score and diagonal
    fn find_trace_op(&self, s: i32, k: i32, trace: &[Vec<(i32, TraceOp)>]) -> Option<TraceOp> {
        if s <= 0 || s as usize >= trace.len() {
            return None;
        }
        for &(tk, op) in &trace[s as usize] {
            if tk == k {
                return Some(op);
            }
        }
        None
    }

    /// Calculate row before extension, considering all wavefronts
    fn row_before_extend_full(
        &self,
        s: i32,
        k: i32,
        wf_m: &[Wavefront],
        wf_i: &[Wavefront],
        wf_d: &[Wavefront],
        x: i32,
        _o: i32,
        _e: i32,
    ) -> i32 {
        let mut best = i32::MIN / 2;

        // From mismatch: M[s-x][k] + 1
        if s >= x {
            best = max(best, wf_m[(s - x) as usize].get(k) + 1);
        }

        // From insertion: I[s][k] (which equals wf_i value at this score)
        if (s as usize) < wf_i.len() {
            best = max(best, wf_i[s as usize].get(k));
        }

        // From deletion: D[s][k]
        if (s as usize) < wf_d.len() {
            best = max(best, wf_d[s as usize].get(k));
        }

        best
    }
}

/// Convenience function for quick alignment with default parameters
pub fn align(query: &[u8], reference: &[u8]) -> Option<Alignment> {
    metrics::histogram!("align_ref_len").record(reference.len() as f64);
    metrics::histogram!("align_query_len").record(query.len() as f64);
    let start = std::time::Instant::now();
    let result = WfAligner::new(AlignParams::default()).align(query, reference);
    metrics::histogram!("wf_align_time_us").record(start.elapsed().as_micros() as f64);
    result
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
        let result = align(b"ACGT", b"ACTT").unwrap();
        assert_eq!(result.score, 4); // mismatch penalty
        assert_eq!(result.cigar_string(), "2=1X1=");
    }

    #[test]
    fn test_single_insertion() {
        let result = align(b"ACGT", b"ACT").unwrap();
        // query has extra G
        assert!(result.score > 0);
        assert!(result.cigar_string().contains('I'));
    }

    #[test]
    fn test_single_deletion() {
        let result = align(b"ACT", b"ACGT").unwrap();
        // query missing G
        assert!(result.score > 0);
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
            cigar: vec![
                CigarOp::Match(3),
                CigarOp::Del(2),
                CigarOp::Match(5),
            ],
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
        let long_ref = b"ACGTACGTACXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXACGTACGTAC";

        let long_result = context_aware_score(&long_gap, long_ref, query, &params);

        // Long gap should cost less than 10x short gap (due to sublinear)
        assert!(long_result.score < short_result.score * 5);
    }
}
