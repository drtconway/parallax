//! WFA2 aligner with two-piece affine gap penalties.
//!
//! This wraps the WFA2-lib C library (via `biodiff-wfa2-sys`) to provide
//! global alignment using the gap-affine-2piece distance metric.
//!
//! The two-piece model uses:
//!   cost(gap of length k) = min(o1 + e1·k, o2 + e2·k)
//!
//! where (o1,e1) are the "short gap" penalties and (o2,e2) are the "long gap"
//! penalties. This concave cost function strongly favours consolidating
//! multiple small indels into a single large indel — exactly what's needed
//! for tandem repeat regions where reads carry extra copies.
//!
//! The WFA2 aligner object is reusable across multiple `align()` calls
//! (it resets internal state on each call).

use biodiff_wfa2_sys::*;
use std::os::raw::c_char;

use super::{Alignment, Kind, Op};
use crate::scores::DivergenceScore;

/// Error types for WFA2 operations.
#[derive(Debug, Clone)]
pub enum Wfa2Error {
    /// Alignment did not complete (status code from WFA2).
    Failed(i32),
}

impl std::fmt::Display for Wfa2Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Wfa2Error::Failed(code) => write!(f, "WFA2 alignment failed with status {}", code),
        }
    }
}

impl std::error::Error for Wfa2Error {}

/// Configuration for the WFA2 aligner.
#[derive(Debug, Clone)]
pub struct Wfa2Config {
    pub mismatch: i32,
    /// Short-gap open penalty (o1, minimap2 default: 4)
    pub gap_open1: i32,
    /// Short-gap extend penalty (e1, minimap2 default: 2)
    pub gap_extend1: i32,
    /// Long-gap open penalty (o2, minimap2 default: 24)
    pub gap_open2: i32,
    /// Long-gap extend penalty (e2, minimap2 default: 1)
    pub gap_extend2: i32,
}

impl Default for Wfa2Config {
    fn default() -> Self {
        Self {
            mismatch: 4,
            gap_open1: 4,
            gap_extend1: 2,
            gap_open2: 24,
            gap_extend2: 1,
        }
    }
}

/// Safe wrapper around the WFA2 C library aligner.
///
/// The underlying `wavefront_aligner_t` is created once and reused
/// across multiple `align()` calls. It is not `Send` or `Sync`
/// because the C library is not thread-safe for a single aligner instance.
pub struct Wfa2Aligner {
    wf_aligner: *mut wavefront_aligner_t,
}

// The WFA2 aligner holds a raw C pointer that is only used from
// a single thread (the per-read alignment pipeline).
// We mark it Send so it can be stored in the Aligner struct which
// is created per-thread.
unsafe impl Send for Wfa2Aligner {}

impl Drop for Wfa2Aligner {
    fn drop(&mut self) {
        if !self.wf_aligner.is_null() {
            unsafe {
                wavefront_aligner_delete(self.wf_aligner);
            }
        }
    }
}

impl Wfa2Aligner {
    /// Create a new WFA2 aligner with two-piece affine gap penalties.
    pub fn new(config: &Wfa2Config) -> Self {
        unsafe {
            let mut attrs = wavefront_aligner_attr_default;

            // Two-piece affine gap model
            attrs.distance_metric = distance_metric_t_gap_affine_2p;
            attrs.alignment_scope = alignment_scope_t_compute_alignment;

            // Set 2-piece penalties.
            //
            // WFA2 convention: penalties are positive values representing costs.
            // match_ = 0 means exact-match costs nothing.
            attrs.affine2p_penalties = affine2p_penalties_t {
                match_: 0,
                mismatch: config.mismatch,
                gap_opening1: config.gap_open1,
                gap_extension1: config.gap_extend1,
                gap_opening2: config.gap_open2,
                gap_extension2: config.gap_extend2,
            };

            // Use low memory mode for production workloads.
            attrs.memory_mode = wavefront_memory_t_wavefront_memory_med;

            let aligner = wavefront_aligner_new(&mut attrs as *mut _);
            assert!(!aligner.is_null(), "wavefront_aligner_new returned null");

            // Global end-to-end alignment
            wavefront_aligner_set_alignment_end_to_end(aligner);

            // Use adaptive WF heuristic for large sequences
            wavefront_aligner_set_heuristic_wfadaptive(aligner, 10, 50, 1);

            Self {
                wf_aligner: aligner,
            }
        }
    }

    /// Create with default (minimap2-like) parameters.
    pub fn with_defaults() -> Self {
        Self::new(&Wfa2Config::default())
    }

    /// Perform global alignment of query against reference.
    ///
    /// WFA2 convention: "pattern" = query, "text" = reference.
    /// Returns an `Alignment` with extended CIGAR (=, X, I, D).
    pub fn align(&mut self, query: &[u8], reference: &[u8]) -> Result<Alignment, Wfa2Error> {
        if query.is_empty() && reference.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::ZERO,
                cigar: Vec::new(),
            });
        }

        // Handle edge case: one sequence empty → pure insertion or deletion
        if query.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::ZERO,
                cigar: vec![Op::new(Kind::Deletion, reference.len())],
            });
        }
        if reference.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::ZERO,
                cigar: vec![Op::new(Kind::Insertion, query.len())],
            });
        }

        unsafe {
            // WFA2 convention: pattern = reference, text = query.
            // This way WFA2 'I' = insertion in query (SAM Insertion)
            // and 'D' = deletion from query (SAM Deletion).
            let status = wavefront_align(
                self.wf_aligner,
                reference.as_ptr() as *const c_char,
                reference.len() as i32,
                query.as_ptr() as *const c_char,
                query.len() as i32,
            );

            if status != 0 {
                return Err(Wfa2Error::Failed(status));
            }

            // Extract CIGAR from the aligner
            let cigar_ptr = (*self.wf_aligner).cigar;
            if cigar_ptr.is_null() {
                return Err(Wfa2Error::Failed(-999));
            }

            let cigar = &*cigar_ptr;
            let begin = cigar.begin_offset as usize;
            let end = cigar.end_offset as usize;

            if end <= begin {
                return Ok(Alignment {
                    divergence: DivergenceScore::ZERO,
                    cigar: Vec::new(),
                });
            }

            let ops_raw = std::slice::from_raw_parts(
                cigar.operations.add(begin) as *const u8,
                end - begin,
            );

            // Convert WFA2 per-column CIGAR characters to run-length encoded ops.
            //
            // WFA2 uses: 'M' = match/mismatch, 'X' = mismatch, 'I' = insertion,
            //            'D' = deletion. But in gap_affine_2p mode it actually
            //            emits 'M' for matches and 'X' for mismatches.
            let cigar_ops = rle_encode_wfa2_cigar(ops_raw, query, reference);

            // Compute divergence consistent with block aligner:
            //   score = sum(match_score=2 for =, -mismatch=4 for X, gap penalties for I/D)
            //   divergence = (query_len - score).max(0)
            // Use single-piece affine (block aligner's convention) for divergence,
            // since the number is only used for diagnostics/tests.
            let score = compute_block_aligner_score(ops_raw, query, reference);
            let edit_score = (query.len() as i32 - score).max(0);

            Ok(Alignment {
                divergence: DivergenceScore::new(edit_score as f64),
                cigar: cigar_ops,
            })
        }
    }
}

/// Convert WFA2's per-column CIGAR characters into run-length encoded
/// noodles `Op` values with extended CIGAR (=, X, I, D).
///
/// WFA2 emits one character per alignment column:
/// - 'M': match (could be = or X depending on actual bases)
/// - 'X': mismatch
/// - 'I': insertion (query base, no ref base consumed)
/// - 'D': deletion (ref base, no query base consumed)
///
/// We convert 'M' by checking actual bases to produce = or X.
fn rle_encode_wfa2_cigar(ops: &[u8], query: &[u8], reference: &[u8]) -> Vec<Op> {
    let mut result = Vec::new();
    let mut q_pos = 0usize;
    let mut r_pos = 0usize;

    let mut current_kind: Option<Kind> = None;
    let mut current_len = 0usize;

    for &op_byte in ops {
        let kind = match op_byte {
            b'M' => {
                // Resolve M to = or X by comparing actual bases
                let is_match = q_pos < query.len()
                    && r_pos < reference.len()
                    && query[q_pos].to_ascii_uppercase()
                        == reference[r_pos].to_ascii_uppercase();
                if is_match {
                    Kind::SequenceMatch
                } else {
                    Kind::SequenceMismatch
                }
            }
            b'X' => Kind::SequenceMismatch,
            b'I' => Kind::Insertion,
            b'D' => Kind::Deletion,
            _ => {
                // Unknown op — skip
                continue;
            }
        };

        // Advance sequence positions
        match op_byte {
            b'M' | b'X' => {
                q_pos += 1;
                r_pos += 1;
            }
            b'I' => {
                q_pos += 1;
            }
            b'D' => {
                r_pos += 1;
            }
            _ => {}
        }

        // RLE: extend or start new run
        if current_kind == Some(kind) {
            current_len += 1;
        } else {
            if let Some(k) = current_kind {
                result.push(Op::new(k, current_len));
            }
            current_kind = Some(kind);
            current_len = 1;
        }
    }

    // Flush last run
    if let Some(k) = current_kind {
        result.push(Op::new(k, current_len));
    }

    result
}

/// Compute a score consistent with the block aligner's NW1 scoring matrix:
///   match = +2, mismatch = -4, gap_open = -6, gap_extend = -2
/// This is used to produce a `divergence` value compatible with existing tests.
fn compute_block_aligner_score(ops: &[u8], query: &[u8], reference: &[u8]) -> i32 {
    const MATCH_SCORE: i32 = 2;
    const MISMATCH_PENALTY: i32 = -4;
    const GAP_OPEN: i32 = -6;
    const GAP_EXTEND: i32 = -2;

    let mut score = 0i32;
    let mut q_pos = 0usize;
    let mut r_pos = 0usize;
    let mut in_gap = false; // Track whether previous op was a gap
    let mut prev_gap_type: u8 = 0; // 'I' or 'D'

    for &op in ops {
        match op {
            b'M' => {
                let is_match = q_pos < query.len()
                    && r_pos < reference.len()
                    && query[q_pos].to_ascii_uppercase()
                        == reference[r_pos].to_ascii_uppercase();
                score += if is_match { MATCH_SCORE } else { MISMATCH_PENALTY };
                q_pos += 1;
                r_pos += 1;
                in_gap = false;
            }
            b'X' => {
                score += MISMATCH_PENALTY;
                q_pos += 1;
                r_pos += 1;
                in_gap = false;
            }
            b'I' => {
                if in_gap && prev_gap_type == b'I' {
                    score += GAP_EXTEND;
                } else {
                    score += GAP_OPEN + GAP_EXTEND;
                }
                q_pos += 1;
                in_gap = true;
                prev_gap_type = b'I';
            }
            b'D' => {
                if in_gap && prev_gap_type == b'D' {
                    score += GAP_EXTEND;
                } else {
                    score += GAP_OPEN + GAP_EXTEND;
                }
                r_pos += 1;
                in_gap = true;
                prev_gap_type = b'D';
            }
            _ => {}
        }
    }

    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wfa2_perfect_match() {
        let mut aligner = Wfa2Aligner::with_defaults();
        let result = aligner.align(b"ACGTACGT", b"ACGTACGT").unwrap();
        assert_eq!(result.cigar.len(), 1);
        assert_eq!(result.cigar[0].kind(), Kind::SequenceMatch);
        assert_eq!(result.cigar[0].len(), 8);
    }

    #[test]
    fn test_wfa2_single_mismatch() {
        let mut aligner = Wfa2Aligner::with_defaults();
        let result = aligner.align(b"ACGTACGT", b"ACATACGT").unwrap();
        let cigar_str = result.cigar_string();
        assert!(
            cigar_str.contains('X'),
            "Expected mismatch in CIGAR: {}",
            cigar_str
        );
    }

    #[test]
    fn test_wfa2_insertion() {
        let mut aligner = Wfa2Aligner::with_defaults();
        // Query has extra bases
        let result = aligner.align(b"ACGTTTACGT", b"ACGTACGT").unwrap();
        let cigar_str = result.cigar_string();
        assert!(
            cigar_str.contains('I'),
            "Expected insertion in CIGAR: {}",
            cigar_str
        );
        assert_eq!(result.query_length(), 10);
        assert_eq!(result.reference_consumed(), 8);
    }

    #[test]
    fn test_wfa2_deletion() {
        let mut aligner = Wfa2Aligner::with_defaults();
        // Reference has extra bases
        let result = aligner.align(b"ACGTACGT", b"ACGTTTACGT").unwrap();
        let cigar_str = result.cigar_string();
        assert!(
            cigar_str.contains('D'),
            "Expected deletion in CIGAR: {}",
            cigar_str
        );
        assert_eq!(result.query_length(), 8);
        assert_eq!(result.reference_consumed(), 10);
    }

    #[test]
    fn test_wfa2_empty_query() {
        let mut aligner = Wfa2Aligner::with_defaults();
        let result = aligner.align(b"", b"ACGT").unwrap();
        assert_eq!(result.cigar.len(), 1);
        assert_eq!(result.cigar[0].kind(), Kind::Deletion);
        assert_eq!(result.cigar[0].len(), 4);
    }

    #[test]
    fn test_wfa2_empty_reference() {
        let mut aligner = Wfa2Aligner::with_defaults();
        let result = aligner.align(b"ACGT", b"").unwrap();
        assert_eq!(result.cigar.len(), 1);
        assert_eq!(result.cigar[0].kind(), Kind::Insertion);
        assert_eq!(result.cigar[0].len(), 4);
    }

    #[test]
    fn test_wfa2_both_empty() {
        let mut aligner = Wfa2Aligner::with_defaults();
        let result = aligner.align(b"", b"").unwrap();
        assert!(result.cigar.is_empty());
    }

    #[test]
    fn test_wfa2_large_insertion_favoured() {
        // This tests the two-piece model: a read with ~40bp extra should
        // produce a single insertion rather than fragmenting it into small
        // indels scattered across matching regions.
        let mut aligner = Wfa2Aligner::with_defaults();

        // Reference: 50bp flanking sequence on each side
        let flank_l = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        let flank_r = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";

        // Query: same flanks with 40bp insert
        let insert = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; // 42bp
        let mut query = Vec::new();
        query.extend_from_slice(flank_l);
        query.extend_from_slice(insert);
        query.extend_from_slice(flank_r);

        let mut reference = Vec::new();
        reference.extend_from_slice(flank_l);
        reference.extend_from_slice(flank_r);

        let result = aligner.align(&query, &reference).unwrap();

        // Count insertions: should be one big insertion, not many small ones
        let insertion_ops: Vec<_> = result
            .cigar
            .iter()
            .filter(|op| op.kind() == Kind::Insertion)
            .collect();

        assert!(
            insertion_ops.len() <= 2,
            "Expected at most 2 insertion ops (consolidation), got {} in CIGAR: {}",
            insertion_ops.len(),
            result.cigar_string()
        );

        let total_ins: usize = insertion_ops.iter().map(|op| op.len()).sum();
        assert_eq!(
            total_ins,
            insert.len(),
            "Total insertion length should match insert size"
        );
    }
}
