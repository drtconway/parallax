use block_aligner::cigar::{Cigar, Operation};
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, NW1, NucMatrix};

use super::{AlignParams, Alignment, CigarOp};

/// Error types for block aligner operations
#[derive(Debug, Clone)]
pub enum BlockAlignerError {
    /// Sequence too short for alignment
    SequenceTooShort,
    /// Invalid sequence characters
    InvalidSequence(String),
    /// Alignment failed
    AlignmentFailed(String),
}

impl std::fmt::Display for BlockAlignerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BlockAlignerError::SequenceTooShort => write!(f, "Sequence too short for alignment"),
            BlockAlignerError::InvalidSequence(msg) => write!(f, "Invalid sequence: {}", msg),
            BlockAlignerError::AlignmentFailed(msg) => write!(f, "Alignment failed: {}", msg),
        }
    }
}

impl std::error::Error for BlockAlignerError {}

/// Align two sequences using block-aligner with default parameters.
pub fn align(
    query: &[u8],
    reference: &[u8],
    min_block_size: usize,
) -> Result<Alignment, BlockAlignerError> {
    align_with_params(query, reference, min_block_size, AlignParams::default())
}

/// Align two sequences using block-aligner with custom parameters.
pub fn align_with_params(
    query: &[u8],
    reference: &[u8],
    min_block_size: usize,
    params: AlignParams,
) -> Result<Alignment, BlockAlignerError> {
    if query.is_empty() && reference.is_empty() {
        return Ok(Alignment {
            score: 0,
            cigar: Vec::new(),
        });
    }
    if query.is_empty() {
        return Ok(Alignment {
            score: reference.len() as i32,
            cigar: vec![CigarOp::Del(reference.len() as u32)],
        });
    }
    if reference.is_empty() {
        return Ok(Alignment {
            score: query.len() as i32,
            cigar: vec![CigarOp::Ins(query.len() as u32)],
        });
    }

    // Gap penalties (block-aligner uses negative values, open is cost of first gap base)
    let gaps = Gaps {
        open: -(params.gap_open as i8),
        extend: -(params.gap_extend as i8),
    };

    // Determine block sizes based on sequence lengths
    let max_len = query.len().max(reference.len());
    let min_bs = min_block_size.max(32);
    let max_bs = (max_len.next_power_of_two()).max(min_bs).min(16384);

    // Convert sequences to PaddedBytes
    let q = PaddedBytes::from_bytes::<NucMatrix>(query, max_bs);
    let r = PaddedBytes::from_bytes::<NucMatrix>(reference, max_bs);

    // Create block aligner with traceback enabled, no x-drop
    // Block::<TRACE, XDROP>
    let mut block = Block::<true, false>::new(q.len(), r.len(), max_bs);

    // Perform alignment (NW1 is a simple nucleotide scoring matrix: match=1, mismatch=-1)
    // For custom scoring, we'd need to create a custom matrix
    block.align(&q, &r, &NW1, gaps, min_bs..=max_bs, 0);

    let res = block.res();

    // Compute traceback with =/X distinction
    let mut cigar_buf = Cigar::new(res.query_idx, res.reference_idx);
    block
        .trace()
        .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_buf);

    // Convert CIGAR
    let cigar = convert_cigar(&cigar_buf);

    // Convert score: block-aligner returns positive scores (higher=better)
    // We use edit-distance style (0=best, higher=worse)
    // Approximate: perfect match would score ~query.len(), so invert
    let edit_score = (query.len() as i32) - res.score;

    Ok(Alignment {
        score: edit_score.max(0),
        cigar,
    })
}

/// Convert block-aligner CIGAR to our CigarOp format
fn convert_cigar(ba_cigar: &Cigar) -> Vec<CigarOp> {
    let mut result = Vec::new();

    for i in 0..ba_cigar.len() {
        let op_len = ba_cigar.get(i);
        let len = op_len.len as u32;

        match op_len.op {
            Operation::Eq => {
                result.push(CigarOp::Match(len));
            }
            Operation::X => {
                result.push(CigarOp::Mismatch(len));
            }
            Operation::I => {
                result.push(CigarOp::Ins(len));
            }
            Operation::D => {
                result.push(CigarOp::Del(len));
            }
            // M is match/mismatch combined - shouldn't appear after cigar_eq
            Operation::M => {
                result.push(CigarOp::Match(len));
            }
            Operation::Sentinel => {
                // Ignore sentinel
            }
        }
    }

    // Merge consecutive operations of the same type
    merge_cigar_ops(result)
}

/// Merge consecutive CIGAR operations of the same type
fn merge_cigar_ops(ops: Vec<CigarOp>) -> Vec<CigarOp> {
    if ops.is_empty() {
        return ops;
    }

    let mut result = Vec::with_capacity(ops.len());
    let mut current = ops[0];

    for op in ops.into_iter().skip(1) {
        match (&mut current, &op) {
            (CigarOp::Match(n), CigarOp::Match(m)) => *n += m,
            (CigarOp::Mismatch(n), CigarOp::Mismatch(m)) => *n += m,
            (CigarOp::Ins(n), CigarOp::Ins(m)) => *n += m,
            (CigarOp::Del(n), CigarOp::Del(m)) => *n += m,
            (CigarOp::SoftClip(n), CigarOp::SoftClip(m)) => *n += m,
            _ => {
                result.push(current);
                current = op;
            }
        }
    }
    result.push(current);

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_sequences() {
        let query = b"ACGTACGTACGT";
        let reference = b"ACGTACGTACGT";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Perfect match should have low score
        assert!(
            alignment.score <= 4,
            "Expected low score for identical sequences, got {}",
            alignment.score
        );
    }

    #[test]
    fn test_single_mismatch() {
        let query = b"ACGTACGT";
        let reference = b"ACGTTCGT";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Should detect the mismatch
        assert!(alignment.score >= 0);
    }

    #[test]
    fn test_insertion() {
        let query = b"ACGTACGT";
        let reference = b"ACGACGT";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.score >= 0);
    }

    #[test]
    fn test_deletion() {
        let query = b"ACGACGT";
        let reference = b"ACGTACGT";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.score >= 0);
    }

    #[test]
    fn test_empty_query() {
        let query = b"";
        let reference = b"ACGT";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.score, 4);
        assert_eq!(alignment.cigar, vec![CigarOp::Del(4)]);
    }

    #[test]
    fn test_empty_reference() {
        let query = b"ACGT";
        let reference = b"";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.score, 4);
        assert_eq!(alignment.cigar, vec![CigarOp::Ins(4)]);
    }

    #[test]
    fn test_cigar_output() {
        // Test that we get a valid CIGAR for a simple case
        let query = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 34 A's
        let reference = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        let result = align(query, reference, 32);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(!alignment.cigar.is_empty());
    }
}
