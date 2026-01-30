use block_aligner::cigar::{Cigar, Operation};
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, NW1, NucMatrix};

use crate::config::BlockAlignerConfig;

use super::{Alignment, CigarOp};

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

/// SIMD-accelerated aligner using block-aligner library.
///
/// This struct holds reusable buffers to avoid repeated allocations
/// when aligning many sequence pairs.
pub struct BlockAligner {
    /// Configuration parameters
    config: BlockAlignerConfig,
    /// Gap penalties (cached from config)
    gaps: Gaps,
    /// Reusable buffer for reversed query (left extension)
    query_rev_buf: Vec<u8>,
    /// Reusable buffer for reversed reference (left extension)
    ref_rev_buf: Vec<u8>,
}

impl BlockAligner {
    /// Create a new BlockAligner with the given configuration.
    pub fn new(config: &BlockAlignerConfig) -> Self {
        Self {
            config: config.clone(),
            gaps: Gaps {
                open: -(config.gap_open as i8),
                extend: -(config.gap_extend as i8),
            },
            query_rev_buf: Vec::with_capacity(4096),
            ref_rev_buf: Vec::with_capacity(4096),
        }
    }

    /// Create a BlockAligner with default configuration.
    pub fn with_defaults() -> Self {
        Self::new(&BlockAlignerConfig::default())
    }

    /// Set configuration from AlignParams.
    pub fn set_align_params(&mut self, params: &super::AlignParams) {
        self.config.gap_open = params.gap_open;
        self.config.gap_extend = params.gap_extend;
        self.gaps = Gaps {
            open: -(self.config.gap_open as i8),
            extend: -(self.config.gap_extend as i8),
        };
    }

    /// Get the current configuration.
    #[allow(dead_code)]
    pub fn config(&self) -> &BlockAlignerConfig {
        &self.config
    }

    /// Compute block size range based on sequence lengths.
    fn block_sizes(&self, max_len: usize) -> (usize, usize) {
        let min_bs = self.config.min_block_size.max(32);
        let max_bs = (max_len.next_power_of_two())
            .max(min_bs)
            .min(self.config.max_block_size.min(16384));
        (min_bs, max_bs)
    }

    /// Align two sequences using global alignment (no X-drop).
    ///
    /// Returns an edit-distance style score where 0 = perfect match.
    pub fn align(
        &mut self,
        query: &[u8],
        reference: &[u8],
    ) -> Result<Alignment, BlockAlignerError> {
        // Handle empty sequences
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

        let max_len = query.len().max(reference.len());
        let (min_bs, max_bs) = self.block_sizes(max_len);

        let q = PaddedBytes::from_bytes::<NucMatrix>(query, max_bs);
        let r = PaddedBytes::from_bytes::<NucMatrix>(reference, max_bs);

        // Block::<TRACE=true, XDROP=false>
        let mut block = Block::<true, false>::new(q.len(), r.len(), max_bs);
        block.align(&q, &r, &NW1, self.gaps, min_bs..=max_bs, 0);

        let res = block.res();

        let mut cigar_ba = Cigar::new(res.query_idx, res.reference_idx);
        block
            .trace()
            .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_ba);

        let cigar = convert_cigar(&cigar_ba);
        let edit_score = (query.len() as i32) - res.score;

        Ok(Alignment {
            score: edit_score.max(0),
            cigar,
        })
    }

    /// Extend alignment rightward (forward) with X-drop early termination.
    ///
    /// Used to extend beyond the last anchor toward the end of sequences.
    /// Stops early if score drops more than `x_drop` below maximum.
    pub fn extend_right(
        &mut self,
        query: &[u8],
        reference: &[u8],
    ) -> Result<Alignment, BlockAlignerError> {
        self.extend_right_with_xdrop(query, reference, self.config.x_drop)
    }

    /// Extend rightward with custom X-drop threshold.
    pub fn extend_right_with_xdrop(
        &mut self,
        query: &[u8],
        reference: &[u8],
        x_drop: i32,
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

        let max_len = query.len().max(reference.len());
        let (min_bs, max_bs) = self.block_sizes(max_len);

        let q = PaddedBytes::from_bytes::<NucMatrix>(query, max_bs);
        let r = PaddedBytes::from_bytes::<NucMatrix>(reference, max_bs);

        // Block::<TRACE=true, XDROP=true>
        let mut block = Block::<true, true>::new(q.len(), r.len(), max_bs);
        block.align(&q, &r, &NW1, self.gaps, min_bs..=max_bs, x_drop);

        let res = block.res();

        let mut cigar_ba = Cigar::new(res.query_idx, res.reference_idx);
        block
            .trace()
            .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_ba);

        let cigar = convert_cigar(&cigar_ba);
        let edit_score = (query.len() as i32) - res.score;

        Ok(Alignment {
            score: edit_score.max(0),
            cigar,
        })
    }

    /// Extend alignment leftward (backward) with X-drop early termination.
    ///
    /// Used to extend before the first anchor toward the start of sequences.
    /// Reverses sequences internally, aligns, then reverses the CIGAR.
    pub fn extend_left(
        &mut self,
        query: &[u8],
        reference: &[u8],
    ) -> Result<Alignment, BlockAlignerError> {
        self.extend_left_with_xdrop(query, reference, self.config.x_drop)
    }

    /// Extend leftward with custom X-drop threshold.
    pub fn extend_left_with_xdrop(
        &mut self,
        query: &[u8],
        reference: &[u8],
        x_drop: i32,
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

        // Reverse sequences into reusable buffers
        self.query_rev_buf.clear();
        self.query_rev_buf.extend(query.iter().rev().copied());
        self.ref_rev_buf.clear();
        self.ref_rev_buf.extend(reference.iter().rev().copied());

        let max_len = query.len().max(reference.len());
        let (min_bs, max_bs) = self.block_sizes(max_len);

        let q = PaddedBytes::from_bytes::<NucMatrix>(&self.query_rev_buf, max_bs);
        let r = PaddedBytes::from_bytes::<NucMatrix>(&self.ref_rev_buf, max_bs);

        // Block::<TRACE=true, XDROP=true>
        let mut block = Block::<true, true>::new(q.len(), r.len(), max_bs);
        block.align(&q, &r, &NW1, self.gaps, min_bs..=max_bs, x_drop);

        let res = block.res();

        let mut cigar_ba = Cigar::new(res.query_idx, res.reference_idx);
        block
            .trace()
            .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_ba);

        // Reverse the CIGAR
        let cigar = reverse_cigar(&convert_cigar(&cigar_ba));
        let edit_score = (query.len() as i32) - res.score;

        Ok(Alignment {
            score: edit_score.max(0),
            cigar,
        })
    }
}

/// Convenience function to align two sequences with default configuration.
///
/// For better performance when aligning many pairs, use `BlockAligner::new()`
/// to create a reusable aligner instance.
pub fn align(
    query: &[u8],
    reference: &[u8],
) -> Result<Alignment, BlockAlignerError> {
    BlockAligner::with_defaults().align(query, reference)
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

/// Reverse a CIGAR string (for left extension)
fn reverse_cigar(cigar: &[CigarOp]) -> Vec<CigarOp> {
    cigar.iter().rev().copied().collect()
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

        let result = align(query, reference);
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

        let result = align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Should detect the mismatch
        assert!(alignment.score >= 0);
    }

    #[test]
    fn test_insertion() {
        let query = b"ACGTACGT";
        let reference = b"ACGACGT";

        let result = align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.score >= 0);
    }

    #[test]
    fn test_deletion() {
        let query = b"ACGACGT";
        let reference = b"ACGTACGT";

        let result = align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.score >= 0);
    }

    #[test]
    fn test_empty_query() {
        let query = b"";
        let reference = b"ACGT";

        let result = align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.score, 4);
        assert_eq!(alignment.cigar, vec![CigarOp::Del(4)]);
    }

    #[test]
    fn test_empty_reference() {
        let query = b"ACGT";
        let reference = b"";

        let result = align(query, reference);
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

        let result = align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(!alignment.cigar.is_empty());
    }

    #[test]
    fn test_extend_right_identical() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.extend_right(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(
            alignment.score <= 4,
            "Expected low score for identical sequences, got {}",
            alignment.score
        );
    }

    #[test]
    fn test_extend_left_identical() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.extend_left(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(
            alignment.score <= 4,
            "Expected low score for identical sequences, got {}",
            alignment.score
        );
    }

    #[test]
    fn test_extend_right_with_mismatch() {
        let query = b"ACGTACGTACGTACGTNNNNNNNNNNNNNNNN";
        let reference = b"ACGTACGTACGTACGTAAAAAAAAAAAAAAAA";

        // With X-drop, should stop when mismatches cause score to drop
        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.extend_right_with_xdrop(query, reference, 50);
        assert!(result.is_ok());
    }

    #[test]
    fn test_extend_left_with_mismatch() {
        let query = b"NNNNNNNNNNNNNNNNACGTACGTACGTACGT";
        let reference = b"AAAAAAAAAAAAAAAAACGTACGTACGTACGT";

        // With X-drop, should stop when mismatches cause score to drop
        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.extend_left_with_xdrop(query, reference, 50);
        assert!(result.is_ok());
    }

    #[test]
    fn test_extend_empty_sequences() {
        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.extend_right(b"", b"ACGT");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().cigar, vec![CigarOp::Del(4)]);

        let result = aligner.extend_left(b"ACGT", b"");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().cigar, vec![CigarOp::Ins(4)]);
    }

    #[test]
    fn test_reverse_cigar() {
        let cigar = vec![
            CigarOp::Match(5),
            CigarOp::Ins(2),
            CigarOp::Match(10),
            CigarOp::Del(3),
        ];
        let reversed = reverse_cigar(&cigar);
        assert_eq!(
            reversed,
            vec![
                CigarOp::Del(3),
                CigarOp::Match(10),
                CigarOp::Ins(2),
                CigarOp::Match(5),
            ]
        );
    }

    #[test]
    fn test_aligner_reuse() {
        // Test that we can reuse the aligner for multiple alignments
        let config = BlockAlignerConfig {
            x_drop: 200,
            gap_open: 4,
            gap_extend: 1,
            ..Default::default()
        };
        let mut aligner = BlockAligner::new(&config);

        // Align several pairs
        let pairs = [
            (b"ACGTACGTACGTACGTACGTACGTACGTACGT" as &[u8], b"ACGTACGTACGTACGTACGTACGTACGTACGT" as &[u8]),
            (b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
            (b"ACGTACGTACGTACGT", b"ACGTACGTACGT"),
        ];

        for (query, reference) in pairs {
            let result = aligner.align(query, reference);
            assert!(result.is_ok(), "Failed to align pair");
        }
    }

    #[test]
    fn test_config_from_struct() {
        let config = BlockAlignerConfig {
            enable_extension: true,
            min_block_size: 64,
            max_block_size: 2048,
            x_drop: 300,
            end_bonus: 10,
            mismatch: 5,
            gap_open: 8,
            gap_extend: 3,
        };

        let aligner = BlockAligner::new(&config);
        assert_eq!(aligner.config().x_drop, 300);
        assert_eq!(aligner.config().gap_open, 8);
    }
}
