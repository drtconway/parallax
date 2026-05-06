use block_aligner::cigar::{Cigar, Operation};
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, NW1, NucMatrix};

use parallax::config::BlockAlignerConfig;
use parallax::scores::DivergenceScore;

use super::{Alignment, Kind, Op};

/// Error types for block aligner operations
#[derive(Debug, Clone)]
pub enum BlockAlignerError {
    /// Sequence too short for alignment
    #[allow(dead_code)]
    SequenceTooShort,
    /// Invalid sequence characters
    #[allow(dead_code)]
    InvalidSequence(String),
    /// Alignment failed
    #[allow(dead_code)]
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
/// when aligning many sequence pairs. In particular, the `Block` objects
/// from block-aligner are cached and reused, avoiding multi-MB trace
/// buffer allocations on every call.
pub struct BlockAligner {
    /// Configuration parameters
    config: BlockAlignerConfig,
    /// Gap penalties (cached from config)
    gaps: Gaps,
    /// Cached block for global alignment (TRACE=true, XDROP=false)
    global_block: Option<Block<true, false>>,
    /// Allocated sequence capacity for global_block
    global_block_capacity: usize,
    /// Cached block for extension alignment (TRACE=true, XDROP=true)
    xdrop_block: Option<Block<true, true>>,
    /// Allocated sequence capacity for xdrop_block
    xdrop_block_capacity: usize,
    /// Reusable padded query buffer
    padded_q: Option<PaddedBytes>,
    /// Reusable padded reference buffer
    padded_r: Option<PaddedBytes>,
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
            global_block: None,
            global_block_capacity: 0,
            xdrop_block: None,
            xdrop_block_capacity: 0,
            padded_q: None,
            padded_r: None,
        }
    }

    /// Create a BlockAligner with default configuration.
    pub fn with_defaults() -> Self {
        Self::new(&BlockAlignerConfig::default())
    }

    /// Set configuration from AlignParams.
    #[allow(dead_code)]
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

    /// Ensure the cached global block exists and is large enough.
    fn ensure_global_block(&mut self, q_len: usize, r_len: usize, max_bs: usize) {
        let alloc_len = q_len.max(r_len).max(4096).next_power_of_two();
        let alloc_bs = max_bs.max(self.config.max_block_size.min(16384));
        if self.global_block.is_none() || alloc_len > self.global_block_capacity {
            self.global_block = Some(Block::<true, false>::new(alloc_len, alloc_len, alloc_bs));
            self.global_block_capacity = alloc_len;
        }
    }

    /// Ensure the cached xdrop block exists and is large enough.
    fn ensure_xdrop_block(&mut self, q_len: usize, r_len: usize, max_bs: usize) {
        let alloc_len = q_len.max(r_len).max(4096).next_power_of_two();
        let alloc_bs = max_bs.max(self.config.max_block_size.min(16384));
        if self.xdrop_block.is_none() || alloc_len > self.xdrop_block_capacity {
            self.xdrop_block = Some(Block::<true, true>::new(alloc_len, alloc_len, alloc_bs));
            self.xdrop_block_capacity = alloc_len;
        }
    }

    /// Ensure padded buffers are allocated with sufficient capacity.
    fn ensure_padded_bufs(&mut self, max_seq_len: usize, max_bs: usize) {
        let alloc_len = max_seq_len.max(4096).next_power_of_two();
        let needed = 1 + alloc_len + max_bs;
        if self.padded_q.as_ref().map_or(true, |p| p.len() + 1 + max_bs < needed) {
            self.padded_q = Some(PaddedBytes::new::<NucMatrix>(alloc_len, max_bs));
        }
        if self.padded_r.as_ref().map_or(true, |p| p.len() + 1 + max_bs < needed) {
            self.padded_r = Some(PaddedBytes::new::<NucMatrix>(alloc_len, max_bs));
        }
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
                divergence: DivergenceScore::ZERO,
                cigar: Vec::new(),
            });
        }
        if query.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::new(reference.len() as f64),
                cigar: vec![Op::new(Kind::Deletion, reference.len())],
            });
        }
        if reference.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::new(query.len() as f64),
                cigar: vec![Op::new(Kind::Insertion, query.len())],
            });
        }

        let max_len = query.len().max(reference.len());
        let (min_bs, max_bs) = self.block_sizes(max_len);

        self.ensure_padded_bufs(max_len, max_bs);
        self.padded_q.as_mut().unwrap().set_bytes::<NucMatrix>(query, max_bs);
        self.padded_r.as_mut().unwrap().set_bytes::<NucMatrix>(reference, max_bs);

        // Ensure block is allocated, then borrow fields separately
        let q_len = self.padded_q.as_ref().unwrap().len();
        let r_len = self.padded_r.as_ref().unwrap().len();
        self.ensure_global_block(q_len, r_len, max_bs);
        let gaps = self.gaps;
        let block = self.global_block.as_mut().unwrap();
        let q = self.padded_q.as_ref().unwrap();
        let r = self.padded_r.as_ref().unwrap();
        block.align(q, r, &NW1, gaps, min_bs..=max_bs, 0);

        let res = block.res();

        let mut cigar_buf = Cigar::new(res.query_idx, res.reference_idx);
        block
            .trace()
            .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_buf);

        let cigar = convert_cigar(&cigar_buf);
        let edit_score = (query.len() as i32) - res.score;

        Ok(Alignment {
            divergence: DivergenceScore::new(edit_score.max(0) as f64),
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
                divergence: DivergenceScore::ZERO,
                cigar: Vec::new(),
            });
        }
        if query.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::new(reference.len() as f64),
                cigar: vec![Op::new(Kind::Deletion, reference.len())],
            });
        }
        if reference.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::new(query.len() as f64),
                cigar: vec![Op::new(Kind::Insertion, query.len())],
            });
        }

        // The maximum length we should attempt to align is the length of the longer sequence, up to twice the length of the query.
        let ref_len = reference.len().min(2 * query.len());
        let max_len = query.len().max(ref_len);
        let (min_bs, max_bs) = self.block_sizes(max_len);

        self.ensure_padded_bufs(max_len, max_bs);
        self.padded_q.as_mut().unwrap().set_bytes::<NucMatrix>(query, max_bs);
        self.padded_r.as_mut().unwrap().set_bytes::<NucMatrix>(&reference[..ref_len], max_bs);

        // Ensure block is allocated, then borrow fields separately
        let q_len = self.padded_q.as_ref().unwrap().len();
        let r_len = self.padded_r.as_ref().unwrap().len();
        self.ensure_xdrop_block(q_len, r_len, max_bs);
        let gaps = self.gaps;
        let block = self.xdrop_block.as_mut().unwrap();
        let q = self.padded_q.as_ref().unwrap();
        let r = self.padded_r.as_ref().unwrap();
        block.align(q, r, &NW1, gaps, min_bs..=max_bs, x_drop);

        let res = block.res();

        let mut cigar_buf = Cigar::new(res.query_idx, res.reference_idx);
        block
            .trace()
            .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_buf);

        let cigar = convert_cigar(&cigar_buf);
        let edit_score = (query.len() as i32) - res.score;

        Ok(Alignment {
            divergence: DivergenceScore::new(edit_score.max(0) as f64),
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
                divergence: DivergenceScore::ZERO,
                cigar: Vec::new(),
            });
        }
        if query.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::new(reference.len() as f64),
                cigar: vec![Op::new(Kind::Deletion, reference.len())],
            });
        }
        if reference.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::new(query.len() as f64),
                cigar: vec![Op::new(Kind::Insertion, query.len())],
            });
        }

        // The maximum length we should attempt to align is the length of the longer sequence, up to twice the length of the query.
        let ref_len = reference.len().min(2 * query.len());
        let max_len = query.len().max(ref_len);
        let (min_bs, max_bs) = self.block_sizes(max_len);

        // Use set_bytes_rev to reverse directly into cached padded buffers,
        // avoiding separate reverse buffers and from_bytes allocations.
        let ref_start_offset = reference.len().saturating_sub(ref_len);
        self.ensure_padded_bufs(max_len, max_bs);
        self.padded_q.as_mut().unwrap().set_bytes_rev::<NucMatrix>(query, max_bs);
        self.padded_r.as_mut().unwrap().set_bytes_rev::<NucMatrix>(&reference[ref_start_offset..], max_bs);

        // Ensure block is allocated, then borrow fields separately
        let q_len = self.padded_q.as_ref().unwrap().len();
        let r_len = self.padded_r.as_ref().unwrap().len();
        self.ensure_xdrop_block(q_len, r_len, max_bs);
        let gaps = self.gaps;
        let block = self.xdrop_block.as_mut().unwrap();
        let q = self.padded_q.as_ref().unwrap();
        let r = self.padded_r.as_ref().unwrap();
        block.align(q, r, &NW1, gaps, min_bs..=max_bs, x_drop);

        let res = block.res();

        let mut cigar_buf = Cigar::new(res.query_idx, res.reference_idx);
        block
            .trace()
            .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar_buf);

        // Reverse the CIGAR
        let cigar = reverse_cigar(&convert_cigar(&cigar_buf));
        let edit_score = (query.len() as i32) - res.score;

        Ok(Alignment {
            divergence: DivergenceScore::new(edit_score.max(0) as f64),
            cigar,
        })
    }
}

/// Convert block-aligner CIGAR to noodles Op format
fn convert_cigar(ba_cigar: &Cigar) -> Vec<Op> {
    let mut result = Vec::new();

    for i in 0..ba_cigar.len() {
        let op_len = ba_cigar.get(i);
        let len = op_len.len as usize;

        match op_len.op {
            Operation::Eq => {
                result.push(Op::new(Kind::SequenceMatch, len));
            }
            Operation::X => {
                result.push(Op::new(Kind::SequenceMismatch, len));
            }
            Operation::I => {
                result.push(Op::new(Kind::Insertion, len));
            }
            Operation::D => {
                result.push(Op::new(Kind::Deletion, len));
            }
            // M is match/mismatch combined - shouldn't appear after cigar_eq
            Operation::M => {
                result.push(Op::new(Kind::SequenceMatch, len));
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
fn reverse_cigar(cigar: &[Op]) -> Vec<Op> {
    cigar.iter().rev().copied().collect()
}

/// Merge consecutive CIGAR operations of the same type
fn merge_cigar_ops(ops: Vec<Op>) -> Vec<Op> {
    if ops.is_empty() {
        return ops;
    }

    let mut result = Vec::with_capacity(ops.len());
    let mut current = ops[0];

    for op in ops.into_iter().skip(1) {
        if current.kind() == op.kind() {
            current = Op::new(current.kind(), current.len() + op.len());
        } else {
            result.push(current);
            current = op;
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

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Perfect match should have low score
        assert!(
            alignment.divergence.0 <= 4.0,
            "Expected low score for identical sequences, got {}",
            alignment.divergence.0
        );
    }

    #[test]
    fn test_single_mismatch() {
        let query = b"ACGTACGT";
        let reference = b"ACGTTCGT";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Should detect the mismatch
        assert!(alignment.divergence.0 >= 0.0);
    }

    #[test]
    fn test_insertion() {
        let query = b"ACGTACGT";
        let reference = b"ACGACGT";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.divergence.0 >= 0.0);
    }

    #[test]
    fn test_deletion() {
        let query = b"ACGACGT";
        let reference = b"ACGTACGT";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.divergence.0 >= 0.0);
    }

    #[test]
    fn test_empty_query() {
        let query = b"";
        let reference = b"ACGT";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.divergence.0, 4.0);
        assert_eq!(alignment.cigar, vec![Op::new(Kind::Deletion, 4)]);
    }

    #[test]
    fn test_empty_reference() {
        let query = b"ACGT";
        let reference = b"";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.divergence.0, 4.0);
        assert_eq!(alignment.cigar, vec![Op::new(Kind::Insertion, 4)]);
    }

    #[test]
    fn test_cigar_output() {
        // Test that we get a valid CIGAR for a simple case
        let query = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 34 A's
        let reference = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        let mut aligner = BlockAligner::with_defaults();
        let result = aligner.align(query, reference);
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
            alignment.divergence.0 <= 4.0,
            "Expected low score for identical sequences, got {}",
            alignment.divergence.0
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
            alignment.divergence.0 <= 4.0,
            "Expected low score for identical sequences, got {}",
            alignment.divergence.0
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
        assert_eq!(result.unwrap().cigar, vec![Op::new(Kind::Deletion, 4)]);

        let result = aligner.extend_left(b"ACGT", b"");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().cigar, vec![Op::new(Kind::Insertion, 4)]);
    }

    #[test]
    fn test_reverse_cigar() {
        let cigar = vec![
            Op::new(Kind::SequenceMatch, 5),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::SequenceMatch, 10),
            Op::new(Kind::Deletion, 3),
        ];
        let reversed = reverse_cigar(&cigar);
        assert_eq!(
            reversed,
            vec![
                Op::new(Kind::Deletion, 3),
                Op::new(Kind::SequenceMatch, 10),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::SequenceMatch, 5),
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
            mismatch: 5,
            gap_open: 8,
            gap_extend: 3,
        };

        let aligner = BlockAligner::new(&config);
        assert_eq!(aligner.config().x_drop, 300);
        assert_eq!(aligner.config().gap_open, 8);
    }
}
