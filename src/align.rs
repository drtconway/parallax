use std::collections::HashMap;
use std::sync::OnceLock;

pub use noodles::sam::alignment::record::cigar::Op;
pub use noodles::sam::alignment::record::cigar::op::Kind;

use parallax::config;
use parallax::scores::{DivergenceScore, QualityScore};
use parallax::utils::telemetry::RecorderExt;
use parallax::utils::telemetry::histogram::HistogramRecorder;
use parallax::utils::telemetry::summary::SimpleSummaryRecorder;

pub mod block;
#[cfg(feature = "wfa2")]
pub mod wfa2;

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
    /// Gap open penalty (positive value) — short gaps (piece 1)
    pub gap_open: i32,
    /// Gap extend penalty (positive value) — short gaps (piece 1)
    pub gap_extend: i32,
    /// Gap open penalty — long gaps (piece 2, two-piece affine)
    pub gap_open2: i32,
    /// Gap extend penalty — long gaps (piece 2, two-piece affine)
    pub gap_extend2: i32,
}

impl Default for AlignParams {
    fn default() -> Self {
        let cfg = config::get();
        Self {
            match_score: cfg.alignment.match_score,
            mismatch: cfg.alignment.mismatch,
            gap_open: cfg.alignment.gap_open,
            gap_extend: cfg.alignment.gap_extend,
            gap_open2: cfg.alignment.gap_open2,
            gap_extend2: cfg.alignment.gap_extend2,
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
            Kind::Insertion | Kind::Deletion => {
                // Two-piece affine: use the cheaper of the two pieces
                let piece1 = self.gap_open as f64 + self.gap_extend as f64 * n;
                let piece2 = self.gap_open2 as f64 + self.gap_extend2 as f64 * n;
                piece1.min(piece2)
            }
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
            Kind::Insertion | Kind::Deletion => {
                // Two-piece affine: cost = max(piece1, piece2)
                // piece1 favours short gaps, piece2 favours long gaps
                let piece1 = -self.gap_open - self.gap_extend * n;
                let piece2 = -self.gap_open2 - self.gap_extend2 * n;
                piece1.max(piece2)
            }
            Kind::SoftClip | Kind::HardClip => 0,
            _ => 0,
        } as f64)
            .into()
    }
}

/// High-level aligner abstraction that wraps the underlying alignment engine.
///
/// Uses WFA2 with two-piece affine gap penalties as the primary alignment
/// engine for gap filling, falling back to block-aligner if WFA2 fails.
/// Extension alignment (left/right with X-drop) still uses block-aligner.
///
/// The aligner maintains a cache of small alignments to speed up repeated calls with the same inputs, which can happen during gap filling.
pub struct DpAligner {
    #[cfg(feature = "wfa2")]
    wfa2: wfa2::Wfa2Aligner,
    inner: block::BlockAligner,
    pub indel_shifter: IndelShifter,
    cache: HashMap<u64, std::result::Result<Alignment, AlignmentError>>,
}

impl DpAligner {
    /// Create a new Aligner from explicit configuration.
    pub fn from_config(
        align_cfg: &crate::config::AlignmentConfig,
        block_cfg: &crate::config::BlockAlignerConfig,
    ) -> Self {
        let _ = align_cfg;
        Self {
            #[cfg(feature = "wfa2")]
            wfa2: wfa2::Wfa2Aligner::new(&wfa2::Wfa2Config {
                mismatch: align_cfg.mismatch,
                gap_open1: align_cfg.gap_open,
                gap_extend1: align_cfg.gap_extend,
                gap_open2: align_cfg.gap_open2,
                gap_extend2: align_cfg.gap_extend2,
            }),
            inner: block::BlockAligner::new(block_cfg),
            indel_shifter: IndelShifter::new(),
            cache: HashMap::new(),
        }
    }

    /// Create an Aligner with explicit configuration.
    #[allow(dead_code)]
    pub fn with_config(config: &crate::config::BlockAlignerConfig) -> Self {
        Self {
            #[cfg(feature = "wfa2")]
            wfa2: wfa2::Wfa2Aligner::with_defaults(),
            inner: block::BlockAligner::new(config),
            indel_shifter: IndelShifter::new(),
            cache: HashMap::new(),
        }
    }

    /// Create an Aligner with default configuration (no global config required).
    ///
    /// Useful for tests where the global config may not be initialized.
    pub fn with_defaults() -> Self {
        Self {
            #[cfg(feature = "wfa2")]
            wfa2: wfa2::Wfa2Aligner::with_defaults(),
            inner: block::BlockAligner::with_defaults(),
            indel_shifter: IndelShifter::new(),
            cache: HashMap::new(),
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
        query_length_recorder().record(query.len());
        ref_length_recorder().record(reference.len());

        // Handle empty sequences directly — no aligner required.
        if query.is_empty() && reference.is_empty() {
            return Ok(Alignment {
                divergence: DivergenceScore::ZERO,
                cigar: Vec::new(),
            });
        }
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

        let key = self.cache_key(query, reference);

        if let Some(cached) = key.and_then(|k| self.cache.get(&k)) {
            return cached.clone();
        }

        let start = std::time::Instant::now();

        // Try WFA2 first (two-piece affine gap penalties), if compiled in.
        #[cfg(feature = "wfa2")]
        let result = {
            let wfa2_start = std::time::Instant::now();
            match self.wfa2.align(query, reference) {
                Ok(mut aln) => {
                    aln.normalize();
                    metrics::histogram!("align_wfa2_us")
                        .record(wfa2_start.elapsed().as_micros() as f64);
                    Ok(aln)
                }
                Err(_wfa2_err) => {
                    // WFA2 failed (heuristic cutoff or other error); record wasted WFA2 time
                    // separately so we can see the fallback rate and cost.
                    metrics::histogram!("align_wfa2_failed_us")
                        .record(wfa2_start.elapsed().as_micros() as f64);
                    let fallback_start = std::time::Instant::now();
                    let r = self
                        .inner
                        .align(query, reference)
                        .map_err(AlignmentError::BlockError)
                        .map(|mut aln| {
                            aln.normalize();
                            aln
                        });
                    metrics::histogram!("align_fallback_us")
                        .record(fallback_start.elapsed().as_micros() as f64);
                    r
                }
            }
        };

        // Without wfa2 feature, use block-aligner directly.
        #[cfg(not(feature = "wfa2"))]
        let result = self
            .inner
            .align(query, reference)
            .map_err(AlignmentError::BlockError)
            .map(|mut aln| {
                aln.normalize();
                aln
            });

        if let Some(k) = key {
            self.cache.insert(k, result.clone());
        }

        let elapsed = start.elapsed().as_secs_f64();
        align_time_recorder().record(elapsed);

        result
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

    /// Compute a cache key for the given query and reference sequences.
    ///
    /// Both sequences must be less than 8 bases, and contain only A/C/G/T (case-insensitive).
    /// The key is a 64-bit integer with the following layout:
    ///     bytes   meaning
    ///     -----   -------
    ///     0       query length (must be < 8) and ref length (must be < 8), packed as (ref_len << 4) | query_len
    ///     1-2     query bases, 2 bits per base (A=00, C=01, G=10, T=11)
    ///     3-4     reference bases, same encoding as query
    fn cache_key(&self, query: &[u8], reference: &[u8]) -> Option<u64> {
        if query.len() >= 8 || reference.len() >= 8 {
            return None;
        }

        let mut key = ((reference.len() as u64) << 4) | (query.len() as u64);
        for &b in query {
            key <<= 2;
            key |= match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None, // Invalid base for cache key
            };
        }
        for &b in reference {
            key <<= 2;
            key |= match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None, // Invalid base for cache key
            };
        }
        Some(key)
    }
}

impl Default for DpAligner {
    fn default() -> Self {
        Self::with_defaults()
    }
}

#[derive(Debug, Clone)]
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

/// Reusable workspace for left-aligning indels.
///
/// Holds the scratch buffers (`cols`, `q_pos`, `r_pos`, `new_cigar`) that
/// `left_align_indels` needs, so they can be allocated once and reused
/// across many alignments instead of being freshly allocated each time.
#[derive(Clone, Debug, Default)]
pub struct IndelShifter {
    /// Per-column expanded CIGAR representation (b'=', b'X', b'I', b'D', b'S', b'H').
    cols: Vec<u8>,
    /// Query bases consumed before each column.
    q_pos: Vec<usize>,
    /// Reference bases consumed before each column.
    r_pos: Vec<usize>,
    /// RLE-encoded CIGAR rebuilt after shifting.
    new_cigar: Vec<Op>,
}

impl IndelShifter {
    /// Create a new, empty `IndelShifter`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Left-align all indels in `alignment` through matching/mismatching sequence.
    ///
    /// See [`Alignment::left_align_indels`] for the algorithm description.
    pub fn left_align_indels(&mut self, alignment: &mut Alignment, query: &[u8], reference: &[u8]) {
        if alignment.cigar.is_empty() {
            return;
        }

        // Expand CIGAR to per-column representation.
        self.cols.clear();
        for &op in &alignment.cigar {
            let ch = match op.kind() {
                Kind::SequenceMatch => b'=',
                Kind::SequenceMismatch => b'X',
                Kind::Insertion => b'I',
                Kind::Deletion => b'D',
                Kind::SoftClip => b'S',
                Kind::HardClip => b'H',
                _ => continue,
            };
            self.cols.resize(self.cols.len() + op.len(), ch);
        }

        let n = self.cols.len();
        if n == 0 {
            return;
        }

        // Build position arrays:
        //   q_pos[i] = query bases consumed before column i
        //   r_pos[i] = reference bases consumed before column i
        self.q_pos.clear();
        self.q_pos.reserve(n + 1);
        self.q_pos.push(0);
        self.r_pos.clear();
        self.r_pos.reserve(n + 1);
        self.r_pos.push(0);
        let mut q = 0usize;
        let mut r = 0usize;
        for &ch in &self.cols {
            match ch {
                b'D' | b'H' => {}
                _ => q += 1,
            }
            match ch {
                b'I' | b'S' | b'H' => {}
                _ => r += 1,
            }
            self.q_pos.push(q);
            self.r_pos.push(r);
        }

        // Process each indel block left-to-right, shifting it as far left
        // as possible through preceding aligned (= or X) columns.
        let mut any_shifted = false;
        let mut i = 0;
        while i < n {
            if self.cols[i] != b'I' && self.cols[i] != b'D' {
                i += 1;
                continue;
            }

            let indel_type = self.cols[i];
            let start = i;
            while i < n && self.cols[i] == indel_type {
                i += 1;
            }
            let end = i;
            let k = end - start;

            // Compute maximum left shift distance.
            let mut shift = 0usize;
            if indel_type == b'I' {
                let q = self.q_pos[start]; // query pos at first insertion column
                while shift < start
                    && (self.cols[start - 1 - shift] == b'='
                        || self.cols[start - 1 - shift] == b'X')
                    && q >= shift + 1
                    && q + k >= shift + 1
                    && query[q - 1 - shift].to_ascii_uppercase()
                        == query[q + k - 1 - shift].to_ascii_uppercase()
                {
                    shift += 1;
                }
            } else {
                let r = self.r_pos[start]; // ref pos at first deletion column
                while shift < start
                    && (self.cols[start - 1 - shift] == b'='
                        || self.cols[start - 1 - shift] == b'X')
                    && r >= shift + 1
                    && r + k >= shift + 1
                    && reference[r - 1 - shift].to_ascii_uppercase()
                        == reference[r + k - 1 - shift].to_ascii_uppercase()
                {
                    shift += 1;
                }
            }

            if shift > 0 {
                any_shifted = true;
                // Rearrange: move the indel block left by `shift` positions.
                // The `shift` aligned columns that were before the indel
                // move to after it.  Their =/X status is preserved because
                // the shift condition guarantees base equality at each step.
                //
                // rotate_left(shift) on [displaced(shift) | indel(k)]
                // yields [indel(k) | displaced(shift)] — no temp buffer needed.
                let base = start - shift;
                self.cols[base..end].rotate_left(shift);

                // Recompute positions only in the rearranged range.
                // Beyond `end` the same columns exist in the same order,
                // so their prefix sums are unchanged.
                let mut q = self.q_pos[base];
                let mut r = self.r_pos[base];
                for c in base..end {
                    match self.cols[c] {
                        b'D' | b'H' => {}
                        _ => q += 1,
                    }
                    match self.cols[c] {
                        b'I' | b'S' | b'H' => {}
                        _ => r += 1,
                    }
                    self.q_pos[c + 1] = q;
                    self.r_pos[c + 1] = r;
                }
            }
        }

        if !any_shifted {
            return;
        }

        // RLE-encode back to CIGAR ops.
        self.new_cigar.clear();
        let mut j = 0;
        while j < n {
            let kind = match self.cols[j] {
                b'=' => Kind::SequenceMatch,
                b'X' => Kind::SequenceMismatch,
                b'I' => Kind::Insertion,
                b'D' => Kind::Deletion,
                b'S' => Kind::SoftClip,
                b'H' => Kind::HardClip,
                _ => {
                    j += 1;
                    continue;
                }
            };
            let mut len = 1;
            while j + len < n && self.cols[j + len] == self.cols[j] {
                len += 1;
            }
            self.new_cigar.push(Op::new(kind, len));
            j += len;
        }

        std::mem::swap(&mut alignment.cigar, &mut self.new_cigar);
    }
}

/// Alignment result
#[derive(Clone, Default, Debug)]
pub struct Alignment {
    /// Edit distance (lower is better)
    #[allow(dead_code)]
    pub divergence: DivergenceScore,
    /// CIGAR operations
    pub cigar: Vec<Op>,
}

impl Alignment {
    /// Create a new perfect match alignment
    pub fn from_perfect_match(length: usize) -> Self {
        Self {
            divergence: DivergenceScore::ZERO,
            cigar: if length > 0 {
                vec![Op::new(Kind::SequenceMatch, length)]
            } else {
                vec![]
            },
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

    /// Produce a "short version" cigar that compresses all the
    /// matches, mismatches, insertions and deletions.
    pub fn summary_cigar_string(&self) -> String {
        let mut left_clip = 0;
        let mut left_clip_kind = 'S';
        let mut right_clip = 0;
        let mut right_clip_kind = 'S';
        let mut match_count = 0;
        let mut indel_count: isize = 0;
        for (i, op) in self.cigar.iter().enumerate() {
            match op.kind() {
                Kind::Match => {
                    match_count += op.len();
                }
                Kind::Insertion => {
                    indel_count += op.len() as isize;
                }
                Kind::Deletion => {
                    indel_count -= op.len() as isize;
                }
                Kind::Skip => {}
                Kind::SoftClip => {
                    if i == 0 {
                        left_clip = op.len();
                        left_clip_kind = 'S';
                    } else {
                        right_clip = op.len();
                        right_clip_kind = 'S';
                    }
                }
                Kind::HardClip => {
                    if i == 0 {
                        left_clip = op.len();
                        left_clip_kind = 'H';
                    } else {
                        right_clip = op.len();
                        right_clip_kind = 'H';
                    }
                }
                Kind::Pad => {}
                Kind::SequenceMatch => {
                    match_count += op.len();
                }
                Kind::SequenceMismatch => {
                    match_count += op.len();
                }
            }
        }
        
        let mut parts = vec![];
        if left_clip > 0 {
            parts.push(format!("{}{}", left_clip, left_clip_kind));
        }
        if match_count > 0 {
            parts.push(format!("{}M", match_count));
        }
        if indel_count > 0 {
            parts.push(format!("{}I", indel_count));
        } else if indel_count < 0 {
            parts.push(format!("{}D", -indel_count));
        }
        if right_clip > 0 {
            parts.push(format!("{}{}", right_clip, right_clip_kind));
        }
        parts.join("")
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
            if op.len() == 0 {
                continue;
            }
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

    /// Compute the Levenshtein Edit Distance
    pub fn edit_distance(&self) -> usize {
        let mut d = 0;
        for op in self.cigar.iter() {
            if matches!(
                op.kind(),
                Kind::Deletion | Kind::Insertion | Kind::SequenceMismatch
            ) {
                d += op.len();
            }
        }
        d
    }

    /// Left-align all indels through matching/mismatching sequence.
    ///
    /// Convenience wrapper that creates a temporary [`IndelShifter`].
    /// If you are calling this in a loop, prefer creating an `IndelShifter`
    /// once and calling [`IndelShifter::left_align_indels`] to reuse buffers.
    #[allow(dead_code)]
    pub fn left_align_indels(&mut self, query: &[u8], reference: &[u8]) {
        IndelShifter::new().left_align_indels(self, query, reference);
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
                return Err(format!("CIGAR op {} has zero length: {:?}", i, op));
            }
            if i > 0 && self.cigar[i].kind() == self.cigar[i - 1].kind() {
                return Err(format!(
                    "CIGAR ops {} and {} are adjacent and of same type: {:?}, {:?}",
                    i - 1,
                    i,
                    self.cigar[i - 1],
                    self.cigar[i]
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
                                op_idx,
                                n,
                                ref_pos,
                                reference.len(),
                                i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): query_pos {} exceeds query length {} at offset {}",
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
                                op_idx,
                                n,
                                ref_pos,
                                reference.len(),
                                i
                            ));
                        }
                        if query_pos >= query.len() {
                            return Err(format!(
                                "CIGAR op {} ({}{ch}): query_pos {} exceeds query length {} at offset {}",
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
                            op_idx,
                            n,
                            query_pos,
                            n,
                            query.len()
                        ));
                    }
                    query_pos = new_pos;
                }
                Kind::Deletion => {
                    let new_pos = ref_pos + n;
                    if new_pos > reference.len() {
                        return Err(format!(
                            "CIGAR op {} ({}{ch}): ref_pos {} + {} exceeds reference length {}",
                            op_idx,
                            n,
                            ref_pos,
                            n,
                            reference.len()
                        ));
                    }
                    ref_pos = new_pos;
                }
                Kind::SoftClip => {
                    let new_pos = query_pos + n;
                    if new_pos > query.len() {
                        return Err(format!(
                            "CIGAR op {} ({}{ch}): query_pos {} + {} exceeds query length {}",
                            op_idx,
                            n,
                            query_pos,
                            n,
                            query.len()
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

    /// Split the alignment at a ref-space position, returning `(left, right)`.
    ///
    /// `ref_pos` is the number of reference bases that go into the left half.
    /// Any op that straddles the boundary is split by kind — both halves receive
    /// an op of the same kind with their respective lengths.
    /// The divergence field is zeroed on both halves (caller re-scores if needed).
    ///
    /// # Panics
    /// Panics if `ref_pos` exceeds the total reference bases consumed by the CIGAR.
    pub fn split_at_ref_pos(&self, ref_pos: usize) -> (Alignment, Alignment) {
        let mut left: Vec<Op> = Vec::new();
        let mut right: Vec<Op> = Vec::new();
        let mut ref_remaining = ref_pos;

        for &op in &self.cigar {
            if ref_remaining == 0 {
                right.push(op);
                continue;
            }
            let n = op.len();
            if op.kind().consumes_reference() {
                if n <= ref_remaining {
                    left.push(op);
                    ref_remaining -= n;
                } else {
                    left.push(Op::new(op.kind(), ref_remaining));
                    right.push(Op::new(op.kind(), n - ref_remaining));
                    ref_remaining = 0;
                }
            } else {
                // Insertions and other non-ref-consuming ops before the boundary
                // belong to the left side.
                left.push(op);
            }
        }

        assert!(
            ref_remaining == 0,
            "ref_pos {ref_pos} exceeds total CIGAR ref bases (CIGAR: {})",
            self.cigar_string()
        );

        (
            Alignment {
                divergence: DivergenceScore::ZERO,
                cigar: left,
            },
            Alignment {
                divergence: DivergenceScore::ZERO,
                cigar: right,
            },
        )
    }

    /// Compute a quality score from the CIGAR (higher is better).
    pub fn quality(&self, params: &AlignParams) -> QualityScore {
        let mut score = 0.0;
        for &op in &self.cigar {
            score += params.quality(op).0;
        }
        QualityScore::new(score)
    }

    pub fn concat(alignments: &[Alignment]) -> Alignment {
        let mut total_divergence = DivergenceScore::ZERO;
        let mut combined_cigar: Vec<Op> = Vec::new();

        for aln in alignments {
            total_divergence = DivergenceScore::new(total_divergence.0 + aln.divergence.0);
            for &op in &aln.cigar {
                if op.len() == 0 {
                    continue;
                }
                if let Some(last) = combined_cigar.last_mut() {
                    if last.kind() == op.kind() {
                        *last = Op::new(op.kind(), last.len() + op.len());
                        continue;
                    }
                }
                combined_cigar.push(op);
            }
        }

        Alignment {
            divergence: total_divergence,
            cigar: combined_cigar,
        }
    }

    /// Reverse the CIGAR operation order (for reverse-strand assembly where ref and
    /// query run in opposite directions).
    pub fn reversed(&self) -> Self {
        let mut cigar = self.cigar.clone();
        cigar.reverse();
        Alignment { divergence: self.divergence, cigar }
    }

    /// Fraction of aligned bases that are exact matches: matches / max(ref_bases, read_bases).
    /// Returns 1.0 for an empty CIGAR.
    pub fn identity(&self) -> f64 {
        let mut matches = 0usize;
        let mut ref_bases = 0usize;
        let mut read_bases = 0usize;
        for op in &self.cigar {
            match op.kind() {
                Kind::SequenceMatch => {
                    matches += op.len();
                    ref_bases += op.len();
                    read_bases += op.len();
                }
                Kind::SequenceMismatch => {
                    ref_bases += op.len();
                    read_bases += op.len();
                }
                Kind::Deletion => ref_bases += op.len(),
                Kind::Insertion => read_bases += op.len(),
                _ => {}
            }
        }
        let denom = ref_bases.max(read_bases);
        if denom == 0 {
            1.0
        } else {
            matches as f64 / denom as f64
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

    /// Build a summary CIGAR for the SA tag.
    ///
    /// Condenses `=`/`X` into `M`, totals insertions and deletions, and
    /// adds soft-clip ops for the unaligned portions of the read.
    /// For reverse-strand segments the clips are flipped to SAM orientation.
    pub fn summary_cigar(
        &self,
        read_start: usize,
        read_end: usize,
        read_len: usize,
        is_reverse: bool,
    ) -> String {
        let mut matches: usize = 0;
        let mut insertions: usize = 0;
        let mut deletions: usize = 0;
        for op in &self.cigar {
            match op.kind() {
                Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match => {
                    matches += op.len();
                }
                Kind::Insertion => insertions += op.len(),
                Kind::Deletion => deletions += op.len(),
                _ => {}
            }
        }

        // Clips in SAM orientation (rc coordinates for reverse strand).
        let (left_clip, right_clip) = if is_reverse {
            (read_len - read_end, read_start)
        } else {
            (read_start, read_len - read_end)
        };

        let mut s = String::new();
        if left_clip > 0 {
            s.push_str(&format!("{}S", left_clip));
        }
        s.push_str(&format!("{}M", matches));
        if insertions > 0 {
            s.push_str(&format!("{}I", insertions));
        }
        if deletions > 0 {
            s.push_str(&format!("{}D", deletions));
        }
        if right_clip > 0 {
            s.push_str(&format!("{}S", right_clip));
        }
        s
    }
}

impl From<Vec<Op>> for Alignment {
    fn from(cigar: Vec<Op>) -> Self {
        let mut divergence = 0usize;
        for &op in &cigar {
            match op.kind() {
                Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion => divergence += op.len(),
                _ => {}
            }
        }
        Self {
            divergence: DivergenceScore::new(divergence as f64),
            cigar,
        }
    }
}

fn query_length_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("qry_len"))
}

fn ref_length_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("ref_len"))
}

fn align_time_recorder() -> &'static SimpleSummaryRecorder {
    static RECORDER: OnceLock<&'static SimpleSummaryRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| SimpleSummaryRecorder::new_registered("aln_time"))
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
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"ACGT", b"ACGT").unwrap();
        assert_eq!(result.divergence.0, 0.0);
        assert_eq!(result.cigar_string(), "4=");
    }

    #[test]
    fn test_single_mismatch() {
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"ACGT", b"ACTT").unwrap();
        assert_eq!(result.divergence.0, 2.0); // mismatch penalty
        assert_eq!(result.cigar_string(), "2=1X1=");
    }

    #[test]
    fn test_single_insertion() {
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"ACGT", b"ACT").unwrap();
        // query has extra G
        assert!(result.divergence.0 > 0.0);
        assert!(result.cigar_string().contains('I'));
    }

    #[test]
    fn test_single_deletion() {
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"ACT", b"ACGT").unwrap();
        // query missing G
        assert!(result.divergence.0 > 0.0);
        assert!(result.cigar_string().contains('D'));
    }

    #[test]
    fn test_empty() {
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"", b"").unwrap();
        assert_eq!(result.divergence.0, 0.0);
        assert!(result.cigar.is_empty());
    }

    #[test]
    fn test_query_empty() {
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"", b"ACGT").unwrap();
        assert_eq!(result.cigar_string(), "4D");
    }

    #[test]
    fn test_reference_empty() {
        let mut aligner = DpAligner::with_defaults();
        let result = aligner.align(b"ACGT", b"").unwrap();
        assert_eq!(result.cigar_string(), "4I");
    }

    #[test]
    fn test_longer_sequences() {
        let mut aligner = DpAligner::with_defaults();
        let query = b"ACGTACGTACGT";
        let reference = b"ACGTACGTACGT";
        let result = aligner.align(query, reference).unwrap();
        assert_eq!(result.divergence.0, 0.0);
        assert_eq!(result.cigar_string(), "12=");
    }

    #[test]
    fn test_with_gaps() {
        let mut aligner = DpAligner::with_defaults();
        let query = b"ACGTACGT";
        let reference = b"ACGTTTTACGT";
        let result = aligner.align(query, reference).unwrap();
        assert!(result.divergence.0 > 0.0);
        // Should have a deletion in the middle
        println!("CIGAR: {}", result.cigar_string());
    }

    #[test]
    fn test_from_data_1() {
        let mut aligner = DpAligner::with_defaults();
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
        let mut aligner = DpAligner::with_defaults();
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
            println!("{:?} -> {:.2} (running total: {:.2})", op, q, total);
        }
        let alignment = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        assert_eq!(alignment.quality(&params).0, total);
    }

    #[test]
    fn test_left_align_insertion_in_repeat() {
        // In a dinucleotide repeat, an insertion can be placed anywhere.
        // Left-alignment should push it to the leftmost position.
        //
        // Reference: ACACACACAC  (10bp)
        // Query:     ACACACACACAC  (12bp, one extra AC)
        //
        // Starting CIGAR: 6= 2I 4=  (insertion in the middle of the repeat)
        let cigar = parse_cigar("6=2I4=").unwrap();
        let mut aln = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        let query = b"ACACACACACAC";
        let reference = b"ACACACACAC";

        aln.left_align_indels(query, reference);
        let cigar_str = aln.cigar_string();
        // Insertion should be pushed to position 0 (leftmost)
        assert!(
            cigar_str.starts_with("2I"),
            "Expected insertion at left, got: {}",
            cigar_str
        );
        // Total lengths should be preserved
        assert_eq!(aln.query_length(), 12);
        assert_eq!(aln.reference_consumed(), 10);
    }

    #[test]
    fn test_left_align_deletion_in_repeat() {
        // Reference: ACACACACACAC  (12bp)
        // Query:     ACACACACAC    (10bp, missing one AC)
        //
        // Starting CIGAR: 6= 2D 4=  (deletion in the middle)
        let cigar = parse_cigar("6=2D4=").unwrap();
        let mut aln = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        let query = b"ACACACACAC";
        let reference = b"ACACACACACAC";

        aln.left_align_indels(query, reference);
        let cigar_str = aln.cigar_string();
        assert!(
            cigar_str.starts_with("2D"),
            "Expected deletion at left, got: {}",
            cigar_str
        );
        assert_eq!(aln.query_length(), 10);
        assert_eq!(aln.reference_consumed(), 12);
    }

    #[test]
    fn test_left_align_no_shift_when_no_repeat() {
        // Non-repetitive context: insertion cannot shift
        // Reference: ACGTXXXX
        // Query:     ACGTNNXXXX (2bp insertion after ACGT)
        let cigar = parse_cigar("4=2I4=").unwrap();
        let mut aln = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        let query = b"ACGTNNXXXX";
        let reference = b"ACGTXXXX";

        aln.left_align_indels(query, reference);
        assert_eq!(
            aln.cigar_string(),
            "4=2I4=",
            "Should not shift in non-repeat context"
        );
    }

    #[test]
    fn test_left_align_preserves_mismatches() {
        // Insertion adjacent to mismatch, with repeat allowing shift
        // Reference: AAAACCCC
        // Query:     AAAAAACCCC  (2bp AA insertion)
        //
        // CIGAR: 4= 2I 4=  — insertion is between AAAA and CCCC
        let cigar = parse_cigar("4=2I4=").unwrap();
        let mut aln = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        let query = b"AAAAAACCCC";
        let reference = b"AAAACCCC";

        aln.left_align_indels(query, reference);
        let cigar_str = aln.cigar_string();
        // Should shift left through the A's
        assert!(
            cigar_str.starts_with("2I"),
            "Expected insertion shifted to left, got: {}",
            cigar_str
        );
        assert_eq!(aln.query_length(), 10);
        assert_eq!(aln.reference_consumed(), 8);
    }

    #[test]
    fn test_left_align_multiple_indels() {
        // Test with a simple repeat insertion
        let cigar2 = parse_cigar("4=2I2=").unwrap();
        let mut aln2 = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar: cigar2,
        };
        let query2 = b"ACACACAC";
        let reference2 = b"ACACAC";
        aln2.left_align_indels(query2, reference2);
        assert!(
            aln2.cigar_string().starts_with("2I"),
            "Expected insertion at left of repeat, got: {}",
            aln2.cigar_string()
        );
    }

    #[test]
    fn test_left_align_with_soft_clip() {
        // Soft clip at start should not be shifted past
        // CIGAR: 2S 4= 2I 4=
        let cigar = parse_cigar("2S4=2I4=").unwrap();
        let mut aln = Alignment {
            divergence: DivergenceScore::ZERO,
            cigar,
        };
        // query includes 2 soft-clipped + 10 aligned
        let query = b"NNAAAAAACCCC";
        let reference = b"AAAACCCC";

        aln.left_align_indels(query, reference);
        let cigar_str = aln.cigar_string();
        // Should shift through the A matches but not past the soft clip
        assert!(
            cigar_str.starts_with("2S"),
            "Soft clip should be preserved, got: {}",
            cigar_str
        );
        assert_eq!(aln.query_length(), 12);
    }

    fn make_aln(cigar_str: &str) -> Alignment {
        Alignment {
            divergence: DivergenceScore::ZERO,
            cigar: parse_cigar(cigar_str).unwrap(),
        }
    }

    // ── split_at_ref_pos tests ─────────────────────────────────────────────

    #[test]
    fn split_at_op_boundary() {
        let aln = make_aln("5=3X2=");
        let (l, r) = aln.split_at_ref_pos(5);
        assert_eq!(l.cigar_string(), "5=");
        assert_eq!(r.cigar_string(), "3X2=");
    }

    #[test]
    fn split_splits_op() {
        let aln = make_aln("10=");
        let (l, r) = aln.split_at_ref_pos(4);
        assert_eq!(l.cigar_string(), "4=");
        assert_eq!(r.cigar_string(), "6=");
    }

    #[test]
    fn split_across_deletion() {
        // 3= + 2D + 5= → split at ref pos 6 (3 from first match, 2 from del, 1 from second match)
        let aln = make_aln("3=2D5=");
        let (l, r) = aln.split_at_ref_pos(6);
        assert_eq!(l.cigar_string(), "3=2D1=");
        assert_eq!(r.cigar_string(), "4=");
    }

    #[test]
    fn split_at_zero() {
        let aln = make_aln("5=3D2=");
        let (l, r) = aln.split_at_ref_pos(0);
        assert_eq!(l.cigar_string(), "");
        assert_eq!(r.cigar_string(), "5=3D2=");
    }

    #[test]
    fn split_at_end() {
        let aln = make_aln("5=3D2=");
        let (l, r) = aln.split_at_ref_pos(10); // 5 + 3 + 2 = 10 total ref bases
        assert_eq!(l.cigar_string(), "5=3D2=");
        assert_eq!(r.cigar_string(), "");
    }

    #[test]
    fn split_insertion_at_boundary_goes_right() {
        // An insertion sitting exactly at the split boundary goes to the right half,
        // consistent with trim_ref_prefix which leaves boundary insertions on the
        // surviving (right) side.
        let aln = make_aln("3=2I5=");
        let (l, r) = aln.split_at_ref_pos(3);
        assert_eq!(l.cigar_string(), "3=");
        assert_eq!(r.cigar_string(), "2I5=");
    }

    #[test]
    fn split_insertion_before_boundary_goes_left() {
        // An insertion mid-left (before ref_remaining hits 0) stays in the left half.
        let aln = make_aln("2=2I3=");
        let (l, r) = aln.split_at_ref_pos(4); // split after 2= + 2I + 2 of 3=
        assert_eq!(l.cigar_string(), "2=2I2=");
        assert_eq!(r.cigar_string(), "1=");
    }

    #[test]
    fn split_mismatch_op() {
        let aln = make_aln("2=4X3=");
        let (l, r) = aln.split_at_ref_pos(4); // splits the 4X at position 2
        assert_eq!(l.cigar_string(), "2=2X");
        assert_eq!(r.cigar_string(), "2X3=");
    }
}
