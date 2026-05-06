#[cfg(feature = "conventional")]
use std::collections::HashMap;

#[cfg(feature = "conventional")]
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

#[cfg(feature = "conventional")]
use crate::align::AlignParams;
#[cfg(feature = "conventional")]
use crate::scores::QualityScore;
#[cfg(feature = "conventional")]
use crate::config;
#[cfg(feature = "conventional")]
use crate::utils::debug::{self, DebugFile, DebugOutput, DebugTsvWriter, TsvRow};
#[cfg(feature = "conventional")]
use crate::utils::join::Joinable;
#[cfg(feature = "conventional")]
use crate::utils::join::sorted_join;

/// SAM flag constants for debug/text SAM output in seeds.
const FLAG_REVERSE: u16 = 0x10;
#[cfg(feature = "conventional")]
const FLAG_SUPPLEMENTARY: u16 = 0x800;

// ── Debug file statics ──────────────────────────────────────────────────────

/// Debug SAM file with seed chains linked via SA tags.
#[cfg(feature = "conventional")]
pub(crate) static CHAINS_SAM: DebugFile<ChainsSamDebug> = DebugFile::new();

/// Debug TSV file with gaps and potential fills.
#[cfg(feature = "conventional")]
static GAP_FILLS: DebugFile<GapFillsDebug> = DebugFile::new();

/// Debug file with failed gap alignment sequences.
#[cfg(feature = "conventional")]
static GAP_ALIGNMENTS: DebugFile<GapAlignmentsDebug> = DebugFile::new();

// ── Concrete debug types ─────────────────────────────────────────────────────

#[cfg(feature = "conventional")]
pub(crate) struct ChainsSamDebug(DebugTsvWriter);

#[cfg(feature = "conventional")]
impl DebugOutput for ChainsSamDebug {
    type Item<'a> = str;
    fn create() -> Option<Self> {
        let path = &config::get().seeding.debug_chains_sam;
        if path.is_empty() { return None; }
        DebugTsvWriter::open(path, debug::sam_header().as_deref()).ok().map(Self)
    }
    fn append(&self, item: &str) { self.0.append(item); }
    fn finish(&self) { self.0.finish(); }
}

#[cfg(feature = "conventional")]
type GapFillsRow<'a> = (&'a str, usize, usize, usize, usize, i64, QualityScore, &'a str, usize, usize, char);

#[cfg(feature = "conventional")]
struct GapFillsDebug(DebugTsvWriter);

#[cfg(feature = "conventional")]
impl GapFillsDebug {
    const HEADERS: &[&str] = &[
        "read_name", "read_len", "read_start", "read_end", "fill_len",
        "cluster_idx", "aln_score", "chrom_name", "ref_start", "ref_end", "strand",
    ];
    const _CHECK: () = assert!(Self::HEADERS.len() == <GapFillsRow<'static> as TsvRow>::NUM_FIELDS);
}

#[cfg(feature = "conventional")]
impl DebugOutput for GapFillsDebug {
    type Item<'a> = GapFillsRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_gap_fills_tsv;
        if path.is_empty() { return None; }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &GapFillsRow<'_>) { self.0.append_row(item); }
    fn finish(&self) { self.0.finish(); }
}

#[cfg(feature = "conventional")]
struct GapAlignmentsDebug(DebugTsvWriter);

#[cfg(feature = "conventional")]
impl DebugOutput for GapAlignmentsDebug {
    type Item<'a> = str;
    fn create() -> Option<Self> {
        let path = &config::get().seeding.debug_gap_alignments;
        if path.is_empty() { return None; }
        DebugTsvWriter::open(path, None).ok().map(Self)
    }
    fn append(&self, item: &str) { self.0.append(item); }
    fn finish(&self) { self.0.finish(); }
}

/// Row: read_name, gap_cluster, filler_cluster,
///      gap_read_start, gap_read_end, gap_chrom, gap_ref_start, gap_ref_end, gap_strand, gap_score,
///      filler_chrom, filler_ref_start, filler_ref_end, filler_strand, filler_score,
///      ref_concordant, same_chrom, same_strand
#[cfg(feature = "conventional")]
type SplitDecisionRow<'a> = (
    &'a str, usize, usize, usize,
    usize, usize, &'a str, usize, usize, char, f64,
    &'a str, usize, usize, char, f64,
    u8, u8, u8,
);

#[cfg(feature = "conventional")]
pub(crate) struct SplitDecisionsDebug(DebugTsvWriter);

#[cfg(feature = "conventional")]
impl SplitDecisionsDebug {
    const HEADERS: &[&str] = &[
        "read_name", "read_len", "gap_cluster", "filler_cluster",
        "gap_read_start", "gap_read_end", "gap_chrom", "gap_ref_start", "gap_ref_end", "gap_strand", "gap_score",
        "filler_chrom", "filler_ref_start", "filler_ref_end", "filler_strand", "filler_score",
        "ref_concordant", "same_chrom", "same_strand",
    ];
    const _CHECK: () = assert!(Self::HEADERS.len() == <SplitDecisionRow<'static> as TsvRow>::NUM_FIELDS);
}

#[cfg(feature = "conventional")]
static SPLIT_DECISIONS: DebugFile<SplitDecisionsDebug> = DebugFile::new();

#[cfg(feature = "conventional")]
impl DebugOutput for SplitDecisionsDebug {
    type Item<'a> = SplitDecisionRow<'a>;
    fn create() -> Option<Self> {
        let _ = Self::_CHECK;
        let path = &config::get().seeding.debug_split_decisions_tsv;
        if path.is_empty() { return None; }
        let header = Self::HEADERS.join("\t");
        DebugTsvWriter::open(path, Some(&header)).ok().map(Self)
    }
    fn append(&self, item: &SplitDecisionRow<'_>) { self.0.append_row(item); }
    fn finish(&self) { self.0.finish(); }
}

#[cfg(feature = "conventional")]
#[derive(Debug)]
pub enum ClusterError {
    SequenceMismatch {
        read_bases: String,
        ref_bases: String,
        read_pos: usize,
        ref_pos: usize,
    },
    AlignmentMismatch {
        gap_index: usize,
        read_start: usize,
        read_end: usize,
        ref_start: usize,
        ref_end: usize,
        error: String,
    },
}

#[cfg(feature = "conventional")]
impl std::fmt::Display for ClusterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClusterError::SequenceMismatch {
                read_bases,
                ref_bases,
                read_pos,
                ref_pos,
            } => write!(
                f,
                "Seed at read_pos={} ref_pos={} does not match reference: read_bases='{}' ref_bases='{}'",
                read_pos, ref_pos, read_bases, ref_bases
            ),
            ClusterError::AlignmentMismatch {
                gap_index,
                read_start,
                read_end,
                ref_start,
                ref_end,
                error,
            } => write!(
                f,
                "Alignment mismatch at gap {}: read[{}-{}] ref[{}-{}]: {}",
                gap_index, read_start, read_end, ref_start, ref_end, error
            ),
        }
    }
}

#[cfg(feature = "conventional")]
impl std::error::Error for ClusterError {}

/// A seed hit representing a k-mer match between read and reference
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct SeedHit {
    /// Chromosome/contig index in the reference
    pub chrom_id: usize,
    /// Diagonal: ref_pos - read_pos (constant for colinear matches)
    pub diagonal: i64,
    /// Position in the reference sequence
    pub ref_pos: usize,
    /// Position in the read sequence
    pub read_pos: usize,
    /// Initial kmer
    pub kmer: u64,
    /// Uniqueness of the most unique k-mer incorporated in this seed (lower is more unique)
    pub kmer_uniqueness: u32,
    /// Number of times this k-mer appears in the read (1 = unique within read)
    pub read_frequency: u32,
    /// Length of the match (initially k, may be extended)
    pub match_len: usize,
}

impl SeedHit {
    /// Create a new seed hit with read_frequency defaulting to 1.
    pub fn new(
        chrom_id: usize,
        ref_pos: usize,
        read_pos: usize,
        kmer: u64,
        kmer_uniqueness: u32,
        match_len: usize,
    ) -> Self {
        Self {
            chrom_id,
            diagonal: ref_pos as i64 - read_pos as i64,
            ref_pos,
            read_pos,
            kmer,
            kmer_uniqueness,
            read_frequency: 1,
            match_len,
        }
    }

    /// Create a new seed hit with an explicit read_frequency.
    pub fn with_read_frequency(
        chrom_id: usize,
        ref_pos: usize,
        read_pos: usize,
        kmer: u64,
        kmer_uniqueness: u32,
        read_frequency: u32,
        match_len: usize,
    ) -> Self {
        Self {
            chrom_id,
            diagonal: ref_pos as i64 - read_pos as i64,
            ref_pos,
            read_pos,
            kmer,
            kmer_uniqueness,
            read_frequency,
            match_len,
        }
    }

    /// End position in the read
    pub fn read_end(&self) -> usize {
        self.read_pos + self.match_len
    }

    /// End position in the reference
    pub fn ref_end(&self) -> usize {
        self.ref_pos + self.match_len
    }
}

impl SeedHit {
    #[allow(dead_code)]
    pub fn fwd_read_range(&self, read_len: usize, is_reverse: bool) -> (usize, usize) {
        if is_reverse {
            (read_len - self.read_end(), read_len - self.read_pos)
        } else {
            (self.read_pos, self.read_end())
        }
    }

    /// Validate that this seed is an exact match between read and reference.
    /// Returns true if all bases match exactly, false otherwise.
    /// If there are mismatches, logs them for debugging.
    #[allow(dead_code)]
    pub fn validate_exact_match(&self, read_seq: &[u8], ref_seq: &[u8], context: &str) -> bool {
        let read_end = self.read_pos + self.match_len;
        let ref_end = self.ref_pos + self.match_len;

        // Bounds check
        if read_end > read_seq.len() {
            log::warn!(
                "[{}] Seed read range [{}, {}) exceeds read length {}",
                context,
                self.read_pos,
                read_end,
                read_seq.len()
            );
            return false;
        }
        if ref_end > ref_seq.len() {
            log::warn!(
                "[{}] Seed ref range [{}, {}) exceeds ref length {}",
                context,
                self.ref_pos,
                ref_end,
                ref_seq.len()
            );
            return false;
        }

        let read_slice = &read_seq[self.read_pos..read_end];
        let ref_slice = &ref_seq[self.ref_pos..ref_end];

        let mut mismatches = Vec::new();
        for (i, (&r, &q)) in ref_slice.iter().zip(read_slice.iter()).enumerate() {
            if r != q {
                mismatches.push((i, r as char, q as char));
            }
        }

        if !mismatches.is_empty() {
            log::warn!(
                "[{}] Seed at ref_pos={} read_pos={} match_len={} has {} mismatches:",
                context,
                self.ref_pos,
                self.read_pos,
                self.match_len,
                mismatches.len()
            );
            for (i, ref_base, read_base) in mismatches.iter().take(10) {
                log::warn!("  Position {}: ref={} read={}", i, ref_base, read_base);
            }
            return false;
        }

        true
    }

    /// Attempt to extend the seed hit if the new k-mer overlaps the current match
    /// or return a new seed hit if the new k-mer does not overlap.
    ///
    /// Two k-mers can be merged only if they overlap - i.e., the new k-mer starts
    /// within the current match region. This ensures no unverified gaps exist.
    pub fn extend(
        &mut self,
        chrom_id: usize,
        chrom_pos: usize,
        read_pos: usize,
        kmer: u64,
        kmer_uniqueness: u32,
        read_frequency: u32,
        k: usize,
    ) -> Option<SeedHit> {
        if chrom_id == self.chrom_id
            && chrom_pos >= self.ref_pos
            && read_pos >= self.read_pos
            && chrom_pos - self.ref_pos == read_pos - self.read_pos
            // Only merge if the new k-mer starts within our current match region
            // (i.e., it overlaps). This prevents creating gaps with unverified bases.
            && chrom_pos - self.ref_pos <= self.match_len
        {
            // Overlaps or extends current match
            let new_end = (chrom_pos - self.ref_pos) + k;
            if new_end > self.match_len {
                self.match_len = new_end;
            }
            if kmer_uniqueness < self.kmer_uniqueness {
                self.kmer_uniqueness = kmer_uniqueness;
            }
            if read_frequency > self.read_frequency {
                self.read_frequency = read_frequency;
            }
            None
        } else {
            // Does not overlap - return new seed hit
            Some(SeedHit::with_read_frequency(
                chrom_id,
                chrom_pos,
                read_pos,
                kmer,
                kmer_uniqueness,
                read_frequency,
                k,
            ))
        }
    }

    /// Format as SAM line string for debug output with proper hard clips.
    ///
    /// # Arguments
    /// * `read_id` - Read name
    /// * `chrom_name` - Chromosome name (not ID)
    /// * `is_reverse` - Whether this seed is on the reverse strand
    /// * `strand_seq` - The sequence for this strand (already rev-comped if reverse)
    /// * `strand_qual` - Quality scores for this strand (already reversed if reverse)
    pub fn to_sam_line(
        &self,
        read_id: &str,
        chrom_name: &str,
        is_reverse: bool,
        strand_seq: &[u8],
        strand_qual: &[u8],
    ) -> String {
        let read_len = strand_seq.len();
        let flag = if is_reverse {
            FLAG_REVERSE
        } else {
            0u16
        };

        // Build CIGAR with hard clips
        // For forward strand: read_pos bases before, match_len bases aligned, rest after
        // For reverse strand: coordinates are already in forward read space after revcomp
        let hclip_start = self.read_pos;
        let hclip_end = read_len.saturating_sub(self.read_pos + self.match_len);

        let cigar = match (hclip_start > 0, hclip_end > 0) {
            (true, true) => format!("{}H{}={}H", hclip_start, self.match_len, hclip_end),
            (true, false) => format!("{}H{}=", hclip_start, self.match_len),
            (false, true) => format!("{}={}H", self.match_len, hclip_end),
            (false, false) => format!("{}=", self.match_len),
        };

        let u = self.kmer_uniqueness;
        let mapq = 60 / u as u8;

        // Extract the aligned portion of the sequence
        let seq_slice = &strand_seq[self.read_pos..self.read_pos + self.match_len];
        let seq_str = String::from_utf8_lossy(seq_slice);

        // Extract the aligned portion of the quality scores, or use * if not available
        let qual_str = {
            let qual_slice = &strand_qual[self.read_pos..self.read_pos + self.match_len];
            // Convert Phred+33 quality scores to ASCII
            qual_slice.iter().map(|&q| q as char).collect::<String>()
        };

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            read_id,
            flag,
            chrom_name,
            self.ref_pos + 1, // 1-based
            mapq,
            cigar,
            seq_str,
            qual_str
        )
    }

    /// Extend the seed match bidirectionally as far as sequences match exactly.
    /// This is the minimap2-style seed extension that maximizes anchor length.
    ///
    /// Returns the number of bases extended (backward + forward).
    pub fn extend_exact(&mut self, read_seq: &[u8], ref_seq: &[u8]) -> usize {
        let mut extended = 0;

        // Extend backward
        let mut back_ext = 0;
        while self.read_pos > back_ext
            && self.ref_pos > back_ext
            && read_seq[self.read_pos - back_ext - 1] == ref_seq[self.ref_pos - back_ext - 1]
        {
            back_ext += 1;
        }
        if back_ext > 0 {
            self.read_pos -= back_ext;
            self.ref_pos -= back_ext;
            self.match_len += back_ext;
            self.diagonal = self.ref_pos as i64 - self.read_pos as i64;
            extended += back_ext;
        }

        // Extend forward
        let mut fwd_ext = 0;
        while self.ref_pos + self.match_len + fwd_ext < ref_seq.len()
            && self.read_pos + self.match_len + fwd_ext < read_seq.len()
            && read_seq[self.read_pos + self.match_len + fwd_ext]
                == ref_seq[self.ref_pos + self.match_len + fwd_ext]
        {
            fwd_ext += 1;
        }
        if fwd_ext > 0 {
            self.match_len += fwd_ext;
            extended += fwd_ext;
        }

        extended
    }
}

#[cfg(feature = "conventional")]
pub mod seed_cluster;

/// A gap in one cluster that is filled by another cluster
#[cfg(feature = "conventional")]
#[derive(Debug, Clone)]
pub struct GapFill {
    /// Index of the cluster containing the gap
    pub cluster_idx: usize,
    /// Gap range in forward-strand read coordinates
    #[allow(dead_code)]
    pub gap: (usize, usize),
    /// Index of the seed before this gap in the cluster's chain
    pub gap_seed_idx: usize,
    /// Index of the cluster that fills this gap
    #[allow(dead_code)]
    pub filler_idx: usize,
    /// Fraction of the gap covered by the filler
    #[allow(dead_code)]
    pub coverage: f64,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
#[cfg(feature = "conventional")]
struct ClusterGap {
    gap_start: usize,
    gap_end: usize,
    cluster_idx: usize,
    gap_seed_idx: usize,
    gap_score: QualityScore,
}

#[cfg(feature = "conventional")]
impl ClusterGap {
    fn new(
        gap_start: usize,
        gap_end: usize,
        cluster_idx: usize,
        gap_seed_idx: usize,
        gap_score: QualityScore,
    ) -> Self {
        Self {
            gap_start,
            gap_end,
            cluster_idx,
            gap_seed_idx,
            gap_score,
        }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.gap_end - self.gap_start
    }
}

#[cfg(feature = "conventional")]
impl Joinable for ClusterGap {
    fn range(&self) -> (usize, usize) {
        (self.gap_start, self.gap_end)
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
#[cfg(feature = "conventional")]
struct ClusterAlignment {
    cluster_start: usize,
    cluster_end: usize,
    cluster_idx: usize,
    cluster_score: QualityScore,
    identity_pct: f64,
}

#[cfg(feature = "conventional")]
impl ClusterAlignment {
    fn new(
        cluster_start: usize,
        cluster_end: usize,
        cluster_idx: usize,
        cluster_score: QualityScore,
        identity_pct: f64,
    ) -> Self {
        Self {
            cluster_start,
            cluster_end,
            cluster_idx,
            cluster_score,
            identity_pct,
        }
    }

    #[allow(dead_code)]
    fn len(&self) -> usize {
        self.cluster_end - self.cluster_start
    }
}

#[cfg(feature = "conventional")]
impl Joinable for ClusterAlignment {
    fn range(&self) -> (usize, usize) {
        (self.cluster_start, self.cluster_end)
    }
}

/// Analyze all gaps across all clusters to find chimeric breakpoints
///
/// Returns a list of gaps that are filled by other clusters, indicating
/// potential chimeric read breakpoints.
#[cfg(feature = "conventional")]
pub fn analyze_gap_fills<'a>(
    read_name: &str,
    clusters: &[seed_cluster::SeedCluster],
    read_len: usize,
    min_gap: usize,
    min_fill: usize,
    tolerance: usize,
    params: &AlignParams,
    chrom_name: &'a dyn Fn(usize) -> &'a str,
) -> Vec<GapFill> {
    if clusters.len() < 2 {
        return Vec::new();
    }
    let start = std::time::Instant::now();

    let dump = GAP_FILLS.is_enabled();

    let mut gaps: Vec<ClusterGap> = Vec::new();
    let mut alignments: Vec<ClusterAlignment> = Vec::new();
    for (cluster_idx, cluster) in clusters.iter().enumerate() {
        // Add alignment region for this cluster
        let (fwd_start, fwd_end) = cluster.fwd_read_range(read_len);
        if fwd_end > fwd_start && fwd_end - fwd_start >= min_fill {
            let score = cluster.quality(params);
            let identity = cluster.total_identity();
            let total_read_length = fwd_end - fwd_start;
            let identity_pct = if total_read_length > 0 {
                100.0 * identity as f64 / total_read_length as f64
            } else {
                0.0
            };
            if identity_pct > 90.0 {
                log::debug!(
                    "Cluster {} alignment: read ({},{}) ref ({}-{}) length {} score {:.2} identity {:.2}%",
                    cluster_idx,
                    fwd_start,
                    fwd_end,
                    cluster.ref_start(),
                    cluster.ref_end(),
                    total_read_length,
                    score,
                    identity_pct
                );
                alignments.push(ClusterAlignment::new(
                    fwd_start,
                    fwd_end,
                    cluster_idx,
                    score,
                    identity_pct,
                ));
                if dump {
                    let chrom_name = format!("chr{}", cluster.chrom_id + 1); // it's a lie, but close enough for debugging
                    let strand = if cluster.is_reverse { '-' } else { '+' };
                    let _fill_len = total_read_length;
                    GAP_FILLS.append(&(
                        read_name, read_len, fwd_start, fwd_end, total_read_length,
                        -(cluster_idx as i64 + 1), score, chrom_name.as_str(),
                        cluster.ref_start(), cluster.ref_end(), strand,
                    ));
                }
            }
        }

        // Find gaps in this cluster
        // Seeds are stored in strand order, so for reverse strand the forward coords
        // are in reverse order. We need to compute gaps in forward coordinate space.
        // Skip clusters that haven't been gap-aligned (they can still act as fillers
        // via the alignment region above, but their gaps can't be scored).
        if cluster.chain.len() < 2 || cluster.gap_alignments.is_empty() {
            continue;
        }
        let n = cluster.chain.len();
        for i in 0..n - 1 {
            let seed1 = &cluster.chain[i];
            let seed1_range = seed1.fwd_read_range(read_len, cluster.is_reverse);
            let seed2 = &cluster.chain[i + 1];
            let seed2_range = seed2.fwd_read_range(read_len, cluster.is_reverse);

            // Compute gap in forward coordinates
            // For forward strand: gap is from seed1.end to seed2.start
            // For reverse strand: seeds are in reverse fwd order, so gap is from seed2.end to seed1.start
            let (gap_start, gap_end) = if cluster.is_reverse {
                (seed2_range.1, seed1_range.0)
            } else {
                (seed1_range.1, seed2_range.0)
            };

            if gap_end > gap_start && gap_end - gap_start >= min_gap {
                let aln = cluster.gap_alignment(i).unwrap();
                let aln_score = cluster
                    .gap_alignment(i)
                    .map(|a| a.quality(params).0)
                    .unwrap_or(0.0);
                let ref_start = seed1.ref_end();
                let ref_end = seed2.ref_pos;
                let identity_pct = if aln.query_length() > 0 {
                    100.0 * aln.total_identity() as f64 / aln.query_length() as f64
                } else {
                    0.0
                };
                log::debug!(
                    "Cluster {} gap between seeds {} and {}: read ({},{}) ref ({}-{}) length {} score {:.2} identity {:.2}% | CIGAR: {}",
                    cluster_idx,
                    i,
                    i + 1,
                    gap_start,
                    gap_end,
                    ref_start,
                    ref_end,
                    gap_end - gap_start,
                    aln_score,
                    identity_pct,
                    aln.cigar_string()
                );
                gaps.push(ClusterGap::new(
                    gap_start,
                    gap_end,
                    cluster_idx,
                    i,
                    QualityScore::new(aln_score),
                ));
                if dump {
                    let chrom_name = format!("chr{}", cluster.chrom_id + 1); // it's a lie, but close enough for debugging
                    let strand = if cluster.is_reverse { '-' } else { '+' };
                    let fill_len = gap_end - gap_start;
                    let (ref_start, ref_end) = if cluster.is_reverse {
                        (seed2.ref_end(), seed1.ref_pos)
                    } else {
                        (seed1.ref_end(), seed2.ref_pos)
                    };
                    GAP_FILLS.append(&(
                        read_name, read_len, gap_start, gap_end, fill_len,
                        (cluster_idx + 1) as i64, QualityScore::new(aln_score), chrom_name.as_str(),
                        ref_start, ref_end, strand,
                    ));
                }
            }
        }
    }
    gaps.sort_by_key(|gap| (gap.gap_start, gap.gap_end, gap.cluster_idx));
    alignments.sort_by_key(|aln| (aln.cluster_start, aln.cluster_end, aln.cluster_idx));

    // Use sorted_join to find gap-alignment pairs from different clusters
    let mut pairs = sorted_join(&gaps, &alignments, tolerance, |gap, aln| {
        let gap_len = (gap.gap_end - gap.gap_start) as isize;
        let aln_len = (aln.cluster_end - aln.cluster_start) as isize;

        if gap.cluster_idx == aln.cluster_idx {
            return false;
        }
        if !gap.gap_score.is_worse_than(QualityScore::ZERO) {
            return false;
        }
        if !aln.cluster_score.is_better_than(QualityScore::ZERO) {
            return false;
        }

        // Reject ref-concordant pairs: if the filler's reference region overlaps
        // the gap cluster's reference region on the same chrom and strand, the
        // filler is likely aliasing within a repeat, not evidence of an SV.
        let gap_cluster = &clusters[gap.cluster_idx];
        let filler_cluster = &clusters[aln.cluster_idx];
        if gap_cluster.chrom_id == filler_cluster.chrom_id
            && gap_cluster.is_reverse == filler_cluster.is_reverse
            && filler_cluster.ref_start() < gap_cluster.ref_end()
            && gap_cluster.ref_start() < filler_cluster.ref_end()
        {
            return false;
        }

        let overlap_start = gap.gap_start.max(aln.cluster_start);
        let overlap_end = gap.gap_end.min(aln.cluster_end);
        if overlap_end <= overlap_start {
            return false;
        }
        let overlap = overlap_end - overlap_start;
        let gap_overlap_ratio = overlap as f64 / gap_len as f64;
        let fill_overlap_ratio = aln_len as f64 / overlap as f64;
        let ratio = gap_overlap_ratio * fill_overlap_ratio;
        if ratio >= 0.5 && ratio <= 1.1 {
            log::debug!(
                "Considering gap fill: gap cluster {} gap ({},{}) length {} score {:.2} | aln cluster {} aln ({},{}) length {} score {:.2} identity {:.2}% | gap ratio {:.2} fill ratio {:.2} | ratio {:.2}",
                gap.cluster_idx,
                gap.gap_start,
                gap.gap_end,
                gap_len,
                gap.gap_score,
                aln.cluster_idx,
                aln.cluster_start,
                aln.cluster_end,
                aln_len,
                aln.cluster_score,
                aln.identity_pct,
                gap_overlap_ratio,
                fill_overlap_ratio,
                ratio
            );
        }
        // Must be from different clusters
        gap.cluster_idx != aln.cluster_idx
            && gap.gap_score.is_worse_than(QualityScore::ZERO)
            && aln.cluster_score.is_better_than(QualityScore::ZERO)
            && ratio >= 0.5
            && ratio <= 1.1
    });

    pairs.sort_by_key(|(gap_idx, aln_idx)| {
        let gap = &gaps[*gap_idx];
        let aln = &alignments[*aln_idx];
        OrderedFloat(gap.gap_score.0 - aln.cluster_score.0)
    });

    let mut selected = HashMap::new();
    for (gap_idx, aln_idx) in &pairs {
        if selected.contains_key(gap_idx) {
            continue;
        }
        selected.insert(*gap_idx, *aln_idx);
        let gap = &gaps[*gap_idx];
        let aln = &alignments[*aln_idx];
        log::debug!(
            "Selected gap fill pair: gap cluster {} gap ({},{}) length {} score {:.2} | aln cluster {} aln ({},{}) length {} score {:.2} | score diff {:.2}",
            gap.cluster_idx,
            gap.gap_start,
            gap.gap_end,
            gap.len(),
            gap.gap_score,
            aln.cluster_idx,
            aln.cluster_start,
            aln.cluster_end,
            aln.len(),
            aln.cluster_score,
            aln.cluster_score.0 - gap.gap_score.0
        );
    }
    let pairs: Vec<(usize, usize)> = selected.into_iter().collect();

    // Convert pairs to GapFill entries with coverage calculation
    let fills: Vec<GapFill> = pairs
        .into_iter()
        .map(|(gap_idx, aln_idx)| {
            let gap = &gaps[gap_idx];
            let aln = &alignments[aln_idx];
            let gap_len = gap.gap_end - gap.gap_start;
            let aln_len = aln.cluster_end - aln.cluster_start;
            let ratio = aln_len as f64 / gap_len as f64;
            let diff = gap_len as isize - aln_len as isize;

            let qual = gap.gap_score.0 - aln.cluster_score.0;

            // Look up reference coordinates for gap and filler clusters
            let gap_cluster = &clusters[gap.cluster_idx];
            let filler_cluster = &clusters[aln.cluster_idx];

            let gap_chrom = chrom_name(gap_cluster.chrom_id);
            let filler_chrom = chrom_name(filler_cluster.chrom_id);

            let same_chrom = gap_cluster.chrom_id == filler_cluster.chrom_id;
            let same_strand = gap_cluster.is_reverse == filler_cluster.is_reverse;

            // Ref-concordance: does the filler's ref region overlap the gap's ref region?
            // For real SVs, the filler will map elsewhere; for spurious splits, it overlaps.
            let gap_ref_start = gap_cluster.ref_start();
            let gap_ref_end = gap_cluster.ref_end();
            let filler_ref_start = filler_cluster.ref_start();
            let filler_ref_end = filler_cluster.ref_end();
            let ref_concordant = same_chrom
                && same_strand
                && filler_ref_start < gap_ref_end
                && gap_ref_start < filler_ref_end;

            log::debug!(
                "Gap fill: cluster {} gap ({},{} - {}bp, score {:.2}) filled by cluster {} aln ({},{} - {}bp, score {:.2}) | ratio {:.2}, diff {}, qual {:.2} | ref_concordant={}, same_chrom={}, same_strand={}",
                gap.cluster_idx,
                gap.gap_start,
                gap.gap_end,
                gap_len,
                gap.gap_score,
                aln.cluster_idx,
                aln.cluster_start,
                aln.cluster_end,
                aln_len,
                aln.cluster_score,
                ratio,
                diff,
                qual,
                ref_concordant,
                same_chrom,
                same_strand,
            );

            // Emit split-decision diagnostic row
            if SPLIT_DECISIONS.is_enabled() {
                let gap_strand = if gap_cluster.is_reverse { '-' } else { '+' };
                let filler_strand = if filler_cluster.is_reverse { '-' } else { '+' };
                SPLIT_DECISIONS.append(&(
                    read_name, read_len, gap.cluster_idx, aln.cluster_idx,
                    gap.gap_start, gap.gap_end, gap_chrom, gap_ref_start, gap_ref_end, gap_strand, gap.gap_score.0,
                    filler_chrom, filler_ref_start, filler_ref_end, filler_strand, aln.cluster_score.0,
                    ref_concordant as u8, same_chrom as u8, same_strand as u8,
                ));
            }

            // Calculate coverage: fraction of gap covered by the alignment
            let gap_len = gap.gap_end - gap.gap_start;
            let overlap_start = gap
                .gap_start
                .max(aln.cluster_start.saturating_sub(tolerance));
            let overlap_end = gap.gap_end.min(aln.cluster_end.saturating_add(tolerance));
            let overlap = overlap_end.saturating_sub(overlap_start);
            let coverage = overlap as f64 / gap_len as f64;

                GapFill {
                    cluster_idx: gap.cluster_idx,
                    gap: (gap.gap_start, gap.gap_end),
                    gap_seed_idx: gap.gap_seed_idx,
                    filler_idx: aln.cluster_idx,
                    coverage,
                }
        })
        .collect();

    metrics::histogram!("analysis_gap_fills").record(start.elapsed().as_micros() as f64);

    fills
}

#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Read {
    pub name: String,
    pub seq: String,
    pub qual: String,
}

impl Read {
    #[allow(dead_code)]
    pub fn new(name: &str, seq: &[u8], qual: &[u8]) -> Self {
        Self {
            name: name.to_string(),
            seq: String::from_utf8_lossy(seq).to_string(),
            qual: String::from_utf8_lossy(qual).to_string(),
        }
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedSaver {
    pub read: Read,
    pub is_reverse: bool,
    pub seeds: Vec<SeedHit>,
}

#[cfg(all(test, feature = "conventional"))]
mod tests {
    use super::*;

    #[test]
    fn test_format_seed_diagram() {
        // Create a chain with 3 seeds and gaps between them
        let chain = vec![
            SeedHit::new(0, 100, 0, 0, 1, 65),   // 65bp seed at ref 100, read 0
            SeedHit::new(0, 176, 73, 0, 1, 111), // 111bp seed at ref 176, read 73 (gap: qry=8, ref=11)
            SeedHit::new(0, 304, 203, 0, 1, 44), // 44bp seed at ref 304, read 203 (gap: qry=19, ref=17)
        ];

        let cluster = seed_cluster::SeedCluster::new(chain, false, 10).unwrap();
        let (qry, ref_line) = cluster.format_seed_diagram();

        println!("{}", qry);
        println!("{}", ref_line);

        // Check that seeds are aligned (same width on both lines)
        assert!(qry.contains("[- 65bp -]"));
        assert!(ref_line.contains("[- 65bp -]"));
        assert!(qry.contains("[- 111bp -]"));
        assert!(ref_line.contains("[- 111bp -]"));
        assert!(qry.contains("[- 44bp -]"));
        assert!(ref_line.contains("[- 44bp -]"));

        // Check that gaps show different values
        assert!(qry.contains("8bp"));
        assert!(ref_line.contains("11bp"));
        assert!(qry.contains("19bp"));
        assert!(ref_line.contains("17bp"));

        // Both lines should have the same length (seeds aligned)
        assert_eq!(qry.len(), ref_line.len());
    }

    #[test]
    fn test_format_seed_diagram_single_seed() {
        let chain = vec![SeedHit::new(0, 100, 0, 0, 1, 50)];
        let cluster = seed_cluster::SeedCluster::new(chain, false, 10).unwrap();
        let (qry, ref_line) = cluster.format_seed_diagram();

        assert_eq!(qry, "QRY: [- 50bp -]");
        assert_eq!(ref_line, "REF: [- 50bp -]");
    }

    #[test]
    fn test_format_seed_diagram_empty() {
        // An empty chain returns None from SeedCluster::new
        // but we can test the method directly by calling format_seed_diagram
        // on a cluster we construct manually
        let cluster = seed_cluster::SeedCluster {
            read_start: 0,
            read_end: 0,
            chain: vec![],
            is_reverse: false,
            chrom_id: 0,
            gap_alignments: vec![],
            split_fill_tags: vec![],
        };
        let (qry, ref_line) = cluster.format_seed_diagram();

        assert_eq!(qry, "QRY:");
        assert_eq!(ref_line, "REF:");
    }
}
