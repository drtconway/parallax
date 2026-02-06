//! Global configuration system for parallax.
//!
//! Uses `confique` for self-documenting TOML configuration files.
//! Configuration is loaded once at startup and accessible globally via `config::get()`.
#![allow(dead_code)]
use std::sync::OnceLock;

use confique::Config;

static CONFIG: OnceLock<ParallaxConfig> = OnceLock::new();

/// Root configuration for parallax aligner.
#[derive(Config, Debug, Clone)]
pub struct ParallaxConfig {
    /// Alignment scoring parameters
    #[config(nested)]
    pub alignment: AlignmentConfig,

    /// Seeding and chaining parameters
    #[config(nested)]
    pub seeding: SeedingConfig,

    /// Alignment filtering thresholds
    #[config(nested)]
    pub filtering: FilteringConfig,

    /// Classification parameters for primary/secondary/supplementary
    #[config(nested)]
    pub classification: ClassificationConfig,

    /// Block aligner parameters for SIMD-accelerated gap filling
    #[config(nested)]
    pub block_aligner: BlockAlignerConfig,
}

/// Alignment scoring parameters.
///
/// These control the WFA (Wavefront Alignment) algorithm and
/// context-aware scoring for homopolymers and STRs.
#[derive(Config, Debug, Clone)]
pub struct AlignmentConfig {
    /// Match score (positive value)
    #[config(default = 2)]
    pub match_score: i32,

    /// Mismatch penalty (positive value, higher = more penalty)
    #[config(default = 4)]
    pub mismatch: i32,

    /// Gap open penalty (first base of gap)
    #[config(default = 4)]
    pub gap_open: i32,

    /// Gap extend penalty for short gaps (linear portion)
    #[config(default = 2)]
    pub gap_extend: i32,

    /// Gap length threshold where sublinear scaling kicks in
    #[config(default = 10)]
    pub sublinear_threshold: u32,

    /// Sublinear coefficient for long gaps.
    /// Penalty = gap_open + gap_extend * threshold + sublinear_coef * log2(len - threshold + 1)
    #[config(default = 4.0)]
    pub sublinear_coef: f64,

    /// Minimum homopolymer length to trigger reduced penalty
    #[config(default = 4)]
    pub homopolymer_min_len: usize,

    /// Penalty multiplier for gaps in homopolymer context (0.0 - 1.0)
    #[config(default = 0.5)]
    pub homopolymer_discount: f64,

    /// Minimum repeat unit count for STR discount (e.g., 3 means ATATAT or CAGCAGCAG)
    #[config(default = 3)]
    pub str_min_repeats: usize,

    /// Penalty multiplier for gaps in STR context (0.0 - 1.0)
    #[config(default = 0.6)]
    pub str_discount: f64,

    /// X-drop threshold for WFA pruning.
    /// Diagonals that fall more than this distance behind the best diagonal
    /// are pruned to accelerate alignment of divergent sequences.
    /// Set to 0 to disable pruning.
    #[config(default = 400)]
    pub x_drop: i32,

    /// Maximum wavefront band width before giving up.
    /// If the wavefront spans more diagonals than this, alignment fails early.
    /// This prevents runaway on highly divergent or unrelated sequences.
    /// Set to 0 to disable this limit.
    #[config(default = 2000)]
    pub max_band_width: i32,
}

/// Seeding and chaining parameters.
///
/// These control how k-mer seeds are collected and clustered into chains.
#[derive(Config, Debug, Clone)]
pub struct SeedingConfig {
    /// Maximum occurrences for a seed to be used (filters highly repetitive k-mers)
    #[config(default = 50)]
    pub max_seed_occurrences: usize,

    /// Minimum distance between seed clusters for DBSCAN clustering
    #[config(default = 100)]
    pub min_seed_cluster_distance: i64,

    /// Variance coefficient for DBSCAN clustering.
    /// Max variance = (read_len * variance_coef)^2
    #[config(default = 0.01)]
    pub variance_coef: f64,

    /// Minimum total seed length for a single-seed chain to be considered
    #[config(default = 50)]
    pub min_single_seed_length: usize,

    /// Minimum gap size (bp) to consider for chimeric splitting.
    /// Gaps smaller than this are bridged with WFA instead.
    #[config(default = 100)]
    pub min_gap_for_split: usize,

    /// Maximum gap length (bp) to attempt WFA alignment on.
    /// Gaps larger than this in either read or reference are not aligned;
    /// instead, an insertion and/or deletion is emitted directly.
    /// This prevents very slow WFA calls on large introns or structural variants.
    #[config(default = 5000)]
    pub max_gap_length: usize,

    /// Tolerance (bp) for matching cluster ranges to gaps.
    /// Allows slight overlaps when detecting gap fills.
    #[config(default = 25)]
    pub gap_fill_tolerance: usize,

    /// Minimum fraction of gap that must be covered by another cluster
    /// to trigger a split (0.0-1.0).
    #[config(default = 0.5)]
    pub min_gap_fill_coverage: f64,

    /// Linear coefficient (α) for gap penalty in chain scoring.
    /// Gap penalty = α * gap_len + β * log2(gap_len)
    /// Higher values penalize gaps more heavily.
    #[config(default = 0.01)]
    pub gap_penalty_linear: f64,

    /// Logarithmic coefficient (β) for gap penalty in chain scoring.
    /// Gap penalty = α * gap_len + β * log2(gap_len)
    /// This penalizes the "surprise" of a gap - even small gaps incur some penalty.
    #[config(default = 0.5)]
    pub gap_penalty_log: f64,

    /// Minimum chain score for a cluster to be considered for alignment.
    /// Clusters scoring below this threshold are discarded.
    /// Set to 0 to disable (keep all clusters).
    #[config(default = 50.0)]
    pub min_chain_score: f64,

    /// Enable anchor-first chaining algorithm.
    /// When true, long seeds are chained first to establish the diagonal,
    /// then short seeds are added only if they're consistent with that diagonal.
    /// This helps avoid spurious indels caused by off-diagonal short seeds.
    #[config(default = false)]
    pub use_anchor_first_chaining: bool,

    /// Minimum seed length to be considered an "anchor" in anchor-first chaining.
    /// Seeds shorter than this are "fillers" that can only be added if they match
    /// the diagonal established by anchors.
    #[config(default = 40)]
    pub anchor_min_length: usize,

    /// Maximum diagonal deviation for filler seeds in anchor-first chaining.
    /// Fillers must have a diagonal within this many bases of the expected diagonal
    /// (interpolated from surrounding anchors).
    #[config(default = 50)]
    pub filler_diagonal_tolerance: i64,

    /// Maximum number of chains to keep per chromosome.
    /// Chains are extracted in descending score order, so this keeps the best ones.
    /// Set to 0 to disable the limit (keep all chains).
    #[config(default = 20)]
    pub max_chains_per_chrom: usize,

    /// Path to write debug SAM file with extended seeds (before clustering).
    /// Useful for visualizing seed placement in IGV alongside final alignments.
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_seeds_sam: String,

    /// Path to write debug TSV file with candidate alignments (before classification).
    /// Columns: read_name, read_start, read_end, read_len, chrom, ref_start, ref_end, strand, score
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_seeds_tsv: String,

    /// Path to write debug JSON file with seeds.
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_seeds_json: String,
    
    /// Path to write debug TSV file with seeds grouped into clusters (before chaining).
    /// Columns: read_name, cluster_id, read_start, read_end, read_len, chrom, ref_start, ref_end, strand, match_len
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_clusters_tsv: String,

    /// Path to write debug TSV file with seed chains/clusters (after chaining, before alignment).
    /// Columns: read_name, cluster_id, read_start, read_end, read_len, chrom, ref_start, ref_end,
    ///          strand, num_seeds, seed_length, coverage, density
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_chains_tsv: String,

    /// Path to write debug SAM file with seed chains linked via SA tags.
    /// Each chain is output with seeds as supplementary alignments, allowing
    /// visualization of chaining in IGV. Leave empty to disable.
    #[config(default = "")]
    pub debug_chains_sam: String,

    /// Path to write debug TSV file with gaps and potential fills.
    /// Columns: read_name, read_len, gap_start, gap_end, fill_len, cluster_id, aln_score,
    ///          chrom, ref_start, ref_end, strand
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_gap_fills_tsv: String,

    /// Path to write debug file with failed alignment strings for gap fills.
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_gap_alignments: String,
}

/// Alignment filtering thresholds.
///
/// These control which alignments are kept vs filtered as low quality.
#[derive(Config, Debug, Clone)]
pub struct FilteringConfig {
    /// Minimum alignment identity (matches / aligned_length) for a valid alignment
    #[config(default = 0.5)]
    pub min_identity: f64,

    /// Maximum context-aware score per aligned base (higher = more errors allowed)
    #[config(default = 0.3)]
    pub max_score_per_base: f64,

    /// Minimum fraction of read covered for a valid alignment
    #[config(default = 0.1)]
    pub min_read_coverage: f64,

    /// Minimum aligned length (bp) - alignments meeting this bypass coverage check.
    /// This handles chimeric reads where a small portion aligns elsewhere.
    #[config(default = 50)]
    pub min_aligned_length: u32,

    /// Reference overlap threshold for deduplicating candidate alignments.
    /// Two alignments are considered duplicates if their reference overlap
    /// fraction (for both) exceeds this threshold. The one with better
    /// read coverage is kept. Set to 1.0 to disable deduplication.
    #[config(default = 0.8)]
    pub ref_overlap_threshold: f64,
}

/// Classification parameters for primary/secondary/supplementary alignments.
#[derive(Config, Debug, Clone)]
pub struct ClassificationConfig {
    /// Overlap threshold for clustering alignments by read region.
    /// Two alignments are in the same cluster if both have > this fraction overlap.
    #[config(default = 0.5)]
    pub overlap_threshold: f64,

    /// Gap open penalty for set scoring.
    /// Applied once per break between alignments in a set (affine gap model).
    /// Higher values penalize fragmented alignments more.
    #[config(default = 50)]
    pub set_gap_open: i64,

    /// Gap extend penalty for set scoring.
    /// Applied per uncovered base in the read (affine gap model).
    #[config(default = 2)]
    pub set_gap_extend: i64,

    /// Use information-theoretic scoring for alignment ranking.
    /// When true, uses N*log2(N) scoring for match runs which rewards
    /// contiguous matches over scattered ones.
    #[config(default = false)]
    pub use_information_score: bool,

    /// Mismatch penalty for information-theoretic scoring.
    #[config(default = 4.0)]
    pub info_mismatch_penalty: f64,

    /// Gap open penalty for information-theoretic scoring.
    #[config(default = 6.0)]
    pub info_gap_open: f64,

    /// Gap extend penalty for information-theoretic scoring.
    #[config(default = 1.0)]
    pub info_gap_extend: f64,
}

/// Block aligner configuration for SIMD-accelerated alignment.
///
/// These parameters control the block-aligner library used for
/// filling gaps between seeds and extending alignments.
#[derive(Config, Debug, Clone)]
pub struct BlockAlignerConfig {
    /// Enable extension alignment at the ends of seed chains.
    /// When true, aligns the unaligned portions of the read beyond the first/last seeds.
    /// When false, these regions are soft-clipped.
    #[config(default = true)]
    pub enable_extension: bool,

    /// Minimum block size for SIMD alignment (must be power of 2, >= 32).
    /// Smaller values allow finer-grained alignment but may be slower.
    #[config(default = 32)]
    pub min_block_size: usize,

    /// Maximum block size for SIMD alignment (must be power of 2, <= 16384).
    /// Larger values handle longer sequences but use more memory.
    #[config(default = 4096)]
    pub max_block_size: usize,

    /// X-drop threshold for extension alignment.
    /// Alignment stops when score drops this far below the maximum.
    /// Higher values = more aggressive extension, lower = more conservative.
    /// Use ~400 for long reads (PacBio/ONT), ~200 for short reads (Illumina).
    #[config(default = 400)]
    pub x_drop: i32,

    /// Bonus score for reaching the end of a sequence during extension.
    /// Encourages alignments to extend fully rather than stopping early.
    #[config(default = 5)]
    pub end_bonus: i32,

    /// Mismatch penalty (positive value).
    /// Used for scoring matrix construction.
    #[config(default = 4)]
    pub mismatch: i32,

    /// Gap open penalty (positive value, cost of starting a gap).
    #[config(default = 6)]
    pub gap_open: i32,

    /// Gap extend penalty (positive value, cost per additional gap base).
    #[config(default = 2)]
    pub gap_extend: i32,
}

impl Default for BlockAlignerConfig {
    fn default() -> Self {
        Self {
            enable_extension: false,
            min_block_size: 32,
            max_block_size: 4096,
            x_drop: 400,
            end_bonus: 5,
            mismatch: 4,
            gap_open: 6,
            gap_extend: 2,
        }
    }
}

/// Initialize the global configuration (call once at startup).
///
/// # Panics
/// Panics if called more than once.
pub fn init(config: ParallaxConfig) {
    CONFIG
        .set(config)
        .expect("Config already initialized - init() called twice");
}

/// Load configuration from a TOML file, falling back to defaults.
///
/// If `path` is None, returns default configuration.
pub fn load(path: Option<&std::path::Path>) -> Result<ParallaxConfig, confique::Error> {
    match path {
        Some(p) => ParallaxConfig::builder().file(p).load(),
        None => ParallaxConfig::builder().load(),
    }
}

/// Get a reference to the global configuration.
///
/// Returns the initialized config, or default config if not explicitly initialized.
/// This allows tests to run without explicit initialization.
pub fn get() -> &'static ParallaxConfig {
    CONFIG.get_or_init(|| {
        ParallaxConfig::builder()
            .load()
            .expect("Failed to load default config")
    })
}

/// Generate a TOML template with all parameters and their documentation.
pub fn generate_template() -> String {
    confique::toml::template::<ParallaxConfig>(confique::toml::FormatOptions::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = ParallaxConfig::builder().load().unwrap();
        assert_eq!(config.alignment.mismatch, 4);
        assert_eq!(config.seeding.max_seed_occurrences, 50);
        assert_eq!(config.filtering.min_identity, 0.5);
    }

    #[test]
    fn test_template_generation() {
        let template = generate_template();
        assert!(template.contains("mismatch"));
        assert!(template.contains("max_seed_occurrences"));
        assert!(template.contains("min_identity"));
    }
}
