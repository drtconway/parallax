//! Global configuration system for parallax.
//!
//! Uses `confique` for self-documenting TOML configuration files.
//! Configuration is loaded once at startup and accessible globally via `config::get()`.
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

    /// Block aligner parameters for SIMD-accelerated gap filling
    #[config(nested)]
    pub block_aligner: BlockAlignerConfig,

    /// Metrics / histogram parameters
    #[config(nested)]
    pub metrics: MetricsConfig,
}

/// Alignment scoring parameters.
///
/// These control the block aligner used for gap filling.
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

    /// Gap open penalty for long gaps (second piece in two-piece affine model).
    /// Used by the WFA2 aligner. Set to 0 to disable two-piece scoring.
    #[config(default = 24)]
    pub gap_open2: i32,

    /// Gap extend penalty for long gaps (second piece in two-piece affine model).
    /// Used by the WFA2 aligner.
    #[config(default = 1)]
    pub gap_extend2: i32,

    /// X-drop threshold for block aligner pruning.
    /// Alignments that fall more than this distance behind the best score
    /// are pruned to accelerate alignment of divergent sequences.
    /// Set to 0 to disable pruning.
    #[config(default = 400)]
    pub x_drop: i32,

    /// Maximum band width before giving up (unused by block aligner, retained
    /// for compatibility with the attic WFA aligner).
    /// Set to 0 to disable this limit.
    #[config(default = 2000)]
    pub max_band_width: i32,
}

/// Seeding and chaining parameters.
///
/// These control how k-mer seeds are collected and clustered into chains.
#[derive(Config, Debug, Clone)]
pub struct SeedingConfig {
    /// Length of the syncmer k-mers used to build the index.  Must match the K
    /// parameter the index was built with.  Used in the edge penalty to ensure
    /// that a heavily-overlapped seed still contributes at least one full k-mer's
    /// worth of weight, since the underlying index guarantees that anchor exists.
    #[config(default = 20)]
    pub kmer_length: usize,

    /// Maximum occurrences for a seed to be used (filters highly repetitive k-mers)
    #[config(default = 500)]
    pub max_seed_occurrences: usize,

    /// Seeds with occurrence count above this threshold are deferred during
    /// initial collection and only rescued into gaps where no low-frequency
    /// seeds fall. Seeds at or below this threshold are collected immediately.
    /// Must be less than `max_seed_occurrences`. Set to 0 to disable rescue
    /// (all seeds up to `max_seed_occurrences` are used directly).
    #[config(default = 10)]
    pub mid_seed_occurrences: usize,

    /// Minimum gap (in read bp) between adjacent seeds that triggers rescue
    /// of deferred high-frequency seeds. Also controls the rescue rate: at
    /// most one seed is rescued per `rescue_spacing` bp of gap.
    #[config(default = 500)]
    pub rescue_spacing: usize,

    /// Minimum total seed length for a single-seed chain to be considered
    #[config(default = 75)]
    pub min_single_seed_length: usize,

    /// Maximum number of seeds in a segment for it to be treated as an excursion
    /// candidate by `ExcursionSegmentFilter`. Set to 0 to disable the filter.
    #[config(default = 3)]
    pub excursion_max_seeds: usize,

    /// Maximum reference span (bp) of an excursion candidate segment.
    #[config(default = 100)]
    pub excursion_max_ref_span: usize,

    /// Minimum gap size (bp) to consider for chimeric splitting.
    /// Gaps smaller than this are bridged with block aligner instead.
    #[config(default = 100)]
    pub min_gap_for_split: usize,

    /// Tolerance (bp) for matching cluster ranges to gaps.
    /// Allows slight overlaps when detecting gap fills.
    #[config(default = 25)]
    pub gap_fill_tolerance: usize,

    /// Maximum read-space overlap (bp) between two clusters to still merge
    /// them as a deletion. Small overlaps arise from microhomology at
    /// breakpoints. Set to 0 to require strictly abutting clusters.
    #[config(default = 50)]
    pub del_merge_max_read_overlap: usize,

    /// Maximum read-space gap (bp) between two clusters to still merge
    /// them as a deletion. For a genuine deletion the read is continuous,
    /// so any gap is due to seed placement granularity.
    #[config(default = 50)]
    pub del_merge_max_read_gap: usize,

    /// Maximum reference-space gap (bp) to bridge by merging clusters.
    /// Larger gaps are kept as separate supplementary alignments.
    #[config(default = 100000)]
    pub del_merge_max_ref_gap: usize,

    /// Threshold for filtering misplaced seeds that would require simultaneous
    /// insertion and deletion during gap alignment. Specifically, the threshold
    /// is on `2 * min(accumulated_insertion, accumulated_deletion)` across
    /// nearby long-gap transitions. This is minimap2's `mm_filter_bad_seeds()`
    /// heuristic. Set to 0 to disable. Default 40.
    #[config(default = 40)]
    pub misplaced_seed_threshold: i64,

    /// Maximum ratio of total diagonal shift to reference span within a
    /// sliding window before seeds are considered "jittery" and removed.
    /// This detects regions where seeds from different repeat copies cause
    /// the diagonal to bounce around — the gap-fill DP aligner produces
    /// better results without these misleading anchors.
    /// Set to 0.0 to disable. Default 0.15 (15 bp shift per 100 bp span).
    #[config(default = 0.15)]
    pub jitter_density_threshold: f64,

    /// Number of inter-seed gaps in the sliding window for jitter detection.
    /// The window must contain at least 2 gaps to measure density, so the
    /// minimum effective value is 2. Default 4.
    #[config(default = 4)]
    pub jitter_window: usize,

    /// Path to write debug SAM file with seed-level weighted interval
    /// scheduling (WIS) results. Each seed is a separate record, tagged
    /// with XE (explanation index) and XS (segment index within explanation).
    /// Group by XS in IGV to see how seeds form segments.
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_wis_sam: String,

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

    /// Path to write debug TSV file with SV-break spans (seed runs flanked by
    /// SV breaks on both sides, where the chain returns to a colinear seed).
    /// Columns: read_name, anchor_before_read_start, anchor_before_read_end,
    ///          chrom, anchor_before_ref_start, anchor_before_ref_end,
    ///          num_sv_breaks, anchor_after_read_start, anchor_after_read_end,
    ///          anchor_after_ref_start, anchor_after_ref_end,
    ///          strand, read_gap, ref_gap
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_sv_spans_tsv: String,

    /// Path to write a FASTQ file containing reads that fail alignment validation.
    /// Each failing read is written once (even if multiple segments fail).
    /// Leave empty to disable.
    #[config(default = "")]
    pub debug_failed_reads_fastq: String,

    /// Fixed penalty applied when chaining two seeds across a structural
    /// variant boundary (different chromosome, different strand, or
    /// non-colinear reference order).  Higher values make the chaining DP
    /// less willing to bridge SVs within a single chain.
    #[config(default = 200.0)]
    pub sv_penalty: f64,

    /// Maximum reference-space distance between two seeds (on the same
    /// chromosome and strand) for a backward ref jump to be treated as a
    /// tandem repeat traversal rather than a genuine SV.
    ///
    /// When the next seed's ref position steps backward (non-colinear), the
    /// DP normally applies `sv_penalty`.  But if both seeds land within
    /// this many bp of each other on the reference, the backward jump is most
    /// likely a tandem repeat expansion: the read contains extra copies of a
    /// short repeat unit that are all anchored to the same narrow ref window.
    /// In that case a logarithmic penalty proportional to the size of the
    /// backward step is used instead of `sv_penalty`, allowing the expansion
    /// seeds to chain through naturally without leaving a large unanchored
    /// read gap.
    ///
    /// Set to 0 to disable (all backward jumps on the same chrom/strand use
    /// `sv_penalty`).
    #[config(default = 400)]
    pub repeat_expansion_max_ref_window: usize,

    /// Fixed additive penalty applied on top of the normal gap cost when
    /// chaining across a tandem repeat traversal (a backward ref jump within
    /// `repeat_expansion_max_ref_window`).  This ensures a backward step of
    /// size `d` always costs more than a forward gap of the same size by
    /// exactly this amount.  Set to 0 to make forward and backward gaps
    /// equally expensive.
    #[config(default = 120.0)]
    pub repeat_expansion_penalty: f64,

    /// Deviation threshold (bp) above which a linear penalty term is added
    /// to the chaining gap cost.  Below this the cost is purely logarithmic,
    /// which is cheap for small insertions and deletions.  Above it a linear component kicks in,
    /// making repeat-copy hops and large gaps increasingly expensive.
    /// Set to 0 to disable (pure logarithmic, legacy behaviour).
    /// The SNV/SV boundary (~50 bp) is a natural choice.
    #[config(default = 50.0)]
    pub gap_linear_threshold: f64,

    /// Scaling factor for the linear penalty term applied to deviation above
    /// `gap_linear_threshold`.  The full gap penalty is:
    ///   ln(1 + min(deviation, threshold)) + k * max(deviation - threshold, 0)
    /// Higher values suppress tandem-repeat copy hops more aggressively but
    /// also make the DP less willing to chain across large genuine deletions.
    /// Set to 0.0 to disable (pure logarithmic, legacy behaviour).
    #[config(default = 0.15)]
    pub gap_linear_scale: f64,

    /// Quadratic scaling factor for the read-gap cost component of the edge
    /// penalty.  The read-gap cost becomes:
    ///   read_gap + read_gap_quad_scale * read_gap²
    /// This makes long unanchored stretches of read disproportionately
    /// expensive, discouraging the DP from chaining seeds that leave large
    /// portions of the read unexplained.  Any sequence genuinely present in
    /// the read but absent from the local reference region is better
    /// represented as an SV breakpoint than as a colinear gap.
    /// Set to 0.0 to disable (linear read-gap cost, legacy behaviour).
    #[config(default = 0.025)]
    pub read_gap_quad_scale: f64,

    /// Maximum identity ratio (worse/better) for merging a ref-overlapping
    /// segment pair by converting the poorly-aligning overlap into an INS.
    ///
    /// When two adjacent segments have a reference overlap, the alignment
    /// identity of each segment within that overlap region is computed as:
    ///   matches / max(read_bases_in_overlap, ref_bases_in_overlap)
    /// A merge fires when the worse identity is below this fraction of the
    /// better identity.  The worse-aligning end is then trimmed to the overlap
    /// boundary, the trimmed query bases become an INS, and the two segments
    /// are merged into one.
    ///
    /// 0.0 disables merging; 1.0 merges any overlapping pair regardless of
    /// quality.  The default 0.5 merges when one side has less than half
    /// of the identity of the other in the overlap region.
    #[config(default = 0.5)]
    pub overlap_merge_max_identity_ratio: f64,

    /// Minimum overlap size (bp) below which segments are always merged,
    /// regardless of identity.  Small overlaps are likely due to seed
    /// placement granularity rather than genuine ambiguity, so forcing a
    /// merge avoids spurious supplementary alignments.
    /// Set to 0 to disable (always require identity check).
    #[config(default = 4000)]
    pub overlap_merge_min_forced: usize,

    /// Use batched prefetching for seed lookups.
    ///
    /// When true, syncmer k-mers are collected into a batch buffer first,
    /// then looked up with software-pipelined prefetching to hide memory
    /// latency in the multi-GB hash tables. This can significantly improve
    /// throughput on large indices.
    #[config(default = true)]
    pub batch_prefetch: bool,

    /// Use collinearity-weighted seed scoring and isolated-seed pre-pruning.
    ///
    /// When true, a collinearity weight c(x) = Σ_y 1/(1 + (diag_x - diag_y)²)
    /// is computed for each seed (summed over seeds on the same chrom/strand
    /// within a diagonal window), and the seed weight becomes:
    ///   length * collinearity / sqrt(kmer_frequency)
    /// Seeds with no colinear neighbours (collinearity ≈ 1.0) are pruned before
    /// the DP, dramatically reducing its O(n²) cost on repetitive reads.
    /// The edge penalty also switches to a linear read-gap model with overlap
    /// truncation instead of a hard read-overlap cutoff.
    ///
    /// When false (default), the original weight formula and edge penalty are used.
    #[config(default = true)]
    pub use_collinearity_weights: bool,

    /// Use the break-count DP (v3) for chaining.
    ///
    /// When true, the chaining DP tracks the number of SV breaks taken so far
    /// as part of its state: `dp[seed][k]` is the best score of a chain ending
    /// at `seed` having made exactly `k` SV breaks.  Each successive SV break
    /// costs one additional `sv_penalty` on top of the base cost, making the
    /// total penalty for `k` breaks `k*(k+1)/2 * sv_penalty` — quadratic in k.
    /// This gives the DP optimal substructure for the objective of maximising
    /// anchored read length while strongly discouraging unnecessary SV breaks,
    /// and prevents segmental-duplication seeds from accumulating enough weight
    /// to justify a double jump away from and back to the primary diagonal.
    ///
    /// Requires `use_collinearity_weights = true`; has no effect otherwise.
    /// The maximum number of SV breaks tracked is controlled by
    /// `max_sv_breaks`.
    #[config(default = true)]
    pub use_break_count_dp: bool,

    /// Maximum number of SV breaks to track in the break-count DP (v3).
    /// The DP state space is O(n * max_sv_breaks), so keeping this small
    /// (4–6) is important for performance.  Chains requiring more breaks than
    /// this are still found but treated as if they used exactly `max_sv_breaks`
    /// breaks (i.e. the penalty stops escalating beyond this point).
    /// Only used when `use_break_count_dp = true`.
    #[config(default = 4)]
    pub max_sv_breaks: usize,

    /// Diagonal window (bp) for collinearity weight computation.
    /// Seeds further apart than this on the diagonal are treated as unrelated.
    /// Only used when `use_collinearity_weights = true`.
    #[config(default = 50.0)]
    pub collinearity_diagonal_cutoff: f64,

    /// Maximum deviation (bp) between ref gap and read gap for a Continuation
    /// edge under the collinearity model. Only used when
    /// `use_collinearity_weights = true`.
    #[config(default = 1000.0)]
    pub collinearity_max_gap_deviation: f64,
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
    pub min_aligned_length: usize,
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
    #[config(default = false)]
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
            mismatch: 4,
            gap_open: 6,
            gap_extend: 2,
        }
    }
}

/// Initialize the global configuration (call once at startup).
///
/// Metrics and histogram configuration.
#[derive(Config, Debug, Clone)]
pub struct MetricsConfig {
    /// Path to write the metrics summary TSV file.
    #[config(default = "parallax-stats.tsv")]
    pub stats_path: String,

    /// Logging interval for progress output (seconds).
    #[config(default = 30.0)]
    pub logging_interval: f64,
}

/// # Panics
/// Panics if called more than once.
pub fn init(config: ParallaxConfig) {
    CONFIG
        .set(config)
        .expect("Config already initialized - init() called twice");
}

/// Load configuration from a TOML file, falling back to defaults.
///
/// Resolution order:
/// 1. Explicit `path` if supplied
/// 2. `.parallax.toml` in the current directory, if it exists
/// 3. `$HOME/.parallax.toml`, if it exists
/// 4. Built-in defaults
pub fn load(path: Option<&std::path::Path>) -> Result<ParallaxConfig, confique::Error> {
    if let Some(p) = path {
        return ParallaxConfig::builder().file(p).load();
    }

    let cwd_config = std::path::Path::new(".parallax.toml");
    if cwd_config.exists() {
        log::info!("Using config: {}", cwd_config.display());
        return ParallaxConfig::builder().file(cwd_config).load();
    }

    if let Some(home) = std::env::var_os("HOME") {
        let home_config = std::path::PathBuf::from(home).join(".parallax.toml");
        if home_config.exists() {
            log::info!("Using config: {}", home_config.display());
            return ParallaxConfig::builder().file(home_config).load();
        }
    }

    ParallaxConfig::builder().load()
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
        assert_eq!(config.seeding.max_seed_occurrences, 500);
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
