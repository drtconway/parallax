//! Unified debug/tracing file infrastructure.
//!
//! Provides a centralized system for writing debug output files during alignment.
//! All debug files are registered at startup and can be written to from anywhere
//! in the codebase.
//!
//! # Usage
//!
//! ```ignore
//! // In main/initialization:
//! debug::init(&config, &reference)?;
//!
//! // Throughout the code:
//! debug::write(DebugFile::Chains, &format!("{}\t{}\t...", read_name, cluster_id));
//!
//! // At shutdown:
//! debug::flush_all();
//! ```

use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    sync::{Mutex, OnceLock},
};

use crate::config::ParallaxConfig;
use crate::reference::InMemoryReference;

/// Identifies a specific debug output file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DebugFile {
    /// Debug SAM file with extended seeds (before clustering)
    Seeds,
    /// TSV file with candidate seeds
    SeedsTsv,
    /// TSV file with seeds grouped into clusters (before chaining)
    ClustersTsv,
    /// TSV file with seed chains/clusters (after chaining, before alignment)
    Chains,
    /// SAM file with seed chains linked via SA tags (for IGV visualization)
    ChainsSam,
    /// A TSV file with gaps and potential fills
    GapFills,
    /// Failed alignment strings
    GapAlignments,
}

impl DebugFile {
    /// Get all variants for iteration
    #[allow(dead_code)]
    pub fn all() -> &'static [DebugFile] {
        &[
            DebugFile::Seeds,
            DebugFile::SeedsTsv,
            DebugFile::ClustersTsv,
            DebugFile::Chains,
            DebugFile::ChainsSam,
            DebugFile::GapFills,
            DebugFile::GapAlignments,
        ]
    }

    /// Human-readable name for logging
    pub fn name(&self) -> &'static str {
        match self {
            DebugFile::Seeds => "seeds SAM",
            DebugFile::SeedsTsv => "seeds TSV",
            DebugFile::ClustersTsv => "clusters TSV",
            DebugFile::Chains => "chains TSV",
            DebugFile::ChainsSam => "chains SAM",
            DebugFile::GapFills => "gap fills TSV",
            DebugFile::GapAlignments => "gap alignments",
        }
    }
}

/// Holds an open debug file writer
struct DebugFileWriter {
    writer: BufWriter<File>,
    #[allow(dead_code)]
    path: String,
}

/// Global registry of debug file writers
struct DebugRegistry {
    files: HashMap<DebugFile, DebugFileWriter>,
}

impl DebugRegistry {
    fn new() -> Self {
        Self {
            files: HashMap::new(),
        }
    }

    fn register(
        &mut self,
        kind: DebugFile,
        path: &str,
        header: Option<&str>,
    ) -> std::io::Result<()> {
        if path.is_empty() {
            return Ok(());
        }

        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        if let Some(h) = header {
            writeln!(writer, "{}", h)?;
        }

        self.files.insert(
            kind,
            DebugFileWriter {
                writer,
                path: path.to_string(),
            },
        );

        log::info!("Debug {} output enabled: {}", kind.name(), path);
        Ok(())
    }

    fn write(&mut self, kind: DebugFile, line: &str) {
        if let Some(file) = self.files.get_mut(&kind) {
            let _ = writeln!(file.writer, "{}", line);
        }
    }

    fn write_raw(&mut self, kind: DebugFile, data: &str) {
        if let Some(file) = self.files.get_mut(&kind) {
            let _ = write!(file.writer, "{}", data);
        }
    }

    fn is_enabled(&self, kind: DebugFile) -> bool {
        self.files.contains_key(&kind)
    }

    fn flush(&mut self, kind: DebugFile) {
        if let Some(file) = self.files.get_mut(&kind) {
            let _ = file.writer.flush();
        }
    }

    fn flush_all(&mut self) {
        for file in self.files.values_mut() {
            let _ = file.writer.flush();
        }
    }
}

/// Global debug registry, initialized once
static DEBUG_REGISTRY: OnceLock<Mutex<DebugRegistry>> = OnceLock::new();

/// Get or create the global debug registry
fn registry() -> &'static Mutex<DebugRegistry> {
    DEBUG_REGISTRY.get_or_init(|| Mutex::new(DebugRegistry::new()))
}

/// Register a debug file with the global registry.
///
/// This should be called during initialization for each debug file type.
/// If `path` is empty, the file is not registered (disabled).
///
/// # Arguments
/// * `kind` - The type of debug file
/// * `path` - Output file path (empty = disabled)
/// * `header` - Optional header line to write at the start
pub fn register(kind: DebugFile, path: &str, header: Option<&str>) -> std::io::Result<()> {
    let mut reg = registry().lock().unwrap();
    reg.register(kind, path, header)
}

/// Write a line to a debug file.
///
/// If the file was not registered or is disabled, this is a no-op.
/// The line will have a newline appended automatically.
pub fn write(kind: DebugFile, line: &str) {
    if let Ok(mut reg) = registry().try_lock() {
        reg.write(kind, line);
    }
}

/// Write raw data to a debug file without appending a newline.
///
/// Useful for multi-line output or when you need precise control over formatting.
#[allow(dead_code)]
pub fn write_raw(kind: DebugFile, data: &str) {
    if let Ok(mut reg) = registry().try_lock() {
        reg.write_raw(kind, data);
    }
}

/// Check if a debug file is enabled.
///
/// Useful for avoiding expensive formatting when debug output is disabled.
pub fn is_enabled(kind: DebugFile) -> bool {
    if let Ok(reg) = registry().try_lock() {
        reg.is_enabled(kind)
    } else {
        false
    }
}

/// Flush a specific debug file to disk.
#[allow(dead_code)]
pub fn flush(kind: DebugFile) {
    if let Ok(mut reg) = registry().try_lock() {
        reg.flush(kind);
    }
}

/// Flush all debug files to disk.
///
/// Should be called at program shutdown to ensure all data is written.
pub fn flush_all() {
    if let Ok(mut reg) = registry().try_lock() {
        reg.flush_all();
    }
}

/// Standard TSV headers for each debug file type
pub mod headers {
    /// Header for alignments TSV
    pub const ALIGNMENTS: &str =
        "read_name\tread_start\tread_end\tread_len\tchrom\tref_start\tref_end\tstrand\tscore";

    /// Header for clusters TSV (seeds with cluster index, before chaining)
    pub const CLUSTERS: &str =
        "read_name\tcluster_id\tread_start\tread_end\tread_len\tchrom\tref_start\tref_end\tstrand\tmatch_len";

    /// Header for chains TSV (one row per seed and one row per gap)
    pub const CHAINS: &str = "read_name\tcluster_id\trow_type\tread_start\tread_end\tread_width\tref_start\tref_end\tref_width\tchrom\tstrand\tuniqueness";

    pub const GAP_FILLS: &str = "read_name\tread_len\tread_start\tread_end\tfill_len\tcluster_idx\taln_score\tchrom_name\tref_start\tref_end\tstrand";

    //pub const GAP_ALIGNMENTS: &str = "source\tsequence";
}

/// Initialize all debug files from configuration.
///
/// This is the main entry point - call once at startup.
///
/// # Arguments
/// * `config` - The parallax configuration containing debug file paths
/// * `reference` - The reference genome (needed for SAM header)
pub fn init(config: &ParallaxConfig, reference: &InMemoryReference) -> std::io::Result<()> {
    let sam_header = if !config.seeding.debug_seeds_sam.is_empty()
        || !config.seeding.debug_chains_sam.is_empty()
    {
        Some(build_sam_header(reference.chromosomes()))
    } else {
        None
    };

    register(
        DebugFile::Seeds,
        &config.seeding.debug_seeds_sam,
        sam_header.as_deref(),
    )?;
    register(
        DebugFile::SeedsTsv,
        &config.seeding.debug_seeds_tsv,
        Some(headers::ALIGNMENTS),
    )?;
    register(
        DebugFile::ClustersTsv,
        &config.seeding.debug_clusters_tsv,
        Some(headers::CLUSTERS),
    )?;
    register(
        DebugFile::Chains,
        &config.seeding.debug_chains_tsv,
        Some(headers::CHAINS),
    )?;
    register(
        DebugFile::ChainsSam,
        &config.seeding.debug_chains_sam,
        sam_header.as_deref(),
    )?;
    register(
        DebugFile::GapFills,
        &config.seeding.debug_gap_fills_tsv,
        Some(headers::GAP_FILLS),
    )?;
    register(
        DebugFile::GapAlignments,
        &config.seeding.debug_gap_alignments,
        None,
    )?;
    Ok(())
}

/// Helper to build SAM header from reference chromosomes.
///
/// # Arguments
/// * `chromosomes` - Iterator of (name, length) pairs
pub fn build_sam_header<'a>(chromosomes: impl Iterator<Item = (&'a str, u64)>) -> String {
    let mut header = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, length) in chromosomes {
        header.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, length));
    }
    header.push_str("@PG\tID:parallax\tPN:parallax\tVN:0.1.0\tCL:debug");
    header
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_debug_file_disabled() {
        // Empty path should not create file
        register(DebugFile::Seeds, "", None).unwrap();
        assert!(!is_enabled(DebugFile::Seeds));
    }

    #[test]
    fn test_debug_file_write() {
        let temp = NamedTempFile::new().unwrap();
        let path = temp.path().to_str().unwrap();

        register(DebugFile::SeedsTsv, path, Some("header")).unwrap();
        assert!(is_enabled(DebugFile::SeedsTsv));

        write(DebugFile::SeedsTsv, "test line");
        flush(DebugFile::SeedsTsv);

        let mut content = String::new();
        std::fs::File::open(path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();

        assert!(content.contains("header"));
        assert!(content.contains("test line"));
    }
}
