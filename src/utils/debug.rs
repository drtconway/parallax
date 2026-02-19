//! Debug/diagnostic file infrastructure.
//!
//! Provides a trait-based system for writing debug output files during alignment.
//! Each debug file is a static singleton in the module where it's used, lazily
//! initialized on first access from configuration.
//!
//! # Architecture
//!
//! - [`DebugOutput`] — trait with `create() -> Option<Self>` for self-describing
//!   debug output destinations (TSV, SAM/BAM, text)
//! - [`DebugFile<T>`] — lazy wrapper that calls `T::create()` on first use
//! - Concrete types live in the modules that use them, each encapsulating
//!   their own config lookup, header, and file path
//!
//! # Usage
//!
//! ```ignore
//! // 1. Define a static in the module that uses it:
//! static CHAINS_TSV: DebugFile<ChainsTsvDebug> = DebugFile::new();
//!
//! // 2. Use it (lazy init happens automatically on first access):
//! if CHAINS_TSV.is_enabled() {
//!     CHAINS_TSV.append(&format!("{}\t{}", read_name, cluster_id));
//! }
//!
//! // 3. At shutdown (once, from main):
//! DebugFile::finish_all();
//! ```

use std::{
    fmt,
    fs::File,
    io::{self, BufWriter, Write},
    sync::{Mutex, OnceLock},
};

// ── Finisher registry ────────────────────────────────────────────────────────

/// Object-safe trait for finalizing debug output files.
trait Finishable: Send + Sync {
    fn finish(&self);
}

impl<T: DebugOutput> Finishable for DebugFile<T> {
    fn finish(&self) {
        if let Some(w) = self.inner.get().and_then(|o| o.as_ref()) {
            w.finish();
        }
    }
}

/// Global registry of all debug files that need to be finished at shutdown.
static FINISHER_REGISTRY: Mutex<Vec<&'static dyn Finishable>> = Mutex::new(Vec::new());

fn register_finisher(f: &'static dyn Finishable) {
    if let Ok(mut reg) = FINISHER_REGISTRY.lock() {
        reg.push(f);
    }
}

// ── Global chromosome info (for SAM headers) ────────────────────────────────

/// Chromosome names and lengths, set once at startup for SAM debug headers.
static CHROMOSOMES: OnceLock<Vec<(String, u64)>> = OnceLock::new();

/// Store the reference chromosome info for debug SAM header generation.
///
/// Call once at startup after the reference is loaded. Idempotent — only the
/// first call has effect.
pub fn set_reference_info<'a>(chroms: impl Iterator<Item = (&'a str, u64)>) {
    let _ = CHROMOSOMES.set(chroms.map(|(n, l)| (n.to_string(), l)).collect());
}

/// Build a SAM header from the stored chromosome info, or `None` if not yet set.
pub fn sam_header() -> Option<String> {
    CHROMOSOMES.get().map(|chroms| build_sam_header(chroms.iter().map(|(n, l)| (n.as_str(), *l))))
}

// ── Trait ────────────────────────────────────────────────────────────────────

/// Trait for debug output destinations.
///
/// Each implementation encapsulates its own config lookup and file creation
/// logic in [`create()`](Self::create). Thread safety is the implementation's
/// responsibility (typically via an internal `Mutex`).
pub trait DebugOutput: Send + Sync + Sized {
    /// The type of item this output accepts.
    ///
    /// Use a specific tuple for TSV outputs (compile-time field count safety),
    /// or `str` for free-form text/SAM outputs.
    type Item<'a>: ?Sized;

    /// Attempt to create this debug output.
    ///
    /// Returns `Some(Self)` if the output is enabled (config path is non-empty
    /// and the file can be opened), or `None` if disabled.
    fn create() -> Option<Self>;

    /// Write a single item to the output.
    fn append(&self, item: &Self::Item<'_>);

    /// Flush and finalize the output.
    fn finish(&self);
}

// ── DebugFile wrapper ────────────────────────────────────────────────────────

/// A lazily-initialized debug output file.
///
/// Wraps an `Option<T>` in a [`OnceLock`], calling [`T::create()`](DebugOutput::create)
/// on first access. If `create()` returns `None`, all subsequent operations are
/// no-ops. Auto-registers with the finisher registry when enabled.
///
/// This is intended to be used as a `static` in the module that owns the debug output.
pub struct DebugFile<T: DebugOutput> {
    inner: OnceLock<Option<T>>,
}

impl<T: DebugOutput> DebugFile<T> {
    /// Create an empty (not yet initialized) debug file slot.
    pub const fn new() -> Self {
        Self {
            inner: OnceLock::new(),
        }
    }

    /// Lazily initialize, returning a reference to the writer if enabled.
    fn get_or_init(&'static self) -> Option<&'static T> {
        self.inner
            .get_or_init(|| {
                let result = T::create();
                if result.is_some() {
                    register_finisher(self);
                }
                result
            })
            .as_ref()
    }

    /// Check if this debug file is enabled.
    ///
    /// Triggers lazy initialization on first call.
    /// Use this to guard expensive formatting before calling [`append()`](Self::append).
    pub fn is_enabled(&'static self) -> bool {
        self.get_or_init().is_some()
    }

    /// Write an item. No-op if disabled.
    ///
    /// Triggers lazy initialization on first call.
    pub fn append<'a>(&'static self, item: &T::Item<'a>) {
        if let Some(w) = self.get_or_init() {
            w.append(item);
        }
    }

    /// Flush and finalize all registered debug files.
    ///
    /// Call once at shutdown to ensure all buffered data is written.
    pub fn finish_all() {
        if let Ok(reg) = FINISHER_REGISTRY.lock() {
            for f in reg.iter() {
                f.finish();
            }
        }
    }
}

// ── TsvRow trait ─────────────────────────────────────────────────────────

/// A type that can be written as a tab-separated row.
///
/// Implemented for tuples of [`Display`](fmt::Display) types via macro.
/// Use with [`DebugTsvWriter::append_row`] for type-safe TSV output.
pub trait TsvRow {
    /// Number of fields in this row (must match the header count).
    const NUM_FIELDS: usize;

    /// Write tab-separated fields followed by a newline.
    fn write_row(&self, w: &mut dyn Write) -> io::Result<()>;
}

macro_rules! impl_tsv_row {
    ($n:expr; $first_idx:tt: $first_ty:ident $(, $idx:tt: $T:ident)*) => {
        impl<$first_ty: fmt::Display $(, $T: fmt::Display)*> TsvRow for ($first_ty, $($T,)*) {
            const NUM_FIELDS: usize = $n;
            fn write_row(&self, w: &mut dyn Write) -> io::Result<()> {
                write!(w, "{}", self.$first_idx)?;
                $(write!(w, "\t{}", self.$idx)?;)*
                writeln!(w)
            }
        }
    };
}

impl_tsv_row!( 1; 0: A);
impl_tsv_row!( 2; 0: A, 1: B);
impl_tsv_row!( 3; 0: A, 1: B, 2: C);
impl_tsv_row!( 4; 0: A, 1: B, 2: C, 3: D);
impl_tsv_row!( 5; 0: A, 1: B, 2: C, 3: D, 4: E);
impl_tsv_row!( 6; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F);
impl_tsv_row!( 7; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F, 6: G);
impl_tsv_row!( 8; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F, 6: G, 7: H);
impl_tsv_row!( 9; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F, 6: G, 7: H, 8: I);
impl_tsv_row!(10; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F, 6: G, 7: H, 8: I, 9: J);
impl_tsv_row!(11; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F, 6: G, 7: H, 8: I, 9: J, 10: K);
impl_tsv_row!(12; 0: A, 1: B, 2: C, 3: D, 4: E, 5: F, 6: G, 7: H, 8: I, 9: J, 10: K, 11: L);

// ── DebugTsvWriter ────────────────────────────────────────────────────────

/// A text/TSV line writer with internal mutex.
///
/// Used by concrete debug types as their internal storage. Not a `DebugOutput`
/// itself — the concrete types wrap this and provide their own `create()`.
pub struct DebugTsvWriter {
    writer: Mutex<BufWriter<File>>,
}

impl DebugTsvWriter {
    pub fn open(path: &str, header: Option<&str>) -> io::Result<Self> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        if let Some(h) = header {
            writeln!(writer, "{}", h)?;
        }
        log::info!("Debug output enabled: {}", path);
        Ok(Self {
            writer: Mutex::new(writer),
        })
    }

    pub fn append(&self, line: &str) {
        if let Ok(mut w) = self.writer.lock() {
            let _ = writeln!(w, "{}", line);
        }
    }

    pub fn finish(&self) {
        if let Ok(mut w) = self.writer.lock() {
            let _ = w.flush();
        }
    }

    /// Write a typed row as tab-separated values.
    pub fn append_row(&self, row: &impl TsvRow) {
        if let Ok(mut w) = self.writer.lock() {
            let _ = row.write_row(&mut *w);
        }
    }
}

// ── SAM header helper ────────────────────────────────────────────────────────

/// Build a text SAM header from reference chromosomes.
pub fn build_sam_header<'a>(chromosomes: impl Iterator<Item = (&'a str, u64)>) -> String {
    let mut header = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, length) in chromosomes {
        header.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, length));
    }
    let version = format!("{}+{}", env!("CARGO_PKG_VERSION"), env!("GIT_VERSION"));
    header.push_str(&format!(
        "@PG\tID:parallax\tPN:parallax\tVN:{}\tCL:debug",
        version
    ));
    header
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    /// A test-only debug type that is always disabled.
    struct DisabledDebug;
    impl DebugOutput for DisabledDebug {
        type Item<'a> = str;
        fn create() -> Option<Self> { None }
        fn append(&self, _: &str) {}
        fn finish(&self) {}
    }

    static TEST_DISABLED: DebugFile<DisabledDebug> = DebugFile::new();

    #[test]
    fn test_debug_file_disabled() {
        assert!(!TEST_DISABLED.is_enabled());
        TEST_DISABLED.append("should be ignored");
    }

    #[test]
    fn test_tsv_writer_directly() {
        let temp = NamedTempFile::new().unwrap();
        let path = temp.path().to_str().unwrap();

        let writer = DebugTsvWriter::open(path, Some("col1\tcol2")).unwrap();
        writer.append("value1\tvalue2");
        writer.finish();

        let mut content = String::new();
        std::fs::File::open(path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();

        assert!(content.contains("col1\tcol2"));
        assert!(content.contains("value1\tvalue2"));
    }
}
