use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::record_buf::RecordBuf;

/// Output format for alignment records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Sam,
    Bam,
    Cram,
}

impl OutputFormat {
    /// Detect format from file extension.
    pub fn from_path(path: &Path) -> Option<Self> {
        match path.extension()?.to_str()? {
            "sam" => Some(Self::Sam),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            _ => None,
        }
    }
}

impl std::str::FromStr for OutputFormat {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s.to_lowercase().as_str() {
            "sam" => Ok(Self::Sam),
            "bam" => Ok(Self::Bam),
            "cram" => Ok(Self::Cram),
            _ => Err(format!(
                "unknown output format: '{}' (expected sam, bam, or cram)",
                s
            )),
        }
    }
}

impl std::fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Sam => write!(f, "SAM"),
            Self::Bam => write!(f, "BAM"),
            Self::Cram => write!(f, "CRAM"),
        }
    }
}

/// Internal enum wrapping format-specific noodles writers.
enum FormatWriter {
    Sam(sam::io::Writer<BufWriter<Box<dyn Write + Send>>>),
    Bam(noodles::bam::io::Writer<noodles::bgzf::io::Writer<Box<dyn Write + Send>>>),
    Cram(noodles::cram::io::Writer<Box<dyn Write + Send>>),
}

/// Builder for creating an AlignmentWriter with proper headers.
///
/// Accumulates header information (@SQ, @RG, @PG) before creating the writer.
/// When `build()` is called, headers are written in the appropriate format.
pub struct AlignmentWriterBuilder {
    output: Box<dyn Write + Send>,
    format: OutputFormat,
    reference_repository: fasta::Repository,
    contigs: Vec<(String, usize)>,
    read_group: Option<String>,
    command_line: Option<String>,
}

impl AlignmentWriterBuilder {
    /// Create a new builder that will write to the given output in the specified format.
    ///
    /// The `reference_repository` is required for CRAM output (the CRAM encoder
    /// needs the reference sequences). For SAM/BAM it is stored but unused.
    pub fn new(
        output: Box<dyn Write + Send>,
        format: OutputFormat,
        reference_repository: fasta::Repository,
    ) -> Self {
        Self {
            output,
            format,
            reference_repository,
            contigs: Vec::new(),
            read_group: None,
            command_line: None,
        }
    }

    /// Add a contig (reference sequence) to the header.
    pub fn add_contig(mut self, name: &str, length: usize) -> Self {
        self.contigs.push((name.to_string(), length));
        self
    }

    /// Add multiple contigs from an iterator of (name, length) pairs.
    pub fn add_contigs<'a, I>(mut self, contigs: I) -> Self
    where
        I: IntoIterator<Item = (&'a str, u64)>,
    {
        for (name, length) in contigs {
            self = self.add_contig(name, length as usize);
        }
        self
    }

    /// Set the read group header line.
    ///
    /// The `rg_line` should be a complete @RG line (e.g., from ReadGroup::to_header_line()).
    pub fn read_group(mut self, rg_line: Option<String>) -> Self {
        self.read_group = rg_line;
        self
    }

    /// Set the command line for the @PG header.
    pub fn command_line(mut self, cmd: &str) -> Self {
        self.command_line = Some(cmd.to_string());
        self
    }

    /// Build the AlignmentWriter, writing all headers to the output.
    ///
    /// Constructs SAM header text, parses it into a noodles Header, then
    /// writes it using the format-specific writer. This ensures @HD, @SQ, @RG,
    /// and @PG records are present in all output formats.
    pub fn build(self) -> std::io::Result<AlignmentWriter> {
        // Construct SAM header text, then parse into a noodles Header.
        let mut header_text = String::new();

        // @HD
        header_text.push_str("@HD\tVN:1.6\tSO:unsorted\n");

        // @SQ
        for (name, length) in &self.contigs {
            header_text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, length));
        }

        // @RG
        if let Some(ref rg) = self.read_group {
            header_text.push_str(rg);
            header_text.push('\n');
        }

        // @PG
        let version = format!("{}+{}", env!("CARGO_PKG_VERSION"), env!("GIT_VERSION"));
        if let Some(ref cmd) = self.command_line {
            header_text.push_str(&format!(
                "@PG\tID:parallax\tPN:parallax\tVN:{}\tCL:{}\n",
                version, cmd
            ));
        } else {
            header_text.push_str(&format!("@PG\tID:parallax\tPN:parallax\tVN:{}\n", version));
        }

        let header: sam::Header = header_text.parse().map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("failed to parse SAM header: {}", e),
            )
        })?;

        // Create format-specific writer and write header
        let format_writer = match self.format {
            OutputFormat::Sam => {
                let mut w = sam::io::Writer::new(BufWriter::new(self.output));
                w.write_header(&header)?;
                FormatWriter::Sam(w)
            }
            OutputFormat::Bam => {
                let mut w = noodles::bam::io::Writer::new(self.output);
                w.write_header(&header)?;
                FormatWriter::Bam(w)
            }
            OutputFormat::Cram => {
                let mut w = noodles::cram::io::writer::Builder::default()
                    .set_reference_sequence_repository(self.reference_repository)
                    .build_from_writer(self.output);
                w.write_header(&header)?;
                FormatWriter::Cram(w)
            }
        };

        let now = std::time::Instant::now();
        let empty = ProgressSnapshot { records: 0, bases: 0, time: now };
        Ok(AlignmentWriter {
            header,
            inner: Mutex::new(format_writer),
            counter: std::sync::atomic::AtomicUsize::new(0),
            bases_written: std::sync::atomic::AtomicU64::new(0),
            start_time: now,
            recent: Mutex::new([empty; AlignmentWriter::WINDOW]),
            recent_pos: std::sync::atomic::AtomicUsize::new(0),
        })
    }
}

/// Snapshot of cumulative progress at a single report point.
#[derive(Clone, Copy)]
struct ProgressSnapshot {
    records: usize,
    bases: u64,
    time: std::time::Instant,
}

/// Thread-safe alignment writer supporting SAM, BAM, and CRAM output.
///
/// Each write operation is atomic — the entire record is written while
/// holding a lock, preventing interleaved output from different threads.
///
/// Create using `AlignmentWriterBuilder` to ensure headers are written first.
pub struct AlignmentWriter {
    header: sam::Header,
    inner: Mutex<FormatWriter>,
    counter: std::sync::atomic::AtomicUsize,
    bases_written: std::sync::atomic::AtomicU64,
    start_time: std::time::Instant,
    // Ring buffer of recent progress snapshots for windowed rate calculation.
    recent: Mutex<[ProgressSnapshot; AlignmentWriter::WINDOW]>,
    recent_pos: std::sync::atomic::AtomicUsize,
}

impl AlignmentWriter {
    /// Number of progress snapshots kept for windowed rate calculation.
    const WINDOW: usize = 8;

    /// Create a builder for constructing an AlignmentWriter with headers.
    pub fn builder(
        output: Box<dyn Write + Send>,
        format: OutputFormat,
        reference_repository: fasta::Repository,
    ) -> AlignmentWriterBuilder {
        AlignmentWriterBuilder::new(output, format, reference_repository)
    }

    /// Returns a reference to the SAM header.
    #[allow(dead_code)]
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Write an alignment record atomically.
    ///
    /// The record is written through the format-specific noodles writer
    /// while holding the lock, ensuring thread safety.
    pub fn write_record(&self, record: &RecordBuf) -> std::io::Result<()> {
        use noodles::sam::alignment::io::Write as _;
        let read_len = record.sequence().len() as u64;
        let mut inner = self.inner.lock().unwrap();
        match &mut *inner {
            FormatWriter::Sam(w) => w.write_alignment_record(&self.header, record),
            FormatWriter::Bam(w) => w.write_alignment_record(&self.header, record),
            FormatWriter::Cram(w) => w.write_alignment_record(&self.header, record),
        }?;

        let n = self.counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let b = self.bases_written.fetch_add(read_len, std::sync::atomic::Ordering::Relaxed);
        if n & 1023 == 0 {
            let now = std::time::Instant::now();
            let snap = ProgressSnapshot { records: n, bases: b, time: now };
            let pos = self.recent_pos.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let (recent_records, recent_bases, recent_secs) = {
                let mut ring = self.recent.lock().unwrap();
                let oldest = ring[pos % Self::WINDOW];
                ring[(pos + 1) % Self::WINDOW] = snap;
                let dr = n.saturating_sub(oldest.records);
                let db = b.saturating_sub(oldest.bases);
                let dt = (now - oldest.time).as_secs_f64().max(1e-6);
                (dr, db, dt)
            };
            let elapsed = (now - self.start_time).as_secs_f64();
            log::info!(
                "Written {} records in {:.0}s [{:.0} rec/s, {:.0} kbp/s]",
                n,
                elapsed,
                recent_records as f64 / recent_secs,
                recent_bases as f64 / 1000.0 / recent_secs,
            );
        }
        Ok(())
    }

    /// Finish the output stream, writing any pending data and format-specific
    /// EOF markers.
    ///
    /// For BAM, this flushes pending data and writes the BGZF EOF block.
    /// For CRAM, this flushes pending containers and writes the EOF container.
    /// For SAM, this simply flushes the buffer.
    ///
    /// Must be called before dropping the writer to ensure valid output.
    pub fn finish(&self) -> std::io::Result<()> {
        let mut inner = self.inner.lock().unwrap();
        match &mut *inner {
            FormatWriter::Sam(w) => w.get_mut().flush(),
            FormatWriter::Bam(w) => w.try_finish(),
            FormatWriter::Cram(w) => w.try_finish(&self.header),
        }?;

        let n = self.counter.load(std::sync::atomic::Ordering::Relaxed);
        let b = self.bases_written.load(std::sync::atomic::Ordering::Relaxed);
        let elapsed = self.start_time.elapsed().as_secs_f64();
        log::info!(
            "Written {} records in {:.0}s [{:.0} rec/s, {:.0} kbp/s]",
            n,
            elapsed,
            n as f64 / elapsed,
            b as f64 / 1000.0 / elapsed,
        );
        Ok(())
    }
}
