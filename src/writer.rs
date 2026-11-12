use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::record_buf::RecordBuf;
use parallax::config;
use parallax::utils::progress::{RateProgress, RateProgressConfig};

pub trait RecordWriter: Send + Sync {
    fn write_record(&self, record: &RecordBuf) -> std::io::Result<()>;
    fn finish(&self) -> std::io::Result<()>;
}

impl<T: RecordWriter> RecordWriter for std::sync::Arc<T> {
    fn write_record(&self, record: &RecordBuf) -> std::io::Result<()> {
        (**self).write_record(record)
    }

    fn finish(&self) -> std::io::Result<()> {
        (**self).finish()
    }
}

pub mod bam_writer;
pub mod sam_writer;
pub mod sorting_writer;

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

        let base_interval = config::get().metrics.logging_interval;

        let progress = RateProgress::with_config(
            RateProgressConfig::default()
                .with_item("alignments")
                .with_unit("bp")
                .with_interval(base_interval * 0.99),
        );
        Ok(AlignmentWriter {
            header,
            inner: Mutex::new(AlignmentWriterInner {
                writer: format_writer,
                progress,
            }),
        })
    }
}

/// Thread-safe alignment writer supporting SAM, BAM, and CRAM output.
///
/// Each write operation is atomic — the entire record is written while
/// holding a lock, preventing interleaved output from different threads.
///
/// Create using `AlignmentWriterBuilder` to ensure headers are written first.
pub struct AlignmentWriter {
    header: sam::Header,
    inner: Mutex<AlignmentWriterInner>,
}

impl AlignmentWriter {
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
}

impl RecordWriter for AlignmentWriter {
    fn write_record(&self, record: &RecordBuf) -> std::io::Result<()> {
        use noodles::sam::alignment::io::Write as _;
        let read_len = record.sequence().len();
        let mut inner = self.inner.lock().unwrap();
        match &mut *&mut inner.writer {
            FormatWriter::Sam(w) => w.write_alignment_record(&self.header, record),
            FormatWriter::Bam(w) => w.write_alignment_record(&self.header, record),
            FormatWriter::Cram(w) => w.write_alignment_record(&self.header, record),
        }?;

        inner.progress.record(read_len as u64);

        Ok(())
    }

    fn finish(&self) -> std::io::Result<()> {
        let mut inner = self.inner.lock().unwrap();
        match &mut *&mut inner.writer {
            FormatWriter::Sam(w) => w.get_mut().flush(),
            FormatWriter::Bam(w) => w.try_finish(),
            FormatWriter::Cram(w) => w.try_finish(&self.header),
        }?;

        inner.progress.finish();

        Ok(())
    }
}

struct AlignmentWriterInner {
    writer: FormatWriter,
    progress: RateProgress,
}
