use std::io::{BufWriter, Write};
use std::sync::Mutex;

/// Builder for creating an AlignmentWriter with proper SAM headers.
///
/// Accumulates header information (@SQ, @RG, @PG) before creating the writer.
/// When `build()` is called, all headers are written in the correct order,
/// and the resulting writer only handles alignment records.
pub struct AlignmentWriterBuilder<W: Write> {
    writer: BufWriter<W>,
    contigs: Vec<(String, usize)>,
    read_group: Option<String>,
    command_line: Option<String>,
}

impl<W: Write> AlignmentWriterBuilder<W> {
    /// Create a new builder that will write to the given output.
    pub fn new(writer: W) -> Self {
        Self {
            writer: BufWriter::new(writer),
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
    /// Headers are written in SAM-standard order: @HD, @SQ, @RG, @PG.
    pub fn build(mut self) -> std::io::Result<AlignmentWriter<W>> {
        // @HD - Header line
        writeln!(self.writer, "@HD\tVN:1.6\tSO:unsorted")?;

        // @SQ - Sequence dictionary
        for (name, length) in &self.contigs {
            writeln!(self.writer, "@SQ\tSN:{}\tLN:{}", name, length)?;
        }

        // @RG - Read group (if provided)
        if let Some(ref rg) = self.read_group {
            writeln!(self.writer, "{}", rg)?;
        }

        // @PG - Program record
        if let Some(ref cmd) = self.command_line {
            writeln!(
                self.writer,
                "@PG\tID:parallax\tPN:parallax\tVN:{}\tCL:{}",
                env!("CARGO_PKG_VERSION"),
                cmd
            )?;
        } else {
            writeln!(
                self.writer,
                "@PG\tID:parallax\tPN:parallax\tVN:{}",
                env!("CARGO_PKG_VERSION")
            )?;
        }

        Ok(AlignmentWriter {
            inner: Mutex::new(WriterInner {
                writer: self.writer,
            }),
        })
    }
}

/// Thread-safe alignment writer that supports multiple threads writing concurrently.
///
/// Each write operation is atomic - the entire SAM record is written as a single
/// unit to prevent interleaved output from different threads.
///
/// Create using `AlignmentWriterBuilder` to ensure headers are written first.
pub struct AlignmentWriter<W: Write> {
    inner: Mutex<WriterInner<W>>,
}

struct WriterInner<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> AlignmentWriter<W> {
    /// Create a builder for constructing an AlignmentWriter with headers.
    pub fn builder(writer: W) -> AlignmentWriterBuilder<W> {
        AlignmentWriterBuilder::new(writer)
    }

    pub fn flush(&self) -> std::io::Result<()> {
        let mut inner = self.inner.lock().unwrap();
        inner.writer.flush()
    }

    /// Write an alignment record atomically.
    ///
    /// The entire record is formatted into a buffer first, then written
    /// as a single operation to prevent interleaving with other threads.
    pub fn write_alignment(
        &self,
        qname: &str,
        flag: u16,
        rname: &str,
        pos: usize,
        mapq: u8,
        cigar: &str,
        rnext: &str,
        pnext: usize,
        tlen: isize,
        seq: &str,
        qual: &str,
        tags: &[(String, String)],
    ) -> std::io::Result<()> {
        // Pre-format the entire line to ensure atomic write
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            qname, flag, rname, pos + 1, mapq, cigar, rnext, pnext + 1, tlen, seq, qual
        );
        for (tag, value) in tags {
            line.push('\t');
            line.push_str(tag);
            line.push(':');
            line.push_str(value);
        }
        line.push('\n');

        // Write atomically while holding the lock
        let mut inner = self.inner.lock().unwrap();
        inner.writer.write_all(line.as_bytes())
    }
}

