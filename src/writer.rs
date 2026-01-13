use std::io::{BufWriter, Write};
use std::sync::Mutex;

/// Thread-safe alignment writer that supports multiple threads writing concurrently.
///
/// Each write operation is atomic - the entire SAM record is written as a single
/// unit to prevent interleaved output from different threads.
pub struct AlignmentWriter<W: Write> {
    inner: Mutex<WriterInner<W>>,
}

struct WriterInner<W: Write> {
    writer: BufWriter<W>,
    header_written: bool,
}

impl<W: Write> AlignmentWriter<W> {
    /// Create a new AlignmentWriter that writes to the given writer.
    pub fn new(writer: W) -> Self {
        Self {
            inner: Mutex::new(WriterInner {
                writer: BufWriter::new(writer),
                header_written: false,
            }),
        }
    }

    pub fn flush(&self) -> std::io::Result<()> {
        let mut inner = self.inner.lock().unwrap();
        inner.writer.flush()
    }

    pub fn write_contig_header(&self, name: &str, length: usize) -> std::io::Result<()> {
        let mut inner = self.inner.lock().unwrap();
        if !inner.header_written {
            writeln!(inner.writer, "@HD\tVN:1.0\tSO:unsorted")?;
            inner.header_written = true;
        }
        writeln!(inner.writer, "@SQ\tSN:{}\tLN:{}", name, length)
    }

    pub fn write_command_header(&self, command: &str) -> std::io::Result<()> {
        let mut inner = self.inner.lock().unwrap();
        if !inner.header_written {
            writeln!(inner.writer, "@HD\tVN:1.0\tSO:unsorted")?;
            inner.header_written = true;
        }
        writeln!(inner.writer, "@PG\tID:parallax\tPN:parallax\tCL:{}", command)
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
        if !inner.header_written {
            writeln!(inner.writer, "@HD\tVN:1.0\tSO:unsorted")?;
            inner.header_written = true;
        }
        inner.writer.write_all(line.as_bytes())
    }
}

