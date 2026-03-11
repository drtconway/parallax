//! Compute per-read alignment coverage from a SAM or BAM file.
//!
//! For each read, uses the primary alignment and any supplementary alignments
//! listed in the SA tag to determine which portions of the read are covered.
//! Secondary alignments (FLAG 0x100) are ignored entirely.
//!
//! Output (TSV to stdout):
//!   read_name  type  read_start  read_end  length
//!
//! One row is emitted per contiguous segment of the read:
//! - `mapped`:    bases covered by ≥1 primary/supplementary alignment
//! - `left`:      unmapped bases at the 5' end (before the first mapped region)
//! - `interior`:  unmapped bases flanked by mapped regions on both sides
//! - `right`:     unmapped bases at the 3' end (after the last mapped region)
//!
//! Coordinates are 0-based half-open [read_start, read_end) in forward-read space.
//!
//! Usage:
//!     cargo run --example read_coverage -- --input alignments.bam

use std::{
    fs::File,
    io::{self, BufReader, Write},
};

use clap::Parser;
use noodles::{
    bam::Record as BamRecord,
    sam::alignment::{
        record::{Cigar, Flags, cigar::op::Kind, data::field::Tag},
        record_buf::{RecordBuf, data::field::Value},
    },
};

#[derive(Parser, Debug)]
#[command(name = "read-coverage")]
#[command(about = "Compute per-read coverage from primary and supplementary alignments")]
struct Args {
    /// Input SAM or BAM file ('-' for stdin as SAM)
    #[arg(short, long)]
    input: String,
}

fn parse_cigar(s: &str) -> Vec<(Kind, usize)> {
    let mut ops = Vec::new();
    let mut len: usize = 0;
    for ch in s.chars() {
        if let Some(d) = ch.to_digit(10) {
            len = len * 10 + d as usize;
        } else {
            let op = match ch {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                _ => {
                    len = 0;
                    continue;
                }
            };
            if len > 0 {
                ops.push((op, len));
            }
            len = 0;
        }
    }
    ops
}

// ── Interval computation ─────────────────────────────────────────────────────

/// Given a CIGAR and strand ('+'/'-'), return the half-open interval
/// `[qstart, qend)` of the read *in forward-read* coordinates that this
/// alignment covers.
///
/// "Covered" means the read bases are not soft/hard-clipped — they are
/// aligned to (or inserted relative to) the reference.
///
/// * `read_len` is the full read length (including any hard-clipped bases)
///   needed only for minus-strand flipping.
fn cigar_read_interval(ops: &[(Kind, usize)], is_reverse: bool, read_len: usize) -> (usize, usize) {
    let mut qpos = 0;
    let mut qstart = None;
    let mut qend = None;

    for &(op, len) in ops {
        let consumes_query = matches!(
            op,
            Kind::Match
                | Kind::Insertion
                | Kind::SoftClip
                | Kind::HardClip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        );
        if consumes_query {
            if !matches!(op, Kind::SoftClip | Kind::HardClip) {
                // This region is covered (not clipped)
                if qstart.is_none() {
                    qstart = Some(qpos);
                }
                qend = Some(qpos + len);
            }
            qpos += len;
        }
    }

    let (s, e) = match (qstart, qend) {
        (Some(s), Some(e)) => (s, e),
        _ => return (0, 0),
    };

    if is_reverse {
        // Flip to forward-read coordinates: [read_len - e, read_len - s)
        (read_len - e, read_len - s)
    } else {
        (s, e)
    }
}

// ── SA tag parsing ────────────────────────────────────────────────────────────

/// Parse the SA tag value and return a list of read intervals covered by
/// supplementary alignments, in forward-read coordinates.
///
/// SA tag format: `(rname,pos,strand,CIGAR,mapq,NM;)*`
fn parse_sa_tag(sa: &str, read_len: usize) -> Vec<(usize, usize)> {
    let mut intervals = Vec::new();
    for entry in sa.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let fields: Vec<&str> = entry.split(',').collect();
        if fields.len() < 4 {
            continue;
        }
        let strand = fields[2];
        let cigar_str = fields[3];
        let is_reverse = strand == "-";
        let ops = parse_cigar(cigar_str);
        if ops.is_empty() {
            continue;
        }
        let interval = cigar_read_interval(&ops, is_reverse, read_len);
        if interval.0 < interval.1 {
            intervals.push(interval);
        }
    }
    intervals
}

// ── Interval utilities ────────────────────────────────────────────────────────

/// Merge a list of (start, end) intervals into a sorted, non-overlapping list.
fn merge_intervals(mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    intervals.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    merged
}

/// Compute coverage statistics from merged intervals.
///
/// Returns `(mapped, terminal_unmapped, interior_unmapped)`.
#[cfg_attr(not(test), allow(dead_code))]
fn coverage_stats(intervals: &[(usize, usize)], read_len: usize) -> (usize, usize, usize) {
    if intervals.is_empty() {
        // Fully unmapped: all bases are "terminal"
        return (0, read_len, 0);
    }

    let mapped: usize = intervals.iter().map(|(s, e)| e - s).sum();
    let first_start = intervals[0].0;
    let last_end = intervals.last().unwrap().1;

    // Bases before the first mapped region + bases after the last mapped region
    let terminal = first_start + (read_len - last_end);

    // Total unmapped = read_len - mapped; subtract terminal to get interior gaps
    let interior = read_len - mapped - terminal;

    (mapped, terminal, interior)
}

// ── Record processing ─────────────────────────────────────────────────────────

static SA_TAG: std::sync::OnceLock<Tag> = std::sync::OnceLock::new();

fn sa_tag() -> Tag {
    *SA_TAG.get_or_init(|| Tag::try_from(*b"SA").expect("invalid tag"))
}

/// Extract coverage intervals from a primary alignment record.
///
/// Returns `None` if the record is secondary (should be skipped) or if the
/// record is supplementary (should be skipped — we get supp info from SA tags
/// on the primary).
///
/// For an unmapped read, returns an empty interval list with the read length.
fn process_record(record: &RecordBuf) -> Option<(String, usize, Vec<(usize, usize)>)> {
    let flags = record.flags();

    // Skip secondary alignments outright
    if flags.contains(Flags::SECONDARY) {
        return None;
    }
    // Skip supplementary alignments — their coverage is captured via SA tags
    // on the primary record
    if flags.contains(Flags::SUPPLEMENTARY) {
        return None;
    }

    let read_name = record
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|| "*".to_string());

    // Determine read length from the sequence field; fall back to CIGAR if
    // the sequence is '*' (not stored, e.g. in some supplementary records).
    let seq_len = record.sequence().len();
    let cigar_ops = record
        .cigar()
        .iter()
        .map(|op| op.unwrap())
        .map(|op| (op.kind(), op.len() as usize))
        .collect::<Vec<_>>();

    // Compute hard-clip total from CIGAR to get the full read length.
    let hard_clips: usize = cigar_ops
        .iter()
        .filter(|(op, _)| *op == Kind::HardClip)
        .map(|(_, l)| l)
        .sum();
    let read_len = seq_len + hard_clips;

    if read_len == 0 {
        return Some((read_name, 0, vec![]));
    }

    let mut intervals: Vec<(usize, usize)> = Vec::new();

    // Primary alignment interval
    if !flags.contains(Flags::UNMAPPED) && !cigar_ops.is_empty() {
        let is_reverse = flags.contains(Flags::REVERSE_COMPLEMENTED);
        let (qstart, qend) = cigar_read_interval(&cigar_ops, is_reverse, read_len);
        if qstart < qend {
            intervals.push((qstart, qend));
        }
    }

    // Supplementary alignment intervals from SA tag
    if let Some(value) = record.data().get(&sa_tag()) {
        let sa_str = match value {
            Value::String(s) => String::from_utf8_lossy(s).to_string(),
            _ => String::new(),
        };
        if !sa_str.is_empty() {
            intervals.extend(parse_sa_tag(&sa_str, read_len));
        }
    }

    Some((read_name, read_len, intervals))
}

/// Identical logic to `process_record` but operates on a lazy `bam::Record`.
///
/// The key difference: `bam::Record::data().get()` scans auxiliary fields
/// linearly and returns the *first* match without ever validating for
/// duplicate tags.  This means a BAM record that has a repeated tag (e.g.
/// two `s1` fields, which is technically invalid) will not cause an error;
/// we simply extract the first SA tag value we encounter.
fn process_bam_record(record: &BamRecord) -> Option<(String, usize, Vec<(usize, usize)>)> {
    let flags = record.flags();

    if flags.contains(Flags::SECONDARY) {
        return None;
    }
    if flags.contains(Flags::SUPPLEMENTARY) {
        return None;
    }

    let read_name = record
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|| "*".to_string());

    let seq_len = record.sequence().len();
    let cigar_ops = record
        .cigar()
        .iter()
        .map(|op| op.unwrap())
        .map(|op| (op.kind(), op.len() as usize))
        .collect::<Vec<_>>();

    let hard_clips: usize = cigar_ops
        .iter()
        .filter(|(op, _)| *op == Kind::HardClip)
        .map(|(_, l)| l)
        .sum();
    let read_len = seq_len + hard_clips;

    if read_len == 0 {
        return Some((read_name, 0, vec![]));
    }

    let mut intervals: Vec<(usize, usize)> = Vec::new();

    if !flags.contains(Flags::UNMAPPED) && !cigar_ops.is_empty() {
        let is_reverse = flags.contains(Flags::REVERSE_COMPLEMENTED);
        let (qstart, qend) = cigar_read_interval(&cigar_ops, is_reverse, read_len);
        if qstart < qend {
            intervals.push((qstart, qend));
        }
    }

    if let Some(Ok(value)) = record.data().get(&sa_tag()) {
        use noodles::sam::alignment::record::data::field::Value as V;
        let sa_str = match value {
            V::String(s) => String::from_utf8_lossy(s).to_string(),
            _ => String::new(),
        };
        if !sa_str.is_empty() {
            intervals.extend(parse_sa_tag(&sa_str, read_len));
        }
    }

    Some((read_name, read_len, intervals))
}

// ── I/O helpers ───────────────────────────────────────────────────────────────

fn detect_format(path: &str) -> &'static str {
    let lower = path.to_lowercase();
    if lower.ends_with(".bam") {
        "bam"
    } else if lower.ends_with(".cram") {
        "cram"
    } else {
        "sam"
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> io::Result<()> {
    let args = Args::parse();

    let out = io::stdout();
    let mut out = io::BufWriter::new(out.lock());

    // Write TSV header
    writeln!(out, "read_name\ttype\tread_start\tread_end\tlength")?;

    let format = detect_format(&args.input);

    match format {
        "bam" => {
            let file = File::open(&args.input)
                .map_err(|e| io::Error::other(format!("Cannot open {}: {e}", args.input)))?;
            let mut reader = noodles::bam::io::Reader::new(file);
            // read_header() advances past the header bytes; we don't need the
            // result for lazy record reading.
            reader.read_header()?;
            for result in reader.records() {
                let record = result?;
                if let Some((read_name, read_len, intervals)) = process_bam_record(&record) {
                    emit_segments(&read_name, read_len, intervals, &mut out)?;
                }
            }
        }
        _ => {
            // SAM (default) — also handles stdin if path is "-"
            let rdr: Box<dyn io::Read> =
                if args.input == "-" {
                    Box::new(io::stdin())
                } else {
                    Box::new(File::open(&args.input).map_err(|e| {
                        io::Error::other(format!("Cannot open {}: {e}", args.input))
                    })?)
                };
            let buf = BufReader::new(rdr);
            let mut reader = noodles::sam::io::Reader::new(buf);
            let header = reader.read_header()?;
            let mut record = RecordBuf::default();
            loop {
                let bytes_read = reader.read_record_buf(&header, &mut record)?;
                if bytes_read == 0 {
                    break;
                }
                emit_record(&record, &mut out)?;
            }
        }
    }

    out.flush()?;
    Ok(())
}

fn emit_segments(
    read_name: &str,
    read_len: usize,
    intervals: Vec<(usize, usize)>,
    out: &mut impl Write,
) -> io::Result<()> {
    let merged = merge_intervals(intervals);

    if merged.is_empty() {
        // Fully unmapped — the entire read is a "left" segment
        writeln!(out, "{read_name}\tleft\t0\t{read_len}\t{read_len}")?;
        return Ok(());
    }

    let mut pos = 0usize;
    for (i, &(start, end)) in merged.iter().enumerate() {
        // Unmapped gap before this mapped segment
        if start > pos {
            let label = if i == 0 { "left" } else { "interior" };
            let len = start - pos;
            if len > 0 {
                writeln!(out, "{read_name}\t{label}\t{pos}\t{start}\t{len}")?;
            }
        }
        // Mapped segment
        let len = end - start;
        writeln!(out, "{read_name}\tmapped\t{start}\t{end}\t{len}")?;
        pos = end;
    }
    // Unmapped tail
    if pos < read_len {
        let len = read_len - pos;
        writeln!(out, "{read_name}\tright\t{pos}\t{read_len}\t{len}")?;
    }

    Ok(())
}

fn emit_record(record: &RecordBuf, out: &mut impl Write) -> io::Result<()> {
    if let Some((read_name, read_len, intervals)) = process_record(record) {
        emit_segments(&read_name, read_len, intervals, out)?;
    }
    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar_simple() {
        let ops = parse_cigar("10S80M10S");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (Kind::SoftClip, 10));
        assert_eq!(ops[1], (Kind::Match, 80));
        assert_eq!(ops[2], (Kind::SoftClip, 10));
    }

    #[test]
    fn test_parse_cigar_hard_clips() {
        let ops = parse_cigar("30H10S50M10H");
        assert_eq!(ops[0], (Kind::HardClip, 30));
        assert_eq!(ops[1], (Kind::SoftClip, 10));
        assert_eq!(ops[2], (Kind::Match, 50));
        assert_eq!(ops[3], (Kind::HardClip, 10));
    }

    #[test]
    fn test_interval_forward_soft_clipped() {
        // 10S 80M 10S on a 100bp read: mapped region is [10, 90)
        let ops = parse_cigar("10S80M10S");
        let (s, e) = cigar_read_interval(&ops, false, 100);
        assert_eq!((s, e), (10, 90));
    }

    #[test]
    fn test_interval_reverse_soft_clipped() {
        // Minus strand 10S 80M 10S on a 100bp read.
        // In RC coords: [10, 90). Flipped: [100-90, 100-10) = [10, 90).
        let ops = parse_cigar("10S80M10S");
        let (s, e) = cigar_read_interval(&ops, true, 100);
        assert_eq!((s, e), (10, 90));
    }

    #[test]
    fn test_interval_hard_clipped_supplementary_forward() {
        // A supplementary on a 100bp read; 30H 10S 50M 10H
        // total_query = 30+10+50+10 = 100; leading=40, trailing=10
        // qstart=40, qend=90; forward: [40, 90)
        let ops = parse_cigar("30H10S50M10H");
        let (s, e) = cigar_read_interval(&ops, false, 100);
        assert_eq!((s, e), (40, 90));
    }

    #[test]
    fn test_interval_hard_clipped_supplementary_reverse() {
        // Same CIGAR on minus strand: flip [40, 90) → [100-90, 100-40) = [10, 60)
        let ops = parse_cigar("30H10S50M10H");
        let (s, e) = cigar_read_interval(&ops, true, 100);
        assert_eq!((s, e), (10, 60));
    }

    #[test]
    fn test_merge_intervals() {
        let iv = vec![(0, 30), (20, 50), (60, 80)];
        let merged = merge_intervals(iv);
        assert_eq!(merged, vec![(0, 50), (60, 80)]);
    }

    #[test]
    fn test_coverage_stats_fully_mapped() {
        let iv = vec![(0, 100)];
        let (mapped, terminal, interior) = coverage_stats(&iv, 100);
        assert_eq!((mapped, terminal, interior), (100, 0, 0));
    }

    #[test]
    fn test_coverage_stats_terminal_gaps() {
        // Mapped region covers [10, 90) — 10 terminal on each end, 0 interior
        let iv = vec![(10, 90)];
        let (mapped, terminal, interior) = coverage_stats(&iv, 100);
        assert_eq!((mapped, terminal, interior), (80, 20, 0));
    }

    #[test]
    fn test_coverage_stats_interior_gap() {
        // Two mapped regions [10, 40) and [60, 90) — interior gap [40, 60)
        let iv = vec![(10, 40), (60, 90)];
        let (mapped, terminal, interior) = coverage_stats(&iv, 100);
        assert_eq!((mapped, terminal, interior), (60, 20, 20));
    }

    #[test]
    fn test_coverage_stats_fully_unmapped() {
        let iv: Vec<(usize, usize)> = vec![];
        let (mapped, terminal, interior) = coverage_stats(&iv, 100);
        assert_eq!((mapped, terminal, interior), (0, 100, 0));
    }

    #[test]
    fn test_parse_sa_tag() {
        // A 200bp read; supplementary: plus strand, CIGAR 100H80M20H → [100, 180)
        let sa = "chr1,12345,+,100H80M20H,60,2;";
        let intervals = parse_sa_tag(sa, 200);
        assert_eq!(intervals, vec![(100, 180)]);
    }

    #[test]
    fn test_parse_sa_tag_reverse() {
        // 200bp read; minus strand, CIGAR 10H80M110H
        // total_query = 200; leading=10, trailing=110; qstart=10, qend=90
        // reverse → [200-90, 200-10) = [110, 190)
        let sa = "chr1,12345,-,10H80M110H,60,2;";
        let intervals = parse_sa_tag(sa, 200);
        assert_eq!(intervals, vec![(110, 190)]);
    }
}
