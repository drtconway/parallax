//! Validate SAM/BAM/CRAM alignment files.
//!
//! For each mapped alignment record, checks:
//! - The reference sequence name exists in the supplied reference FASTA.
//! - The alignment start and end coordinates are within the chromosome bounds.
//! - The CIGAR read-length equals the sequence length stored in the record.
//! - The CIGAR reference span does not extend beyond the chromosome end.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Args;
use noodles::sam::alignment::io::Read as AlignmentRead;

use crate::error::Result;
use crate::reference::InMemoryReference;

/// Arguments for the `validate` subcommand.
#[derive(Args, Debug, Clone)]
pub struct ValidateArgs {
    /// Path to the input SAM file (BAM/CRAM support coming later)
    pub input: PathBuf,

    /// Path to the reference FASTA used for coordinate bounds checking
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Emit a line for every valid record, not just errors
    #[arg(short = 'v', long)]
    pub verbose: bool,

    /// Maximum allowed mismatch rate over 'M' bases across the whole file (0.0–1.0)
    #[arg(long, default_value = "0.01")]
    pub max_mismatch_rate: f64,
}

/// Counts returned by `validate_record` for accumulation across the file.
struct RecordStats {
    problems: usize,
    m_matches: u64,
    m_mismatches: u64,
}

/// A validation problem found in a single record.
#[derive(Debug)]
struct Problem {
    record_name: String,
    kind: ProblemKind,
}

#[derive(Debug)]
enum ProblemKind {
    UnknownReference(String),
    StartOutOfRange { start: usize, chrom_len: u64 },
    EndOutOfRange { end: usize, chrom_len: u64 },
    CigarReadLengthMismatch { cigar_len: usize, seq_len: usize },
    /// A `=` op where the bases differ.
    SequenceMatchMismatch { ref_pos: usize, read_pos: usize, ref_base: u8, read_base: u8 },
    SequenceMatchMismatches { first_ref_pos: usize, count: usize },
    /// A `X` op where the bases are the same.
    SequenceMismatchMatch { ref_pos: usize, read_pos: usize, base: u8 },
    SequenceMismatchMatches { first_ref_pos: usize, count: usize },
}

impl std::fmt::Display for ProblemKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProblemKind::UnknownReference(name) => {
                write!(f, "reference sequence '{}' not found in FASTA", name)
            }
            ProblemKind::StartOutOfRange { start, chrom_len } => {
                write!(f, "alignment start {} exceeds chromosome length {}", start, chrom_len)
            }
            ProblemKind::EndOutOfRange { end, chrom_len } => {
                write!(f, "alignment end {} exceeds chromosome length {}", end, chrom_len)
            }
            ProblemKind::CigarReadLengthMismatch { cigar_len, seq_len } => {
                write!(f, "CIGAR read length {} != sequence length {}", cigar_len, seq_len)
            }
            ProblemKind::SequenceMatchMismatch { ref_pos, read_pos, ref_base, read_base } => {
                write!(
                    f,
                    "'=' op mismatch at ref pos {}, read pos {}: ref={} read={}",
                    ref_pos, read_pos,
                    *ref_base as char, *read_base as char,
                )
            }
            ProblemKind::SequenceMatchMismatches { first_ref_pos, count } => {
                write!(
                    f,
                    "... and {} more '=' mismatches (first at ref pos {})",
                    count, first_ref_pos,
                )
            }
            ProblemKind::SequenceMismatchMatch { ref_pos, read_pos, base } => {
                write!(
                    f,
                    "'X' op at ref pos {}, read pos {} but both bases are {}",
                    ref_pos, read_pos, *base as char,
                )
            }
            ProblemKind::SequenceMismatchMatches { first_ref_pos, count } => {
                write!(
                    f,
                    "... and {} more 'X' false mismatches (first at ref pos {})",
                    count, first_ref_pos,
                )
            }
        }
    }
}

pub fn run(args: ValidateArgs) -> Result<()> {
    let reference = InMemoryReference::load(&args.reference, false)?;

    // Build name -> (index, length) map from the reference.
    let ref_map: HashMap<String, (usize, u64)> = reference
        .chromosomes()
        .enumerate()
        .map(|(i, (name, len))| (name.to_string(), (i, len)))
        .collect();

    let input_path = &args.input;
    let ext = input_path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    let mut problems: u64 = 0;
    let mut records_checked: u64 = 0;
    let mut total_m_matches: u64 = 0;
    let mut total_m_mismatches: u64 = 0;

    match ext.as_str() {
        "sam" => {
            let file = std::fs::File::open(input_path)?;
            let mut reader = noodles::sam::io::Reader::new(BufReader::new(file));
            let header = reader.read_alignment_header()?;
            for result in reader.alignment_records(&header) {
                let record = result?;
                let stats = validate_record(
                    record.as_ref(),
                    &header,
                    &ref_map,
                    &reference,
                    args.verbose,
                )?;
                problems += stats.problems as u64;
                total_m_matches += stats.m_matches;
                total_m_mismatches += stats.m_mismatches;
                records_checked += 1;
            }
        }
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "unsupported input format '{}'; only SAM is supported so far",
                    ext
                ),
            )
            .into());
        }
    }

    let total_m_bases = total_m_matches + total_m_mismatches;
    let mismatch_rate = if total_m_bases > 0 {
        total_m_mismatches as f64 / total_m_bases as f64
    } else {
        0.0
    };

    log::info!(
        "Checked {} records: {} problem(s) found; M-op mismatch rate {:.4}% ({}/{} bases)",
        records_checked,
        problems,
        mismatch_rate * 100.0,
        total_m_mismatches,
        total_m_bases,
    );

    if mismatch_rate > args.max_mismatch_rate {
        eprintln!(
            "ERROR: M-op mismatch rate {:.4}% exceeds maximum {:.4}%",
            mismatch_rate * 100.0,
            args.max_mismatch_rate * 100.0,
        );
        problems += 1;
    }

    if problems > 0 {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("{} validation error(s) found", problems),
        )
        .into())
    } else {
        Ok(())
    }
}

/// Validate a single alignment record.
fn validate_record(
    record: &dyn noodles::sam::alignment::Record,
    header: &noodles::sam::Header,
    ref_map: &HashMap<String, (usize, u64)>,
    reference: &InMemoryReference,
    verbose: bool,
) -> std::io::Result<RecordStats> {
    let flags = record.flags()?;
    if flags.is_unmapped() {
        return Ok(RecordStats { problems: 0, m_matches: 0, m_mismatches: 0 });
    }

    let name = record
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|| "<unnamed>".to_string());

    let mut count = 0;
    let mut m_matches = 0u64;
    let mut m_mismatches = 0u64;

    // Resolve reference sequence name from the SAM header.
    let (ref_name, chrom_len) = match record.reference_sequence(header) {
        Some(Ok((ref_name, map))) => {
            let ref_name_str = ref_name.to_string();
            // Use length from our FASTA reference map, falling back to the SAM header length.
            let chrom_len = ref_map
                .get(&ref_name_str)
                .map(|&(_, len)| len)
                .or_else(|| Some(map.length().get() as u64));
            (ref_name_str, chrom_len)
        }
        Some(Err(e)) => return Err(e),
        None => {
            // No reference sequence; record is effectively unmapped.
            return Ok(RecordStats { problems: 0, m_matches: 0, m_mismatches: 0 });
        }
    };

    // Check the reference name is known in our FASTA.
    if !ref_map.contains_key(&ref_name) {
        emit_problem(
            &Problem {
                record_name: name.clone(),
                kind: ProblemKind::UnknownReference(ref_name.clone()),
            },
        );
        count += 1;
        // Can't do coordinate checks without length.
        return Ok(RecordStats { problems: count, m_matches: 0, m_mismatches: 0 });
    }

    let chrom_len = chrom_len.unwrap_or(u64::MAX);

    // Alignment start (1-based, inclusive → convert to 0-based).
    let start_1based = match record.alignment_start() {
        Some(Ok(pos)) => usize::from(pos),
        Some(Err(e)) => return Err(e),
        None => return Ok(RecordStats { problems: 0, m_matches: 0, m_mismatches: 0 }),
    };
    let start_0based = start_1based - 1;

    if start_0based as u64 >= chrom_len {
        emit_problem(&Problem {
            record_name: name.clone(),
            kind: ProblemKind::StartOutOfRange {
                start: start_0based,
                chrom_len,
            },
        });
        count += 1;
    }

    // CIGAR validation: read-length consistency, reference span, and '='/'X' base matching.
    let cigar = record.cigar();
    let cigar_read_len = cigar.read_length()?;
    let seq_len = record.sequence().len();
    // seq_len check must happen before we consume the sequence iterator below.

    // Sequence length of 0 with a non-empty CIGAR is allowed (sequence elided).
    if seq_len > 0 && cigar_read_len != seq_len {
        emit_problem(&Problem {
            record_name: name.clone(),
            kind: ProblemKind::CigarReadLengthMismatch {
                cigar_len: cigar_read_len,
                seq_len,
            },
        });
        count += 1;
    }

    // Walk the CIGAR once: check reference span and validate '='/'X' ops.
    let chrom_idx = ref_map[&ref_name].0;
    let ref_seq = reference.sequence(chrom_idx);
    // Collect read sequence upfront so both sides are plain slices.
    let read_seq: Vec<u8> = record.sequence().iter().map(|b| b.to_ascii_uppercase()).collect();

    // Report the first MAX_REPORTED per category individually, then summarise the rest.
    const MAX_REPORTED: usize = 10;
    let mut eq_mismatch_count = 0usize;
    let mut eq_first_suppressed: Option<usize> = None;
    let mut x_match_count = 0usize;
    let mut x_first_suppressed: Option<usize> = None;

    let mut read_pos = 0usize;
    let mut ref_pos = start_0based;

    for result in cigar.iter() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        let op = result?;
        let op_len = op.len();
        match op.kind() {
            Kind::Match => {
                let ref_chunk = &ref_seq[ref_pos..ref_pos + op_len];
                let read_chunk = &read_seq[read_pos..read_pos + op_len];
                for (&ref_base, &read_base) in ref_chunk.iter().zip(read_chunk) {
                    if ref_base == read_base {
                        m_matches += 1;
                    } else {
                        m_mismatches += 1;
                    }
                }
                read_pos += op_len;
                ref_pos += op_len;
            }
            Kind::SequenceMismatch => {
                let ref_chunk = &ref_seq[ref_pos..ref_pos + op_len];
                let read_chunk = &read_seq[read_pos..read_pos + op_len];
                for (i, (&ref_base, &read_base)) in ref_chunk.iter().zip(read_chunk).enumerate() {
                    if ref_base == read_base {
                        x_match_count += 1;
                        if x_match_count <= MAX_REPORTED {
                            emit_problem(&Problem {
                                record_name: name.clone(),
                                kind: ProblemKind::SequenceMismatchMatch {
                                    ref_pos: ref_pos + i,
                                    read_pos: read_pos + i,
                                    base: ref_base,
                                },
                            });
                            count += 1;
                        } else if x_first_suppressed.is_none() {
                            x_first_suppressed = Some(ref_pos + i);
                        }
                    }
                }
                read_pos += op_len;
                ref_pos += op_len;
            }
            Kind::SequenceMatch => {
                let ref_chunk = &ref_seq[ref_pos..ref_pos + op_len];
                let read_chunk = &read_seq[read_pos..read_pos + op_len];
                for (i, (&ref_base, &read_base)) in ref_chunk.iter().zip(read_chunk).enumerate() {
                    if ref_base != read_base {
                        eq_mismatch_count += 1;
                        if eq_mismatch_count <= MAX_REPORTED {
                            emit_problem(&Problem {
                                record_name: name.clone(),
                                kind: ProblemKind::SequenceMatchMismatch {
                                    ref_pos: ref_pos + i,
                                    read_pos: read_pos + i,
                                    ref_base,
                                    read_base,
                                },
                            });
                            count += 1;
                        } else if eq_first_suppressed.is_none() {
                            eq_first_suppressed = Some(ref_pos + i);
                        }
                    }
                }
                read_pos += op_len;
                ref_pos += op_len;
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += op_len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += op_len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    if let Some(first_rp) = eq_first_suppressed {
        let suppressed = eq_mismatch_count - MAX_REPORTED;
        emit_problem(&Problem {
            record_name: name.clone(),
            kind: ProblemKind::SequenceMatchMismatches {
                first_ref_pos: first_rp,
                count: suppressed,
            },
        });
        count += 1;
    }

    if let Some(first_rp) = x_first_suppressed {
        let suppressed = x_match_count - MAX_REPORTED;
        emit_problem(&Problem {
            record_name: name.clone(),
            kind: ProblemKind::SequenceMismatchMatches {
                first_ref_pos: first_rp,
                count: suppressed,
            },
        });
        count += 1;
    }

    let end_0based = ref_pos; // ref_pos after the walk is start + ref_span
    if end_0based as u64 > chrom_len {
        emit_problem(&Problem {
            record_name: name.clone(),
            kind: ProblemKind::EndOutOfRange {
                end: end_0based,
                chrom_len,
            },
        });
        count += 1;
    }

    if verbose && count == 0 {
        log::debug!("OK: {}", name);
    }

    Ok(RecordStats { problems: count, m_matches, m_mismatches })
}

fn emit_problem(p: &Problem) {
    eprintln!("ERROR [{}]: {}", p.record_name, p.kind);
}
