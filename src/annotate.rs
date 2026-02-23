//! Annotate structural variant VCF records with library sequence identity.
//!
//! Given a library FASTA (mobile elements, viral genomes, etc.), a genome reference,
//! and a VCF of structural variants, screens variant sequences against the library
//! using syncmer lookup and confirms hits with DP alignment.

use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use crate::align::{Aligner, Kind, Op};
use crate::error::ParallaxError;
use crate::index::{self, Index, IndexBuilder};
use crate::kmers::Kmer;
use crate::reference::InMemoryReference;
use crate::utils::hasher::FnvHasher;
use crate::utils::sequence::reverse_complement_into;

/// Configuration for the annotate subcommand.
pub struct AnnotateConfig {
    pub library_fasta: PathBuf,
    pub reference_fasta: Option<PathBuf>,
    pub vcf_path: PathBuf,
    pub output: Option<PathBuf>,
    pub info_field: Option<String>,
    pub min_score: f64,
    pub threads: usize,
}

/// Best library match for a variant.
struct LibraryMatch {
    name: String,
    qual: f64,
    identity: f64,
    cigar: String,
    is_reverse: bool,
}

// ── helpers ──────────────────────────────────────────────────────────────

/// Extract a value from a VCF INFO field by key.
/// INFO is semicolon-separated "KEY=VALUE" pairs.
fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| {
        let (k, v) = kv.split_once('=')?;
        if k == key { Some(v) } else { None }
    })
}

/// Compute sequence identity from an extended CIGAR.
/// Identity = (sequence-matching columns) / (total alignment columns).
fn alignment_identity(cigar: &[Op]) -> f64 {
    let mut matches: u64 = 0;
    let mut columns: u64 = 0;
    for &op in cigar {
        let n = op.len() as u64;
        match op.kind() {
            Kind::SequenceMatch => {
                matches += n;
                columns += n;
            }
            Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion | Kind::Match => {
                columns += n;
            }
            // Soft/hard clips don't count as alignment columns
            _ => {}
        }
    }
    if columns == 0 {
        0.0
    } else {
        matches as f64 / columns as f64
    }
}

/// Convert identity (0..1) to PHRED quality, capped at 60.
fn identity_to_phred(identity: f64) -> f64 {
    if identity >= 1.0 {
        60.0
    } else if identity <= 0.0 {
        0.0
    } else {
        let q = -20.0 * (1.0 - identity).log10();
        q.min(60.0)
    }
}

// ── variant sequence extraction ──────────────────────────────────────────

/// Extract query sequence(s) from a VCF record for screening.
///
/// Returns zero or more byte sequences to screen against the library.
fn extract_variant_sequences(
    chrom: &str,
    pos: usize, // 1-based VCF POS
    ref_allele: &str,
    alt_allele: &str,
    info: &str,
    info_field: Option<&str>,
    reference: Option<&InMemoryReference>,
    ref_chrom_map: &HashMap<String, usize>,
) -> Vec<Vec<u8>> {
    let mut seqs = Vec::new();

    // If the user specified an INFO field (e.g. INSSEQ), try that first.
    if let Some(field_name) = info_field {
        if let Some(val) = info_value(info, field_name) {
            let s: Vec<u8> = val.bytes().filter(|b| b.is_ascii_alphabetic()).collect();
            if s.len() >= 20 {
                seqs.push(s);
                return seqs;
            }
        }
    }

    let svtype = info_value(info, "SVTYPE").unwrap_or("");

    match svtype {
        "INS" => {
            // Inserted sequence: ALT bases minus the anchor base
            if !alt_allele.starts_with('<') && alt_allele.len() > ref_allele.len() {
                let ins_seq: Vec<u8> = alt_allele[ref_allele.len()..].as_bytes().to_vec();
                if ins_seq.len() >= 20 {
                    seqs.push(ins_seq);
                }
            }
            // Also check INSSEQ in INFO
            if let Some(val) = info_value(info, "INSSEQ") {
                let s: Vec<u8> = val.bytes().filter(|b| b.is_ascii_alphabetic()).collect();
                if s.len() >= 20 {
                    seqs.push(s);
                }
            }
        }
        "DEL" => {
            // Deleted reference sequence
            if let Some(ref_genome) = reference {
                if let Some(&chrom_idx) = ref_chrom_map.get(chrom) {
                    let start = pos; // 0-based start (VCF POS is 1-based, after anchor)
                    let end = if let Some(end_str) = info_value(info, "END") {
                        end_str.parse::<usize>().unwrap_or(start)
                    } else if let Some(svlen_str) = info_value(info, "SVLEN") {
                        let svlen: i64 = svlen_str.parse().unwrap_or(0);
                        (pos as i64 + svlen.abs()) as usize
                    } else {
                        start
                    };
                    if end > start {
                        let chrom_len = ref_genome.chrom_length(chrom_idx) as usize;
                        let clamped_end = end.min(chrom_len);
                        if clamped_end > start {
                            let del_seq = ref_genome.get_seq(chrom_idx, start, clamped_end);
                            if del_seq.len() >= 20 {
                                seqs.push(del_seq.to_vec());
                            }
                        }
                    }
                }
            }
        }
        "DUP" => {
            // Duplicated reference region
            if let Some(ref_genome) = reference {
                if let Some(&chrom_idx) = ref_chrom_map.get(chrom) {
                    let start = pos.saturating_sub(1); // 0-based
                    let end = if let Some(end_str) = info_value(info, "END") {
                        end_str.parse::<usize>().unwrap_or(start)
                    } else {
                        start
                    };
                    if end > start {
                        let chrom_len = ref_genome.chrom_length(chrom_idx) as usize;
                        let clamped_end = end.min(chrom_len);
                        if clamped_end > start {
                            let dup_seq = ref_genome.get_seq(chrom_idx, start, clamped_end);
                            if dup_seq.len() >= 20 {
                                seqs.push(dup_seq.to_vec());
                            }
                        }
                    }
                }
            }
        }
        _ => {
            // For non-symbolic alleles, try the ALT if it's long enough
            if !alt_allele.starts_with('<') && alt_allele.len() > ref_allele.len() {
                let ins_seq: Vec<u8> = alt_allele[ref_allele.len()..].as_bytes().to_vec();
                if ins_seq.len() >= 20 {
                    seqs.push(ins_seq);
                }
            }
        }
    }

    seqs
}

// ── syncmer screening ────────────────────────────────────────────────────

/// Minimum fraction of syncmers that must hit a library sequence to be
/// considered a candidate (avoids spurious alignments).
const MIN_HIT_FRACTION: f64 = 0.15;

/// Maximum number of candidate library sequences to keep per strand.
const MAX_CANDIDATES_PER_STRAND: usize = 5;

/// Find the top library sequence candidates for one strand of the query.
///
/// Generates forward-only syncmers, tallies hits per library sequence,
/// normalises by the smaller of the query and library syncmer counts
/// (to allow partial matches while penalising spurious hits against
/// large library sequences), then filters by minimum hit fraction.
///
/// Returns the top candidates as `(chrom_idx, normalised_score)`
/// where score is `hits / min(query_syncmers, lib_syncmers)`.
fn find_strand_candidates(
    strand_seq: &[u8],
    library_index: &Index<20, 15>,
    library: &InMemoryReference,
    label: &str,
    strand_name: &str,
) -> Vec<(usize, f64)> {
    let mut hit_counts: HashMap<usize, u32> = HashMap::new();
    let mut total_syncmers: u32 = 0;

    Kmer::<20>::kmerize_open_syncmers_fwd::<15, FnvHasher, _, _>(
        strand_seq,
        [(); 15],
        |_pos, kmer| {
            total_syncmers += 1;
            library_index.with(&kmer, |_count, loci| {
                for &loc in loci {
                    let (chrom_idx, _pos) = index::decode_locus(loc);
                    *hit_counts.entry(chrom_idx).or_insert(0) += 1;
                }
            });
        },
    );

    if total_syncmers == 0 {
        return Vec::new();
    }

    // Compute normalised scores: hits / min(query_syncmers, lib_syncmers).
    // Using the minimum allows partial matches (e.g. a short query matching
    // part of a long library sequence) while suppressing spurious low-count
    // hits against large library sequences.
    let mut scored: Vec<(usize, f64, u32)> = hit_counts
        .into_iter()
        .map(|(chrom_idx, hits)| {
            let lib_syncmers = library_index.chrom_info(chrom_idx).syncmer_count;
            let denominator = if lib_syncmers > 0 {
                (total_syncmers as u64).min(lib_syncmers) as f64
            } else {
                // Fallback for legacy indexes without syncmer_count
                total_syncmers as f64
            };
            let score = hits as f64 / denominator;
            (chrom_idx, score, hits)
        })
        .collect();

    if log::log_enabled!(log::Level::Debug) {
        for &(chrom_idx, score, hits) in &scored {
            log::debug!(
                "For {} ({}), candidate {} ({} hits, score {:.3})",
                label,
                strand_name,
                library.chrom_name(chrom_idx),
                hits,
                score
            );
        }
    }

    // Filter and rank by normalised score
    scored.retain(|&(_, score, _)| score >= MIN_HIT_FRACTION);
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    scored.truncate(MAX_CANDIDATES_PER_STRAND);

    log::debug!(
        "For {} ({}), {} query syncmers, {} candidates (min score {:.2})",
        label,
        strand_name,
        total_syncmers,
        scored.len(),
        MIN_HIT_FRACTION
    );

    scored.into_iter().map(|(idx, score, _)| (idx, score)).collect()
}

/// Screen a query sequence against the library index on both strands,
/// then confirm top candidates with DP alignment.
///
/// Each strand independently selects its top candidates, which are then
/// combined for alignment. This ensures that a strong match on one strand
/// is not masked by noise on the other.
fn screen_and_align(
    query_id: &str,
    query: &[u8],
    library_index: &Index<20, 15>,
    library: &InMemoryReference,
    aligner: &mut Aligner,
    rc_buf: &mut Vec<u8>,
) -> Option<LibraryMatch> {
    if query.len() < 20 {
        return None;
    }

    // Find top candidates per strand
    let fwd_candidates = find_strand_candidates(query, library_index, library, query_id, "+");

    reverse_complement_into(query, rc_buf);
    let rev_candidates = find_strand_candidates(rc_buf, library_index, library, query_id, "-");

    // Combine: tag each with its strand, then align
    let mut best: Option<LibraryMatch> = None;

    for &(chrom_idx, _score) in &fwd_candidates {
        if let Some(m) = try_align(query_id, query, chrom_idx, false, library, aligner) {
            let dominated = best.as_ref().is_some_and(|b| b.qual >= m.qual);
            if !dominated {
                best = Some(m);
            }
        }
    }

    for &(chrom_idx, _score) in &rev_candidates {
        if let Some(m) = try_align(query_id, rc_buf, chrom_idx, true, library, aligner) {
            let dominated = best.as_ref().is_some_and(|b| b.qual >= m.qual);
            if !dominated {
                best = Some(m);
            }
        }
    }

    best
}

/// Align a strand-specific query against a single library sequence candidate.
fn try_align(
    query_id: &str,
    strand_seq: &[u8],
    chrom_idx: usize,
    is_reverse: bool,
    library: &InMemoryReference,
    aligner: &mut Aligner,
) -> Option<LibraryMatch> {
    let lib_seq = library.sequence(chrom_idx);
    let lib_name = library.chrom_name(chrom_idx).to_string();

    let numerator = lib_seq.len().min(strand_seq.len());
    let denominator = lib_seq.len().max(strand_seq.len());

    if denominator == 0 || (numerator as f64 / denominator as f64) < 0.5 {
        log::debug!(
            "Skipping alignment of {} ({} strand) to {} due to length mismatch ({:.2})",
            query_id,
            if is_reverse { "reverse" } else { "forward" },
            lib_name,
            numerator as f64 / denominator as f64
        );
        return None;
    }

    log::debug!(
        "Aligning {} ({} strand, {}bp) to library {} ({}bp)",
        query_id,
        if is_reverse { "reverse" } else { "forward" },
        strand_seq.len(),
        lib_name,
        lib_seq.len()
    );

    let aln = aligner.align(strand_seq, lib_seq)?;
    let identity = alignment_identity(&aln.cigar);
    let qual = identity_to_phred(identity);
    let strand_ch = if is_reverse { '-' } else { '+' };

    log::debug!(
        "Aligned {} ({}) to library {}: identity {:.1}%, qual {:.1}",
        query_id,
        strand_ch,
        lib_name,
        100.0 * identity,
        qual
    );

    Some(LibraryMatch {
        name: lib_name,
        qual,
        identity,
        cigar: aln.cigar_string(),
        is_reverse,
    })
}

// ── main entry point ─────────────────────────────────────────────────────

/// Run the annotate subcommand.
pub fn run(config: AnnotateConfig) -> Result<(), ParallaxError> {
    // 1. Load library sequences and build syncmer index
    log::info!(
        "Loading library sequences from {}",
        config.library_fasta.display()
    );
    let library = InMemoryReference::load(&config.library_fasta, false)?;

    log::info!(
        "Building library index ({} sequences)",
        library.num_chroms()
    );
    let library_index: Index<20, 15> = IndexBuilder::build_parallel(&library, None, config.threads);

    // 2. Optionally load genome reference (needed for DEL/DUP sequence extraction)
    let reference = if let Some(ref ref_path) = config.reference_fasta {
        log::info!("Loading genome reference from {}", ref_path.display());
        Some(InMemoryReference::load(ref_path, false)?)
    } else {
        None
    };

    // Build name→index map for reference contigs
    let ref_chrom_map: HashMap<String, usize> = reference
        .as_ref()
        .map(|r| {
            (0..r.num_chroms())
                .map(|i| (r.chrom_name(i).to_string(), i))
                .collect()
        })
        .unwrap_or_default();

    // 3. Open VCF input (handles plain text, gzip, bgzf via niffler)
    let vcf_reader: Box<dyn BufRead> = if config.vcf_path.to_str() == Some("-") {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let (reader, _compression) = niffler::from_path(&config.vcf_path)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        Box::new(BufReader::new(reader))
    };

    // 4. Open output
    let mut out: Box<dyn Write> = if let Some(ref path) = config.output {
        Box::new(BufWriter::new(std::fs::File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    // 5. Create aligner (uses defaults since we don't need the full config)
    let mut aligner = Aligner::with_defaults();
    let mut rc_buf: Vec<u8> = Vec::new();

    let mut n_records: u64 = 0;
    let mut n_annotated: u64 = 0;

    for line_result in vcf_reader.lines() {
        let line = line_result?;

        if line.starts_with('#') {
            // Inject INFO definitions before the #CHROM line
            if line.starts_with("#CHROM") {
                writeln!(
                    out,
                    "##INFO=<ID=LIBRARY,Number=1,Type=String,\
                     Description=\"Best matching library sequence name\">"
                )?;
                writeln!(
                    out,
                    "##INFO=<ID=LIBRARY_STRAND,Number=1,Type=String,\
                     Description=\"Strand of library match (+ or -)\">"
                )?;
                writeln!(
                    out,
                    "##INFO=<ID=LIBRARY_QUAL,Number=1,Type=Float,\
                     Description=\"Library match quality (PHRED-scaled identity)\">"
                )?;
                writeln!(
                    out,
                    "##INFO=<ID=LIBRARY_CIGAR,Number=1,Type=String,\
                     Description=\"CIGAR of alignment to library sequence\">"
                )?;
            }
            writeln!(out, "{}", line)?;
            continue;
        }

        n_records += 1;

        // Parse VCF data line (tab-separated)
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            writeln!(out, "{}", line)?;
            continue;
        }

        let chrom = fields[0];
        let pos: usize = fields[1].parse().unwrap_or(0);
        let id = fields[2];
        let ref_allele = fields[3];
        let alt_allele = fields[4];
        let info = fields[7];

        // Extract variant sequence(s) to screen
        let query_seqs = extract_variant_sequences(
            chrom,
            pos,
            ref_allele,
            alt_allele,
            info,
            config.info_field.as_deref(),
            reference.as_ref(),
            &ref_chrom_map,
        );

        // Screen each query against library index and align best candidate
        let best_match = query_seqs
            .iter()
            .flat_map(|seq| {
                screen_and_align(id, seq, &library_index, &library, &mut aligner, &mut rc_buf)
            })
            .max_by(|a, b| {
                a.qual
                    .partial_cmp(&b.qual)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some(ref m) = best_match {

            log::info!(
                "Best library match for {}: {} (strand {}, qual {:.1}) ({:.1}% identity)",
                id,
                m.name,
                if m.is_reverse { '-' } else { '+' },
                m.qual,
                m.identity * 100.0
            );

            if m.qual >= config.min_score {
                // Write annotated record: append to INFO column
                let strand_ch = if m.is_reverse { '-' } else { '+' };
                let annotation = format!(
                    "LIBRARY={};LIBRARY_STRAND={};LIBRARY_QUAL={:.1};LIBRARY_CIGAR={}",
                    m.name, strand_ch, m.qual, m.cigar
                );
                let new_info = if info == "." {
                    annotation
                } else {
                    format!("{};{}", info, annotation)
                };

                for (i, f) in fields.iter().enumerate() {
                    if i > 0 {
                        write!(out, "\t")?;
                    }
                    if i == 7 {
                        write!(out, "{}", new_info)?;
                    } else {
                        write!(out, "{}", f)?;
                    }
                }
                writeln!(out)?;

                n_annotated += 1;
                continue;
            }
        }

        // No qualifying match — write record unchanged
        writeln!(out, "{}", line)?;
    }

    log::info!(
        "Processed {} records, annotated {} ({:.1}%)",
        n_records,
        n_annotated,
        if n_records > 0 {
            100.0 * n_annotated as f64 / n_records as f64
        } else {
            0.0
        }
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_info_value() {
        let info = "SVTYPE=INS;END=12345;INSSEQ=ACGTACGT";
        assert_eq!(info_value(info, "SVTYPE"), Some("INS"));
        assert_eq!(info_value(info, "END"), Some("12345"));
        assert_eq!(info_value(info, "INSSEQ"), Some("ACGTACGT"));
        assert_eq!(info_value(info, "MISSING"), None);
    }

    #[test]
    fn test_info_value_dot() {
        assert_eq!(info_value(".", "SVTYPE"), None);
    }

    #[test]
    fn test_alignment_identity_perfect() {
        let cigar = vec![Op::new(Kind::SequenceMatch, 100)];
        assert!((alignment_identity(&cigar) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_alignment_identity_with_mismatches() {
        let cigar = vec![
            Op::new(Kind::SequenceMatch, 90),
            Op::new(Kind::SequenceMismatch, 10),
        ];
        assert!((alignment_identity(&cigar) - 0.9).abs() < 1e-9);
    }

    #[test]
    fn test_alignment_identity_with_indels() {
        let cigar = vec![
            Op::new(Kind::SequenceMatch, 90),
            Op::new(Kind::Insertion, 5),
            Op::new(Kind::Deletion, 5),
        ];
        assert!((alignment_identity(&cigar) - 0.9).abs() < 1e-9);
    }

    #[test]
    fn test_alignment_identity_empty() {
        let cigar: Vec<Op> = vec![];
        assert!((alignment_identity(&cigar) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_identity_to_phred() {
        // Perfect identity → 60 (cap)
        assert!((identity_to_phred(1.0) - 60.0).abs() < 1e-9);

        // Zero identity → 0
        assert!((identity_to_phred(0.0) - 0.0).abs() < 1e-9);

        // 90% identity → Q10
        assert!((identity_to_phred(0.9) - 10.0).abs() < 0.01);

        // 99% identity → Q20
        assert!((identity_to_phred(0.99) - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_extract_insertion_from_alt() {
        let ref_chrom_map = HashMap::new();
        let seqs = extract_variant_sequences(
            "chr1",
            100,
            "A",
            "AACGTACGTACGTACGTACGTACGT", // 24-base insertion
            "SVTYPE=INS",
            None,
            None,
            &ref_chrom_map,
        );
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0], b"ACGTACGTACGTACGTACGTACGT");
    }

    #[test]
    fn test_extract_short_insertion_skipped() {
        let ref_chrom_map = HashMap::new();
        let seqs = extract_variant_sequences(
            "chr1",
            100,
            "A",
            "AACGT", // Only 4 bases — too short
            "SVTYPE=INS",
            None,
            None,
            &ref_chrom_map,
        );
        assert!(seqs.is_empty());
    }

    #[test]
    fn test_extract_from_info_field() {
        let ref_chrom_map = HashMap::new();
        let long_seq = "A".repeat(30);
        let info = format!("SVTYPE=INS;MYSEQ={}", long_seq);
        let seqs = extract_variant_sequences(
            "chr1",
            100,
            "A",
            "<INS>",
            &info,
            Some("MYSEQ"),
            None,
            &ref_chrom_map,
        );
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].len(), 30);
    }
}
