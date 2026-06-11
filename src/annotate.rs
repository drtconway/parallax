//! Annotate structural variant VCF records with library sequence identity.
//!
//! Given a library FASTA (mobile elements, viral genomes, etc.), a genome reference,
//! and a VCF of structural variants, screens variant sequences against the library
//! using syncmer lookup and confirms hits with DP alignment.

use std::collections::HashMap;
use std::io::{self, BufReader, BufWriter};
use std::path::PathBuf;

use noodles::vcf;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use noodles::vcf::header::record::value::{Map, map::Info as InfoMap};
use noodles::vcf::variant::io::Write as VcfWrite;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use parallax::index::{Index, PackedLocus};

use crate::align::{DpAligner, Kind, Op};
use parallax::error::ParallaxError;
use parallax::reference::InMemoryReference;
use parallax::utils::hasher::FnvHasher;
use parallax::utils::sequence::reverse_complement_into;
use parallax::{
    index::fwd_index::{FwdIndex, FwdIndexBuilder},
    kmers::Kmer,
};

/// Configuration for the annotate subcommand.
pub struct AnnotateConfig {
    pub library_fasta: PathBuf,
    pub index_path: Option<PathBuf>,
    pub reference_fasta: Option<PathBuf>,
    pub vcf_path: PathBuf,
    pub output: Option<PathBuf>,
    pub info_field: Option<String>,
    pub min_score: f64,
    pub emit_cigar: bool,
    pub portable: bool,
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

/// Extract a string value from a RecordBuf INFO field.
fn get_info_str<'a>(info: &'a vcf::variant::record_buf::Info, key: &str) -> Option<&'a str> {
    match info.get(key)? {
        Some(InfoValue::String(s)) => Some(s.as_str()),
        _ => None,
    }
}

/// Extract an integer value from a RecordBuf INFO field.
fn get_info_int(info: &vcf::variant::record_buf::Info, key: &str) -> Option<i32> {
    match info.get(key)? {
        Some(InfoValue::Integer(i)) => Some(*i),
        _ => None,
    }
}

/// Format the record IDs for logging (semicolon-separated, or ".").
fn record_id_string(record: &vcf::variant::RecordBuf) -> String {
    let ids = record.ids().as_ref();
    if ids.is_empty() {
        ".".to_string()
    } else {
        ids.iter().map(String::as_str).collect::<Vec<_>>().join(";")
    }
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

/// Extract query sequence(s) from a VCF RecordBuf for screening.
///
/// Returns zero or more byte sequences to screen against the library.
fn extract_variant_sequences(
    record: &vcf::variant::RecordBuf,
    info_field: Option<&str>,
    reference: Option<&InMemoryReference>,
    ref_chrom_map: &HashMap<String, usize>,
) -> Vec<Vec<u8>> {
    let chrom = record.reference_sequence_name();
    let pos = record.variant_start().map(|p| usize::from(p)).unwrap_or(0);
    let ref_allele = record.reference_bases();
    let alt_alleles: &[String] = record.alternate_bases().as_ref();
    let alt_allele = alt_alleles.first().map(|s| s.as_str()).unwrap_or("");
    let info = record.info();

    let mut seqs = Vec::new();

    // If the user specified an INFO field (e.g. INSSEQ), try that first.
    if let Some(field_name) = info_field {
        if let Some(val) = get_info_str(info, field_name) {
            let s: Vec<u8> = val.bytes().filter(|b| b.is_ascii_alphabetic()).collect();
            if s.len() >= 20 {
                seqs.push(s);
                return seqs;
            }
        }
    }

    let svtype = get_info_str(info, "SVTYPE").unwrap_or("");

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
            if let Some(val) = get_info_str(info, "INSSEQ") {
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
                    let end = if let Some(end_val) = get_info_int(info, "END") {
                        end_val as usize
                    } else if let Some(svlen_val) = get_info_int(info, "SVLEN") {
                        (pos as i64 + (svlen_val as i64).abs()) as usize
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
                    let end = if let Some(end_val) = get_info_int(info, "END") {
                        end_val as usize
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
    library_index: &FwdIndex<20, 15>,
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
                for loc in loci {
                    let (chrom_idx, _pos, _strand) = loc.unpack();
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

    scored
        .into_iter()
        .map(|(idx, score, _)| (idx, score))
        .collect()
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
    library_index: &FwdIndex<20, 15>,
    library: &InMemoryReference,
    aligner: &mut DpAligner,
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
    aligner: &mut DpAligner,
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
    // 1. Load library sequences and either load or build syncmer index
    log::info!(
        "Loading library sequences from {}",
        config.library_fasta.display()
    );
    let library = InMemoryReference::load(&config.library_fasta, false)?;

    let library_index: FwdIndex<20, 15> = if let Some(ref index_path) = config.index_path {
        if index_path.join("chrom_info.json").exists() {
            log::info!("Loading library index from {}", index_path.display());
            if config.portable {
                FwdIndex::load(index_path)?
            } else {
                FwdIndex::load_feather(index_path)?
            }
        } else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!(
                    "Index not found at {}. Use 'parallax index' to build it first.",
                    index_path.display()
                ),
            )
            .into());
        }
    } else {
        log::info!(
            "Building library index ({} sequences)",
            library.num_chroms()
        );
        FwdIndexBuilder::build_parallel(&library, None, config.threads)
    };

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
    let vcf_input: Box<dyn io::BufRead> = if config.vcf_path.to_str() == Some("-") {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let (reader, _compression) = niffler::from_path(&config.vcf_path)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        Box::new(BufReader::new(reader))
    };

    let mut reader = vcf::io::Reader::new(vcf_input);
    let mut header = reader.read_header()?;

    // Add INFO field definitions to the header
    header.infos_mut().insert(
        String::from("LIBRARY"),
        Map::<InfoMap>::new(
            Number::Count(1),
            Type::String,
            "Best matching library sequence name",
        ),
    );
    header.infos_mut().insert(
        String::from("LIBRARY_STRAND"),
        Map::<InfoMap>::new(
            Number::Count(1),
            Type::String,
            "Strand of library match (+ or -)",
        ),
    );
    header.infos_mut().insert(
        String::from("LIBRARY_IDENTITY"),
        Map::<InfoMap>::new(
            Number::Count(1),
            Type::Integer,
            "Library percent match identity",
        ),
    );
    header.infos_mut().insert(
        String::from("LIBRARY_QUAL"),
        Map::<InfoMap>::new(
            Number::Count(1),
            Type::Float,
            "Library match quality (PHRED-scaled identity)",
        ),
    );
    if config.emit_cigar {
        header.infos_mut().insert(
            String::from("LIBRARY_CIGAR"),
            Map::<InfoMap>::new(
                Number::Count(1),
                Type::String,
                "CIGAR of alignment to library sequence",
            ),
        );
    }

    // 4. Open output and write header
    let vcf_output: Box<dyn io::Write> = if let Some(ref path) = config.output {
        Box::new(BufWriter::new(std::fs::File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut writer = vcf::io::Writer::new(vcf_output);
    writer.write_header(&header)?;

    // 5. Process records
    let mut aligner = DpAligner::with_defaults();
    let mut rc_buf: Vec<u8> = Vec::new();

    let mut n_records: u64 = 0;
    let mut n_annotated: u64 = 0;

    let mut record = vcf::variant::RecordBuf::default();
    while reader.read_record_buf(&header, &mut record)? != 0 {
        n_records += 1;

        let id = record_id_string(&record);

        // Extract variant sequence(s) to screen
        let query_seqs = extract_variant_sequences(
            &record,
            config.info_field.as_deref(),
            reference.as_ref(),
            &ref_chrom_map,
        );

        // Screen each query against library index and align best candidate
        let best_match = query_seqs
            .iter()
            .flat_map(|seq| {
                screen_and_align(
                    &id,
                    seq,
                    &library_index,
                    &library,
                    &mut aligner,
                    &mut rc_buf,
                )
            })
            .max_by(|a, b| {
                a.qual
                    .partial_cmp(&b.qual)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some(ref m) = best_match {
            if m.qual >= config.min_score {
                log::info!(
                    "Match for {}: {} (strand {}, qual {:.1}) ({:.1}% identity)",
                    id,
                    m.name,
                    if m.is_reverse { '-' } else { '+' },
                    m.qual,
                    m.identity * 100.0
                );

                let strand_str = if m.is_reverse { "-" } else { "+" };
                record.info_mut().insert(
                    String::from("LIBRARY"),
                    Some(InfoValue::String(m.name.clone())),
                );
                record.info_mut().insert(
                    String::from("LIBRARY_STRAND"),
                    Some(InfoValue::String(strand_str.to_string())),
                );
                record.info_mut().insert(
                    String::from("LIBRARY_IDENTITY"),
                    Some(InfoValue::Integer((m.identity * 100.0) as i32)),
                );
                record.info_mut().insert(
                    String::from("LIBRARY_QUAL"),
                    Some(InfoValue::Float(m.qual as f32)),
                );
                if config.emit_cigar {
                    record.info_mut().insert(
                        String::from("LIBRARY_CIGAR"),
                        Some(InfoValue::String(m.cigar.clone())),
                    );
                }
                n_annotated += 1;
            }
        }

        writer.write_variant_record(&header, &record)?;
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
    fn test_get_info_str() {
        let info: vcf::variant::record_buf::Info = [
            (
                String::from("SVTYPE"),
                Some(InfoValue::String("INS".into())),
            ),
            (String::from("END"), Some(InfoValue::Integer(12345))),
            (
                String::from("INSSEQ"),
                Some(InfoValue::String("ACGTACGT".into())),
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(get_info_str(&info, "SVTYPE"), Some("INS"));
        assert_eq!(get_info_str(&info, "INSSEQ"), Some("ACGTACGT"));
        assert_eq!(get_info_str(&info, "MISSING"), None);
        // Integer field should return None for get_info_str
        assert_eq!(get_info_str(&info, "END"), None);
    }

    #[test]
    fn test_get_info_int() {
        let info: vcf::variant::record_buf::Info = [
            (String::from("END"), Some(InfoValue::Integer(12345))),
            (
                String::from("SVTYPE"),
                Some(InfoValue::String("INS".into())),
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(get_info_int(&info, "END"), Some(12345));
        assert_eq!(get_info_int(&info, "MISSING"), None);
        // String field should return None for get_info_int
        assert_eq!(get_info_int(&info, "SVTYPE"), None);
    }

    #[test]
    fn test_get_info_empty() {
        let info: vcf::variant::record_buf::Info = Default::default();
        assert_eq!(get_info_str(&info, "SVTYPE"), None);
        assert_eq!(get_info_int(&info, "END"), None);
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

        // 90% identity → Q20  (-20 * log10(0.1))
        assert!((identity_to_phred(0.9) - 20.0).abs() < 0.01);

        // 99% identity → Q40  (-20 * log10(0.01))
        assert!((identity_to_phred(0.99) - 40.0).abs() < 0.01);
    }

    #[test]
    fn test_extract_insertion_from_alt() {
        use noodles::core::Position;
        use vcf::variant::record_buf::AlternateBases;

        let ref_chrom_map = HashMap::new();
        let info: vcf::variant::record_buf::Info = [(
            String::from("SVTYPE"),
            Some(InfoValue::String("INS".into())),
        )]
        .into_iter()
        .collect();
        let record = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::new(100).expect("valid position"))
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec![
                "AACGTACGTACGTACGTACGTACGT".to_string(),
            ]))
            .set_info(info)
            .build();
        let seqs = extract_variant_sequences(&record, None, None, &ref_chrom_map);
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0], b"ACGTACGTACGTACGTACGTACGT");
    }

    #[test]
    fn test_extract_short_insertion_skipped() {
        use noodles::core::Position;
        use vcf::variant::record_buf::AlternateBases;

        let ref_chrom_map = HashMap::new();
        let info: vcf::variant::record_buf::Info = [(
            String::from("SVTYPE"),
            Some(InfoValue::String("INS".into())),
        )]
        .into_iter()
        .collect();
        let record = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::new(100).expect("valid position"))
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec!["AACGT".to_string()]))
            .set_info(info)
            .build();
        let seqs = extract_variant_sequences(&record, None, None, &ref_chrom_map);
        assert!(seqs.is_empty());
    }

    #[test]
    fn test_extract_from_info_field() {
        use noodles::core::Position;
        use vcf::variant::record_buf::AlternateBases;

        let ref_chrom_map = HashMap::new();
        let long_seq = "A".repeat(30);
        let info: vcf::variant::record_buf::Info = [
            (
                String::from("SVTYPE"),
                Some(InfoValue::String("INS".into())),
            ),
            (String::from("MYSEQ"), Some(InfoValue::String(long_seq))),
        ]
        .into_iter()
        .collect();
        let record = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::new(100).expect("valid position"))
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec!["<INS>".to_string()]))
            .set_info(info)
            .build();
        let seqs = extract_variant_sequences(&record, Some("MYSEQ"), None, &ref_chrom_map);
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].len(), 30);
    }
}
