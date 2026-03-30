//! Simulate long reads from a reference genome.
//!
//! Samples reads with:
//! - Random positions from the reference
//! - Random strand (forward or reverse complement)
//! - Read lengths from a normal distribution
//! - Optional substitution errors at a configurable rate
//! - Optional structural variants applied via a VCF file
//!
//! Usage:
//!   cargo run --release --example simulate_reads -- \\
//!     --reference genome.fa \\
//!     --output reads.fq \\
//!     --num-reads 1000 \\
//!     --mean-length 15000 \\
//!     --std-dev 3000 \\
//!     --error-rate 0.01 \\
//!     --primary-only
//!
//! With structural variants:
//!   cargo run --release --example simulate_reads -- \\
//!     --reference genome.fa \\
//!     --vcf variants.vcf \\
//!     --output reads.fq \\
//!     --num-reads 1000
//!
//! When --vcf is given, reads are biased to overlap variant regions unless
//! --global-sampling is also specified.
//!
//! Use --primary-only to exclude ALT scaffolds, random, decoy, and unplaced contigs.
//!
//! Read names encode the expected alignment on the original reference:
//!   sim_00000001:chr1_1000_2000_+
//! With structural variants, reads may map to multiple segments:
//!   sim_00000001:chr1_1000_1500_+,chr1_3000_3500_+
//! With --include-errors-in-name, errors are appended:
//!   sim_00000001:chr1_1000_2000_+:A50T,G100C
//! Where A50T means position 50 (1-based in read) was changed from A to T.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use parallax::utils::rope::Rope;
use rand::prelude::*;

#[derive(Parser, Debug)]
#[command(name = "simulate_reads")]
#[command(about = "Simulate long reads from a reference genome")]
struct Args {
    /// Input reference FASTA file
    #[arg(short, long)]
    reference: PathBuf,

    /// Output FASTQ file
    #[arg(short, long)]
    output: PathBuf,

    /// Number of reads to generate
    #[arg(short, long, default_value = "1000")]
    num_reads: usize,

    /// Mean read length
    #[arg(short, long, default_value = "15000")]
    mean_length: f64,

    /// Standard deviation of read length
    #[arg(short = 'd', long = "std-dev", default_value = "3000.0")]
    std_dev: f64,

    /// Minimum read length
    #[arg(long, default_value = "500")]
    min_length: usize,

    /// Maximum read length (0 = no limit)
    #[arg(long, default_value = "0")]
    max_length: usize,

    /// Random seed (optional)
    #[arg(short, long)]
    seed: Option<u64>,

    /// Base quality score (Phred)
    #[arg(short, long, default_value = "20")]
    quality: u8,

    /// Substitution error rate (0.0 to 1.0)
    #[arg(short, long, default_value = "0.0")]
    error_rate: f64,

    /// Include introduced errors in read name
    #[arg(long, default_value = "false")]
    include_errors_in_name: bool,

    /// Only use primary chromosomes (exclude ALT, random, decoy, unplaced contigs)
    #[arg(short = 'p', long)]
    primary_only: bool,

    /// VCF file containing structural variants to apply to the reference.
    /// Supported SVTYPE values: DEL, DUP, INV. INS is not supported.
    /// Overlapping variants on the same contig are forbidden.
    #[arg(long)]
    vcf: Option<PathBuf>,

    /// When --vcf is given, sample reads uniformly across the whole genome
    /// instead of biasing towards variant regions.
    #[arg(long)]
    global_sampling: bool,
}

/// Complement a single nucleotide
#[inline]
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Reverse complement a sequence
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// A substitution error: position (1-based), original base, new base
struct Substitution {
    pos: usize,   // 1-based position in read
    from: u8,     // original base
    to: u8,       // substituted base
}

impl std::fmt::Display for Substitution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}{}", self.from as char, self.pos, self.to as char)
    }
}

/// Introduce substitution errors into a sequence
/// Returns the mutated sequence and a list of substitutions made
fn introduce_errors(seq: &[u8], error_rate: f64, rng: &mut StdRng) -> (Vec<u8>, Vec<Substitution>) {
    if error_rate <= 0.0 {
        return (seq.to_vec(), Vec::new());
    }

    let mut result = seq.to_vec();
    let mut substitutions = Vec::new();

    for (i, base) in result.iter_mut().enumerate() {
        // Skip non-ACGT bases
        if !BASES.contains(base) {
            continue;
        }

        let r: f64 = rng.random();
        if r < error_rate {
            let original = *base;
            // Pick a different base
            loop {
                let new_base = BASES[rng.random_range(0..4)];
                if new_base != original {
                    substitutions.push(Substitution {
                        pos: i + 1,  // 1-based
                        from: original,
                        to: new_base,
                    });
                    *base = new_base;
                    break;
                }
            }
        }
    }

    (result, substitutions)
}

/// A chromosome/contig with its sequence
struct Contig {
    name: String,
    sequence: Vec<u8>,
}

/// Check if a contig name represents a primary chromosome.
/// Primary chromosomes are chr1-22, chrX, chrY, chrM without any suffix.
fn is_primary_contig(name: &str) -> bool {
    // Non-primary indicators
    if name.contains("_alt") {
        return false;
    }
    if name.contains("_random") {
        return false;
    }
    if name.contains("_decoy") {
        return false;
    }
    if name.starts_with("chrUn") {
        return false;
    }
    if name.starts_with("HLA-") {
        return false;
    }
    
    // Primary: starts with "chr" and has no underscores
    if name.starts_with("chr") && !name.contains('_') {
        return true;
    }
    
    // Also accept numeric chromosomes (1, 2, ..., X, Y, MT)
    if name.parse::<u32>().is_ok() {
        return true;
    }
    if name == "X" || name == "Y" || name == "MT" {
        return true;
    }
    
    false
}

// ── Structural variant types ─────────────────────────────────────────────────

#[derive(Clone, Debug)]
enum SvType {
    Del,
    Dup,
    Inv,
}

#[derive(Clone, Debug)]
struct StructuralVariant {
    chrom: String,
    /// 0-based start of affected region.
    start: usize,
    /// 0-based exclusive end of affected region.
    end: usize,
    svtype: SvType,
}

/// A variant's location in the modified-contig coordinate space.
struct VariantRegion {
    start: usize,
    end: usize,
}

/// A block mapping modified-sequence coordinates back to original-reference
/// coordinates.  Used to translate read positions into expected alignments.
struct MappingBlock {
    /// Start in modified-sequence coords (0-based).
    modified_start: usize,
    /// Exclusive end in modified-sequence coords.
    modified_end: usize,
    /// Start in original-reference coords (0-based).
    ref_start: usize,
    /// Exclusive end in original-reference coords.
    ref_end: usize,
    /// True if this block is an inversion (reverse complement).
    inverted: bool,
}

/// A contig whose sequence may have been modified by structural variants.
struct ModifiedContig {
    name: String,
    sequence: Vec<u8>,
    variant_regions: Vec<VariantRegion>,
    mapping_blocks: Vec<MappingBlock>,
}

/// Index into the flat list of all variant regions across all contigs.
struct SamplingTarget {
    contig_idx: usize,
    region_start: usize,
    region_end: usize,
}

/// An expected alignment segment in original-reference coordinates.
struct ExpectedSegment {
    chrom: String,
    /// 1-based start position.
    start: usize,
    /// 1-based inclusive end position.
    end: usize,
    strand: char,
}

impl std::fmt::Display for ExpectedSegment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}_{}_{}_{}", self.chrom, self.start, self.end, self.strand)
    }
}

/// Map a read's modified-coordinate range back to original-reference segments.
fn compute_expected_segments(
    chrom: &str,
    read_start: usize,
    read_end: usize,
    is_reverse: bool,
    blocks: &[MappingBlock],
) -> Vec<ExpectedSegment> {
    let mut segments = Vec::new();
    for block in blocks {
        if block.modified_end <= read_start || block.modified_start >= read_end {
            continue;
        }
        let inter_start = read_start.max(block.modified_start);
        let inter_end = read_end.min(block.modified_end);
        let offset_start = inter_start - block.modified_start;
        let offset_end = inter_end - block.modified_start;

        let (ref_seg_start, ref_seg_end) = if block.inverted {
            (block.ref_end - offset_end, block.ref_end - offset_start)
        } else {
            (block.ref_start + offset_start, block.ref_start + offset_end)
        };

        let strand = if is_reverse ^ block.inverted { '-' } else { '+' };
        segments.push(ExpectedSegment {
            chrom: chrom.to_string(),
            start: ref_seg_start + 1, // 1-based
            end: ref_seg_end,         // 1-based inclusive (= 0-based exclusive)
            strand,
        });
    }
    segments
}

// ── VCF loading ──────────────────────────────────────────────────────────────

/// Extract a value from a VCF INFO field by key (e.g. "SVTYPE", "END").
fn parse_info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|field| {
        let (k, v) = field.split_once('=')?;
        if k == key { Some(v) } else { None }
    })
}

/// Load structural variants from a VCF file.
///
/// Only DEL, DUP, and INV are accepted. Returns an error on INS or if
/// required INFO fields (SVTYPE, END) are missing.
fn load_vcf(path: &PathBuf) -> std::io::Result<Vec<StructuralVariant>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut variants = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = fields[0].to_string();
        let pos: usize = fields[1]
            .parse()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Invalid POS: {e}")))?;
        let info = fields[7];

        let svtype_str = parse_info_value(info, "SVTYPE").ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Missing SVTYPE in INFO at {chrom}:{pos}"),
            )
        })?;

        let svtype = match svtype_str {
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "INS" => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("INS variants are not supported (at {chrom}:{pos})"),
                ))
            }
            other => {
                log::warn!("Skipping unsupported SVTYPE {other} at {chrom}:{pos}");
                continue;
            }
        };

        let end: usize = parse_info_value(info, "END")
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Missing END in INFO at {chrom}:{pos}"),
                )
            })?
            .parse()
            .map_err(|e| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, format!("Invalid END: {e}"))
            })?;

        // VCF POS is 1-based and is the anchor base preceding the event.
        // The affected region starts at position POS+1 (1-based) = POS (0-based).
        // VCF INFO:END is the 1-based position of the last affected base,
        // which equals the 0-based exclusive end.
        variants.push(StructuralVariant {
            chrom,
            start: pos,
            end,
            svtype,
        });
    }

    Ok(variants)
}

/// Group variants by chromosome, sort by position, and check that no two
/// variants on the same contig overlap.
fn group_and_validate_variants(
    variants: Vec<StructuralVariant>,
) -> std::io::Result<HashMap<String, Vec<StructuralVariant>>> {
    let mut by_chrom: HashMap<String, Vec<StructuralVariant>> = HashMap::new();
    for var in variants {
        by_chrom.entry(var.chrom.clone()).or_default().push(var);
    }

    for (chrom, vars) in by_chrom.iter_mut() {
        vars.sort_by_key(|v| (v.start, v.end));
        for i in 1..vars.len() {
            if vars[i].start < vars[i - 1].end {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "Overlapping variants on {chrom}: [{}, {}) and [{}, {})",
                        vars[i - 1].start,
                        vars[i - 1].end,
                        vars[i].start,
                        vars[i].end
                    ),
                ));
            }
        }
    }

    Ok(by_chrom)
}

// ── Rope-based contig modification ───────────────────────────────────────────

/// Apply structural variants to a contig using a [`Rope`] and materialise the
/// result.  Returns the modified sequence and the position of each variant in
/// the modified coordinate space.
///
/// `variants` must be sorted by `start` (left-to-right) with no overlaps.
/// Internally the function iterates left-to-right, building the rope by
/// concatenating unchanged and modified segments so that original-reference
/// coordinates remain valid throughout.
fn apply_variants(
    sequence: &[u8],
    variants: &[StructuralVariant],
) -> (Vec<u8>, Vec<VariantRegion>, Vec<MappingBlock>) {
    // Pre-compute reverse-complement buffers for INV variants (the rope
    // borrows these alongside the original sequence).
    let inv_buffers: Vec<Vec<u8>> = variants
        .iter()
        .filter(|v| matches!(v.svtype, SvType::Inv))
        .map(|v| reverse_complement(&sequence[v.start..v.end]))
        .collect();

    let mut pieces: Vec<Rope<[u8]>> = Vec::new();
    let mut regions = Vec::new();
    let mut blocks = Vec::new();
    let mut cursor = 0usize;
    let mut modified_pos = 0usize;
    let mut inv_idx = 0usize;

    for var in variants {
        // Unchanged gap before this variant.
        if var.start > cursor {
            let gap = var.start - cursor;
            pieces.push(Rope::from(&sequence[cursor..var.start]));
            blocks.push(MappingBlock {
                modified_start: modified_pos,
                modified_end: modified_pos + gap,
                ref_start: cursor,
                ref_end: var.start,
                inverted: false,
            });
            modified_pos += gap;
        }

        let region_len = var.end - var.start;

        match var.svtype {
            SvType::Del => {
                // Nothing added — the region is removed.
                regions.push(VariantRegion {
                    start: modified_pos,
                    end: modified_pos,
                });
            }
            SvType::Dup => {
                // Two copies of the region.
                pieces.push(Rope::from(&sequence[var.start..var.end]));
                pieces.push(Rope::from(&sequence[var.start..var.end]));
                regions.push(VariantRegion {
                    start: modified_pos,
                    end: modified_pos + 2 * region_len,
                });
                blocks.push(MappingBlock {
                    modified_start: modified_pos,
                    modified_end: modified_pos + region_len,
                    ref_start: var.start,
                    ref_end: var.end,
                    inverted: false,
                });
                blocks.push(MappingBlock {
                    modified_start: modified_pos + region_len,
                    modified_end: modified_pos + 2 * region_len,
                    ref_start: var.start,
                    ref_end: var.end,
                    inverted: false,
                });
                modified_pos += 2 * region_len;
            }
            SvType::Inv => {
                // Reverse-complemented region from the pre-computed buffer.
                pieces.push(Rope::from(inv_buffers[inv_idx].as_slice()));
                inv_idx += 1;
                regions.push(VariantRegion {
                    start: modified_pos,
                    end: modified_pos + region_len,
                });
                blocks.push(MappingBlock {
                    modified_start: modified_pos,
                    modified_end: modified_pos + region_len,
                    ref_start: var.start,
                    ref_end: var.end,
                    inverted: true,
                });
                modified_pos += region_len;
            }
        }

        cursor = var.end;
    }

    // Trailing unchanged sequence.
    if cursor < sequence.len() {
        let tail = sequence.len() - cursor;
        pieces.push(Rope::from(&sequence[cursor..]));
        blocks.push(MappingBlock {
            modified_start: modified_pos,
            modified_end: modified_pos + tail,
            ref_start: cursor,
            ref_end: sequence.len(),
            inverted: false,
        });
    }

    let rope = if pieces.is_empty() {
        Rope::from(sequence)
    } else {
        pieces.into_iter().reduce(|a, b| a + b).unwrap()
    };

    (Vec::from(rope), regions, blocks)
}

/// Build modified contigs by applying the variants grouped by chromosome.
fn build_modified_contigs(
    contigs: Vec<Contig>,
    variants_by_chrom: &HashMap<String, Vec<StructuralVariant>>,
) -> Vec<ModifiedContig> {
    contigs
        .into_iter()
        .map(|c| {
            if let Some(vars) = variants_by_chrom.get(&c.name) {
                let (seq, regions, blocks) = apply_variants(&c.sequence, vars);
                log::info!(
                    "Applied {} variant(s) to {} ({} bp -> {} bp)",
                    vars.len(),
                    c.name,
                    c.sequence.len(),
                    seq.len(),
                );
                ModifiedContig {
                    name: c.name,
                    sequence: seq,
                    variant_regions: regions,
                    mapping_blocks: blocks,
                }
            } else {
                let len = c.sequence.len();
                ModifiedContig {
                    name: c.name,
                    sequence: c.sequence,
                    variant_regions: Vec::new(),
                    mapping_blocks: vec![MappingBlock {
                        modified_start: 0,
                        modified_end: len,
                        ref_start: 0,
                        ref_end: len,
                        inverted: false,
                    }],
                }
            }
        })
        .collect()
}

/// Load all contigs from a FASTA file
fn load_reference(path: &PathBuf, primary_only: bool) -> std::io::Result<Vec<Contig>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reader = noodles::fasta::io::Reader::new(reader);

    let mut contigs = Vec::new();

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).to_string();
        let sequence: Vec<u8> = record.sequence().as_ref().to_vec();

        // Skip very short contigs
        if sequence.len() < 1000 {
            log::warn!("Skipping short contig: {} ({} bp)", name, sequence.len());
            continue;
        }

        // Skip non-primary contigs if primary_only is set
        if primary_only && !is_primary_contig(&name) {
            log::info!("Skipping non-primary contig: {}", name);
            continue;
        }

        log::info!("Loaded contig: {} ({} bp)", name, sequence.len());
        contigs.push(Contig { name, sequence });
    }

    Ok(contigs)
}

struct NormalDistribution {
    mean: f64,
    std_dev: f64,
    cache: Vec<f64>,
}

impl NormalDistribution {
    fn new(mean: f64, std_dev: f64) -> Self {
        Self {
            mean,
            std_dev,
            cache: Vec::new(),
        }
    }

    fn sample(&mut self, rng: &mut StdRng) -> f64 {
        if self.cache.len() == 0 {
            let u1: f64 = rng.random_range(0.0..1.0);
            let u2: f64 = rng.random_range(0.0..1.0);

            let v = (-2.0 * u1.ln()).sqrt();
            let w = 2.0 * std::f64::consts::PI * u2;
            let z0 = v * w.cos();
            let z1 = v * w.sin();

            self.cache.push(self.mean + z1 * self.std_dev);
            self.cache.push(self.mean + z0 * self.std_dev);
        }
        self.cache.pop().unwrap()
    }
}

/// Sample a read length, respecting min/max bounds
fn sample_length(
    rng: &mut StdRng,
    dist: &mut NormalDistribution,
    min_length: usize,
    max_length: usize,
) -> usize {
    let min_length = min_length as f64;
    let max_length = max_length as f64;
    loop {
        let len = dist.sample(rng).round();
        if len >= min_length && (max_length == 0.0 || len <= max_length) {
            return len as usize;
        }
    }
}

/// Write a FASTQ record
fn write_fastq_record<W: Write>(
    writer: &mut W,
    name: &str,
    sequence: &[u8],
    quality: u8,
) -> std::io::Result<()> {
    // Header
    writeln!(writer, "@{}", name)?;

    // Sequence
    writer.write_all(sequence)?;
    writeln!(writer)?;

    // Plus line
    writeln!(writer, "+")?;

    // Quality (constant quality score for all bases)
    let qual_char = (quality + 33) as char;
    for _ in 0..sequence.len() {
        write!(writer, "{}", qual_char)?;
    }
    writeln!(writer)?;

    Ok(())
}

fn main() -> std::io::Result<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    // Initialize RNG
    let mut rng: StdRng = match args.seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_rng(&mut rand::rng()),
    };

    // Load reference
    log::info!(
        "Loading reference from {:?}{}",
        args.reference,
        if args.primary_only { " (primary contigs only)" } else { "" }
    );
    let contigs = load_reference(&args.reference, args.primary_only)?;

    if contigs.is_empty() {
        return Err(std::io::Error::other("No contigs found in reference"));
    }

    log::info!(
        "Loaded {} contigs, total size: {} bp",
        contigs.len(),
        contigs.iter().map(|c| c.sequence.len()).sum::<usize>()
    );

    // Optionally load and apply structural variants.
    let has_vcf = args.vcf.is_some();
    let contigs = if let Some(vcf_path) = &args.vcf {
        log::info!("Loading structural variants from {:?}", vcf_path);
        let variants = load_vcf(vcf_path)?;
        log::info!("Loaded {} variant record(s)", variants.len());
        let by_chrom = group_and_validate_variants(variants)?;
        build_modified_contigs(contigs, &by_chrom)
    } else {
        contigs
            .into_iter()
            .map(|c| {
                let len = c.sequence.len();
                ModifiedContig {
                    name: c.name,
                    sequence: c.sequence,
                    variant_regions: Vec::new(),
                    mapping_blocks: vec![MappingBlock {
                        modified_start: 0,
                        modified_end: len,
                        ref_start: 0,
                        ref_end: len,
                        inverted: false,
                    }],
                }
            })
            .collect()
    };

    // Biased sampling: when a VCF is provided and --global-sampling is NOT set,
    // reads are drawn so they overlap at least one variant region.
    let biased = has_vcf && !args.global_sampling;

    let sampling_targets: Vec<SamplingTarget> = contigs
        .iter()
        .enumerate()
        .flat_map(|(idx, c)| {
            c.variant_regions.iter().map(move |r| SamplingTarget {
                contig_idx: idx,
                region_start: r.start,
                region_end: r.end,
            })
        })
        .collect();

    if biased && sampling_targets.is_empty() {
        return Err(std::io::Error::other(
            "--vcf was given but no applicable variants were found on any loaded contig",
        ));
    }

    // Calculate total genome size for weighted sampling (global mode).
    let total_size: usize = contigs.iter().map(|c| c.sequence.len()).sum();
    let contig_weights: Vec<f64> = contigs
        .iter()
        .map(|c| c.sequence.len() as f64 / total_size as f64)
        .collect();

    log::info!(
        "Modified genome size: {} bp, {} variant target(s), sampling: {}",
        total_size,
        sampling_targets.len(),
        if biased { "biased" } else { "global" }
    );

    // Create length distribution
    let mut length_dist = NormalDistribution::new(args.mean_length, args.std_dev);

    // Open output file
    let file = File::create(&args.output)?;
    let mut writer = BufWriter::new(file);

    log::info!("Generating {} reads...", args.num_reads);

    let mut total_bases = 0usize;
    let mut forward_count = 0usize;
    let mut reverse_count = 0usize;

    for i in 0..args.num_reads {
        // Sample read length (needed early for biased position sampling).
        let read_len = sample_length(&mut rng, &mut length_dist, args.min_length, args.max_length);

        // Choose contig and start position.
        let (contig_idx, start) = if biased {
            // Pick a random variant target and sample a position that overlaps it.
            let target = &sampling_targets[rng.random_range(0..sampling_targets.len())];
            let contig = &contigs[target.contig_idx];

            if read_len >= contig.sequence.len() {
                continue;
            }

            // Treat zero-width DEL breakpoints as width-1 so the overlap
            // window is non-empty.
            let window_end = target.region_end.max(target.region_start + 1);

            let valid_start = target.region_start.saturating_sub(read_len - 1);
            let valid_end = (window_end - 1).min(contig.sequence.len().saturating_sub(read_len));

            if valid_start > valid_end {
                continue;
            }

            (target.contig_idx, rng.random_range(valid_start..=valid_end))
        } else {
            // Global sampling: weighted by contig length.
            let contig_idx = {
                let r: f64 = rng.random();
                let mut cumulative = 0.0;
                let mut idx = 0;
                for (j, &w) in contig_weights.iter().enumerate() {
                    cumulative += w;
                    if r < cumulative {
                        idx = j;
                        break;
                    }
                }
                idx
            };

            let contig = &contigs[contig_idx];
            if read_len >= contig.sequence.len() {
                continue;
            }

            let max_start = contig.sequence.len() - read_len;
            (contig_idx, rng.random_range(0..=max_start))
        };

        let contig = &contigs[contig_idx];

        // Sample strand
        let is_reverse: bool = rng.random();

        // Extract sequence
        let seq = &contig.sequence[start..start + read_len];
        let seq = if is_reverse {
            reverse_count += 1;
            reverse_complement(seq)
        } else {
            forward_count += 1;
            seq.to_vec()
        };

        // Introduce errors
        let (seq, substitutions) = introduce_errors(&seq, args.error_rate, &mut rng);

        // Generate read name with expected original-reference segments.
        let segments = compute_expected_segments(
            &contig.name,
            start,
            start + read_len,
            is_reverse,
            &contig.mapping_blocks,
        );
        let segments_str = segments
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let errors_str = if !args.include_errors_in_name || substitutions.is_empty() {
            String::new()
        } else {
            format!(
                ":{}",
                substitutions
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        let read_name = format!("sim_{:08}:{}{}", i, segments_str, errors_str);

        // Write FASTQ record
        write_fastq_record(&mut writer, &read_name, &seq, args.quality)?;

        total_bases += read_len;

        if (i + 1) % 10000 == 0 {
            log::info!("Generated {} reads...", i + 1);
        }
    }

    writer.flush()?;

    let mean_len = total_bases as f64 / args.num_reads as f64;
    log::info!(
        "Done! Generated {} reads ({} forward, {} reverse)",
        args.num_reads,
        forward_count,
        reverse_count
    );
    log::info!(
        "Total bases: {}, mean length: {:.1}",
        total_bases,
        mean_len
    );

    Ok(())
}
