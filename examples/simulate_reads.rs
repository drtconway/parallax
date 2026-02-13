//! Simulate long reads from a reference genome.
//!
//! Samples reads with:
//! - Random positions from the reference
//! - Random strand (forward or reverse complement)
//! - Read lengths from a normal distribution
//! - Optional substitution errors at a configurable rate
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
//! Use --primary-only to exclude ALT scaffolds, random, decoy, and unplaced contigs.
//!
//! Read names encode the true origin:
//!   sim_00000001:chr1_1000_2000_+
//! With --include-errors-in-name, errors are appended:
//!   sim_00000001:chr1_1000_2000_+:A50T,G100C
//! Where A50T means position 50 (1-based in read) was changed from A to T.

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
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

    // Calculate total genome size for weighted sampling
    let total_size: usize = contigs.iter().map(|c| c.sequence.len()).sum();
    let contig_weights: Vec<f64> = contigs
        .iter()
        .map(|c| c.sequence.len() as f64 / total_size as f64)
        .collect();

    log::info!(
        "Loaded {} contigs, total size: {} bp",
        contigs.len(),
        total_size
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
        // Sample contig (weighted by length)
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

        // Sample read length
        let read_len = sample_length(&mut rng, &mut length_dist, args.min_length, args.max_length);

        // Ensure we can fit the read
        if read_len >= contig.sequence.len() {
            continue;
        }

        // Sample start position
        let max_start = contig.sequence.len() - read_len;
        let start: usize = rng.random_range(0..=max_start);

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

        // Generate read name
        let strand_char = if is_reverse { '-' } else { '+' };
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
        let read_name = format!(
            "sim_{:08}:{}_{}_{}_{}{}",
            i,
            contig.name,
            start + 1,        // 1-based position
            start + read_len, // inclusive end position
            strand_char,
            errors_str
        );

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
