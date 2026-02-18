//! Extract syncmers from a reference genome and produce a BED file.
//!
//! This extracts open syncmers (K=20, S=15) using FnvHasher — the same
//! parameters used by the parallax indexer — and outputs a BED file with:
//!   - chrom, start, end (0-based, half-open)
//!   - syncmer sequence as the name
//!   - 20 * log2(frequency) as the score (rounded to integer)
//!
//! Usage:
//!     cargo run --release --example syncmer_bed -- \
//!         --reference genome.fa \
//!         --output syncmers.bed \
//!         --primary-only

use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use parallax::kmers::Kmer;
use parallax::reference::ChromInfo;
use parallax::utils::hasher::FnvHasher;

const K: usize = 21;
const S: usize = 13;

#[derive(Parser, Debug)]
#[command(name = "syncmer-bed")]
#[command(about = "Extract syncmers from a reference genome and produce a BED file")]
struct Args {
    /// Input reference FASTA file
    #[arg(short, long)]
    reference: PathBuf,

    /// Output BED file (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Only use primary chromosomes (exclude ALT, random, decoy, unplaced contigs)
    #[arg(short = 'p', long)]
    primary_only: bool,
}

/// Infer niffler compression format from file extension.
fn format_from_ext(path: &PathBuf) -> niffler::Format {
    match path.extension().and_then(|e| e.to_str()) {
        Some("gz") => niffler::Format::Gzip,
        Some("bz2") => niffler::Format::Bzip,
        Some("xz") => niffler::Format::Lzma,
        Some("zst") => niffler::Format::Zstd,
        _ => niffler::Format::No,
    }
}

/// Load FASTA records, returning (chrom_name, sequence) pairs.
/// Supports compressed input (gzip, bzip2, xz, zstd) via niffler.
fn load_chromosomes(
    path: &PathBuf,
    primary_only: bool,
) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    use noodles::fasta;

    let (reader, format) = niffler::from_path(path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    if format != niffler::Format::No {
        log::info!("Detected {:?} compression on input", format);
    }
    let mut reader = fasta::io::Reader::new(BufReader::new(reader));
    let mut chroms = Vec::new();

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let description = record
            .description()
            .map(|d| String::from_utf8_lossy(d).into_owned());

        if primary_only {
            let info =
                ChromInfo::from_header(&name, description.as_deref().unwrap_or(""));
            if !info.is_primary() {
                continue;
            }
        }

        let seq: Vec<u8> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|&b| b.to_ascii_uppercase())
            .collect();

        chroms.push((name, seq));
    }

    Ok(chroms)
}

fn main() -> std::io::Result<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    // ── Pass 1: count k-mer frequencies ────────────────────────────────
    log::info!("Pass 1: counting syncmer frequencies");
    let chroms = load_chromosomes(&args.reference, args.primary_only)?;

    let mut kmer_counts: HashMap<u64, u32> = HashMap::new();
    for (name, seq) in &chroms {
        let mut n = 0u64;
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            seq.as_slice(),
            [(); S],
            |_pos, kmer| {
                *kmer_counts.entry(kmer.0).or_insert(0) += 1;
                n += 1;
            },
        );
        log::info!("  {} - {} syncmers", name, n);
    }

    let total: u64 = kmer_counts.values().map(|&c| c as u64).sum();
    let unique = kmer_counts.len();
    log::info!(
        "Pass 1 complete: {} total syncmers, {} unique k-mers",
        total,
        unique
    );

    // ── Pass 2: emit BED ───────────────────────────────────────────────
    log::info!("Pass 2: writing BED file");

    let writer: Box<dyn Write> = match &args.output {
        Some(path) => {
            let fmt = format_from_ext(path);
            if fmt != niffler::Format::No {
                log::info!("Writing {:?} compressed output", fmt);
            }
            niffler::to_path(path, fmt, niffler::Level::Six)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?
        }
        None => Box::new(std::io::stdout().lock()),
    };
    let mut out = BufWriter::new(writer);

    for (name, seq) in &chroms {
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
            seq.as_slice(),
            [(); S],
            |pos, kmer| {
                let freq = kmer_counts.get(&kmer.0).copied().unwrap_or(1);
                let kmer_str = kmer.to_string();
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    name,
                    pos,
                    pos + K,
                    kmer_str,
                    freq,
                )
                .expect("write failed");
            },
        );
    }

    out.flush()?;
    log::info!("Done");

    Ok(())
}
