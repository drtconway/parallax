use std::{collections::HashMap, path::PathBuf};

use clap::Parser;
use parallax::kmers::Kmer;
use parallax::utils::hasher::{FnvHasher, Hasher, IdentityHasher};

#[derive(Parser, Debug)]
#[command(name = "syncmer-histogram")]
#[command(about = "Compute syncmer k-mer frequency histograms comparing IdentityHasher vs FnvHasher")]
struct Args {
    /// Input FASTA file
    input: PathBuf,
}

const K: usize = 20;
const S: usize = 15;

/// Build k-mer frequency histogram for syncmers using the given hasher.
/// Returns (histogram, total_syncmers, unique_kmers)
fn build_histogram<H: Hasher>(fasta_path: &PathBuf) -> std::io::Result<(HashMap<u64, u64>, u64, usize)> {
    let mut reader = noodles::fasta::io::reader::Builder::default()
        .build_from_path(fasta_path)?;

    let mut kmer_counts: HashMap<u64, u64> = HashMap::new();

    for record in reader.records() {
        let record = record?;
        let seq = record.sequence().as_ref();

        Kmer::<K>::kmerize_open_syncmers_fwd::<S, H, _, _>(seq, [(); S], |_pos, kmer| {
            *kmer_counts.entry(kmer.0).or_insert(0) += 1;
        });
    }

    // Build frequency histogram (count -> how many k-mers have that count)
    let mut histogram: HashMap<u64, u64> = HashMap::new();
    for &count in kmer_counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }

    let total: u64 = kmer_counts.values().sum();
    let unique = kmer_counts.len();

    Ok((histogram, total, unique))
}

fn main() -> std::io::Result<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    // First pass: IdentityHasher
    log::info!("Pass 1: Computing syncmer histogram with IdentityHasher...");
    let (identity_histogram, identity_total, identity_unique) =
        build_histogram::<IdentityHasher>(&args.input)?;
    log::info!(
        "IdentityHasher: {} total syncmers, {} unique k-mers",
        identity_total,
        identity_unique
    );

    // Second pass: FnvHasher
    log::info!("Pass 2: Computing syncmer histogram with FnvHasher...");
    let (fnv_histogram, fnv_total, fnv_unique) =
        build_histogram::<FnvHasher>(&args.input)?;
    log::info!(
        "FnvHasher: {} total syncmers, {} unique k-mers",
        fnv_total,
        fnv_unique
    );

    // Collect all frequencies present in either histogram
    let mut all_freqs: Vec<u64> = identity_histogram
        .keys()
        .chain(fnv_histogram.keys())
        .copied()
        .collect();
    all_freqs.sort();
    all_freqs.dedup();

    // Output TSV header
    println!("frequency\tidentity_count\tfnv_count");

    // Output collated histogram
    for freq in all_freqs {
        let identity_count = identity_histogram.get(&freq).copied().unwrap_or(0);
        let fnv_count = fnv_histogram.get(&freq).copied().unwrap_or(0);
        println!("{}\t{}\t{}", freq, identity_count, fnv_count);
    }

    Ok(())
}
