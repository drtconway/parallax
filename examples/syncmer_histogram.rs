use std::{collections::HashMap, path::PathBuf};

use clap::Parser;
use parallax::kmers::Kmer;
use parallax::utils::hasher::{FnvHasher, Hasher, IdentityHasher, Splitmix64Hasher};

#[derive(Parser, Debug)]
#[command(name = "syncmer-histogram")]
#[command(
    about = "Compute syncmer k-mer frequency histograms comparing IdentityHasher vs FnvHasher vs Splitmix64Hasher"
)]
struct Args {
    /// Value of k for k-mers
    #[arg(short, long, default_value = "20")]
    k: usize,

    /// Value of s for syncmers
    #[arg(short, long, default_value = "15")]
    s: usize,

    /// Input FASTA file
    input: PathBuf,
}

/// Build k-mer frequency histogram for syncmers using the given hasher.
/// Returns (histogram, total_syncmers, unique_kmers)
fn build_histogram<const K: usize, const S: usize, H: Hasher>(
    fasta_path: &PathBuf,
) -> std::io::Result<(HashMap<u64, u64>, u64, usize)> {
    log::info!(
        "Pass: Computing syncmer histogram with K={} S={} H={}",
        K,
        S,
        H::NAME
    );
    let mut reader = noodles::fasta::io::reader::Builder::default().build_from_path(fasta_path)?;

    let mut kmer_counts: HashMap<u64, u64> = HashMap::new();

    let now = std::time::Instant::now();

    for record in reader.records() {
        let record = record?;
        let seq = record.sequence().as_ref();

        Kmer::<K>::kmerize_open_syncmers_fwd::<S, H, _, _>(seq, [(); S], |_pos, kmer| {
            *kmer_counts.entry(kmer.0).or_insert(0) += 1;
        });
    }

    let elapsed = now.elapsed();

    // Build frequency histogram (count -> how many k-mers have that count)
    let mut histogram: HashMap<u64, u64> = HashMap::new();
    for &count in kmer_counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }

    let total: u64 = kmer_counts.values().sum();
    let unique = kmer_counts.len();

    log::info!(
        "{}: {} total syncmers, {} unique k-mers, computed in {:.2?}",
        H::NAME,
        total,
        unique,
        elapsed
    );

    Ok((histogram, total, unique))
}

fn establish_s<const K: usize, H: Hasher>(
    s: usize,
    fasta_path: &PathBuf,
) -> std::io::Result<(HashMap<u64, u64>, u64, usize)> {
    match s {
        10 => build_histogram::<K, 10, H>(fasta_path),
        11 => build_histogram::<K, 11, H>(fasta_path),
        12 => build_histogram::<K, 12, H>(fasta_path),
        13 => build_histogram::<K, 13, H>(fasta_path),
        14 => build_histogram::<K, 14, H>(fasta_path),
        15 => build_histogram::<K, 15, H>(fasta_path),
        16 => build_histogram::<K, 16, H>(fasta_path),
        17 => build_histogram::<K, 17, H>(fasta_path),
        18 => build_histogram::<K, 18, H>(fasta_path),
        19 => build_histogram::<K, 19, H>(fasta_path),
        20 => build_histogram::<K, 20, H>(fasta_path),
        21 => build_histogram::<K, 21, H>(fasta_path),
        22 => build_histogram::<K, 22, H>(fasta_path),
        23 => build_histogram::<K, 23, H>(fasta_path),
        24 => build_histogram::<K, 24, H>(fasta_path),
        25 => build_histogram::<K, 25, H>(fasta_path),
        26 => build_histogram::<K, 26, H>(fasta_path),
        27 => build_histogram::<K, 27, H>(fasta_path),
        28 => build_histogram::<K, 28, H>(fasta_path),
        29 => build_histogram::<K, 29, H>(fasta_path),
        30 => build_histogram::<K, 30, H>(fasta_path),
        _ => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Unsupported S value",
        )),
    }
}

fn establish_k<H: Hasher>(
    k: usize,
    s: usize,
    fasta_path: &PathBuf,
) -> std::io::Result<(HashMap<u64, u64>, u64, usize)> {
    match k {
        15 => establish_s::<15, H>(s, fasta_path),
        16 => establish_s::<16, H>(s, fasta_path),
        17 => establish_s::<17, H>(s, fasta_path),
        18 => establish_s::<18, H>(s, fasta_path),
        19 => establish_s::<19, H>(s, fasta_path),
        20 => establish_s::<20, H>(s, fasta_path),
        21 => establish_s::<21, H>(s, fasta_path),
        22 => establish_s::<22, H>(s, fasta_path),
        23 => establish_s::<23, H>(s, fasta_path),
        24 => establish_s::<24, H>(s, fasta_path),
        25 => establish_s::<25, H>(s, fasta_path),
        26 => establish_s::<26, H>(s, fasta_path),
        27 => establish_s::<27, H>(s, fasta_path),
        28 => establish_s::<28, H>(s, fasta_path),
        29 => establish_s::<29, H>(s, fasta_path),
        30 => establish_s::<30, H>(s, fasta_path),
        _ => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Unsupported K value",
        )),
    }
}

fn main() -> std::io::Result<()> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    if args.s >= args.k {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "S must be less than K",
        ));
    }

    let (identity_histogram, _, _) =
        establish_k::<IdentityHasher>(args.k, args.s, &args.input)?;
    let (fnv_histogram, _, _) =
        establish_k::<FnvHasher>(args.k, args.s, &args.input)?;
    let (splitmix_histogram, _, _) =
        establish_k::<Splitmix64Hasher>(args.k, args.s, &args.input)?;
        
    // Collect all frequencies present in either histogram
    let mut all_freqs: Vec<u64> = identity_histogram
        .keys()
        .chain(fnv_histogram.keys())
        .chain(splitmix_histogram.keys())
        .copied()
        .collect();
    all_freqs.sort();
    all_freqs.dedup();

    // Output TSV header
    println!("frequency\tidentity_count\tfnv_count\tsplitmix_count");

    // Output collated histogram
    for freq in all_freqs {
        let identity_count = identity_histogram.get(&freq).copied().unwrap_or(0);
        let fnv_count = fnv_histogram.get(&freq).copied().unwrap_or(0);
        let splitmix_count = splitmix_histogram.get(&freq).copied().unwrap_or(0);
        println!("{}\t{}\t{}\t{}", freq, identity_count, fnv_count, splitmix_count);
    }

    Ok(())
}
