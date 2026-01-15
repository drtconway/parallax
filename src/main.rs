
use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod align;
mod index;
mod metrics;
mod reads;
mod reference;
mod writer;
mod kmers;
mod error;
mod utils;

#[derive(Parser)]
#[command(name = "parallax")]
#[command(about = "Sequence indexing and alignment utilities", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create an index from a reference FASTA
    Align {
        /// Path to reference FASTA
        fasta: PathBuf,

        /// Path to FASTQ file with reads to process
        fastq: PathBuf,

        /// Path to output SAM file
        sam: Option<PathBuf>,

        /// Path to index directory (to save/load prebuilt index)
        #[arg(short = 'x', long)]
        index: Option<PathBuf>,

        /// Number of threads to use for alignment
        #[arg(short = 't', long, default_value = "4")]
        threads: usize,
    },
}

fn inner_main(cli: Cli) -> Result<(), error::ParallaxError> {
    match cli.command {
        Commands::Align { fasta, fastq, sam, index, threads } => {
            // Load reference into memory first
            let reference = reference::InMemoryReference::load(&fasta)?;
            
            // Either load or build the index
            // Check for chrom_info.json to verify index exists (not just empty directory)
            let idx: index::Index<20, 15> = if let Some(ref index_path) = index {
                if index_path.join("chrom_info.json").exists() {
                    log::info!("Loading index from {}", index_path.display());
                    index::Index::load(index_path)?
                } else {
                    log::info!("Building index from {}", fasta.display());
                    let idx = index::IndexBuilder::build_parallel(&reference, threads);
                    log::info!("Saving index to {}", index_path.display());
                    idx.save(index_path)?;
                    idx
                }
            } else {
                log::info!("Building index from {}", fasta.display());
                index::IndexBuilder::build_parallel(&reference, threads)
            };
            log::info!("Finished indexing {}", fasta.display());
            
            reads::process_reads_parallel(&idx, &reference, fastq.to_str().unwrap(), sam.as_ref().map(|p| p.to_str().unwrap()), threads)?;
        }
    }

    Ok(())
}

fn main() {
    env_logger::builder()
    .filter_level(log::LevelFilter::Info)
    .init();

    // Install metrics recorder
    let metrics_handle = metrics::SummaryRecorder::install()
        .expect("Failed to install metrics recorder");

    let cli = Cli::parse();

    let result = inner_main(cli);

    // Print metrics summary before exiting
    metrics_handle.print_summary();

    if let Err(err) = result {
        eprintln!("Error: {}", err);
        let mut err: &dyn std::error::Error = &err;
        while let Some(source) = err.source() {
            eprintln!("Caused by: {}", source);
            err = source;
        }
        std::process::exit(1);
    }
}
