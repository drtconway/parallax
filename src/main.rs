
use std::{fs::File, io::BufReader, path::PathBuf};

use clap::{Parser, Subcommand};
use noodles::fasta;

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

        /// Number of threads to use for alignment
        #[arg(short = 't', long, default_value = "4")]
        threads: usize,
    },
}

fn inner_main(cli: Cli) -> Result<(), error::ParallaxError> {
    match cli.command {
        Commands::Align { fasta, fastq, sam, threads } => {
            log::info!("Indexing reference from {}", fasta.display());
            let reader = File::open(&fasta).map(BufReader::new)?;
            let reader = fasta::io::Reader::new(reader);
            let index: index::Index<20, 15> = index::Index::try_from(reader)?;
            log::info!("Finished indexing {}", fasta.display());
            
            // Load reference into memory for efficient parallel access
            let reference = reference::InMemoryReference::load(&fasta)?;
            
            reads::process_reads_parallel(&index, &reference, fastq.to_str().unwrap(), sam.as_ref().map(|p| p.to_str().unwrap()), threads)?;
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
