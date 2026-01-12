
use std::{fs::File, io::BufReader, path::PathBuf};

use clap::{Parser, Subcommand};
use noodles::fasta;

mod align;
mod index;
mod reads;
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
    Index {
        /// Path to reference FASTA
        fasta: PathBuf,

        /// Path to FASTQ file with reads to process
        fastq: PathBuf
    },
}

fn inner_main(cli: Cli) -> Result<(), error::ParallaxError> {
    match cli.command {
        Commands::Index { fasta, fastq } => {
            log::info!("Indexing reference from {}", fasta.display());
            let reader = File::open(&fasta).map(BufReader::new)?;
            let reader = fasta::io::Reader::new(reader);
            let index: index::Index<20, 15> = index::Index::try_from(reader)?;
            log::info!("Finished indexing {}", fasta.display());
            reads::process_reads(&index, fastq.to_str().unwrap())?;
        }
    }

    Ok(())
}

fn main() {
    env_logger::builder()
    .filter_level(log::LevelFilter::Info)
    .init();

    let cli = Cli::parse();

    if let Err(err) = inner_main(cli) {
        eprintln!("Error: {}", err);
        let mut err: &dyn std::error::Error = &err;
        while let Some(source) = err.source() {
            eprintln!("Caused by: {}", source);
            err = source;
        }
        std::process::exit(1);
    }
}
