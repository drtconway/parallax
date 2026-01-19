
use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod align;
mod config;
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
    /// Align reads to a reference genome
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

        /// Path to BED file with regions of interest (only index these regions)
        #[arg(short = 'b', long)]
        bed: Option<PathBuf>,

        /// Number of threads to use for alignment
        #[arg(short = 't', long, default_value = "4")]
        threads: usize,

        /// Only use primary chromosomes (exclude ALT, unlocalized, unplaced contigs)
        #[arg(short = 'p', long)]
        primary_only: bool,

        /// Path to configuration file (TOML format)
        #[arg(short = 'c', long)]
        config: Option<PathBuf>,
    },

    /// Generate a template configuration file with documented defaults
    GenerateConfig {
        /// Output path for the config file (use - for stdout)
        #[arg(default_value = "parallax.toml")]
        output: PathBuf,
    },
}

fn inner_main(cli: Cli) -> Result<(), error::ParallaxError> {
    match cli.command {
        Commands::GenerateConfig { output } => {
            let template = config::generate_template();
            if output.as_os_str() == "-" {
                println!("{}", template);
            } else {
                std::fs::write(&output, &template)?;
                log::info!("Generated config template: {}", output.display());
            }
        }

        Commands::Align { fasta, fastq, sam, index, bed, threads, primary_only, config: config_path } => {
            // Load and initialize configuration
            let cfg = config::load(config_path.as_deref())
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
            config::init(cfg);

            // Load reference into memory first
            let reference = reference::InMemoryReference::load(&fasta, primary_only)?;
            
            // Either load or build the index
            // Check for chrom_info.json to verify index exists (not just empty directory)
            // Load BED regions if provided
            let bed_regions = if let Some(ref bed_path) = bed {
                Some(index::load_bed_regions(bed_path)?)
            } else {
                None
            };

            let idx: index::Index<20, 15> = if let Some(ref index_path) = index {
                if index_path.join("chrom_info.json").exists() {
                    log::info!("Loading index from {}", index_path.display());
                    index::Index::load(index_path)?
                } else {
                    log::info!("Building index from {}", fasta.display());
                    let idx = index::IndexBuilder::build_parallel(&reference, bed_regions.as_ref(), threads);
                    log::info!("Saving index to {}", index_path.display());
                    idx.save(index_path)?;
                    idx
                }
            } else {
                log::info!("Building index from {}", fasta.display());
                index::IndexBuilder::build_parallel(&reference, bed_regions.as_ref(), threads)
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
