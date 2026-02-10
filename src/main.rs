
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

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

/// Read group information for SAM/BAM output.
///
/// Read groups identify subsets of reads sharing common properties
/// (e.g., same sequencing run, library, or sample).
#[derive(Args, Debug, Clone, Default)]
pub struct ReadGroup {
    /// Read group identifier (required if any RG field is set)
    #[arg(long = "rg-id")]
    pub id: Option<String>,

    /// Sample name
    #[arg(long = "rg-sm")]
    pub sample: Option<String>,

    /// Library identifier
    #[arg(long = "rg-lb")]
    pub library: Option<String>,

    /// Platform/technology (e.g., ILLUMINA, PACBIO, ONT)
    #[arg(long = "rg-pl")]
    pub platform: Option<String>,

    /// Platform unit (e.g., flowcell-barcode.lane)
    #[arg(long = "rg-pu")]
    pub platform_unit: Option<String>,

    /// Sequencing center
    #[arg(long = "rg-cn")]
    pub center: Option<String>,

    /// Description
    #[arg(long = "rg-ds")]
    pub description: Option<String>,

    /// Run date (ISO 8601)
    #[arg(long = "rg-dt")]
    pub date: Option<String>,

    /// Platform model
    #[arg(long = "rg-pm")]
    pub platform_model: Option<String>,
}

impl ReadGroup {
    /// Returns true if any read group field is set.
    pub fn is_set(&self) -> bool {
        self.id.is_some()
            || self.sample.is_some()
            || self.library.is_some()
            || self.platform.is_some()
            || self.platform_unit.is_some()
            || self.center.is_some()
            || self.description.is_some()
            || self.date.is_some()
            || self.platform_model.is_some()
    }

    /// Format as SAM @RG header line. Returns None if no fields are set.
    pub fn to_header_line(&self) -> Option<String> {
        if !self.is_set() {
            return None;
        }

        let mut parts = vec!["@RG".to_string()];

        // ID is required - default to "1" if not specified but other fields are
        let id = self.id.as_deref().unwrap_or("1");
        parts.push(format!("ID:{}", id));

        if let Some(ref sm) = self.sample {
            parts.push(format!("SM:{}", sm));
        }
        if let Some(ref lb) = self.library {
            parts.push(format!("LB:{}", lb));
        }
        if let Some(ref pl) = self.platform {
            parts.push(format!("PL:{}", pl));
        }
        if let Some(ref pu) = self.platform_unit {
            parts.push(format!("PU:{}", pu));
        }
        if let Some(ref cn) = self.center {
            parts.push(format!("CN:{}", cn));
        }
        if let Some(ref ds) = self.description {
            parts.push(format!("DS:{}", ds));
        }
        if let Some(ref dt) = self.date {
            parts.push(format!("DT:{}", dt));
        }
        if let Some(ref pm) = self.platform_model {
            parts.push(format!("PM:{}", pm));
        }

        Some(parts.join("\t"))
    }

    /// Get the read group ID for per-read RG:Z: tags.
    pub fn id(&self) -> Option<&str> {
        if self.is_set() {
            Some(self.id.as_deref().unwrap_or("1"))
        } else {
            None
        }
    }
}

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

        /// Use portable index format (slower I/O, smaller files)
        #[arg(long)]
        portable: bool,

        /// Path to configuration file (TOML format)
        #[arg(short = 'c', long)]
        config: Option<PathBuf>,

        /// Read group information
        #[command(flatten)]
        read_group: ReadGroup,
    },

    /// Generate a template configuration file with documented defaults
    GenerateConfig {
        /// Output path for the config file (use - for stdout)
        #[arg(default_value = "parallax.toml")]
        output: PathBuf,
    },
}

fn inner_main(cli: Cli, command_line: &str) -> Result<(), error::ParallaxError> {
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

        Commands::Align { fasta, fastq, sam, index, bed, threads, primary_only, portable, config: config_path, read_group } => {
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
                    if portable {
                        index::Index::load(index_path)?
                    } else {
                        index::Index::load_feather(index_path)?
                    }
                } else {
                    log::info!("Building index from {}", fasta.display());
                    let idx = index::IndexBuilder::build_parallel(&reference, bed_regions.as_ref(), threads);
                    log::info!("Saving index to {}", index_path.display());
                    if portable {
                        idx.save(index_path)?;
                    } else {
                        idx.save_feather(index_path)?;
                    }
                    idx
                }
            } else {
                log::info!("Building index from {}", fasta.display());
                index::IndexBuilder::build_parallel(&reference, bed_regions.as_ref(), threads)
            };
            log::info!("Finished indexing {}", fasta.display());
            
            let rg_header = read_group.to_header_line();
            reads::process_reads_parallel(&idx, &reference, fastq.to_str().unwrap(), sam.as_ref().map(|p| p.to_str().unwrap()), threads, command_line, rg_header.as_deref())?;
        }
    }

    Ok(())
}

fn main() {
    // Capture command line before clap consumes it
    let command_line: String = std::env::args().collect::<Vec<_>>().join(" ");

    env_logger::builder()
    .filter_level(log::LevelFilter::Info)
    .init();

    // Install metrics recorder
    let metrics_handle = metrics::SummaryRecorder::install()
        .expect("Failed to install metrics recorder");

    let cli = Cli::parse();

    let result = inner_main(cli, &command_line);

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
