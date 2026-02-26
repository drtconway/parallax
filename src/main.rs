
use std::path::PathBuf;
use std::io::Write;

use clap::{Args, Parser, Subcommand};

mod align;
mod annotate;
mod cluster;
mod config;
mod index;
mod metrics;
mod reads;
mod reference;
mod scores;
mod writer;
mod kmers;
mod error;
mod utils;

use writer::OutputFormat;

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

/// Full version string including git hash
const VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    "+",
    env!("GIT_VERSION")
);

/// Index building options shared between `index` and `align` commands.
#[derive(Args, Debug, Clone)]
pub struct IndexOptions {
    /// Path to BED file with regions of interest (only index these regions)
    #[arg(short = 'b', long)]
    pub bed: Option<PathBuf>,

    /// Number of threads to use
    #[arg(short = 't', long, default_value = "4")]
    pub threads: usize,

    /// Only use primary chromosomes (exclude ALT, unlocalized, unplaced contigs)
    #[arg(short = 'p', long)]
    pub primary_only: bool,

    /// Use portable index format (slower I/O, smaller files)
    #[arg(long)]
    pub portable: bool,
}

#[derive(Parser)]
#[command(name = "parallax")]
#[command(version = VERSION)]
#[command(about = "Sequence indexing and alignment utilities", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build an index for a reference genome
    Index {
        /// Path to reference FASTA
        fasta: PathBuf,

        /// Path to output index directory
        #[arg(short = 'o', long)]
        output: PathBuf,

        /// Index options
        #[command(flatten)]
        options: IndexOptions,
    },

    /// Align reads to a reference genome
    Align {
        /// Path to reference FASTA
        fasta: PathBuf,

        /// Path to reads file (FASTQ, FASTQ.gz, or unaligned BAM)
        reads: PathBuf,

        /// Path to output alignment file (SAM/BAM/CRAM)
        output: Option<PathBuf>,

        /// Output format (sam, bam, cram). If omitted, inferred from output
        /// file extension, or defaults to SAM.
        #[arg(short = 'O', long)]
        output_format: Option<OutputFormat>,

        /// Path to index directory (to load prebuilt index)
        #[arg(short = 'x', long)]
        index: Option<PathBuf>,

        /// Index options (used if building index on-the-fly)
        #[command(flatten)]
        index_options: IndexOptions,

        /// Path to configuration file (TOML format)
        #[arg(short = 'c', long)]
        config: Option<PathBuf>,

        /// Read group information
        #[command(flatten)]
        read_group: ReadGroup,
    },

    /// Annotate structural variant VCF records with library sequence identity
    Annotate {
        /// Path to library FASTA (mobile elements, viral genomes, etc.)
        library: PathBuf,

        /// Input VCF file (plain or bgzipped; use - for stdin)
        vcf: PathBuf,

        /// Path to a pre-built index directory (skips index building)
        #[arg(short = 'i', long)]
        index: Option<PathBuf>,

        /// Path to genome reference FASTA (needed for DEL/DUP annotation)
        #[arg(short = 'r', long)]
        reference: Option<PathBuf>,

        /// Output VCF path (default: stdout)
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,

        /// INFO field to use as query sequence instead of the ALT allele
        #[arg(long)]
        info_field: Option<String>,

        /// Minimum alignment quality (PHRED) to report a match
        #[arg(long, default_value = "20.0")]
        min_score: f64,

        /// Include CIGAR string in output
        #[arg(long)]
        emit_cigar: bool,

        /// Use portable index format (Parquet instead of Feather)
        #[arg(long)]
        portable: bool,

        /// Number of threads
        #[arg(short = 't', long, default_value = "4")]
        threads: usize,
    },

    /// Cluster repeat element instances and produce representative sequences
    Cluster {
        #[command(flatten)]
        args: cluster::ClusterArgs,
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
        Commands::Cluster { args } => {
            cluster::run(args)?;
        }

        Commands::GenerateConfig { output } => {
            let template = config::generate_template();
            if output.as_os_str() == "-" {
                println!("{}", template);
            } else {
                std::fs::write(&output, &template)?;
                log::info!("Generated config template: {}", output.display());
            }
        }

        Commands::Index { fasta, output, options } => {
            log::info!("Building index from {}", fasta.display());

            // Load reference
            let reference = reference::InMemoryReference::load(&fasta, options.primary_only)?;

            // Load BED regions if provided
            let bed_regions = if let Some(ref bed_path) = options.bed {
                Some(index::load_bed_regions(bed_path)?)
            } else {
                None
            };

            // Build index
            let idx: index::Index<20, 15> = index::IndexBuilder::build_parallel(
                &reference,
                bed_regions.as_ref(),
                options.threads,
            );

            // Save index
            log::info!("Saving index to {}", output.display());
            if options.portable {
                idx.save(&output)?;
            } else {
                idx.save_feather(&output)?;
            }
            log::info!("Index complete");
        }

        Commands::Annotate { library, vcf, index, reference, output, info_field, min_score, emit_cigar, portable, threads } => {
            annotate::run(annotate::AnnotateConfig {
                library_fasta: library,
                index_path: index,
                reference_fasta: reference,
                vcf_path: vcf,
                output,
                info_field,
                min_score,
                emit_cigar,
                portable,
                threads,
            })?;
        }

        Commands::Align { fasta, reads, output, output_format, index, index_options, config: config_path, read_group } => {
            // Load and initialize configuration
            let cfg = config::load(config_path.as_deref())
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
            config::init(cfg);

            // Determine output format: explicit flag > extension > SAM default
            let fmt = output_format
                .or_else(|| output.as_ref().and_then(|p| OutputFormat::from_path(p)))
                .unwrap_or(OutputFormat::Sam);
            log::info!("Output format: {}", fmt);

            // Load reference into memory first
            let reference = reference::InMemoryReference::load(&fasta, index_options.primary_only)?;
            
            // Either load or build the index
            let idx: index::Index<20, 15> = if let Some(ref index_path) = index {
                if index_path.join("chrom_info.json").exists() {
                    log::info!("Loading index from {}", index_path.display());
                    if index_options.portable {
                        index::Index::load(index_path)?
                    } else {
                        index::Index::load_feather(index_path)?
                    }
                } else {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Index not found at {}. Use 'parallax index' to build it first.", index_path.display())
                    ).into());
                }
            } else {
                // Build index on-the-fly
                log::info!("Building index from {}", fasta.display());
                let bed_regions = if let Some(ref bed_path) = index_options.bed {
                    Some(index::load_bed_regions(bed_path)?)
                } else {
                    None
                };
                index::IndexBuilder::build_parallel(&reference, bed_regions.as_ref(), index_options.threads)
            };
            log::info!("Finished indexing {}", fasta.display());

            // Validate that the index and reference are compatible
            idx.validate_reference(&reference)?;

            let rg_header = read_group.to_header_line();
            reads::process_reads_parallel(&idx, &reference, reads.to_str().unwrap(), output.as_ref().map(|p| p.to_str().unwrap()), index_options.threads, command_line, rg_header.as_deref(), fmt)?;
        }
    }

    Ok(())
}

fn main() {
    // Capture command line before clap consumes it
    let command_line: String = std::env::args().collect::<Vec<_>>().join(" ");

    env_logger::builder()
    .filter_level(log::LevelFilter::Info)
    .format(|buf, record| {
        writeln!(
            buf,
            "[{} {:5}] {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            record.level(),
            record.args()
        )
    })
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
