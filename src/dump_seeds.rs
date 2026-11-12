use std::collections::{HashMap, HashSet};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::Args;
use parallax::error::Result;
use parallax::index::{Index, IndexHit};
use parallax::kmers::Kmer;
use parallax::reference::InMemoryReference;
use zip::write::SimpleFileOptions;

#[derive(Args, Debug, Clone)]
pub struct DumpSeedsArgs {
    /// Path to reference FASTA
    pub fasta: PathBuf,

    /// Path to reads file (FASTQ or FASTQ.gz)
    pub reads: PathBuf,

    /// Path to index directory
    #[arg(short = 'x', long)]
    pub index: Option<PathBuf>,

    /// Path to configuration file (TOML format)
    #[arg(short = 'c', long)]
    pub config: Option<PathBuf>,

    /// Only dump seeds for these read names (may be repeated)
    #[arg(long = "read", value_name = "NAME")]
    pub read_names: Vec<String>,

    /// Maximum k-mer occurrences (filters highly repetitive k-mers)
    #[arg(long, default_value = "50")]
    pub max_seed_occurrences: usize,

    /// Write each read's seeds as a separate compressed TSV entry inside this zip file
    #[arg(long, value_name = "PATH")]
    pub zip: Option<PathBuf>,
}

const HEADER: &str = "read_id\tstrand\tkmer\tchrom\tref_pos\tread_pos\tkmer_multiplicity\n";

pub fn run(args: DumpSeedsArgs) -> Result<()> {
    let reference = InMemoryReference::load(&args.fasta, false)?;

    let idx: std::sync::Arc<dyn Index> = if let Some(ref index_path) = args.index {
        parallax::index::load_index(index_path)?
    } else {
        log::info!("Building index from {}", args.fasta.display());
        let built: parallax::index::fwd_index::FwdIndex<20, 15> =
            parallax::index::fwd_index::FwdIndexBuilder::build_parallel(&reference, None, 4);
        std::sync::Arc::new(built)
    };

    let filter: HashSet<String> = args.read_names.into_iter().collect();

    let (reader, _) = niffler::from_path(&args.reads)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    let reader = io::BufReader::new(reader);
    let mut fastq = noodles::fastq::io::Reader::new(reader);

    let zip_options = SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);

    let mut zip_writer = args.zip.as_ref().map(|path| {
        let file = std::fs::File::create(path)
            .unwrap_or_else(|e| panic!("failed to create zip file {}: {e}", path.display()));
        zip::ZipWriter::new(file)
    });

    let mut stdout_out;
    let mut stdout_writer: Option<BufWriter<io::StdoutLock<'_>>> = if zip_writer.is_none() {
        let stdout = io::stdout();
        stdout_out = BufWriter::new(stdout.lock());
        writeln!(stdout_out, "{}", HEADER.trim_end())?;
        Some(stdout_out)
    } else {
        None
    };

    let mut loci_buf = Vec::new();

    let mut read_count = 0;
    for record in fastq.records() {
        read_count += 1;

        let record = record?;
        let read_id = String::from_utf8_lossy(record.name()).into_owned();

        if !filter.is_empty() && !filter.contains(&read_id) {
            continue;
        }

        let seq: &[u8] = record.sequence().as_ref();
        let seq = seq.to_vec();

        // Accumulate all rows for this read so we can write them as one zip entry.
        let mut read_rows: Vec<u8> = Vec::new();

        for is_reverse in [false, true] {
            let strand_seq = if is_reverse {
                let mut buf = Vec::new();
                parallax::utils::sequence::reverse_complement_into(&seq, &mut buf);
                buf
            } else {
                seq.clone()
            };

            let strand_label = if is_reverse { "-" } else { "+" };

            let mut atoms: Vec<(u64, usize, usize, usize, u32, usize)> = Vec::new();

            idx.find_seeds(&strand_seq, &mut |hit: IndexHit<'_>| {
                let IndexHit { query_pos, seed_kmer, loci, k } = hit;
                let hit_count = loci.len();
                if hit_count > args.max_seed_occurrences {
                    return;
                }
                idx.unpack_loci(loci, &mut loci_buf);
                for &(chrom_id, ref_pos) in &loci_buf {
                    atoms.push((seed_kmer, chrom_id, ref_pos, query_pos, hit_count as u32, k));
                }
            });

            let mut read_freq: HashMap<u64, u32> = HashMap::new();
            for &(kmer, _, _, _, _, _) in &atoms {
                *read_freq.entry(kmer).or_insert(0) += 1;
            }

            for (kmer, chrom_id, ref_pos, strand_local_pos, _, _) in atoms {
                let kmer_str = Kmer::<20>::from(kmer).to_string();
                let chrom_name = reference.chrom_name(chrom_id);
                let multiplicity = read_freq.get(&kmer).copied().unwrap_or(1);
                let row = format!(
                    "{read_id}\t{strand_label}\t{kmer_str}\t{chrom_name}\t{ref_pos}\t{strand_local_pos}\t{multiplicity}\n"
                );
                if let Some(ref mut w) = stdout_writer {
                    w.write_all(row.as_bytes())?;
                } else {
                    read_rows.extend_from_slice(row.as_bytes());
                }
            }
        }

        if let Some(ref mut zw) = zip_writer {
            let entry_name = format!("{read_id}.tsv");
            zw.start_file(&entry_name, zip_options)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            zw.write_all(HEADER.as_bytes())
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            zw.write_all(&read_rows)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        }

        if read_count % 100 == 0 {
            log::info!("Processed {} reads", read_count);
        }
    }

    if let Some(ref mut w) = stdout_writer {
        w.flush()?;
    }
    if let Some(zw) = zip_writer {
        zw.finish()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    }

    Ok(())
}
