use std::collections::{HashMap, HashSet};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::Args;
use parallax::error::Result;
use parallax::index::{Index, IndexHit};
use parallax::kmers::Kmer;
use parallax::reference::InMemoryReference;

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
    #[arg(long, default_value = "500")]
    pub max_seed_occurrences: usize,
}

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

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    writeln!(
        out,
        "read_id\tstrand\tkmer\tchrom\tref_pos\tread_pos\tkmer_multiplicity"
    )?;

    let (reader, _) = niffler::from_path(&args.reads)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    let reader = io::BufReader::new(reader);
    let mut fastq = noodles::fastq::io::Reader::new(reader);

    let mut loci_buf = Vec::new();

    for record in fastq.records() {
        let record = record?;
        let read_id = String::from_utf8_lossy(record.name()).into_owned();

        if !filter.is_empty() && !filter.contains(&read_id) {
            continue;
        }

        let seq: &[u8] = record.sequence().as_ref();
        let seq = seq.to_vec();

        for is_reverse in [false, true] {
            let strand_seq = if is_reverse {
                let mut buf = Vec::new();
                parallax::utils::sequence::reverse_complement_into(&seq, &mut buf);
                buf
            } else {
                seq.clone()
            };

            let strand_label = if is_reverse { "-" } else { "+" };

            // Collect raw per-kmer hits before any merging.
            // (kmer, chrom_id, ref_pos, read_pos, hit_count)
            let mut atoms: Vec<(u64, usize, usize, usize, u32)> = Vec::new();

            idx.find_seeds(&strand_seq, &mut |hit: IndexHit<'_>| {
                let IndexHit { query_pos, seed_kmer, loci, .. } = hit;
                let hit_count = loci.len();
                if hit_count > args.max_seed_occurrences {
                    return;
                }
                idx.unpack_loci(loci, &mut loci_buf);
                for &(chrom_id, ref_pos) in &loci_buf {
                    atoms.push((seed_kmer, chrom_id, ref_pos, query_pos, hit_count as u32));
                }
            });

            // kmer_multiplicity: how many times this kmer appears in this read's hit list.
            // Matches the read_frequency computed in gather_seeds_batched Phase 1b.
            let mut read_freq: HashMap<u64, u32> = HashMap::new();
            for &(kmer, _, _, _, _) in &atoms {
                *read_freq.entry(kmer).or_insert(0) += 1;
            }

            for (kmer, chrom_id, ref_pos, read_pos, _) in atoms {
                let kmer_str = Kmer::<20>::from(kmer).to_string();
                let chrom_name = reference.chrom_name(chrom_id);
                let multiplicity = read_freq.get(&kmer).copied().unwrap_or(1);
                writeln!(
                    out,
                    "{read_id}\t{strand_label}\t{kmer_str}\t{chrom_name}\t{ref_pos}\t{read_pos}\t{multiplicity}"
                )?;
            }
        }
    }

    out.flush()?;
    Ok(())
}
