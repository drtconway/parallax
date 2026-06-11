//! Select reads that have at least one unique seed in any of the provided BED regions.
//!
//! For each read in the input FASTQ, syncmer k-mers are extracted and looked up
//! in the index. A read is selected if any k-mer maps to exactly one reference
//! locus (a unique seed) and that locus falls within one of the BED intervals.
//! Both forward and reverse-complement strands are checked.

use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Args;
use noodles::fastq;

use parallax::index::{Index, PackedLocus};
use parallax::{
    error::Result,
    index::{self, BedRegions, fwd_index::FwdIndex},
    kmers::Kmer,
    utils::hasher::FnvHasher,
    utils::sequence::reverse_complement_into,
};

/// Arguments for the `select` subcommand.
#[derive(Args, Debug, Clone)]
pub struct SelectArgs {
    /// Path to prebuilt index directory
    pub index: PathBuf,

    /// BED file with regions of interest
    pub bed: PathBuf,

    /// Input FASTQ file (plain or gzip/bzip2/xz compressed)
    pub input: PathBuf,

    /// Output FASTQ file (uncompressed); use - for stdout
    #[arg(short = 'o', long, default_value = "-")]
    pub output: PathBuf,

    /// Use portable index format (Parquet instead of Feather)
    #[arg(long)]
    pub portable: bool,
}

/// Returns true if `pos` falls within any interval in the sorted, merged list.
fn in_regions(intervals: &[(usize, usize)], pos: usize) -> bool {
    let idx = intervals.partition_point(|&(start, _)| start <= pos);
    if idx == 0 {
        return false;
    }
    pos < intervals[idx - 1].1
}

/// Returns Some(region_idx) if pos lies within the region.
fn which_region(intervals: &[(usize, usize)], pos: usize) -> Option<usize> {
    let idx = intervals.partition_point(|&(start, _)| start <= pos);
    if idx > 0 && pos < intervals[idx - 1].1 {
        return Some(idx - 1);
    }
    None
}

/// Holds the index and per-call reusable buffers for read selection.
struct Selector<'a, const K: usize, const S: usize> {
    index: &'a FwdIndex<K, S>,
    chrom_names: Vec<String>,
    regions: BedRegions,
    /// Batch buffer reused across calls: (read_pos, kmer_value)
    kmer_batch: Vec<(usize, u64)>,
    /// Reverse-complement scratch buffer
    rc_buf: Vec<u8>,
}

impl<'a, const K: usize, const S: usize> Selector<'a, K, S> {
    fn new(index: &'a FwdIndex<K, S>, regions: BedRegions) -> Self {
        let chrom_names = index
            .all_chrom_info()
            .iter()
            .map(|c| c.name.clone())
            .collect();
        Selector {
            index,
            chrom_names,
            regions,
            kmer_batch: Vec::new(),
            rc_buf: Vec::new(),
        }
    }

    /// Fill `batch` from `seq` then run a prefetch-pipelined batch lookup.
    /// Returns true if any k-mer is uniquely placed in the target regions.
    ///
    /// Takes `batch` as an explicit `&mut` so callers can pass `self.kmer_batch`
    /// independently of any other field borrow (e.g. `self.rc_buf`).
    fn strand_has_hit(
        index: &FwdIndex<K, S>,
        chrom_names: &[String],
        regions: &BedRegions,
        batch: &mut Vec<(usize, u64)>,
        seq: &[u8],
    ) -> bool {
        batch.clear();
        Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(seq, [(); S], |pos, kmer| {
            batch.push((pos, kmer.0));
        });

        const BRIEF_MODE: bool = false;

        if BRIEF_MODE {
            let mut found = false;
            index.lookup_batch(batch, |_read_pos, _kmer_val, hit_count, loci| {
                if found || hit_count != 1 {
                    return;
                }
                let (chrom_id, ref_pos, _) = loci[0].unpack();
                if chrom_id >= chrom_names.len() {
                    return;
                }
                if let Some(intervals) = regions.get(&chrom_names[chrom_id]) {
                    if in_regions(intervals, ref_pos) {
                        found = true;
                    }
                }
            });
            found
        } else {
            let mut hits: Vec<Vec<f64>> = Vec::new();
            index.lookup_batch(batch, |_read_pos, kmer_val, hit_count, loci| {
                if hit_count > 10 {
                    return;
                }
                let w = 1.0 / (hit_count as f64);
                for loc in loci.iter() {
                    let (chrom_id, ref_pos, _) = loc.unpack();
                    if chrom_id >= chrom_names.len() {
                        return;
                    }
                    while hits.len() <= chrom_id {
                        hits.push(vec![]);
                    }
                    if let Some(intervals) = regions.get(&chrom_names[chrom_id]) {
                        if let Some(idx) = which_region(intervals, ref_pos) {
                            let chrom_hits = &mut hits[chrom_id];
                            while chrom_hits.len() <= idx {
                                chrom_hits.push(0.0);
                            }
                            chrom_hits[idx] += w;
                            let x = Kmer::<K>::from(kmer_val);
                            log::debug!("hit {}:{}-{} with {}", &chrom_names[chrom_id], intervals[idx].0, intervals[idx].1, x.to_string());
                        }
                    }
                }
            });
            let mut found = false;
            for (chrom_id, chrom_hits) in hits.iter().enumerate() {
                for (idx, count) in chrom_hits.iter().enumerate() {
                    if *count < 2.0 {
                        continue;
                    }
                    found = true;
                    let ivl = regions.get(&chrom_names[chrom_id]).unwrap()[idx];
                    log::debug!("hits on {}:{}-{}\t{:.1}", &chrom_names[chrom_id], ivl.0, ivl.1, count);
                }
            }
            if found {
                log::debug!(".");
            }
            found
        }
    }

    /// Returns `(fwd_hit, rev_hit)`. Skips the RC check if the forward strand
    /// already matches.
    fn check(&mut self, seq: &[u8]) -> (bool, bool) {
        let fwd = Self::strand_has_hit(
            self.index, &self.chrom_names, &self.regions, &mut self.kmer_batch, seq,
        );
        reverse_complement_into(seq, &mut self.rc_buf);
        let rev = Self::strand_has_hit(
            self.index, &self.chrom_names, &self.regions, &mut self.kmer_batch, &self.rc_buf,
        );
        (fwd, rev)
    }
}

pub fn run(args: SelectArgs) -> Result<()> {
    log::info!("Loading index from {}", args.index.display());
    let idx: FwdIndex<20, 15> = if args.portable {
        FwdIndex::load(&args.index)?
    } else {
        FwdIndex::load_feather(&args.index)?
    };

    let regions = index::load_bed_regions(&args.bed)?;
    let mut selector = Selector::<20, 15>::new(&idx, regions);

    log::info!("Reading reads from {}", args.input.display());
    let (decompressed, compression) = niffler::from_path(&args.input)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    if compression != niffler::Format::No {
        log::info!("Detected {:?} compression on input", compression);
    }
    let mut reader = fastq::io::Reader::new(std::io::BufReader::new(decompressed));

    let output: Box<dyn Write> = if args.output.as_os_str() == "-" {
        Box::new(std::io::stdout())
    } else {
        log::info!("Writing selected reads to {}", args.output.display());
        Box::new(std::fs::File::create(&args.output)?)
    };
    let mut writer = fastq::io::Writer::new(BufWriter::new(output));

    let mut total = 0usize;
    let mut selected = 0usize;
    let mut fwd_hit_count = 0usize;
    let mut rev_hit_count = 0usize;

    for result in reader.records() {
        let record = result?;
        total += 1;

        let seq: &[u8] = record.sequence().as_ref();
        let (fwd_hit, rev_hit) = selector.check(seq);

        if fwd_hit {
            fwd_hit_count += 1;
        }
        if rev_hit {
            rev_hit_count += 1;
        }
        if fwd_hit || rev_hit {
            writer.write_record(&record)?;
            selected += 1;
        }
        if total & 1023 == 0 {
            log::info!(
                "scanned {} records, with {} hits ({} fwd, {} rev)",
                total, selected, fwd_hit_count, rev_hit_count
            );
        }
    }

    log::info!(
        "Selected {}/{} reads with unique seeds in target regions",
        selected, total
    );

    Ok(())
}
