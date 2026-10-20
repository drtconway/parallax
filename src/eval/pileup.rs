//! Per-chromosome pileup summary for truth-free alignment evaluation.
//!
//! Walks a coordinate-sorted, indexed BAM over one chromosome and counts the
//! A/C/G/T/deletion alleles (and insertions) at every reference position, for
//! two classes of alignment: all primary and supplementary alignments, and
//! those with MAPQ >= `min_mapq`.  Each position is folded into histograms
//! rather than written out, so the output stays small whatever the depth.
//!
//! Outputs, all prefixed with `--output-prefix`:
//! - `.hist.tsv`: histograms of depth, and (at depth >= `min_depth`) of major,
//!   minor and insertion allele fractions, by MAPQ class, stratum and depth bin.
//! - `.totals.tsv`: per MAPQ class and stratum totals (positions, depth,
//!   reference/mismatch/deletion/insertion counts, non-reference majors,
//!   candidate heterozygous sites).
//! - `.sites.bed`: positions whose minor allele fraction is at least `het_af`.
//!   On a haploid sample aligned to its own assembly (CHM13) these are almost
//!   all misplaced paralogous reads.
//! - `.windows.tsv`: mean depth in fixed windows, to find spikes and holes.

use std::collections::VecDeque;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::Args;
use noodles::core::Region;
use noodles::sam::alignment::record::cigar::op::Kind;

use crate::error::Result;
use crate::reference::Reference;

/// Arguments for the `eval-pileup` subcommand.
#[derive(Args, Debug, Clone)]
pub struct PileupArgs {
    /// Coordinate-sorted, indexed BAM file
    pub bam: PathBuf,

    /// Reference FASTA (with .fai)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Chromosome to summarise
    #[arg(long)]
    pub chrom: String,

    /// Strata TSV (name<TAB>bed path, no header).  A position belongs to the
    /// first stratum listed that contains it, or to "other".
    #[arg(long)]
    pub strata: Option<PathBuf>,

    /// MAPQ threshold for the high-confidence class
    #[arg(long, default_value = "20")]
    pub min_mapq: u8,

    /// Minimum depth for allele fraction histograms and candidate sites
    #[arg(long, default_value = "10")]
    pub min_depth: u32,

    /// Minimum minor allele fraction for a candidate heterozygous site
    #[arg(long, default_value = "0.2")]
    pub het_af: f64,

    /// Window size for the depth track
    #[arg(long, default_value = "10000")]
    pub window: u64,

    /// Prefix for output files
    #[arg(short = 'o', long)]
    pub output_prefix: String,
}

const NUM_CLASSES: usize = 2;
const NUM_ALLELES: usize = 5;
const DEL: usize = 4;
const ALLELE_CHARS: [u8; NUM_ALLELES] = *b"ACGT-";
/// Depths above this share the last depth histogram bin.
const MAX_DEPTH: usize = 1000;
const AF_BINS: usize = 100;
/// Depth bins for allele fraction histograms: [d, 2d), [2d, 4d), ... with
/// d = `min_depth`; the last bin is open-ended.
const NUM_DEPTH_BINS: usize = 5;

fn allele_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[derive(Clone, Copy, Default, Debug)]
struct Counts {
    alleles: [u32; NUM_ALLELES],
    /// Insertions immediately after this position.
    ins: u32,
}

type Column = [Counts; NUM_CLASSES];

/// Stratum membership for one chromosome, queried at non-decreasing positions.
pub struct Strata {
    names: Vec<String>,
    /// Sorted, merged half-open intervals per stratum.
    intervals: Vec<Vec<(u64, u64)>>,
    cursors: Vec<usize>,
}

impl Strata {
    pub fn none() -> Self {
        Self::from_intervals(Vec::new(), Vec::new())
    }

    fn from_intervals(names: Vec<String>, mut intervals: Vec<Vec<(u64, u64)>>) -> Self {
        for ivs in intervals.iter_mut() {
            ivs.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ivs.len());
            for &(s, e) in ivs.iter() {
                match merged.last_mut() {
                    Some(last) if s <= last.1 => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            *ivs = merged;
        }
        let cursors = vec![0; names.len()];
        Strata { names, intervals, cursors }
    }

    /// Load the BEDs listed in `tsv`, keeping only intervals on `chrom`.
    pub fn load(tsv: &Path, chrom: &str) -> Result<Self> {
        let mut names = Vec::new();
        let mut intervals = Vec::new();
        for line in BufReader::new(std::fs::File::open(tsv)?).lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split('\t');
            let (Some(name), Some(bed)) = (fields.next(), fields.next()) else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("bad strata line (want name<TAB>bed): {}", line),
                )
                .into());
            };
            if name == "name" && bed == "bed" {
                continue; // header
            }
            names.push(name.to_string());
            intervals.push(read_bed(Path::new(bed), chrom)?);
        }
        Ok(Self::from_intervals(names, intervals))
    }

    /// Number of strata including the catch-all "other".
    pub fn len(&self) -> usize {
        self.names.len() + 1
    }

    pub fn name(&self, i: usize) -> &str {
        self.names.get(i).map(String::as_str).unwrap_or("other")
    }

    /// Stratum of `pos`.  Successive calls must have non-decreasing `pos`.
    fn classify(&mut self, pos: u64) -> usize {
        for (i, ivs) in self.intervals.iter().enumerate() {
            let c = &mut self.cursors[i];
            while *c < ivs.len() && ivs[*c].1 <= pos {
                *c += 1;
            }
            if *c < ivs.len() && ivs[*c].0 <= pos {
                return i;
            }
        }
        self.names.len()
    }
}

fn read_bed(path: &Path, chrom: &str) -> Result<Vec<(u64, u64)>> {
    let (reader, _format) = niffler::from_path(path).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidInput, format!("{}: {}", path.display(), e))
    })?;
    let mut intervals = Vec::new();
    for line in BufReader::new(reader).lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let mut fields = line.split('\t');
        if fields.next() != Some(chrom) {
            continue;
        }
        let parse = |f: Option<&str>| f.and_then(|s| s.trim().parse::<u64>().ok());
        if let (Some(s), Some(e)) = (parse(fields.next()), parse(fields.next())) {
            if s < e {
                intervals.push((s, e));
            }
        }
    }
    Ok(intervals)
}

/// Per MAPQ class and stratum totals.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Totals {
    /// Positions with a non-N reference base.
    pub positions: u64,
    /// Positions with depth >= `min_depth`.
    pub called_positions: u64,
    pub depth_sum: u64,
    pub ref_sum: u64,
    pub mismatch_sum: u64,
    pub del_sum: u64,
    pub ins_sum: u64,
    /// Called positions whose major allele is not the reference base.
    pub major_nonref: u64,
    /// Called positions whose minor allele fraction is >= `het_af`.
    pub candidate_sites: u64,
}

/// Histograms and totals, indexed by [class][stratum][...].
pub struct Stats {
    num_strata: usize,
    depth: Vec<u64>,
    major_af: Vec<u64>,
    minor_af: Vec<u64>,
    ins_af: Vec<u64>,
    totals: Vec<Totals>,
}

impl Stats {
    fn new(num_strata: usize) -> Self {
        let cs = NUM_CLASSES * num_strata;
        Stats {
            num_strata,
            depth: vec![0; cs * (MAX_DEPTH + 1)],
            major_af: vec![0; cs * NUM_DEPTH_BINS * AF_BINS],
            minor_af: vec![0; cs * NUM_DEPTH_BINS * AF_BINS],
            ins_af: vec![0; cs * NUM_DEPTH_BINS * AF_BINS],
            totals: vec![Totals::default(); cs],
        }
    }

    fn cs(&self, class: usize, stratum: usize) -> usize {
        class * self.num_strata + stratum
    }

    fn af_index(&self, class: usize, stratum: usize, depth_bin: usize, af: f64) -> usize {
        let bin = ((af * AF_BINS as f64) as usize).min(AF_BINS - 1);
        (self.cs(class, stratum) * NUM_DEPTH_BINS + depth_bin) * AF_BINS + bin
    }

    pub fn totals(&self, class: usize, stratum: usize) -> Totals {
        self.totals[self.cs(class, stratum)]
    }
}

/// Thresholds controlling how positions are summarised.
#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub min_depth: u32,
    pub het_af: f64,
    pub window: u64,
}

/// Streaming pileup over one chromosome.  Alignments must be added in order
/// of start position; columns are summarised as soon as no later alignment
/// can reach them.
pub struct Pileup<W: Write> {
    chrom: String,
    reference: Vec<u8>,
    params: Params,
    class_names: [String; NUM_CLASSES],
    strata: Strata,
    /// Columns for positions `next`, `next + 1`, ...
    buf: VecDeque<Column>,
    next: u64,
    stats: Stats,
    window_start: u64,
    window_positions: u64,
    window_depth: [u64; NUM_CLASSES],
    sites: W,
    windows: W,
}

impl<W: Write> Pileup<W> {
    pub fn new(
        chrom: String,
        reference: Vec<u8>,
        params: Params,
        class_names: [String; NUM_CLASSES],
        strata: Strata,
        mut sites: W,
        mut windows: W,
    ) -> io::Result<Self> {
        writeln!(
            sites,
            "#chrom\tstart\tend\tmapq_class\tstratum\tdepth\tref\tmajor\tmajor_count\tminor\tminor_count"
        )?;
        writeln!(windows, "chrom\tstart\tend\tpositions\t{}\t{}",
            format_args!("mean_depth_{}", class_names[0]),
            format_args!("mean_depth_{}", class_names[1]))?;
        let stats = Stats::new(strata.len());
        Ok(Pileup {
            chrom,
            reference,
            params,
            class_names,
            strata,
            buf: VecDeque::new(),
            next: 0,
            stats,
            window_start: 0,
            window_positions: 0,
            window_depth: [0; NUM_CLASSES],
            sites,
            windows,
        })
    }

    /// Add one alignment starting at 0-based `start`.  `classes[c]` says
    /// whether it counts towards MAPQ class `c`.
    pub fn add(
        &mut self,
        start: u64,
        ops: &[(Kind, usize)],
        seq: &[u8],
        classes: [bool; NUM_CLASSES],
    ) -> io::Result<()> {
        if start < self.next {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("alignments not coordinate sorted at {}:{}", self.chrom, start + 1),
            ));
        }
        self.flush_to(start)?;

        let span: u64 = ops
            .iter()
            .filter(|(k, _)| k.consumes_reference())
            .map(|&(_, n)| n as u64)
            .sum();
        let needed = (start + span - self.next) as usize;
        if self.buf.len() < needed {
            self.buf.resize(needed, Column::default());
        }

        let mut r = (start - self.next) as usize;
        let mut q = 0usize;
        for &(kind, n) in ops {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for i in 0..n {
                        if let Some(a) = seq.get(q + i).copied().and_then(allele_index) {
                            let col = &mut self.buf[r + i];
                            for c in 0..NUM_CLASSES {
                                if classes[c] {
                                    col[c].alleles[a] += 1;
                                }
                            }
                        }
                    }
                    r += n;
                    q += n;
                }
                Kind::Deletion => {
                    for i in 0..n {
                        let col = &mut self.buf[r + i];
                        for c in 0..NUM_CLASSES {
                            if classes[c] {
                                col[c].alleles[DEL] += 1;
                            }
                        }
                    }
                    r += n;
                }
                Kind::Insertion => {
                    if r > (start - self.next) as usize {
                        let col = &mut self.buf[r - 1];
                        for c in 0..NUM_CLASSES {
                            if classes[c] {
                                col[c].ins += 1;
                            }
                        }
                    }
                    q += n;
                }
                Kind::Skip => r += n,
                Kind::SoftClip => q += n,
                Kind::HardClip | Kind::Pad => {}
            }
        }
        Ok(())
    }

    /// Summarise every position before `pos`.
    fn flush_to(&mut self, pos: u64) -> io::Result<()> {
        let end = pos.min(self.reference.len() as u64);
        while self.next < end {
            let col = self.buf.pop_front().unwrap_or_default();
            self.emit(self.next, &col)?;
            self.next += 1;
        }
        Ok(())
    }

    fn emit(&mut self, pos: u64, col: &Column) -> io::Result<()> {
        if pos >= self.window_start + self.params.window {
            self.write_window()?;
            self.window_start = pos - pos % self.params.window;
        }
        let ref_base = self.reference[pos as usize];
        let Some(ref_idx) = allele_index(ref_base) else {
            return Ok(());
        };
        let stratum = self.strata.classify(pos);
        self.window_positions += 1;

        for c in 0..NUM_CLASSES {
            let counts = &col[c];
            let depth: u32 = counts.alleles.iter().sum();
            self.window_depth[c] += depth as u64;
            let cs = self.stats.cs(c, stratum);
            self.stats.depth[cs * (MAX_DEPTH + 1) + (depth as usize).min(MAX_DEPTH)] += 1;

            let ref_count = counts.alleles[ref_idx];
            let del = counts.alleles[DEL];
            let t = &mut self.stats.totals[cs];
            t.positions += 1;
            t.depth_sum += depth as u64;
            t.ref_sum += ref_count as u64;
            t.del_sum += del as u64;
            t.mismatch_sum += (depth - ref_count - del) as u64;
            t.ins_sum += counts.ins as u64;

            if depth < self.params.min_depth || depth == 0 {
                continue;
            }
            t.called_positions += 1;

            // Top two alleles; ties go to the reference, then to lower index.
            let mut order: [usize; NUM_ALLELES] = [0, 1, 2, 3, 4];
            order.sort_by_key(|&a| (std::cmp::Reverse(counts.alleles[a]), a != ref_idx, a));
            let (major, minor) = (order[0], order[1]);
            let major_af = counts.alleles[major] as f64 / depth as f64;
            let minor_af = counts.alleles[minor] as f64 / depth as f64;
            let ins_af = (counts.ins as f64 / depth as f64).min(1.0);
            if major != ref_idx {
                t.major_nonref += 1;
            }
            let is_candidate = minor_af >= self.params.het_af;
            if is_candidate {
                t.candidate_sites += 1;
            }

            let ratio = depth / self.params.min_depth.max(1);
            let depth_bin = (ratio.max(1).ilog2() as usize).min(NUM_DEPTH_BINS - 1);
            let i = self.stats.af_index(c, stratum, depth_bin, major_af);
            self.stats.major_af[i] += 1;
            let i = self.stats.af_index(c, stratum, depth_bin, minor_af);
            self.stats.minor_af[i] += 1;
            let i = self.stats.af_index(c, stratum, depth_bin, ins_af);
            self.stats.ins_af[i] += 1;

            if is_candidate {
                writeln!(
                    self.sites,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    self.chrom,
                    pos,
                    pos + 1,
                    self.class_names[c],
                    self.strata.name(stratum),
                    depth,
                    ref_base.to_ascii_uppercase() as char,
                    ALLELE_CHARS[major] as char,
                    counts.alleles[major],
                    ALLELE_CHARS[minor] as char,
                    counts.alleles[minor],
                )?;
            }
        }
        Ok(())
    }

    fn write_window(&mut self) -> io::Result<()> {
        let end = (self.window_start + self.params.window).min(self.reference.len() as u64);
        if end > self.window_start {
            let mean = |d: u64| {
                if self.window_positions > 0 {
                    d as f64 / self.window_positions as f64
                } else {
                    0.0
                }
            };
            writeln!(
                self.windows,
                "{}\t{}\t{}\t{}\t{:.3}\t{:.3}",
                self.chrom,
                self.window_start,
                end,
                self.window_positions,
                mean(self.window_depth[0]),
                mean(self.window_depth[1]),
            )?;
        }
        self.window_positions = 0;
        self.window_depth = [0; NUM_CLASSES];
        Ok(())
    }

    /// Summarise the remaining positions and return the statistics and sinks.
    pub fn finish(mut self) -> io::Result<(Stats, Strata, W, W)> {
        self.flush_to(self.reference.len() as u64)?;
        self.write_window()?;
        self.sites.flush()?;
        self.windows.flush()?;
        Ok((self.stats, self.strata, self.sites, self.windows))
    }
}

fn write_hist(
    out: &mut impl Write,
    chrom: &str,
    class_names: &[String; NUM_CLASSES],
    strata: &Strata,
    stats: &Stats,
    min_depth: u32,
) -> io::Result<()> {
    writeln!(out, "chrom\tmapq_class\tstratum\tkind\tdepth_bin\tbin\tcount")?;
    for c in 0..NUM_CLASSES {
        for s in 0..stats.num_strata {
            let cs = stats.cs(c, s);
            let base = cs * (MAX_DEPTH + 1);
            for d in 0..=MAX_DEPTH {
                let n = stats.depth[base + d];
                if n > 0 {
                    writeln!(out, "{}\t{}\t{}\tdepth\tNA\t{}\t{}",
                        chrom, class_names[c], strata.name(s), d, n)?;
                }
            }
            for (kind, hist) in [
                ("major_af", &stats.major_af),
                ("minor_af", &stats.minor_af),
                ("ins_af", &stats.ins_af),
            ] {
                for db in 0..NUM_DEPTH_BINS {
                    let depth_lo = min_depth.max(1) << db;
                    for b in 0..AF_BINS {
                        let n = hist[(cs * NUM_DEPTH_BINS + db) * AF_BINS + b];
                        if n > 0 {
                            writeln!(out, "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
                                chrom, class_names[c], strata.name(s), kind,
                                depth_lo, b as f64 / AF_BINS as f64, n)?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn write_totals(
    out: &mut impl Write,
    chrom: &str,
    class_names: &[String; NUM_CLASSES],
    strata: &Strata,
    stats: &Stats,
) -> io::Result<()> {
    writeln!(
        out,
        "chrom\tmapq_class\tstratum\tpositions\tcalled_positions\tdepth_sum\tref_sum\tmismatch_sum\tdel_sum\tins_sum\tmajor_nonref\tcandidate_sites"
    )?;
    for c in 0..NUM_CLASSES {
        for s in 0..stats.num_strata {
            let t = stats.totals(c, s);
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chrom, class_names[c], strata.name(s), t.positions, t.called_positions,
                t.depth_sum, t.ref_sum, t.mismatch_sum, t.del_sum, t.ins_sum,
                t.major_nonref, t.candidate_sites,
            )?;
        }
    }
    Ok(())
}

fn create(path: String) -> io::Result<BufWriter<std::fs::File>> {
    Ok(BufWriter::new(std::fs::File::create(path)?))
}

pub fn run(args: PileupArgs) -> Result<()> {
    let mut reference = Reference::open(&args.reference)?;
    let Some((chrom_idx, chrom_len)) = reference
        .chromosomes()
        .enumerate()
        .find(|(_, (name, _))| *name == args.chrom)
        .map(|(i, (_, len))| (i, len))
    else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("chromosome {} not in {}", args.chrom, args.reference.display()),
        )
        .into());
    };
    let mut ref_seq = Vec::new();
    reference.get_seq_into(chrom_idx, 0, chrom_len as usize, &mut ref_seq)?;
    ref_seq.make_ascii_uppercase();

    let strata = match &args.strata {
        Some(path) => Strata::load(path, &args.chrom)?,
        None => Strata::none(),
    };
    let class_names = ["all".to_string(), format!("q{}", args.min_mapq)];
    let params = Params {
        min_depth: args.min_depth,
        het_af: args.het_af,
        window: args.window.max(1),
    };
    let prefix = &args.output_prefix;
    let mut pileup = Pileup::new(
        args.chrom.clone(),
        ref_seq,
        params,
        class_names.clone(),
        strata,
        create(format!("{}.sites.bed", prefix))?,
        create(format!("{}.windows.tsv", prefix))?,
    )?;

    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(&args.bam)?;
    let header = reader.read_header()?;
    let mut records: u64 = 0;
    if header.reference_sequences().contains_key(args.chrom.as_bytes()) {
        let region: Region = args.chrom.parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("{}: {}", args.chrom, e))
        })?;
        let mut ops: Vec<(Kind, usize)> = Vec::new();
        let mut seq: Vec<u8> = Vec::new();
        for result in reader.query(&header, &region)?.records() {
            let record = result?;
            let flags = record.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_qc_fail()
                || flags.is_duplicate()
            {
                continue;
            }
            let Some(start) = record.alignment_start() else {
                continue;
            };
            let start = usize::from(start?) as u64 - 1;
            // A missing MAPQ (255) does not count as confident.
            let confident = record
                .mapping_quality()
                .is_some_and(|q| u8::from(q) >= args.min_mapq);
            ops.clear();
            for op in record.cigar().iter() {
                let op = op?;
                ops.push((op.kind(), op.len()));
            }
            seq.clear();
            seq.extend(record.sequence().iter());
            pileup.add(start, &ops, &seq, [true, confident])?;
            records += 1;
        }
    } else {
        log::warn!("{} not in BAM header; reporting zero depth", args.chrom);
    }

    let (stats, strata, _, _) = pileup.finish()?;
    let mut hist = create(format!("{}.hist.tsv", prefix))?;
    write_hist(&mut hist, &args.chrom, &class_names, &strata, &stats, args.min_depth)?;
    hist.flush()?;
    let mut totals = create(format!("{}.totals.tsv", prefix))?;
    write_totals(&mut totals, &args.chrom, &class_names, &strata, &stats)?;
    totals.flush()?;
    log::info!("{}: {} alignments summarised", args.chrom, records);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pileup(reference: &[u8], min_depth: u32, strata: Strata) -> Pileup<Vec<u8>> {
        Pileup::new(
            "chrT".to_string(),
            reference.to_vec(),
            Params { min_depth, het_af: 0.2, window: 4 },
            ["all".to_string(), "q20".to_string()],
            strata,
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn counts_alleles_deletions_and_insertions() {
        //               0123456789
        let reference = b"ACGTACGTAC";
        let mut p = pileup(reference, 2, Strata::none());
        // Exact match over [0, 6).
        p.add(0, &[(Kind::SequenceMatch, 6)], b"ACGTAC", [true, true]).unwrap();
        // Mismatch at 2 (G->T), deletion of 3, insertion after 4; low MAPQ.
        p.add(
            1,
            &[
                (Kind::SoftClip, 2),
                (Kind::SequenceMatch, 1),
                (Kind::SequenceMismatch, 1),
                (Kind::Deletion, 1),
                (Kind::SequenceMatch, 1),
                (Kind::Insertion, 2),
                (Kind::SequenceMatch, 1),
            ],
            b"NNCTAGGC",
            [true, false],
        )
        .unwrap();
        let (stats, _, sites, _) = p.finish().unwrap();
        let all = stats.totals(0, 0);
        assert_eq!(all.positions, 10);
        assert_eq!(all.depth_sum, 6 + 5);
        assert_eq!(all.mismatch_sum, 1);
        assert_eq!(all.del_sum, 1);
        assert_eq!(all.ins_sum, 1);
        // Depth >= 2 at positions 1..=5; minor AF 0.5 at 2 and 3.
        assert_eq!(all.called_positions, 5);
        assert_eq!(all.candidate_sites, 2);
        let q20 = stats.totals(1, 0);
        assert_eq!(q20.depth_sum, 6);
        assert_eq!(q20.candidate_sites, 0);
        let sites = String::from_utf8(sites).unwrap();
        assert!(sites.contains("chrT\t2\t3\tall\tother\t2\tG\tG\t1\tT\t1"), "{}", sites);
        assert!(sites.contains("chrT\t3\t4\tall\tother\t2\tT\tT\t1\t-\t1"), "{}", sites);
    }

    #[test]
    fn rejects_unsorted_input() {
        let mut p = pileup(b"ACGTACGTAC", 1, Strata::none());
        p.add(5, &[(Kind::SequenceMatch, 2)], b"CG", [true, true]).unwrap();
        assert!(p.add(2, &[(Kind::SequenceMatch, 2)], b"GT", [true, true]).is_err());
    }

    #[test]
    fn strata_and_windows() {
        let strata = Strata::from_intervals(
            vec!["rep".to_string()],
            vec![vec![(6, 8), (2, 4), (3, 5)]],
        );
        let mut p = pileup(b"ACGTNNGTAC", 1, strata);
        p.add(0, &[(Kind::SequenceMatch, 10)], b"ACGTNNGTAC", [true, true]).unwrap();
        let (stats, strata, _, windows) = p.finish().unwrap();
        assert_eq!(strata.name(0), "rep");
        assert_eq!(strata.name(1), "other");
        // rep = [2, 5) + [6, 8), N at 4 and 5 are skipped.
        assert_eq!(stats.totals(0, 0).positions, 4);
        assert_eq!(stats.totals(0, 1).positions, 4);
        let windows = String::from_utf8(windows).unwrap();
        let rows: Vec<&str> = windows.lines().skip(1).collect();
        assert_eq!(rows, vec![
            "chrT\t0\t4\t4\t1.000\t1.000",
            "chrT\t4\t8\t2\t1.000\t1.000",
            "chrT\t8\t10\t2\t1.000\t1.000",
        ]);
    }
}
