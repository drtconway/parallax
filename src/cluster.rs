//! Cluster genomic repeat element instances and produce representative sequences.
//!
//! This subcommand reads a BED file of repeat element annotations
//! (produced by `scripts/ucsc_repeats_to_bed.py`) together with a reference
//! FASTA, extracts the genomic sequence of each instance (reverse-complemented
//! for minus-strand entries), and clusters instances within each
//! `(repClass, repFamily, repName)` group by sequence similarity using
//! dense k-mer frequency vectors and cosine similarity.
//!
//! For each cluster the longest (most complete) instance is written to the
//! output FASTA so that the resulting catalogue can be used directly by the
//! `annotate` subcommand.
//!
//! ## Algorithm
//!
//! 1. For each strand-resolved instance, build a dense `u8` frequency vector
//!    over the 4^K k-mer space (K=5 → 1024 entries, 1 KiB per instance).
//! 2. Precompute L2 norms for each vector.
//! 3. Within each `(repClass, repFamily, repName)` group, compute all-pairs
//!    cosine similarity and link pairs above `--min-cosine`.
//! 4. Form connected-component clusters via UnionFind; emit the longest
//!    member of each cluster.

#![allow(dead_code)]

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;

use crate::align::{Aligner, Kind, Op};
use crate::error::Result;
use crate::kmers::Kmer;
use crate::reference::InMemoryReference;
use crate::utils::sequence::reverse_complement_into;
use crate::utils::union_find_2::UnionFind;

// ---------------------------------------------------------------------------
// K-mer / MinHash parameters
// ---------------------------------------------------------------------------

/// K-mer length for frequency vectors.  Smaller K produces a denser,
/// more informative vector at the cost of k-mer specificity.  K=5
/// gives 4^5 = 1024 distinct k-mers – compact enough for u8 vectors
/// (1 KiB each) while still discriminating repeat families well.
const K: usize = 5;

/// Total number of distinct k-mers: 4^K.
const NUM_KMERS: usize = 1 << (2 * K);

const DEFAULT_MIN_COSINE: f32 = 0.75;
const DEFAULT_MIN_IDENTITY: f64 = 0.85;

// ---------------------------------------------------------------------------
// CLI args
// ---------------------------------------------------------------------------

/// Cluster repeat element instances and produce representative sequences.
#[derive(Args, Debug, Clone)]
pub struct ClusterArgs {
    /// Reference FASTA file (plain or compressed; must be readable by niffler)
    pub reference: PathBuf,

    /// Input BED file produced by scripts/ucsc_repeats_to_bed.py
    /// (columns: chrom, start, end, repName|repClass|repFamily, score, strand)
    pub bed: PathBuf,

    /// Output FASTA file (default: stdout)
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Number of threads
    #[arg(short = 't', long, default_value = "4")]
    pub threads: usize,

    /// Minimum cosine similarity to consider a pair for alignment verification
    #[arg(long, default_value_t = DEFAULT_MIN_COSINE)]
    pub min_cosine: f32,

    /// Minimum alignment identity fraction to merge two instances
    #[arg(long, default_value_t = DEFAULT_MIN_IDENTITY)]
    pub min_identity: f64,

    /// Emit diagnostic histograms (cosine, cluster-size, inter-cluster identity)
    #[arg(long)]
    pub stats: bool,
}

// ---------------------------------------------------------------------------
// BED record
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct BedRecord {
    chrom: String,
    start: usize, // 0-based, half-open
    end: usize,
    rep_name: String,
    rep_class: String,
    rep_family: String,
    strand: char,
}

impl BedRecord {
    fn genomic_len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    fn group_key(&self) -> (&str, &str, &str) {
        (&self.rep_class, &self.rep_family, &self.rep_name)
    }
}

fn parse_bed(path: &PathBuf) -> std::io::Result<Vec<BedRecord>> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.splitn(7, '\t').collect();
        if fields.len() < 6 {
            log::warn!(
                "BED line {}: expected ≥6 fields, got {}, skipping",
                lineno + 1,
                fields.len()
            );
            continue;
        }
        let chrom = fields[0].to_string();
        let start: usize = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => {
                log::warn!("BED line {}: bad start field, skipping", lineno + 1);
                continue;
            }
        };
        let end: usize = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => {
                log::warn!("BED line {}: bad end field, skipping", lineno + 1);
                continue;
            }
        };
        let name = fields[3];
        let strand = fields[5].chars().next().unwrap_or('+');

        // Name format from ucsc_repeats_to_bed.py: "repName|repClass|repFamily"
        let mut parts = name.splitn(3, '|');
        let rep_name = match parts.next() {
            Some(s) => s.to_string(),
            None => {
                log::warn!(
                    "BED line {}: could not parse name '{}', skipping",
                    lineno + 1,
                    name
                );
                continue;
            }
        };
        let rep_class = match parts.next() {
            Some(s) => s.to_string(),
            None => {
                log::warn!(
                    "BED line {}: could not parse name '{}', skipping",
                    lineno + 1,
                    name
                );
                continue;
            }
        };
        let rep_family = match parts.next() {
            Some(s) => s.to_string(),
            None => {
                log::warn!(
                    "BED line {}: could not parse name '{}', skipping",
                    lineno + 1,
                    name
                );
                continue;
            }
        };

        records.push(BedRecord {
            chrom,
            start,
            end,
            rep_name,
            rep_class,
            rep_family,
            strand,
        });
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// K-mer frequency vector
// ---------------------------------------------------------------------------

/// Build a dense frequency vector over the 4^K k-mer space.
///
/// Counts are accumulated in a caller-supplied `u16` buffer (one per
/// worker thread) and then scaled to `u8`:
///   - if max ≤ 255: direct truncation (lossless);
///   - otherwise: linear rescale so that max maps to 255.
fn kmer_freq_vector(seq: &[u8], counts: &mut [u16; NUM_KMERS]) -> Vec<u8> {
    counts.iter_mut().for_each(|c| *c = 0);
    Kmer::<K>::kmerize_fwd(seq, |_pos, kmer| {
        let idx = kmer.0 as usize;
        counts[idx] = counts[idx].saturating_add(1);
    });
    let max_count = counts.iter().copied().max().unwrap_or(0);
    if max_count == 0 {
        return vec![0u8; NUM_KMERS];
    }
    if max_count <= 255 {
        counts.iter().map(|&c| c as u8).collect()
    } else {
        let scale = max_count as f32 / 255.0;
        counts.iter().map(|&c| (c as f32 / scale) as u8).collect()
    }
}

// ---------------------------------------------------------------------------
// Cosine similarity on dense u8 frequency vectors
// ---------------------------------------------------------------------------

/// Euclidean (L2) norm of a `u8` frequency vector.
///
/// The inner loop accumulates into `u32` so that 1024 × 255² = 66 million
/// fits comfortably.  The simple sum-of-squares loop auto-vectorises well
/// on x86-64 and aarch64.
#[inline]
fn l2_norm(v: &[u8]) -> f32 {
    let sum_sq: u32 = v.iter().map(|&x| (x as u32) * (x as u32)).sum();
    (sum_sq as f32).sqrt()
}

/// Cosine similarity between two dense `u8` frequency vectors whose L2
/// norms have been precomputed.
///
/// Returns 0.0 when either vector is the zero vector.  The dot-product
/// loop is written to mirror `l2_norm` for auto-vectorisation.
#[inline]
fn cosine_similarity(a: &[u8], b: &[u8], norm_a: f32, norm_b: f32) -> f32 {
    debug_assert_eq!(a.len(), b.len());
    let denom = norm_a * norm_b;
    if denom == 0.0 {
        return 0.0;
    }
    let dot: u32 = a.iter().zip(b.iter()).map(|(&x, &y)| x as u32 * y as u32).sum();
    dot as f32 / denom
}

// ---------------------------------------------------------------------------
// Alignment identity & edit distance from extended CIGAR
// ---------------------------------------------------------------------------

/// Compute identity fraction and edit distance from an extended CIGAR.
///
/// Identity = (sequence-match columns) / (total alignment columns).
/// Edit distance = mismatches + insertions + deletions (by column count).
fn alignment_stats(cigar: &[Op]) -> (f64, usize) {
    let mut matches = 0usize;
    let mut columns = 0usize;
    let mut edits = 0usize;
    for &op in cigar {
        let n = op.len();
        match op.kind() {
            Kind::SequenceMatch => {
                matches += n;
                columns += n;
            }
            Kind::SequenceMismatch => {
                edits += n;
                columns += n;
            }
            Kind::Insertion | Kind::Deletion => {
                edits += n;
                columns += n;
            }
            Kind::Match => {
                // Ambiguous M — count as column but not match or edit
                columns += n;
            }
            _ => {}
        }
    }
    let identity = if columns == 0 {
        0.0
    } else {
        matches as f64 / columns as f64
    };
    (identity, edits)
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

pub fn run(args: ClusterArgs) -> Result<()> {
    // 1. Load reference
    log::info!("Loading reference from {} ...", args.reference.display());
    let reference = Arc::new(InMemoryReference::load(&args.reference, false)?);
    let chrom_index: Arc<HashMap<String, usize>> = Arc::new(
        (0..reference.num_chroms())
            .map(|i| (reference.chrom_name(i).to_string(), i))
            .collect(),
    );

    // 2. Parse BED
    log::info!("Parsing BED from {} ...", args.bed.display());
    let records: Arc<Vec<BedRecord>> = Arc::new(parse_bed(&args.bed)?);
    let n = records.len();
    log::info!("  {} records", n);

    // 3. Compute k-mer frequency vectors in parallel
    log::info!(
        "Computing k-mer frequency vectors (K={K}, {NUM_KMERS} k-mers, {} threads) ...",
        args.threads
    );

    // Work channel: main → workers (record index)
    let (tx_work, rx_work) = crossbeam::channel::bounded::<usize>(args.threads * 8);
    // Result channel: workers → main (index, frequency vector, L2 norm, cached sequence)
    let (tx_result, rx_result) =
        crossbeam::channel::unbounded::<(usize, Vec<u8>, f32, Vec<u8>)>();

    let workers: Vec<_> = (0..args.threads)
        .map(|_| {
            let rx = rx_work.clone();
            let tx = tx_result.clone();
            let reference = Arc::clone(&reference);
            let records = Arc::clone(&records);
            let chrom_index = Arc::clone(&chrom_index);
            std::thread::spawn(move || {
                let mut rc_buf = Vec::new();
                let mut counts = [0u16; NUM_KMERS];
                for idx in rx {
                    let rec = &records[idx];
                    let ci = match chrom_index.get(&rec.chrom) {
                        Some(&ci) => ci,
                        None => {
                            log::warn!(
                                "Chromosome '{}' not found in reference; skipping record {}",
                                rec.chrom,
                                idx
                            );
                            continue;
                        }
                    };
                    let raw = reference.get_seq(ci, rec.start, rec.end);
                    if raw.iter().any(|&b| b == b'N' || b == b'n') {
                        continue;
                    }
                    let (freq, seq) = if rec.strand == '-' {
                        reverse_complement_into(raw, &mut rc_buf);
                        (kmer_freq_vector(&rc_buf, &mut counts), rc_buf.clone())
                    } else {
                        (kmer_freq_vector(raw, &mut counts), raw.to_vec())
                    };
                    let norm = l2_norm(&freq);
                    tx.send((idx, freq, norm, seq)).ok();
                }
                // tx dropped here → workers signal completion
            })
        })
        .collect();

    // Drop the main thread's copies of the channel endpoints so that the
    // result channel closes when all workers finish.
    drop(tx_result);
    drop(rx_work);

    // Feed work items
    for i in 0..n {
        tx_work.send(i).unwrap();
    }
    drop(tx_work); // signal workers: no more work

    // Collect results; records for which the chromosome was missing retain the
    // pre-filled empty entry (zero vector, zero norm, empty seq) and will form singletons.
    let mut per_record: Vec<(Vec<u8>, f32, Vec<u8>)> = (0..n)
        .map(|_| (vec![0u8; NUM_KMERS], 0.0f32, Vec::new()))
        .collect();
    let mut vecs_done = 0usize;
    for (idx, freq, norm, seq) in rx_result {
        per_record[idx] = (freq, norm, seq);
        vecs_done += 1;
        if vecs_done % 100_000 == 0 {
            log::info!("  {vecs_done}/{n} frequency vectors computed ...");
        }
    }
    for w in workers {
        w.join().unwrap();
    }

    // 4. Group records by (repClass, repFamily, repName)
    log::info!("Grouping into element families ...");
    let mut groups: HashMap<(String, String, String), Vec<usize>> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        groups
            .entry((
                rec.rep_class.clone(),
                rec.rep_family.clone(),
                rec.rep_name.clone(),
            ))
            .or_default()
            .push(i);
    }
    log::info!("  {} groups", groups.len());

    // 5. All-pairs cosine clustering within each group
    log::info!(
        "Running cosine clustering (min_cosine={:.2}, min_identity={:.2}) ...",
        args.min_cosine,
        args.min_identity
    );
    let mut aligner = Aligner::with_defaults();
    let uf = UnionFind::new(n);

    let mut total_candidates = 0usize;
    let mut total_cosine_pass = 0usize;
    let mut total_already_merged = 0usize;
    let mut total_aligned = 0usize;
    let mut total_linked = 0usize;

    // Cosine histogram buckets: [0,0.05), [0.05,0.10), ..., [0.95,1.0], [1.0]
    // Indexed as floor(c * 20), capped at 20.
    let mut cosine_hist = if args.stats { vec![0usize; 21] } else { Vec::new() };

    let mut sorted_groups: Vec<_> = groups.iter().collect();
    sorted_groups.sort_by_key(|(k, _)| *k);

    let mut groups_done = 0usize;
    let num_groups = sorted_groups.iter().filter(|(_, v)| v.len() >= 2).count();
    for ((rc, rf, rn), indices) in &sorted_groups {
        if indices.len() < 2 {
            continue; // singleton group, nothing to do
        }

        groups_done += 1;
        log::info!(
            "  [{groups_done}/{num_groups}] {rc}/{rf}/{rn}: {} instances",
            indices.len()
        );

        // ── Norm diagnostics ────────────────────────────────────────────
        if log::log_enabled!(log::Level::Debug) {
            let norms: Vec<f32> = indices.iter().map(|&i| per_record[i].1).collect();
            let non_zero: Vec<f32> = norms.iter().copied().filter(|&n| n > 0.0).collect();
            let mean = if non_zero.is_empty() {
                0.0
            } else {
                non_zero.iter().sum::<f32>() / non_zero.len() as f32
            };
            log::debug!(
                "Group {rc}/{rf}/{rn}: {} instances, {}/{} non-zero norms, mean_norm={mean:.1}",
                indices.len(),
                non_zero.len(),
                indices.len()
            );
        }

        // ── All-pairs cosine similarity ─────────────────────────────────
        // Filter to indices with non-zero frequency vectors.
        let active: Vec<usize> = indices
            .iter()
            .copied()
            .filter(|&i| per_record[i].1 > 0.0)
            .collect();

        let num_pairs = active.len() * (active.len().saturating_sub(1)) / 2;
        if num_pairs > 1_000_000 {
            log::warn!(
                "    {rc}/{rf}/{rn}: {} active instances → {} pairs (large group!)",
                active.len(),
                num_pairs
            );
        }

        let mut percentage_done = 0usize;
        let mut group_candidates = 0usize;
        let mut group_linked = 0usize;
        for (ai, &a) in active.iter().enumerate() {
            let (ref freq_a, norm_a, _) = per_record[a];
            for &b in &active[ai + 1..] {
                let (ref freq_b, norm_b, _) = per_record[b];
                group_candidates += 1;
                let new_percentage = (group_candidates * 100) / num_pairs;
                if new_percentage >= percentage_done + 1 {
                    percentage_done = new_percentage;
                    log::info!(
                        "    {rc}/{rf}/{rn}\t{percentage_done}%"
                    );
                }
                let c = cosine_similarity(freq_a, freq_b, norm_a, norm_b);
                if args.stats {
                    let hist_bucket = ((c * 20.0) as usize).min(20);
                    cosine_hist[hist_bucket] += 1;
                }
                if c >= args.min_cosine {
                    total_cosine_pass += 1;

                    // Skip alignment if already in the same partition
                    if uf.find(a) == uf.find(b) {
                        total_already_merged += 1;
                        // Already connected, no need to align
                    } else {
                        // Alignment verification
                        total_aligned += 1;
                        let seq_a = &per_record[a].2;
                        let seq_b = &per_record[b].2;
                        if let Some(aln) = aligner.align(seq_a, seq_b) {
                            let (identity, _edit_dist) = alignment_stats(&aln.cigar);
                            if identity >= args.min_identity {
                                group_linked += 1;
                                uf.union(a, b);
                            }
                        }
                    }
                }
                if log::log_enabled!(log::Level::Debug) {
                    let ra = &records[a];
                    let rb = &records[b];
                    log::debug!(
                        "  pair {}:{}-{} vs {}:{}-{} cos={:.3}{}",
                        ra.chrom,
                        ra.start,
                        ra.end,
                        rb.chrom,
                        rb.start,
                        rb.end,
                        c,
                        if c >= args.min_cosine {
                            " [MERGED]"
                        } else {
                            ""
                        }
                    );
                }
            }
        }

        log::info!("    {group_candidates} pairs, {group_linked} merged");

        total_candidates += group_candidates;
        total_linked += group_linked;
    }

    // Print cosine histogram
    if args.stats {
        log::info!("Pairwise cosine distribution:");
        log::info!("  {:>12}  {:>8}  {}", "range", "count", "bar");
        for (i, &count) in cosine_hist.iter().enumerate() {
            let lo = i as f32 * 0.05;
            let hi = if i == 20 { 1.0f32 } else { lo + 0.05 };
            let bar_len = if total_candidates > 0 {
                ((count as f64 / total_candidates as f64) * 200.0).ceil() as usize
            } else {
                0
            };
            let bar = "#".repeat(bar_len);
            log::info!("  [{lo:.2}, {hi:.2})  {count:>8}  {bar}");
        }
    }

    log::info!(
        "  {} pairs checked, {} cosine pass, {} already merged, {} aligned, {} linked (cosine \u{2265} {:.2}, identity \u{2265} {:.2})",
        total_candidates,
        total_cosine_pass,
        total_already_merged,
        total_aligned,
        total_linked,
        args.min_cosine,
        args.min_identity
    );

    // 6. Select best representative per cluster (longest sequence = most complete copy)
    // and compute cluster sizes in a single pass.
    log::info!("Selecting cluster representatives ...");
    let mut cluster_best: HashMap<usize, (usize, usize)> = HashMap::new(); // root → (best_idx, best_len)
    let mut cluster_size: HashMap<usize, usize> = HashMap::new(); // root → member count

    for i in 0..n {
        let root = uf.find(i);
        let len = records[i].genomic_len();
        *cluster_size.entry(root).or_insert(0) += 1;
        let entry = cluster_best.entry(root).or_insert((i, 0));
        if len > entry.1 {
            *entry = (i, len);
        }
    }

    let num_clusters = cluster_best.len();
    let num_singletons = cluster_size.values().filter(|&&sz| sz == 1).count();
    log::info!("  {} clusters ({} singletons)", num_clusters, num_singletons);

    // Cluster size histogram (log2-scale buckets: 1, 2, 3-4, 5-8, 9-16, ...)
    if args.stats {
        let mut size_hist: Vec<(String, usize)> = Vec::new();
        let mut sizes: Vec<usize> = cluster_size.values().copied().collect();
        sizes.sort_unstable();
        // Build log2-scale buckets
        let mut lo = 1usize;
        loop {
            let hi = if lo <= 2 { lo } else { (lo * 2) - 1 };
            let count = sizes.iter().filter(|&&s| s >= lo && s <= hi).count();
            if count > 0 {
                let label = if lo == hi {
                    format!("{lo}")
                } else {
                    format!("{lo}-{hi}")
                };
                size_hist.push((label, count));
            }
            if hi >= *sizes.last().unwrap_or(&1) {
                break;
            }
            lo = hi + 1;
        }
        let max_count = size_hist.iter().map(|(_, c)| *c).max().unwrap_or(1);
        log::info!("Cluster size distribution:");
        log::info!("  {:>12}  {:>8}  {}", "size", "clusters", "bar");
        for (label, count) in &size_hist {
            let bar_len = (*count as f64 / max_count as f64 * 40.0).ceil() as usize;
            let bar = "#".repeat(bar_len);
            log::info!("  {:>12}  {:>8}  {}", label, count, bar);
        }
    }

    // 7. Inter-cluster identity: pairwise alignment of representative sequences
    //    within each (repClass, repFamily, repName) group.
    if args.stats {
        log::info!("Computing inter-cluster pairwise identity for representative sequences ...");
        let mut rep_groups: HashMap<(&str, &str, &str), Vec<usize>> = HashMap::new();
        for &(idx, _) in cluster_best.values() {
            // Skip records with no cached sequence
            if per_record[idx].2.is_empty() {
                continue;
            }
            let rec = &records[idx];
            rep_groups
                .entry(rec.group_key())
                .or_default()
                .push(idx);
        }

        // Identity histogram: [0,1%), [1,2%), ..., [99,100%], i.e. 101 buckets
        // but we'll use 5%-wide buckets for readability: [0,5%), [5,10%), ..., [95,100%], [100%]
        let mut ident_hist = vec![0usize; 21];
        let mut total_rep_pairs = 0usize;

        for ((_rc, _rf, _rn), reps) in &rep_groups {
            if reps.len() < 2 {
                continue;
            }
            for (ai, &a) in reps.iter().enumerate() {
                let seq_a = &per_record[a].2;
                for &b in &reps[ai + 1..] {
                    let seq_b = &per_record[b].2;
                    total_rep_pairs += 1;
                    if let Some(aln) = aligner.align(seq_a, seq_b) {
                        let (identity, _) = alignment_stats(&aln.cigar);
                        let pct = identity * 100.0;
                        let bucket = ((pct / 5.0) as usize).min(20);
                        ident_hist[bucket] += 1;
                    }
                }
            }
        }

        log::info!("Inter-cluster representative identity ({total_rep_pairs} pairs):");
        log::info!("  {:>12}  {:>8}  {}", "identity%", "count", "bar");
        let max_count = ident_hist.iter().copied().max().unwrap_or(1);
        for (i, &count) in ident_hist.iter().enumerate() {
            let lo = i as f64 * 5.0;
            let hi = if i == 20 { 100.0 } else { lo + 5.0 };
            let bar_len = if max_count > 0 {
                (count as f64 / max_count as f64 * 40.0).ceil() as usize
            } else {
                0
            };
            let bar = "#".repeat(bar_len);
            log::info!("  [{lo:5.1},{hi:5.1})  {count:>8}  {bar}");
        }
    } // end if args.stats (inter-cluster identity)

    // 8. Write output FASTA, sorted by (repClass, repFamily, repName, chrom, start)
    // for deterministic, human-browsable output.
    let mut rep_indices: Vec<usize> = cluster_best.values().map(|&(idx, _)| idx).collect();
    rep_indices.sort_unstable_by(|&a, &b| {
        let ra = &records[a];
        let rb = &records[b];
        ra.rep_class
            .cmp(&rb.rep_class)
            .then(ra.rep_family.cmp(&rb.rep_family))
            .then(ra.rep_name.cmp(&rb.rep_name))
            .then(ra.chrom.cmp(&rb.chrom))
            .then(ra.start.cmp(&rb.start))
    });

    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut written = 0usize;

    for rep_idx in rep_indices {
        let rec = &records[rep_idx];
        let seq = &per_record[rep_idx].2;

        // Skip records with no cached sequence (N-containing or missing chrom)
        if seq.is_empty() {
            continue;
        }

        let root = uf.find(rep_idx);
        let size = cluster_size.get(&root).copied().unwrap_or(1);

        // FASTA header: >repName|repClass|repFamily chrom:start-end(strand) cluster_size=N
        writeln!(
            out,
            ">{}|{}|{} {}:{}-{}({}) cluster_size={}",
            rec.rep_name,
            rec.rep_class,
            rec.rep_family,
            rec.chrom,
            rec.start,
            rec.end,
            rec.strand,
            size
        )?;

        // Sequence in 60-character lines
        for chunk in seq.chunks(60) {
            out.write_all(chunk)?;
            writeln!(out)?;
        }

        written += 1;
    }

    match &args.output {
        Some(path) => log::info!("Done. Wrote {} sequences to {}.", written, path.display()),
        None => log::info!("Done. Wrote {} sequences.", written),
    }

    Ok(())
}
