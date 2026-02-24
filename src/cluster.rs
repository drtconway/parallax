//! Cluster genomic repeat element instances and produce representative sequences.
//!
//! This subcommand reads a BED file of repeat element annotations
//! (produced by `scripts/ucsc_repeats_to_bed.py`) together with a reference
//! FASTA, extracts the genomic sequence of each instance (reverse-complemented
//! for minus-strand entries), and clusters instances within each
//! `(repClass, repFamily, repName)` group by sequence similarity using open
//! syncmer MinHash LSH.
//!
//! For each cluster the longest (most complete) instance is written to the
//! output FASTA so that the resulting catalogue can be used directly by the
//! `annotate` subcommand.
//!
//! ## Algorithm
//!
//! 1. Extract open syncmer hashes (K=21, S=13, FnvHasher) from each
//!    strand-resolved instance sequence.
//! 2. Build a length-`(bands × rows)` MinHash signature per instance.
//! 3. For each element group, bucket instances by `(band, band_hash)`.
//!    Instances that collide in any band are candidate pairs.
//! 4. Verify candidate pairs with exact set-intersection Jaccard on their
//!    syncmer sets; link pairs above `--min-jaccard`.
//! 5. Form connected-component clusters via UnionFind; emit the longest
//!    member of each cluster.

#![allow(dead_code)]

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;

use crate::error::Result;
use crate::kmers::Kmer;
use crate::reference::InMemoryReference;
use crate::utils::hasher::FnvHasher;
use crate::utils::sequence::reverse_complement_into;
use crate::utils::union_find::UnionFind;

// ---------------------------------------------------------------------------
// Syncmer / MinHash parameters
// ---------------------------------------------------------------------------

/// K-mer length for open syncmers.
const K: usize = 15;
/// S-mer length for open syncmers (must be < K).
const S: usize = 11;

const DEFAULT_BANDS: usize = 20;
const DEFAULT_ROWS: usize = 5;
const DEFAULT_MIN_JACCARD: f32 = 0.25;

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

    /// Number of LSH bands
    #[arg(long, default_value_t = DEFAULT_BANDS)]
    pub bands: usize,

    /// Rows per LSH band  (total MinHash functions = bands × rows)
    #[arg(long, default_value_t = DEFAULT_ROWS)]
    pub rows: usize,

    /// Minimum Jaccard similarity to merge two instances into the same cluster
    #[arg(long, default_value_t = DEFAULT_MIN_JACCARD)]
    pub min_jaccard: f32,
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
            log::warn!("BED line {}: expected ≥6 fields, got {}, skipping", lineno + 1, fields.len());
            continue;
        }
        let chrom = fields[0].to_string();
        let start: usize = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => { log::warn!("BED line {}: bad start field, skipping", lineno + 1); continue; }
        };
        let end: usize = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => { log::warn!("BED line {}: bad end field, skipping", lineno + 1); continue; }
        };
        let name = fields[3];
        let strand = fields[5].chars().next().unwrap_or('+');

        // Name format from ucsc_repeats_to_bed.py: "repName|repClass|repFamily"
        let mut parts = name.splitn(3, '|');
        let rep_name = match parts.next() {
            Some(s) => s.to_string(),
            None => { log::warn!("BED line {}: could not parse name '{}', skipping", lineno + 1, name); continue; }
        };
        let rep_class = match parts.next() {
            Some(s) => s.to_string(),
            None => { log::warn!("BED line {}: could not parse name '{}', skipping", lineno + 1, name); continue; }
        };
        let rep_family = match parts.next() {
            Some(s) => s.to_string(),
            None => { log::warn!("BED line {}: could not parse name '{}', skipping", lineno + 1, name); continue; }
        };

        records.push(BedRecord { chrom, start, end, rep_name, rep_class, rep_family, strand });
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// Syncmer extraction
// ---------------------------------------------------------------------------

/// Collect the raw 64-bit hash of every open syncmer in `seq`.
fn syncmer_hashes(seq: &[u8]) -> Vec<u64> {
    let mut hashes = Vec::new();
    Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
        seq,
        [(); S],
        |_pos, kmer| hashes.push(kmer.0),
    );
    hashes
}

// ---------------------------------------------------------------------------
// MinHash
// ---------------------------------------------------------------------------

/// Parameterised family of h = `n` multiply-xor-shift hash functions.
///
/// Each function `f_i(x) = ((a_i * x) ^ ((a_i * x) >> 33)) * b_i`.
/// The coefficients are generated deterministically from a fixed seed so that
/// results are reproducible across runs.
struct MinHasher {
    a: Vec<u64>,
    b: Vec<u64>,
}

impl MinHasher {
    fn new(n: usize) -> Self {
        // LCG with Knuth's constants for deterministic coefficient generation.
        let mut state: u64 = 0x6c62272e07bb0142;
        let mut next = || -> u64 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state
        };
        // `a` must be odd to be bijective on u64.
        let a = (0..n).map(|_| next() | 1).collect();
        let b = (0..n).map(|_| next()).collect();
        Self { a, b }
    }

    /// Compute the MinHash signature for a sorted, deduplicated slice of syncmer hashes.
    fn signature(&self, syncmers: &[u64]) -> Vec<u64> {
        let n = self.a.len();
        let mut sig = vec![u64::MAX; n];
        for &s in syncmers {
            for i in 0..n {
                let h = self.a[i].wrapping_mul(s);
                let h = h ^ (h >> 33);
                let h = h.wrapping_mul(self.b[i]);
                if h < sig[i] {
                    sig[i] = h;
                }
            }
        }
        sig
    }
}

// ---------------------------------------------------------------------------
// LSH band hashing
// ---------------------------------------------------------------------------

/// Hash the `rows` signature values for `band` using FNV-1a byte-by-byte.
fn band_hash(sig: &[u64], band: usize, rows: usize) -> u64 {
    let start = band * rows;
    let mut h: u64 = 0xcbf2_9ce4_8422_2325; // FNV offset basis
    for &v in &sig[start..start + rows] {
        for byte in v.to_le_bytes() {
            h ^= byte as u64;
            h = h.wrapping_mul(0x0000_0100_0000_01B3); // FNV prime
        }
    }
    h
}

// ---------------------------------------------------------------------------
// Exact Jaccard on sorted, deduplicated sets
// ---------------------------------------------------------------------------

fn jaccard(a: &[u64], b: &[u64]) -> f32 {
    let mut intersection = 0usize;
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => { intersection += 1; i += 1; j += 1; }
        }
    }
    let union = a.len() + b.len() - intersection;
    if union == 0 { 0.0 } else { intersection as f32 / union as f32 }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

pub fn run(args: ClusterArgs) -> Result<()> {
    let num_hashes = args.bands * args.rows;

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

    // 3. Compute syncmer MinHash signatures in parallel
    log::info!(
        "Computing MinHash signatures (K={K}, S={S}, {} hash functions, {} threads) ...",
        num_hashes, args.threads
    );
    let hasher = Arc::new(MinHasher::new(num_hashes));

    // Work channel: main → workers (record index)
    let (tx_work, rx_work) = crossbeam::channel::bounded::<usize>(args.threads * 8);
    // Result channel: workers → main (index, sorted-deduped syncmer hashes, signature)
    let (tx_result, rx_result) =
        crossbeam::channel::unbounded::<(usize, Vec<u64>, Vec<u64>)>();

    let workers: Vec<_> = (0..args.threads)
        .map(|_| {
            let rx = rx_work.clone();
            let tx = tx_result.clone();
            let reference = Arc::clone(&reference);
            let records = Arc::clone(&records);
            let chrom_index = Arc::clone(&chrom_index);
            let hasher = Arc::clone(&hasher);
            std::thread::spawn(move || {
                let mut rc_buf = Vec::new();
                for idx in rx {
                    let rec = &records[idx];
                    let ci = match chrom_index.get(&rec.chrom) {
                        Some(&ci) => ci,
                        None => {
                            log::warn!(
                                "Chromosome '{}' not found in reference; skipping record {}",
                                rec.chrom, idx
                            );
                            continue;
                        }
                    };
                    let raw = reference.get_seq(ci, rec.start, rec.end);
                    let mut hashes = if rec.strand == '-' {
                        reverse_complement_into(raw, &mut rc_buf);
                        syncmer_hashes(&rc_buf)
                    } else {
                        syncmer_hashes(raw)
                    };
                    hashes.sort_unstable();
                    hashes.dedup();
                    let sig = hasher.signature(&hashes);
                    tx.send((idx, hashes, sig)).ok();
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
    // pre-filled empty entry and will form singletons.
    let mut per_record: Vec<(Vec<u64>, Vec<u64>)> =
        (0..n).map(|_| (Vec::new(), vec![u64::MAX; num_hashes])).collect();
    let mut sigs_done = 0usize;
    for (idx, hashes, sig) in rx_result {
        per_record[idx] = (hashes, sig);
        sigs_done += 1;
        if sigs_done % 100_000 == 0 {
            log::info!("  {sigs_done}/{n} signatures computed ...");
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
            .entry((rec.rep_class.clone(), rec.rep_family.clone(), rec.rep_name.clone()))
            .or_default()
            .push(i);
    }
    log::info!("  {} groups", groups.len());

    // 5. LSH clustering within each group
    log::info!(
        "Running LSH clustering (bands={}, rows={}, min_jaccard={:.2}) ...",
        args.bands, args.rows, args.min_jaccard
    );
    let mut uf = UnionFind::new();
    // Register every record so singletons exist in the structure.
    for i in 0..n {
        uf.find(i);
    }

    let mut total_candidates = 0usize;
    let mut total_linked = 0usize;

    // Jaccard histogram buckets: [0,0.05), [0.05,0.10), ..., [0.95,1.0], [1.0]
    // Indexed as floor(j * 20), capped at 20.
    let mut jaccard_hist = vec![0usize; 21];

    let mut sorted_groups: Vec<_> = groups.iter().collect();
    sorted_groups.sort_by_key(|(k, _)| *k);

    let mut groups_done = 0usize;
    let num_groups = sorted_groups.iter().filter(|(_, v)| v.len() >= 2).count();
    for ((rc, rf, rn), indices) in &sorted_groups {
        if indices.len() < 2 {
            continue; // singleton group, nothing to do
        }

        groups_done += 1;
        log::info!("  [{groups_done}/{num_groups}] {rc}/{rf}/{rn}: {} instances", indices.len());

        // ── Syncmer set size diagnostics ────────────────────────────────
        if log::log_enabled!(log::Level::Debug) {
            let sizes: Vec<usize> = indices.iter().map(|&i| per_record[i].0.len()).collect();
            let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
            let min = sizes.iter().copied().min().unwrap_or(0);
            let max = sizes.iter().copied().max().unwrap_or(0);
            let mut s = sizes.clone();
            s.sort_unstable();
            let median = s[s.len() / 2];
            log::debug!(
                "Group {rc}/{rf}/{rn}: {} instances, syncmer set sizes: min={min} median={median} mean={mean:.0} max={max}",
                indices.len()
            );
        }

        // Build per-band LSH buckets: band_index → (band_hash → sorted record indices).
        // Keeping bands separate lets us apply the canonical-band criterion below.
        let mut band_buckets: Vec<HashMap<u64, Vec<usize>>> =
            (0..args.bands).map(|_| HashMap::new()).collect();
        for &i in indices.iter() {
            let sig = &per_record[i].1;
            for (band, bmap) in band_buckets.iter_mut().enumerate() {
                let bh = band_hash(sig, band, args.rows);
                bmap.entry(bh).or_default().push(i);
            }
        }
        // Sort each bucket so indices are ascending; makes pair (a,b) canonical (a<b)
        // without a branch, and lets us slice bucket[ai+1..] for the inner loop.
        for bmap in &mut band_buckets {
            for bucket in bmap.values_mut() {
                bucket.sort_unstable();
            }
        }

        // ── LSH bucket size diagnostics ─────────────────────────────────
        if log::log_enabled!(log::Level::Debug) {
            let non_singleton_buckets: Vec<usize> = band_buckets
                .iter()
                .flat_map(|bmap| bmap.values())
                .filter(|b| b.len() > 1)
                .map(|b| b.len())
                .collect();
            let max_bucket = non_singleton_buckets.iter().copied().max().unwrap_or(0);
            let total_bucket_pairs: usize = non_singleton_buckets
                .iter()
                .map(|&sz| sz * (sz - 1) / 2)
                .sum();
            log::debug!(
                "  LSH: {}/{} buckets non-singleton, max_bucket_size={max_bucket}, \
                 raw candidate pairs (with duplicates)={total_bucket_pairs}",
                non_singleton_buckets.len(),
                band_buckets.iter().map(|bm| bm.len()).sum::<usize>()
            );
        }

        // Enumerate candidate pairs using the canonical-band criterion.
        //
        // A pair (a, b) is processed exactly once: in the lowest-numbered band
        // where they collide.  For a pair appearing in bands {2, 5, 11}, we
        // process it in band 2 and skip it in bands 5 and 11.
        //
        // The check "did they collide in any earlier band?" is pure arithmetic:
        // band_hash(sig_a, j) == band_hash(sig_b, j) for j in 0..band.
        // No HashSet or extra allocation is needed.
        let mut group_candidates = 0usize;
        let mut group_linked = 0usize;
        for (band, bmap) in band_buckets.iter().enumerate() {
            for bucket in bmap.values() {
                if bucket.len() < 2 {
                    continue;
                }
                let mut skip_count = 0usize;
                // bucket is sorted ascending, so a < b always holds.
                for (ai, &a) in bucket.iter().enumerate() {
                    let sig_a = &per_record[a].1;
                    for &b in &bucket[ai + 1..] {
                        let sig_b = &per_record[b].1;
                        // Skip if this pair collided in any earlier band.
                        let is_canonical = (0..band).all(|j| {
                            band_hash(sig_a, j, args.rows) != band_hash(sig_b, j, args.rows)
                        });
                        if !is_canonical {
                            skip_count += 1;
                            continue;
                        }
                        group_candidates += 1;
                        if group_candidates % 100_000 == 0 {
                            log::info!("{rc}/{rf}/{rn}\t{group_candidates} pairs checked ...");
                        }
                        let j = jaccard(&per_record[a].0, &per_record[b].0);
                        let hist_bucket = ((j * 20.0) as usize).min(20);
                        jaccard_hist[hist_bucket] += 1;
                        if j >= args.min_jaccard {
                            group_linked += 1;
                            uf.union(a, b);
                        }
                        if log::log_enabled!(log::Level::Debug) {
                            let ra = &records[a];
                            let rb = &records[b];
                            log::debug!(
                                "  pair {}:{}-{} vs {}:{}-{} J={:.3}{}",
                                ra.chrom, ra.start, ra.end,
                                rb.chrom, rb.start, rb.end,
                                j,
                                if j >= args.min_jaccard { " [MERGED]" } else { "" }
                            );
                        }
                    }
                }
                log::info!(
                    "    band {band}: {}/{} pairs skipped by canonical-band criterion",
                    skip_count,
                    bucket.len() * (bucket.len() - 1) / 2
                );
            }
        }

        log::info!(
            "    {group_candidates} candidate pairs, {group_linked} merged"
        );

        total_candidates += group_candidates;
        total_linked += group_linked;
    }

    // Print Jaccard histogram
    if log::log_enabled!(log::Level::Debug) {
        log::debug!("Pairwise Jaccard distribution of candidate pairs:");
        log::debug!("  {:>12}  {:>8}  {}", "range", "count", "bar");
        for (i, &count) in jaccard_hist.iter().enumerate() {
            let lo = i as f32 * 0.05;
            let hi = if i == 20 { 1.0f32 } else { lo + 0.05 };
            let bar = "#".repeat((count as f64 / (total_candidates as f64 / 40.0)).ceil() as usize);
            log::debug!("  [{lo:.2}, {hi:.2})  {count:>8}  {bar}");
        }
    }

    log::info!(
        "  {} candidate pairs checked, {} merged (Jaccard ≥ {:.2})",
        total_candidates, total_linked, args.min_jaccard
    );

    // 6. Select best representative per cluster (longest sequence = most complete copy)
    // and compute cluster sizes in a single pass.
    log::info!("Selecting cluster representatives ...");
    let mut cluster_best: HashMap<usize, (usize, usize)> = HashMap::new(); // root → (best_idx, best_len)
    let mut cluster_size: HashMap<usize, usize> = HashMap::new();          // root → member count

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
    log::info!("  {} clusters", num_clusters);

    // 7. Write output FASTA, sorted by (repClass, repFamily, repName, chrom, start)
    // for deterministic, human-browsable output.
    let mut rep_indices: Vec<usize> = cluster_best.values().map(|&(idx, _)| idx).collect();
    rep_indices.sort_unstable_by(|&a, &b| {
        let ra = &records[a];
        let rb = &records[b];
        ra.rep_class.cmp(&rb.rep_class)
            .then(ra.rep_family.cmp(&rb.rep_family))
            .then(ra.rep_name.cmp(&rb.rep_name))
            .then(ra.chrom.cmp(&rb.chrom))
            .then(ra.start.cmp(&rb.start))
    });

    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut seq_buf: Vec<u8> = Vec::new();
    let mut written = 0usize;

    for rep_idx in rep_indices {
        let rec = &records[rep_idx];
        let ci = match chrom_index.get(&rec.chrom) {
            Some(&ci) => ci,
            None => continue,
        };

        // Extract strand-resolved sequence into seq_buf
        let raw = reference.get_seq(ci, rec.start, rec.end);
        if rec.strand == '-' {
            reverse_complement_into(raw, &mut seq_buf);
        } else {
            seq_buf.clear();
            seq_buf.extend_from_slice(raw);
        }

        let root = uf.find(rep_idx);
        let size = cluster_size.get(&root).copied().unwrap_or(1);

        // FASTA header: >repName|repClass|repFamily chrom:start-end(strand) cluster_size=N
        writeln!(
            out,
            ">{}|{}|{} {}:{}-{}({}) cluster_size={}",
            rec.rep_name, rec.rep_class, rec.rep_family,
            rec.chrom, rec.start, rec.end, rec.strand,
            size
        )?;

        // Sequence in 60-character lines
        for chunk in seq_buf.chunks(60) {
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
