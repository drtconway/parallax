use std::collections::HashMap;

use crate::error::Result;
use crate::index::Index;
use crate::kmers::Kmer;
use crate::utils::{GroupByTrait, Selection};

pub fn process_reads<const K: usize, const S: usize>(
    index: &Index<K, S>,
    fastq: &str,
) -> Result<()> {
    log::info!("Processing reads from {}", fastq);

    let reader = std::fs::File::open(fastq).map(std::io::BufReader::new)?;
    let mut reader = noodles::fastq::io::Reader::new(reader);

    let mut len_n = 0;
    let mut len_s = 0.0;
    let mut len_s2 = 0.0;

    let mut d_hist: HashMap<u64, usize> = HashMap::new();
    let mut hist: HashMap<usize, usize> = HashMap::new();
    for record in reader.records() {
        let record = record?;
        let seq = record.sequence().as_ref();

        let mut fwd_kmers: Vec<Kmer<K>> = Vec::new();
        let mut fwd_hits: Vec<(usize, i64, usize, usize)> = Vec::new();
        let mut rev_hits: Vec<(usize, i64, usize, usize)> = Vec::new();
        let mut hit_vec: Vec<(usize, usize)> = Vec::new();
        let mut j = 0;
        for (pos, selection) in Kmer::<K>::open_syncmer_iter(seq, [(); S]) {
            j += 1;
            let fwd: Option<Kmer<K>> = match &selection {
                Selection::Left(kmer) => Some(*kmer),
                Selection::Both(kmer, _) => Some(*kmer),
                _ => None,
            };
            if let Some(kmer) = fwd {
                fwd_kmers.push(kmer);
                hit_vec.clear();
                index.with(&kmer, |chrom_id, chrom_pos| {
                    hit_vec.push((chrom_id, chrom_pos));
                });
                if hit_vec.len() == 1 {
                    let (chrom_id, chrom_pos) = hit_vec[0];
                    let d = chrom_pos as i64 - pos as i64;
                    fwd_hits.push((chrom_id, d, chrom_pos, pos));
                }
            }

            let rev: Option<Kmer<K>> = match &selection {
                Selection::Right(kmer) => Some(*kmer),
                Selection::Both(_, kmer) => Some(*kmer),
                _ => None,
            };
            if let Some(kmer) = rev {
                hit_vec.clear();
                index.with(&kmer, |chrom_id, chrom_pos| {
                    hit_vec.push((chrom_id, chrom_pos));
                });
                if hit_vec.len() == 1 {
                    let (chrom_id, chrom_pos) = hit_vec[0];
                    let d = chrom_pos as i64 - pos as i64;
                    rev_hits.push((chrom_id, d, chrom_pos, pos));
                }
            }
        }

        fwd_kmers.sort_unstable();
        fwd_kmers.dedup();
        let n_fwd_kmers = fwd_kmers.len() as f64;

        len_n += 1;
        len_s += n_fwd_kmers;
        len_s2 += n_fwd_kmers * n_fwd_kmers;

        for i in 1..fwd_kmers.len() {
            assert!(fwd_kmers[i].0 > fwd_kmers[i - 1].0);
            let d = (fwd_kmers[i].0 - fwd_kmers[i - 1].0).ilog2() as u64;
            *d_hist.entry(d).or_insert(0) += 1;
        }
        fwd_kmers.clear();

        fwd_hits.sort_unstable();
        rev_hits.sort_unstable();

        for (_chrom_id, hit_list) in fwd_hits
            .as_slice()
            .group_by(|(chrom_id, _, _, _)| *chrom_id)
        {
            *hist.entry(hit_list.len()).or_insert(0) += 1;
        }
        for (_chrom_id, hit_list) in rev_hits
            .as_slice()
            .group_by(|(chrom_id, _, _, _)| *chrom_id)
        {
            *hist.entry(hit_list.len()).or_insert(0) += 1;
        }
    }

    if true {
        let len_mean = (len_s as f64) / (len_n as f64);
        let len_var = (len_s2 as f64) / (len_n as f64) - len_mean * len_mean;
        let len_sd = len_var.sqrt();
        log::info!(
            "Read length stats: n={}, mean={:.2}, sd={:.2}",
            len_n,
            len_mean,
            len_sd
        );
    }

    if false {
        let mut d_hist_vec: Vec<(u64, usize)> =
            d_hist.iter().map(|(&delta, &freq)| (delta, freq)).collect();
        d_hist_vec.sort_by_key(|(delta, _freq)| *delta);
        log::info!("Read delta histogram:");
        for (delta, freq) in d_hist_vec {
            log::info!("  {}\t{}", delta, freq);
        }
    }
    if false {
        let mut n = 0;
        let mut s = 0.0;
        let mut s2 = 0.0;
        for (&delta, &freq) in &d_hist {
            n += freq;
            s += delta as f64 * freq as f64;
            s2 += (delta as f64) * (delta as f64) * freq as f64;
        }
        let mean = s / (n as f64);
        let var = s2 / (n as f64) - mean * mean;
        log::info!(
            "Read delta stats: n={}, mean={:.2}, s2/n = {}, var={:.2}",
            n,
            mean,
            s2 / (n as f64),
            var
        );
        let stddev = var.sqrt();
        log::info!(
            "Read hit stats: n={}, mean={:.2}, stddev={:.2}",
            n,
            mean,
            stddev
        );
    }
    if false {
        let mut hist_vec: Vec<(usize, usize)> = hist.into_iter().collect();
        hist_vec.sort_by_key(|(count, _freq)| *count);
        log::info!("Read hit histogram:");
        for (count, freq) in hist_vec {
            log::info!("  {}\t{}", count, freq);
        }
    }

    Ok(())
}
