use std::{collections::HashMap, path::Path};

use parallax::index::{IndexHit, load_index};

pub fn analyse(_fasta: &Path, path: &Path) -> std::io::Result<()> {
    //let reference =
    //    InMemoryReference::load(&fasta, true).map_err(|err| std::io::Error::other(err))?;
    log::info!("loading index {}", path.display());
    let index = load_index(path)?;
    log::info!("done.");

    let n = index.all_chrom_info().len();

    let mut hist: HashMap<usize, usize> = HashMap::new();

    let mut seeds: Vec<Vec<i32>> = vec![Vec::new(); n];

    let mut i: usize = 0;
    for hit in index.iter() {
        i += 1;
        if i & 0xfffff == 0 {
            log::info!("processed {} k-mers", i);
        }
        let IndexHit {
            query_pos: _,
            seed_kmer: _,
            loci,
            k: _,
            unpack_locus,
        } = hit;
        let n = loci.len();
        *hist.entry(n).or_default() += 1;
        for &locus in loci {
            let (chrom_idx, chrom_pos, strand) = unpack_locus(locus);
            let pos = match strand {
                parallax::index::Strand::Forward => chrom_pos as i32,
                parallax::index::Strand::Reverse => -(chrom_pos as i32),
            };
            seeds[chrom_idx].push(pos);
        }
    }

    println!("freq\tcount");
    let mut frequencies: Vec<usize> = hist.keys().copied().collect();
    frequencies.sort_unstable();

    for f in frequencies.into_iter() {
        let c = hist.get(&f).copied().unwrap_or(0);
        println!("{}\t{}", f, c);
    }
    println!("");

    let mut fwd_gap_hist: HashMap<u32, u32> = HashMap::new();
    let mut rev_gap_hist: HashMap<u32, u32> = HashMap::new();
    let mut any_gap_hist: HashMap<u32, u32> = HashMap::new();

    let mut fwd_seeds = Vec::new();
    let mut rev_seeds = Vec::new();
    for i in 0..n {
        log::info!("sorting seeds for {}", i);
        seeds[i].sort_by_key(|x| (x.abs(), -x.signum()));

        fwd_seeds.clear();
        rev_seeds.clear();
        for &seed_pos in seeds[i].iter() {
            if seed_pos >= 0 {
                fwd_seeds.push(seed_pos as u32);
            } else {
                rev_seeds.push((-seed_pos) as u32);
            }
        }

        let mut j_fwd = 0;
        let mut j_rev = 0;
        while j_fwd < fwd_seeds.len() && j_rev < rev_seeds.len() {
            let fwd_seed = fwd_seeds[j_fwd];
            let rev_seed = rev_seeds[j_rev];
            if fwd_seed < rev_seed {
                j_fwd += 1;
                if j_fwd < fwd_seeds.len() {
                    let next_fwd_seed = fwd_seeds[j_fwd];
                    let fwd_gap = next_fwd_seed - fwd_seed;
                    *fwd_gap_hist.entry(fwd_gap).or_default() += 1;
                    let any_gap = fwd_gap.min(rev_seed - fwd_seed);
                    *any_gap_hist.entry(any_gap).or_default() += 1;
                } else {
                    let any_gap = rev_seed - fwd_seed;
                    *any_gap_hist.entry(any_gap).or_default() += 1;
                }
                continue;
            }
            if rev_seed < fwd_seed {
                j_rev += 1;
                if j_rev < rev_seeds.len() {
                    let next_rev_seed = rev_seeds[j_rev];
                    let rev_gap = next_rev_seed - rev_seed;
                    *rev_gap_hist.entry(rev_gap).or_default() += 1;
                    let any_gap = rev_gap.min(fwd_seed - rev_seed);
                    *any_gap_hist.entry(any_gap).or_default() += 1;
                } else {
                    let any_gap = fwd_seed - rev_seed;
                    *any_gap_hist.entry(any_gap).or_default() += 1;
                }
                continue;
            }
            j_fwd += 1;
            j_rev += 1;
            if j_fwd < fwd_seeds.len() {
                let next_fwd_seed = fwd_seeds[j_fwd];
                let fwd_gap = next_fwd_seed - fwd_seed;
                *fwd_gap_hist.entry(fwd_gap).or_default() += 1;
            }
            if j_rev < rev_seeds.len() {
                let next_rev_seed = rev_seeds[j_rev];
                let rev_gap = next_rev_seed - rev_seed;
                *rev_gap_hist.entry(rev_gap).or_default() += 1;
            }
            let any_gap = 0;
            *any_gap_hist.entry(any_gap).or_default() += 1;
        }

        while j_fwd < fwd_seeds.len() {
            let fwd_seed = fwd_seeds[j_fwd];
            j_fwd += 1;
            if j_fwd < fwd_seeds.len() {
                let next_fwd_seed = fwd_seeds[j_fwd];
                let fwd_gap = next_fwd_seed - fwd_seed;
                *fwd_gap_hist.entry(fwd_gap).or_default() += 1;
                *any_gap_hist.entry(fwd_gap).or_default() += 1;
            }
        }

        while j_rev < rev_seeds.len() {
            let rev_seed = rev_seeds[j_rev];
            j_rev += 1;
            if j_rev < rev_seeds.len() {
                let next_rev_seed = rev_seeds[j_rev];
                let rev_gap = next_rev_seed - rev_seed;
                *rev_gap_hist.entry(rev_gap).or_default() += 1;
                *any_gap_hist.entry(rev_gap).or_default() += 1;
            }
        }
    }

    let mut gaps: Vec<u32> = fwd_gap_hist.keys().chain(rev_gap_hist.keys()).chain(any_gap_hist.keys()).copied().collect();
    gaps.sort_unstable();
    gaps.dedup();

    println!("gap\tfwd_count\trev_count\tany_count");
    for gap in gaps.into_iter() {
        let fwd_count = fwd_gap_hist.get(&gap).copied().unwrap_or(0);
        let rev_count = rev_gap_hist.get(&gap).copied().unwrap_or(0);
        let any_count = any_gap_hist.get(&gap).copied().unwrap_or(0);
        println!("{gap}\t{fwd_count}\t{rev_count}\t{any_count}");
    }

    Ok(())
}
