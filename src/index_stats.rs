use std::{collections::HashMap, path::Path};

use parallax::index::{IndexHit, load_index};

pub fn analyse(_fasta: &Path, path: &Path) -> std::io::Result<()> {
    log::info!("loading index {}", path.display());
    let index = load_index(path)?;
    log::info!("done.");

    let n = index.all_chrom_info().len();

    let mut freq_hist: HashMap<usize, usize> = HashMap::new();
    let mut seeds: Vec<Vec<u32>> = vec![Vec::new(); n];

    let mut loci_buffer: Vec<(usize, usize)> = Vec::new();
    let mut i: usize = 0;
    for hit in index.iter() {
        i += 1;
        if i & 0xfffff == 0 {
            log::info!("processed {} k-mers", i);
        }
        let IndexHit { loci, .. } = hit;
        *freq_hist.entry(loci.len()).or_default() += 1;
        index.unpack_loci(loci, &mut loci_buffer);
        for &(chrom_idx, chrom_pos) in &loci_buffer {
            seeds[chrom_idx].push(chrom_pos as u32);
        }
    }

    println!("freq\tcount");
    let mut frequencies: Vec<usize> = freq_hist.keys().copied().collect();
    frequencies.sort_unstable();
    for f in frequencies {
        println!("{}\t{}", f, freq_hist[&f]);
    }
    println!();

    let mut gap_hist: HashMap<u32, u32> = HashMap::new();

    for i in 0..n {
        log::info!("sorting seeds for chrom {}", i);
        seeds[i].sort_unstable();
        for pair in seeds[i].windows(2) {
            let gap = pair[1] - pair[0];
            *gap_hist.entry(gap).or_default() += 1;
        }
    }

    let mut gaps: Vec<u32> = gap_hist.keys().copied().collect();
    gaps.sort_unstable();

    println!("gap\tcount");
    for gap in gaps {
        println!("{}\t{}", gap, gap_hist[&gap]);
    }

    Ok(())
}
