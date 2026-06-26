use std::{collections::HashMap, path::Path};

use parallax::{
    index::{self, IndexBuilder, IndexHit, load_index},
    reference::InMemoryReference,
};

pub fn analyse(fasta: &Path, path: Option<&Path>) -> std::io::Result<()> {
    let reference =
        InMemoryReference::load(fasta, true).map_err(|err| std::io::Error::other(err))?;

    const K: usize = 20;

    println!("k\ts\tgap\tcount");
    analyze_inner::<K, 8>(&reference)?;
    analyze_inner::<K, 9>(&reference)?;
    analyze_inner::<K, 10>(&reference)?;
    analyze_inner::<K, 11>(&reference)?;
    analyze_inner::<K, 12>(&reference)?;
    analyze_inner::<K, 13>(&reference)?;
    analyze_inner::<K, 14>(&reference)?;
    analyze_inner::<K, 15>(&reference)?;

    if true {
        return Ok(());
    }

    if let Some(path) = path {
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
    }

    Ok(())
}

fn analyze_inner<const K: usize, const S: usize>(
    reference: &InMemoryReference,
) -> std::io::Result<()> {
    let builder = index::asymmetric_index::AsymmetricIndexBuilder::<K, S>::make(&reference);

    let mut gaps: HashMap<u32, u32> = HashMap::new();
    let mut prev_chrom_id = u32::MAX;
    let mut prev_pos = 0u32;
    builder.kmers(&mut |_x, chrom_id, pos| {
        if chrom_id != prev_chrom_id {
            prev_chrom_id = chrom_id;
            prev_pos = pos;
            log::info!(
                "Processing chromosome {}",
                reference.chrom_info(chrom_id as usize).name
            );
        } else {
            let gap = pos - prev_pos;
            *gaps.entry(gap).or_default() += 1;
            prev_pos = pos;
        }
    });

    for (gap, count) in gaps.iter() {
        println!("{}\t{}\t{}\t{}", K, S, gap, count);
    }

    Ok(())
}
