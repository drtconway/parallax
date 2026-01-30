use crate::{
    reads::seeds::SeedCluster,
    utils::debug::{self, DebugFile},
};

pub fn write_clusters_debug(
    clusters: &[SeedCluster],
    read_name: &str,
    chrom_name: &str,
    strand_seq: &[u8],
    strand_qual: &[u8],
    read_len: usize,
    is_reverse: bool,
) {
    if false && debug::is_enabled(DebugFile::ChainsSam) {
        for (cluster_id, cluster) in clusters.iter().enumerate() {
            // Write debug chain SAM with SA tags linking seeds
            cluster.write_chain_sam(
                read_name,
                cluster_id, // cluster index as ID
                chrom_name,
                strand_seq,
                strand_qual,
            );
        }
    }

    // Write debug clusters TSV (seeds with cluster index)
    if debug::is_enabled(DebugFile::ClustersTsv) {
        for (cluster_id, cluster) in clusters.iter().enumerate() {
            let strand = if is_reverse { "-" } else { "+" };
            for hit in cluster.chain.iter() {
                // Convert strand coordinates to forward coordinates
                let (fwd_start, fwd_end) = if is_reverse {
                    (read_len - hit.read_end(), read_len - hit.read_pos)
                } else {
                    (hit.read_pos, hit.read_end())
                };
                debug::write(
                    DebugFile::ClustersTsv,
                    &format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        read_name,
                        cluster_id,
                        fwd_start,
                        fwd_end,
                        strand_seq.len(),
                        chrom_name,
                        hit.ref_pos,
                        hit.ref_end(),
                        strand,
                        hit.match_len,
                    ),
                );
            }
        }
    }
}

pub mod rmq_dp {
    use crate::reads::seeds::{SeedCluster, SeedHit};

    const MAX_DIAGONAL_DIST: i64 = 20000; // max diagonal distance for banded chaining
    const THRESHOLD: i64 = 400; // skip heuristic threshold
    const W: f64 = 2.0; // gap penalty weight

    pub fn collect_chains(seeds: &mut [SeedHit], is_reverse: bool) -> Vec<SeedCluster> {
        // Sort by reference position
        seeds.sort_by_key(|s| s.ref_pos);

        let n = seeds.len();
        let mut f = vec![0i64; n]; // best score ending at i
        let mut pred = vec![-1i32; n]; // predecessor for traceback

        // DP with banding + skip heuristic
        for i in 0..n {
            let q_i = seeds[i].read_pos;
            let r_i = seeds[i].ref_pos;
            let l_i = seeds[i].match_len;
            f[i] = l_i as i64; // base score: just this seed

            // Look at predecessors (with optimizations)
            let mut n_skip = 0;
            let max_skip = 25; // skip heuristic parameter

            for j in (0..i).rev() {
                let q_j = seeds[j].read_pos;
                let r_j = seeds[j].ref_pos;
                let l_j = seeds[j].match_len;

                // Banding: skip if diagonal too far
                let d_i = r_i as i64 - q_i as i64;
                let d_j = r_j as i64 - q_j as i64;
                if (d_i - d_j).abs() > MAX_DIAGONAL_DIST {
                    continue;
                }

                // Colinearity check
                if q_j >= q_i || r_j >= r_i {
                    continue; // not colinear
                }

                // Compute gaps
                let gap_q = q_i - q_j;
                let gap_r = r_i - r_j;

                // Match bonus: non-overlapping portion
                let alpha = l_i.min(l_j).min(gap_q).min(gap_r) as i64;

                // Gap penalty
                let g = (gap_q as i64 - gap_r as i64).abs() as usize; // indel size
                let diag_dev = (d_i - d_j).abs() as f64;
                let beta = if g == 0 {
                    0.0
                } else {
                    0.01 * W * g as f64 + 0.5 * (g as f64).log2()
                } + 0.05 * diag_dev;

                // Score this transition
                let score = f[j] + alpha - beta as i64;

                if score > f[i] {
                    f[i] = score;
                    pred[i] = j as i32;
                    n_skip = 0;
                } else if score > f[i] - THRESHOLD {
                    n_skip += 1;
                    if n_skip > max_skip {
                        break; // Skip heuristic: stop early
                    }
                }
            }
        }

        // Extract chains by backtracking from best endpoints
        //let max_chains = config::get().seeding.max_chains_per_chrom;
        let max_chains: usize = 100;
        extract_chains(&f, &pred, seeds, is_reverse, max_chains)
    }

    const MIN_CHAIN_SCORE: i64 = 75; // minimum score to accept a chain

    fn extract_chains(
        f: &[i64],
        pred: &[i32],
        seeds: &[SeedHit],
        is_reverse: bool,
        max_chains: usize,
    ) -> Vec<SeedCluster> {
        let mut chains = Vec::new();
        let mut used = vec![false; f.len()];

        loop {
            // Check chain limit (0 means no limit)
            if max_chains > 0 && chains.len() >= max_chains {
                break;
            }

            // Find best unused endpoint
            let best = (0..f.len()).filter(|&i| !used[i]).max_by_key(|&i| f[i]);

            let Some(i) = best else { break };
            if f[i] < MIN_CHAIN_SCORE {
                break;
            }
            let score = f[i];

            // Backtrack
            let mut i = i as i32;
            let mut chain = Vec::new();
            while i >= 0 {
                chain.push(seeds[i as usize]);
                used[i as usize] = true;
                i = pred[i as usize];
            }

            chain.reverse();

            let read_start = chain.first().map(|h| h.read_pos).unwrap_or(0);
            let read_end = chain.last().map(|h| h.read_end()).unwrap_or(0);
            let read_span = read_end.saturating_sub(read_start);
            let ref_start = chain.first().map(|h| h.ref_pos).unwrap_or(0);
            let ref_end = chain.last().map(|h| h.ref_end()).unwrap_or(0);
            let ref_span = ref_end.saturating_sub(ref_start);
            let mapping_density = score as f64 / (read_span.min(ref_span) as f64 + 1.0);
            log::debug!(
                "Extracted chain of length {} (score {}, read_span {}, ref_span {}, mapping_density {:.3}, read_start {})",
                chain.len(),
                score,
                read_span,
                ref_span,
                mapping_density,
                read_start
            );
            chains.push(SeedCluster::new(chain, is_reverse, 8).unwrap());
        }

        chains
    }
}

pub mod agglomerative {
    use ordered_float::OrderedFloat;

    use crate::{
        reads::seeds::{SeedCluster, SeedHit},
        utils::GroupByTrait,
    };

    struct DiagonalGroup<'a> {
        diagonal: i64,
        seeds: &'a [SeedHit],
    }

    impl<'a> DiagonalGroup<'a> {
        fn new(diagonal: i64, seeds: &'a [SeedHit]) -> Self {
            Self { diagonal, seeds }
        }

        fn weight(&self) -> f64 {
            chain_weight(self.seeds)
        }

        fn total_length(&self) -> usize {
            self.seeds.iter().map(|s| s.match_len).sum()
        }
    }

    pub fn collect_chains(
        seeds: &mut [SeedHit],
        chrom_name: &str,
        is_reverse: bool,
    ) -> Vec<SeedCluster> {
        let mut diagonals: Vec<DiagonalGroup> = Vec::new();
        for (diagonal, group) in seeds.group_by(|seed| seed.diagonal) {
            diagonals.push(DiagonalGroup::new(diagonal, group));
        }
        diagonals.sort_by_key(|g| OrderedFloat(-g.weight()));

        if log::log_enabled!(log::Level::Debug) {
            let strand = if is_reverse { "-" } else { "+" };
            for diag_group in diagonals.iter() {
                log::debug!(
                    "Diagonal {}:{}({}): {} seeds, weight {:.3}, total length {}",
                    chrom_name,
                    diag_group.diagonal,
                    strand,
                    diag_group.seeds.len(),
                    diag_group.weight(),
                    diag_group.total_length()
                );
            }
        }

        let mut clusters: Vec<SeedCluster> = Vec::new();
        for diag_group in diagonals.iter() {
            let mut merged_into_existing = false;
            for i in 0..clusters.len() {
                let diag = clusters[i].diagonal();
                let diag_dist = (diag - diag_group.diagonal as f64).abs();
                if diag_dist > 1000.0 {
                    continue;
                }
                if let Some(merged) = merge_chains(&clusters[i].chain, diag_group.seeds) {
                    log::debug!(
                        "Merging diagonal {} into existing chain {:.0} (dist {:.1})",
                        diag_group.diagonal,
                        diag,
                        diag_dist
                    );
                    clusters[i] = SeedCluster::new(merged, is_reverse, 8).unwrap();
                    merged_into_existing = true;
                    break;
                }
            }
            if !merged_into_existing && diag_group.weight() >= 250.0 {
                clusters.push(SeedCluster::new(diag_group.seeds.to_vec(), is_reverse, 8).unwrap());
            }
            clusters.sort_by_cached_key(|cluster| OrderedFloat(-chain_weight(&cluster.chain)));
        }

        if log::log_enabled!(log::Level::Info) {
            let max_weight = clusters
                .first()
                .map(|c| chain_weight(&c.chain))
                .unwrap_or(1.0);
            for (cluster_id, cluster) in clusters.iter().enumerate() {
                log::info!(
                    "Final cluster {}: {}:{}-{}({}), {}bp, {} seeds, diagonal {:.0}, weight {:.3} ({:.5})",
                    cluster_id,
                    chrom_name,
                    cluster.ref_start(),
                    cluster.ref_end(),
                    if is_reverse { "-" } else { "+" },
                    cluster.ref_end() - cluster.ref_start(),
                    cluster.chain.len(),
                    cluster.diagonal(),
                    chain_weight(&cluster.chain),
                    chain_weight(&cluster.chain) / max_weight
                );
            }
        }

        clusters
    }

    /// Take two chains if the seeds in the right hand one are
    /// colinear and non-overlapping with those in the left hand one.
    ///
    /// Small overlaps (up to MAX_OVERLAP bases) are resolved by trimming the shorter seed.
    /// If trimming makes a seed shorter than MIN_SEED, the seed is dropped entirely.
    const MAX_OVERLAP: usize = 10;
    const MIN_SEED: usize = 10;

    fn merge_chains(lhs: &[SeedHit], rhs: &[SeedHit]) -> Option<Vec<SeedHit>> {
        let mut merged = Vec::with_capacity(lhs.len() + rhs.len());
        let mut i = 0;
        let mut j = 0;

        while i < lhs.len() && j < rhs.len() {
            let left = &lhs[i];
            let right = &rhs[j];

            // Case 1: left comes entirely before right (no overlap)
            if left.read_end() <= right.read_pos && left.ref_end() <= right.ref_pos {
                merged.push(left.clone());
                i += 1;
                continue;
            }

            // Case 2: right comes entirely before left (no overlap)
            if right.read_end() <= left.read_pos && right.ref_end() <= left.ref_pos {
                merged.push(right.clone());
                j += 1;
                continue;
            }

            // Case 3: They overlap - check if colinear
            let left_before_right = left.read_pos < right.read_pos && left.ref_pos < right.ref_pos;
            let right_before_left = right.read_pos < left.read_pos && right.ref_pos < left.ref_pos;

            if !left_before_right && !right_before_left {
                // Non-colinear: one is before in read but after in ref (or vice versa)
                log::info!(
                    "Non-colinear seeds: left read [{}-{}), ref [{}-{}), right read [{}-{}), ref [{}-{})",
                    left.read_pos,
                    left.read_end(),
                    left.ref_pos,
                    left.ref_end(),
                    right.read_pos,
                    right.read_end(),
                    right.ref_pos,
                    right.ref_end()
                );
                return None;
            }

            // Determine which seed comes first positionally
            let (first, second) = if left_before_right {
                (left, right)
            } else {
                (right, left)
            };

            // Calculate overlap in both dimensions
            let read_overlap = first.read_end().saturating_sub(second.read_pos);
            let ref_overlap = first.ref_end().saturating_sub(second.ref_pos);
            let overlap = read_overlap.max(ref_overlap);

            if overlap > MAX_OVERLAP {
                let min_diag = first.diagonal.min(second.diagonal);
                log::info!(
                    "Overlap too large ({}): first read [{}-{}), ref [{}-{}), second read [{}-{}), ref [{}-{}), uniqueness ({}, {}), diagonals ({}, {})",
                    overlap,
                    first.read_pos,
                    first.read_end(),
                    first.ref_pos,
                    first.ref_end(),
                    second.read_pos,
                    second.read_end(),
                    second.ref_pos,
                    second.ref_end(),
                    first.kmer_uniqueness,
                    second.kmer_uniqueness,
                    first.diagonal - min_diag,
                    second.diagonal - min_diag
                );
                return None;
            }

            // Trim the shorter seed to resolve the overlap
            if first.match_len <= second.match_len {
                // First seed is shorter (or equal) - trim from its end
                let new_len = first.match_len.saturating_sub(overlap);
                if new_len >= MIN_SEED {
                    let mut trimmed = first.clone();
                    trimmed.match_len = new_len;
                    merged.push(trimmed);
                }
                // else: first seed too short after trimming, drop it

                // Advance the index of first seed; second stays for next iteration
                if left_before_right {
                    i += 1;
                } else {
                    j += 1;
                }
            } else {
                // Second seed is shorter - trim from its start
                // Push first seed unchanged
                merged.push(first.clone());

                let new_len = second.match_len.saturating_sub(overlap);
                if new_len >= MIN_SEED {
                    let mut trimmed = second.clone();
                    trimmed.read_pos += overlap;
                    trimmed.ref_pos += overlap;
                    trimmed.match_len = new_len;
                    merged.push(trimmed);
                }
                // else: second seed too short after trimming, drop it

                // Both seeds have been handled
                i += 1;
                j += 1;
            }
        }

        // Drain remaining seeds from lhs
        while i < lhs.len() {
            merged.push(lhs[i].clone());
            i += 1;
        }

        // Drain remaining seeds from rhs
        while j < rhs.len() {
            merged.push(rhs[j].clone());
            j += 1;
        }

        Some(merged)
    }

    fn chain_weight(chain: &[SeedHit]) -> f64 {
        let mut w = 0.0;
        for i in 0..chain.len() {
            let seed = &chain[i];
            let l = seed.match_len as f64;

            // Positive weight for longer seeds
            w += l * l.ln();

            if i > 0 {
                let prev = &chain[i - 1];
                let gap_q = seed.read_pos.saturating_sub(prev.read_end()) as f64;
                let gap_r = seed.ref_pos.saturating_sub(prev.ref_end()) as f64;
                let gap_min = gap_q.min(gap_r);

                // Penalty for gaps between seeds
                if gap_min > 0.0 {
                    w -= 0.5 * gap_min;
                }
            }
        }
        w
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// Helper to create a SeedHit with just the fields relevant to merge_chains
        fn seed(read_pos: usize, ref_pos: usize, match_len: usize) -> SeedHit {
            SeedHit::new(0, ref_pos, read_pos, 0, 0, match_len)
        }

        #[test]
        fn test_merge_empty_chains() {
            let lhs: Vec<SeedHit> = vec![];
            let rhs: Vec<SeedHit> = vec![];
            let result = merge_chains(&lhs, &rhs);
            assert_eq!(result, Some(vec![]));
        }

        #[test]
        fn test_merge_left_empty() {
            let lhs: Vec<SeedHit> = vec![];
            let rhs = vec![seed(10, 100, 5)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 1);
            assert_eq!(merged[0].read_pos, 10);
        }

        #[test]
        fn test_merge_right_empty() {
            let lhs = vec![seed(10, 100, 5)];
            let rhs: Vec<SeedHit> = vec![];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 1);
            assert_eq!(merged[0].read_pos, 10);
        }

        #[test]
        fn test_merge_non_overlapping_lhs_before_rhs() {
            // lhs: [0-5) in read, [100-105) in ref
            // rhs: [10-15) in read, [110-115) in ref
            let lhs = vec![seed(0, 100, 5)];
            let rhs = vec![seed(10, 110, 5)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            assert_eq!(merged[0].read_pos, 0);
            assert_eq!(merged[1].read_pos, 10);
        }

        #[test]
        fn test_merge_non_overlapping_rhs_before_lhs() {
            // rhs comes before lhs in both coordinates
            let lhs = vec![seed(20, 200, 5)];
            let rhs = vec![seed(5, 100, 5)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            assert_eq!(merged[0].read_pos, 5);
            assert_eq!(merged[1].read_pos, 20);
        }

        #[test]
        fn test_merge_interleaved_chains() {
            // Seeds interleave: lhs has positions 0, 20; rhs has 10, 30
            let lhs = vec![seed(0, 100, 5), seed(20, 200, 5)];
            let rhs = vec![seed(10, 150, 5), seed(30, 250, 5)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 4);
            assert_eq!(merged[0].read_pos, 0);
            assert_eq!(merged[1].read_pos, 10);
            assert_eq!(merged[2].read_pos, 20);
            assert_eq!(merged[3].read_pos, 30);
        }

        #[test]
        fn test_merge_large_overlapping_in_read_fails() {
            // Seeds overlap by 5 bases in read coordinates (> MAX_OVERLAP)
            // lhs: read [0-10), rhs: read [5-15) - overlap at [5-10) = 5 bases
            let lhs = vec![seed(0, 100, 10)];
            let rhs = vec![seed(5, 200, 10)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_none());
        }

        #[test]
        fn test_merge_large_overlapping_in_ref_fails() {
            // Seeds overlap by 5 bases in reference coordinates (> MAX_OVERLAP)
            // lhs: ref [100-110), rhs: ref [105-115) - overlap at [105-110) = 5 bases
            let lhs = vec![seed(0, 100, 10)];
            let rhs = vec![seed(20, 105, 10)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_none());
        }

        #[test]
        fn test_merge_non_colinear_fails() {
            // lhs before rhs in read, but after in reference (non-colinear)
            let lhs = vec![seed(0, 200, 5)]; // read [0-5), ref [200-205)
            let rhs = vec![seed(10, 100, 5)]; // read [10-15), ref [100-105)
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_none());
        }

        #[test]
        fn test_merge_adjacent_seeds() {
            // Seeds are exactly adjacent (no gap, no overlap)
            let lhs = vec![seed(0, 100, 10)]; // read [0-10), ref [100-110)
            let rhs = vec![seed(10, 110, 10)]; // read [10-20), ref [110-120)
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
        }

        #[test]
        fn test_merge_multiple_seeds_each_chain() {
            // lhs: 3 seeds, rhs: 2 seeds, all colinear
            let lhs = vec![seed(0, 100, 5), seed(10, 200, 5), seed(40, 500, 5)];
            let rhs = vec![seed(20, 300, 5), seed(30, 400, 5)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 5);
            // Check ordering
            let read_positions: Vec<usize> = merged.iter().map(|s| s.read_pos).collect();
            assert_eq!(read_positions, vec![0, 10, 20, 30, 40]);
        }

        #[test]
        fn test_merge_fails_midway() {
            // First seeds are compatible, but later ones conflict
            let lhs = vec![
                seed(0, 100, 5),   // compatible
                seed(15, 200, 10), // read [15-25), ref [200-210)
            ];
            let rhs = vec![
                seed(10, 150, 5),  // compatible with first lhs
                seed(20, 180, 10), // overlaps with second lhs in read
            ];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_none());
        }

        // ========== New tests for small overlap trimming ==========

        #[test]
        fn test_merge_small_overlap_trims_shorter_seed() {
            // lhs: read [0-15), ref [100-115) - 15 bases
            // rhs: read [14-34), ref [114-134) - 20 bases
            // Overlap: 1 base in read, 1 base in ref (within MAX_OVERLAP)
            // lhs is shorter, so it should be trimmed from the end
            let lhs = vec![seed(0, 100, 15)];
            let rhs = vec![seed(14, 114, 20)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            // First seed should be trimmed to 14 bases
            assert_eq!(merged[0].read_pos, 0);
            assert_eq!(merged[0].match_len, 14);
            assert_eq!(merged[0].read_end(), 14);
            // Second seed unchanged
            assert_eq!(merged[1].read_pos, 14);
            assert_eq!(merged[1].match_len, 20);
        }

        #[test]
        fn test_merge_small_overlap_trims_second_seed_when_shorter() {
            // lhs: read [0-20), ref [100-120) - 20 bases
            // rhs: read [18-33), ref [118-133) - 15 bases
            // Overlap: 2 bases (within MAX_OVERLAP)
            // rhs is shorter, so it should be trimmed from the start
            let lhs = vec![seed(0, 100, 20)];
            let rhs = vec![seed(18, 118, 15)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            // First seed unchanged
            assert_eq!(merged[0].read_pos, 0);
            assert_eq!(merged[0].match_len, 20);
            // Second seed should be trimmed from start: 2 bases removed
            assert_eq!(merged[1].read_pos, 20);
            assert_eq!(merged[1].ref_pos, 120);
            assert_eq!(merged[1].match_len, 13);
        }

        #[test]
        fn test_merge_small_overlap_drops_seed_if_too_short_after_trim() {
            // lhs: read [0-11), ref [100-111) - 11 bases
            // rhs: read [10-25), ref [110-125) - 15 bases
            // Overlap: 1 base
            // lhs is shorter, trimming 1 base leaves 10 bases (exactly MIN_SEED, should keep)
            let lhs = vec![seed(0, 100, 11)];
            let rhs = vec![seed(10, 110, 15)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            assert_eq!(merged[0].match_len, 10); // trimmed to exactly MIN_SEED
        }

        #[test]
        fn test_merge_small_overlap_drops_seed_below_min() {
            // lhs: read [0-10), ref [100-110) - 10 bases
            // rhs: read [9-25), ref [109-125) - 16 bases
            // Overlap: 1 base
            // Trimming lhs would leave 9 bases (< MIN_SEED), so drop it
            let lhs = vec![seed(0, 100, 10)];
            let rhs = vec![seed(9, 109, 16)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 1); // lhs dropped
            assert_eq!(merged[0].read_pos, 9); // only rhs remains
        }

        #[test]
        fn test_merge_overlap_exactly_at_max() {
            // Test with overlap exactly at MAX_OVERLAP (2 bases)
            // lhs: read [0-20), ref [100-120) - 20 bases
            // rhs: read [18-38), ref [118-138) - 20 bases
            // Overlap: exactly 2 bases
            let lhs = vec![seed(0, 100, 20)];
            let rhs = vec![seed(18, 118, 20)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            // First seed trimmed (both equal length, first is trimmed)
            assert_eq!(merged[0].match_len, 18);
        }

        #[test]
        fn test_merge_overlap_just_over_max_fails() {
            // Test with overlap just over MAX_OVERLAP (3 bases)
            // lhs: read [0-20), ref [100-120) - 20 bases
            // rhs: read [17-37), ref [117-137) - 20 bases
            // Overlap: 3 bases (> MAX_OVERLAP)
            let lhs = vec![seed(0, 100, 20)];
            let rhs = vec![seed(17, 117, 20)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_none());
        }

        #[test]
        fn test_merge_rhs_before_lhs_with_small_overlap() {
            // rhs comes before lhs, with small overlap
            // rhs: read [0-15), ref [100-115) - 15 bases
            // lhs: read [14-34), ref [114-134) - 20 bases
            // Overlap: 1 base
            let lhs = vec![seed(14, 114, 20)];
            let rhs = vec![seed(0, 100, 15)];
            let result = merge_chains(&lhs, &rhs);
            assert!(result.is_some());
            let merged = result.unwrap();
            assert_eq!(merged.len(), 2);
            // rhs (first) should be trimmed
            assert_eq!(merged[0].read_pos, 0);
            assert_eq!(merged[0].match_len, 14);
            // lhs (second) unchanged
            assert_eq!(merged[1].read_pos, 14);
            assert_eq!(merged[1].match_len, 20);
        }
    }
}

pub mod layered {
    use std::collections::BTreeSet;

    use ordered_float::OrderedFloat;

    use crate::reads::seeds::{SeedCluster, SeedHit};

    #[derive(Debug, PartialEq, Eq)]
    enum Cluster {
        Singleton(SeedHit),
        Chain(Vec<SeedHit>),
    }

    impl Cluster {
        pub fn weight(&self) -> f64 {
            match self {
                Cluster::Singleton(seed) => {
                    let w = seed.match_len as f64;
                    w * w.ln() / (seed.kmer_uniqueness as f64).sqrt()
                }
                Cluster::Chain(seeds) => {
                    let mut w = 0.0;
                    for i in 0..seeds.len() {
                        let seed = &seeds[i];
                        let l = seed.match_len as f64;
                        w += l * l.ln() / (seed.kmer_uniqueness as f64).sqrt();
                        if i > 0 {
                            let prev = &seeds[i - 1];
                            let gap_q = seed.read_pos.saturating_sub(prev.read_end()) as f64;
                            let gap_r = seed.ref_pos.saturating_sub(prev.ref_end()) as f64;
                            let gap_min = gap_q.min(gap_r);
                            if gap_min > 0.0 {
                                w -= 0.5 * gap_min;
                            }
                        }
                    }
                    w
                }
            }
        }
    }

    impl PartialOrd for Cluster {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for Cluster {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            let self_weight = self.weight();
            let other_weight = other.weight();
            OrderedFloat(self_weight).cmp(&OrderedFloat(other_weight))
        }
    }

    pub fn collect_chains(
        seeds: &mut [SeedHit],
        chrom_name: &str,
        is_reverse: bool,
    ) -> Vec<SeedCluster> {
        let mut clusters: BTreeSet<Cluster> = seeds
            .iter()
            .map(|s| Cluster::Singleton(s.clone()))
            .collect();

        for cluster in clusters.iter() {
            log::info!(
                "Initial cluster: {:?}, weight {:.3}",
                cluster,
                cluster.weight()
            );
        }

        vec![]
    }
}
