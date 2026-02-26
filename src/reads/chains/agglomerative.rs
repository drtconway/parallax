#![allow(dead_code)]
use ordered_float::OrderedFloat;

use crate::{
    reads::seeds::{SeedCluster, SeedHit},
    utils::{GroupByTrait, human::HumanReadable},
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
                "Final cluster {}: {}:{}-{}({}), {}bp, {} seeds, diagonal {:.0}, weight {} ({:.5})",
                cluster_id,
                chrom_name,
                cluster.ref_start(),
                cluster.ref_end(),
                if is_reverse { "-" } else { "+" },
                cluster.ref_end() - cluster.ref_start(),
                cluster.chain.len(),
                cluster.diagonal(),
                chain_weight(&cluster.chain).human(),
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
        // Seeds overlap by 12 bases in read coordinates (> MAX_OVERLAP of 10)
        // lhs: read [0-20), rhs: read [8-28) - overlap at [8-20) = 12 bases
        let lhs = vec![seed(0, 100, 20)];
        let rhs = vec![seed(8, 200, 20)];
        let result = merge_chains(&lhs, &rhs);
        assert!(result.is_none());
    }

    #[test]
    fn test_merge_large_overlapping_in_ref_fails() {
        // Seeds overlap by 12 bases in reference coordinates (> MAX_OVERLAP of 10)
        // lhs: ref [100-120), rhs: ref [108-128) - overlap at [108-120) = 12 bases
        let lhs = vec![seed(0, 100, 20)];
        let rhs = vec![seed(40, 108, 20)];
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
        // Test with overlap exactly at MAX_OVERLAP (10 bases)
        // lhs: read [0-20), ref [100-120) - 20 bases
        // rhs: read [10-30), ref [110-130) - 20 bases
        // Overlap: exactly 10 bases
        let lhs = vec![seed(0, 100, 20)];
        let rhs = vec![seed(10, 110, 20)];
        let result = merge_chains(&lhs, &rhs);
        assert!(result.is_some());
        let merged = result.unwrap();
        assert_eq!(merged.len(), 2);
        // First seed trimmed (both equal length, first is trimmed)
        assert_eq!(merged[0].match_len, 10);
    }

    #[test]
    fn test_merge_overlap_just_over_max_fails() {
        // Test with overlap just over MAX_OVERLAP (11 bases > MAX_OVERLAP of 10)
        // lhs: read [0-20), ref [100-120) - 20 bases
        // rhs: read [9-29), ref [109-129) - 20 bases
        // Overlap: 11 bases (> MAX_OVERLAP)
        let lhs = vec![seed(0, 100, 20)];
        let rhs = vec![seed(9, 109, 20)];
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
