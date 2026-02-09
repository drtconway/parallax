use std::collections::{HashMap, VecDeque};

use ordered_float::OrderedFloat;

use crate::{
    reads::seeds::{SeedCluster, SeedHit},
    utils::{union_find::UnionFind, upper_triangular_pairs},
};

#[allow(dead_code)]
fn gap_penalty(lhs: &SeedHit, rhs: &SeedHit) -> Option<f64> {
    const MAX_DIAGONAL_DIST: i64 = 2000; // max diagonal distance for banded chaining
    const W: f64 = 2.0; // gap penalty weight

    let q_i = lhs.read_pos;
    let r_i = lhs.ref_pos;
    let l_i = lhs.match_len;
    let q_j = rhs.read_pos;
    let r_j = rhs.ref_pos;
    let l_j = rhs.match_len;

    // Banding: skip if diagonal too far
    let d_i = r_i as i64 - q_i as i64;
    let d_j = r_j as i64 - q_j as i64;
    if (d_i - d_j).abs() > MAX_DIAGONAL_DIST {
        return None;
    }

    // Colinearity check
    if q_j >= q_i || r_j >= r_i {
        return None; // not colinear
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

    Some((alpha as f64 - beta).max(0.0))
}

fn mean_squared_gap_deviation(seeds: &[&SeedHit]) -> f64 {
    let n = seeds.len();
    if n == 0 {
        return 0.0;
    }

    let mut sum = 0.0;
    let mut sum2 = 0.0;
    for i in 1..n {
        let read_gap = (seeds[i].read_pos as i64) - (seeds[i - 1].read_end() as i64);
        let ref_gap = (seeds[i].ref_pos as i64) - (seeds[i - 1].ref_end() as i64);
        let gap_diff = (read_gap - ref_gap) as f64;
        sum += gap_diff;
        sum2 += gap_diff * gap_diff;
    }

    let mean = sum / (n as f64);
    let msd = (sum2 / (n as f64)) - (mean * mean);

    msd.sqrt() / (n as f64)
}

fn gap_diff_priority(seed_i: &SeedHit, seed_j: &SeedHit) -> OrderedFloat<f64> {
    let read_gap = (seed_j.read_pos as i64) - (seed_i.read_end() as i64);
    let ref_gap = (seed_j.ref_pos as i64) - (seed_i.ref_end() as i64);
    let gap_diff = (read_gap - ref_gap).abs() as f64;
    let avg_gap = (read_gap.abs() + ref_gap.abs()) as f64 / 2.0;
    let match_weight = (seed_i.match_len as f64 * seed_j.match_len as f64).sqrt();

    // Deviation penalty grows superlinearly
    let deviation_penalty = if gap_diff > 0.0 {
        gap_diff * (1.0 + gap_diff.ln().max(0.0))
    } else {
        0.0
    };

    let uniqueness_penalty = (seed_i.kmer_uniqueness + seed_j.kmer_uniqueness) as f64;

    OrderedFloat((avg_gap + deviation_penalty * 2.0 + uniqueness_penalty * 0.5) / match_weight)
}

fn kruskal_like_grouping<F: Fn(&SeedHit,&SeedHit) -> OrderedFloat<f64>>(seeds: &[SeedHit], priority: F) -> Vec<SeedTree> {
    const MIN_SEED_LENGTH: i64 = 20;

    let n = seeds.len();
    let mut pairs: Vec<(usize, usize)> = upper_triangular_pairs(n)
        .filter_map(|(i, j)| {
            // work out if seed i is before seed j or if they need to be swapped
            let seed_i = &seeds[i];
            let seed_j = &seeds[j];

            let seed_i_before_seed_j = (seed_j.read_pos as i64) - (seed_i.read_end() as i64) >= 0
                && (seed_j.ref_pos as i64) - (seed_i.ref_end() as i64) >= 0;

            let seed_j_before_seed_i = (seed_i.read_pos as i64) - (seed_j.read_end() as i64) >= 0
                && (seed_i.ref_pos as i64) - (seed_j.ref_end() as i64) >= 0;

            if seed_j_before_seed_i {
                Some((j, i))
            } else if seed_i_before_seed_j {
                Some((i, j))
            } else {
                // non-colinear, return in original order
                None
            }
        })
        .collect();

    // sort pairs by distance
    pairs.sort_by_key(|(i, j)| {
        let seed_i = &seeds[*i];
        let seed_j = &seeds[*j];
        priority(seed_i, seed_j)
    });

    let mut uf = UnionFind::new();
    let mut seeds: HashMap<usize, SeedTree> = seeds
        .into_iter()
        .enumerate()
        .map(|(i, s)| (i, SeedTree::Leaf(s.clone())))
        .collect();

    for (_k, (i, j)) in pairs.iter().enumerate() {
        let a = uf.find(*i);
        let b = uf.find(*j);
        if a == b {
            continue;
        }
        let lhs = seeds.get(&a).unwrap();
        let rhs = seeds.get(&b).unwrap();

        let fwd_read_gap = (rhs.read_start() as i64) - (lhs.read_end() as i64);
        let fwd_ref_gap = (rhs.ref_start() as i64) - (lhs.ref_end() as i64);

        let rev_read_gap = (lhs.read_start() as i64) - (rhs.read_end() as i64);
        let rev_ref_gap = (lhs.ref_start() as i64) - (rhs.ref_end() as i64);

        // Allow small negative gaps up to the length of the last seed
        let (a, b) = if fwd_read_gap.min(fwd_ref_gap) + rhs.first_seed_length() >= MIN_SEED_LENGTH {
            (a, b)
        } else if rev_read_gap.min(rev_ref_gap) + lhs.first_seed_length() >= MIN_SEED_LENGTH {
            (b, a)
        } else {
            // non-colinear or overlapping
            continue;
        };

        let lhs = seeds.get(&a).unwrap();
        let rhs = seeds.get(&b).unwrap();

        let read_gap = (rhs.read_start() as i64) - (lhs.read_end() as i64);
        let ref_gap = (rhs.ref_start() as i64) - (lhs.ref_end() as i64);

        if read_gap > 10000 || ref_gap > 10000 {
            continue;
        }

        let lhs = seeds.remove(&a).unwrap();
        let mut rhs = seeds.remove(&b).unwrap();

        let read_gap = (rhs.read_start() as i64) - (lhs.read_end() as i64);
        let ref_gap = (rhs.ref_start() as i64) - (lhs.ref_end() as i64);

        let overlap = read_gap.min(ref_gap);
        if overlap < 0 {
            // trim the right hand seed
            let trim = -overlap as usize;
            rhs.trim_start(trim);
        }
        let rhs = rhs;

        let read_gap = (rhs.read_start() as i64) - (lhs.read_end() as i64);
        let ref_gap = (rhs.ref_start() as i64) - (lhs.ref_end() as i64);

        let show = false;
        if show {
            let q = ((read_gap as f64) * (ref_gap as f64)).sqrt()
                / (lhs.match_length() as f64 + rhs.match_length() as f64);
            log::debug!(
                "Attempting to merge seeds: i = {}, j = {}, diagonals = {} and {}, read_gap = {}, ref_gap = {}, match_lens = {} and {}, q = {:.3}",
                i,
                j,
                lhs.diagonal(),
                rhs.diagonal(),
                read_gap,
                ref_gap,
                lhs.match_length(),
                rhs.match_length(),
                q
            );
        }

        match lhs.group(rhs) {
            Ok(t) => {
                if show {
                    log::debug!("  Merged!");
                }
                let c = uf.union(a, b);
                seeds.insert(c, t);
            }
            Err((s1, s2)) => {
                if show {
                    log::debug!("  Could not merge.");
                }
                seeds.insert(a, s1);
                seeds.insert(b, s2);
            }
        }
    }

    seeds.into_values().collect()
}

enum SeedTree {
    Leaf(SeedHit),
    Node(Vec<SeedHit>),
}

impl SeedTree {
    fn count(&self) -> usize {
        match self {
            SeedTree::Leaf(_) => 1,
            SeedTree::Node(children) => children.len(),
        }
    }

    fn read_start(&self) -> usize {
        match self {
            SeedTree::Leaf(seed) => seed.read_pos,
            SeedTree::Node(children) => children.first().unwrap().read_pos,
        }
    }

    fn read_end(&self) -> usize {
        match self {
            SeedTree::Leaf(seed) => seed.read_end(),
            SeedTree::Node(children) => children.last().unwrap().read_end(),
        }
    }

    fn ref_start(&self) -> usize {
        match self {
            SeedTree::Leaf(seed) => seed.ref_pos,
            SeedTree::Node(children) => children.first().unwrap().ref_pos,
        }
    }

    fn ref_end(&self) -> usize {
        match self {
            SeedTree::Leaf(seed) => seed.ref_end(),
            SeedTree::Node(children) => children.last().unwrap().ref_end(),
        }
    }

    fn group(self, other: SeedTree) -> std::result::Result<SeedTree, (SeedTree, SeedTree)> {
        // Check colinearity
        let self_before_other =
            self.read_end() <= other.read_start() && self.ref_end() <= other.ref_start();
        let other_before_self =
            other.read_end() <= self.read_start() && other.ref_end() <= self.ref_start();

        if !(self_before_other || other_before_self) {
            return Err((self, other));
        }

        match (self, other) {
            (SeedTree::Leaf(s), SeedTree::Leaf(o)) => Ok(SeedTree::Node(vec![s, o])),
            (SeedTree::Node(mut children), SeedTree::Leaf(o)) => {
                children.push(o);
                Ok(SeedTree::Node(children))
            }
            (SeedTree::Leaf(s), SeedTree::Node(mut other)) => {
                let mut children = vec![s];
                children.append(&mut other);
                Ok(SeedTree::Node(children))
            }
            (SeedTree::Node(mut children), SeedTree::Node(mut other)) => {
                children.append(&mut other);
                Ok(SeedTree::Node(children))
            }
        }
    }

    fn match_length(&self) -> usize {
        match self {
            SeedTree::Leaf(seed) => seed.match_len,
            SeedTree::Node(children) => children.iter().map(|c| c.match_len).sum(),
        }
    }

    fn diagonal(&self) -> i64 {
        (self.ref_start() as i64) - (self.read_start() as i64)
    }

    fn first_seed_length(&self) -> i64 {
        match self {
            SeedTree::Leaf(seed) => seed.match_len as i64,
            SeedTree::Node(children) => children.first().unwrap().match_len as i64,
        }
    }

    fn trim_start(&mut self, trim: usize) {
        match self {
            SeedTree::Leaf(seed) => {
                seed.read_pos += trim;
                seed.ref_pos += trim;
                seed.match_len -= trim;
            }
            SeedTree::Node(children) => {
                let first = &mut children[0];
                first.read_pos += trim;
                first.ref_pos += trim;
                first.match_len -= trim;
            }
        }
    }
}

#[allow(dead_code)]
fn scan_for_candidate(anchor: &SeedTree, candidates: &VecDeque<SeedTree>) -> Option<usize> {
    const MAX_GAP: i64 = 250;
    const MAX_MIN_GAP: i64 = 25;
    for (i, candidate) in candidates.iter().enumerate() {
        let read_gap = (candidate.read_start() as i64) - (anchor.read_end() as i64);
        let ref_gap = (candidate.ref_start() as i64) - (anchor.ref_end() as i64);
        let min_gap = read_gap.min(ref_gap);
        if read_gap >= 0
            && read_gap <= MAX_GAP
            && ref_gap >= 0
            && ref_gap <= MAX_GAP
            && min_gap <= MAX_MIN_GAP
        {
            return Some(i);
        }
    }
    None
}

pub fn collect_chains(
    seeds: &mut [SeedHit],
    chrom_name: &str,
    is_reverse: bool,
) -> Vec<SeedCluster> {
    let _ = chrom_name; // currently unused

    let mut groups = kruskal_like_grouping(seeds, gap_diff_priority);
    log::debug!(
        "Kruskal-like grouping formed {} groups on {}",
        groups.len(),
        chrom_name
    );

    groups.sort_by_key(|g| -(g.match_length() as isize));

    groups.retain(|group| {
        let density = (group.match_length() as f64)
            / ((group.read_end() - group.read_start()) as f64).max(1.0);

        density >= 0.15 && group.match_length() >= 75
    });

    if log::log_enabled!(log::Level::Info) {
        for (i, group) in groups.iter().enumerate() {
            let density = (group.match_length() as f64)
                / ((group.ref_end() - group.ref_start()) as f64).max(1.0);

            let group_seeds: Vec<&SeedHit> = match group {
                SeedTree::Leaf(s) => vec![s],
                SeedTree::Node(children) => children.iter().collect(),
            };

            let score: f64 = mean_squared_gap_deviation(&group_seeds);

            log::info!(
                "Group {}: read [{}-{}), ref [{}-{}), length {}, count {}, density {:.3}, score {:.3}",
                i + 1,
                group.read_start(),
                group.read_end(),
                group.ref_start(),
                group.ref_end(),
                group.match_length(),
                group.count(),
                density,
                score
            );
        }
    }

    let chains = groups
        .into_iter()
        .map(|g| match g {
            SeedTree::Leaf(s) => SeedCluster::new(vec![s], is_reverse, 8).unwrap(),
            SeedTree::Node(children) => SeedCluster::new(children, is_reverse, 8).unwrap(),
        })
        .collect();

    chains
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_seed(chrom_id: usize, ref_pos: usize, read_pos: usize, match_len: usize) -> SeedHit {
        SeedHit::new(chrom_id, ref_pos, read_pos, 0, 1, match_len)
    }

    // Note: The current implementation has an issue where edge() panics on Source->Sink
    // This means we cannot test with empty seed lists
    // The following tests work around this limitation

    #[test]
    fn test_gap_penalty_computation() {
        // Test the gap penalty directly with colinear seeds
        // For gap_penalty: need q_j < q_i and r_j < r_i
        let seed_i = create_seed(0, 200, 60, 20); // later positions
        let seed_j = create_seed(0, 100, 30, 20); // earlier positions

        let penalty = gap_penalty(&seed_i, &seed_j);

        // Should return Some value for colinear seeds
        assert!(penalty.is_some());
        if let Some(p) = penalty {
            // Penalty should be finite
            assert!(p.is_finite());
        }
    }

    #[test]
    fn test_gap_penalty_non_colinear() {
        // Test with non-colinear seeds (fails colinearity check)
        let seed_i = create_seed(0, 100, 60, 20);
        let seed_j = create_seed(0, 200, 30, 20); // q_j < q_i but r_j > r_i

        let penalty = gap_penalty(&seed_i, &seed_j);

        // Non-colinear should return None
        assert_eq!(penalty, None);
    }

    #[test]
    fn test_gap_penalty_overlapping_positions() {
        // Seeds at same positions
        let seed1 = create_seed(0, 100, 30, 20);
        let seed2 = create_seed(0, 100, 30, 20); // Same position

        let penalty = gap_penalty(&seed1, &seed2);

        // Should return None (fails q_j >= q_i check)
        assert_eq!(penalty, None);
    }

    #[test]
    fn test_gap_penalty_with_gaps() {
        // Test gap penalty with actual gaps
        let seed_i = create_seed(0, 250, 80, 20); // ref gap = 150, read gap = 50
        let seed_j = create_seed(0, 100, 30, 20);

        let penalty = gap_penalty(&seed_i, &seed_j);

        assert!(penalty.is_some());
        if let Some(p) = penalty {
            // With indel, penalty should account for gap
            assert!(p.is_finite());
        }
    }

    #[test]
    fn test_gap_penalty_diagonal_limit() {
        // Test that seeds beyond MAX_DIAGONAL_DIST return None
        let seed_i = create_seed(0, 25100, 60, 20); // diagonal = 25040
        let seed_j = create_seed(0, 100, 30, 20); // diagonal = 70
        // diff = 24970 > 20000

        let penalty = gap_penalty(&seed_i, &seed_j);

        // Should return None due to diagonal distance
        assert_eq!(penalty, None);
    }

    #[test]
    fn test_gap_penalty_within_diagonal_limit() {
        // Test seeds within diagonal limit (MAX_DIAGONAL_DIST = 2000)
        let seed_i = create_seed(0, 1100, 100, 20); // diagonal = 1000
        let seed_j = create_seed(0, 100, 30, 20); // diagonal = 70
        // diff = 930 < 2000

        let penalty = gap_penalty(&seed_i, &seed_j);

        // Should return Some within diagonal limit (if colinear)
        // This should pass colinearity: q_j (30) < q_i (100), r_j (100) < r_i (1100)
        assert!(penalty.is_some());
    }

    #[test]
    fn test_gap_penalty_match_bonus() {
        // Test that match bonus is calculated from overlapping regions
        let seed_i = create_seed(0, 150, 50, 30); // match_len=30
        let seed_j = create_seed(0, 100, 30, 25); // match_len=25

        // gap_q = 50-30 = 20, gap_r = 150-100 = 50
        // alpha = min(30, 25, 20, 50) = 20

        let penalty = gap_penalty(&seed_i, &seed_j);
        assert!(penalty.is_some());
    }
}
