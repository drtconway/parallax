#![allow(dead_code)]
use crate::reads::seeds::{SeedCluster, SeedHit};

const MAX_DIAGONAL_DIST: i64 = 2000; // max diagonal distance for banded chaining
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
