//! Farthest-first traversal for diverse representative selection.
//!
//! Iteratively selects the point most distant (lowest cosine similarity)
//! from all currently selected representatives, stopping when every
//! remaining point has cosine similarity ≥ `min_cosine` to at least one
//! representative.  This produces a compact set of representatives that
//! covers the diversity of the group.
//!
//! Complexity: O(n·k) cosine evaluations, where k is the number of
//! representatives selected (typically small).

use parallax::utils::human::CommaReadable;

use super::{cosine_similarity, GroupStats};

/// Run farthest-first traversal to select diverse representatives.
///
/// `per_member` is indexed by group-local index and contains
/// `(freq_vector, l2_norm, sequence)` for each member.
///
/// `genomic_lens` gives the genomic length of each member.  The initial
/// seed is the longest sequence (most complete instance).
///
/// Returns the list of group-local indices chosen as representatives,
/// plus accumulated statistics.
pub fn cluster_group(
    per_member: &[(Vec<u8>, f32, Vec<u8>)],
    genomic_lens: &[usize],
    group_label: &str,
    min_cosine: f32,
    max_representatives: usize,
    _collect_stats: bool,
) -> (Vec<usize>, GroupStats) {
    let group_n = per_member.len();
    let mut stats = GroupStats::default();

    // Filter to members with non-zero frequency vectors
    let active: Vec<usize> = (0..group_n)
        .filter(|&i| per_member[i].1 > 0.0)
        .collect();

    if active.is_empty() {
        return (Vec::new(), stats);
    }

    // Seed: pick the longest active member
    let seed = *active
        .iter()
        .max_by_key(|&&i| genomic_lens[i])
        .unwrap();

    let mut representatives: Vec<usize> = vec![seed];

    // dist[j] = max cosine similarity from active[j] to any selected rep.
    // We want to stop when min over all non-reps of dist[j] >= min_cosine,
    // and we want to pick argmin(dist) next (most distant = lowest cosine).
    let mut best_sim: Vec<f32> = active
        .iter()
        .map(|&i| {
            let (ref freq_s, norm_s, _) = per_member[seed];
            let (ref freq_i, norm_i, _) = per_member[i];
            cosine_similarity(freq_s, freq_i, norm_s, norm_i)
        })
        .collect();

    loop {
        // Find the active point with the lowest similarity to any rep
        // (i.e. most distant).  Skip points already selected.
        let mut worst_idx = None;
        let mut worst_sim = f32::MAX;

        for (ai, &global) in active.iter().enumerate() {
            if representatives.contains(&global) {
                continue;
            }
            if best_sim[ai] < worst_sim {
                worst_sim = best_sim[ai];
                worst_idx = Some(ai);
            }
        }

        let Some(ai) = worst_idx else {
            // All active points are already representatives
            break;
        };

        // Stop if we've reached the maximum number of representatives,
        // or if the most distant point is already well-covered.
        if representatives.len() >= max_representatives {
            break;
        }
        if worst_sim >= min_cosine {
            break;
        }

        // Add this point as a new representative
        let new_rep = active[ai];
        representatives.push(new_rep);

        // Update best_sim for all active points
        let (ref freq_r, norm_r, _) = per_member[new_rep];
        for (aj, &j) in active.iter().enumerate() {
            let (ref freq_j, norm_j, _) = per_member[j];
            let sim = cosine_similarity(freq_r, freq_j, norm_r, norm_j);
            if sim > best_sim[aj] {
                best_sim[aj] = sim;
            }
        }
    }

    log::info!(
        "    {} selected {} representatives (min_cosine={:.2})",
        group_label,
        representatives.len().commas(),
        min_cosine
    );

    stats.clusters = representatives.len();
    stats.singletons = 0;
    stats.cluster_sizes = vec![1; representatives.len()];

    (representatives, stats)
}
