use ordered_float::OrderedFloat;
use suffix::SuffixTable;

use crate::{align::anchor::Anchor, utils::GroupByTrait};

pub fn select_lcs_anchors(query: &[u8], reference: &[u8], min_length: usize) -> Vec<Anchor> {
    let anchors = find_lcs_anchors(query, reference, min_length);
    let anchors = remove_redundant_anchors(anchors);
    let mut diagonals_with_weights: Vec<(isize, f64, &[Anchor])> = anchors
        .group_by(|a| a.diagonal())
        .map(|group| {
            let diag = group.0;
            let weight: f64 = group
                .1
                .iter()
                .map(|a| {
                    let l = a.length as f64;
                    l * l.log2()
                })
                .sum();
            (diag, weight, group.1)
        })
        .collect();
    diagonals_with_weights.sort_by_key(|&(_, w, _)| OrderedFloat(-w));
    let mut selected = Vec::new();
    for (_diagonal, _wight, anchors) in diagonals_with_weights.iter() {
        if selected.is_empty() {
            selected.extend(anchors.iter().cloned());
            selected.sort_by(Anchor::order_by_query_pos);
            continue;
        }
        let average_diagonal: f64 =
            selected.iter().map(|a| a.diagonal() as f64).sum::<f64>() / (selected.len() as f64);
        let mut extras = Vec::new();
        for anchor in *anchors {
            let i =
                selected.partition_point(|a: &Anchor| a.query_pos + a.length <= anchor.query_pos);
            let mut query_gap_start = 0;
            let mut ref_gap_start = 0;
            let mut query_gap_end = query.len();
            let mut ref_gap_end = reference.len();
            if i > 0 {
                let prev = &selected[i - 1];
                query_gap_start = prev.query_pos + prev.length;
                ref_gap_start = prev.ref_pos + prev.length;
                // Check it is non-overlapping and colinear with the previous item
                let prev_query_end = prev.query_pos + prev.length;
                let prev_ref_end = prev.ref_pos + prev.length;
                if prev_query_end > anchor.query_pos || prev_ref_end >= anchor.ref_pos {
                    continue;
                }
            }
            if i < selected.len() {
                let next = &selected[i];
                query_gap_end = next.query_pos;
                ref_gap_end = next.ref_pos;
                // Check it is non-overlapping and colinear with the next item
                let anchor_query_end = anchor.query_pos + anchor.length;
                let anchor_ref_end = anchor.ref_pos + anchor.length;
                if anchor_query_end > next.query_pos || anchor_ref_end >= next.ref_pos {
                    continue;
                }
            }

            let query_gap = query_gap_end - query_gap_start;
            let ref_gap = ref_gap_end - ref_gap_start;

            let tolerance = ((query_gap + ref_gap) as f64).sqrt() as usize + 1;
            let delta_diagonal = (anchor.diagonal() as f64 - average_diagonal).abs();
            log::debug!(
                "Evaluating extra anchor: query_pos={}, ref_pos={}, length={}, query_gap={}, ref_gap={}, diagonal={}, delta_diagonal={}, tolerance={}",
                anchor.query_pos,
                anchor.ref_pos,
                anchor.length,
                query_gap,
                ref_gap,
                anchor.diagonal(),
                delta_diagonal,
                tolerance
            );
            if delta_diagonal > tolerance as f64 {
                log::debug!("  Rejected due to diagonal deviation");
                continue;
            }

            extras.push(anchor.clone());
            log::debug!(
                "Selected extra anchor: query_pos={}, ref_pos={}, length={}, diagonal={}",
                anchor.query_pos,
                anchor.ref_pos,
                anchor.length,
                anchor.diagonal()
            );
        }
        selected.extend(extras);
        selected.sort_by(Anchor::order_by_query_pos);
    }
    selected
}

/// Find all common substrings between query and reference that are at least min_length long.
/// Returns a vector of Anchor structs with query position, reference position, and length.
pub fn find_lcs_anchors(query: &[u8], reference: &[u8], min_length: usize) -> Vec<Anchor> {
    if query.is_empty() || reference.is_empty() || min_length == 0 {
        return Vec::new();
    }

    let combined = format!(
        "{}#{}",
        String::from_utf8_lossy(query),
        String::from_utf8_lossy(reference)
    );
    let st = SuffixTable::new(&combined);

    let sa = st.table(); // Suffix Array
    let lcp = st.lcp_lens(); // LCP Array
    let sep_idx = query.len();
    let mut anchors = Vec::new();

    // Linear scan of LCP array to find all cross-string matches
    for i in 1..sa.len() {
        let p1 = sa[i] as usize;
        let p2 = sa[i - 1] as usize;
        let match_len = lcp[i] as usize;

        // Check if suffixes come from different strings and meet min_length
        if match_len >= min_length {
            if p1 < sep_idx && p2 > sep_idx {
                // p1 is in query, p2 is in reference
                let ref_pos = p2 - sep_idx - 1; // adjust for separator
                anchors.push(Anchor::new(p1, ref_pos, match_len));
            } else if p2 < sep_idx && p1 > sep_idx {
                // p2 is in query, p1 is in reference
                let ref_pos = p1 - sep_idx - 1; // adjust for separator
                anchors.push(Anchor::new(p2, ref_pos, match_len));
            }
        }
    }

    // Sort by query position, then by length (descending)
    // Sort by diagonal, then by query_pos, then by length descending
    // This ensures that for each diagonal, we see longer anchors first at each position
    anchors.sort_by(|a, b| {
        a.diagonal()
            .cmp(&b.diagonal())
            .then(a.query_pos.cmp(&b.query_pos))
            .then(b.length.cmp(&a.length))
    });

    anchors
}

/// Remove anchors that are completely contained within another anchor on the same diagonal.
/// An anchor A is contained in anchor B if they're on the same diagonal and
/// A's query range [q, q+len) is a subset of B's query range.
pub fn remove_redundant_anchors(mut anchors: Vec<Anchor>) -> Vec<Anchor> {
    if anchors.len() <= 1 {
        return anchors;
    }

    let mut k = 0;
    let mut i = 0;

    while i < anchors.len() {
        let current = anchors[i].clone();
        let current_diag = current.diagonal();
        let current_end = current.query_pos + current.length;

        // This anchor is not redundant, keep it
        if k != i {
            anchors[k] = current.clone();
        }
        k += 1;

        // Skip all subsequent anchors on the same diagonal that are contained within this one
        let mut j = i + 1;
        while j < anchors.len() {
            let next = &anchors[j];
            if next.diagonal() != current_diag {
                // Different diagonal, stop skipping
                break;
            }
            let next_end = next.query_pos + next.length;
            if next.query_pos >= current.query_pos && next_end <= current_end {
                // next is fully contained in current, skip it
                j += 1;
            } else {
                // next extends beyond current, don't skip
                break;
            }
        }
        i = j;
    }
    anchors.truncate(k);

    anchors
}

#[cfg(test)]
mod tests {

    use super::*;

    /// Returns the longest anchor's substring and length.
    /// Helper for testing.
    fn extract_longest(query: &[u8], anchors: &[Anchor]) -> (usize, String) {
        if anchors.is_empty() {
            return (0, String::new());
        }

        // Find the longest anchor
        let longest = anchors.iter().max_by_key(|a| a.length).unwrap();

        let end = (longest.query_pos + longest.length).min(query.len());
        let substring = String::from_utf8_lossy(&query[longest.query_pos..end]).to_string();
        (longest.length, substring)
    }

    #[test]
    fn test_remove_redundant_anchors() {
        // Anchor at (26, 100, 89) covers query [26..115]
        // Anchor at (27, 101, 88) covers query [27..115] - contained in the first
        // Anchor at (27, 101, 87) covers query [27..114] - also contained
        let anchors = vec![
            Anchor::new(26, 100, 89), // covers [26..115]
            Anchor::new(27, 101, 88), // covers [27..115] - redundant
            Anchor::new(27, 101, 87), // covers [27..114] - redundant
            Anchor::new(28, 102, 87), // covers [28..115] - redundant
        ];

        let result = remove_redundant_anchors(anchors);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].query_pos, 26);
        assert_eq!(result[0].length, 89);
    }

    #[test]
    fn test_remove_redundant_keeps_non_overlapping() {
        // Two non-overlapping anchors on the same diagonal
        let anchors = vec![
            Anchor::new(0, 10, 5),  // covers [0..5]
            Anchor::new(10, 20, 5), // covers [10..15] - not overlapping
        ];

        let result = remove_redundant_anchors(anchors);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_remove_redundant_partial_overlap() {
        // Partially overlapping anchors - neither contains the other
        let anchors = vec![
            Anchor::new(0, 10, 10), // covers [0..10]
            Anchor::new(5, 15, 10), // covers [5..15] - extends beyond
        ];

        let result = remove_redundant_anchors(anchors);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_remove_redundant_different_diagonals() {
        // Same query range but different diagonals - both should be kept
        let anchors = vec![
            Anchor::new(10, 100, 50), // diagonal = 90
            Anchor::new(10, 110, 50), // diagonal = 100, same query range but different ref
        ];

        let result = remove_redundant_anchors(anchors);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_identical_sequences() {
        let query = b"ACGTACGT";
        let reference = b"ACGTACGT";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 8);
        assert_eq!(substr, "ACGTACGT");
    }

    #[test]
    fn test_no_common_substring() {
        let query = b"AAAA";
        let reference = b"TTTT";
        let anchors = find_lcs_anchors(query, reference, 1);
        assert!(anchors.is_empty());
    }

    #[test]
    fn test_common_prefix() {
        let query = b"ACGTXXXX";
        let reference = b"ACGTYYYY";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 4);
        assert_eq!(substr, "ACGT");
    }

    #[test]
    fn test_common_suffix() {
        let query = b"XXXXACGT";
        let reference = b"YYYYACGT";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 4);
        assert_eq!(substr, "ACGT");
    }

    #[test]
    fn test_common_middle() {
        let query = b"XXXACGTYYY";
        let reference = b"AAAACGTBBB";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 4);
        assert_eq!(substr, "ACGT");
    }

    #[test]
    fn test_multiple_common_substrings_finds_all() {
        // Has "ACG" and "TTTTTT" - should find both with min_length=3
        let query = b"ACGXXXTTTTTT";
        let reference = b"YYYACGTTTTTT";
        let anchors = find_lcs_anchors(query, reference, 3);
        // Should have at least 2 anchors (ACG and TTTTTT)
        assert!(
            anchors.len() >= 2,
            "Expected at least 2 anchors, got {}",
            anchors.len()
        );
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 6);
        assert_eq!(substr, "TTTTTT");
    }

    #[test]
    fn test_min_length_filters() {
        let query = b"ACGXXXTTTTTT";
        let reference = b"YYYACGTTTTTT";
        // With min_length=5, should only get TTTTTT
        let anchors = find_lcs_anchors(query, reference, 5);
        assert!(!anchors.is_empty());
        for anchor in &anchors {
            assert!(
                anchor.length >= 5,
                "Anchor length {} is less than min_length 5",
                anchor.length
            );
        }
    }

    #[test]
    fn test_single_base_match() {
        let query = b"A";
        let reference = b"A";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 1);
        assert_eq!(substr, "A");
    }

    #[test]
    fn test_query_is_substring_of_reference() {
        let query = b"ACGT";
        let reference = b"XXXACGTYYY";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 4);
        assert_eq!(substr, "ACGT");
    }

    #[test]
    fn test_reference_is_substring_of_query() {
        let query = b"XXXACGTYYY";
        let reference = b"ACGT";
        let anchors = find_lcs_anchors(query, reference, 1);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 4);
        assert_eq!(substr, "ACGT");
    }

    #[test]
    fn test_long_repeated_sequence() {
        let query = b"ACACACACACAC";
        let reference = b"GTGTGTACACACACACGT";
        let anchors = find_lcs_anchors(query, reference, 4);
        let (len, substr) = extract_longest(query, &anchors);
        assert!(
            len >= 8,
            "Expected at least 8 matching bases, got {}: {}",
            len,
            substr
        );
    }

    #[test]
    fn test_realistic_dna_sequences() {
        // Simulating a small insertion in the query
        let reference = b"ATCGATCGATCGATCG";
        let query = b"ATCGATCGXXXXATCGATCG";
        let anchors = find_lcs_anchors(query, reference, 4);
        let (len, substr) = extract_longest(query, &anchors);
        assert_eq!(len, 8);
        assert_eq!(substr, "ATCGATCG");
    }

    #[test]
    fn test_empty_query() {
        let query = b"";
        let reference = b"ACGT";
        let anchors = find_lcs_anchors(query, reference, 1);
        assert!(anchors.is_empty());
    }

    #[test]
    fn test_empty_reference() {
        let query = b"ACGT";
        let reference = b"";
        let anchors = find_lcs_anchors(query, reference, 1);
        assert!(anchors.is_empty());
    }

    #[test]
    fn test_anchor_positions_are_correct() {
        let query = b"XXXACGTYYY";
        let reference = b"AAACGTBBB";
        let anchors = find_lcs_anchors(query, reference, 4);
        assert!(!anchors.is_empty());
        let anchor = anchors.iter().find(|a| a.length == 4).unwrap();
        assert_eq!(anchor.query_pos, 3); // ACGT starts at position 3 in query
        assert_eq!(anchor.ref_pos, 2); // ACGT starts at position 2 in reference
        assert_eq!(anchor.length, 4);
    }

    #[test]
    fn test_returns_multiple_anchors() {
        // Two distinct matches: AAAA at different positions
        let query = b"AAAAXXXXXBBBB";
        let reference = b"AAAAYYYYYBBBB";
        let anchors = find_lcs_anchors(query, reference, 4);
        // Should find AAAA and BBBB
        assert!(
            anchors.len() >= 2,
            "Expected at least 2 anchors, got {:?}",
            anchors
        );
    }

    #[test]
    fn large_strings_1() {
        let reference = b"TTAGGCAGTGCCCCAGTGGGGACGCTGTGTGGGGTCTCCAATCCCACATTTCCCTTCTGCACTGCTCTAGCAGAGGTTCTCCATGAGGGCTCTTACCCTGCAGCAAACTTCTGCCTGGGCATTCAGGCATTTCTGTACAACCTCTGAAATCTAGGTGGAAGTTCCCAAACCTCAATTCTTGACTTCTGTGCACCCACAGGCTCAACACCACATAGGAGCTGCCAAGGCTTGGGGCTTGCACCCTCTGAAGCCACAGCCTGAGCTGTACTTTGGCTCCTTGTAGCCATGGCTAGAGTGGCTGGGACACAGGACACCAAATCCCTAGGCTGCACACAGCAGGTGGGCCCTAGGCCCTGCCCACAAAACAATTTTTTCCTCCTAGGACTCTGGGCCTGTGATGGGAAGGGCTGCTGTGAAGAAGACCTCTGACATGCCCTGGAGACATTTTCCCTATTGTCTTGGCAATTAACATTTAGCTCCTCATTAATCATGCAAATTTTTGCAGCCAGCTTGAATTTCTCCTCAGAAAATGGGTTTTTCTTTCCTATTACATTGTTAGGCTGCAAAATTTCCAAACTTTTATGCTCTGTTTCCCTTTTAAACTGAATGCTTTTTAACAGCACTCAAGTCACCTCTTGAATGCTTTGCTGCTTAGAAATTTCTTCTGCCAGATGCCCTAAATCATCTCCCTCAAGTTCAAAGTTCCACAGATCTTTAGGGTGGGAGCAAAATGTCACCAGTCTGTTTGCTAAAACATAGCAAGAGTCACCATTTCTCCCC";
        let query = b"CTAGGCAAGGGGATTCTCTGTGGGGGCTCACTCCCCATATTTCCCTTCCACATGGCCAGTAGAGGTTCTCCATGAGGGCTCTGCCCCTGCAGCAAACTTCTGCTTGGACATGCAGGCATTTCCATATGTCCTCTGAAATCTAGGCGGAGGTTCCCAAACCTCAATTCTTGACTCTGTGCACTGACAGGCTCAGCATCACATGAAAATCACAAGGCTTGGGGAGGGCTTATCCTTCAAGCAATGACCTCTTTGAGCCAGCTGGAGCTGAAGCAGCGGGAGTGAGGCACCATGTCCTGAGGCGGCACAGAGCAGGGCAGCCCTGGGCTCAGCCCAGGAAACCATTTTTCCCTACTAGGTGTCTGTGCCTGTGATGGGAGGGGTGGCCATGAAGACCTCTACATGGCCTGGAGACATTTTCTCCATTGCCTTAATGATTAACATTTGGCACATTTTTCAGATACATGTGGCTGTAAATATGTATCAATATCAGTGATATTTGCAGCTGGCTTGAATTTCTCTTCAAACAATGGGTTTTTCTTTAGTATTGCATCATCAGTACTGTAATTTTAAGTTTTTTTGCTTGTGCTTCCTCTTTTCATGCTTTCTTGAGAAATTTCTTCTGCCAGATACCCTAAATCATATCTGTCTCAAGTTCAACGTTCCACAGTCTCCAGGGCAGGGCAAGTCACCAGCACCAGTCTCTTCTGCCAAAGCATGCAAGAGTCACCTTTGCTCCAG";
        let mut anchors = find_lcs_anchors(query, reference, 20);
        anchors.sort_by(Anchor::order_by_length);
        for anchor in &anchors {
            println!(
                "Anchor: query_pos={}, ref_pos={}, length={}, diagonal={}",
                anchor.query_pos,
                anchor.ref_pos,
                anchor.length,
                anchor.diagonal()
            );
        }
        assert_eq!(anchors.len(), 9);
        assert_eq!(anchors[0].diagonal(), 11);
    }

    #[test]
    fn large_strings_2() {
        let reference = b"TTAGGCAGTGCCCCAGTGGGGACGCTGTGTGGGGTCTCCAATCCCACATTTCCCTTCTGCACTGCTCTAGCAGAGGTTCTCCATGAGGGCTCTTACCCTGCAGCAAACTTCTGCCTGGGCATTCAGGCATTTCTGTACAACCTCTGAAATCTAGGTGGAAGTTCCCAAACCTCAATTCTTGACTTCTGTGCACCCACAGGCTCAACACCACATAGGAGCTGCCAAGGCTTGGGGCTTGCACCCTCTGAAGCCACAGCCTGAGCTGTACTTTGGCTCCTTGTAGCCATGGCTAGAGTGGCTGGGACACAGGACACCAAATCCCTAGGCTGCACACAGCAGGTGGGCCCTAGGCCCTGCCCACAAAACAATTTTTTCCTCCTAGGACTCTGGGCCTGTGATGGGAAGGGCTGCTGTGAAGAAGACCTCTGACATGCCCTGGAGACATTTTCCCTATTGTCTTGGCAATTAACATTTAGCTCCTCATTAATCATGCAAATTTTTGCAGCCAGCTTGAATTTCTCCTCAGAAAATGGGTTTTTCTTTCCTATTACATTGTTAGGCTGCAAAATTTCCAAACTTTTATGCTCTGTTTCCCTTTTAAACTGAATGCTTTTTAACAGCACTCAAGTCACCTCTTGAATGCTTTGCTGCTTAGAAATTTCTTCTGCCAGATGCCCTAAATCATCTCCCTCAAGTTCAAAGTTCCACAGATCTTTAGGGTGGGAGCAAAATGTCACCAGTCTGTTTGCTAAAACATAGCAAGAGTCACCATTTCTCCCC";
        let query = b"CTAGGCAAGGGGATTCTCTGTGGGGGCTCACTCCCCATATTTCCCTTCCACATGGCCAGTAGAGGTTCTCCATGAGGGCTCTGCCCCTGCAGCAAACTTCTGCTTGGACATGCAGGCATTTCCATATGTCCTCTGAAATCTAGGCGGAGGTTCCCAAACCTCAATTCTTGACTCTGTGCACTGACAGGCTCAGCATCACATGAAAATCACAAGGCTTGGGGAGGGCTTATCCTTCAAGCAATGACCTCTTTGAGCCAGCTGGAGCTGAAGCAGCGGGAGTGAGGCACCATGTCCTGAGGCGGCACAGAGCAGGGCAGCCCTGGGCTCAGCCCAGGAAACCATTTTTCCCTACTAGGTGTCTGTGCCTGTGATGGGAGGGGTGGCCATGAAGACCTCTACATGGCCTGGAGACATTTTCTCCATTGCCTTAATGATTAACATTTGGCACATTTTTCAGATACATGTGGCTGTAAATATGTATCAATATCAGTGATATTTGCAGCTGGCTTGAATTTCTCTTCAAACAATGGGTTTTTCTTTAGTATTGCATCATCAGTACTGTAATTTTAAGTTTTTTTGCTTGTGCTTCCTCTTTTCATGCTTTCTTGAGAAATTTCTTCTGCCAGATACCCTAAATCATATCTGTCTCAAGTTCAACGTTCCACAGTCTCCAGGGCAGGGCAAGTCACCAGCACCAGTCTCTTCTGCCAAAGCATGCAAGAGTCACCTTTGCTCCAG";
        let mut anchors = find_lcs_anchors(query, reference, 15);
        anchors.sort_by_key(|a| a.diagonal());
        let mut best_diagonal = 0isize;
        let mut best_weight = f64::MIN;
        for ancs in anchors.group_by(|a| a.diagonal()) {
            let diag = ancs.0;
            let weight: f64 = ancs
                .1
                .iter()
                .map(|a| {
                    let l = a.length as f64;
                    l * l.log2()
                })
                .sum();
            if weight > best_weight {
                best_weight = weight;
                best_diagonal = diag;
            }
        }
        anchors.retain(|a| a.diagonal() == best_diagonal);
        for anchor in &anchors {
            println!(
                "Anchor: query_pos={}, ref_pos={}, length={}, diagonal={}",
                anchor.query_pos,
                anchor.ref_pos,
                anchor.length,
                anchor.diagonal()
            );
        }
        assert_eq!(best_diagonal, 11);
        assert_eq!(anchors.len(), 24);
        assert_eq!(anchors[0].diagonal(), 11);
    }

    #[test]
    fn long_string_anchors() {
        let reference = b"ACATAAAAACTAGATGGAAGCATTCTCAGAAACTACTTTGTGATGATTGCATTCGACTCACAGAGTTGAACATTCCTATAGATAGAGCAGGTTGTAAACAATCTTTTTGTAGAATCTGCGATTGGACATTTGGAATGCTTTGAGGCCTACTGTAGTAAAGGAAATAACTTCATCTAAAAACCAAACGGACGCATTCACAGTACAATTCTTAGTGATCATTGGATTGAACTAACAGAGCTGAACATTCCTTTAGATGGAGCAGTTTCCAAACCCACTTTCTGTAGAATCTGCAAGTGGATATTTGGACTTCTCTGAGGATTTCGTTGGAAACGGGATATACTTCCCAGAACTACACGGAGCATTGTGAGAAACTTCTTTGTGATGTTTGCATTCAACTCACAGAGTTGAACCTTGCTTTCATAGTTCAGCTTTCAAACACTCCTTTTGTAGAATCTGCAAGTGGATATTTGGGCCACTTTGTGGCCTTCCTTCGAAACGGGTATATCTTCACATCAAACCTAGACAGAAGCATTCTCAGAATGTTTCCTGTGATGACTGCATTCAACTCACAGAGGTGAACAATCCTGTTGATGCAGCAGTTTTGAAACTCTCTTTCTTTGGATTCTGCAAGTTGATATGTGGACCTCTGTGAAGATTTCGTTGGAAATGGGTTCATCTTCACAGAAAAACTAAACAGAAGCATTCTCAGAAACTGCTTTGTGATGTTTGTGTTCCACTTCAAGAATTGAACTTTCCTCTTGACAGAGCAGCTCTGAAACCCTCTTTTTCTAGAATCTGCAAGTGGACATTTGGAGGGCTTTGAGGCCTGTGGTGGAAAAGGAAAATCTTCCCATAAAAACTAGATGGAAGCATTCTCAGAAACTACTTTGTGATGATTGCATTCGACTCACAGAGTTGAACATTCCTATAGATAGAGCAGGTTGTAAACAATCTTTTTGTAGAATCTGCGATTGGAGATTTGGACTGCTTTGAGGCCTTCTGTAGTAAAGGAAATAACTTCATCTAAAAACCAAACGGAAGCATTCACAGACAATTCTTAGTGATCATTGCATTGAACTAACAGAGCTGAACATTCCTTTAGATGGAGCAGTTTCCAAACCCACTTTCTGTAGAATCTGCAAGTGGATATTTGGACTTCTCTGAGGATTTCGTTGGAAACGGGATAAACTTCCCAGAACTACACGGAAGCATTGTGAGAAACTTCTTTGTGATGTTTGCATTCAACTCACAGAGTTGAACCTTGCTTTGATAGTTCAGCTTTCAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGACCACTTTGTGGCCTTCCTTCGAAAAGGCTATATCTTCACATCAAACCTAGACAGAAGCATTCTCAGAATGTTTCCTGTGATGACTGCATTCAACTCACAGAGGTGAACAATCCTGTTGATGGAGCAGTTTTGAAACTCTCTTTCTTTGGATTCTGCAAGTGGATATGTGGACCTCTGTGAAGATTTCGTTGGAAACGGGTTCATCTTCACAGAAAAACTAAACAGAAGCATTCTCAGAAACTGCTTTTTGATGTTTGTGTTCCACTTCAAGAATTGAACTTTCCTCTTGACAGAGCAGCTCTGAAACCCTCTTTTTCTAGAATCTGCAAGTGGACATTTGGAGGGCTTTGAGGCCTGTGGTGGAAAAGGAAAATCTTCACATAAAAACTAGATGGAAGCATTCTCAGAAACTACTTTGTGATGATTGCATTCGACTCACAGAGTTGAACATTCCTATAGATAGAGCAGGTTGTAAACAATCTTTTTGTAGAATCTGCGATTGGAGATTTGGACTGCTTTGAGGCCTACTGTAGTAAAGGAAATAACTTCATCTAAAAACCAAACGGAAGCATTCACAGACAATTCTTAGTGATCATTGGATTGAACTAACAGAGCTGAACATTCCTTTAGATGGCGCAGTTTCCATACACACTTTCTGTAGAATCTGCAAGTGGATATTTGGACCTCTCTGAGGATTTCGTTGGAAACGGGATAAATTTCCCAGAACTACACGGAAGCATTCTGAGAAACTTCTTTGTGATGTTTGCATTCAACTCACAGAGTTGAACCTTGCTTTCATAGTTCAGCTTTCAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGACCACTTTGTGGCCTTCCTTCGAAACGGGTATATCTTCACATCAAACCTAGACAGAAGCATTCTCAGAATGTTTCCTGGGATGACTGCATTCAACTCACAGAGGTGAACAATCCTGCTGATGGAGCAGTTTTGAAACTCTCTTTCTTTGGATTCTGCAAGTGGATATGTGGACCTCTGTGAAGATTTCGTTGGAAACGGGTTCATCTTCACAGAAAAACTAAACAGGAGCATTCTCAGAAACTGCTTTGTGATGTTTGTGTTCCACTTCAG";
        let query = b"CCATAAAAACTAGATGTTAGCATTCTCAGAAACTACTTTGTGATGATTGCATTCGACTCACAGAGTTGAACATTCCTATAGATAGAGCAGGTTGTAAACAATCTTTTTGTAGAATATGCGATTGGAGATTTGGACTGCTTTGAGGCCTACTGTAGTAAAGGAAATAACTTCATCTAAAAACCAAACGGAAGCATTCACAGACAATTCTTAGTGATCATTGCATTGATCTAACAGAGCTGAACATTCCTTTAGATGGCGCAGTTTCCAAACACACTTTCTGTAGAATCTGCAAGTGGATATTTGGACCTCTCTGAGGATTTCGTTGGAAACGGGATAAACTTCCCAGAACTACACGGAAGCATTCTGAGAAACTTCTTTGTGATGTTTGCATTCAACTCACAGAGTTGAACCTTGCTTTCATAGTTCAGCTTTCAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGACCACTTTGTGGCCTTCCTTCGAAACGGGTATATCTTCACATCAAACCTAGACAGAAGCATTCTCAGAATGTTTCCTGTGATGACTGCATTCAACTCACAGAGGTGAACAATCCTGCTGATGGAGCAGTTTTGAAACTCTCTTTCTTTGGATTCTGCAAGTGGATATGTGGACTTCTGTGAAGATTTCGTTGGAAACGGGTTCATCTTCACAGAAAAACTAAACAGGAGCATTCGCAGAAACTGCTTTGTGATGTTTGTGTTCCACTTCAAGAATTGAACTTTCCTCTTGACAGAGCAGCTCTGAAACCCTCTTTTTCTAGAATCTGCAAGTGGACATTTGGAGGGCTTTGAGGCCTGTGGTGGAAAAGGAAAATCTTCACATAAAAACTAGATGGAAGCATTCTCAGAAACTACTTTGTGATGATTGCATTCGACTCACAGAGTTGAACATTCGTATAGATAGAGCAGGTTGTAAACAATCTTTTTGTAGAATCTGCGATTGGAGATTTGGACTGCTTTGAGGCCTACTGTAGTAAAGGAAATAACTTCATCTAAAAACCAAACGGAAGCATTCACAGACAATTCTTAGTGATCATTGCATTGAACTAACAGAGCTGAACATTCCTTTAGATGGCGCAGTTTCCAAACACACTTTCTGTAGAATCTGCAAGTGGATATTTGGACCTCTCTGAGGATTTCGTTGGAAACGGGATAAACTTCCCAGAACTACACGGAAGCATTGTGAGAAACTTCTTTGTGATGTTTGCATTCAACTCACAGAGTTGAACCTTGCTTTCATAGTTCAGCTTTCAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGACCACTTTGTGGCCTCCTTCGAAACGGGTATATCTTCACATCAAACCTAGACAGAAGCATTCTCAGAATGTTTCCTGTGATGATTGCATTCAACTCACAGAGGTGAACAATCCTGTTGATGGAGCAGTTTTGAAACTCTCTTTCTTTGGATTCTGCAAGTGGATATGTGGACCTCTGTGAAGATTTCGTTGGAAACGGGTTCATCTTCACAGAAAAACTAAACAGAAGCATTCTCAGAAACTGCTTTGTGATGCTTGTGTTCCACTTCAGGAATTGAACTTTCCTCTTGAAAGAGCAGCTCTGAAACCCTCTTTTTCTAGAATCTGCAAGTGGACATTTGGAGGGCTTTGAGGACTGTGGTGGAAAAGGAAAATCTTCCCATAAAAACTAGATGTTAGCATTCTCAGAAACTACTTTGTGATGATTGCATTCGACTCACAGAGTTGAACATTCTTATAGATAGAGCAGGTTGTAAACAATCTTTTTGTAGAATCTGCGATTGGAGATTTGGACTGCTTTGAGGCCTACTGTAGTAAAGGAAATAACTTCATCTAAAAGCCAAACGGAAGCATTCACAGACAATTCTTAGTGATCATTGGATTGAACTAACAGAGCTGAACATTCCTTTAGATGGAGCAGTTTCCAAACACACTTTCTGTAGAATCTGCAAGTGGATATTTGGACTTCTCTGAGGATTTCGTTGGATACGGGATAAACTTCCCAGAACTACACGGAAGCATTCTGAGAAACTTCTTTGTGATGTTTGCATTCAACTCACAGAGTTGAACCTTGCTTTCATAGTTCAGCTTTCAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGACCACTTTGTGGCCTTCCTTCGAAACGGGTATATCTTCACATCAAACCTAGACAGAAGCATTCTCAGAATGTTTCCTGTGATGATTGCATTCAACTCACAGAGGTGAACAATCCTGTTGATGGAGGAGGTTTGAAACTCTCTTTCTTTGGATTCTGCAAGTGGATATGTGGACCTCTGTGAAGATTTCGTTGGAAACGGGTTCATCTTCACAGAAAAACTAAACAGAAGCATTCTCAGAAACTGCTTTGTGATGTTTGTGTTCCACTTCAT";
        let anchors = select_lcs_anchors(query, reference, 15);

        assert_eq!(anchors.len(), 17);
    }
}
