use std::cmp::Reverse;

use crate::{
    align::Anchor,
    kmers::Kmer,
    utils::{longest_common_prefix, range_set::RangeSet},
};

pub fn select_kmer_anchors(read_seq: &[u8], ref_seq: &[u8], min_length: usize) -> Vec<Anchor> {
    let mut finder = KmerAnchorFinder::<16>::new(min_length);
    finder.find_anchors(read_seq, ref_seq)
}

struct KmerAnchorFinder<const K: usize> {
    min_length: usize,
    query_kmers: Vec<u64>,
    query_positions: Vec<usize>,
    query_permutation: Vec<isize>,
    ref_kmers: Vec<u64>,
    ref_positions: Vec<usize>,
    ref_permutation: Vec<isize>,
    tmp_kmers: Vec<u64>,
    tmp_positions: Vec<usize>,
}

impl<const K: usize> KmerAnchorFinder<K> {
    pub fn new(min_length: usize) -> Self {
        Self {
            min_length,
            query_kmers: Vec::new(),
            query_positions: Vec::new(),
            query_permutation: Vec::new(),
            ref_kmers: Vec::new(),
            ref_positions: Vec::new(),
            ref_permutation: Vec::new(),
            tmp_kmers: Vec::new(),
            tmp_positions: Vec::new(),
        }
    }

    pub fn find_anchors(&mut self, read_seq: &[u8], ref_seq: &[u8]) -> Vec<Anchor> {
        Self::extract_and_sort_kmers(
            read_seq,
            &mut self.query_kmers,
            &mut self.query_positions,
            &mut self.query_permutation,
            &mut self.tmp_kmers,
            &mut self.tmp_positions,
        );
        Self::extract_and_sort_kmers(
            ref_seq,
            &mut self.ref_kmers,
            &mut self.ref_positions,
            &mut self.ref_permutation,
            &mut self.tmp_kmers,
            &mut self.tmp_positions,
        );

        let mut diagonals = Vec::new();

        let mut q_idx = 0;
        let mut r_idx = 0;
        while q_idx < self.query_kmers.len() && r_idx < self.ref_kmers.len() {
            let q_kmer = self.query_kmers[q_idx];
            let r_kmer = self.ref_kmers[r_idx];
            if q_kmer < r_kmer {
                q_idx += 1;
                continue;
            }
            if r_kmer < q_kmer {
                r_idx += 1;
                continue;
            }
            // q_kmer == r_kmer

            // First, find the range of equal k-mers in the query
            let q_start = q_idx;
            while q_idx < self.query_kmers.len() && self.query_kmers[q_idx] == q_kmer {
                q_idx += 1;
            }
            let q_end = q_idx;

            // Next, find the range of equal k-mers in the reference
            let r_start = r_idx;
            while r_idx < self.ref_kmers.len() && self.ref_kmers[r_idx] == r_kmer {
                r_idx += 1;
            }
            let r_end = r_idx;

            let q_count = q_end - q_start;
            let r_count = r_end - r_start;

            if q_count * r_count > 25 {
                // Too many hits, skip
                continue;
            }

            for qi in q_start..q_end {
                let q_pos = self.query_positions[qi];
                for ri in r_start..r_end {
                    let r_pos = self.ref_positions[ri];
                    let diag = q_pos as isize - r_pos as isize;
                    diagonals.push((diag, qi, ri));
                }
            }
        }

        let mut seeds = Vec::new();

        let mut diag_index = 0;
        while diag_index < diagonals.len() {
            let diag = diagonals[diag_index].0;
            let start_index = diag_index;
            while diag_index < diagonals.len() && diagonals[diag_index].0 == diag {
                diag_index += 1;
            }
            let end_index = diag_index;

            // Now we have all hits for this diagonal in diagonals[start_index..end_index]
            let hits: &mut [(isize, usize, usize)] = &mut diagonals[start_index..end_index];
            if hits.len() < 2 {
                continue;
            }

            hits.sort_by_key(|(_diag, i, _j)| self.query_positions[*i]);

            // Find overlapping seeds and use them to compose maximal seeds.
            let mut seed_start = 0;
            while seed_start < hits.len() {
                let (_diag, first_qi, first_ri) = hits[seed_start];
                let first_qpos = self.query_positions[first_qi];
                let first_rpos = self.ref_positions[first_ri];
                let mut seed_end = seed_start + 1;
                let mut last_qend = first_qpos + K;
                let mut last_rend = first_rpos + K;
                while seed_end < hits.len() {
                    let (_diag, qi, ri) = hits[seed_end];
                    let qpos = self.query_positions[qi];
                    let rpos = self.ref_positions[ri];
                    if qpos < last_qend && rpos < last_rend {
                        // Overlaps with previous seed
                        last_qend = last_qend.max(qpos + K);
                        last_rend = last_rend.max(rpos + K);
                        seed_end += 1;
                    } else {
                        break;
                    }
                }
                // We have a seed from seed_start..seed_end
                let qpos = self.query_positions[first_qi];
                let rpos = self.ref_positions[first_ri];
                let match_len = last_qend - qpos;
                let extension = longest_common_prefix(
                    &read_seq[qpos + match_len..],
                    &ref_seq[rpos + match_len..],
                );
                let match_len = match_len + extension;
                if match_len >= self.min_length {
                    seeds.push(Anchor::new(qpos, rpos, match_len));
                }
                seed_start = seed_end;
            }
        }

        seeds.sort_by_key(|s| (Reverse(s.length), s.ref_pos));

        let mut wanted: Vec<Anchor> = Vec::new();
        let mut occupied_query = RangeSet::new();
        let mut occupied_ref = RangeSet::new();
        for seed in seeds {
            if seed.length == 0 {
                continue;
            }
            let q_start = seed.query_pos;
            let q_end = seed.query_pos + seed.length;
            let r_start = seed.ref_pos;
            let r_end = seed.ref_pos + seed.length;

            if occupied_query.contains_overlap(q_start, q_end)
                || occupied_ref.contains_overlap(r_start, r_end)
            {
                continue;
            }

            // Check seed is colinear with existing seeds
            let mut colinear = true;
            for existing_seed in &wanted {
                if (seed.query_pos < existing_seed.query_pos
                    && seed.ref_pos > existing_seed.ref_pos)
                    || (seed.query_pos > existing_seed.query_pos
                        && seed.ref_pos < existing_seed.ref_pos)
                {
                    colinear = false;
                    break;
                }
            }
            if !colinear {
                continue;
            }

            // No overlap, take this seed
            occupied_query.add_range(q_start, q_end);
            occupied_ref.add_range(r_start, r_end);
            wanted.push(seed);
        }

        wanted.sort_by(|a, b| a.query_pos.cmp(&b.query_pos));

        wanted
    }

    pub fn extract_and_sort_kmers(
        seq: &[u8],
        kmers: &mut Vec<u64>,
        positions: &mut Vec<usize>,
        permutation: &mut Vec<isize>,
        tmp_kmers: &mut Vec<u64>,
        tmp_positions: &mut Vec<usize>,
    ) {
        kmers.clear();
        positions.clear();
        Kmer::<K>::kmerize_fwd(seq, |pos, kmer| {
            kmers.push(kmer.0);
            positions.push(pos);
        });
        permutation.clear();
        permutation.extend(0..kmers.len() as isize);
        permutation.sort_by_key(|&i| kmers[i as usize]);

        tmp_kmers.clear();
        tmp_positions.clear();
        for i in permutation.iter() {
            tmp_kmers.push(kmers[*i as usize]);
            tmp_positions.push(positions[*i as usize]);
        }
        std::mem::swap(kmers, tmp_kmers);
        std::mem::swap(positions, tmp_positions);
    }
}
