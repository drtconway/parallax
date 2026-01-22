use std::{cmp::Reverse, collections::HashMap};

use crate::{
    align::{AlignParams, Alignment, WfAligner},
    kmers::Kmer,
    utils::{longest_common_prefix, range_set::RangeSet},
};

pub struct MiniAligner<const K: usize> {
    query_kmers: Vec<u64>,
    query_positions: Vec<usize>,
    query_permutation: Vec<isize>,
    ref_kmers: Vec<u64>,
    ref_positions: Vec<usize>,
    ref_permutation: Vec<isize>,
    tmp_kmers: Vec<u64>,
    tmp_positions: Vec<usize>,
}

impl<const K: usize> MiniAligner<K> {
    pub fn align(&mut self, read_seq: &[u8], ref_seq: &[u8]) -> Option<Alignment> {
        self.query_kmers.clear();
        self.query_positions.clear();
        self.query_permutation.clear();
        Kmer::<K>::kmerize_fwd(read_seq, |pos, kmer| {
            self.query_kmers.push(kmer.0);
            self.query_positions.push(pos);
        });
        self.query_permutation
            .extend(0..self.query_kmers.len() as isize);
        self.query_permutation
            .sort_by_key(|&i| self.query_kmers[i as usize]);

        self.tmp_kmers.clear();
        self.tmp_positions.clear();
        for &i in &self.query_permutation {
            self.tmp_kmers.push(self.query_kmers[i as usize]);
            self.tmp_positions.push(self.query_positions[i as usize]);
        }
        std::mem::swap(&mut self.query_kmers, &mut self.tmp_kmers);
        std::mem::swap(&mut self.query_positions, &mut self.tmp_positions);
        for i in 1..self.query_kmers.len() {
            assert!(self.query_kmers[i - 1] <= self.query_kmers[i]);
        }

        self.ref_kmers.clear();
        self.ref_positions.clear();
        self.ref_permutation.clear();
        Kmer::<K>::kmerize_fwd(ref_seq, |pos, kmer| {
            self.ref_kmers.push(kmer.0);
            self.ref_positions.push(pos);
        });
        self.ref_permutation
            .extend(0..self.ref_kmers.len() as isize);
        self.ref_permutation
            .sort_by_key(|&i| self.ref_kmers[i as usize]);

        self.tmp_kmers.clear();
        self.tmp_positions.clear();
        for &i in &self.ref_permutation {
            self.tmp_kmers.push(self.ref_kmers[i as usize]);
            self.tmp_positions.push(self.ref_positions[i as usize]);
        }
        std::mem::swap(&mut self.ref_kmers, &mut self.tmp_kmers);
        std::mem::swap(&mut self.ref_positions, &mut self.tmp_positions);
        for i in 1..self.ref_kmers.len() {
            assert!(self.ref_kmers[i - 1] <= self.ref_kmers[i]);
        }

        let mut diagonals: Vec<(isize, usize, usize)> = Vec::new();

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

        diagonals.sort_unstable();

        let mut mini_seeds = Vec::new();

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

            hits.sort_by_key(|(_diag, i, j)| self.query_positions[*i]);

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
                mini_seeds.push(MiniSeed::new(
                    qpos,
                    rpos,
                    self.query_kmers[first_qi],
                    match_len,
                ));
                seed_start = seed_end;
            }
        }

        mini_seeds.sort_by_key(|s| (Reverse(s.match_len), s.ref_pos));

        let mut wanted: Vec<MiniSeed> = Vec::new();
        let mut occupied_query = RangeSet::new();
        let mut occupied_ref = RangeSet::new();
        for seed in mini_seeds {
            if seed.match_len == 0 {
                continue;
            }
            let q_start = seed.query_pos;
            let q_end = seed.query_pos + seed.match_len;
            let r_start = seed.ref_pos;
            let r_end = seed.ref_pos + seed.match_len;

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

        wanted.sort_by_key(|s| s.query_pos);

        let mut alignments: Vec<Alignment> = Vec::new();
        let n = wanted.len();
        for i in 0..n {
            let seed = &wanted[i];
            if i == 0 {
                // For the first seed, align from start of read/ref
                let q_start = 0;
                let q_end = seed.query_pos;
                let r_start = 0;
                let r_end = seed.ref_pos;
                if q_end > q_start || r_end > r_start {
                    let query = &read_seq[q_start..q_end];
                    let reference = &ref_seq[r_start..r_end];
                    //println!("({}) Aligning initial gap: read[{}:{}], ref[{}:{}]", i, q_start, q_end, r_start, r_end);
                    let aln = WfAligner::new(AlignParams::default()).align(query, reference)?;
                    alignments.push(aln);
                }
                continue;
            }

            // Align up to this seed
            let prev = &wanted[i - 1];
            let q_start = prev.query_pos + prev.match_len;
            let q_end = seed.query_pos;
            let r_start = prev.ref_pos + prev.match_len;
            let r_end = seed.ref_pos;
            if q_end > q_start || r_end > r_start {
                let query = &read_seq[q_start..q_end];
                let reference = &ref_seq[r_start..r_end];
                //println!("({}) Aligning initial gap: read[{}:{}], ref[{}:{}]", i, q_start, q_end, r_start, r_end);
                let aln = WfAligner::new(AlignParams::default()).align(query, reference)?;
                alignments.push(aln);
            }

            // Add the seed match as a perfect alignment
            let seed_query = &read_seq[seed.query_pos..(seed.query_pos + seed.match_len)];
            let seed_ref = &ref_seq[seed.ref_pos..(seed.ref_pos + seed.match_len)];
            assert_eq!(seed_query, seed_ref);
            let seed_aln = Alignment::from_perfect_match(seed_query.len());
            alignments.push(seed_aln);

            if i == n - 1 {
                // For the last seed, align to end of read/ref
                let q_start = seed.query_pos + seed.match_len;
                let q_end = read_seq.len();
                let r_start = seed.ref_pos + seed.match_len;
                let r_end = ref_seq.len();
                if q_end > q_start || r_end > r_start {
                    let query = &read_seq[q_start..q_end];
                    let reference = &ref_seq[r_start..r_end];
                    //println!("({}) Aligning initial gap: read[{}:{}], ref[{}:{}]", i + 1, q_start, q_end, r_start, r_end);
                    let aln = WfAligner::new(AlignParams::default()).align(query, reference)?;
                    alignments.push(aln);
                }
            }
        }

        let final_alignment = Alignment::concat(&alignments);

        Some(final_alignment)
    }
}

struct MiniSeed {
    query_pos: usize,
    ref_pos: usize,
    kmer: u64,
    match_len: usize,
}

impl MiniSeed {
    pub fn new(query_pos: usize, ref_pos: usize, kmer: u64, match_len: usize) -> Self {
        Self {
            query_pos,
            ref_pos,
            kmer,
            match_len,
        }
    }
}

impl Default for MiniAligner<11> {
    fn default() -> Self {
        Self {
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn align_test_1() {
        let mut aligner: MiniAligner<11> = MiniAligner::default();

        let query = b"ATTTTTATATATACTTATATATTTATATATATTTTTATATATACTCATATATTTATATATATTTTATATATACTTATTTATATATATATATTTTTATATATATTTAATTTTTACATATATTTATATTTTTATATATTTATATATTTATATATTTTTATATTTTATATATATGTTTATATATTTATATATTATATATATTTATATATATTTATATATTTATATATTATATATTTATATATATTTATATATTTATATATTATATATATTTATATATTTATATATTTATATATTACATATATTTATATATATTTATATATTTATATATGTTTATATATTTATATATTATATATATTTATATATATTTATATTATATATATACTTATATATTTATATATATTTTTATATATACTTATATATTTATATATATTTTTATATATACTTATATATTTATATATATTTTTATATATACTTATATATATTTTTTATATATTTATATATTTTTATATATATTTAATTTTTAT";
        let reference = b"TTTTATATATACTTATATATTTATATATATTTTTATATATACTCATATATTTATATATATTTTATATATACTTATTTTATATATATATATTTTTATATATATTTAATTTTTAC";
        let result = aligner.align(query, reference);
        let alignment = result.expect("alignment failed");
        let (r, a, q) = alignment.blast_style(reference, query);
        println!("REF: {}", r);
        println!("ALN: {}", a);
        println!("QRY: {}", q);
        assert_eq!(alignment.score, 844);
    }
}
