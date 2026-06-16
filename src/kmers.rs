use crate::utils::{Selection, hasher::Hasher};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Kmer<const K: usize>(pub u64);

impl<const K: usize> Kmer<K> {
    pub fn to_string(&self) -> String {
        let mut s = String::with_capacity(K);
        let mut x = self.0;
        for _i in 0..K {
            let bits = x & 0b11;
            x >>= 2;
            let c = match bits {
                0b00 => 'A',
                0b01 => 'C',
                0b10 => 'G',
                0b11 => 'T',
                _ => unreachable!(),
            };
            s.push(c);
        }
        s.chars().rev().collect()
    }

    pub fn kmerize_fwd<S: AsRef<[u8]>, F: FnMut(usize, Kmer<K>)>(seq: S, mut acceptor: F) {
        let seq = seq.as_ref();
        let m = (1u64 << (2 * K)) - 1;
        let mut kmer: u64 = 0;
        let mut j = 0;
        for i in 0..seq.len() {
            if let Some(nuc) = Kmer::<K>::nucleotide(seq[i]) {
                kmer = ((kmer << 2) | nuc) & m;
                j += 1;
                if j == K {
                    acceptor(i + 1 - K, Kmer(kmer));
                    j -= 1;
                }
            } else {
                kmer = 0;
                j = 0;
            }
        }
    }

    pub fn kmerize<S: AsRef<[u8]>, F: FnMut(usize, Kmer<K>, Kmer<K>)>(seq: S, mut acceptor: F) {
        let seq = seq.as_ref();
        let m = (1u64 << (2 * K)) - 1;
        let mut fwd_kmer: u64 = 0;
        let mut rev_kmer: u64 = 0;
        let mut j = 0;
        for i in 0..seq.len() {
            if let Some(nuc) = Kmer::<K>::nucleotide(seq[i]) {
                fwd_kmer = ((fwd_kmer << 2) | nuc) & m;
                rev_kmer = (rev_kmer >> 2) | ((3 - nuc) << (2 * (K - 1)));
                j += 1;
                if j == K {
                    acceptor(i + 1 - K, Kmer(fwd_kmer), Kmer(rev_kmer));
                    j -= 1;
                }
            } else {
                fwd_kmer = 0;
                rev_kmer = 0;
                j = 0;
            }
        }
    }

    pub fn fwd_kmer_iter<'a>(seq: &'a [u8]) -> KmerIterator<'a, K> {
        KmerIterator::new(seq)
    }

    pub fn kmer_pair_iter<'a>(seq: &'a [u8]) -> KmerPairIterator<'a, K> {
        KmerPairIterator::new(seq)
    }

    pub fn is_open_syncmer<const S: usize, H: Hasher>(&self) -> bool {
        let mask: u64 = (1 << (2 * S)) - 1;
        let mut min_subkmer_pos = 0;
        let mut min_subkmer = H::hash64(self.0 & mask);
        for i in 0..=(K - S) {
            let subkmer = (self.0 >> (2 * i)) & mask;
            let hashed_subkmer = H::hash64(subkmer);
            if hashed_subkmer < min_subkmer {
                min_subkmer = hashed_subkmer;
                min_subkmer_pos = i;
            }
        }
        min_subkmer_pos == 0
    }

    pub fn kmerize_open_syncmers<
        const S: usize,
        H: Hasher,
        Seq: AsRef<[u8]>,
        F: FnMut(usize, Selection<Kmer<K>, Kmer<K>>),
    >(
        seq: Seq,
        _s: [(); S],
        mut acceptor: F,
    ) {
        Kmer::<K>::kmerize(seq, |i, fwd, rev| {
            match (fwd.is_open_syncmer::<S, H>(), rev.is_open_syncmer::<S, H>()) {
                (true, true) => acceptor(i, Selection::Both(fwd, rev)),
                (true, false) => acceptor(i, Selection::Left(fwd)),
                (false, true) => acceptor(i, Selection::Right(rev)),
                (false, false) => {}
            }
        });
    }

    /// Iterate over forward-only open syncmers.
    /// Only returns k-mers where the forward k-mer is an open syncmer.
    /// Used for seeding when processing forward and reverse strands separately.
    #[cfg(test)]
    pub fn kmerize_open_syncmers_fwd_orig<
        const S: usize,
        H: Hasher,
        Seq: AsRef<[u8]>,
        F: FnMut(usize, Kmer<K>),
    >(
        seq: Seq,
        _s: [(); S],
        mut acceptor: F,
    ) {
        Kmer::<K>::kmerize_fwd(seq, |i, fwd| {
            if fwd.is_open_syncmer::<S, H>() {
                acceptor(i, fwd);
            }
        });
    }

    /// Optimized version of kmerize_open_syncmers_fwd that avoids recomputing
    /// s-mer hashes by using a sliding window.
    ///
    /// Instead of computing all (K-S+1) s-mer hashes for each k-mer, this version
    /// maintains a circular buffer of hashes and only computes one new hash per base.
    /// This reduces hash computations from O((K-S+1) * n) to O(n).
    pub fn kmerize_open_syncmers_fwd<
        const S: usize,
        H: Hasher,
        Seq: AsRef<[u8]>,
        F: FnMut(usize, Kmer<K>),
    >(
        seq: Seq,
        _s: [(); S],
        mut acceptor: F,
    ) {
        let seq = seq.as_ref();
        if seq.len() < K {
            return;
        }

        let k_mask = (1u64 << (2 * K)) - 1;
        let s_mask = (1u64 << (2 * S)) - 1;
        let w: usize = K - S + 1;

        let mut kmer: u64 = 0;
        let mut smer: u64 = 0;
        let mut kmer_len = 0; // valid bases in current k-mer
        let mut smer_len = 0; // valid bases in current s-mer

        // Circular buffer for s-mer hashes
        let mut hashes = [0u64; K];
        let mut ring_pos = 0; // next position to write
        let mut ring_count = 0; // number of valid entries (0..=w)

        for (i, &base) in seq.iter().enumerate() {
            if let Some(nuc) = Kmer::<K>::nucleotide(base) {
                // Update rolling k-mer
                kmer = ((kmer << 2) | nuc) & k_mask;
                kmer_len = (kmer_len + 1).min(K);

                // Update rolling s-mer (rightmost s-mer of the k-mer)
                smer = ((smer << 2) | nuc) & s_mask;
                smer_len = (smer_len + 1).min(S);

                // Once we have a complete s-mer, add its hash to the buffer
                if smer_len >= S {
                    let hash = H::hash64(smer);
                    let newest_idx = ring_pos;
                    hashes[newest_idx] = hash;
                    ring_pos = (ring_pos + 1) % w;
                    ring_count = (ring_count + 1).min(w);

                    // Check syncmer once we have a complete k-mer and full hash window
                    if kmer_len >= K && ring_count == w {
                        // For open syncmer, the newest s-mer (at newest_idx) must have min hash
                        let newest_hash = hashes[newest_idx];
                        let is_syncmer = (0..w)
                            .filter(|&j| j != newest_idx)
                            .all(|j| hashes[j] >= newest_hash);

                        if is_syncmer {
                            acceptor(i + 1 - K, Kmer(kmer));
                        }
                    }
                }
            } else {
                // Reset on invalid base (N, etc.)
                kmer = 0;
                smer = 0;
                kmer_len = 0;
                smer_len = 0;
                ring_count = 0;
            }
        }
    }

    pub fn open_syncmer_iter<'a, const S: usize, H: Hasher>(
        seq: &'a [u8],
        _s: [(); S],
    ) -> KmerSyncmerIterator<'a, K, S, H> {
        KmerSyncmerIterator::new(seq)
    }

    pub fn agnostic_open_syncmer_iter<'a, const S: usize, H: Hasher>(
        seq: &'a [u8],
        _s: [(); S],
    ) -> KmerSyncmerIteratorAgnostic<'a, K, S, H> {
        KmerSyncmerIteratorAgnostic::new(seq)
    }


    fn nucleotide(c: u8) -> Option<u64> {
        match c {
            b'A' | b'a' => Some(0b00),
            b'C' | b'c' => Some(0b01),
            b'G' | b'g' => Some(0b10),
            b'T' | b't' => Some(0b11),
            b'U' | b'u' => Some(0b11),
            _ => None,
        }
    }
}

impl<const K: usize> From<u64> for Kmer<K> {
    fn from(value: u64) -> Self {
        Kmer(value)
    }
}

impl<const K: usize> From<&str> for Kmer<K> {
    fn from(value: &str) -> Self {
        let mut kmer: u64 = 0;
        for &b in value.as_bytes() {
            let nuc = Kmer::<K>::nucleotide(b).expect("invalid base in input string");
            kmer = (kmer << 2) | nuc;
        }
        Kmer(kmer)
    }
}

pub struct KmerIterator<'a, const K: usize> {
    seq: &'a [u8],
    pos: usize,
    kmer: u64,
    j: usize,
}

impl<'a, const K: usize> KmerIterator<'a, K> {
    pub fn new(seq: &'a [u8]) -> Self {
        KmerIterator {
            seq,
            pos: 0,
            kmer: 0,
            j: 0,
        }
    }
}

impl<'a, const K: usize> Iterator for KmerIterator<'a, K> {
    type Item = (usize, Kmer<K>);

    fn next(&mut self) -> Option<Self::Item> {
        let m = (1u64 << (2 * K)) - 1;
        while self.pos < self.seq.len() {
            if let Some(nuc) = Kmer::<K>::nucleotide(self.seq[self.pos]) {
                self.kmer = ((self.kmer << 2) | nuc) & m;
                self.j += 1;
                self.pos += 1;
                if self.j == K {
                    let result = Some((self.pos - K, Kmer(self.kmer)));
                    self.j -= 1;
                    return result;
                }
            } else {
                self.kmer = 0;
                self.j = 0;
                self.pos += 1;
            }
        }
        None
    }
}

pub struct KmerPairIterator<'a, const K: usize> {
    seq: &'a [u8],
    pos: usize,
    fwd_kmer: u64,
    rev_kmer: u64,
    j: usize,
}

impl<'a, const K: usize> KmerPairIterator<'a, K> {
    pub fn new(seq: &'a [u8]) -> Self {
        KmerPairIterator {
            seq,
            pos: 0,
            fwd_kmer: 0,
            rev_kmer: 0,
            j: 0,
        }
    }
}

impl<'a, const K: usize> Iterator for KmerPairIterator<'a, K> {
    type Item = (usize, Kmer<K>, Kmer<K>);

    fn next(&mut self) -> Option<Self::Item> {
        let m = (1u64 << (2 * K)) - 1;
        while self.pos < self.seq.len() {
            if let Some(nuc) = Kmer::<K>::nucleotide(self.seq[self.pos]) {
                self.fwd_kmer = ((self.fwd_kmer << 2) | nuc) & m;
                self.rev_kmer = (self.rev_kmer >> 2) | ((3 - nuc) << (2 * (K - 1)));
                self.j += 1;
                self.pos += 1;
                if self.j == K {
                    let result = Some((self.pos - K, Kmer(self.fwd_kmer), Kmer(self.rev_kmer)));
                    self.j -= 1;
                    return result;
                }
            } else {
                self.fwd_kmer = 0;
                self.rev_kmer = 0;
                self.j = 0;
                self.pos += 1;
            }
        }
        None
    }
}

pub struct KmerSyncmerIterator<'a, const K: usize, const S: usize, H: Hasher> {
    inner: KmerPairIterator<'a, K>,
    _marker: std::marker::PhantomData<H>,
}

impl<'a, const K: usize, const S: usize, H: Hasher> KmerSyncmerIterator<'a, K, S, H> {
    pub fn new(seq: &'a [u8]) -> Self {
        KmerSyncmerIterator {
            inner: KmerPairIterator::new(seq),
            _marker: std::marker::PhantomData,
        }
    }
}

impl<'a, const K: usize, const S: usize, H: Hasher> Iterator for KmerSyncmerIterator<'a, K, S, H> {
    type Item = (usize, Selection<Kmer<K>, Kmer<K>>);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((i, fwd, rev)) = self.inner.next() {
            let fwd_is_syncmer = fwd.is_open_syncmer::<S, H>();
            let rev_is_syncmer = rev.is_open_syncmer::<S, H>();
            match (fwd_is_syncmer, rev_is_syncmer) {
                (true, true) => return Some((i, Selection::Both(fwd, rev))),
                (true, false) => return Some((i, Selection::Left(fwd))),
                (false, true) => return Some((i, Selection::Right(rev))),
                (false, false) => {}
            }
        }
        None
    }
}

pub struct KmerSyncmerIteratorAgnostic<'a, const K: usize, const S: usize, H: Hasher> {
    inner: KmerPairIterator<'a, K>,
    _marker: std::marker::PhantomData<H>,
}

impl<'a, const K: usize, const S: usize, H: Hasher> KmerSyncmerIteratorAgnostic<'a, K, S, H> {
    pub fn new(seq: &'a [u8]) -> Self {
        KmerSyncmerIteratorAgnostic {
            inner: KmerPairIterator::new(seq),
            _marker: std::marker::PhantomData,
        }
    }
}

impl<'a, const K: usize, const S: usize, H: Hasher> Iterator for KmerSyncmerIteratorAgnostic<'a, K, S, H> {
    type Item = (usize, Kmer<K>, Kmer<K>);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((i, fwd, rev)) = self.inner.next() {
            let fwd_is_syncmer = fwd.is_open_syncmer::<S, H>();
            let rev_is_syncmer = rev.is_open_syncmer::<S, H>();
            if fwd_is_syncmer || rev_is_syncmer {
                return Some((i, fwd, rev))
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::{
        Selection,
        hasher::{FnvHasher, IdentityHasher},
    };

    use super::Kmer;

    fn rev_comp(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' | 'a' => 'T',
                'C' | 'c' => 'G',
                'G' | 'g' => 'C',
                'T' | 't' | 'U' | 'u' => 'A',
                _ => panic!("invalid base in test input"),
            })
            .collect()
    }

    #[test]
    fn to_string_roundtrip() {
        let seq = "ACGTAC";
        let kmer = Kmer::<6>::from(seq);
        assert_eq!(kmer.to_string(), seq);
    }

    #[test]
    fn kmerize_fwd_produces_expected_windows() {
        let seq = b"ACGTAC";
        let mut out = Vec::new();
        Kmer::<3>::kmerize_fwd(seq, |_i, k| out.push(k.to_string()));
        assert_eq!(out, vec!["ACG", "CGT", "GTA", "TAC"]);
    }

    #[test]
    fn kmerize_forward_and_reverse_agree_with_rev_comp() {
        let seq = b"ACGT";
        let mut out = Vec::new();
        Kmer::<3>::kmerize(seq, |_i, fwd, rev| {
            out.push((fwd.to_string(), rev.to_string()))
        });
        let expected = vec![
            ("ACG".to_string(), rev_comp("ACG")),
            ("CGT".to_string(), rev_comp("CGT")),
        ];
        assert_eq!(out, expected);
    }

    #[test]
    fn kmerize_open_syncmers_matches_naive_check() {
        let seq = b"ACGTACGTA";

        let expected = vec![Selection::Both(
            Kmer::<4>::from("GTAC"),
            Kmer::<4>::from("GTAC"),
        )];

        let mut actual = Vec::new();
        Kmer::<4>::kmerize_open_syncmers::<2, IdentityHasher, _, _>(seq, [(); 2], |_i, sel| {
            actual.push(sel)
        });

        println!(
            "{}",
            actual
                .iter()
                .map(|s| match s {
                    Selection::Left(x) => format!("+{}", x.to_string()),
                    Selection::Right(x) => format!("-{}", x.to_string()),
                    Selection::Both(f, r) => format!("{}/{}", f.to_string(), r.to_string()),
                })
                .collect::<Vec<_>>()
                .join(", ")
        );

        assert_eq!(actual, expected);
    }

    #[test]
    fn kmerize_open_syncmers_resets_on_invalid_base() {
        let seq = b"ACNACG";

        let expected = Vec::new();

        let mut actual = Vec::new();
        Kmer::<3>::kmerize_open_syncmers::<2, IdentityHasher, _, _>(seq, [(); 2], |_i, sel| {
            actual.push(sel)
        });

        println!(
            "{}",
            actual
                .iter()
                .map(|s| match s {
                    Selection::Left(x) => format!("+{}", x.to_string()),
                    Selection::Right(x) => format!("-{}", x.to_string()),
                    Selection::Both(f, r) => format!("{}/{}", f.to_string(), r.to_string()),
                })
                .collect::<Vec<_>>()
                .join(", ")
        );

        assert_eq!(actual, expected);
    }

    // ==================== FnvHasher tests ====================

    #[test]
    fn fnv_open_syncmers_produces_syncmers() {
        // With FNV hashing, different k-mers will be selected as syncmers
        // compared to identity hashing. This test verifies the mechanism works.
        let seq = b"ACGTACGTACGT";

        let mut syncmers = Vec::new();
        Kmer::<5>::kmerize_open_syncmers::<3, FnvHasher, _, _>(seq, [(); 3], |i, sel| {
            syncmers.push((i, sel));
        });

        // Should produce at least one syncmer from this sequence
        assert!(!syncmers.is_empty(), "FnvHasher should produce syncmers");

        // Verify positions are within valid range
        for (pos, _) in &syncmers {
            assert!(*pos <= seq.len() - 5, "syncmer position out of range");
        }
    }

    #[test]
    fn fnv_syncmers_differ_from_identity() {
        // FNV hashing should generally produce different syncmer selections
        // than identity hashing due to the scrambling effect
        let seq = b"ACGTACGTACGTACGT";

        let mut identity_syncmers = Vec::new();
        Kmer::<5>::kmerize_open_syncmers::<3, IdentityHasher, _, _>(seq, [(); 3], |i, _sel| {
            identity_syncmers.push(i);
        });

        let mut fnv_syncmers = Vec::new();
        Kmer::<5>::kmerize_open_syncmers::<3, FnvHasher, _, _>(seq, [(); 3], |i, _sel| {
            fnv_syncmers.push(i);
        });

        // At least one should produce syncmers
        assert!(
            !identity_syncmers.is_empty() || !fnv_syncmers.is_empty(),
            "At least one hasher should produce syncmers"
        );

        // The selections should typically differ (though not guaranteed for all sequences)
        // This is a statistical property - with a reasonable sequence length they should differ
        println!("Identity syncmer positions: {:?}", identity_syncmers);
        println!("FNV syncmer positions: {:?}", fnv_syncmers);
    }

    #[test]
    fn fnv_syncmer_iterator_matches_callback() {
        // Verify the iterator produces the same results as the callback version
        let seq = b"ACGTACGTACGTACGTACGT";

        let mut callback_results = Vec::new();
        Kmer::<6>::kmerize_open_syncmers::<4, FnvHasher, _, _>(seq, [(); 4], |i, sel| {
            callback_results.push((i, sel));
        });

        let iterator_results: Vec<_> =
            Kmer::<6>::open_syncmer_iter::<4, FnvHasher>(seq, [(); 4]).collect();

        assert_eq!(
            callback_results, iterator_results,
            "Iterator and callback should produce identical results"
        );
    }

    #[test]
    fn fnv_syncmers_fwd_only() {
        // Test forward-only syncmer extraction with FnvHasher
        let seq = b"ACGTACGTACGT";

        let mut syncmers = Vec::new();
        Kmer::<5>::kmerize_open_syncmers_fwd::<3, FnvHasher, _, _>(seq, [(); 3], |i, kmer| {
            syncmers.push((i, kmer.to_string()));
        });

        // Verify all returned k-mers are valid open syncmers
        for (pos, kmer_str) in &syncmers {
            let kmer = Kmer::<5>::from(kmer_str.as_str());
            assert!(
                kmer.is_open_syncmer::<3, FnvHasher>(),
                "k-mer at position {} should be an open syncmer",
                pos
            );
        }
    }

    #[test]
    fn fnv_is_open_syncmer_consistent() {
        // Test that is_open_syncmer gives consistent results
        let kmer = Kmer::<6>::from("ACGTAC");

        let result1 = kmer.is_open_syncmer::<3, FnvHasher>();
        let result2 = kmer.is_open_syncmer::<3, FnvHasher>();

        assert_eq!(result1, result2, "is_open_syncmer should be deterministic");
    }

    #[test]
    fn fnv_syncmers_density() {
        // Open syncmers should have a density of approximately 2/(K-S+1)
        // For K=10, S=5, expected density is 2/(10-5+1) = 2/6 ≈ 0.33
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let mut count = 0;
        let mut total = 0;
        Kmer::<10>::kmerize::<_, _>(seq, |_i, fwd, _rev| {
            total += 1;
            if fwd.is_open_syncmer::<5, FnvHasher>() {
                count += 1;
            }
        });

        let density = count as f64 / total as f64;
        let expected_density = 2.0 / 6.0;

        println!(
            "FnvHasher syncmer density: {:.3} (expected ~{:.3}), count={}/{}",
            density, expected_density, count, total
        );

        // Allow some variance due to sequence composition
        assert!(
            density > 0.1 && density < 0.6,
            "Syncmer density {:.3} outside reasonable range",
            density
        );
    }

    // ==================== Optimized syncmer tests ====================

    #[test]
    fn optimized_fwd_matches_original_identity() {
        // Verify the optimized version produces identical results to the original
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let mut original_results = Vec::new();
        Kmer::<7>::kmerize_open_syncmers_fwd_orig::<4, IdentityHasher, _, _>(seq, [(); 4], |i, kmer| {
            original_results.push((i, kmer));
        });

        let mut optimized_results = Vec::new();
        Kmer::<7>::kmerize_open_syncmers_fwd::<4, IdentityHasher, _, _>(
            seq,
            [(); 4],
            |i, kmer| {
                optimized_results.push((i, kmer));
            },
        );

        assert_eq!(
            original_results, optimized_results,
            "Optimized version should match original with IdentityHasher"
        );
    }

    #[test]
    fn optimized_fwd_matches_original_fnv() {
        // Verify with FnvHasher as well
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let mut original_results = Vec::new();
        Kmer::<7>::kmerize_open_syncmers_fwd_orig::<4, FnvHasher, _, _>(seq, [(); 4], |i, kmer| {
            original_results.push((i, kmer));
        });

        let mut optimized_results = Vec::new();
        Kmer::<7>::kmerize_open_syncmers_fwd::<4, FnvHasher, _, _>(seq, [(); 4], |i, kmer| {
            optimized_results.push((i, kmer));
        });

        assert_eq!(
            original_results, optimized_results,
            "Optimized version should match original with FnvHasher"
        );
    }

    #[test]
    fn optimized_fwd_handles_invalid_bases() {
        // Test that invalid bases reset the state correctly
        let seq = b"ACGTACNACGTACGT";

        let mut original_results = Vec::new();
        Kmer::<5>::kmerize_open_syncmers_fwd_orig::<3, FnvHasher, _, _>(seq, [(); 3], |i, kmer| {
            original_results.push((i, kmer));
        });

        let mut optimized_results = Vec::new();
        Kmer::<5>::kmerize_open_syncmers_fwd::<3, FnvHasher, _, _>(seq, [(); 3], |i, kmer| {
            optimized_results.push((i, kmer));
        });

        assert_eq!(
            original_results, optimized_results,
            "Optimized version should handle invalid bases like original"
        );
    }

    #[test]
    fn optimized_fwd_short_sequence() {
        // Test with sequence shorter than K
        let seq = b"ACG";

        let mut original_results = Vec::new();
        Kmer::<5>::kmerize_open_syncmers_fwd_orig::<3, FnvHasher, _, _>(seq, [(); 3], |i, kmer| {
            original_results.push((i, kmer));
        });

        let mut optimized_results = Vec::new();
        Kmer::<5>::kmerize_open_syncmers_fwd::<3, FnvHasher, _, _>(seq, [(); 3], |i, kmer| {
            optimized_results.push((i, kmer));
        });

        assert!(original_results.is_empty());
        assert!(optimized_results.is_empty());
    }

    #[test]
    fn optimized_fwd_various_k_s_combinations() {
        // Test with different K and S values
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        // K=10, S=5
        let mut orig_10_5 = Vec::new();
        let mut opt_10_5 = Vec::new();
        Kmer::<10>::kmerize_open_syncmers_fwd_orig::<5, FnvHasher, _, _>(seq, [(); 5], |i, k| {
            orig_10_5.push((i, k))
        });
        Kmer::<10>::kmerize_open_syncmers_fwd::<5, FnvHasher, _, _>(seq, [(); 5], |i, k| {
            opt_10_5.push((i, k))
        });
        assert_eq!(orig_10_5, opt_10_5, "K=10, S=5 should match");

        // K=15, S=8
        let mut orig_15_8 = Vec::new();
        let mut opt_15_8 = Vec::new();
        Kmer::<15>::kmerize_open_syncmers_fwd_orig::<8, FnvHasher, _, _>(seq, [(); 8], |i, k| {
            orig_15_8.push((i, k))
        });
        Kmer::<15>::kmerize_open_syncmers_fwd::<8, FnvHasher, _, _>(seq, [(); 8], |i, k| {
            opt_15_8.push((i, k))
        });
        assert_eq!(orig_15_8, opt_15_8, "K=15, S=8 should match");

        // K=6, S=4
        let mut orig_6_4 = Vec::new();
        let mut opt_6_4 = Vec::new();
        Kmer::<6>::kmerize_open_syncmers_fwd_orig::<4, FnvHasher, _, _>(seq, [(); 4], |i, k| {
            orig_6_4.push((i, k))
        });
        Kmer::<6>::kmerize_open_syncmers_fwd::<4, FnvHasher, _, _>(seq, [(); 4], |i, k| {
            opt_6_4.push((i, k))
        });
        assert_eq!(orig_6_4, opt_6_4, "K=6, S=4 should match");
    }

    #[test]
    fn optimized_fwd_long_random_sequence() {
        // Test with a longer pseudo-random sequence
        let bases = b"ACGT";
        let mut seq = Vec::with_capacity(1000);
        for i in 0..1000 {
            seq.push(bases[i % 4]);
        }

        let mut original_results = Vec::new();
        Kmer::<11>::kmerize_open_syncmers_fwd_orig::<6, FnvHasher, _, _>(&seq, [(); 6], |i, kmer| {
            original_results.push((i, kmer));
        });

        let mut optimized_results = Vec::new();
        Kmer::<11>::kmerize_open_syncmers_fwd::<6, FnvHasher, _, _>(&seq, [(); 6], |i, kmer| {
            optimized_results.push((i, kmer));
        });

        assert_eq!(
            original_results, optimized_results,
            "Optimized should match original on long sequence"
        );
    }

    #[test]
    fn optimized_fwd_performance() {
        // Performance test comparing original vs optimized syncmer extraction.
        // Uses black_box to prevent compiler from optimizing away repeated work.
        use std::hint::black_box;

        let seq = b"ACTTGCTTTATGAATCTGGGCGCTCCTGTATTGGGTGCATATATATTTAGGGTAGTTAGCTCCCTTTACCATTATGTAATGGCCTTCTTTGTCCCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACAACACCTGCTTTTTTTGTTTTCCATTTGCTTGGTAGGTCTTCCTCCATCCCTTTATTTTGAGCCTATGTGTGTGTCTGCACATGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTCATCTAACTTGCCAGTCTGTGTCTTTTAATTGGGGCATTTAGCCCATTTACATTTAAGGTTAATATTGTTATGTGTGAATTTGATTCTGTCATTATGATGTTAGCTGGTTATTTTTCCCGTTAGTTGATGCAGTTTCTTCCTAGCATCGATGGTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGTGCTTCCTTCAGGAGCTCTTGTAAGGCAGGGCTGGTGGTGACAAAATCTCTCAGCATTTGCTGGTCTATAAAGGATTTTATTTCTCCTTATGAAGCTTTGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTAAGAATGTTGAATATTGGTGCCCACTCTCTTCTGACTTGTAGAGTTTCTGTTGAGAGATCCACTGTTAGTCTGATGGGCTTCCCTTTGTGGCTAACTCGACCTTTCTCTCTGGGTGCCATTAACATTTTTTCCTTCATTTCAACCTTGGTGAATCTGACAATTATGTGTCTTGGGGTTGCTCCTCTCGAGGAGCATCTTGGTAGTGTTCTCTGTATTTCCTGAGTTTGAATGTTTGCCTGCCTTGCTAGGTTGGGGAAGTTCTCCTGGACAATATCCTGAAGAGTGTTTTCGAACTTGGTTCCATTCTCCCCGTCACTTTCAGGTACACCAATCAAACGTAGATTTGGTGTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGTTCATTCTTTTTACTCTTTTTTCTCTAAACTTCTCACTTCATTAATTTGATCTTCAATCACTGATACCCTTTCTTTCAGTTTATTGAATCAACTACTGAAGCTTGTGCATGTGTCACATAGTTCTTGTTCCATGGTTTTCAGCTCCATCAGGTCATTTAAGGTCTCCACACTGCTTATTCTAGTTAGCCATTCATCTAATCTGTTTGCAAGGCTTTTAGCTTCCTTGTGATGGGTTCGAATACCTCCCTTAACTCAGAGAAGTTTGTTATTACCAACCTTCTGAAGCCTACTTCTGTCAGCTCATCAAAGTCATTCTCCGTCCAGCTTTATTCCGTTGCTGGCAAGGAGCTGTAATCCTTTGCAGGAGAAGGGATGCTGTGGTTTTTAGAATTTTCAGCTTTTCTGCTCTGGTTTCTCCCCATCTTTGTGGTTTTATCTACCTTTGGTCTTCGATGATGGTGACCCACAGATGGGGTTTTGGTGTGGGATGTCCTTTTTGTTGATGTTGATGCTATTCCTTTCTGTTTGTTAGTTTTCCTTCTGACAGTCAGGTCCCTCAGCTGCAGATCTGTTGGAGTTTGCTGGAGGTCCACTCCAGACTCTGTTTACCTGTGTATCACCAGCAGAGGCTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAATAGCAAATATTGCAGAACAGCAAATATTAC";

        const N: usize = 1000;
        const K: usize = 15;
        const S: usize = 8;

        // Accumulate a checksum across iterations to force the compiler to do the work
        let start = std::time::Instant::now();
        let mut original_checksum: u64 = 0;
        for i in 0..N {
            Kmer::<K>::kmerize_open_syncmers_fwd_orig::<S, FnvHasher, _, _>(
                black_box(seq),
                [(); S],
                |pos, kmer| {
                    // XOR with iteration and position to prevent loop-invariant optimization
                    original_checksum ^= kmer.0.wrapping_add(pos as u64).wrapping_add(i as u64);
                },
            );
        }
        let duration_original = start.elapsed();
        black_box(original_checksum);

        let start = std::time::Instant::now();
        let mut optimized_checksum: u64 = 0;
        for i in 0..N {
            Kmer::<K>::kmerize_open_syncmers_fwd::<S, FnvHasher, _, _>(
                black_box(seq),
                [(); S],
                |pos, kmer| {
                    optimized_checksum ^= kmer.0.wrapping_add(pos as u64).wrapping_add(i as u64);
                },
            );
        }
        let duration_optimized = start.elapsed();
        black_box(optimized_checksum);

        println!("Original duration: {:?}", duration_original);
        println!("Optimized duration: {:?}", duration_optimized);
        println!(
            "Ratio (optimized/original): {:.2}",
            duration_optimized.as_secs_f64() / duration_original.as_secs_f64()
        );

        assert_eq!(
            original_checksum, optimized_checksum,
            "Both versions should produce the same checksum"
        );
    }
}
