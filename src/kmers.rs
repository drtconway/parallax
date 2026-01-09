#![allow(dead_code)]

use crate::utils::Selection;

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

    pub fn is_open_syncmer<const S: usize>(&self) -> bool {
        let mask: u64 = (1 << (2 * S)) - 1;
        let mut min_subkmer = u64::MAX;
        for i in 0..=(K - S) {
            let subkmer = (self.0 >> (2 * i)) & mask;
            if subkmer < min_subkmer {
                min_subkmer = subkmer;
            }
        }
        let first_subkmer = self.0 & mask;
        first_subkmer == min_subkmer
    }

    pub fn kmerize_open_syncmers<
        const S: usize,
        Seq: AsRef<[u8]>,
        F: FnMut(usize, Selection<Kmer<K>, Kmer<K>>),
    >(
        seq: Seq,
        _s: [(); S],
        mut acceptor: F,
    ) {
        Kmer::<K>::kmerize(seq, |i, fwd, rev| {
            match (fwd.is_open_syncmer::<S>(), rev.is_open_syncmer::<S>()) {
                (true, true) => acceptor(i, Selection::Both(fwd, rev)),
                (true, false) => acceptor(i, Selection::Left(fwd)),
                (false, true) => acceptor(i, Selection::Right(rev)),
                (false, false) => {}
            }
        });
    }

    pub fn open_syncmer_iter<'a, const S: usize>(
        seq: &'a [u8],
        _s: [(); S],
    ) -> KmerSyncmerIterator<'a, K, S> {
        KmerSyncmerIterator::new(seq)
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

pub struct KmerSyncmerIterator<'a, const K: usize, const S: usize> {
    inner: KmerPairIterator<'a, K>,
}

impl<'a, const K: usize, const S: usize> KmerSyncmerIterator<'a, K, S> {
    pub fn new(seq: &'a [u8]) -> Self {
        KmerSyncmerIterator {
            inner: KmerPairIterator::new(seq),
        }
    }
}

impl<'a, const K: usize, const S: usize> Iterator for KmerSyncmerIterator<'a, K, S> {
    type Item = (usize, Selection<Kmer<K>, Kmer<K>>);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((i, fwd, rev)) = self.inner.next() {
            let fwd_is_syncmer = fwd.is_open_syncmer::<S>();
            let rev_is_syncmer = rev.is_open_syncmer::<S>();
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

#[cfg(test)]
mod tests {
    use crate::utils::Selection;

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
        Kmer::<4>::kmerize_open_syncmers(seq, [(); 2], |_i, sel| actual.push(sel));

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
        Kmer::<3>::kmerize_open_syncmers(seq, [(); 2], |_i, sel| actual.push(sel));

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
}
