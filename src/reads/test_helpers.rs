use std::io::{BufRead, BufReader};
use crate::reads::compound::AtomicSeed;

pub struct RawSeedRow {
    pub read_id: String,
    pub strand: bool,
    pub kmer: String,
    pub chrom: String,
    pub ref_pos: u32,
    pub read_pos: u32,
    pub kmer_multiplicity: u32,
}

pub fn load_seeds(path: &str) -> Vec<RawSeedRow> {
    let file = std::fs::File::open(path).unwrap_or_else(|_| panic!("cannot open {path}"));
    let (reader, _) = niffler::get_reader(Box::new(file)).expect("niffler open");
    let mut lines = BufReader::new(reader).lines();
    lines.next(); // skip header
    let mut rows = Vec::new();
    for line in lines {
        let line = line.expect("line read");
        let cols: Vec<&str> = line.splitn(7, '\t').collect();
        assert_eq!(cols.len(), 7, "expected 7 columns, got {}: {line}", cols.len());
        rows.push(RawSeedRow {
            read_id: cols[0].to_owned(),
            strand: cols[1] == "-",
            kmer: cols[2].to_owned(),
            chrom: cols[3].to_owned(),
            ref_pos: cols[4].parse().unwrap(),
            read_pos: cols[5].parse().unwrap(),
            kmer_multiplicity: cols[6].parse().unwrap(),
        });
    }
    rows
}

pub fn rows_to_atomic(rows: &[RawSeedRow], k: usize) -> (Vec<AtomicSeed>, Vec<String>) {
    let read_len = rows.iter().map(|r| r.read_pos + k as u32).max().unwrap_or(0);
    let mut chrom_names: Vec<String> = Vec::new();
    let mut chrom_id_of = |name: &str| -> u32 {
        if let Some(i) = chrom_names.iter().position(|n| n == name) {
            return i as u32;
        }
        let i = chrom_names.len() as u32;
        chrom_names.push(name.to_owned());
        i
    };
    let atoms: Vec<AtomicSeed> = rows.iter().map(|r| {
        let chrom_id = chrom_id_of(&r.chrom);
        let kmer_u64 = kmer_str_to_u64(&r.kmer, k);
        AtomicSeed::new(r.read_pos, read_len, k, chrom_id, r.ref_pos, r.strand, kmer_u64, r.kmer_multiplicity)
    }).collect();
    (atoms, chrom_names)
}

pub fn kmer_str_to_u64(s: &str, k: usize) -> u64 {
    let mut v: u64 = 0;
    for (i, b) in s.bytes().enumerate() {
        if i >= k { break; }
        let bits: u64 = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        };
        v = (v << 2) | bits;
    }
    v
}
