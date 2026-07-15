use std::io::{BufRead, BufReader};

use crate::reads::compound::{
    chain_seeds, AtomicSeed, ChainingDPScheme, ChainResult, CompoundSeed, DPConfig, EdgeType,
    FullDPScheme, GapComputable, Seed, SeedCollection, Weighted,
};

// ── Fixture loading ────────────────────────────────────────────────────────────

struct RawSeedRow {
    _read_id: String,
    strand: bool, // true = reverse
    kmer: String,
    chrom: String,
    ref_pos: u32,
    read_pos: u32,
    kmer_multiplicity: u32,
}

fn load_seeds(path: &str) -> Vec<RawSeedRow> {
    let file = std::fs::File::open(path).unwrap_or_else(|_| panic!("cannot open {path}"));
    let (reader, _) = niffler::get_reader(Box::new(file)).expect("niffler open");
    let mut lines = BufReader::new(reader).lines();

    // skip header
    lines.next();

    let mut rows = Vec::new();
    for line in lines {
        let line = line.expect("line read");
        let cols: Vec<&str> = line.splitn(7, '\t').collect();
        assert_eq!(
            cols.len(),
            7,
            "expected 7 columns, got {}: {line}",
            cols.len()
        );
        rows.push(RawSeedRow {
            _read_id: cols[0].to_owned(),
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

/// Build a fake chrom_id mapping from whatever chroms appear in the fixture.
/// Returns (rows_as_AtomicSeed, chrom_names).
fn rows_to_atomic(rows: &[RawSeedRow], k: usize) -> (Vec<AtomicSeed>, Vec<String>) {
    // Stable chrom ordering in first-appearance order.
    let mut chrom_names: Vec<String> = Vec::new();
    let mut chrom_id_of = |name: &str| -> u32 {
        if let Some(i) = chrom_names.iter().position(|n| n == name) {
            return i as u32;
        }
        let i = chrom_names.len() as u32;
        chrom_names.push(name.to_owned());
        i
    };

    let atoms: Vec<AtomicSeed> = rows
        .iter()
        .map(|r| {
            let chrom_id = chrom_id_of(&r.chrom);
            // kmer string → u64 by treating first 8 bytes as big-endian (just for identity)
            let kmer_u64 = kmer_str_to_u64(&r.kmer, k);
            AtomicSeed::new(
                r.read_pos,
                chrom_id,
                r.ref_pos,
                r.strand,
                kmer_u64,
                r.kmer_multiplicity,
            )
        })
        .collect();

    (atoms, chrom_names)
}

/// Encode a DNA kmer string as a 2-bit u64 (A=0, C=1, G=2, T=3).
/// Truncates or zero-pads to k bases.
fn kmer_str_to_u64(s: &str, k: usize) -> u64 {
    let mut v: u64 = 0;
    for (i, b) in s.bytes().enumerate() {
        if i >= k {
            break;
        }
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

fn fixture_path() -> &'static str {
    "tests/data/SRR29147690.1001-seeds.tsv.gz"
}

// ── Unit tests for AtomicSeed::gap_to ─────────────────────────────────────────

fn atomic(read_pos: u32, ref_pos: u32, is_reverse: bool, mult: u32) -> AtomicSeed {
    AtomicSeed::new(read_pos, 0, ref_pos, is_reverse, 0, mult)
}

#[test]
fn atomic_gap_no_overlap_fwd() {
    let k = 5;
    // lhs ends at read 10, rhs starts at read 15 → read_gap = 5
    // lhs ref ends at ref 110, rhs ref starts at ref 115 → ref_gap = 5
    let lhs = atomic(5, 100, false, 1);
    let rhs = atomic(15, 115, false, 1);
    let gap = lhs.gap_to(&rhs, k).unwrap();
    assert_eq!(gap.read_gap, 5);
    assert_eq!(gap.ref_gap, 10); // lhs ref end=105, rhs ref start=115 → 10
    assert_eq!(gap.weight_trimmed, 0.0);
}

#[test]
fn atomic_gap_no_overlap_rev() {
    let k = 5;
    // Reverse strand: ref_pos increases with read_pos (5' end in read-traversal direction).
    // lhs read=5, ref=100; rhs read=15, ref=110
    // read_gap = 15 - (5+5) = 5; ref_gap = 110 - (100+5) = 5
    let lhs = atomic(5, 100, true, 1);
    let rhs = atomic(15, 110, true, 1);
    let gap = lhs.gap_to(&rhs, k).unwrap();
    assert_eq!(gap.read_gap, 5);
    assert_eq!(gap.ref_gap, 5);
    assert_eq!(gap.weight_trimmed, 0.0);
}

#[test]
fn atomic_gap_cross_strand_is_sv_break() {
    let k = 5;
    // lhs read=0, rhs read=10: read_gap = 10 - (0 + 5) = 5
    let lhs = atomic(0, 100, false, 1);
    let rhs = atomic(10, 110, true, 1);
    let gap = lhs.gap_to(&rhs, k).unwrap();
    assert_eq!(gap.ref_gap, i64::MIN);
    assert_eq!(gap.read_gap, 5);
}

#[test]
fn atomic_gap_partial_overlap_same_diagonal() {
    let k = 5;
    // lhs at read=0, ref=100; rhs at read=3, ref=103 — same diagonal (100), 2-base overlap
    // read_gap = 3 - (0 + 5) = -2; ref_gap = 103 - (100 + 5) = -2
    let lhs = atomic(0, 100, false, 1);
    let rhs = atomic(3, 103, false, 1);
    let gap = lhs.gap_to(&rhs, k).unwrap();
    assert_eq!(gap.read_gap, -2);
    assert_eq!(gap.ref_gap, -2);
    assert_eq!(gap.weight_trimmed, 0.0);
}

#[test]
fn atomic_gap_partial_overlap_different_diagonal_is_none() {
    let k = 5;
    // Different diagonals with overlap → None (no valid split between two k-mers)
    let lhs = atomic(0, 100, false, 1);
    let rhs = atomic(3, 110, false, 1);
    assert!(lhs.gap_to(&rhs, k).is_none());
}

#[test]
fn atomic_gap_fully_consumed_is_none() {
    let k = 5;
    // rhs starts inside lhs and ends inside lhs → fully consumed
    let lhs = atomic(0, 100, false, 1);
    let rhs = atomic(1, 101, false, 1); // overlap = 4 bases >= k=5? No, overlap = 4 < 5
    // overlap = k - (rhs.read_pos - lhs.read_pos) = 5 - 1 = 4 < k → not fully consumed
    assert!(lhs.gap_to(&rhs, k).is_some());

    // Now make rhs fully consumed: rhs.read_pos = 0, same as lhs — overlap = k
    let rhs2 = atomic(0, 100, false, 1);
    assert!(lhs.gap_to(&rhs2, k).is_none());
}

// ── Unit tests for AtomicSeed weight ──────────────────────────────────────────

#[test]
fn atom_weight_unique() {
    let a = atomic(0, 100, false, 1);
    assert!((a.weight() - 1.0).abs() < 1e-9);
}

#[test]
fn atom_weight_mult4() {
    let a = atomic(0, 100, false, 4);
    assert!((a.weight() - 0.5).abs() < 1e-9);
}

// ── Unit tests for CompoundSeed weight ────────────────────────────────────────

fn make_compound(atoms: &[AtomicSeed]) -> CompoundSeed<'_> {
    CompoundSeed::new(atoms)
}

#[test]
fn compound_weight_is_sum_of_atoms() {
    let atoms = vec![
        atomic(0, 100, false, 1),
        atomic(5, 105, false, 4),
        atomic(10, 110, false, 9),
    ];
    let cs = make_compound(&atoms);
    let expected = 1.0 + 0.5 + (1.0 / 3.0);
    assert!((cs.weight() - expected).abs() < 1e-9);
}

// ── Unit tests for SeedCollection::compound_seeds ────────────────────────────

#[test]
fn compound_seeds_merges_overlapping_on_same_diagonal() {
    let k = 5;
    // Three atoms on the same diagonal (ref - read = 100), each overlapping the next.
    let seeds = vec![
        atomic(0, 100, false, 1), // read 0..5,  ref 100..105
        atomic(3, 103, false, 1), // read 3..8,  ref 103..108  — overlaps prev by 2
        atomic(6, 106, false, 1), // read 6..11, ref 106..111  — overlaps prev by 2
    ];
    let col = SeedCollection::new(k, seeds);
    let compounds = col.compound_seeds();
    assert_eq!(
        compounds.len(),
        1,
        "all three should merge into one compound seed"
    );
    assert_eq!(compounds[0].atoms().len(), 3);
}

#[test]
fn compound_seeds_splits_on_positive_read_gap() {
    let k = 5;
    // Two atoms on the same diagonal but with a gap between them.
    let seeds = vec![
        atomic(0, 100, false, 1),  // read 0..5
        atomic(10, 110, false, 1), // read 10..15 — gap of 5, same diagonal
    ];
    let col = SeedCollection::new(k, seeds);
    let compounds = col.compound_seeds();
    assert_eq!(
        compounds.len(),
        2,
        "gap should split into two compound seeds"
    );
}

#[test]
fn compound_seeds_splits_on_different_diagonal() {
    let k = 5;
    // Two atoms that abut in read space but are on different diagonals.
    let seeds = vec![
        atomic(0, 100, false, 1), // diagonal = 100
        atomic(5, 110, false, 1), // diagonal = 105, read_gap = 0 but ref_gap = 5 ≠ read_gap
    ];
    let col = SeedCollection::new(k, seeds);
    let compounds = col.compound_seeds();
    assert_eq!(compounds.len(), 2, "different diagonals should not merge");
}

#[test]
fn compound_seeds_splits_on_cross_strand() {
    let k = 5;
    let seeds = vec![atomic(0, 100, false, 1), atomic(5, 110, true, 1)];
    let col = SeedCollection::new(k, seeds);
    let compounds = col.compound_seeds();
    assert_eq!(compounds.len(), 2);
}

#[test]
fn compound_seeds_singleton_when_no_merges() {
    let k = 5;
    let seeds = vec![atomic(0, 100, false, 1)];
    let col = SeedCollection::new(k, seeds);
    let compounds = col.compound_seeds();
    assert_eq!(compounds.len(), 1);
    assert_eq!(compounds[0].atoms().len(), 1);
}

#[test]
fn compound_seeds_empty_collection() {
    let col = SeedCollection::new(5, vec![]);
    assert!(col.compound_seeds().is_empty());
}

// ── Unit tests for CompoundSeed::gap_to overlap trimming ──────────────────────

#[test]
fn compound_gap_no_overlap() {
    let k = 5;
    // lhs: atoms at read 0, 5; rhs: atoms at read 15, 20
    let lhs_atoms = vec![atomic(0, 100, false, 1), atomic(5, 105, false, 1)];
    let rhs_atoms = vec![atomic(15, 115, false, 1), atomic(20, 120, false, 1)];
    let lhs = make_compound(&lhs_atoms);
    let rhs = make_compound(&rhs_atoms);
    // lhs ends at read 10 (last atom read_pos=5, +k=5 → 10); rhs starts at 15 → gap=5
    let gap = lhs.gap_to(&rhs, k).unwrap();
    assert_eq!(gap.read_gap, 5);
    assert_eq!(gap.ref_gap, 5);
    assert_eq!(gap.weight_trimmed, 0.0);
}

#[test]
fn compound_gap_overlap_trims_lighter_atoms() {
    let k = 5;
    // lhs: two atoms; last one (read=8, ref=108) is unique (mult=1, weight=1.0)
    // rhs: first atom (read=10, ref=110) overlaps, also unique (mult=1)
    // Overlap region: read 10..13 (lhs ends at 13, rhs starts at 10)
    // lhs_overlap atom: read_pos=8, read_pos+k=13 > rhs.read_pos=10 ✓
    // rhs_overlap atom: read_pos=10 < lhs_end=13 ✓
    // Candidate splits: {13 (lhs atom end), 10 (rhs atom start)}
    // At split=10: lhs_trimmed=1.0 (atom at 8 ends at 13 > 10), rhs_trimmed=0.0 → total=1.0
    // At split=13: lhs_trimmed=0.0 (no lhs atom ends after 13), rhs_trimmed=1.0 (atom at 10 < 13) → total=1.0
    // Both equal — picks first encountered (split=10)
    let lhs_atoms = vec![atomic(0, 100, false, 4), atomic(8, 108, false, 1)];
    let rhs_atoms = vec![atomic(10, 110, false, 1), atomic(18, 118, false, 4)];
    let lhs = make_compound(&lhs_atoms);
    let rhs = make_compound(&rhs_atoms);
    let gap = lhs.gap_to(&rhs, k).unwrap();
    // Either split is equally good (weight_trimmed = 1.0)
    assert!((gap.weight_trimmed - 1.0).abs() < 1e-9);
}

#[test]
fn compound_gap_overlap_prefers_trimming_repetitive_atom() {
    let k = 5;
    // lhs last atom: read=8, mult=16 (weight=0.25) — repetitive
    // rhs first atom: read=10, mult=1 (weight=1.0) — unique
    // At split=13 (lhs atom end retained, rhs atom trimmed): weight_trimmed = 1.0
    // At split=10 (lhs atom trimmed, rhs atom retained): weight_trimmed = 0.25
    // Optimal: trim the lhs repetitive atom → split at 10
    let lhs_atoms = vec![atomic(0, 100, false, 1), atomic(8, 108, false, 16)];
    let rhs_atoms = vec![atomic(10, 110, false, 1), atomic(18, 118, false, 1)];
    let lhs = make_compound(&lhs_atoms);
    let rhs = make_compound(&rhs_atoms);
    let gap = lhs.gap_to(&rhs, k).unwrap();
    assert!((gap.weight_trimmed - 0.25).abs() < 1e-9);
}

// ── Integration test: load fixture, build SeedCollection, prune ───────────────

#[test]
fn prune_isolated_seeds_chr13_primary_survives() {
    let path = fixture_path();
    if !std::path::Path::new(path).exists() {
        eprintln!("fixture not found at {path}, skipping");
        return;
    }

    let k = 20;
    let rows = load_seeds(path);

    // Keep only chr13 seeds.
    let chr13: Vec<RawSeedRow> = rows
        .into_iter()
        .filter(|r| r.chrom == "chr13")
        .collect();
    let (atoms, _) = rows_to_atomic(&chr13, k);

    let collection = SeedCollection::new(k, atoms.clone());
    let before = collection.hits.len();

    let scheme = FullDPScheme::new(DPConfig::default());
    let pruned = collection.prune_isolated_seeds(&scheme);

    // The primary diagonal is densely covered — expect most seeds to survive.
    assert!(
        pruned.hits.len() > before / 2,
        "too many chr13 fwd seeds pruned: {before} → {}",
        pruned.hits.len()
    );

    let seeds = pruned.compound_seeds();

    if false {
        for seed in seeds.iter() {
            let strand = if seed.is_reverse() { "-" } else { "+" };
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                seed.diagonal(),
                seed.read_start(),
                seed.read_end(k),
                seed.ref_start(),
                seed.ref_end(k),
                strand,
                seed.weight(),
                seed.atoms().len(),
                seed.length(k)
            );
        }
    }

    let scheme = FullDPScheme::new(DPConfig::default());
    let ChainResult { score, chain, edge_types } = chain_seeds(&seeds, k, &scheme).unwrap();

    if false {
        println!("Chained {} compound seeds with score {score:.1}", chain.len());
        for (i, &idx) in chain.iter().enumerate() {
            let seed = &seeds[idx];
            let edge_type = if i == 0 { String::from("START") } else { format!("{:?}", edge_types[i - 1]) };
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}",
                i,
                seed.read_start(),
                seed.read_end(k),
                seed.ref_start(),
                seed.ref_end(k),
                if seed.is_reverse() { "-" } else { "+" },
                seed.weight(),
                seed.atoms().len(),
                seed.length(k),
                edge_type
            );
        }
    }

    //assert!(false);
}

#[test]
fn prune_isolated_seeds_removes_scattered_noise() {
    // Build a small synthetic set: one pair of colinear seeds and one isolated seed
    // on a totally different chromosome. The isolated seed should be pruned.
    let k = 5;
    let seeds = vec![
        // Colinear pair on chrom 0
        AtomicSeed::new(0, 0, 100, false, 0, 1),
        AtomicSeed::new(10, 0, 110, false, 0, 1),
        // Isolated seed on chrom 1 — no neighbour within max_neighbour_gap
        AtomicSeed::new(5, 1, 200, false, 0, 1),
    ];

    let collection = SeedCollection::new(k, seeds);
    let scheme = FullDPScheme::new(DPConfig::default());
    let pruned = collection.prune_isolated_seeds(&scheme);

    assert_eq!(
        pruned.hits.len(),
        2,
        "isolated chrom-1 seed should be pruned"
    );
    assert!(pruned.hits.iter().all(|s| s.chrom_id == 0));
}

// ── FullDPScheme edge_penalty ─────────────────────────────────────────────────

#[test]
fn edge_penalty_continuation() {
    let k = 5;
    let scheme = FullDPScheme::new(DPConfig::default());
    let lhs = atomic(0, 100, false, 1);
    let rhs = atomic(10, 110, false, 1); // read_gap=5, ref_gap=5, deviation=0
    let (penalty, edge_type) = scheme.edge_penalty(&lhs, &rhs, k).unwrap();
    assert_eq!(edge_type, EdgeType::Continuation);
    assert!(
        penalty < 200.0,
        "continuation should be cheap, got {penalty}"
    );
}

#[test]
fn edge_penalty_sv_break_cross_strand() {
    let k = 5;
    let scheme = FullDPScheme::new(DPConfig::default());
    let lhs = atomic(0, 100, false, 1);
    let rhs = atomic(10, 110, true, 1);
    let (penalty, edge_type) = scheme.edge_penalty(&lhs, &rhs, k).unwrap();
    assert_eq!(edge_type, EdgeType::SvBreak);
    assert!(penalty >= 200.0);
}

#[test]
fn edge_penalty_sv_break_large_diagonal_shift() {
    let k = 5;
    let scheme = FullDPScheme::new(DPConfig::default());
    let lhs = atomic(0, 100, false, 1);
    // read_gap = 5, ref_gap = 5005 → deviation = 5000 >> max_gap_deviation=1000
    let rhs = atomic(10, 5110, false, 1);
    let (penalty, edge_type) = scheme.edge_penalty(&lhs, &rhs, k).unwrap();
    assert_eq!(edge_type, EdgeType::SvBreak);
    assert!(penalty >= 200.0);
}

// ── chain_seeds ───────────────────────────────────────────────────────────────

#[test]
fn chain_empty_returns_none() {
    let seeds: Vec<AtomicSeed> = vec![];
    let scheme = FullDPScheme::new(DPConfig::default());
    assert!(chain_seeds(&seeds, 5, &scheme).is_none());
}

#[test]
fn chain_single_seed() {
    let k = 5;
    let scheme = FullDPScheme::new(DPConfig::default());
    let seeds = vec![atomic(0, 100, false, 1)];
    let ChainResult { score, chain, edge_types } = chain_seeds(&seeds, k, &scheme).unwrap();
    assert_eq!(chain, vec![0]);
    assert!(edge_types.is_empty());
    assert!((score - 1.0).abs() < 1e-9);
}

#[test]
fn chain_collinear_beats_isolated() {
    let k = 5;
    let scheme = FullDPScheme::new(DPConfig::default());
    // Three collinear seeds that abut (read_gap = 0, zero penalty) so chaining
    // definitely pays.  Plus one off-diagonal seed interleaved in read space.
    // read_pos order after sorting: 0, 5, 10, 15
    let seeds = vec![
        atomic(0,  100,  false, 1),   // collinear, diagonal 100
        atomic(5,  5000, false, 1),   // off-diagonal (isolated)
        atomic(5,  105,  false, 1),   // collinear, diagonal 100, abutting prev
        atomic(10, 110,  false, 1),   // collinear, diagonal 100, abutting prev
    ];
    let mut seeds = seeds;
    seeds.sort_by_key(|s| (s.read_pos(), s.ref_pos()));

    let ChainResult { chain, edge_types, .. } = chain_seeds(&seeds, k, &scheme).unwrap();

    // The winning chain should be the three collinear seeds with no SV breaks.
    assert_eq!(chain.len(), 3, "chain: {chain:?}");
    assert!(edge_types.iter().all(|&e| e == EdgeType::Continuation));
}

#[test]
fn chain_sv_break_recorded() {
    let k = 5;
    // Seeds abut (read_gap = 0) so there is no read-gap cost, only the SV
    // penalty from the large diagonal shift.  Use a tiny sv_penalty so the
    // chain is worth forming and the test is about edge classification.
    let cfg = DPConfig { sv_penalty: 0.1, ..DPConfig::default() };
    let scheme = FullDPScheme::new(cfg);
    let seeds = vec![
        atomic(0, 100,  false, 1),
        atomic(5, 5110, false, 1),  // abutting, large diagonal shift → SV break
    ];
    let ChainResult { chain, edge_types, .. } = chain_seeds(&seeds, k, &scheme).unwrap();
    assert_eq!(chain.len(), 2);
    assert_eq!(edge_types[0], EdgeType::SvBreak);
}

#[test]
fn chain_fixture_chr13_rev_is_single_segment() {
    let path = fixture_path();
    if !std::path::Path::new(path).exists() {
        eprintln!("fixture not found at {path}, skipping");
        return;
    }

    let k = 20;
    let rows = load_seeds(path);

    // Take only chr13 reverse-strand seeds — the primary alignment diagonal.
    let chr13_rev: Vec<RawSeedRow> = rows
        .into_iter()
        .filter(|r| r.chrom == "chr13" && r.strand)
        .collect();
    let (atoms, _) = rows_to_atomic(&chr13_rev, k);

    let scheme = FullDPScheme::new(DPConfig::default());

    // Prune isolated seeds, then build compound seeds, then chain.
    let collection = SeedCollection::new(k, atoms);
    let collection = collection.prune_isolated_seeds(&scheme);
    let mut compounds = collection.compound_seeds();

    // chain_seeds requires read_pos order.
    compounds.sort_by_key(|s| s.read_pos());

    let n_compounds = compounds.len();
    let read_span = compounds.last().map(|s| s.read_end(k)).unwrap_or(0);

    let ChainResult { chain, edge_types, score } =
        chain_seeds(&compounds, k, &scheme).unwrap();

    let first = &compounds[chain[0]];
    let last  = &compounds[*chain.last().unwrap()];
    eprintln!(
        "compounds={n_compounds}  chain_len={}  score={score:.1}  \
         read_span=0..{read_span}  chain_read={}..{}",
        chain.len(), first.read_pos(), last.read_end(k),
    );

    // The primary chr13 reverse alignment should chain as a single continuous
    // segment — no SV breaks.
    let sv_breaks = edge_types.iter().filter(|&&e| e == EdgeType::SvBreak).count();
    assert_eq!(
        sv_breaks, 0,
        "expected no SV breaks in chr13 rev chain, got {sv_breaks} \
         (chain len={}, score={score:.1})",
        chain.len()
    );
}
