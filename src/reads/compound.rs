use std::iter::repeat;

use noodles::sam::alignment::{
    RecordBuf,
    record::{
        Flags,
        cigar::{Op, op::Kind},
        data::field::Tag,
    },
    record_buf::{Cigar, Data, data::field::Value},
};
use parallax::{
    reference::InMemoryReference,
    utils::{piecewise::Piecewise, sequence::complement},
};

use crate::reads::builder::build_record;

pub trait Seed {
    fn length(&self, k: usize) -> usize;

    fn diagonal(&self) -> i64 {
        self.ref_start() as i64 - self.read_start() as i64
    }

    fn read_start(&self) -> u32;

    fn read_end(&self, k: usize) -> u32;

    fn chrom_id(&self) -> u32;

    fn is_reverse(&self) -> bool;

    fn ref_start(&self) -> u32;

    fn ref_end(&self, k: usize) -> u32;

    fn multiplicity(&self) -> u32;

    /// Verified exact-match extension beyond the seed's k-mer boundaries, in
    /// read coordinates.  Returns `(left_ext, right_ext)` where `left_ext` is
    /// the number of matching bases to the left of `read_start()` and
    /// `right_ext` is the number to the right of `read_end(k)`.
    fn extension(&self) -> (u32, u32) {
        (0, 0)
    }

    fn to_string(&self, k: usize) -> String;
}

pub trait Weighted: Seed {
    fn weight(&self) -> f64;
}

/// The result of resolving the gap between two seeds, including the optimal
/// split of any overlapping atoms.
pub struct GapResult {
    /// Reference-space gap between the effective end of lhs and effective start
    /// of rhs after trimming.  Positive = gap, negative = overlap within
    /// tolerance, semantics match `edge_penalty_v2`.
    pub ref_gap: i64,
    /// Read-space gap (0 when the seeds overlap in read coordinates).
    pub read_gap: i64,
}

/// Implemented by seed types that can compute the optimal gap to a following
/// seed, resolving any read-coordinate overlap by finding the split point that
/// minimises total weight trimmed (accounting for per-atom frequencies).
pub trait GapComputable: Weighted {
    /// Storage for the trimmed seeds produced by `gap_to` when two seeds
    /// Compute the gap from `self` (lhs) to `other` (rhs).
    ///
    /// Returns `None` if `other` is fully consumed by `self` and should be
    /// skipped by the DP.
    ///
    /// Otherwise returns `Some((gap, lhs_retained, rhs_retained))` where:
    ///
    /// - `gap` describes the reference- and read-space gap after resolving any
    ///   overlap.
    /// - `lhs_retained` is `Some(s)` if lhs was trimmed, holding only the
    ///   atoms that were kept.  `None` means lhs is used in full.
    /// - `rhs_retained` is the symmetric value for rhs.
    fn gap_to(&self, other: &Self, k: usize) -> Option<(GapResult, Option<Self>, Option<Self>)>
    where
        Self: Sized;
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EdgeType {
    Continuation,
    SvBreak,
    Repeat,
}

/// Configuration shared by all `ChainingDPScheme` implementations.
pub struct DPConfig {
    /// SV break penalty used by `chain_seeds` (flat DP).
    pub sv_penalty: f64,
    pub repeat_expansion_penalty: f64,
    pub repeat_expansion_max_ref_window: u32,
    /// Cost of read-space gap (bases of uncovered read).  Three slopes: low for
    /// short gaps, steeper in the middle, steepest above `threshold_hi` to
    /// strongly penalise long uncovered stretches.
    pub read_gap_cost: Piecewise<2>,
    /// Cost of ref-space deviation `|ref_gap - read_gap|` along a Continuation
    /// edge.  Two slopes: high for small deviations, lower above `threshold` to
    /// flatten the cost for large indel-like deviations.
    pub ref_dev_cost: Piecewise<1>,
    pub max_gap_deviation: f64,
    pub ref_overlap_tolerance: i64,
}

impl Default for DPConfig {
    fn default() -> Self {
        Self {
            sv_penalty: 30.0,
            repeat_expansion_penalty: 20.0,
            repeat_expansion_max_ref_window: 400,
            read_gap_cost: Piecewise {
                base: 0.02,
                breakpoints: [15.0, 60.0],
                slope_increments: [0.03, 0.15],
            },
            ref_dev_cost: Piecewise {
                base: 0.01,
                breakpoints: [50.0],
                slope_increments: [-0.009],
            },
            max_gap_deviation: 1000.0,
            ref_overlap_tolerance: 10,
        }
    }
}

impl DPConfig {
    /// Config tuned for `chain_seeds_multi`.  Uses a lower `sv_penalty` so that
    /// a well-supported seed on a different strand/chrom can justify one or two
    /// breaks; the escalating cost structure in `chain_seeds_multi` then
    /// discourages long chains of spurious breaks.
    pub fn default_multi() -> Self {
        Self {
            sv_penalty: 6.0,
            repeat_expansion_penalty: 4.0,
            ..Self::default()
        }
    }
}

/// A chaining DP scheme.  `edge_penalty` is generic over any `S: GapComputable`
/// so a single `FullDPScheme` instance handles both `AtomicSeed`s and
/// `CompoundSeed`s without dynamic dispatch or duplication.
///
/// The DP recurrence using this trait is:
///   `dp[i] = max over j: dp[j] + rhs.weight() - edge_penalty(lhs, rhs)`
///
/// where `edge_penalty` already incorporates any weight lost to trimming at
/// the lhs/rhs boundary (see `GapComputable::gap_to`).
pub trait ChainingDPScheme {
    /// Returns `None` if `rhs` is fully consumed by `lhs`.  Otherwise returns
    /// `Some((total_penalty, edge_type))` where `total_penalty` already
    /// incorporates the trim cost; the caller scores the edge as
    /// `dp[lhs] + rhs.weight() - total_penalty`.
    fn edge_penalty<S: GapComputable>(&self, lhs: &S, rhs: &S, k: usize)
    -> Option<(f64, EdgeType)>;

    /// Maximum read-space gap (in bases) within which a non-SV neighbour can
    /// exist.  Used by `prune_isolated_seeds` to bound the forward scan.
    fn max_neighbour_gap(&self) -> i64;
}

// ── Seed types ────────────────────────────────────────────────────────────────

pub struct GenericSeed {
    read_pos: u32,
    chrom_id: u32,
    ref_pos: u32,
    length: u32,
    is_reverse: bool,
    multiplicity: u32,
    weight: f32,
}

impl GenericSeed {
    pub fn from_seed<S: Seed + Weighted>(seed: &S, k: usize) -> Self {
        GenericSeed {
            read_pos: seed.read_start(),
            chrom_id: seed.chrom_id(),
            ref_pos: seed.ref_start(),
            length: seed.length(k) as u32,
            is_reverse: seed.is_reverse(),
            multiplicity: seed.multiplicity(),
            weight: seed.weight() as f32,
        }
    }

    /// Trim `bases` off the left (read-space) side of the seed.
    /// For forward strand, advances `ref_pos`; for reverse, shrinks `length`
    /// (the high ref end is the read-left side on reverse strand).
    /// Returns any leftover that couldn't be absorbed because the seed reached zero length.
    pub fn trim_left(&mut self, bases: u32) -> u32 {
        let trimmed = bases.min(self.length);
        self.read_pos += trimmed;
        if !self.is_reverse {
            self.ref_pos += trimmed;
        }
        self.length -= trimmed;
        bases - trimmed
    }

    /// Trim `bases` off the right (read-space) side of the seed.
    /// For reverse strand, advances `ref_pos` (the low ref end is the read-right side);
    /// for forward strand, just shrinks `length`.
    /// Returns any leftover that couldn't be absorbed because the seed reached zero length.
    pub fn trim_right(&mut self, bases: u32) -> u32 {
        let trimmed = bases.min(self.length);
        if self.is_reverse {
            self.ref_pos += trimmed;
        }
        self.length -= trimmed;
        bases - trimmed
    }
}

impl Seed for GenericSeed {
    fn is_reverse(&self) -> bool {
        self.is_reverse
    }

    fn length(&self, _k: usize) -> usize {
        self.length as usize
    }

    fn multiplicity(&self) -> u32 {
        self.multiplicity
    }

    fn to_string(&self, _k: usize) -> String {
        repeat('N').take(self.length as usize).collect::<String>()
    }

    fn read_start(&self) -> u32 {
        self.read_pos
    }

    fn read_end(&self, _k: usize) -> u32 {
        self.read_pos + self.length
    }

    fn chrom_id(&self) -> u32 {
        self.chrom_id
    }

    fn ref_start(&self) -> u32 {
        self.ref_pos
    }

    fn ref_end(&self, _k: usize) -> u32 {
        self.ref_pos + self.length
    }
}

impl Weighted for GenericSeed {
    fn weight(&self) -> f64 {
        self.weight as f64
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomicSeed {
    /// Position in the original (forward-sense) read, regardless of strand.
    /// For reverse-strand seeds this is `read_len - strand_local_pos - k`.
    read_pos: u32,
    chrom_id: u32,
    /// Positive = forward (1-based ref pos), negative = reverse (1-based ref pos).
    ref_pos: i32,
    kmer: u64,
    kmer_multiplicity: u32,
}

impl AtomicSeed {
    /// Construct an `AtomicSeed`.
    ///
    /// `strand_local_pos` is the k-mer start position in the strand-local sequence
    /// (position in the RC sequence for reverse-strand seeds).  `read_len` and `k`
    /// are used to convert reverse-strand positions to forward-read coordinates:
    /// `read_pos = read_len - strand_local_pos - k`.
    pub fn new(
        strand_local_pos: u32,
        read_len: u32,
        k: usize,
        chrom_id: u32,
        ref_pos: u32,
        is_reverse: bool,
        kmer: u64,
        kmer_multiplicity: u32,
    ) -> Self {
        let read_pos = if is_reverse {
            read_len - strand_local_pos - k as u32
        } else {
            strand_local_pos
        };
        let ref_pos_stored = (ref_pos + 1) as i32;
        let ref_pos_stored = if is_reverse {
            -ref_pos_stored
        } else {
            ref_pos_stored
        };
        Self {
            read_pos,
            chrom_id,
            ref_pos: ref_pos_stored,
            kmer,
            kmer_multiplicity,
        }
    }

    pub fn ref_pos_and_strand(&self) -> (u32, bool) {
        if self.ref_pos < 0 {
            (-self.ref_pos as u32 - 1, true)
        } else {
            (self.ref_pos as u32 - 1, false)
        }
    }

    pub fn kmer(&self) -> u64 {
        self.kmer
    }
    pub fn kmer_multiplicity(&self) -> u32 {
        self.kmer_multiplicity
    }

    pub fn atom_weight(&self) -> f64 {
        1.0 / (self.kmer_multiplicity as f64).sqrt()
    }

    /// Signed reference gap from the end of `self` to the start of `other`,
    /// in the direction of read traversal.  Positive = gap; negative = overlap.
    /// Both atoms must be on the same strand.
    pub fn ref_gap_to(&self, other: &Self, k: usize) -> i64 {
        if self.is_reverse() {
            // ref decreases as read_pos increases; lhs has larger ref than rhs.
            self.ref_start() as i64 - (other.ref_end(k) as i64)
        } else {
            other.ref_start() as i64 - (self.ref_end(k) as i64)
        }
    }
}

impl Seed for AtomicSeed {
    fn read_start(&self) -> u32 {
        self.read_pos
    }

    fn read_end(&self, k: usize) -> u32 {
        self.read_pos + k as u32
    }

    fn chrom_id(&self) -> u32 {
        self.chrom_id
    }

    fn is_reverse(&self) -> bool {
        self.ref_pos_and_strand().1
    }

    fn ref_start(&self) -> u32 {
        self.ref_pos_and_strand().0
    }

    fn ref_end(&self, k: usize) -> u32 {
        self.ref_start() + k as u32
    }

    fn length(&self, k: usize) -> usize {
        k
    }

    fn multiplicity(&self) -> u32 {
        self.kmer_multiplicity
    }

    /// Diagonal index — constant along a gapless match.
    /// Forward: `ref_pos - read_pos` (ref increases with read_pos).
    /// Reverse: `ref_pos + read_pos` (ref decreases as read_pos increases).
    fn diagonal(&self) -> i64 {
        if self.is_reverse() {
            self.ref_start() as i64 + self.read_start() as i64
        } else {
            self.ref_start() as i64 - self.read_start() as i64
        }
    }

    fn to_string(&self, k: usize) -> String {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        let mut bases: [char; 32] = ['N'; 32];
        let mut x = self.kmer;
        for i in 0..k {
            bases[i] = BASES[(x & 3) as usize];
            x >>= 2;
        }
        bases[0..k].into_iter().rev().collect::<String>()
    }
}

impl Weighted for AtomicSeed {
    fn weight(&self) -> f64 {
        self.atom_weight()
    }
}

impl GapComputable for AtomicSeed {
    fn gap_to(&self, other: &Self, k: usize) -> Option<(GapResult, Option<Self>, Option<Self>)> {
        // In forward-read coordinates, lhs.read_pos < rhs.read_pos always.
        // read_gap > 0: gap between seeds; < 0: overlap.
        let read_gap = other.read_pos as i64 - (self.read_pos as i64 + k as i64);

        // Cross-chrom or cross-strand: SV break sentinel.
        if self.chrom_id != other.chrom_id || self.is_reverse() != other.is_reverse() {
            return Some((
                GapResult {
                    ref_gap: i64::MIN,
                    read_gap,
                },
                None,
                None,
            ));
        }

        let ref_gap = self.ref_gap_to(other, k);

        if read_gap >= 0 {
            Some((GapResult { ref_gap, read_gap }, None, None))
        } else if (-read_gap) as usize >= k {
            None
        } else {
            // Partial overlap: only valid on the same diagonal.
            if self.diagonal() != other.diagonal() {
                return None;
            }
            Some((GapResult { ref_gap, read_gap }, None, None))
        }
    }
}

/// A collection of `AtomicSeed`s sorted by (chrom, diagonal, read_pos).
pub struct SeedCollection {
    pub k: usize,
    pub hits: Vec<AtomicSeed>,
}

impl SeedCollection {
    pub fn new(k: usize, hits: Vec<AtomicSeed>) -> Self {
        let mut hits = hits;
        hits.sort_by_key(|hit| (hit.chrom_id, hit.diagonal(), hit.read_pos));
        Self { k, hits }
    }

    /// Remove seeds that have no colinear neighbour — i.e. every reachable
    /// neighbour within `max_neighbour_gap` read bases would require an SV
    /// break edge.  Such seeds can only ever appear as isolated single-seed
    /// segments and contribute nothing to a multi-seed chain.
    ///
    /// Seeds are sorted by read_pos for the scan.  For each seed we look
    /// forward until the read gap exceeds `max_neighbour_gap`; if we find at
    /// least one neighbour whose edge is not an SvBreak, both seeds are marked
    /// as having a neighbour.  Seeds with no neighbour are dropped.
    /// Merge adjacent atomic seeds on the same diagonal into `CompoundSeed`s.
    ///
    /// Iterates the hits in `(chrom, diagonal, read_pos)` order and groups
    /// consecutive atoms that overlap or abut on the same diagonal
    /// (`-k < read_gap <= 0` and `ref_gap == read_gap`).  Any positive read
    /// gap, a cross-chrom/strand result, or a diagonal mismatch starts a new
    /// compound seed.
    ///
    /// This is functionally equivalent to the seed-extension merge pass in the
    /// old pipeline but preserves per-atom frequencies for weight computation.
    pub fn compound_seeds(&self) -> Vec<CompoundSeed<'_>> {
        let k = self.k;
        if self.hits.is_empty() {
            return Vec::new();
        }

        let mut result = Vec::new();
        let mut group_start = 0usize;

        for i in 1..self.hits.len() {
            let prev = &self.hits[i - 1];
            let curr = &self.hits[i];

            let merge = matches!(
                prev.gap_to(curr, k),
                Some((GapResult { read_gap, ref_gap }, _, _))
                    if ref_gap != i64::MIN && read_gap <= 0 && ref_gap == read_gap
            );

            if !merge {
                result.push(CompoundSeed::new(&self.hits[group_start..i]));
                group_start = i;
            }
        }
        result.push(CompoundSeed::new(&self.hits[group_start..]));
        result
    }

    pub fn prune_isolated_seeds(mut self, scheme: &impl ChainingDPScheme) -> Self {
        let k = self.k;
        let n = self.hits.len();

        // Work in read-position order for the forward scan.
        let mut by_read: Vec<usize> = (0..n).collect();
        by_read.sort_unstable_by_key(|&i| self.hits[i].read_pos);

        let mut has_neighbour = vec![false; n];

        for rank in 0..n {
            let i = by_read[rank];
            if has_neighbour[i] {
                continue;
            }
            let si = &self.hits[i];
            let si_read_end = si.read_end(k);

            for rank2 in rank + 1..n {
                let j = by_read[rank2];
                let sj = &self.hits[j];

                // Once the read gap exceeds max_neighbour_gap no further
                // neighbour can produce a non-SV edge.
                let read_gap = sj.read_start() as i64 - si_read_end as i64;
                if read_gap > scheme.max_neighbour_gap() {
                    break;
                }

                if let Some((_, edge_type)) = scheme.edge_penalty(si, sj, k) {
                    if edge_type != EdgeType::SvBreak {
                        has_neighbour[i] = true;
                        has_neighbour[j] = true;
                        break;
                    }
                }
            }
        }

        // Retain only seeds that have at least one non-SV neighbour, preserving
        // the original (chrom, diagonal, read_pos) sort order.
        let hits = std::mem::take(&mut self.hits);
        self.hits = hits
            .into_iter()
            .enumerate()
            .filter(|(i, _)| has_neighbour[*i])
            .map(|(_, s)| s)
            .collect();

        self
    }
}

/// A merged seed: a contiguous run of `AtomicSeed`s on the same diagonal,
/// borrowed from a `SeedCollection`.  No heap allocation.
pub struct CompoundSeed<'a> {
    atoms: &'a [AtomicSeed],
    left_ext: u32,
    right_ext: u32,
    min_multiplicity: u32,
    weight: f64,
}

impl<'a> CompoundSeed<'a> {
    pub fn new(atoms: &'a [AtomicSeed]) -> Self {
        debug_assert!(!atoms.is_empty());
        let min_multiplicity = atoms.iter().map(|a| a.multiplicity()).min().unwrap();
        let weight = atoms.iter().map(|a| a.atom_weight()).sum();
        Self {
            atoms,
            left_ext: 0,
            right_ext: 0,
            min_multiplicity,
            weight,
        }
    }

    pub fn with_extensions(atoms: &'a [AtomicSeed], left_ext: u32, right_ext: u32) -> Self {
        let mut s = Self::new(atoms);
        s = s.set_extensions(left_ext, right_ext);
        s
    }

    pub fn set_extensions(self, left_ext: u32, right_ext: u32) -> Self {
        let atom_weight: f64 = self.atoms.iter().map(|a| a.atom_weight()).sum();
        let ext_weight =
            ((left_ext + right_ext) as f64).sqrt() / (self.min_multiplicity as f64).sqrt();
        Self {
            left_ext,
            right_ext,
            weight: atom_weight + ext_weight,
            ..self
        }
    }

    pub fn atoms(&self) -> &[AtomicSeed] {
        self.atoms
    }
}

impl<'a> Seed for CompoundSeed<'a> {
    fn read_start(&self) -> u32 {
        self.atoms[0].read_start()
    }

    fn read_end(&self, k: usize) -> u32 {
        let n = self.atoms.len();
        self.atoms[n - 1].read_end(k)
    }

    fn chrom_id(&self) -> u32 {
        self.atoms[0].chrom_id()
    }

    fn is_reverse(&self) -> bool {
        self.atoms[0].is_reverse()
    }

    fn diagonal(&self) -> i64 {
        self.atoms[0].diagonal()
    }

    fn length(&self, k: usize) -> usize {
        let n = self.atoms.len();
        self.atoms[n - 1].read_end(k) as usize - self.atoms[0].read_start() as usize
    }

    fn ref_start(&self) -> u32 {
        if self.is_reverse() {
            self.atoms.last().unwrap().ref_start()
        } else {
            self.atoms[0].ref_start()
        }
    }

    fn ref_end(&self, k: usize) -> u32 {
        if self.is_reverse() {
            self.atoms[0].ref_end(k)
        } else {
            self.atoms.last().unwrap().ref_end(k)
        }
    }

    fn multiplicity(&self) -> u32 {
        self.min_multiplicity
    }

    fn extension(&self) -> (u32, u32) {
        (self.left_ext, self.right_ext)
    }

    fn to_string(&self, k: usize) -> String {
        let mut result = String::new();
        let mut covered_up_to: u32 = 0;
        for atom in self.atoms {
            let atom_str = atom.to_string(k);
            let atom_start = atom.read_start();
            let skip = if atom_start < covered_up_to {
                (covered_up_to - atom_start) as usize
            } else {
                0
            };
            result.push_str(&atom_str[skip..]);
            covered_up_to = atom_start + k as u32;
        }
        result
    }
}

impl<'a> Weighted for CompoundSeed<'a> {
    fn weight(&self) -> f64 {
        self.weight
    }
}

impl<'a> GapComputable for CompoundSeed<'a> {
    fn gap_to(&self, other: &Self, k: usize) -> Option<(GapResult, Option<Self>, Option<Self>)> {
        let self_read_end = self.read_end(k);
        let read_gap_signed = other.read_start() as i64 - self_read_end as i64;

        let sv_break =
            self.chrom_id() != other.chrom_id() || self.is_reverse() != other.is_reverse();

        if sv_break && read_gap_signed >= 0 {
            return Some((
                GapResult {
                    ref_gap: i64::MIN,
                    read_gap: read_gap_signed,
                },
                Some(CompoundSeed::new(self.atoms).set_extensions(self.left_ext, self.right_ext)),
                Some(
                    CompoundSeed::new(other.atoms).set_extensions(other.left_ext, other.right_ext),
                ),
            ));
        }

        // Fully consumed (same-strand only).
        if !sv_break && -read_gap_signed >= other.length(k) as i64 {
            return None;
        }

        if read_gap_signed >= 0 {
            // No atom-level read overlap.  Compute the raw gap from atom boundaries.
            let lhs_last = self.atoms.last().unwrap();
            let rhs_first = &other.atoms[0];
            let read_gap = rhs_first.read_start() as i64 - lhs_last.read_end(k) as i64;
            let ref_gap = lhs_last.ref_gap_to(rhs_first, k);

            let lhs_ext_right = self_read_end as i64 + self.right_ext as i64;
            let rhs_ext_left = other.read_start() as i64 - other.left_ext as i64;
            let ext_read_overlap = lhs_ext_right - rhs_ext_left;

            let ext_ref_overlap = if !sv_break {
                -ref_gap + self.right_ext as i64 + other.left_ext as i64
            } else {
                0
            };

            let overlap = ext_read_overlap.max(ext_ref_overlap).max(0);

            if overlap == 0 {
                return Some((
                    GapResult { ref_gap, read_gap },
                    Some(
                        CompoundSeed::new(self.atoms).set_extensions(self.left_ext, self.right_ext),
                    ),
                    Some(
                        CompoundSeed::new(other.atoms)
                            .set_extensions(other.left_ext, other.right_ext),
                    ),
                ));
            }

            let lhs_avail = self.right_ext as i64;
            let rhs_avail = other.left_ext as i64;

            if overlap <= lhs_avail + rhs_avail {
                let (lhs_trim, rhs_trim) = if self.weight <= other.weight {
                    let lt = overlap.min(lhs_avail);
                    (lt, overlap - lt)
                } else {
                    let rt = overlap.min(rhs_avail);
                    (overlap - rt, rt)
                };
                let new_lhs_right = self.right_ext - lhs_trim as u32;
                let new_rhs_left = other.left_ext - rhs_trim as u32;
                return Some((
                    GapResult {
                        ref_gap: 0,
                        read_gap: 0,
                    },
                    Some(
                        CompoundSeed::new(self.atoms).set_extensions(self.left_ext, new_lhs_right),
                    ),
                    Some(
                        CompoundSeed::new(other.atoms)
                            .set_extensions(new_rhs_left, other.right_ext),
                    ),
                ));
            }

            // Overlap exceeds combined extension margin — return raw gap.
            return Some((
                GapResult { ref_gap, read_gap },
                Some(CompoundSeed::new(self.atoms).set_extensions(self.left_ext, self.right_ext)),
                Some(
                    CompoundSeed::new(other.atoms).set_extensions(other.left_ext, other.right_ext),
                ),
            ));
        }

        // ── Atom-level read overlap ───────────────────────────────────────────────
        let lhs_ext_end = self_read_end as i64 + self.right_ext as i64;
        let rhs_ext_start = other.read_start() as i64 - other.left_ext as i64;
        let total_overlap = lhs_ext_end - rhs_ext_start;

        if total_overlap <= 0 {
            let lhs_last = self.atoms.last().unwrap();
            let rhs_first = &other.atoms[0];
            let read_gap = rhs_first.read_start() as i64 - lhs_last.read_end(k) as i64;
            let ref_gap = lhs_last.ref_gap_to(rhs_first, k);
            return Some((
                GapResult { ref_gap, read_gap },
                Some(CompoundSeed::new(self.atoms).set_extensions(self.left_ext, self.right_ext)),
                Some(
                    CompoundSeed::new(other.atoms).set_extensions(other.left_ext, other.right_ext),
                ),
            ));
        }

        let lhs_ext_avail = self.right_ext as i64;
        let rhs_ext_avail = other.left_ext as i64;
        let ext_capacity = lhs_ext_avail + rhs_ext_avail;

        let (ext_lhs_trim, ext_rhs_trim) = if ext_capacity >= total_overlap {
            if self.weight <= other.weight {
                let lt = total_overlap.min(lhs_ext_avail);
                (lt, total_overlap - lt)
            } else {
                let rt = total_overlap.min(rhs_ext_avail);
                (total_overlap - rt, rt)
            }
        } else {
            (lhs_ext_avail, rhs_ext_avail)
        };

        let atom_overlap = (total_overlap - ext_lhs_trim - ext_rhs_trim).max(0);
        let new_lhs_right = self.right_ext - ext_lhs_trim as u32;
        let new_rhs_left = other.left_ext - ext_rhs_trim as u32;

        if atom_overlap == 0 {
            return Some((
                GapResult {
                    ref_gap: 0,
                    read_gap: 0,
                },
                Some(CompoundSeed::new(self.atoms).set_extensions(self.left_ext, new_lhs_right)),
                Some(CompoundSeed::new(other.atoms).set_extensions(new_rhs_left, other.right_ext)),
            ));
        }

        // Atom-level overlap remains — run the split-point optimisation.
        let lhs_overlap_start = self
            .atoms
            .partition_point(|a| a.read_end(k) <= other.read_start());
        let lhs_overlap = &self.atoms[lhs_overlap_start..];
        let rhs_overlap_end = other
            .atoms
            .partition_point(|a| a.read_start() < self_read_end);
        let rhs_overlap = &other.atoms[..rhs_overlap_end];

        let mut candidates: Vec<u32> = Vec::with_capacity(lhs_overlap.len() + rhs_overlap.len());
        for a in lhs_overlap {
            candidates.push(a.read_end(k));
        }
        for a in rhs_overlap {
            candidates.push(a.read_start());
        }
        candidates.sort_unstable();
        candidates.dedup();

        let mut best: Option<(f64, u32, GapResult)> = None;

        for &p in &candidates {
            let lhs_atom_trimmed: f64 = lhs_overlap
                .iter()
                .filter(|a| a.read_end(k) > p)
                .map(|a| a.atom_weight())
                .sum();
            let rhs_atom_trimmed: f64 = rhs_overlap
                .iter()
                .filter(|a| a.read_start() < p)
                .map(|a| a.atom_weight())
                .sum();
            let total_trimmed = lhs_atom_trimmed + rhs_atom_trimmed;

            let lhs_effective_last = self.atoms.iter().rev().find(|a| a.read_end(k) <= p);
            let rhs_effective_first = other.atoms.iter().find(|a| a.read_end(k) >= p);

            let (ref_gap, read_gap) = if sv_break {
                let read_gap = match (lhs_effective_last, rhs_effective_first) {
                    (Some(l), Some(r)) => r.read_start() as i64 - (l.read_end(k)) as i64,
                    (None, Some(r)) => r.read_start() as i64 - self.atoms[0].read_start() as i64,
                    (Some(l), None) => {
                        other.atoms.last().unwrap().read_start() as i64 - l.read_end(k) as i64
                    }
                    (None, None) => continue,
                };
                (i64::MIN, read_gap)
            } else {
                match (lhs_effective_last, rhs_effective_first) {
                    (Some(l), Some(r)) => (
                        l.ref_gap_to(r, k),
                        r.read_start() as i64 - l.read_end(k) as i64,
                    ),
                    (None, Some(r)) => {
                        let l = &self.atoms[0];
                        (
                            l.ref_gap_to(r, k),
                            r.read_start() as i64 - l.read_end(k) as i64,
                        )
                    }
                    (Some(l), None) => {
                        let r = other.atoms.last().unwrap();
                        (
                            l.ref_gap_to(r, k),
                            r.read_start() as i64 - l.read_end(k) as i64,
                        )
                    }
                    (None, None) => continue,
                }
            };

            if best.as_ref().map_or(true, |&(bt, ..)| total_trimmed < bt) {
                best = Some((total_trimmed, p, GapResult { ref_gap, read_gap }));
            }
        }

        best.map(|(_, best_p, gap)| {
            let lhs_split = self
                .atoms
                .partition_point(|a| a.read_end(k) <= best_p);
            let rhs_split = other.atoms.partition_point(|a| a.read_start() < best_p);

            let lhs_out = if lhs_split < self.atoms.len() && lhs_split > 0 {
                CompoundSeed::new(&self.atoms[..lhs_split])
                    .set_extensions(self.left_ext, new_lhs_right)
            } else {
                CompoundSeed::new(self.atoms).set_extensions(self.left_ext, new_lhs_right)
            };
            let rhs_out = if rhs_split > 0 && rhs_split < other.atoms.len() {
                CompoundSeed::new(&other.atoms[rhs_split..])
                    .set_extensions(new_rhs_left, other.right_ext)
            } else {
                CompoundSeed::new(other.atoms).set_extensions(new_rhs_left, other.right_ext)
            };
            (gap, Some(lhs_out), Some(rhs_out))
        })
    }
}

// ── Seed extension ────────────────────────────────────────────────────────────

/// Compute the exact-match extension of `seed` in read coordinates.
///
/// Walks outward from each boundary of the seed, comparing the strand-local
/// query sequence against the reference, stopping at the first mismatch or
/// after `max_ext` bases.
///
/// Returns `(left_ext, right_ext)` where `left_ext` is the number of matching
/// bases to the left of `seed.read_start()` and `right_ext` is the number to
/// the right of `seed.read_end(k)`, both in read coordinates.
pub fn compute_extensions(
    seed: &impl Seed,
    k: usize,
    query: &[u8],
    query_rc: &[u8],
    reference: &InMemoryReference,
    max_ext: u32,
) -> (u32, u32) {
    let chrom = seed.chrom_id() as usize;
    let is_rev = seed.is_reverse();
    let ref_seq = reference.sequence(chrom);

    // For a forward seed the strand-local query is `query` and both query and
    // ref run in the same direction.  The seed occupies:
    //   query[read_start .. read_end]  ↔  ref[ref_start .. ref_end]
    // Extending forward (toward higher indices in both) gives right_ext;
    // extending backward gives left_ext.  left_ext == right_ext in read coords.
    //
    // For a reverse seed the strand-local query is `query_rc`.  In strand-local
    // coordinates the seed runs from low to high, mirroring the ref:
    //   query_rc[sl_start .. sl_end]  ↔  ref[ref_start .. ref_end]
    // where sl_start = query_len - read_end, sl_end = query_len - read_start.
    // Extending forward in strand-local coords (higher sl, lower ref) corresponds
    // to extending LEFT in read coords; extending backward (lower sl, lower ref)
    // corresponds to extending RIGHT in read coords.

    let (strand_seq, sl_fwd_start, sl_bwd_start, ref_fwd_start, ref_bwd_start) = if is_rev {
        let sl_end = query.len() - seed.read_start() as usize; // strand-local end of seed
        let sl_start = query.len() - seed.read_end(k) as usize; // strand-local start of seed
        // Forward in strand-local = left in read = toward HIGHER ref (since reverse reads
        // map high-ref to low-read-pos). ref_fwd_start is the HIGH ref coord (ref_end).
        // Backward in strand-local = right in read = toward LOWER ref. ref_bwd_start = ref_start.
        (
            query_rc,
            sl_end,
            sl_start,
            seed.ref_end(k) as usize,
            seed.ref_start() as usize,
        )
    } else {
        let read_end = seed.read_end(k) as usize;
        let read_start = seed.read_start() as usize;
        (
            query,
            read_end,
            read_start,
            seed.ref_end(k) as usize,
            seed.ref_start() as usize,
        )
    };

    // Extend forward in strand-local coords.
    let fwd_max = max_ext
        .min((strand_seq.len() - sl_fwd_start) as u32)
        .min((ref_seq.len() - ref_fwd_start) as u32) as usize;
    let fwd_ext = (0..fwd_max)
        .take_while(|&i| strand_seq[sl_fwd_start + i] == ref_seq[ref_fwd_start + i])
        .count() as u32;

    // Extend backward in strand-local coords (both indices count down from their starts).
    // Cap to avoid underflow: can't extend further back than the smaller of the two starts.
    let bwd_max = max_ext.min(sl_bwd_start as u32).min(ref_bwd_start as u32) as usize;
    let bwd_ext = (1..=bwd_max)
        .take_while(|&i| strand_seq[sl_bwd_start - i] == ref_seq[ref_bwd_start - i])
        .count() as u32;

    // Map strand-local fwd/bwd back to read-coordinate left/right.
    if is_rev {
        // fwd in strand-local = left in read; bwd = right in read
        (fwd_ext, bwd_ext)
    } else {
        // fwd in strand-local = right in read; bwd = left in read
        (bwd_ext, fwd_ext)
    }
}

// ── Chain bounds ─────────────────────────────────────────────────────────────

/// Compute the effective boundaries of each seed in a chain as `GenericSeed`s,
/// after extension trimming and atom-split trimming have been resolved.
///
/// Each output seed has `read_pos` = effective left read edge and `length` =
/// effective span, encoding both boundaries.  A single forward pass calls
/// `gap_to` on each consecutive pair of original seeds; the trimmed right edge
/// of lhs and trimmed left edge of rhs are written back into the copies
/// immediately, so each seed's right boundary is settled before it is used as
/// lhs in the next iteration.
///
/// Seeds trimmed to zero length (squeezed from both sides) are left at zero —
/// they contribute a zero-length match and keep the gap alignments contiguous.
pub fn compute_chain_bounds<'a>(
    seeds: &'a [CompoundSeed<'a>],
    chain: &[usize],
    edge_types: &[EdgeType],
    k: usize,
) -> Vec<GenericSeed> {
    let n = chain.len();
    if n == 0 {
        return Vec::new();
    }

    // Initialise each copy with the full extended boundaries:
    //   read_pos  = read_start - left_ext  (effective left read edge)
    //   length    = seed.length(k) + left_ext + right_ext
    //   ref_pos   = genomic low coord of extended seed (ref_start() - left_ext for fwd,
    //               ref_start() - right_ext for rev since right_ext extends toward lower ref)
    let mut copies: Vec<GenericSeed> = chain
        .iter()
        .map(|&i| {
            let s = &seeds[i];
            let (left_ext, right_ext) = s.extension();
            let read_pos = s.read_start().saturating_sub(left_ext);
            let length = s.length(k) as u32 + left_ext + right_ext;
            // ref_pos is always the low ref coord of the extended seed.
            // fwd: ref_start (= low coord) minus left_ext extends left.
            // rev: ref_start (= low coord of last atom) minus right_ext extends right-in-read
            //      which is the low-ref direction.
            let ref_pos = if s.is_reverse() {
                s.ref_start().saturating_sub(right_ext)
            } else {
                s.ref_start().saturating_sub(left_ext)
            };
            GenericSeed {
                read_pos,
                chrom_id: s.chrom_id(),
                ref_pos,
                length,
                is_reverse: s.is_reverse(),
                multiplicity: s.multiplicity(),
                weight: s.weight() as f32,
            }
        })
        .collect();

    for rhs_rank in 1..n {
        let lhs_rank = rhs_rank - 1;

        let is_continuation = edge_types[lhs_rank] == EdgeType::Continuation;

        // ── Read-space (all pairs): trim lhs right first, spill to rhs left ─────
        let lhs_right = copies[lhs_rank].read_pos + copies[lhs_rank].length;
        let rhs_left = copies[rhs_rank].read_pos;
        let read_overlap = lhs_right.saturating_sub(rhs_left);
        let spill = copies[lhs_rank].trim_right(read_overlap);
        let remaining = copies[rhs_rank].trim_left(spill);
        debug_assert_eq!(
            remaining,
            0,
            "read overlap exceeded combined seed lengths at lhs_rank={lhs_rank} rhs_rank={rhs_rank}\n\
             read_overlap={read_overlap} spill={spill}\n\
             lhs copy: read_pos={} length={} ref_pos={}\n\
             rhs copy: read_pos={} length={} ref_pos={}",
            copies[lhs_rank].read_pos,
            copies[lhs_rank].length,
            copies[lhs_rank].ref_pos,
            copies[rhs_rank].read_pos,
            copies[rhs_rank].length,
            copies[rhs_rank].ref_pos,
        );

        if !is_continuation {
            continue;
        }

        // ── Ref-space (continuation pairs only): same priority ─────────────────
        let ref_overlap = if !copies[lhs_rank].is_reverse {
            // Forward: lhs ref right = ref_pos + length; rhs ref left = ref_pos.
            (copies[lhs_rank].ref_pos + copies[lhs_rank].length)
                .saturating_sub(copies[rhs_rank].ref_pos)
        } else {
            // Reverse: lhs ref left = ref_pos; rhs ref right = ref_pos + length.
            (copies[rhs_rank].ref_pos + copies[rhs_rank].length)
                .saturating_sub(copies[lhs_rank].ref_pos)
        };
        let spill = copies[lhs_rank].trim_right(ref_overlap);
        let remaining = copies[rhs_rank].trim_left(spill);
        assert_eq!(
            remaining,
            0,
            "ref overlap exceeded combined seed lengths at lhs_rank={lhs_rank} rhs_rank={rhs_rank}\n\
             ref_overlap={ref_overlap} spill={spill}\n\
             lhs copy after read trim: read_pos={} length={} ref_pos={} is_reverse={}\n\
             rhs copy after read trim: read_pos={} length={} ref_pos={} is_reverse={}",
            copies[lhs_rank].read_pos,
            copies[lhs_rank].length,
            copies[lhs_rank].ref_pos,
            copies[lhs_rank].is_reverse,
            copies[rhs_rank].read_pos,
            copies[rhs_rank].length,
            copies[rhs_rank].ref_pos,
            copies[rhs_rank].is_reverse,
        );
    }

    copies
}

// ── isolated seed pruning ─────────────────────────────────────────────────────

pub fn prune_isolated_seeds<S: Seed + GapComputable>(
    seeds: &mut Vec<S>,
    k: usize,
    max_read_gap: i64,
    max_ref_gap: i64,
    min_isolated_seed_len: usize,
) {
    let n = seeds.len();
    let mut ordering: Vec<usize> = (0..n).collect();
    ordering.sort_unstable_by_key(|i| {
        let seed = &seeds[*i];
        seed.read_start()
    });
    let mut retain: Vec<bool> = vec![false; n];

    // First mark all seeds we are certain to keep.
    for i in 0..n {
        let j = ordering[i];
        if seeds[j].length(k) > min_isolated_seed_len {
            retain[j] = true;
        }
    }

    for i in 0..n {
        if retain[ordering[i]] {
            continue;
        }
        let lhs = &seeds[ordering[i]];
        let mut neighbour_found = false;
        for j in (i + 1)..n {
            if neighbour_found && retain[ordering[j]] {
                continue;
            }
            let rhs = &seeds[ordering[j]];
            if let Some((g, _, _)) = lhs.gap_to(rhs, k) {
                if g.read_gap > max_read_gap {
                    break;
                }
                if g.ref_gap == i64::MIN {
                    continue;
                }
                if g.ref_gap.abs() <= max_ref_gap {
                    retain[ordering[j]] = true;
                    neighbour_found = true;
                }
            }
        }
        if neighbour_found {
            retain[ordering[i]] = true;
        }
    }
    let mut i = 0;
    seeds.retain(|_| {
        let keep = retain[i];
        i += 1;
        keep
    });
}

pub enum TagValue {
    Str(String),
    Int(i64),
    Flt(f64),
}

/// Build a `RecordBuf` for a single seed, suitable for writing via a `RecordWriter`.
///
/// The record uses hard clips for the non-seed read portions and a single `=<len>`
/// CIGAR op for the seed itself.  For reverse-strand seeds the sequence is
/// reverse-complemented (matching what IGV expects for FLAG=0x10 records).
pub fn seed_to_record<S: Seed + Weighted>(
    name: &str,
    k: usize,
    read_len: usize,
    seed: &S,
    seq: &[u8],  // forward-strand query bases for seed's [read_start, read_end)
    qual: &[u8], // Phred+33 quality for the same range
    tags: Vec<(String, TagValue)>,
) -> RecordBuf {
    let len = seed.length(k);
    let read_left = seed.read_start() as usize;
    let read_right = read_len - seed.read_end(k) as usize;
    let (left_clip, right_clip) = if seed.is_reverse() {
        (read_right, read_left)
    } else {
        (read_left, read_right)
    };

    let mut cigar_ops: Vec<Op> = Vec::with_capacity(3);
    if left_clip > 0 {
        cigar_ops.push(Op::new(Kind::HardClip, left_clip));
    }
    cigar_ops.push(Op::new(Kind::SequenceMatch, len));
    if right_clip > 0 {
        cigar_ops.push(Op::new(Kind::HardClip, right_clip));
    }
    let cigar: Cigar = cigar_ops.iter().copied().collect();

    let out_seq: Vec<u8> = if seed.is_reverse() {
        seq.iter().rev().map(|&b| complement(b)).collect()
    } else {
        seq.to_vec()
    };
    let out_qual: Vec<u8> = if seed.is_reverse() {
        qual.iter().rev().copied().collect()
    } else {
        qual.to_vec()
    };

    let mapq = (seed.weight().floor() as u8).min(254);
    let mut flags = Flags::empty();
    if seed.is_reverse() {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }

    let mut data_tags: Vec<(Tag, Value)> = Vec::with_capacity(tags.len());
    for (key, value) in tags {
        let bytes = key.as_bytes();
        if bytes.len() == 2 {
            let tag = Tag::from([bytes[0], bytes[1]]);
            let v = match value {
                TagValue::Str(s) => Value::from(s.as_str()),
                TagValue::Int(i) => Value::from(i as i32),
                TagValue::Flt(f) => Value::from(f as f32),
            };
            data_tags.push((tag, v));
        }
    }
    let data: Data = data_tags.into_iter().collect();

    build_record(
        name,
        flags,
        seed.chrom_id() as usize,
        (seed.ref_start() + 1) as usize,
        mapq,
        cigar,
        None,
        None,
        &out_seq,
        &out_qual,
        data,
    )
}

#[cfg(test)]
mod prune_tests {
    use super::*;

    fn make_atomic(read_pos: u32, ref_pos: u32, chrom_id: u32, k: usize) -> AtomicSeed {
        AtomicSeed::new(
            read_pos,
            read_pos + k as u32 * 10,
            k,
            chrom_id,
            ref_pos,
            false,
            0,
            1,
        )
    }

    #[test]
    fn keeps_neighbouring_pair() {
        let k = 10;
        // Two seeds with read gap 0 and ref gap 0 — clear neighbours.
        let mut seeds = vec![make_atomic(0, 1000, 0, k), make_atomic(10, 1010, 0, k)];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn removes_isolated_seed() {
        let k = 10;
        // Three seeds: first and third are neighbours, second is isolated
        // (too far in ref from any neighbour within read gap).
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(10, 5000, 0, k), // ref gap of 3990 from seed 0 — isolated
            make_atomic(20, 1020, 0, k),
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 2);
        assert_eq!(seeds[0].ref_start(), 1000);
        assert_eq!(seeds[1].ref_start(), 1020);
    }

    #[test]
    fn removes_all_if_all_isolated() {
        let k = 10;
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(100, 2000, 0, k), // read gap 90 > max_read_gap 20
            make_atomic(200, 3000, 0, k),
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 0);
    }

    #[test]
    fn cross_chrom_seeds_not_neighbours() {
        let k = 10;
        // Same read positions, but different chroms — ref_gap sentinel, not neighbours.
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(10, 1010, 1, k), // chrom 1
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 0);
    }

    #[test]
    fn empty_input() {
        let k = 10;
        let mut seeds: Vec<AtomicSeed> = vec![];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 0);
    }

    // Tests for min_isolated_seed_len behaviour.
    // AtomicSeed::length(k) == k, so with k=10, length==10.

    #[test]
    fn long_seed_retained_without_neighbour() {
        // min_isolated_seed_len = 9 → length(10) > 9, so seed kept unconditionally.
        let k = 10;
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),   // isolated — far from any neighbour
            make_atomic(500, 9000, 0, k), // read gap 490 >> max_read_gap 20
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 9);
        // Both seeds qualify as "long" (length 10 > 9), retained unconditionally.
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn seed_at_boundary_not_retained() {
        // min_isolated_seed_len = 10 → condition is length > 10, which is false for
        // AtomicSeed with k=10. Isolated seeds should still be pruned.
        let k = 10;
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(500, 9000, 0, k), // read gap 490 >> max_read_gap 20
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 10);
        assert_eq!(seeds.len(), 0);
    }

    #[test]
    fn mix_long_and_short_seeds() {
        // min_isolated_seed_len = 9 → all seeds with k=10 are "long" and retained.
        // But with a stricter threshold (100), isolated short seeds are pruned.
        let k = 10;
        // Seed 0 and seed 2 are neighbours; seed 1 is far in ref.
        // With min_isolated_seed_len=100, seed 1 is short — must earn its place
        // via a neighbour, but has none within max_ref_gap=200.
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(10, 5000, 0, k), // ref gap 3990 — isolated
            make_atomic(20, 1020, 0, k),
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 2);
        assert_eq!(seeds[0].ref_start(), 1000);
        assert_eq!(seeds[1].ref_start(), 1020);

        // Now lower the threshold so the isolated seed is "long" and kept.
        let mut seeds2 = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(10, 5000, 0, k),
            make_atomic(20, 1020, 0, k),
        ];
        prune_isolated_seeds(&mut seeds2, k, 20, 200, 9);
        assert_eq!(seeds2.len(), 3);
    }

    #[test]
    fn large_threshold_forces_neighbour_check() {
        // min_isolated_seed_len = 100 >> k=10, so no seed is unconditionally long.
        // All seeds must earn their place via a neighbour.
        let k = 10;
        let mut seeds = vec![
            make_atomic(0, 1000, 0, k),
            make_atomic(10, 1010, 0, k),
            make_atomic(500, 9000, 0, k), // read gap 490 >> max_read_gap 20 — isolated
        ];
        prune_isolated_seeds(&mut seeds, k, 20, 200, 100);
        assert_eq!(seeds.len(), 2);
        assert_eq!(seeds[0].ref_start(), 1000);
        assert_eq!(seeds[1].ref_start(), 1010);
    }
}

#[cfg(test)]
mod ref_extent_tests {
    use super::*;

    // read_len large enough that forward/reverse positions don't wrap
    const READ_LEN: u32 = 10000;

    fn fwd_atom(strand_local_pos: u32, ref_pos: u32, k: usize) -> AtomicSeed {
        AtomicSeed::new(strand_local_pos, READ_LEN, k, 0, ref_pos, false, 0, 1)
    }

    fn rev_atom(strand_local_pos: u32, ref_pos: u32, k: usize) -> AtomicSeed {
        AtomicSeed::new(strand_local_pos, READ_LEN, k, 0, ref_pos, true, 0, 1)
    }

    #[test]
    fn forward_single_atom() {
        let k = 10;
        // atom at ref 1000, covers [1000, 1010)
        let atoms = [fwd_atom(50, 1000, k)];
        let seed = CompoundSeed::new(&atoms);
        assert_eq!(seed.ref_start(), 1000);
        assert_eq!(seed.ref_end(k), 1010);
    }

    #[test]
    fn forward_multi_atom() {
        let k = 10;
        // atoms at ref 1000, 1005, 1010 (overlapping by 5) — span [1000, 1020)
        let atoms = [
            fwd_atom(50, 1000, k),
            fwd_atom(55, 1005, k),
            fwd_atom(60, 1010, k),
        ];
        let seed = CompoundSeed::new(&atoms);
        assert_eq!(seed.ref_start(), 1000);
        assert_eq!(seed.ref_end(k), 1020);
    }

    #[test]
    fn reverse_single_atom() {
        let k = 10;
        // reverse atom: strand_local_pos=50, ref_pos=2000 covers [2000, 2010)
        let collection = SeedCollection::new(k, vec![rev_atom(50, 2000, k)]);
        let compounds = collection.compound_seeds();
        assert_eq!(compounds.len(), 1);
        let seed = &compounds[0];
        assert_eq!(seed.ref_start(), 2000);
        assert_eq!(seed.ref_end(k), 2010);
    }

    #[test]
    fn reverse_multi_atom() {
        let k = 10;
        // Reverse strand: diagonal = ref_pos + read_pos = const.
        // read_pos = READ_LEN - strand_local - k, so as strand_local increases,
        // read_pos decreases and ref_pos increases to keep the diagonal constant.
        //   strand_local=60 → read_pos=9930, ref=2010, diag=11940; ref span [2010,2020)
        //   strand_local=50 → read_pos=9940, ref=2000, diag=11940; ref span [2000,2010)
        // Atoms sorted by ascending read_pos: [9930, 9940]. Full ref span: [2000, 2020).
        let a0 = rev_atom(60, 2010, k); // read_pos=9930, ref span [2010,2020)
        let a1 = rev_atom(50, 2000, k); // read_pos=9940, ref span [2000,2010)
        // a0 has lower read_pos → goes first; a1 has higher read_pos → goes last
        let atoms = [a0, a1];
        let seed = CompoundSeed::new(&atoms);
        // ref_start = atoms.last().ref_start() = a1.ref_pos() = 2000 ✓
        // ref_end   = atoms[0].ref_end(k)      = a0.ref_pos()+k = 2010+10 = 2020 ✓
        assert_eq!(seed.ref_start(), 2000);
        assert_eq!(seed.ref_end(k), 2020);
    }
}

// ── DP scheme ─────────────────────────────────────────────────────────────────

pub struct FullDPScheme {
    pub cfg: DPConfig,
}

impl FullDPScheme {
    pub fn new(cfg: DPConfig) -> Self {
        Self { cfg }
    }

    fn read_gap_cost(&self, read_gap: i64) -> f64 {
        self.cfg.read_gap_cost.eval(read_gap as f64)
    }

    fn ref_dev_cost(&self, deviation: f64) -> f64 {
        self.cfg.ref_dev_cost.eval(deviation)
    }

    fn classify_gap(&self, gap: &GapResult, weight_trimmed: f64) -> (f64, EdgeType) {
        let c = &self.cfg;
        let read_gap_cost = self.read_gap_cost(gap.read_gap);

        if gap.ref_gap >= -c.ref_overlap_tolerance {
            let deviation = (gap.ref_gap - gap.read_gap).unsigned_abs() as f64;
            if deviation > c.max_gap_deviation {
                return (
                    read_gap_cost + c.sv_penalty + weight_trimmed,
                    EdgeType::SvBreak,
                );
            }
            (
                read_gap_cost + self.ref_dev_cost(deviation) + weight_trimmed,
                EdgeType::Continuation,
            )
        } else {
            let deviation = (gap.ref_gap - gap.read_gap).unsigned_abs() as f64;
            if c.repeat_expansion_max_ref_window > 0
                && deviation <= c.repeat_expansion_max_ref_window as f64
            {
                (
                    read_gap_cost
                        + c.repeat_expansion_penalty
                        + self.ref_dev_cost(deviation)
                        + weight_trimmed,
                    EdgeType::Repeat,
                )
            } else {
                (
                    read_gap_cost + c.sv_penalty + weight_trimmed,
                    EdgeType::SvBreak,
                )
            }
        }
    }
}

impl ChainingDPScheme for FullDPScheme {
    fn max_neighbour_gap(&self) -> i64 {
        self.cfg.max_gap_deviation as i64
    }

    fn edge_penalty<S: GapComputable>(
        &self,
        lhs: &S,
        rhs: &S,
        k: usize,
    ) -> Option<(f64, EdgeType)> {
        let (gap, lhs_eff, rhs_eff) = lhs.gap_to(rhs, k)?;
        if gap.ref_gap == i64::MIN {
            let read_gap_cost = self.read_gap_cost(gap.read_gap);
            return Some((read_gap_cost + self.cfg.sv_penalty, EdgeType::SvBreak));
        }
        let lhs_eff_w = lhs_eff.as_ref().unwrap_or(lhs).weight();
        let rhs_eff_w = rhs_eff.as_ref().unwrap_or(rhs).weight();
        let weight_trimmed = (lhs.weight() - lhs_eff_w) + (rhs.weight() - rhs_eff_w);
        Some(self.classify_gap(&gap, weight_trimmed))
    }
}

// ── Chaining DP ───────────────────────────────────────────────────────────────

/// The result of a chaining DP run.
pub struct ChainResult {
    /// DP score of the best chain.
    pub score: f64,
    /// Indices into the input slice for the seeds in the best chain, in order.
    pub chain: Vec<usize>,
    /// Edge type for each consecutive pair in `chain` (length = chain.len() - 1).
    pub edge_types: Vec<EdgeType>,
}

/// Run the chaining DP over `seeds` using `scheme`.
///
/// Seeds must be sorted by `read_pos` (ascending) before calling — the DP
/// scans forward in read coordinates only.
///
/// Returns `None` if `seeds` is empty.  Otherwise returns the best-scoring
/// chain with its traceback.
pub fn chain_seeds<S, Scheme>(seeds: &[S], k: usize, scheme: &Scheme) -> Option<ChainResult>
where
    S: GapComputable,
    Scheme: ChainingDPScheme,
{
    let n = seeds.len();
    if n == 0 {
        return None;
    }

    // Sort by read_pos so seeds are processed in read-coordinate order.
    // The returned chain indices map back to the original `seeds` slice.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| seeds[i].read_start());

    let mut dp = vec![0.0f64; n];
    let mut prev = vec![usize::MAX; n];
    let mut best_edge = vec![EdgeType::Continuation; n];

    for rank in 0..n {
        let i = order[rank];
        dp[i] = seeds[i].weight();

        // Consider all earlier seeds as potential predecessors.
        for r in (0..rank).rev() {
            let j = order[r];
            if let Some((penalty, edge_type)) = scheme.edge_penalty(&seeds[j], &seeds[i], k) {
                let candidate = dp[j] + seeds[i].weight() - penalty;
                if candidate > dp[i] {
                    dp[i] = candidate;
                    prev[i] = j;
                    best_edge[i] = edge_type;
                }
            }
        }
    }

    // Find the highest-scoring endpoint.
    let best_end = (0..n).max_by(|&a, &b| dp[a].partial_cmp(&dp[b]).unwrap())?;

    // Traceback using original indices.
    let mut chain = Vec::new();
    let mut cur = best_end;
    loop {
        chain.push(cur);
        let p = prev[cur];
        if p == usize::MAX {
            break;
        }
        cur = p;
    }
    chain.reverse();

    let edge_types = chain.windows(2).map(|w| best_edge[w[1]]).collect();

    Some(ChainResult {
        score: dp[best_end],
        chain,
        edge_types,
    })
}

/// Run the chaining DP with escalating SV-break costs.
///
/// State is `dp[i][b]` — best score ending at seed `i` having taken exactly
/// `b` SV breaks (capped at `max_sv_breaks`).  The k-th break costs
/// `k * scheme.sv_penalty()`, so the second break costs twice the first, etc.
/// This strongly discourages short spurious detours that require two breaks.
///
/// Pass a `scheme` built from `DPConfig::default_multi()` (low sv_penalty) so
/// that well-supported seeds can justify one or two breaks without the flat
/// chain's high single-break cost applying.
///
/// Seeds are sorted internally by `read_pos`; returned `chain` indices refer to
/// the original `seeds` slice.
///
/// Returns `None` if `seeds` is empty.
pub fn chain_seeds_multi<S, Scheme>(
    seeds: &[S],
    k: usize,
    scheme: &Scheme,
    max_sv_breaks: usize,
) -> Option<ChainResult>
where
    S: GapComputable,
    Scheme: ChainingDPScheme + SvPenalty,
{
    let n = seeds.len();
    if n == 0 {
        return None;
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&i| seeds[i].read_start());

    // dp[i][b]: best score ending at seeds[i] with b SV breaks taken (b capped at max_sv_breaks).
    let inf = f64::NEG_INFINITY;
    let mut dp = vec![vec![inf; max_sv_breaks + 1]; n];
    let mut prev = vec![vec![(usize::MAX, 0usize); max_sv_breaks + 1]; n];
    let mut prev_edge = vec![vec![EdgeType::Continuation; max_sv_breaks + 1]; n];

    for rank in 0..n {
        let i = order[rank];
        dp[i][0] = seeds[i].weight();

        for r in (0..rank).rev() {
            let j = order[r];
            let Some((base_penalty, edge_type)) = scheme.edge_penalty(&seeds[j], &seeds[i], k)
            else {
                continue;
            };

            let is_sv = edge_type == EdgeType::SvBreak;
            // For SV edges, edge_penalty() bakes in the flat sv_penalty.  We
            // replace that with the escalating cost, so strip it from the base.
            let gap_only_penalty = if is_sv {
                base_penalty - scheme.sv_penalty()
            } else {
                base_penalty
            };

            for b in 0..=max_sv_breaks {
                if dp[j][b] == inf {
                    continue;
                }
                let (b_new, sv_cost) = if is_sv {
                    let b_new = (b + 1).min(max_sv_breaks);
                    // Escalating cost: the b_new-th break costs b_new * sv_penalty.
                    (b_new, b_new as f64 * scheme.sv_penalty())
                } else {
                    (b, 0.0)
                };

                let candidate = dp[j][b] + seeds[i].weight() - gap_only_penalty - sv_cost;
                if candidate > dp[i][b_new] {
                    dp[i][b_new] = candidate;
                    prev[i][b_new] = (j, b);
                    prev_edge[i][b_new] = edge_type;
                }
            }
        }
    }

    // Best (seed, break_count) globally.
    let (best_i, best_b) = (0..n)
        .flat_map(|i| (0..=max_sv_breaks).map(move |b| (i, b)))
        .filter(|&(i, b)| dp[i][b] != inf)
        .max_by(|&(ai, ab), &(bi, bb)| {
            dp[ai][ab]
                .partial_cmp(&dp[bi][bb])
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;

    // Traceback: walk prev pointers, collecting seed indices and edges.
    let mut chain = Vec::new();
    let mut edge_types = Vec::new();
    let mut cur_i = best_i;
    let mut cur_b = best_b;
    loop {
        chain.push(cur_i);
        let (p_i, p_b) = prev[cur_i][cur_b];
        if p_i == usize::MAX {
            break;
        }
        edge_types.push(prev_edge[cur_i][cur_b]);
        cur_b = p_b;
        cur_i = p_i;
    }
    chain.reverse();
    edge_types.reverse();

    Some(ChainResult {
        score: dp[best_i][best_b],
        chain,
        edge_types,
    })
}

/// Provides the sv_penalty value used by `chain_seeds_multi` for escalating
/// break costs.  Implemented by `FullDPScheme`.
pub trait SvPenalty {
    fn sv_penalty(&self) -> f64;
}

impl SvPenalty for FullDPScheme {
    fn sv_penalty(&self) -> f64 {
        self.cfg.sv_penalty
    }
}

#[cfg(test)]
#[path = "compound_tests.rs"]
mod tests;
