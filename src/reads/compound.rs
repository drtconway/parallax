use parallax::utils::piecewise::Piecewise;

pub trait Seed {
    fn read_pos(&self) -> u32;

    fn chrom_id(&self) -> u32;

    fn ref_pos(&self) -> u32;

    fn is_reverse(&self) -> bool;

    fn length(&self, k: usize) -> usize;

    fn diagonal(&self) -> i64 {
        self.ref_pos() as i64 - self.read_pos() as i64
    }

    fn read_start(&self) -> u32 {
        self.read_pos()
    }

    fn read_end(&self, k: usize) -> u32 {
        self.read_pos() + self.length(k) as u32
    }

    fn ref_start(&self) -> u32 {
        self.ref_pos()
    }

    fn ref_end(&self, k: usize) -> u32 {
        self.ref_pos() + self.length(k) as u32
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
    /// overlap in read coordinates.  The memento is caller-owned and must
    /// outlive the references returned by `gap_to`.
    ///
    /// `AtomicSeed` uses `()` — it never produces overlap that requires
    /// trimming.  `CompoundSeed<'a>` uses
    /// `(Option<CompoundSeed<'a>>, Option<CompoundSeed<'a>>)`, storing the
    /// retained portions of lhs and rhs respectively as sub-slices of the
    /// originals (no allocation).
    type Memento: Default;

    /// Compute the gap from `self` (lhs) to `other` (rhs).
    ///
    /// Returns `None` if `other` is fully consumed by `self` and should be
    /// skipped by the DP.
    ///
    /// Otherwise returns `Some((gap, lhs_retained, rhs_retained))` where:
    ///
    /// - `gap` describes the reference- and read-space gap after resolving any
    ///   overlap.
    /// - `lhs_retained` is `Some(s)` if lhs was trimmed, where `s` borrows
    ///   from `memento` and holds only the atoms that were kept.  `None` means
    ///   lhs is used in full.
    /// - `rhs_retained` is the symmetric value for rhs.
    ///
    /// When there is no read-coordinate overlap both retained values are `None`
    /// and the memento is left at its default.  When there is overlap, the
    /// implementation finds the split point that minimises total trimmed weight
    /// and stores the retained sub-seeds in `memento` so their lifetimes are
    /// tied to the caller's stack frame.
    ///
    /// Callers compute the weight lost to trimming as
    /// `lhs.weight() - lhs_retained.map(|s| s.weight()).unwrap_or(lhs.weight())`
    /// and symmetrically for rhs.
    fn gap_to<'a>(
        &'a self,
        other: &'a Self,
        k: usize,
        memento: &'a mut Self::Memento,
    ) -> Option<(GapResult, Option<&'a Self>, Option<&'a Self>)>;
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
            self.ref_pos() as i64 - (other.ref_pos() as i64 + k as i64)
        } else {
            other.ref_pos() as i64 - (self.ref_pos() as i64 + k as i64)
        }
    }
}

impl Seed for AtomicSeed {
    fn read_pos(&self) -> u32 {
        self.read_pos
    }
    fn chrom_id(&self) -> u32 {
        self.chrom_id
    }
    fn ref_pos(&self) -> u32 {
        self.ref_pos_and_strand().0
    }
    fn is_reverse(&self) -> bool {
        self.ref_pos_and_strand().1
    }
    fn length(&self, _k: usize) -> usize {
        _k
    }

    /// Diagonal index — constant along a gapless match.
    /// Forward: `ref_pos - read_pos` (ref increases with read_pos).
    /// Reverse: `ref_pos + read_pos` (ref decreases as read_pos increases).
    fn diagonal(&self) -> i64 {
        if self.is_reverse() {
            self.ref_pos() as i64 + self.read_pos() as i64
        } else {
            self.ref_pos() as i64 - self.read_pos() as i64
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
    type Memento = ();

    fn gap_to<'a>(
        &'a self,
        other: &'a Self,
        k: usize,
        _memento: &'a mut (),
    ) -> Option<(GapResult, Option<&'a Self>, Option<&'a Self>)> {
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
                prev.gap_to(curr, k, &mut Default::default()),
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
            let si_read_end = si.read_pos() + k as u32;

            for rank2 in rank + 1..n {
                let j = by_read[rank2];
                let sj = &self.hits[j];

                // Once the read gap exceeds max_neighbour_gap no further
                // neighbour can produce a non-SV edge.
                let read_gap = sj.read_pos() as i64 - si_read_end as i64;
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
}

impl<'a> CompoundSeed<'a> {
    pub fn new(atoms: &'a [AtomicSeed]) -> Self {
        debug_assert!(!atoms.is_empty());
        Self { atoms }
    }

    pub fn atoms(&self) -> &[AtomicSeed] {
        self.atoms
    }
}

impl<'a> Seed for CompoundSeed<'a> {
    fn read_pos(&self) -> u32 {
        self.atoms[0].read_pos()
    }
    fn chrom_id(&self) -> u32 {
        self.atoms[0].chrom_id()
    }
    fn ref_pos(&self) -> u32 {
        self.atoms[0].ref_pos()
    }
    fn is_reverse(&self) -> bool {
        self.atoms[0].is_reverse()
    }

    fn diagonal(&self) -> i64 {
        self.atoms[0].diagonal()
    }

    fn length(&self, k: usize) -> usize {
        let last = self.atoms.last().unwrap();
        (last.read_pos() + k as u32 - self.atoms[0].read_pos()) as usize
    }
    
    fn to_string(&self, k: usize) -> String {
        let mut result = String::new();
        let mut covered_up_to: u32 = 0;
        for atom in self.atoms {
            let atom_str = atom.to_string(k);
            let atom_start = atom.read_pos();
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
        self.atoms.iter().map(|a| a.atom_weight()).sum()
    }
}

impl<'a> GapComputable for CompoundSeed<'a> {
    type Memento = (Option<CompoundSeed<'a>>, Option<CompoundSeed<'a>>); // (lhs_trimmed, rhs_trimmed)

    fn gap_to<'b>(
        &'b self,
        other: &'b CompoundSeed<'a>,
        k: usize,
        memento: &'b mut Self::Memento,
    ) -> Option<(GapResult, Option<&'b CompoundSeed<'a>>, Option<&'b CompoundSeed<'a>>)> {
        let self_read_end = self.read_pos() + self.length(k) as u32;
        let read_gap_signed = other.read_pos() as i64 - self_read_end as i64;

        let sv_break = self.chrom_id() != other.chrom_id() || self.is_reverse() != other.is_reverse();

        if sv_break && read_gap_signed >= 0 {
            // Cross-chrom/strand with no read overlap: no trimming needed.
            return Some((
                GapResult {
                    ref_gap: i64::MIN,
                    read_gap: read_gap_signed,
                },
                None,
                None,
            ));
        }

        // Fully consumed (same-strand only — cross-strand overlaps are handled below).
        if !sv_break && -read_gap_signed >= other.length(k) as i64 {
            return None;
        }

        if read_gap_signed >= 0 {
            // No overlap: gap is between the last atom of lhs and first atom of rhs.
            let lhs_last = self.atoms.last().unwrap();
            let rhs_first = &other.atoms[0];
            let lhs_end = lhs_last.read_pos() + k as u32;
            let read_gap = rhs_first.read_pos() as i64 - lhs_end as i64;
            let ref_gap = lhs_last.ref_gap_to(rhs_first, k);
            return Some((GapResult { ref_gap, read_gap }, None, None));
        }

        // Overlap: find atoms from each side in the overlap region.
        let lhs_overlap_start = self
            .atoms
            .partition_point(|a| a.read_pos() + k as u32 <= other.read_pos());
        let lhs_overlap = &self.atoms[lhs_overlap_start..];
        let rhs_overlap_end = other
            .atoms
            .partition_point(|a| a.read_pos() < self_read_end);
        let rhs_overlap = &other.atoms[..rhs_overlap_end];

        // Candidate split points p: lhs retains atoms with read_pos+k <= p,
        // rhs retains atoms with read_pos >= p.
        let mut candidates: Vec<u32> = Vec::with_capacity(lhs_overlap.len() + rhs_overlap.len());
        for a in lhs_overlap {
            candidates.push(a.read_pos() + k as u32);
        }
        for a in rhs_overlap {
            candidates.push(a.read_pos());
        }
        candidates.sort_unstable();
        candidates.dedup();

        let mut best: Option<(f64, u32, GapResult)> = None; // (total_trimmed, split_p, gap)

        for &p in &candidates {
            let lhs_trimmed: f64 = lhs_overlap
                .iter()
                .filter(|a| a.read_pos() + k as u32 > p)
                .map(|a| a.atom_weight())
                .sum();
            let rhs_trimmed: f64 = rhs_overlap
                .iter()
                .filter(|a| a.read_pos() < p)
                .map(|a| a.atom_weight())
                .sum();
            let total_trimmed = lhs_trimmed + rhs_trimmed;

            let lhs_effective_last = self
                .atoms
                .iter()
                .rev()
                .find(|a| a.read_pos() + k as u32 <= p);
            let rhs_effective_first = other.atoms.iter().find(|a| a.read_pos() >= p);

            let (ref_gap, read_gap) = if sv_break {
                // For SV breaks, ref_gap is always the sentinel; only read_gap matters.
                let read_gap = match (lhs_effective_last, rhs_effective_first) {
                    (Some(l), Some(r)) => r.read_pos() as i64 - (l.read_pos() + k as u32) as i64,
                    (None, Some(r)) => r.read_pos() as i64 - self.atoms[0].read_pos() as i64,
                    (Some(l), None) => {
                        other.atoms.last().unwrap().read_pos() as i64
                            - (l.read_pos() + k as u32) as i64
                    }
                    (None, None) => continue,
                };
                (i64::MIN, read_gap)
            } else {
                match (lhs_effective_last, rhs_effective_first) {
                    (Some(l), Some(r)) => {
                        let read_gap = r.read_pos() as i64 - (l.read_pos() + k as u32) as i64;
                        (l.ref_gap_to(r, k), read_gap)
                    }
                    (None, Some(r)) => {
                        let l = &self.atoms[0];
                        (
                            l.ref_gap_to(r, k),
                            r.read_pos() as i64 - l.read_pos() as i64,
                        )
                    }
                    (Some(l), None) => {
                        let r = other.atoms.last().unwrap();
                        (
                            l.ref_gap_to(r, k),
                            r.read_pos() as i64 - (l.read_pos() + k as u32) as i64,
                        )
                    }
                    (None, None) => continue,
                }
            };

            let is_better = best
                .as_ref()
                .map_or(true, |&(best_t, ..)| total_trimmed < best_t);
            if is_better {
                best = Some((total_trimmed, p, GapResult { ref_gap, read_gap }));
            }
        }

        best.map(|(_, best_p, gap)| {
            // lhs retains atoms with read_pos + k <= best_p (a prefix of self.atoms).
            let lhs_split = self.atoms.partition_point(|a| a.read_pos() + k as u32 <= best_p);
            // rhs retains atoms with read_pos >= best_p (a suffix of other.atoms).
            let rhs_split = other.atoms.partition_point(|a| a.read_pos() < best_p);

            // Store the retained (trimmed-down) seed only when trimming actually
            // removed some atoms.  None means the full seed is used unchanged.
            *memento = (
                if lhs_split < self.atoms.len() && lhs_split > 0 {
                    Some(CompoundSeed::new(&self.atoms[..lhs_split]))
                } else {
                    None
                },
                if rhs_split > 0 && rhs_split < other.atoms.len() {
                    Some(CompoundSeed::new(&other.atoms[rhs_split..]))
                } else {
                    None
                },
            );
            (gap, memento.0.as_ref(), memento.1.as_ref())
        })
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
        let mut memento = S::Memento::default();
        let (gap, lhs_trimmed, rhs_trimmed) = lhs.gap_to(rhs, k, &mut memento)?;
        if gap.ref_gap == i64::MIN {
            let read_gap_cost = self.read_gap_cost(gap.read_gap);
            return Some((read_gap_cost + self.cfg.sv_penalty, EdgeType::SvBreak));
        }
        let weight_trimmed =
            lhs_trimmed.map(|s| lhs.weight() - s.weight()).unwrap_or(0.0)
            + rhs_trimmed.map(|s| rhs.weight() - s.weight()).unwrap_or(0.0);
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
    order.sort_unstable_by_key(|&i| seeds[i].read_pos());

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
    order.sort_unstable_by_key(|&i| seeds[i].read_pos());

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
