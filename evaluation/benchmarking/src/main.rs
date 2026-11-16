//! Replay benchmark for parallax's bits-scored chaining DP (`best_chain`).
//!
//! Reads a dump of real `best_chain` inputs (written by a parallax build with
//! `instrumentation.patch` applied, via PARALLAX_DUMP_DP) and runs several
//! implementations over them, checking that they agree and timing them.  Gap
//! alignments are replayed from the scores logged in the dump.
//!
//!   A  orig      copy of segments.rs `best_chain` (f64, AoS)
//!   D  skip      A + skip cost evaluation for predecessors the prefix-max
//!                bound already rules out (exact: identical results to A)
//!   B  blocked   SoA, i32 coords, f32 scores, polynomial log2, 8-wide
//!                branch-free blocks written for autovectorisation
//!   C  fearless  same algorithm as B with explicit fearless_simd f32x8/i32x8
//!
//! B and C come in two flavours: without (…0) and with (…1) a block-level
//! version of D's skip.  All variants skip gap alignments that could not win
//! even if perfect (`gap_bits_bound`), as `best_chain` does.
//!
//! Other modes, selected by environment variable:
//!
//!   DPBENCH_FX=1              look gaps up in an Fx-hashed map instead of
//!                             SipHash (parallax's gap cache uses SipHash)
//!   DPBENCH_GAPSPAN=1         how many gap lookups span another seed
//!   DPBENCH_GAPSTUDY=<log>    gap sizes, cheap proxies for the alignment
//!                             score, and their effect on the DP; <log> is
//!                             the PARALLAX_GAPLOG file from the same run

use std::collections::HashMap;
use std::hint::black_box;
use std::time::Instant;

mod gapstudy;
mod simd;

#[derive(Clone, Copy, Debug)]
pub struct Ext {
    pub read_start: i64,
    pub read_end: i64,
    pub entry: i64,
    pub exit: i64,
    pub diagonal: i64,
    pub chrom_id: u32,
    pub is_reverse: bool,
}

pub struct Call {
    /// Index of the read this call belongs to (format-2 dumps only).
    pub read_idx: u32,
    pub ext: Vec<Ext>,
    /// Gap-alignment results observed in the real run, keyed by (q, r) rank.
    /// std SipHash, as parallax's `gap_cache` uses.
    pub gaps: HashMap<(usize, usize), Option<f64>>,
    /// The same, with an Fx-style multiplicative hash.
    pub gaps_fx: HashMap<(usize, usize), Option<f64>, std::hash::BuildHasherDefault<FxHasher>>,
}

/// rustc-hash's FxHasher algorithm (inlined to avoid a dependency).
#[derive(Default, Clone, Copy)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.hash = (self.hash.rotate_left(5) ^ w).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

impl std::hash::Hasher for FxHasher {
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut b = [0u8; 8];
            b[..chunk.len()].copy_from_slice(chunk);
            self.add(u64::from_le_bytes(b));
        }
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    fn finish(&self) -> u64 {
        self.hash
    }
}

fn load(path: &str, limit: usize) -> Vec<Call> {
    let data = std::fs::read(path).expect("read dump");
    let mut p = 0usize;
    let u32_at = |p: &mut usize| {
        let v = u32::from_le_bytes(data[*p..*p + 4].try_into().unwrap());
        *p += 4;
        v
    };
    let i64_at = |p: &mut usize| {
        let v = i64::from_le_bytes(data[*p..*p + 8].try_into().unwrap());
        *p += 8;
        v
    };
    let mut calls = Vec::new();
    // Format 2 starts with "DPD2" and prefixes each record with the read index.
    let v2 = data.len() >= 4 && &data[0..4] == b"DPD2";
    if v2 {
        p = 4;
    }
    while p < data.len() && calls.len() < limit {
        let read_idx = if v2 { u32_at(&mut p) } else { 0 };
        let n = u32_at(&mut p) as usize;
        let mut ext = Vec::with_capacity(n);
        for _ in 0..n {
            let read_start = i64_at(&mut p);
            let read_end = i64_at(&mut p);
            let entry = i64_at(&mut p);
            let exit = i64_at(&mut p);
            let diagonal = i64_at(&mut p);
            let chrom_id = u32_at(&mut p);
            let is_reverse = u32_at(&mut p) != 0;
            ext.push(Ext { read_start, read_end, entry, exit, diagonal, chrom_id, is_reverse });
        }
        let g = u32_at(&mut p) as usize;
        let mut gaps = HashMap::with_capacity(g);
        for _ in 0..g {
            let r = u32_at(&mut p) as usize;
            let q = u32_at(&mut p) as usize;
            let b = f64::from_bits(i64_at(&mut p) as u64);
            gaps.insert((q, r), if b.is_nan() { None } else { Some(b) });
        }
        let gaps_fx = gaps.iter().map(|(&k, &v)| (k, v)).collect();
        calls.push(Call { read_idx, ext, gaps, gaps_fx });
    }
    calls
}

// ── Configuration (ChainingConfig defaults) ─────────────────────────────────

#[derive(Clone, Copy)]
pub struct Cfg {
    pub match_bits: f64,
    pub mismatch_bits: f64,
    pub gap_open_bits: f64,
    pub gap_extend_bits: f64,
    pub sv_break_bits: f64,
    pub gap_align_min: usize,
    pub gap_align_max: usize,
}

pub const CFG: Cfg = Cfg {
    match_bits: 2.0,
    mismatch_bits: 8.0,
    gap_open_bits: 6.0,
    gap_extend_bits: 1.0,
    sv_break_bits: 35.0,
    gap_align_min: 20,
    gap_align_max: 2000,
};

const GENOME_BITS: f64 = 31.5;
const GAP_ALIGN_PREDECESSORS: usize = 3;

/// Upper bound on a gap alignment's bits, as segments.rs `gap_bits_bound`:
/// every base of the shorter side matches, plus one indel for the difference.
pub fn gap_bits_bound(read_gap: i64, ref_gap: i64, cfg: &Cfg) -> f64 {
    let diff = read_gap.abs_diff(ref_gap);
    let indel = if diff == 0 { 0.0 } else { cfg.gap_open_bits + cfg.gap_extend_bits * diff as f64 };
    cfg.match_bits * read_gap.min(ref_gap) as f64 - indel
}

/// The same configuration, pre-converted for the f32 kernels.
#[derive(Clone, Copy)]
pub struct CfgF32 {
    pub mb: f32,
    pub mismatch: f32,
    pub open: f32,
    pub extend: f32,
    pub sv: f32,
    pub sv1: f32,
    pub svg: f32,
    pub gmin: i32,
    pub gmax: i32,
}

impl CfgF32 {
    pub fn from(c: &Cfg) -> Self {
        Self {
            mb: c.match_bits as f32,
            mismatch: c.mismatch_bits as f32,
            open: c.gap_open_bits as f32,
            extend: c.gap_extend_bits as f32,
            sv: c.sv_break_bits as f32,
            sv1: (c.sv_break_bits + 1.0) as f32,
            svg: (c.sv_break_bits + GENOME_BITS) as f32,
            gmin: c.gap_align_min as i32,
            gmax: c.gap_align_max as i32,
        }
    }

    /// `gap_bits_bound` in f32.
    pub fn gap_bound(&self, read_gap: i32, ref_gap: i32) -> f32 {
        let diff = read_gap.abs_diff(ref_gap);
        let indel = if diff == 0 { 0.0 } else { self.open + self.extend * diff as f32 };
        self.mb * read_gap.min(ref_gap) as f32 - indel
    }
}

// ── A: copy of segments.rs best_chain ────────────────────────────────────────

fn best_chain_orig(
    ext: &[Ext],
    order: &[usize],
    mask: &[bool],
    cfg: &Cfg,
    gap: &mut (impl FnMut(usize, usize) -> Option<f64> + ?Sized),
) -> Option<(f64, Vec<usize>)> {
    let order: Vec<usize> = order.iter().copied().filter(|&s| mask[s]).collect();
    let n = order.len();
    let mut dp = vec![0.0f64; n];
    let mut prev = vec![usize::MAX; n];
    // prefix_max[q] = max(dp[0..=q]).
    let mut prefix_max = vec![0.0f64; n];

    for r in 0..n {
        let i = &ext[order[r]];
        let coverage = cfg.match_bits * (i.read_end - i.read_start) as f64;
        dp[r] = coverage;
        let mut nearest: Vec<usize> = Vec::with_capacity(GAP_ALIGN_PREDECESSORS);
        for q in (0..r).rev() {
            let j = &ext[order[q]];
            let in_gap_window = nearest.len() < GAP_ALIGN_PREDECESSORS
                && j.read_end >= i.read_start - cfg.gap_align_max as i64;
            if !in_gap_window && prefix_max[q] + coverage <= dp[r] {
                break;
            }
            if j.read_end >= i.read_end {
                continue;
            }
            let new_bases = i.read_end - i.read_start.max(j.read_end);
            let read_gap = i.read_start - j.read_end;
            let cost = if j.chrom_id != i.chrom_id {
                cfg.sv_break_bits + GENOME_BITS
            } else if j.is_reverse != i.is_reverse {
                cfg.sv_break_bits + 1.0 + (1.0 + (i.entry - j.exit).abs() as f64).log2()
            } else {
                let shift = if i.is_reverse {
                    j.diagonal - i.diagonal
                } else {
                    i.diagonal - j.diagonal
                }
                .abs();
                let indel = match (shift, read_gap > 0) {
                    (0, true) => cfg.mismatch_bits,
                    (0, false) => 0.0,
                    _ => cfg.gap_open_bits + cfg.gap_extend_bits * shift as f64,
                };
                let ref_gap = if i.is_reverse { j.exit - i.entry } else { i.entry - j.exit };
                if nearest.len() < GAP_ALIGN_PREDECESSORS
                    && read_gap >= cfg.gap_align_min as i64
                    && read_gap <= cfg.gap_align_max as i64
                    && ref_gap >= 0
                    && ref_gap <= cfg.gap_align_max as i64
                {
                    nearest.push(q);
                }
                indel.min(cfg.sv_break_bits + (1.0 + shift as f64).log2())
            };
            let candidate = dp[q] + cfg.match_bits * new_bases as f64 - cost;
            if candidate > dp[r] {
                dp[r] = candidate;
                prev[r] = q;
            }
        }
        for &q in &nearest {
            let j = &ext[order[q]];
            let read_gap = i.read_start - j.read_end;
            let ref_gap = if i.is_reverse { j.exit - i.entry } else { i.entry - j.exit };
            if dp[q] + coverage + gap_bits_bound(read_gap, ref_gap, cfg) <= dp[r] {
                continue;
            }
            let Some(bits) = gap(order[q], order[r]) else {
                continue;
            };
            let candidate = dp[q] + coverage + bits;
            if candidate > dp[r] {
                dp[r] = candidate;
                prev[r] = q;
            }
        }
        prefix_max[r] = if r == 0 { dp[r] } else { prefix_max[r - 1].max(dp[r]) };
    }

    let best = (0..n).max_by(|&a, &b| dp[a].total_cmp(&dp[b]))?;
    let mut ranks = Vec::new();
    let mut cur = best;
    loop {
        ranks.push(cur);
        if prev[cur] == usize::MAX {
            break;
        }
        cur = prev[cur];
    }
    ranks.reverse();
    Some((dp[best], ranks.iter().map(|&r| order[r]).collect()))
}

// ── D: original + exact skip of bound-excluded predecessors ─────────────────

fn best_chain_skip(
    ext: &[Ext],
    order: &[usize],
    mask: &[bool],
    cfg: &Cfg,
    gap: &mut (impl FnMut(usize, usize) -> Option<f64> + ?Sized),
) -> Option<(f64, Vec<usize>)> {
    let order: Vec<usize> = order.iter().copied().filter(|&s| mask[s]).collect();
    let n = order.len();
    let mut dp = vec![0.0f64; n];
    let mut prev = vec![usize::MAX; n];
    let mut prefix_max = vec![0.0f64; n];
    let gmin = cfg.gap_align_min as i64;
    let gmax = cfg.gap_align_max as i64;

    for r in 0..n {
        let i = &ext[order[r]];
        let coverage = cfg.match_bits * (i.read_end - i.read_start) as f64;
        dp[r] = coverage;
        let mut nearest = [0usize; GAP_ALIGN_PREDECESSORS];
        let mut nn = 0;
        let mut q = r;
        // Phase 1: full evaluation until the bound rules out every earlier seed.
        while q > 0 {
            q -= 1;
            if prefix_max[q] + coverage <= dp[r] {
                q += 1;
                break;
            }
            let j = &ext[order[q]];
            if j.read_end >= i.read_end {
                continue;
            }
            let new_bases = i.read_end - i.read_start.max(j.read_end);
            let read_gap = i.read_start - j.read_end;
            let cost = if j.chrom_id != i.chrom_id {
                cfg.sv_break_bits + GENOME_BITS
            } else if j.is_reverse != i.is_reverse {
                cfg.sv_break_bits + 1.0 + (1.0 + (i.entry - j.exit).abs() as f64).log2()
            } else {
                let shift = (i.diagonal - j.diagonal).abs();
                let indel = match (shift, read_gap > 0) {
                    (0, true) => cfg.mismatch_bits,
                    (0, false) => 0.0,
                    _ => cfg.gap_open_bits + cfg.gap_extend_bits * shift as f64,
                };
                let ref_gap = if i.is_reverse { j.exit - i.entry } else { i.entry - j.exit };
                if nn < GAP_ALIGN_PREDECESSORS
                    && read_gap >= gmin
                    && read_gap <= gmax
                    && ref_gap >= 0
                    && ref_gap <= gmax
                {
                    nearest[nn] = q;
                    nn += 1;
                }
                indel.min(cfg.sv_break_bits + (1.0 + shift as f64).log2())
            };
            let candidate = dp[q] + cfg.match_bits * new_bases as f64 - cost;
            if candidate > dp[r] {
                dp[r] = candidate;
                prev[r] = q;
            }
        }
        // Phase 2: no earlier seed can win; only look for gap-alignment
        // candidates still inside the window.
        let win_lo = i.read_start - gmax;
        while q > 0 && nn < GAP_ALIGN_PREDECESSORS {
            q -= 1;
            let j = &ext[order[q]];
            if j.read_end < win_lo {
                break;
            }
            if j.read_end >= i.read_end || j.chrom_id != i.chrom_id || j.is_reverse != i.is_reverse {
                continue;
            }
            let read_gap = i.read_start - j.read_end;
            let ref_gap = if i.is_reverse { j.exit - i.entry } else { i.entry - j.exit };
            if read_gap >= gmin && read_gap <= gmax && ref_gap >= 0 && ref_gap <= gmax {
                nearest[nn] = q;
                nn += 1;
            }
        }
        for &q in &nearest[..nn] {
            // Could this gap alignment possibly win?
            let j = &ext[order[q]];
            let read_gap = i.read_start - j.read_end;
            let ref_gap = if i.is_reverse { j.exit - i.entry } else { i.entry - j.exit };
            let prunable = dp[q] + coverage + gap_bits_bound(read_gap, ref_gap, cfg) <= dp[r];
            // Does another seed on the same chrom/strand lie wholly inside
            // this gap, collinear in both read and reference?
            let spanned = SPAN_MODE.with(|m| *m.borrow())
                && order.iter().any(|&m| {
                    let e = &ext[m];
                    m != order[q]
                        && m != order[r]
                        && e.chrom_id == i.chrom_id
                        && e.is_reverse == i.is_reverse
                        && e.read_start >= j.read_end
                        && e.read_end <= i.read_start
                        && if i.is_reverse {
                            e.entry <= j.exit && e.exit >= i.entry
                        } else {
                            e.entry >= j.exit && e.exit <= i.entry
                        }
                });
            GAP_STATS.with(|s| {
                let mut s = s.borrow_mut();
                s.0 += 1;
                s.1 += prunable as u64;
                s.3 += spanned as u64;
            });
            if prunable {
                continue;
            }
            let Some(bits) = gap(order[q], order[r]) else {
                continue;
            };
            let candidate = dp[q] + coverage + bits;
            if candidate > dp[r] {
                dp[r] = candidate;
                prev[r] = q;
                GAP_STATS.with(|s| {
                    let mut s = s.borrow_mut();
                    s.2 += 1;
                    s.4 += spanned as u64;
                });
            }
        }
        prefix_max[r] = if r == 0 { dp[r] } else { prefix_max[r - 1].max(dp[r]) };
    }

    let best = (0..n).max_by(|&a, &b| dp[a].total_cmp(&dp[b]))?;
    let mut ranks = Vec::new();
    let mut cur = best;
    loop {
        ranks.push(cur);
        if prev[cur] == usize::MAX {
            break;
        }
        cur = prev[cur];
    }
    ranks.reverse();
    Some((dp[best], ranks.iter().map(|&r| order[r]).collect()))
}

thread_local! {
    /// (gap lookups, lookups an upper bound rules out, lookups that improved dp[r],
    ///  lookups whose gap contains another same-locus seed, improving lookups that did)
    pub static GAP_STATS: std::cell::RefCell<(u64, u64, u64, u64, u64)> = const { std::cell::RefCell::new((0, 0, 0, 0, 0)) };
    /// Compute the (slow) `spanned` statistic; off during timed runs.
    pub static SPAN_MODE: std::cell::RefCell<bool> = const { std::cell::RefCell::new(false) };
}

// ── Shared SoA layout for B and C ───────────────────────────────────────────

/// Lanes per block.
pub const W: usize = 8;
/// Padding before the first seed, so a block's loads never go out of bounds.
pub const PAD: usize = W;

/// Struct-of-arrays copy of the selected seeds in rank order, with `PAD`
/// sentinel entries in front (never valid predecessors).
#[derive(Default)]
pub struct Soa {
    pub n: usize,
    pub idx: Vec<usize>,
    pub rs: Vec<i32>,
    pub re: Vec<i32>,
    pub entry: Vec<i32>,
    pub exit: Vec<i32>,
    pub diag: Vec<i32>,
    /// chrom_id * 2 + is_reverse
    pub key: Vec<i32>,
    pub dp: Vec<f32>,
    pub pm: Vec<f32>,
    pub prev: Vec<u32>,
}

impl Soa {
    /// The marshalling step: gather the masked seeds into SoA buffers.
    pub fn fill(&mut self, ext: &[Ext], order: &[usize], mask: &[bool]) {
        self.idx.clear();
        self.idx.extend(order.iter().copied().filter(|&s| mask[s]));
        let n = self.idx.len();
        self.n = n;
        let len = PAD + n;
        for v in [&mut self.rs, &mut self.re, &mut self.entry, &mut self.exit, &mut self.diag, &mut self.key] {
            v.clear();
            v.resize(len, 0);
        }
        self.dp.clear();
        self.dp.resize(len, f32::NEG_INFINITY);
        self.pm.clear();
        self.pm.resize(len, f32::NEG_INFINITY);
        self.prev.clear();
        self.prev.resize(len, u32::MAX);
        for p in 0..PAD {
            // read_end = i32::MAX makes a sentinel fail `read_end < i.read_end`.
            self.re[p] = i32::MAX;
            self.key[p] = -2;
        }
        for (r, &s) in self.idx.iter().enumerate() {
            let e = &ext[s];
            let p = PAD + r;
            self.rs[p] = e.read_start as i32;
            self.re[p] = e.read_end as i32;
            self.entry[p] = e.entry as i32;
            self.exit[p] = e.exit as i32;
            self.diag[p] = e.diagonal as i32;
            self.key[p] = (e.chrom_id as i32) * 2 + e.is_reverse as i32;
        }
    }

    /// Best endpoint and traceback, as in the original.
    pub fn traceback(&self) -> Option<(f32, Vec<usize>)> {
        let n = self.n;
        let dp = &self.dp[PAD..PAD + n];
        let best = (0..n).max_by(|&a, &b| dp[a].total_cmp(&dp[b]))?;
        let mut ranks = Vec::new();
        let mut cur = best;
        loop {
            ranks.push(cur);
            let p = self.prev[PAD + cur];
            if p == u32::MAX {
                break;
            }
            cur = p as usize;
        }
        ranks.reverse();
        Some((dp[best], ranks.iter().map(|&r| self.idx[r]).collect()))
    }
}

/// log2(1 + t) on [0, 1] as t + t(1 - t)(a0 + a1 t + a2 t² + a3 t³);
/// least-squares fit, max error 1.8e-5, exact at t = 0 and t = 1.
pub const L2A: [f32; 4] = [0.441_835_87, -0.266_536_28, 0.147_018_77, -0.044_243_74];

/// log2(1 + x) for x ≥ 0, via exponent extraction and the polynomial above.
#[inline(always)]
pub fn log2_1p(x: f32) -> f32 {
    let y = 1.0 + x;
    let bits = y.to_bits();
    let e = ((bits >> 23) as i32 - 127) as f32;
    let t = f32::from_bits((bits & 0x007f_ffff) | 0x3f80_0000) - 1.0;
    let p = t + t * (1.0 - t) * (L2A[0] + t * (L2A[1] + t * (L2A[2] + t * L2A[3])));
    e + p
}

// ── B: blocked SoA kernel written for autovectorisation ─────────────────────

pub fn best_chain_blocked<const SKIP: bool>(
    ext: &[Ext],
    order: &[usize],
    mask: &[bool],
    cfg: &CfgF32,
    soa: &mut Soa,
    gap: &mut (impl FnMut(usize, usize) -> Option<f64> + ?Sized),
) -> Option<(f32, Vec<usize>)> {
    soa.fill(ext, order, mask);
    let n = soa.n;
    for r in 0..n {
        let ri = PAD + r;
        let (irs, ire, ient, idiag, ikey) = (soa.rs[ri], soa.re[ri], soa.entry[ri], soa.diag[ri], soa.key[ri]);
        let ichrom = ikey >> 1;
        let irev = ikey & 1 == 1;
        let coverage = cfg.mb * (ire - irs) as f32;
        let mut best = coverage;
        let mut best_q = u32::MAX;
        let mut nearest = [0usize; GAP_ALIGN_PREDECESSORS];
        let mut nn = 0usize;
        let win_lo = irs - cfg.gmax;
        let mut hi = ri;
        while hi > PAD {
            let top = hi - 1;
            let in_win = nn < GAP_ALIGN_PREDECESSORS && soa.re[top] >= win_lo;
            let hopeless = soa.pm[top] + coverage <= best;
            if hopeless && !in_win {
                break;
            }
            let lo = hi - W;
            let re = &soa.re[lo..hi];
            let ex = &soa.exit[lo..hi];
            let ky = &soa.key[lo..hi];
            let mut elig = [false; W];
            if !(SKIP && hopeless) {
                let dg = &soa.diag[lo..hi];
                let dpv = &soa.dp[lo..hi];
                let mut cand = [0f32; W];
                for l in 0..W {
                    let (rej, exj, dgj, kj, dpj) = (re[l], ex[l], dg[l], ky[l], dpv[l]);
                    let valid = rej < ire;
                    let same_strand = kj == ikey;
                    let same_chrom = (kj >> 1) == ichrom;
                    let d = ient.wrapping_sub(exj);
                    let dist = d.wrapping_abs();
                    let shift = idiag.wrapping_sub(dgj).wrapping_abs();
                    let x = if same_strand { shift } else { dist };
                    let lg = log2_1p(x as f32);
                    let read_gap = irs.wrapping_sub(rej);
                    let indel_nz = cfg.open + cfg.extend * shift as f32;
                    let indel_z = if read_gap > 0 { cfg.mismatch } else { 0.0 };
                    let indel = if shift == 0 { indel_z } else { indel_nz };
                    let c_same = indel.min(cfg.sv + lg);
                    let c_strand = cfg.sv1 + lg;
                    let cost = if same_strand { c_same } else if same_chrom { c_strand } else { cfg.svg };
                    let new_bases = ire.wrapping_sub(irs.max(rej));
                    let c = dpj + cfg.mb * new_bases as f32 - cost;
                    cand[l] = if valid { c } else { f32::NEG_INFINITY };
                    let ref_gap = if irev { d.wrapping_neg() } else { d };
                    elig[l] = valid
                        & same_strand
                        & (read_gap >= cfg.gmin)
                        & (read_gap <= cfg.gmax)
                        & (ref_gap >= 0)
                        & (ref_gap <= cfg.gmax);
                }
                // Highest lane wins ties, matching the original's descending scan.
                let mut bm = f32::NEG_INFINITY;
                let mut bl = 0usize;
                for l in (0..W).rev() {
                    if cand[l] > bm {
                        bm = cand[l];
                        bl = l;
                    }
                }
                if bm > best {
                    best = bm;
                    best_q = (lo + bl - PAD) as u32;
                }
            } else {
                for l in 0..W {
                    let (rej, exj, kj) = (re[l], ex[l], ky[l]);
                    let d = ient.wrapping_sub(exj);
                    let read_gap = irs.wrapping_sub(rej);
                    let ref_gap = if irev { d.wrapping_neg() } else { d };
                    elig[l] = (rej < ire)
                        & (kj == ikey)
                        & (read_gap >= cfg.gmin)
                        & (read_gap <= cfg.gmax)
                        & (ref_gap >= 0)
                        & (ref_gap <= cfg.gmax);
                }
            }
            if nn < GAP_ALIGN_PREDECESSORS {
                for l in (0..W).rev() {
                    if elig[l] {
                        nearest[nn] = lo + l - PAD;
                        nn += 1;
                        if nn == GAP_ALIGN_PREDECESSORS {
                            break;
                        }
                    }
                }
            }
            hi = lo;
        }
        for &q in &nearest[..nn] {
            let qi = PAD + q;
            let read_gap = irs - soa.re[qi];
            let ref_gap = if irev { soa.exit[qi] - ient } else { ient - soa.exit[qi] };
            if soa.dp[qi] + coverage + cfg.gap_bound(read_gap, ref_gap) <= best {
                continue;
            }
            let Some(bits) = gap(soa.idx[q], soa.idx[r]) else {
                continue;
            };
            let candidate = soa.dp[qi] + coverage + bits as f32;
            if candidate > best {
                best = candidate;
                best_q = q as u32;
            }
        }
        soa.dp[ri] = best;
        soa.prev[ri] = best_q;
        soa.pm[ri] = soa.pm[ri - 1].max(best);
    }
    soa.traceback()
}

// ── Harness ──────────────────────────────────────────────────────────────────

struct Outcome {
    score: f64,
    chain: Vec<usize>,
}

fn run_variant(
    name: &str,
    calls: &[Call],
    reps: usize,
    fx: bool,
    mut f: impl FnMut(usize, &Call, &mut dyn FnMut(usize, usize) -> Option<f64>) -> Option<(f64, Vec<usize>)>,
) -> (f64, Vec<Option<Outcome>>, u64) {
    let mut best_secs = f64::INFINITY;
    let mut outcomes = Vec::new();
    let mut misses = 0u64;
    for rep in 0..reps {
        let mut out = Vec::with_capacity(calls.len());
        let t0 = Instant::now();
        for (k, call) in calls.iter().enumerate() {
            let mut gap = |j: usize, i: usize| -> Option<f64> {
                let hit = if fx { call.gaps_fx.get(&(j, i)) } else { call.gaps.get(&(j, i)) };
                match hit {
                    Some(v) => *v,
                    None => {
                        misses += 1;
                        None
                    }
                }
            };
            let res = f(k, call, &mut gap);
            out.push(black_box(res).map(|(score, chain)| Outcome { score, chain }));
        }
        let secs = t0.elapsed().as_secs_f64();
        best_secs = best_secs.min(secs);
        if rep == 0 {
            outcomes = out;
        }
    }
    eprintln!("  ran {name}");
    (best_secs, outcomes, misses)
}

fn compare(name: &str, base: &[Option<Outcome>], other: &[Option<Outcome>]) -> String {
    let mut chain_diff = 0usize;
    let mut max_score_diff = 0f64;
    let mut max_rel = 0f64;
    for (a, b) in base.iter().zip(other) {
        match (a, b) {
            (Some(a), Some(b)) => {
                if a.chain != b.chain {
                    chain_diff += 1;
                }
                let d = (a.score - b.score).abs();
                max_score_diff = max_score_diff.max(d);
                max_rel = max_rel.max(d / a.score.abs().max(1.0));
            }
            (None, None) => {}
            _ => chain_diff += 1,
        }
    }
    format!(
        "{name}: chains differing {}/{} ({:.3}%), max |score diff| {:.4} bits (rel {:.1e})",
        chain_diff,
        base.len(),
        100.0 * chain_diff as f64 / base.len() as f64,
        max_score_diff,
        max_rel
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        // Without a dump, just check the SIMD kernel's building blocks.
        simd::self_test();
        eprintln!("usage: dpbench <dump> [reps] [limit]   (see README.md)");
        return;
    };
    let reps: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(5);
    let limit: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    simd::self_test();

    let calls = load(path, limit);
    let total_n: usize = calls.iter().map(|c| c.ext.len()).sum();
    let max_coord = calls
        .iter()
        .flat_map(|c| c.ext.iter())
        .map(|e| e.read_end.max(e.entry).max(e.exit).max(e.diagonal.abs()))
        .max()
        .unwrap_or(0);
    eprintln!(
        "loaded {} calls, {} seeds total (mean n {:.1}); max |coord| {} (i32 max {})",
        calls.len(),
        total_n,
        total_n as f64 / calls.len() as f64,
        max_coord,
        i32::MAX
    );
    let level = fearless_simd::Level::new();
    eprintln!("fearless_simd level: {:?}", level);
    if let Ok(gaplog) = std::env::var("DPBENCH_GAPSTUDY") {
        gapstudy::run(&calls, &gaplog);
        return;
    }
    if std::env::var("DPBENCH_GAPSPAN").is_ok() {
        // One untimed pass of D (exact semantics) collecting gap statistics.
        SPAN_MODE.with(|m| *m.borrow_mut() = true);
        for c in &calls {
            let order: Vec<usize> = (0..c.ext.len()).collect();
            let mask = vec![true; c.ext.len()];
            let mut gap = |j: usize, i: usize| c.gaps.get(&(j, i)).copied().flatten();
            black_box(best_chain_skip(&c.ext, &order, &mask, &CFG, &mut gap));
        }
        GAP_STATS.with(|s| {
            let (n, prunable, improved, spanned, improved_spanned) = *s.borrow();
            let pct = |a: u64, b: u64| 100.0 * a as f64 / b.max(1) as f64;
            println!("gap lookups: {n}");
            println!("  improved dp[r]:                         {:5.1}%", pct(improved, n));
            println!("  ruled out by simple upper bound:        {:5.1}%", pct(prunable, n));
            println!("  gap contains another same-locus seed:   {:5.1}%", pct(spanned, n));
            println!("    of those, improved dp[r]:             {:5.1}%", pct(improved_spanned, spanned));
            println!("  truly seedless gaps:                    {:5.1}%", pct(n - spanned, n));
            println!("    of those, improved dp[r]:             {:5.1}%", pct(improved - improved_spanned, n - spanned));
            println!("  share of all improvements from spanned gaps: {:5.1}%", pct(improved_spanned, improved));
        });
        return;
    }
    let fx = std::env::var("DPBENCH_FX").is_ok();
    println!("gap lookups via {}", if fx { "FxHash map" } else { "std SipHash map (as parallax's gap_cache)" });

    let cfg = CFG;
    let cfg32 = CfgF32::from(&cfg);
    let all_order: Vec<Vec<usize>> = calls.iter().map(|c| (0..c.ext.len()).collect()).collect();
    let all_mask: Vec<Vec<bool>> = calls.iter().map(|c| vec![true; c.ext.len()]).collect();

    let (t_a, o_a, m_a) = run_variant("A orig", &calls, reps, fx, |k, c, gap| {
        best_chain_orig(&c.ext, &all_order[k], &all_mask[k], &cfg, gap)
    });
    let (t_d, o_d, m_d) = run_variant("D skip", &calls, reps, fx, |k, c, gap| {
        best_chain_skip(&c.ext, &all_order[k], &all_mask[k], &cfg, gap)
    });
    let mut soa = Soa::default();
    let (t_b0, o_b0, m_b0) = run_variant("B0 blocked", &calls, reps, fx, |k, c, gap| {
        best_chain_blocked::<false>(&c.ext, &all_order[k], &all_mask[k], &cfg32, &mut soa, gap)
            .map(|(s, ch)| (s as f64, ch))
    });
    let (t_b1, o_b1, m_b1) = run_variant("B1 blocked+skip", &calls, reps, fx, |k, c, gap| {
        best_chain_blocked::<true>(&c.ext, &all_order[k], &all_mask[k], &cfg32, &mut soa, gap)
            .map(|(s, ch)| (s as f64, ch))
    });
    let (t_c0, o_c0, m_c0) = run_variant("C0 fearless", &calls, reps, fx, |k, c, gap| {
        simd::best_chain_fearless(level, false, &c.ext, &all_order[k], &all_mask[k], &cfg32, &mut soa, gap)
            .map(|(s, ch)| (s as f64, ch))
    });
    let (t_c1, o_c1, m_c1) = run_variant("C1 fearless+skip", &calls, reps, fx, |k, c, gap| {
        simd::best_chain_fearless(level, true, &c.ext, &all_order[k], &all_mask[k], &cfg32, &mut soa, gap)
            .map(|(s, ch)| (s as f64, ch))
    });

    // Marshalling alone: the per-call SoA fill (mask filter + 6 i32 arrays).
    let mut t_fill = f64::INFINITY;
    for _ in 0..reps {
        let t0 = Instant::now();
        for (k, c) in calls.iter().enumerate() {
            soa.fill(&c.ext, &all_order[k], &all_mask[k]);
            black_box(&soa);
        }
        t_fill = t_fill.min(t0.elapsed().as_secs_f64());
    }
    // Gap lookups alone: the same (q, r) lookups the DP makes, with nothing else.
    let mut t_look = f64::INFINITY;
    for _ in 0..reps {
        let t0 = Instant::now();
        let mut hits = 0usize;
        for c in &calls {
            if fx {
                for k in c.gaps.keys() {
                    hits += black_box(c.gaps_fx.get(k)).is_some() as usize;
                }
            } else {
                let keys: Vec<_> = c.gaps.keys().copied().collect();
                for k in &keys {
                    hits += black_box(c.gaps.get(k)).is_some() as usize;
                }
            }
        }
        black_box(hits);
        t_look = t_look.min(t0.elapsed().as_secs_f64());
    }

    println!("calls {}  mean n {:.1}  reps {} (best of)", calls.len(), total_n as f64 / calls.len() as f64, reps);
    println!(
        "SoA marshalling alone: {:.2} us/call ({:.1}% of C1);  distinct gap lookups alone: {:.2} us/call",
        t_fill * 1e6 / calls.len() as f64,
        100.0 * t_fill / t_c1,
        t_look * 1e6 / calls.len() as f64
    );
    println!("{:<22} {:>10} {:>12} {:>9}  {}", "variant", "total ms", "us/call", "speedup", "gap-log misses");
    for (name, t, m) in [
        ("A  orig (f64)", t_a, m_a),
        ("D  orig+skip (f64)", t_d, m_d),
        ("B0 blocked f32", t_b0, m_b0),
        ("B1 blocked f32+skip", t_b1, m_b1),
        ("C0 fearless f32x8", t_c0, m_c0),
        ("C1 fearless+skip", t_c1, m_c1),
    ] {
        println!(
            "{:<22} {:>10.1} {:>12.2} {:>8.2}x  {}",
            name,
            t * 1e3,
            t * 1e6 / calls.len() as f64,
            t_a / t,
            m / reps as u64
        );
    }
    println!("{}", compare("D  vs A", &o_a, &o_d));
    println!("{}", compare("B0 vs A", &o_a, &o_b0));
    println!("{}", compare("B1 vs A", &o_a, &o_b1));
    println!("{}", compare("C0 vs A", &o_a, &o_c0));
    println!("{}", compare("C1 vs A", &o_a, &o_c1));
    println!("{}", compare("C0 vs B0", &o_b0, &o_c0));
    println!("{}", compare("C1 vs B1", &o_b1, &o_c1));
    GAP_STATS.with(|s| {
        let s = s.borrow();
        println!(
            "gap lookups per rep {}: {:.1}% provably unable to win (upper bound), {:.1}% actually improved dp[r]",
            s.0 / reps as u64,
            100.0 * s.1 as f64 / s.0.max(1) as f64,
            100.0 * s.2 as f64 / s.0.max(1) as f64
        );
    });
}
