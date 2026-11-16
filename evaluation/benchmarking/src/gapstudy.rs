//! Gap-alignment study: the sizes of the gaps the chaining DP aligns, and how
//! well cheap proxies predict their score in bits.
//!
//! Inputs: a format-2 `best_chain` dump (PARALLAX_DUMP_DP) and the gap log
//! (PARALLAX_GAPLOG) from the same instrumented run.

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use crate::{CFG, Call, Cfg};

pub struct GapRec {
    read_idx: u32,
    qs: u32,
    qe: u32,
    rs: u32,
    re: u32,
    chrom: u32,
    rev: bool,
    ok: bool,
    bits: f64,
    ns: u64,
    m: u32,
    x: u32,
    ins: u32,
    del: u32,
    ops: u32,
    q: Vec<u8>,
    r: Vec<u8>,
}

type Key = (u32, u32, u32, u32, u32, u32, bool);

impl GapRec {
    fn key(&self) -> Key {
        (self.read_idx, self.qs, self.qe, self.rs, self.re, self.chrom, self.rev)
    }
    fn qg(&self) -> usize {
        self.q.len()
    }
    fn rg(&self) -> usize {
        self.r.len()
    }
}

fn load_gaplog(path: &str) -> Vec<GapRec> {
    let data = std::fs::read(path).expect("read gap log");
    let mut p = 0usize;
    let rd32 = |p: &mut usize| {
        let v = u32::from_le_bytes(data[*p..*p + 4].try_into().unwrap());
        *p += 4;
        v
    };
    let rd64 = |p: &mut usize| {
        let v = u64::from_le_bytes(data[*p..*p + 8].try_into().unwrap());
        *p += 8;
        v
    };
    let mut out = Vec::new();
    while p < data.len() {
        let read_idx = rd32(&mut p);
        let (qs, qe, rs, re) = (rd32(&mut p), rd32(&mut p), rd32(&mut p), rd32(&mut p));
        let chrom = rd32(&mut p);
        let rev = rd32(&mut p) != 0;
        let ok = rd32(&mut p) != 0;
        let bits = f64::from_bits(rd64(&mut p));
        let ns = rd64(&mut p);
        let (m, x, ins, del, ops) = (rd32(&mut p), rd32(&mut p), rd32(&mut p), rd32(&mut p), rd32(&mut p));
        let ql = rd32(&mut p) as usize;
        let rl = rd32(&mut p) as usize;
        let q = data[p..p + ql].to_vec();
        p += ql;
        let r = data[p..p + rl].to_vec();
        p += rl;
        out.push(GapRec { read_idx, qs, qe, rs, re, chrom, rev, ok, bits, ns, m, x, ins, del, ops, q, r });
    }
    out
}

// ── Proxies ──────────────────────────────────────────────────────────────────

#[inline]
fn enc(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// k-mer code at each start position (u32::MAX where the k-mer contains an N).
fn codes(s: &[u8], k: usize, out: &mut Vec<u32>) {
    out.clear();
    if s.len() < k {
        return;
    }
    let mask = (1u32 << (2 * k)) - 1;
    let (mut v, mut valid) = (0u32, 0usize);
    for (i, &b) in s.iter().enumerate() {
        let c = enc(b);
        if c > 3 {
            v = 0;
            valid = 0;
        } else {
            v = ((v << 2) | c as u32) & mask;
            valid += 1;
        }
        if i + 1 >= k {
            out.push(if valid >= k { v } else { u32::MAX });
        }
    }
}

#[derive(Default)]
struct Scratch {
    cq: Vec<u32>,
    cr: Vec<u32>,
    cnt: Vec<u16>,
    touched: Vec<u32>,
}

/// k-mer spectrum containment: shared k-mers (multiset) / min(#k-mers).
fn spectrum(q: &[u8], r: &[u8], k: usize, s: &mut Scratch) -> f64 {
    codes(q, k, &mut s.cq);
    codes(r, k, &mut s.cr);
    let size = 1usize << (2 * k);
    if s.cnt.len() < size {
        s.cnt.resize(size, 0);
    }
    let nq = s.cq.iter().filter(|&&c| c != u32::MAX).count();
    let nr = s.cr.iter().filter(|&&c| c != u32::MAX).count();
    if nq == 0 || nr == 0 {
        return 0.0;
    }
    s.touched.clear();
    for &c in &s.cr {
        if c != u32::MAX {
            if s.cnt[c as usize] == 0 {
                s.touched.push(c);
            }
            s.cnt[c as usize] += 1;
        }
    }
    let mut shared = 0usize;
    for &c in &s.cq {
        if c != u32::MAX && s.cnt[c as usize] > 0 {
            s.cnt[c as usize] -= 1;
            shared += 1;
        }
    }
    for &c in &s.touched {
        s.cnt[c as usize] = 0;
    }
    shared as f64 / nq.min(nr) as f64
}

/// Order-aware k-mer score: fraction of Q's k-mers found in R within ±w of
/// the left-anchored (j = i) or right-anchored (j = i + rg - qg) diagonal,
/// so a single indel anywhere is tolerated.
fn banded(q: &[u8], r: &[u8], k: usize, w: i64, s: &mut Scratch) -> f64 {
    codes(q, k, &mut s.cq);
    codes(r, k, &mut s.cr);
    let (nq, nr) = (s.cq.len() as i64, s.cr.len() as i64);
    if nq == 0 || nr == 0 {
        return 0.0;
    }
    let d = r.len() as i64 - q.len() as i64;
    let (mut valid, mut hit) = (0usize, 0usize);
    for i in 0..nq {
        let c = s.cq[i as usize];
        if c == u32::MAX {
            continue;
        }
        valid += 1;
        let mut found = false;
        for center in [i, i + d] {
            let (lo, hi) = ((center - w).max(0), (center + w).min(nr - 1));
            let mut j = lo;
            while j <= hi {
                if s.cr[j as usize] == c {
                    found = true;
                    break;
                }
                j += 1;
            }
            if found {
                break;
            }
        }
        hit += found as usize;
    }
    if valid == 0 { 0.0 } else { hit as f64 / valid as f64 }
}

fn indel_bits(n: usize, cfg: &Cfg) -> f64 {
    if n == 0 { 0.0 } else { cfg.gap_open_bits + cfg.gap_extend_bits * n as f64 }
}

/// Bits of the best alignment with at most one indel (placed optimally),
/// scored exactly as `alignment_bits` would.  O(qg + rg).
fn one_indel_bits(q: &[u8], r: &[u8], cfg: &Cfg) -> f64 {
    let (a, b) = (q.len(), r.len());
    if a == 0 || b == 0 {
        return -indel_bits(a.max(b), cfg);
    }
    let sc = |x: u8, y: u8| {
        let (ex, ey) = (enc(x), enc(y));
        if ex == ey && ex < 4 { cfg.match_bits } else { -cfg.mismatch_bits }
    };
    let m = a.min(b);
    let e = a.max(b) - m;
    let left = |t: usize| sc(q[t], r[t]);
    let right = |t: usize| if b >= a { sc(q[t], r[t + e]) } else { sc(q[t + e], r[t]) };
    let total_right: f64 = (0..m).map(right).sum();
    let (mut best, mut pl, mut sr) = (total_right, 0.0, total_right);
    for t in 0..m {
        pl += left(t);
        sr -= right(t);
        best = best.max(pl + sr);
    }
    best - indel_bits(e, cfg)
}

/// Bits estimate from a k-mer match fraction `f`, given the background
/// fraction `b0` for unrelated sequence: k-mer survival → identity p →
/// match_bits per identical base, mismatch_bits per other, one indel.
fn est_from_fraction(f: f64, b0: f64, k: usize, qg: usize, rg: usize, cfg: &Cfg) -> f64 {
    let surv = ((f - b0) / (1.0 - b0)).clamp(0.0, 1.0);
    let p = surv.powf(1.0 / k as f64);
    let lmin = qg.min(rg) as f64;
    lmin * (cfg.match_bits * p - cfg.mismatch_bits * (1.0 - p)) - indel_bits(qg.abs_diff(rg), cfg)
}

/// Upper bound on any alignment's bits: every base of the shorter side
/// matches, plus one indel for the length difference.
fn upper_bound(qg: usize, rg: usize, cfg: &Cfg) -> f64 {
    cfg.match_bits * qg.min(rg) as f64 - indel_bits(qg.abs_diff(rg), cfg)
}

// ── Helpers ──────────────────────────────────────────────────────────────────

const EDGES: [usize; 7] = [0, 32, 64, 128, 256, 512, 1024];
const LABELS: [&str; 7] = ["<32", "32-63", "64-127", "128-255", "256-511", "512-1023", "1024-2000"];

fn lbin(l: usize) -> usize {
    EDGES.iter().rposition(|&e| l >= e).unwrap()
}

fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.total_cmp(b));
    v[v.len() / 2]
}

/// Mann-Whitney AUC of `score` separating positives from negatives.
fn auc(pairs: &mut [(f64, bool)]) -> f64 {
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (mut rank_sum, mut npos) = (0.0, 0usize);
    let mut i = 0;
    while i < pairs.len() {
        let mut j = i;
        while j < pairs.len() && pairs[j].0 == pairs[i].0 {
            j += 1;
        }
        let avg = (i + j + 1) as f64 / 2.0; // ranks are 1-based
        for p in &pairs[i..j] {
            if p.1 {
                rank_sum += avg;
                npos += 1;
            }
        }
        i = j;
    }
    let nneg = pairs.len() - npos;
    if npos == 0 || nneg == 0 {
        return f64::NAN;
    }
    (rank_sum - (npos * (npos + 1)) as f64 / 2.0) / (npos as f64 * nneg as f64)
}

fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (x, y) in xs.iter().zip(ys) {
        sxy += (x - mx) * (y - my);
        sxx += (x - mx) * (x - mx);
        syy += (y - my) * (y - my);
    }
    sxy / (sxx * syy).sqrt()
}

// ── DP replay with a pluggable gap scorer ───────────────────────────────────

/// `best_chain` (f64, exact) over a dumped call; `gapf(q, r, need)` scores
/// the gap, where the gap candidate beats the current best iff bits > need.
fn replay(call: &Call, cfg: &Cfg, gapf: &mut impl FnMut(usize, usize, f64) -> Option<f64>) -> Option<(f64, Vec<usize>)> {
    let ext = &call.ext;
    let n = ext.len();
    let mut dp = vec![0.0f64; n];
    let mut prev = vec![usize::MAX; n];
    let mut prefix_max = vec![0.0f64; n];
    for r in 0..n {
        let i = &ext[r];
        let coverage = cfg.match_bits * (i.read_end - i.read_start) as f64;
        dp[r] = coverage;
        let mut nearest: Vec<usize> = Vec::with_capacity(3);
        for q in (0..r).rev() {
            let j = &ext[q];
            let in_gap_window = nearest.len() < 3 && j.read_end >= i.read_start - cfg.gap_align_max as i64;
            if !in_gap_window && prefix_max[q] + coverage <= dp[r] {
                break;
            }
            if j.read_end >= i.read_end {
                continue;
            }
            let new_bases = i.read_end - i.read_start.max(j.read_end);
            let read_gap = i.read_start - j.read_end;
            let cost = if j.chrom_id != i.chrom_id {
                cfg.sv_break_bits + 31.5
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
                if nearest.len() < 3
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
            let need = dp[r] - dp[q] - coverage;
            let Some(bits) = gapf(q, r, need) else {
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
    Some((dp[best], ranks))
}

/// The records a chain would be emitted as: maximal runs of `Continuation`
/// edges per `bits_edge_type`, as (chrom, strand, read start, read end, diagonal).
fn records(call: &Call, chain: &[usize], max_bridge: i64) -> Vec<(u32, bool, i64, i64, i64)> {
    let e = &call.ext;
    let mut out = Vec::new();
    let mut start = chain[0];
    for w in chain.windows(2) {
        let (l, r) = (&e[w[0]], &e[w[1]]);
        let read_gap = r.read_start - l.read_end;
        let ref_gap = if l.is_reverse { l.exit - r.entry } else { r.entry - l.exit };
        let cont = l.chrom_id == r.chrom_id
            && l.is_reverse == r.is_reverse
            && ref_gap >= -(r.read_end - r.read_start)
            && read_gap <= max_bridge
            && (ref_gap - read_gap).abs() <= max_bridge;
        if !cont {
            let s = &e[start];
            out.push((s.chrom_id, s.is_reverse, s.read_start, l.read_end, s.diagonal));
            start = w[1];
        }
    }
    let (s, last) = (&e[start], &e[*chain.last().unwrap()]);
    out.push((s.chrom_id, s.is_reverse, s.read_start, last.read_end, s.diagonal));
    out
}

/// Same records within `tol` bp on read boundaries and diagonal.
fn same_records(a: &[(u32, bool, i64, i64, i64)], b: &[(u32, bool, i64, i64, i64)], tol: i64) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(x, y)| {
            x.0 == y.0 && x.1 == y.1 && (x.2 - y.2).abs() <= tol && (x.3 - y.3).abs() <= tol && (x.4 - y.4).abs() <= tol
        })
}

/// Does another seed lie wholly inside the gap from q to r, collinear in
/// read and reference on the same chrom/strand?
fn spans_seed(call: &Call, q: usize, r: usize) -> bool {
    let (j, i) = (&call.ext[q], &call.ext[r]);
    call.ext.iter().enumerate().any(|(m, e)| {
        m != q
            && m != r
            && e.chrom_id == i.chrom_id
            && e.is_reverse == i.is_reverse
            && e.read_start >= j.read_end
            && e.read_end <= i.read_start
            && if i.is_reverse { e.entry <= j.exit && e.exit >= i.entry } else { e.entry >= j.exit && e.exit <= i.entry }
    })
}

fn key_for(call: &Call, q: usize, r: usize) -> Key {
    let (j, i) = (&call.ext[q], &call.ext[r]);
    let (lo, hi) = if i.is_reverse { (i.entry, j.exit) } else { (j.exit, i.entry) };
    (call.read_idx, j.read_end as u32, i.read_start as u32, lo as u32, hi as u32, i.chrom_id, i.is_reverse)
}

// ── Main study ───────────────────────────────────────────────────────────────

const NPROXY: usize = 4;
const PROXY_NAMES: [&str; NPROXY] = ["4-mer spectrum", "4-mer banded", "8-mer spectrum", "1-indel ungapped"];

pub fn run(calls: &[Call], gaplog: &str) {
    let cfg = CFG;
    let recs = load_gaplog(gaplog);
    let n = recs.len();
    let reads: HashSet<u32> = recs.iter().map(|g| g.read_idx).collect();
    let total_ns: u64 = recs.iter().map(|g| g.ns).sum();
    println!("=== 1. Gap alignments made by the DP ===");
    println!(
        "{} alignments over {} reads ({:.0}/read); total {:.2} s, mean {:.1} us",
        n,
        reads.len(),
        n as f64 / reads.len() as f64,
        total_ns as f64 / 1e9,
        total_ns as f64 / 1e3 / n as f64
    );

    // Size profile, by the longer side.
    println!("\nby max(read gap, ref gap):");
    println!("{:>10} {:>7} {:>7} {:>9} {:>11} {:>9} {:>9}", "bp", "%n", "%time", "mean us", "median bits", "%bits>0", "identity");
    let mut by: Vec<Vec<usize>> = vec![Vec::new(); 7];
    for (k, g) in recs.iter().enumerate() {
        by[lbin(g.qg().max(g.rg()))].push(k);
    }
    for b in 0..7 {
        let ix = &by[b];
        if ix.is_empty() {
            continue;
        }
        let ns: u64 = ix.iter().map(|&k| recs[k].ns).sum();
        let mut bits: Vec<f64> = ix.iter().filter(|&&k| recs[k].ok).map(|&k| recs[k].bits).collect();
        let pos = bits.iter().filter(|&&v| v > 0.0).count();
        let (mut am, mut at) = (0u64, 0u64);
        for &k in ix {
            let g = &recs[k];
            am += g.m as u64;
            at += (g.m + g.x + g.ins + g.del) as u64;
        }
        println!(
            "{:>10} {:>6.1}% {:>6.1}% {:>9.1} {:>11.1} {:>8.1}% {:>8.1}%",
            LABELS[b],
            100.0 * ix.len() as f64 / n as f64,
            100.0 * ns as f64 / total_ns as f64,
            ns as f64 / 1e3 / ix.len() as f64,
            median(&mut bits),
            100.0 * pos as f64 / ix.len() as f64,
            100.0 * am as f64 / at.max(1) as f64
        );
    }
    let marg = |f: &dyn Fn(&GapRec) -> usize, title: &str| {
        let mut c = [0usize; 7];
        for g in &recs {
            c[lbin(f(g))] += 1;
        }
        let parts: Vec<String> = (0..7).map(|b| format!("{} {:.1}%", LABELS[b], 100.0 * c[b] as f64 / n as f64)).collect();
        println!("{title}: {}", parts.join(" | "));
    };
    println!();
    marg(&|g| g.qg(), "read gap");
    marg(&|g| g.rg(), "ref gap ");
    let diff_edges = [0usize, 1, 10, 50, 200, 1000];
    let mut dc = [0usize; 6];
    for g in &recs {
        dc[diff_edges.iter().rposition(|&e| g.qg().abs_diff(g.rg()) >= e).unwrap()] += 1;
    }
    println!(
        "|read gap - ref gap|: 0 {:.1}% | 1-9 {:.1}% | 10-49 {:.1}% | 50-199 {:.1}% | 200-999 {:.1}% | 1000+ {:.1}%",
        100.0 * dc[0] as f64 / n as f64,
        100.0 * dc[1] as f64 / n as f64,
        100.0 * dc[2] as f64 / n as f64,
        100.0 * dc[3] as f64 / n as f64,
        100.0 * dc[4] as f64 / n as f64,
        100.0 * dc[5] as f64 / n as f64
    );
    let ref_empty = recs.iter().filter(|g| g.rg() == 0).count();
    println!("ref gap empty (pure insertion): {:.1}%", 100.0 * ref_empty as f64 / n as f64);
    let mut ident: Vec<f64> = recs
        .iter()
        .filter(|g| g.ok && g.m + g.x + g.ins + g.del > 0)
        .map(|g| g.m as f64 / (g.m + g.x + g.ins + g.del) as f64)
        .collect();
    let multi_indel = recs.iter().filter(|g| g.ok && g.ops >= 2).count();
    println!(
        "median per-alignment identity {:.2}%;  alignments with >= 2 indels: {:.1}%",
        100.0 * median(&mut ident),
        100.0 * multi_indel as f64 / n as f64
    );
    // The cache is keyed by seed pair; how many alignments repeat identical
    // gap coordinates within a read?
    let distinct: HashSet<Key> = recs.iter().map(|g| g.key()).collect();
    let mut seen: HashSet<Key> = HashSet::new();
    let dup_ns: u64 = recs.iter().filter(|g| !seen.insert(g.key())).map(|g| g.ns).sum();
    println!(
        "distinct gap coordinates: {} of {} alignments ({:.1}% are repeats of an identical gap in the same read, {:.1}% of alignment time)",
        distinct.len(),
        n,
        100.0 * (n - distinct.len()) as f64 / n as f64,
        100.0 * dup_ns as f64 / total_ns as f64
    );

    // ── 2. Proxies ──
    println!("\n=== 2. Cheap proxies vs the real alignment's bits ===");
    let mut s = Scratch::default();
    let proxy = |p: usize, q: &[u8], r: &[u8], s: &mut Scratch| -> f64 {
        match p {
            0 => spectrum(q, r, 4, s),
            1 => banded(q, r, 4, 4, s),
            2 => spectrum(q, r, 8, s),
            _ => one_indel_bits(q, r, &cfg),
        }
    };
    // Raw proxy values for the real pairs, and timing.
    let mut raw = vec![vec![0f64; n]; NPROXY];
    let mut proxy_us = [0f64; NPROXY];
    for p in 0..NPROXY {
        let t0 = Instant::now();
        for (k, g) in recs.iter().enumerate() {
            raw[p][k] = proxy(p, &g.q, &g.r, &mut s);
        }
        proxy_us[p] = t0.elapsed().as_secs_f64() * 1e6 / n as f64;
    }
    // Decoys: each read gap against the ref gap of an alignment from another
    // read with a similar ref-gap length, to see what unrelated sequence scores.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&k| recs[k].rg());
    let mut decoy_of = vec![usize::MAX; n];
    for t in 0..n {
        let k = order[t];
        for dt in 1..50 {
            let c = order[(t + dt) % n];
            if recs[c].read_idx != recs[k].read_idx {
                decoy_of[k] = c;
                break;
            }
        }
    }
    let mut draw = vec![vec![f64::NAN; n]; NPROXY];
    for p in 0..NPROXY {
        for k in 0..n {
            if decoy_of[k] != usize::MAX {
                draw[p][k] = proxy(p, &recs[k].q, &recs[decoy_of[k]].r, &mut s);
            }
        }
    }
    // Background match fraction per length bin (by min side), from decoys.
    let kof = [4usize, 4, 8, 0];
    let mut b0 = vec![[0f64; 7]; NPROXY];
    for p in 0..3 {
        let mut sum = [0f64; 7];
        let mut cnt = [0usize; 7];
        for k in 0..n {
            if !draw[p][k].is_nan() {
                let b = lbin(recs[k].qg().min(recs[decoy_of[k]].rg()));
                sum[b] += draw[p][k];
                cnt[b] += 1;
            }
        }
        for b in 0..7 {
            b0[p][b] = if cnt[b] > 0 { sum[b] / cnt[b] as f64 } else { 0.0 };
        }
    }
    let to_est = |p: usize, f: f64, qg: usize, rg: usize| -> f64 {
        if p == 3 { f } else { est_from_fraction(f, b0[p][lbin(qg.min(rg))], kof[p], qg, rg, &cfg) }
    };
    let est: Vec<Vec<f64>> = (0..NPROXY)
        .map(|p| (0..n).map(|k| to_est(p, raw[p][k], recs[k].qg(), recs[k].rg())).collect())
        .collect();
    let dest: Vec<Vec<f64>> = (0..NPROXY)
        .map(|p| {
            (0..n)
                .map(|k| if decoy_of[k] == usize::MAX { f64::NAN } else { to_est(p, draw[p][k], recs[k].qg(), recs[decoy_of[k]].rg()) })
                .collect()
        })
        .collect();

    println!("\nbackground (unrelated-pair) match fraction by min side:");
    println!("{:>18} {}", "", LABELS.iter().map(|l| format!("{:>9}", l)).collect::<String>());
    for p in 0..3 {
        println!("{:>18} {}", PROXY_NAMES[p], (0..7).map(|b| format!("{:>9.3}", b0[p][b])).collect::<String>());
    }

    println!("\naccuracy (est = proxy turned into bits; real pairs with an alignment):");
    println!(
        "{:>18} {:>8} {:>8} {:>10} {:>10} {:>12} {:>13}",
        "proxy", "us/gap", "pearson", "med|err|", "AUC bits>0", "decoys est>0", "real>0 est>0"
    );
    let okix: Vec<usize> = (0..n).filter(|&k| recs[k].ok).collect();
    let truth: Vec<f64> = okix.iter().map(|&k| recs[k].bits).collect();
    for p in 0..NPROXY {
        let e: Vec<f64> = okix.iter().map(|&k| est[p][k]).collect();
        let mut errs: Vec<f64> = e.iter().zip(&truth).map(|(a, b)| (a - b).abs()).collect();
        let mut pairs: Vec<(f64, bool)> = e.iter().zip(&truth).map(|(&a, &b)| (a, b > 0.0)).collect();
        let dec_pos = (0..n).filter(|&k| !dest[p][k].is_nan() && dest[p][k] > 0.0).count();
        let dec_n = (0..n).filter(|&k| !dest[p][k].is_nan()).count();
        let real_pos: Vec<usize> = okix.iter().copied().filter(|&k| recs[k].bits > 0.0).collect();
        let sens = real_pos.iter().filter(|&&k| est[p][k] > 0.0).count();
        println!(
            "{:>18} {:>8.2} {:>8.3} {:>10.1} {:>10.3} {:>11.1}% {:>12.1}%",
            PROXY_NAMES[p],
            proxy_us[p],
            pearson(&e, &truth),
            median(&mut errs),
            auc(&mut pairs),
            100.0 * dec_pos as f64 / dec_n.max(1) as f64,
            100.0 * sens as f64 / real_pos.len().max(1) as f64
        );
    }
    println!("(real alignment: {:.1} us/gap on average)", total_ns as f64 / 1e3 / n as f64);

    println!("\nAUC (bits>0 vs <=0) by max side:");
    println!("{:>18} {}", "", LABELS.iter().map(|l| format!("{:>9}", l)).collect::<String>());
    for p in 0..NPROXY {
        let row: String = (0..7)
            .map(|b| {
                let mut pairs: Vec<(f64, bool)> =
                    by[b].iter().filter(|&&k| recs[k].ok).map(|&k| (est[p][k], recs[k].bits > 0.0)).collect();
                format!("{:>9.3}", auc(&mut pairs))
            })
            .collect();
        println!("{:>18} {}", PROXY_NAMES[p], row);
    }

    // ── 3. DP impact ──
    println!("\n=== 3. Effect on the DP (replaying {} best_chain calls) ===", calls.len());
    let index: HashMap<Key, usize> = recs.iter().enumerate().map(|(k, g)| (g.key(), k)).collect();
    let mut missing = 0usize;
    let max_bridge = 120i64;
    let mut run_strategy = |name: &str,
                            decide: &mut dyn FnMut(&Call, usize, usize, usize, f64) -> Decision,
                            base: Option<&Vec<Option<(f64, Vec<usize>)>>>|
     -> Vec<Option<(f64, Vec<usize>)>> {
        let mut aligned: HashSet<usize> = HashSet::new();
        let mut lookups = 0usize;
        let mut out = Vec::with_capacity(calls.len());
        for call in calls {
            let mut gapf = |q: usize, r: usize, need: f64| -> Option<f64> {
                lookups += 1;
                let Some(&k) = index.get(&key_for(call, q, r)) else {
                    missing += 1;
                    return None;
                };
                match decide(call, q, r, k, need) {
                    Decision::Skip => None,
                    Decision::Value(v) => Some(v),
                    Decision::Align => {
                        aligned.insert(k);
                        if recs[k].ok { Some(recs[k].bits) } else { None }
                    }
                }
            };
            out.push(replay(call, &cfg, &mut gapf));
        }
        // Alignment time a coordinate-keyed cache would spend: each distinct
        // gap aligned at most once per read.
        let used_ns: u64 = aligned.iter().map(|&k| recs[k].ns).sum();
        // With today's seed-pair-keyed cache, every logged alignment of a gap
        // that is still looked up is repeated.
        let aligned_keys: HashSet<Key> = aligned.iter().map(|&k| recs[k].key()).collect();
        let used_ns_seedpair: u64 = recs.iter().filter(|g| aligned_keys.contains(&g.key())).map(|g| g.ns).sum();
        let n_seedpair = recs.iter().filter(|g| aligned_keys.contains(&g.key())).count();
        eprintln!(
            "SEEDPAIR {name}: alignments {} of {} ({:.1}% fewer), alignment time {:.1}% less",
            n_seedpair,
            n,
            100.0 * (1.0 - n_seedpair as f64 / n as f64),
            100.0 * (1.0 - used_ns_seedpair as f64 / total_ns as f64)
        );
        let (mut chains, mut recs_chg, mut maxd) = (0usize, 0usize, 0f64);
        if let Some(b) = base {
            for ((x, y), call) in b.iter().zip(&out).zip(calls) {
                match (x, y) {
                    (Some(x), Some(y)) => {
                        if x.1 != y.1 {
                            chains += 1;
                            if !same_records(&records(call, &x.1, max_bridge), &records(call, &y.1, max_bridge), 50) {
                                recs_chg += 1;
                            }
                        }
                        maxd = maxd.max((x.0 - y.0).abs());
                    }
                    (None, None) => {}
                    _ => {
                        chains += 1;
                        recs_chg += 1;
                    }
                }
            }
        }
        println!(
            "{:<44} {:>8.1}% {:>9.2}% {:>9.2}% {:>10.1}   ({} lookups)",
            name,
            100.0 * (1.0 - used_ns as f64 / total_ns as f64),
            100.0 * chains as f64 / calls.len() as f64,
            100.0 * recs_chg as f64 / calls.len() as f64,
            maxd,
            lookups
        );
        out
    };
    println!("time saved = vs today's alignment time; 'records chg' = emitted records differ by > 50 bp or in locus");
    println!("{:<44} {:>9} {:>10} {:>10} {:>11}", "strategy", "time saved", "chains chg", "records chg", "max |dscore|");
    let exact = run_strategy("exact, cache keyed by gap coordinates", &mut |_, _, _, _, _| Decision::Align, None);
    let ub = |k: usize| upper_bound(recs[k].qg(), recs[k].rg(), &cfg);
    run_strategy("  + upper-bound screen (exact)", &mut |_, _, _, k, need| {
        if ub(k) <= need { Decision::Skip } else { Decision::Align }
    }, Some(&exact));
    run_strategy("  + skip gaps spanning a seed", &mut |c, q, r, k, need| {
        if ub(k) <= need || spans_seed(c, q, r) { Decision::Skip } else { Decision::Align }
    }, Some(&exact));
    run_strategy("  + align only read gaps <= max_bridge", &mut |_, _, _, k, need| {
        if ub(k) <= need || recs[k].qg() as i64 > max_bridge { Decision::Skip } else { Decision::Align }
    }, Some(&exact));
    {
        let e8 = &est[2];
        run_strategy("  + 8-mer estimate for read gaps > max_bridge", &mut |_, _, _, k, need| {
            if ub(k) <= need {
                Decision::Skip
            } else if recs[k].qg() as i64 > max_bridge {
                Decision::Value(e8[k])
            } else {
                Decision::Align
            }
        }, Some(&exact));
    }
    for p in 0..NPROXY {
        let ep = &est[p];
        run_strategy(&format!("{}: replace alignment", PROXY_NAMES[p]), &mut |_, _, _, k, _| Decision::Value(ep[k]), Some(&exact));
        for margin in [10.0, 30.0] {
            run_strategy(&format!("{}: bound + screen, margin {margin}", PROXY_NAMES[p]), &mut |_, _, _, k, need| {
                if ub(k) <= need || ep[k] + margin <= need { Decision::Skip } else { Decision::Align }
            }, Some(&exact));
        }
    }
    // A build with gap_bits_bound never aligns the gaps it rules out, so they
    // are missing from its gap log; they score as no alignment here, which
    // cannot change a decision since such a gap could not win anyway.
    println!("(lookups without a logged alignment, i.e. ruled out by gap_bits_bound in the run: {missing})");
}

enum Decision {
    Skip,
    Align,
    Value(f64),
}
