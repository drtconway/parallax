//! Variant C: the blocked kernel with explicit fearless_simd vectors.

use fearless_simd::prelude::*;
use fearless_simd::{Level, dispatch, f32x8, i32x8, u32x8};
use fearless_simd_macros::simd;

use crate::{CfgF32, Ext, GAP_ALIGN_PREDECESSORS, L2A, PAD, Soa, W, log2_1p};

/// log2(1 + x) per lane; bit-identical to the scalar `log2_1p`.
#[inline(always)]
fn log2_1p_v<S: Simd>(simd: S, x: f32x8<S>) -> f32x8<S> {
    let y = x + 1.0;
    let bits: u32x8<S> = y.to_bits();
    let e_i: i32x8<S> = (bits >> 23).bitcast::<i32x8<S>>() - 127;
    let e = f32x8::float_from(e_i);
    let t = f32x8::from_bits((bits & 0x007f_ffff) | 0x3f80_0000) - 1.0;
    let one = f32x8::splat(simd, 1.0);
    let poly = L2A_V(simd, t);
    let p = t + t * (one - t) * poly;
    e + p
}

#[allow(non_snake_case)]
#[inline(always)]
fn L2A_V<S: Simd>(simd: S, t: f32x8<S>) -> f32x8<S> {
    // L2A[0] + t * (L2A[1] + t * (L2A[2] + t * L2A[3])), same operation order as scalar.
    let a3 = f32x8::splat(simd, L2A[3]);
    let a2 = f32x8::splat(simd, L2A[2]);
    let a1 = f32x8::splat(simd, L2A[1]);
    let a0 = f32x8::splat(simd, L2A[0]);
    a0 + t * (a1 + t * (a2 + t * a3))
}

pub fn self_test() {
    let level = Level::new();
    dispatch!(level, simd => self_test_k(simd));
}

#[simd]
fn self_test_k<S: Simd>(simd: S) {
    let a = i32x8::from_slice(simd, &[0, 1, 2, 3, 4, 5, 6, 7]);
    let m = a.simd_lt(i32x8::splat(simd, 3)).to_bitmask();
    assert_eq!(m, 0b111, "bitmask lane order");
    let mut worst = 0f64;
    for base in [0u32, 1, 7, 30, 100, 1000, 4095, 123_456, 1 << 20, 1 << 24, 1 << 28] {
        let xs: [f32; 8] = core::array::from_fn(|l| (base + l as u32 * 3) as f32);
        let v = log2_1p_v(simd, f32x8::from_slice(simd, &xs));
        for l in 0..W {
            let got = v.as_slice()[l];
            assert_eq!(got.to_bits(), log2_1p(xs[l]).to_bits(), "simd/scalar log2 mismatch at {}", xs[l]);
            let exact = (1.0 + xs[l] as f64).log2();
            worst = worst.max((got as f64 - exact).abs());
        }
    }
    assert!(worst < 1e-4, "log2 approximation error {worst}");
    eprintln!("self-test ok: bitmask lane order, log2_1p simd==scalar, max err vs exact {worst:.2e}");
}

pub fn best_chain_fearless(
    level: Level,
    skip: bool,
    ext: &[Ext],
    order: &[usize],
    mask: &[bool],
    cfg: &CfgF32,
    soa: &mut Soa,
    gap: &mut dyn FnMut(usize, usize) -> Option<f64>,
) -> Option<(f32, Vec<usize>)> {
    soa.fill(ext, order, mask);
    dispatch!(level, simd => kernel(simd, skip, cfg, soa, gap));
    soa.traceback()
}

#[simd]
fn kernel<S: Simd>(
    simd: S,
    skip: bool,
    cfg: &CfgF32,
    soa: &mut Soa,
    gap: &mut dyn FnMut(usize, usize) -> Option<f64>,
) {
    let n = soa.n;
    let zero_i = i32x8::splat(simd, 0);
    let zero_f = f32x8::splat(simd, 0.0);
    let neg_inf = f32x8::splat(simd, f32::NEG_INFINITY);
    let v_mb = f32x8::splat(simd, cfg.mb);
    let v_mismatch = f32x8::splat(simd, cfg.mismatch);
    let v_open = f32x8::splat(simd, cfg.open);
    let v_extend = f32x8::splat(simd, cfg.extend);
    let v_sv = f32x8::splat(simd, cfg.sv);
    let v_sv1 = f32x8::splat(simd, cfg.sv1);
    let v_svg = f32x8::splat(simd, cfg.svg);
    let v_gmin = i32x8::splat(simd, cfg.gmin);
    let v_gmax = i32x8::splat(simd, cfg.gmax);

    for r in 0..n {
        let ri = PAD + r;
        let (irs, ire, ient, idiag, ikey) = (soa.rs[ri], soa.re[ri], soa.entry[ri], soa.diag[ri], soa.key[ri]);
        let irev = ikey & 1 == 1;
        let coverage = cfg.mb * (ire - irs) as f32;
        let v_irs = i32x8::splat(simd, irs);
        let v_ire = i32x8::splat(simd, ire);
        let v_ient = i32x8::splat(simd, ient);
        let v_idiag = i32x8::splat(simd, idiag);
        let v_ikey = i32x8::splat(simd, ikey);
        let v_ichrom = i32x8::splat(simd, ikey >> 1);
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
            let rej = i32x8::from_slice(simd, &soa.re[lo..hi]);
            let exj = i32x8::from_slice(simd, &soa.exit[lo..hi]);
            let kj = i32x8::from_slice(simd, &soa.key[lo..hi]);
            let valid = rej.simd_lt(v_ire);
            let same_strand = kj.simd_eq(v_ikey);
            let d = v_ient - exj;
            let read_gap = v_irs - rej;
            if !(skip && hopeless) {
                let dgj = i32x8::from_slice(simd, &soa.diag[lo..hi]);
                let dpj = f32x8::from_slice(simd, &soa.dp[lo..hi]);
                let same_chrom = (kj >> 1).simd_eq(v_ichrom);
                let dist = d.abs();
                let shift = (v_idiag - dgj).abs();
                let x = same_strand.select(shift, dist);
                let lg = log2_1p_v(simd, f32x8::float_from(x));
                let indel_nz = v_open + v_extend * f32x8::float_from(shift);
                let indel_z = read_gap.simd_gt(zero_i).select(v_mismatch, zero_f);
                let indel = shift.simd_eq(zero_i).select(indel_z, indel_nz);
                let c_same = indel.min(v_sv + lg);
                let c_strand = v_sv1 + lg;
                let cost = same_strand.select(c_same, same_chrom.select(c_strand, v_svg));
                let new_bases = v_ire - v_irs.max(rej);
                let c = dpj + v_mb * f32x8::float_from(new_bases) - cost;
                let cand = valid.select(c, neg_inf);
                let bm = cand.reduce_max();
                if bm > best {
                    // Highest lane among the maxima, matching the descending scan.
                    let m = cand.simd_eq(f32x8::splat(simd, bm)).to_bitmask();
                    let l = 63 - m.leading_zeros() as usize;
                    best = bm;
                    best_q = (lo + l - PAD) as u32;
                }
            }
            if nn < GAP_ALIGN_PREDECESSORS {
                let ref_gap = if irev { -d } else { d };
                let elig = valid
                    & same_strand
                    & read_gap.simd_ge(v_gmin)
                    & read_gap.simd_le(v_gmax)
                    & ref_gap.simd_ge(zero_i)
                    & ref_gap.simd_le(v_gmax);
                let mut bits = elig.to_bitmask();
                while bits != 0 && nn < GAP_ALIGN_PREDECESSORS {
                    let l = 63 - bits.leading_zeros() as usize;
                    nearest[nn] = lo + l - PAD;
                    nn += 1;
                    bits &= !(1u64 << l);
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
}
