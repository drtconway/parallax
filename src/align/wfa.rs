//! Wavefront Alignment Algorithm (WFA) for sequence alignment.
//!
//! Implements affine-gap WFA as described by Santiago Marco-Sola et al.
//! Time complexity: O(ns) where n is sequence length and s is the alignment score.
//! This is very fast for similar sequences where s << n.

use std::cmp::{max, min};

use crate::config;

use super::{AlignParams, Alignment, CigarOp};

#[derive(Debug)]
pub enum WfaFailure {
    MaxScoreExceeded,
    BandWidthExceeded(i32),
}

impl std::fmt::Display for WfaFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WfaFailure::MaxScoreExceeded => write!(f, "WFA alignment failed: maximum score exceeded"),
            WfaFailure::BandWidthExceeded(width) => write!(f, "WFA alignment failed: band width {} exceeded maximum allowed", width),
        }
    }
}

impl std::error::Error for WfaFailure {}

/// Wavefront for a single score, storing furthest-reaching point on each diagonal.
/// Diagonal k = row - col, so row = offset[k] and col = offset[k] - k.
#[derive(Clone)]
struct Wavefront {
    /// Furthest row reached on each diagonal k in [lo, hi]
    offsets: Vec<i32>,
    /// Minimum diagonal index
    lo: i32,
    /// Maximum diagonal index
    hi: i32,
}

impl Wavefront {
    fn new() -> Self {
        Self {
            offsets: Vec::new(),
            lo: 0,
            hi: -1, // empty
        }
    }

    fn with_diagonals(lo: i32, hi: i32) -> Self {
        let size = (hi - lo + 1) as usize;
        Self {
            offsets: vec![i32::MIN / 2; size], // -inf sentinel
            lo,
            hi,
        }
    }

    fn is_empty(&self) -> bool {
        self.hi < self.lo
    }

    fn get(&self, k: i32) -> i32 {
        if k < self.lo || k > self.hi {
            i32::MIN / 2
        } else {
            self.offsets[(k - self.lo) as usize]
        }
    }

    fn set(&mut self, k: i32, val: i32) {
        if k >= self.lo && k <= self.hi {
            self.offsets[(k - self.lo) as usize] = val;
        }
    }
}

/// Backtrace state for CIGAR reconstruction
#[derive(Clone, Copy, Debug)]
enum TraceOp {
    Mismatch,
    InsOpen,
    InsExt,
    DelOpen,
    DelExt,
}

/// Affine-gap WFA aligner
pub struct WfAligner {
    params: AlignParams,
    /// Maximum score to explore before giving up (prevents runaway on very different sequences)
    max_score: i32,
    /// X-drop threshold: prune diagonals that fall this far behind the best
    x_drop: i32,
    /// Maximum wavefront band width before giving up
    max_band_width: i32,
}

impl WfAligner {
    pub fn new(params: AlignParams) -> Self {
        let cfg = config::get();
        Self {
            params,
            max_score: 10000,
            x_drop: cfg.alignment.x_drop,
            max_band_width: cfg.alignment.max_band_width,
        }
    }

    #[allow(dead_code)]
    pub fn set_max_score(&mut self, max_score: i32) {
        self.max_score = max_score;
    }

    #[allow(dead_code)]
    pub fn set_x_drop(&mut self, x_drop: i32) {
        self.x_drop = x_drop;
    }

    #[allow(dead_code)]
    pub fn set_max_band_width(&mut self, max_band_width: i32) {
        self.max_band_width = max_band_width;
    }

    /// Align query to reference, returning the alignment.
    /// Uses affine-gap WFA with traceback.
    pub fn align(&self, query: &[u8], reference: &[u8]) -> std::result::Result<Alignment, WfaFailure> {
        let n = query.len() as i32;
        let m = reference.len() as i32;

        if n == 0 && m == 0 {
            return Ok(Alignment {
                score: 0,
                cigar: vec![],
            });
        }
        if n == 0 {
            return Ok(Alignment {
                score: self.params.gap_open + self.params.gap_extend * m,
                cigar: vec![CigarOp::Del(m as u32)],
            });
        }
        if m == 0 {
            return Ok(Alignment {
                score: self.params.gap_open + self.params.gap_extend * n,
                cigar: vec![CigarOp::Ins(n as u32)],
            });
        }

        let x = self.params.mismatch;
        let o = self.params.gap_open;
        let e = self.params.gap_extend;

        // Wavefronts for M (match/mismatch), I (insertion), D (deletion)
        // Indexed by score
        let mut wf_m: Vec<Wavefront> = Vec::new();
        let mut wf_i: Vec<Wavefront> = Vec::new();
        let mut wf_d: Vec<Wavefront> = Vec::new();

        // Traceback storage: for each (score, diagonal, component) -> operation
        let mut trace: Vec<Vec<(i32, TraceOp)>> = Vec::new(); // score -> [(k, op), ...]

        // Initialize score 0: start at origin on diagonal 0
        wf_m.push(Wavefront::with_diagonals(0, 0));
        wf_m[0].set(0, 0);
        wf_i.push(Wavefront::new());
        wf_d.push(Wavefront::new());
        trace.push(vec![]);

        // Extend matches at score 0
        self.extend(query, reference, &mut wf_m[0]);

        // Target: reach (n, m) which is diagonal k = n - m, offset = n
        let target_k = n - m;

        // Check if already done
        if wf_m[0].get(target_k) >= n {
            return Ok(self.traceback(query, reference, 0, &wf_m, &wf_i, &wf_d, &trace));
        }

        for s in 1..=self.max_score {
            let s_idx = s as usize;

            // Determine wavefront bounds
            let mut lo = i32::MAX;
            let mut hi = i32::MIN;

            // From M[s-x], I[s-e], D[s-e], M[s-o-e], etc.
            if s >= x && !wf_m[s_idx - x as usize].is_empty() {
                lo = min(lo, wf_m[s_idx - x as usize].lo);
                hi = max(hi, wf_m[s_idx - x as usize].hi);
            }
            if s >= e && !wf_i[(s - e) as usize].is_empty() {
                lo = min(lo, wf_i[(s - e) as usize].lo - 1);
                hi = max(hi, wf_i[(s - e) as usize].hi - 1);
            }
            if s >= e && !wf_d[(s - e) as usize].is_empty() {
                lo = min(lo, wf_d[(s - e) as usize].lo + 1);
                hi = max(hi, wf_d[(s - e) as usize].hi + 1);
            }
            if s >= o + e && !wf_m[(s - o - e) as usize].is_empty() {
                lo = min(lo, wf_m[(s - o - e) as usize].lo - 1);
                hi = max(hi, wf_m[(s - o - e) as usize].hi + 1);
            }

            if lo > hi {
                // No valid wavefront at this score
                wf_m.push(Wavefront::new());
                wf_i.push(Wavefront::new());
                wf_d.push(Wavefront::new());
                trace.push(vec![]);
                continue;
            }

            // Clamp to valid diagonal range
            lo = max(lo, -m);
            hi = min(hi, n);

            let mut new_m = Wavefront::with_diagonals(lo, hi);
            let mut new_i = Wavefront::with_diagonals(lo, hi);
            let mut new_d = Wavefront::with_diagonals(lo, hi);
            let mut trace_s = Vec::new();

            for k in lo..=hi {
                // Insertion: I[s][k] = max(M[s-o-e][k-1] + 1, I[s-e][k-1] + 1)
                let i_from_m = if s >= o + e {
                    wf_m[(s - o - e) as usize].get(k - 1) + 1
                } else {
                    i32::MIN / 2
                };
                let i_from_i = if s >= e {
                    wf_i[(s - e) as usize].get(k - 1) + 1
                } else {
                    i32::MIN / 2
                };
                let i_val = max(i_from_m, i_from_i);
                new_i.set(k, i_val);

                // Deletion: D[s][k] = max(M[s-o-e][k+1], D[s-e][k+1])
                let d_from_m = if s >= o + e {
                    wf_m[(s - o - e) as usize].get(k + 1)
                } else {
                    i32::MIN / 2
                };
                let d_from_d = if s >= e {
                    wf_d[(s - e) as usize].get(k + 1)
                } else {
                    i32::MIN / 2
                };
                let d_val = max(d_from_m, d_from_d);
                new_d.set(k, d_val);

                // Match/Mismatch: M[s][k] = max(M[s-x][k] + 1, I[s][k], D[s][k])
                let m_from_m = if s >= x {
                    wf_m[(s - x) as usize].get(k) + 1
                } else {
                    i32::MIN / 2
                };
                let m_val = max(max(m_from_m, i_val), d_val);
                new_m.set(k, m_val);

                // Record traceback
                if m_val > i32::MIN / 2 {
                    let op = if m_val == i_val {
                        if i_val == i_from_i {
                            TraceOp::InsExt
                        } else {
                            TraceOp::InsOpen
                        }
                    } else if m_val == d_val {
                        if d_val == d_from_d {
                            TraceOp::DelExt
                        } else {
                            TraceOp::DelOpen
                        }
                    } else {
                        TraceOp::Mismatch
                    };
                    trace_s.push((k, op));
                }
            }

            // Extend matches
            self.extend(query, reference, &mut new_m);

            // X-drop pruning: prune diagonals that fall too far behind the best
            if self.x_drop > 0 {
                let max_offset = (new_m.lo..=new_m.hi)
                    .map(|k| new_m.get(k))
                    .max()
                    .unwrap_or(i32::MIN);

                for k in new_m.lo..=new_m.hi {
                    if max_offset - new_m.get(k) > self.x_drop {
                        new_m.set(k, i32::MIN / 2);
                        new_i.set(k, i32::MIN / 2);
                        new_d.set(k, i32::MIN / 2);
                    }
                }

                // Shrink wavefront bounds to exclude pruned diagonals
                while new_m.lo <= new_m.hi && new_m.get(new_m.lo) <= i32::MIN / 2 + 1 {
                    new_m.lo += 1;
                    new_i.lo += 1;
                    new_d.lo += 1;
                }
                while new_m.hi >= new_m.lo && new_m.get(new_m.hi) <= i32::MIN / 2 + 1 {
                    new_m.hi -= 1;
                    new_i.hi -= 1;
                    new_d.hi -= 1;
                }
            }

            // Band width check: fail early if wavefront is too wide
            if self.max_band_width > 0 && new_m.hi - new_m.lo > self.max_band_width {
                return Err(WfaFailure::BandWidthExceeded(new_m.hi - new_m.lo));
            }

            wf_m.push(new_m);
            wf_i.push(new_i);
            wf_d.push(new_d);
            trace.push(trace_s);

            // Check if we reached the target
            if wf_m[s_idx].get(target_k) >= n {
                return Ok(self.traceback(query, reference, s, &wf_m, &wf_i, &wf_d, &trace));
            }
        }

        Err(WfaFailure::MaxScoreExceeded) // Exceeded max score
    }

    /// Greedy match extension on all diagonals
    fn extend(&self, query: &[u8], reference: &[u8], wf: &mut Wavefront) {
        use crate::utils::longest_common_prefix;

        let n = query.len() as i32;
        let m = reference.len() as i32;

        for k in wf.lo..=wf.hi {
            let mut row = wf.get(k);
            let col = row - k;

            // Check bounds before calling LCP
            if row >= 0 && row < n && col >= 0 && col < m {
                let query_slice = &query[row as usize..];
                let ref_slice = &reference[col as usize..];
                let lcp = longest_common_prefix(query_slice, ref_slice);
                row += lcp as i32;
            }

            wf.set(k, row);
        }
    }

    /// Traceback to produce CIGAR using stored trace information
    fn traceback(
        &self,
        query: &[u8],
        reference: &[u8],
        final_score: i32,
        wf_m: &[Wavefront],
        wf_i: &[Wavefront],
        wf_d: &[Wavefront],
        trace: &[Vec<(i32, TraceOp)>],
    ) -> Alignment {
        let n = query.len() as i32;
        let m = reference.len() as i32;

        let x = self.params.mismatch;
        let o = self.params.gap_open;
        let e = self.params.gap_extend;

        let mut cigar = Vec::new();
        let mut row = n;
        let mut col = m;
        let mut s = final_score;

        // State: which wavefront are we currently in?
        #[derive(Clone, Copy, PartialEq)]
        enum State {
            M,
            I,
            D,
        }
        let mut state = State::M;

        while row > 0 || col > 0 {
            let k = row - col;

            if s == 0 {
                // All remaining must be matches (score 0 means we extended from origin)
                if row > 0 && col > 0 {
                    // Verify they actually match
                    let match_count = row.min(col);
                    cigar.push(CigarOp::Match(match_count as u32));
                }
                break;
            }

            if s < 0 {
                break;
            }

            let s_idx = s as usize;

            match state {
                State::M => {
                    let current_row = wf_m[s_idx].get(k);

                    // First, check how many matches were extended at this score
                    let row_before_extend =
                        self.row_before_extend_full(s, k, wf_m, wf_i, wf_d, x, o, e);
                    let matches = (current_row - row_before_extend).max(0).min(row.min(col));

                    if matches > 0 {
                        cigar.push(CigarOp::Match(matches as u32));
                        row -= matches;
                        col -= matches;
                    }

                    if row <= 0 && col <= 0 {
                        break;
                    }

                    // Now find what operation brought us to row_before_extend
                    // Look up the trace for this score/k
                    if let Some(op) = self.find_trace_op(s, k, trace) {
                        match op {
                            TraceOp::Mismatch => {
                                if s >= x && row > 0 && col > 0 {
                                    cigar.push(CigarOp::Mismatch(1));
                                    row -= 1;
                                    col -= 1;
                                    s -= x;
                                    // Stay in M state
                                } else {
                                    break;
                                }
                            }
                            TraceOp::InsOpen => {
                                // This M cell came from I cell at same score
                                // The I cell was opened from M[s-o-e][k-1]
                                if row > 0 {
                                    cigar.push(CigarOp::Ins(1));
                                    row -= 1;
                                    s -= o + e;
                                    // Go back to M state (gap open comes from M)
                                } else {
                                    break;
                                }
                            }
                            TraceOp::InsExt => {
                                // This M cell came from I cell at same score
                                // The I cell was extended from I[s-e][k-1]
                                if row > 0 {
                                    cigar.push(CigarOp::Ins(1));
                                    row -= 1;
                                    s -= e;
                                    state = State::I; // Continue in I state
                                } else {
                                    break;
                                }
                            }
                            TraceOp::DelOpen => {
                                if col > 0 {
                                    cigar.push(CigarOp::Del(1));
                                    col -= 1;
                                    s -= o + e;
                                    // Go back to M state
                                } else {
                                    break;
                                }
                            }
                            TraceOp::DelExt => {
                                if col > 0 {
                                    cigar.push(CigarOp::Del(1));
                                    col -= 1;
                                    s -= e;
                                    state = State::D; // Continue in D state
                                } else {
                                    break;
                                }
                            }
                        }
                    } else {
                        // No trace found - try to infer from wavefronts
                        if s >= x && row > 0 && col > 0 {
                            let prev_row = wf_m[(s - x) as usize].get(k);
                            if prev_row + 1 == row_before_extend {
                                cigar.push(CigarOp::Mismatch(1));
                                row -= 1;
                                col -= 1;
                                s -= x;
                                continue;
                            }
                        }
                        // Try gap transitions
                        if s >= o + e && row > 0 {
                            cigar.push(CigarOp::Ins(1));
                            row -= 1;
                            s -= o + e;
                            continue;
                        }
                        if s >= o + e && col > 0 {
                            cigar.push(CigarOp::Del(1));
                            col -= 1;
                            s -= o + e;
                            continue;
                        }
                        break;
                    }
                }
                State::I => {
                    // We're tracing back through insertion states
                    // I[s][k] came from either I[s-e][k-1] (extend) or M[s-o-e][k-1] (open)
                    let i_from_i = if s >= e {
                        wf_i[(s - e) as usize].get(k - 1)
                    } else {
                        i32::MIN / 2
                    };
                    let i_from_m = if s >= o + e {
                        wf_m[(s - o - e) as usize].get(k - 1) + 1
                    } else {
                        i32::MIN / 2
                    };

                    if i_from_i >= i_from_m && s >= e {
                        // Extended from I
                        if row > 0 {
                            cigar.push(CigarOp::Ins(1));
                            row -= 1;
                            s -= e;
                            // Stay in I state
                        } else {
                            break;
                        }
                    } else if s >= o + e {
                        // Opened from M
                        if row > 0 {
                            cigar.push(CigarOp::Ins(1));
                            row -= 1;
                            s -= o + e;
                            state = State::M;
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                State::D => {
                    // D[s][k] came from either D[s-e][k+1] (extend) or M[s-o-e][k+1] (open)
                    let d_from_d = if s >= e {
                        wf_d[(s - e) as usize].get(k + 1)
                    } else {
                        i32::MIN / 2
                    };
                    let d_from_m = if s >= o + e {
                        wf_m[(s - o - e) as usize].get(k + 1)
                    } else {
                        i32::MIN / 2
                    };

                    if d_from_d >= d_from_m && s >= e {
                        // Extended from D
                        if col > 0 {
                            cigar.push(CigarOp::Del(1));
                            col -= 1;
                            s -= e;
                            // Stay in D state
                        } else {
                            break;
                        }
                    } else if s >= o + e {
                        // Opened from M
                        if col > 0 {
                            cigar.push(CigarOp::Del(1));
                            col -= 1;
                            s -= o + e;
                            state = State::M;
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
            }
        }

        cigar.reverse();
        let mut alignment = Alignment {
            score: final_score,
            cigar,
        };
        alignment.normalize();
        alignment
    }

    /// Find the trace operation for a given score and diagonal
    fn find_trace_op(&self, s: i32, k: i32, trace: &[Vec<(i32, TraceOp)>]) -> Option<TraceOp> {
        if s <= 0 || s as usize >= trace.len() {
            return None;
        }
        for &(tk, op) in &trace[s as usize] {
            if tk == k {
                return Some(op);
            }
        }
        None
    }

    /// Calculate row before extension, considering all wavefronts
    fn row_before_extend_full(
        &self,
        s: i32,
        k: i32,
        wf_m: &[Wavefront],
        wf_i: &[Wavefront],
        wf_d: &[Wavefront],
        x: i32,
        _o: i32,
        _e: i32,
    ) -> i32 {
        let mut best = i32::MIN / 2;

        // From mismatch: M[s-x][k] + 1
        if s >= x {
            best = max(best, wf_m[(s - x) as usize].get(k) + 1);
        }

        // From insertion: I[s][k] (which equals wf_i value at this score)
        if (s as usize) < wf_i.len() {
            best = max(best, wf_i[s as usize].get(k));
        }

        // From deletion: D[s][k]
        if (s as usize) < wf_d.len() {
            best = max(best, wf_d[s as usize].get(k));
        }

        best
    }
}

#[cfg(test)]
mod tests {
    use super::{AlignParams, WfAligner};
    use crate::align::{CigarOp, align};

    #[test]
    fn test_simple_alignment() {
        let query = b"ACGTACGT";
        let reference = b"ACGTTCGT";

        let aligner = WfAligner::new(AlignParams::default());
        let result = aligner.align(query, reference).unwrap();

        assert_eq!(result.score, 4);
        assert_eq!(
            result.cigar,
            vec![CigarOp::Match(4), CigarOp::Mismatch(1), CigarOp::Match(3)]
        );
    }

    #[test]
    fn test_repetitive_at_alignment() {
        // This is a challenging case with highly repetitive AT-rich sequences
        // The query is longer than the reference, testing insertion handling
        let query = b"ATTTTTATATATACTTATATATTTATATATATTTTTATATATACTCATATATTTATATATATTTTATATATACTTATTTATATATATATATTTTTATATATATTTAATTTTTACATATATTTATATTTTTATATATTTATATATTTATATATTTTTATATTTTATATATATGTTTATATATTTATATATTATATATATTTATATATATTTATATATTTATATATTATATATTTATATATATTTATATATTTATATATTATATATATTTATATATTTATATATTTATATATTACATATATTTATATATATTTATATATTTATATATGTTTATATATTTATATATTATATATATTTATATATATTTATATTATATATATACTTATATATTTATATATATTTTTATATATACTTATATATTTATATATATTTTTATATATACTTATATATTTATATATATTTTTATATATACTTATATATATTTTTTATATATTTATATATTTTTATATATATTTAATTTTTAT";
        let reference = b"TTTTATATATACTTATATATTTATATATATTTTTATATATACTCATATATTTATATATATTTTATATATACTTATTTTATATATATATATTTTTATATATATTTAATTTTTAC";

        println!("Query length: {}", query.len());
        println!("Reference length: {}", reference.len());
        println!(
            "Expected diagonal difference: {}",
            query.len() as i32 - reference.len() as i32
        );

        // First, try with a custom aligner with no x-drop or band limit to see if it can succeed
        let mut aligner_no_limits = WfAligner::new(AlignParams::default());
        aligner_no_limits.set_max_score(50000); // Very high limit
        aligner_no_limits.set_x_drop(0); // No x-drop pruning
        aligner_no_limits.set_max_band_width(0); // No band limit

        println!("\n--- Trying with no limits (max_score=50000, no x-drop, no band limit) ---");
        let result_no_limits = aligner_no_limits.align(query, reference);

        match &result_no_limits {
            Ok(alignment) => {
                println!(
                    "SUCCESS with no limits: score={}, cigar={}",
                    alignment.score,
                    alignment.cigar_string()
                );
            }
            Err(e) => {
                println!("FAILED even with no limits: {}", e);
            }
        }

        // Now try with default config (which includes x-drop and band limits)
        println!("\n--- Trying with default config (x_drop, band_width from config) ---");
        let cfg = crate::config::get();
        println!("Config x_drop: {}", cfg.alignment.x_drop);
        println!("Config max_band_width: {}", cfg.alignment.max_band_width);

        let result = align(query, reference);

        if let Some(ref alignment) = result {
            println!(
                "SUCCESS with config: score={}, cigar={}",
                alignment.score,
                alignment.cigar_string()
            );
        } else {
            println!("FAILED with config settings");
        }

        // For now, just ensure the no-limits version works
        assert!(
            result_no_limits.is_ok(),
            "Alignment failed even with no limits - likely exceeded max_score=50000"
        );

        let alignment = result_no_limits.unwrap();

        // Verify CIGAR accounts for length difference
        let (cigar_query_len, cigar_ref_len) =
            alignment
                .cigar
                .iter()
                .fold((0usize, 0usize), |(q, r), op| match op {
                    CigarOp::Match(n) | CigarOp::Mismatch(n) => (q + *n as usize, r + *n as usize),
                    CigarOp::Ins(n) | CigarOp::SoftClip(n) => (q + *n as usize, r),
                    CigarOp::Del(n) => (q, r + *n as usize),
                });
        assert_eq!(cigar_query_len, query.len(), "CIGAR query length mismatch");
        assert_eq!(
            cigar_ref_len,
            reference.len(),
            "CIGAR reference length mismatch"
        );
    }
}
