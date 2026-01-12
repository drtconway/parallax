//! Wavefront Alignment Algorithm (WFA) for sequence alignment.
//!
//! Implements affine-gap WFA as described by Santiago Marco-Sola et al.
//! Time complexity: O(ns) where n is sequence length and s is the alignment score.
//! This is very fast for similar sequences where s << n.

use std::cmp::{max, min};

/// Alignment scoring parameters
#[derive(Clone, Copy, Debug)]
pub struct AlignParams {
    /// Mismatch penalty (positive value)
    pub mismatch: i32,
    /// Gap open penalty (positive value)
    pub gap_open: i32,
    /// Gap extend penalty (positive value)
    pub gap_extend: i32,
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            mismatch: 4,
            gap_open: 6,
            gap_extend: 2,
        }
    }
}

/// CIGAR operation
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),    // '='
    Mismatch(u32), // 'X'
    Ins(u32),      // 'I' - insertion in query
    Del(u32),      // 'D' - deletion in query
    SoftClip(u32), // 'S' - soft clipped bases
}

impl CigarOp {
    /// Format as SAM-style CIGAR string
    pub fn to_string(&self) -> String {
        match self {
            CigarOp::Match(n) => format!("{}=", n),
            CigarOp::Mismatch(n) => format!("{}X", n),
            CigarOp::Ins(n) => format!("{}I", n),
            CigarOp::Del(n) => format!("{}D", n),
            CigarOp::SoftClip(n) => format!("{}S", n),
        }
    }
}

/// Alignment result
#[derive(Clone, Debug)]
pub struct Alignment {
    /// Alignment score (edit distance style - lower is better)
    pub score: i32,
    /// CIGAR operations
    pub cigar: Vec<CigarOp>,
}

impl Alignment {
    /// Format CIGAR as string
    pub fn cigar_string(&self) -> String {
        self.cigar.iter().map(|op| op.to_string()).collect()
    }

    /// Compute the query (read) length consumed by this CIGAR.
    /// This is the sum of M, I, S, =, X operations.
    /// For valid SAM, this must equal the length of the SEQ field.
    pub fn query_length(&self) -> u64 {
        self.cigar.iter().map(|op| match op {
            CigarOp::Match(n) => *n as u64,
            CigarOp::Mismatch(n) => *n as u64,
            CigarOp::Ins(n) => *n as u64,
            CigarOp::SoftClip(n) => *n as u64,
            CigarOp::Del(_) => 0,
        }).sum()
    }

    /// Compute the reference span consumed by this CIGAR.
    /// This is the sum of M, D, N, =, X operations.
    pub fn reference_span(&self) -> u64 {
        self.cigar.iter().map(|op| match op {
            CigarOp::Match(n) => *n as u64,
            CigarOp::Mismatch(n) => *n as u64,
            CigarOp::Del(n) => *n as u64,
            CigarOp::Ins(_) => 0,
            CigarOp::SoftClip(_) => 0,
        }).sum()
    }

    /// Merge adjacent operations of same type
    pub fn normalize(&mut self) {
        if self.cigar.is_empty() {
            return;
        }
        let mut merged = Vec::with_capacity(self.cigar.len());
        for op in self.cigar.drain(..) {
            if let Some(last) = merged.last_mut() {
                match (last, op) {
                    (CigarOp::Match(n), CigarOp::Match(m)) => *n += m,
                    (CigarOp::Mismatch(n), CigarOp::Mismatch(m)) => *n += m,
                    (CigarOp::Ins(n), CigarOp::Ins(m)) => *n += m,
                    (CigarOp::Del(n), CigarOp::Del(m)) => *n += m,
                    (CigarOp::SoftClip(n), CigarOp::SoftClip(m)) => *n += m,
                    _ => merged.push(op),
                }
            } else {
                merged.push(op);
            }
        }
        self.cigar = merged;
    }
}

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
    Match,
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
}

impl WfAligner {
    pub fn new(params: AlignParams) -> Self {
        Self {
            params,
            max_score: 10000,
        }
    }

    pub fn with_max_score(mut self, max_score: i32) -> Self {
        self.max_score = max_score;
        self
    }

    /// Align query to reference, returning the alignment.
    /// Uses affine-gap WFA with traceback.
    pub fn align(&self, query: &[u8], reference: &[u8]) -> Option<Alignment> {
        let n = query.len() as i32;
        let m = reference.len() as i32;

        if n == 0 && m == 0 {
            return Some(Alignment {
                score: 0,
                cigar: vec![],
            });
        }
        if n == 0 {
            return Some(Alignment {
                score: self.params.gap_open + self.params.gap_extend * m,
                cigar: vec![CigarOp::Del(m as u32)],
            });
        }
        if m == 0 {
            return Some(Alignment {
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
            return Some(self.traceback(query, reference, 0, &wf_m, &wf_i, &wf_d, &trace));
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

            wf_m.push(new_m);
            wf_i.push(new_i);
            wf_d.push(new_d);
            trace.push(trace_s);

            // Check if we reached the target
            if wf_m[s_idx].get(target_k) >= n {
                return Some(self.traceback(query, reference, s, &wf_m, &wf_i, &wf_d, &trace));
            }
        }

        None // Exceeded max score
    }

    /// Greedy match extension on all diagonals
    fn extend(&self, query: &[u8], reference: &[u8], wf: &mut Wavefront) {
        let n = query.len() as i32;
        let m = reference.len() as i32;

        for k in wf.lo..=wf.hi {
            let mut row = wf.get(k);
            let mut col = row - k;

            while row < n && col >= 0 && col < m {
                // Case-insensitive comparison (FASTA may have lowercase repeat-masked regions)
                if query[row as usize].to_ascii_uppercase() == reference[col as usize].to_ascii_uppercase() {
                    row += 1;
                    col += 1;
                } else {
                    break;
                }
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
                    let row_before_extend = self.row_before_extend_full(s, k, wf_m, wf_i, wf_d, x, o, e);
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
                            TraceOp::Match => {
                                // Match operations are handled through extension, not trace
                                // This shouldn't happen, but if it does, treat as mismatch fallback
                                if s >= x && row > 0 && col > 0 {
                                    cigar.push(CigarOp::Mismatch(1));
                                    row -= 1;
                                    col -= 1;
                                    s -= x;
                                } else {
                                    break;
                                }
                            }
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
                    let i_from_i = if s >= e { wf_i[(s - e) as usize].get(k - 1) } else { i32::MIN / 2 };
                    let i_from_m = if s >= o + e { wf_m[(s - o - e) as usize].get(k - 1) + 1 } else { i32::MIN / 2 };

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
                    let d_from_d = if s >= e { wf_d[(s - e) as usize].get(k + 1) } else { i32::MIN / 2 };
                    let d_from_m = if s >= o + e { wf_m[(s - o - e) as usize].get(k + 1) } else { i32::MIN / 2 };

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

/// Convenience function for quick alignment with default parameters
pub fn align(query: &[u8], reference: &[u8]) -> Option<Alignment> {
    WfAligner::new(AlignParams::default()).align(query, reference)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical() {
        let result = align(b"ACGT", b"ACGT").unwrap();
        assert_eq!(result.score, 0);
        assert_eq!(result.cigar_string(), "4=");
    }

    #[test]
    fn test_single_mismatch() {
        let result = align(b"ACGT", b"ACTT").unwrap();
        assert_eq!(result.score, 4); // mismatch penalty
        assert_eq!(result.cigar_string(), "2=1X1=");
    }

    #[test]
    fn test_single_insertion() {
        let result = align(b"ACGT", b"ACT").unwrap();
        // query has extra G
        assert!(result.score > 0);
        assert!(result.cigar_string().contains('I'));
    }

    #[test]
    fn test_single_deletion() {
        let result = align(b"ACT", b"ACGT").unwrap();
        // query missing G
        assert!(result.score > 0);
        assert!(result.cigar_string().contains('D'));
    }

    #[test]
    fn test_empty() {
        let result = align(b"", b"").unwrap();
        assert_eq!(result.score, 0);
        assert!(result.cigar.is_empty());
    }

    #[test]
    fn test_query_empty() {
        let result = align(b"", b"ACGT").unwrap();
        assert_eq!(result.cigar_string(), "4D");
    }

    #[test]
    fn test_reference_empty() {
        let result = align(b"ACGT", b"").unwrap();
        assert_eq!(result.cigar_string(), "4I");
    }

    #[test]
    fn test_longer_sequences() {
        let query = b"ACGTACGTACGT";
        let reference = b"ACGTACGTACGT";
        let result = align(query, reference).unwrap();
        assert_eq!(result.score, 0);
        assert_eq!(result.cigar_string(), "12=");
    }

    #[test]
    fn test_with_gaps() {
        let query = b"ACGTACGT";
        let reference = b"ACGTTTTACGT";
        let result = align(query, reference).unwrap();
        assert!(result.score > 0);
        // Should have a deletion in the middle
        println!("CIGAR: {}", result.cigar_string());
    }
}
