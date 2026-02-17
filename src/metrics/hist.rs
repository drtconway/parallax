//! Binned histogram with adaptive warmup.
//!
//! Collects the first [`WARMUP_COUNT`] values to learn the data range, then
//! allocates ~[`NUM_WARMUP_BINS`] fixed-width bins spanning that range.
//! Values arriving after warmup that fall outside the initial range are
//! accommodated by growing the bin arrays on demand.
//!
//! This gives O(1) recording and cache-friendly sequential quantile scans,
//! with accuracy proportional to `1 / NUM_WARMUP_BINS`.

/// Number of values to buffer before establishing the bin layout.
const WARMUP_COUNT: usize = 1000;

/// Number of bins to span the warmup range.
const NUM_WARMUP_BINS: usize = 100;

/// A histogram that learns its scale from the first [`WARMUP_COUNT`] values.
///
/// After warmup, `bin_index(v) = floor((v - offset) * scale)` where
/// `offset` = warmup min and `scale` = `NUM_WARMUP_BINS / range`.
/// Bins with index ≥ 0 live in `pos`, bins with index < 0 in `neg`
/// (stored as `neg[(-index) - 1]`).
#[derive(Clone)]
pub struct BinnedHistogram {
    warmup_buffer: Option<Vec<f64>>,
    bins: Bins,
}

#[derive(Clone)]
struct Bins {
    /// Reciprocal of bin width: bin_index = floor((value - offset) * scale)
    scale: f64,
    /// Value corresponding to bin 0's lower edge (warmup min)
    offset: f64,
    /// Counts for non-negative bin indices [0..)
    pos: Vec<u64>,
    /// Counts for negative bin indices; neg[i] = count for bin_index -(i+1)
    neg: Vec<u64>,
    /// Total count of values recorded (post-warmup only; pre-warmup values
    /// are counted when the buffer is flushed).
    count: u64,
}

impl BinnedHistogram {
    pub fn new() -> Self {
        Self {
            warmup_buffer: Some(Vec::with_capacity(WARMUP_COUNT)),
            bins: Bins {
                scale: 1.0,
                offset: 0.0,
                pos: Vec::new(),
                neg: Vec::new(),
                count: 0,
            },
        }
    }

    /// Record a single value.
    #[inline]
    pub fn add(&mut self, value: f64) {
        if let Some(buf) = &mut self.warmup_buffer {
            buf.push(value);
            if buf.len() >= WARMUP_COUNT {
                self.flush_warmup();
            }
        } else {
            self.bins.record(value);
        }
    }

    /// Ensure the warmup buffer has been flushed (e.g. before querying quantiles).
    fn ensure_flushed(&mut self) {
        if self.warmup_buffer.is_some() {
            self.flush_warmup();
        }
    }

    fn flush_warmup(&mut self) {
        if let Some(buffer) = self.warmup_buffer.take() {
            self.bins.init_from_warmup(buffer);
        }
    }

    /// Estimate the value at quantile `q` (0.0 – 1.0).
    pub fn quantile(&mut self, q: f64) -> Option<f64> {
        self.ensure_flushed();
        self.bins.quantile(q)
    }

    /// Iterate over all non-zero bins in value order, yielding `(midpoint, count)` pairs.
    ///
    /// Flushes any pending warmup buffer first.
    pub fn bins(&mut self) -> impl Iterator<Item = (f64, u64)> + '_ {
        self.ensure_flushed();
        self.bins.iter_bins()
    }

    #[cfg(test)]    
    pub fn count(&self) -> u64 {
        match &self.warmup_buffer {
            Some(buf) => buf.len() as u64,
            None => self.bins.count,
        }
    }
}

impl Default for BinnedHistogram {
    fn default() -> Self {
        Self::new()
    }
}

// ── Bins implementation ─────────────────────────────────────────────────

impl Bins {
    /// Establish the bin layout from warmup data and replay all values.
    fn init_from_warmup(&mut self, values: Vec<f64>) {
        if values.is_empty() {
            return;
        }

        let min_val = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max_val = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let range = max_val - min_val;

        self.offset = min_val;
        self.scale = if range > 0.0 {
            NUM_WARMUP_BINS as f64 / range
        } else {
            1.0
        };

        // Pre-allocate for the warmup range (indices 0..=NUM_WARMUP_BINS)
        self.pos = vec![0u64; NUM_WARMUP_BINS + 1];

        for v in values {
            self.record(v);
        }
    }

    #[inline]
    fn bin_index(&self, value: f64) -> i64 {
        ((value - self.offset) * self.scale).floor() as i64
    }

    /// Recover the midpoint value of a bin from its index.
    #[inline]
    fn bin_midpoint(&self, index: i64) -> f64 {
        (index as f64 + 0.5) / self.scale + self.offset
    }

    #[inline]
    fn record(&mut self, value: f64) {
        let idx = self.bin_index(value);
        if idx >= 0 {
            let i = idx as usize;
            if i >= self.pos.len() {
                self.pos.resize((i + 1).next_power_of_two(), 0);
            }
            self.pos[i] += 1;
        } else {
            let i = (-idx - 1) as usize;
            if i >= self.neg.len() {
                self.neg.resize((i + 1).next_power_of_two(), 0);
            }
            self.neg[i] += 1;
        }
        self.count += 1;
    }

    fn quantile(&self, q: f64) -> Option<f64> {
        if self.count == 0 {
            return None;
        }
        debug_assert!((0.0..=1.0).contains(&q));

        let rank = (q * self.count as f64).ceil().max(1.0) as u64;
        let mut running = 0u64;

        // Walk negative bins: highest index first → most-negative values first
        for i in (0..self.neg.len()).rev() {
            let c = self.neg[i];
            if c == 0 {
                continue;
            }
            running += c;
            if running >= rank {
                let bin_key = -(i as i64) - 1;
                return Some(self.bin_midpoint(bin_key));
            }
        }

        // Walk positive bins: lowest index first → smallest values first
        for (i, &c) in self.pos.iter().enumerate() {
            if c == 0 {
                continue;
            }
            running += c;
            if running >= rank {
                return Some(self.bin_midpoint(i as i64));
            }
        }

        // Fallback: return midpoint of last occupied positive bin
        None
    }

    /// Iterate all non-zero bins in ascending value order: `(midpoint, count)`.
    fn iter_bins(&self) -> impl Iterator<Item = (f64, u64)> + '_ {
        // Negative bins: highest array index = most-negative values first
        let neg_iter = (0..self.neg.len())
            .rev()
            .filter(|&i| self.neg[i] > 0)
            .map(|i| {
                let bin_key = -(i as i64) - 1;
                (self.bin_midpoint(bin_key), self.neg[i])
            });
        // Positive bins: ascending
        let pos_iter = self.pos.iter().enumerate().filter(|&(_, &c)| c > 0).map(
            |(i, &c)| (self.bin_midpoint(i as i64), c),
        );
        neg_iter.chain(pos_iter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let mut h = BinnedHistogram::new();
        assert_eq!(h.count(), 0);
        assert_eq!(h.quantile(0.5), None);
    }

    #[test]
    fn test_single_value() {
        let mut h = BinnedHistogram::new();
        h.add(42.0);
        assert_eq!(h.count(), 1);
        let q = h.quantile(0.5).unwrap();
        // Single value in warmup-flushed bin → should be close to 42
        assert!((q - 42.0).abs() < 1.0, "q={}", q);
    }

    #[test]
    fn test_uniform() {
        let mut h = BinnedHistogram::new();
        for i in 1..=2000 {
            h.add(i as f64);
        }
        assert_eq!(h.count(), 2000);

        let median = h.quantile(0.5).unwrap();
        let err = (median - 1000.0).abs() / 1000.0;
        assert!(err < 0.05, "median={}, err={}", median, err);

        let p90 = h.quantile(0.9).unwrap();
        let err = (p90 - 1800.0).abs() / 1800.0;
        assert!(err < 0.05, "p90={}, err={}", p90, err);
    }

    #[test]
    fn test_negative_values() {
        let mut h = BinnedHistogram::new();
        for i in -500..=500 {
            h.add(i as f64);
        }
        assert_eq!(h.count(), 1001);
        let median = h.quantile(0.5).unwrap();
        assert!(median.abs() < 20.0, "median={}", median);
    }

    #[test]
    fn test_values_outside_warmup_range() {
        let mut h = BinnedHistogram::new();
        // Warmup on 0..1000
        for i in 0..WARMUP_COUNT {
            h.add(i as f64);
        }
        // Now add values far outside
        h.add(-100.0);
        h.add(5000.0);
        assert_eq!(h.count(), WARMUP_COUNT as u64 + 2);
        // Should not panic, and quantile should still work
        let _ = h.quantile(0.5).unwrap();
    }

    #[test]
    fn test_bins_iterator() {
        let mut h = BinnedHistogram::new();
        for i in 1..=2000 {
            h.add(i as f64);
        }
        let bins: Vec<_> = h.bins().collect();
        // All bins should have positive counts
        assert!(bins.iter().all(|&(_, c)| c > 0));
        // Midpoints should be in ascending order
        for w in bins.windows(2) {
            assert!(w[0].0 < w[1].0, "{} >= {}", w[0].0, w[1].0);
        }
        // Total count must equal records added
        let total: u64 = bins.iter().map(|&(_, c)| c).sum();
        assert_eq!(total, 2000);
    }

    #[test]
    fn test_constant_value() {
        let mut h = BinnedHistogram::new();
        for _ in 0..2000 {
            h.add(7.0);
        }
        let q = h.quantile(0.5).unwrap();
        assert!((q - 7.0).abs() < 1.0, "q={}", q);
    }
}