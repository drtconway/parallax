#![allow(dead_code)]
//! DDSketch implementation for quantile estimation, backed by flat arrays.
//!
//! DDSketch (Distributed Distribution Sketch) is a data structure that provides
//! accurate quantile estimation with a relative error guarantee. This implementation
//! is based on the paper "DDSketch: A Fast and Fully-Mergeable Quantile Sketch with
//! Relative-Error Guarantees" by Masson, Rim, and Lee.
//!
//! The key property is that for any quantile q, the estimated value v̂ satisfies:
//! |v̂ - v| ≤ α * v, where α is the relative accuracy parameter.
//!
//! Bucket counts are stored in contiguous `Vec<u64>` arrays rather than `BTreeMap`s,
//! giving O(1) amortized recording and cache-friendly sequential access during
//! quantile queries.

/// A DDSketch for estimating quantiles with relative error guarantees.
///
/// The sketch uses a logarithmic mapping to assign values to buckets,
/// which provides the relative error guarantee. Bucket counts are stored
/// in flat arrays indexed by `(bucket_key - offset)`, enabling direct
/// indexing instead of tree lookups.
#[derive(Clone)]
pub struct DDSketch {
    /// Relative accuracy parameter (0 < alpha < 1)
    alpha: f64,
    /// Precomputed gamma = (1 + alpha) / (1 - alpha)
    gamma: f64,
    /// Precomputed ln(gamma) for bucket index calculation
    ln_gamma: f64,
    /// Counts for positive-value buckets. Index i → bucket key (i as i32 + pos_offset).
    pos_buckets: Vec<u64>,
    /// Bucket key corresponding to pos_buckets[0].
    pos_offset: i32,
    /// Counts for negative-value buckets (keyed on absolute value).
    neg_buckets: Vec<u64>,
    /// Bucket key corresponding to neg_buckets[0].
    neg_offset: i32,
    /// Count of zero values
    zero_count: u64,
    /// Total count of all values
    count: u64,
    /// Minimum value seen
    min: f64,
    /// Maximum value seen
    max: f64,
    /// Sum of all values
    sum: f64,
}

/// Initial allocation size for bucket arrays (centered around first key).
const INITIAL_CAPACITY: usize = 64;

impl DDSketch {
    /// Create a new DDSketch with the given relative accuracy.
    ///
    /// # Arguments
    /// * `alpha` - Relative accuracy parameter (0 < alpha < 1).
    ///   Smaller values give more accurate quantiles but use more memory.
    ///   Common values: 0.01 (1% error), 0.005 (0.5% error)
    ///
    /// # Panics
    /// Panics if alpha is not in (0, 1).
    pub fn new(alpha: f64) -> Self {
        assert!(alpha > 0.0 && alpha < 1.0, "alpha must be in (0, 1)");

        let gamma = (1.0 + alpha) / (1.0 - alpha);
        let ln_gamma = gamma.ln();

        Self {
            alpha,
            gamma,
            ln_gamma,
            pos_buckets: Vec::new(),
            pos_offset: 0,
            neg_buckets: Vec::new(),
            neg_offset: 0,
            zero_count: 0,
            count: 0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            sum: 0.0,
        }
    }

    /// Create a DDSketch with default 1% relative accuracy.
    pub fn with_default_accuracy() -> Self {
        Self::new(0.01)
    }

    /// Map a value to its bucket index.
    #[inline]
    fn bucket_index(&self, value: f64) -> i32 {
        // Using the logarithmic mapping: index = ceil(log_gamma(value))
        (value.ln() / self.ln_gamma).ceil() as i32
    }

    /// Get the representative value for a bucket index.
    fn bucket_value(&self, index: i32) -> f64 {
        // The representative value is 2 * gamma^index / (1 + gamma)
        // This is the geometric midpoint of the bucket
        2.0 * self.gamma.powi(index) / (1.0 + self.gamma)
    }

    /// Ensure the bucket array can hold `bucket_key`, returning the array index.
    ///
    /// On first insertion, allocates `INITIAL_CAPACITY` slots centered around `bucket_key`.
    /// Subsequent out-of-range keys trigger geometric growth via [`Self::grow`].
    #[inline]
    fn ensure_index(buckets: &mut Vec<u64>, offset: &mut i32, bucket_key: i32) -> usize {
        if !buckets.is_empty() {
            let rel = bucket_key - *offset;
            if rel >= 0 && (rel as usize) < buckets.len() {
                return rel as usize;
            }
            return Self::grow(buckets, offset, bucket_key);
        }
        // First insertion: center initial allocation around this key
        *offset = bucket_key - (INITIAL_CAPACITY as i32 / 2);
        buckets.resize(INITIAL_CAPACITY, 0);
        (bucket_key - *offset) as usize
    }

    /// Cold path: grow the bucket array to include `bucket_key`.
    #[cold]
    fn grow(buckets: &mut Vec<u64>, offset: &mut i32, bucket_key: i32) -> usize {
        let old_len = buckets.len();
        if bucket_key < *offset {
            // Extend below: shift existing data rightward
            let deficit = (*offset - bucket_key) as usize;
            let new_len = (old_len + deficit).next_power_of_two().max(INITIAL_CAPACITY);
            let shift = new_len - old_len;
            buckets.resize(new_len, 0);
            buckets.copy_within(0..old_len, shift);
            buckets[..shift].fill(0);
            *offset -= shift as i32;
        } else {
            // Extend above
            let needed = (bucket_key - *offset) as usize + 1;
            let new_len = needed.next_power_of_two().max(INITIAL_CAPACITY);
            buckets.resize(new_len, 0);
        }
        (bucket_key - *offset) as usize
    }

    /// Add a value to the sketch.
    #[inline]
    pub fn add(&mut self, value: f64) {
        self.count += 1;
        self.sum += value;
        self.min = self.min.min(value);
        self.max = self.max.max(value);

        if value > 0.0 {
            let key = self.bucket_index(value);
            let idx = Self::ensure_index(&mut self.pos_buckets, &mut self.pos_offset, key);
            self.pos_buckets[idx] += 1;
        } else if value < 0.0 {
            let key = self.bucket_index(-value);
            let idx = Self::ensure_index(&mut self.neg_buckets, &mut self.neg_offset, key);
            self.neg_buckets[idx] += 1;
        } else {
            self.zero_count += 1;
        }
    }

    /// Get the total count of values added.
    pub fn count(&self) -> u64 {
        self.count
    }

    /// Get the minimum value.
    pub fn min(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.min)
        }
    }

    /// Get the maximum value.
    pub fn max(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.max)
        }
    }

    /// Get the sum of all values.
    pub fn sum(&self) -> f64 {
        self.sum
    }

    /// Get the mean of all values.
    pub fn mean(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.sum / self.count as f64)
        }
    }

    /// Estimate the value at the given quantile.
    ///
    /// # Arguments
    /// * `q` - Quantile to estimate, in [0, 1]. For example, 0.5 for median,
    ///   0.95 for 95th percentile.
    ///
    /// # Returns
    /// The estimated value at the quantile, or None if the sketch is empty.
    pub fn quantile(&self, q: f64) -> Option<f64> {
        if self.count == 0 {
            return None;
        }

        assert!(q >= 0.0 && q <= 1.0, "quantile must be in [0, 1]");

        // Handle edge cases
        if q == 0.0 {
            return Some(self.min);
        }
        if q == 1.0 {
            return Some(self.max);
        }

        // Target rank (1-indexed)
        let rank = (q * self.count as f64).ceil() as u64;
        let mut running_count = 0u64;

        // Negative buckets: descending bucket key = most negative values first.
        // Highest array index = highest bucket key = largest absolute value.
        for i in (0..self.neg_buckets.len()).rev() {
            let c = self.neg_buckets[i];
            if c == 0 {
                continue;
            }
            running_count += c;
            if running_count >= rank {
                return Some(-self.bucket_value(i as i32 + self.neg_offset));
            }
        }

        // Then zeros
        running_count += self.zero_count;
        if running_count >= rank {
            return Some(0.0);
        }

        // Positive buckets: ascending bucket key = smallest positive values first.
        for (i, &c) in self.pos_buckets.iter().enumerate() {
            if c == 0 {
                continue;
            }
            running_count += c;
            if running_count >= rank {
                return Some(self.bucket_value(i as i32 + self.pos_offset));
            }
        }

        // Should not reach here, but return max as fallback
        Some(self.max)
    }

    /// Merge another DDSketch into this one.
    ///
    /// # Panics
    /// Panics if the sketches have different alpha values.
    pub fn merge(&mut self, other: &DDSketch) {
        assert!(
            (self.alpha - other.alpha).abs() < f64::EPSILON,
            "Cannot merge sketches with different alpha values"
        );

        if other.count == 0 {
            return;
        }

        self.count += other.count;
        self.sum += other.sum;
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
        self.zero_count += other.zero_count;

        Self::merge_buckets(
            &mut self.pos_buckets,
            &mut self.pos_offset,
            &other.pos_buckets,
            other.pos_offset,
        );
        Self::merge_buckets(
            &mut self.neg_buckets,
            &mut self.neg_offset,
            &other.neg_buckets,
            other.neg_offset,
        );
    }

    /// Merge source bucket array into destination, aligning by bucket key.
    fn merge_buckets(
        dst: &mut Vec<u64>,
        dst_offset: &mut i32,
        src: &[u64],
        src_offset: i32,
    ) {
        for (i, &count) in src.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let key = i as i32 + src_offset;
            let idx = Self::ensure_index(dst, dst_offset, key);
            dst[idx] += count;
        }
    }

    /// Get the number of buckets currently in use (non-zero count).
    pub fn bucket_count(&self) -> usize {
        self.pos_buckets.iter().filter(|&&c| c > 0).count()
            + self.neg_buckets.iter().filter(|&&c| c > 0).count()
    }
}

impl Default for DDSketch {
    fn default() -> Self {
        Self::with_default_accuracy()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_sketch() {
        let sketch = DDSketch::new(0.01);
        assert_eq!(sketch.count(), 0);
        assert_eq!(sketch.quantile(0.5), None);
        assert_eq!(sketch.min(), None);
        assert_eq!(sketch.max(), None);
    }

    #[test]
    fn test_single_value() {
        let mut sketch = DDSketch::new(0.01);
        sketch.add(42.0);

        assert_eq!(sketch.count(), 1);
        assert_eq!(sketch.min(), Some(42.0));
        assert_eq!(sketch.max(), Some(42.0));
        assert_eq!(sketch.sum(), 42.0);

        // All quantiles should return approximately 42
        let q = sketch.quantile(0.5).unwrap();
        assert!((q - 42.0).abs() / 42.0 <= 0.01);
    }

    #[test]
    fn test_uniform_distribution() {
        let mut sketch = DDSketch::new(0.01);

        // Add values 1 to 1000
        for i in 1..=1000 {
            sketch.add(i as f64);
        }

        assert_eq!(sketch.count(), 1000);
        assert_eq!(sketch.min(), Some(1.0));
        assert_eq!(sketch.max(), Some(1000.0));

        // Check median (should be around 500)
        let median = sketch.quantile(0.5).unwrap();
        let expected = 500.0;
        let relative_error = (median - expected).abs() / expected;
        assert!(
            relative_error <= 0.02,
            "median={}, expected={}, error={}",
            median,
            expected,
            relative_error
        );

        // Check p99 (should be around 990)
        let p99 = sketch.quantile(0.99).unwrap();
        let expected = 990.0;
        let relative_error = (p99 - expected).abs() / expected;
        assert!(
            relative_error <= 0.02,
            "p99={}, expected={}, error={}",
            p99,
            expected,
            relative_error
        );
    }

    #[test]
    fn test_negative_values() {
        let mut sketch = DDSketch::new(0.01);

        for i in -100..=100 {
            sketch.add(i as f64);
        }

        assert_eq!(sketch.count(), 201);
        assert_eq!(sketch.min(), Some(-100.0));
        assert_eq!(sketch.max(), Some(100.0));

        // Median should be around 0
        let median = sketch.quantile(0.5).unwrap();
        assert!(median.abs() <= 5.0, "median={}", median);
    }

    #[test]
    fn test_merge() {
        let mut sketch1 = DDSketch::new(0.01);
        let mut sketch2 = DDSketch::new(0.01);

        for i in 1..=500 {
            sketch1.add(i as f64);
        }
        for i in 501..=1000 {
            sketch2.add(i as f64);
        }

        sketch1.merge(&sketch2);

        assert_eq!(sketch1.count(), 1000);
        assert_eq!(sketch1.min(), Some(1.0));
        assert_eq!(sketch1.max(), Some(1000.0));
    }
}
