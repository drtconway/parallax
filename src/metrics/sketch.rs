#![allow(dead_code)]
//! DDSketch implementation for quantile estimation.
//!
//! DDSketch (Distributed Distribution Sketch) is a data structure that provides
//! accurate quantile estimation with a relative error guarantee. This implementation
//! is based on the paper "DDSketch: A Fast and Fully-Mergeable Quantile Sketch with
//! Relative-Error Guarantees" by Masson, Rim, and Lee.
//!
//! The key property is that for any quantile q, the estimated value v̂ satisfies:
//! |v̂ - v| ≤ α * v, where α is the relative accuracy parameter.

use std::collections::BTreeMap;

/// A DDSketch for estimating quantiles with relative error guarantees.
///
/// The sketch uses a logarithmic mapping to assign values to buckets,
/// which provides the relative error guarantee.
#[derive(Clone)]
pub struct DDSketch {
    /// Relative accuracy parameter (0 < alpha < 1)
    alpha: f64,
    /// Precomputed gamma = (1 + alpha) / (1 - alpha)
    gamma: f64,
    /// Precomputed ln(gamma) for bucket index calculation
    ln_gamma: f64,
    /// Buckets for positive values: index -> count
    positive_buckets: BTreeMap<i32, u64>,
    /// Buckets for negative values: index -> count
    negative_buckets: BTreeMap<i32, u64>,
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
            positive_buckets: BTreeMap::new(),
            negative_buckets: BTreeMap::new(),
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

    /// Add a value to the sketch.
    pub fn add(&mut self, value: f64) {
        self.count += 1;
        self.sum += value;
        self.min = self.min.min(value);
        self.max = self.max.max(value);

        if value > 0.0 {
            let index = self.bucket_index(value);
            *self.positive_buckets.entry(index).or_insert(0) += 1;
        } else if value < 0.0 {
            let index = self.bucket_index(-value);
            *self.negative_buckets.entry(index).or_insert(0) += 1;
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

        // First, go through negative buckets (in descending order of absolute value)
        for (&index, &bucket_count) in self.negative_buckets.iter().rev() {
            running_count += bucket_count;
            if running_count >= rank {
                return Some(-self.bucket_value(index));
            }
        }

        // Then zeros
        running_count += self.zero_count;
        if running_count >= rank {
            return Some(0.0);
        }

        // Finally, positive buckets (in ascending order)
        for (&index, &bucket_count) in self.positive_buckets.iter() {
            running_count += bucket_count;
            if running_count >= rank {
                return Some(self.bucket_value(index));
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

        for (&index, &count) in &other.positive_buckets {
            *self.positive_buckets.entry(index).or_insert(0) += count;
        }

        for (&index, &count) in &other.negative_buckets {
            *self.negative_buckets.entry(index).or_insert(0) += count;
        }
    }

    /// Get the number of buckets currently in use.
    pub fn bucket_count(&self) -> usize {
        self.positive_buckets.len() + self.negative_buckets.len()
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
        assert!(relative_error <= 0.02, "median={}, expected={}, error={}", median, expected, relative_error);

        // Check p99 (should be around 990)
        let p99 = sketch.quantile(0.99).unwrap();
        let expected = 990.0;
        let relative_error = (p99 - expected).abs() / expected;
        assert!(relative_error <= 0.02, "p99={}, expected={}, error={}", p99, expected, relative_error);
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
