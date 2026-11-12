//! Type-safe score wrappers for alignment scoring.
//!
//! These newtypes prevent accidental mixing of scores with different semantics:
//! - `DivergenceScore`: Lower is better (edit distance, gap penalties)
//! - `QualityScore`: Higher is better (alignment scores, ranking)

use std::cmp::Ordering;

/// A divergence/distance metric where **lower is better**.
///
/// Used for:
/// - Edit distance
/// - Gap penalties  
/// - Mismatch counts
///
/// Always non-negative (clamped at 0.0).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DivergenceScore(pub f64);

impl DivergenceScore {
    /// Zero divergence (perfect match).
    pub const ZERO: Self = Self(0.0);

    /// Infinite divergence (no alignment possible).
    pub const INFINITY: Self = Self(f64::INFINITY);

    /// Create a new divergence score, clamping negative values to 0.
    pub fn new(value: f64) -> Self {
        if value < 0.0 {
            panic!("DivergenceScore cannot be negative: {}", value);
        }
        Self(value)
    }

    /// Get the underlying value.
    pub fn value(self) -> f64 {
        self.0
    }

    /// Returns true if this score is better than (less than) other.
    pub fn is_better_than(self, other: Self) -> bool {
        self.0 < other.0
    }

    /// Returns true if this score is worse than (greater than) other.
    pub fn is_worse_than(self, other: Self) -> bool {
        self.0 > other.0
    }
}

impl PartialOrd for DivergenceScore {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl From<u32> for DivergenceScore {
    fn from(value: u32) -> Self {
        Self::new(value as f64)
    }
}

impl From<i32> for DivergenceScore {
    fn from(value: i32) -> Self {
        Self::new(value as f64)
    }
}

impl From<f64> for DivergenceScore {
    fn from(value: f64) -> Self {
        Self::new(value)
    }
}

impl std::fmt::Display for DivergenceScore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.2}d", self.0)
    }
}

/// A quality/alignment score where **higher is better**.
///
/// Used for:
/// - Alignment scores (matches - penalties)
/// - Ranking scores
/// - Cluster quality
///
/// Can be negative (when penalties exceed matches).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct QualityScore(pub f64);

impl QualityScore {
    /// Zero quality score.
    pub const ZERO: Self = Self(0.0);

    /// Negative infinity (worst possible).
    pub const NEG_INFINITY: Self = Self(f64::NEG_INFINITY);

    /// Positive infinity (best possible).
    pub const INFINITY: Self = Self(f64::INFINITY);

    /// Create a new quality score.
    pub fn new(value: f64) -> Self {
        Self(value)
    }

    /// Get the underlying value.
    pub fn value(self) -> f64 {
        self.0
    }

    /// Returns true if this score is better than (greater than) other.
    pub fn is_better_than(self, other: Self) -> bool {
        self.0 > other.0
    }

    /// Returns true if this score is worse than (less than) other.
    pub fn is_worse_than(self, other: Self) -> bool {
        self.0 < other.0
    }
}

impl PartialOrd for QualityScore {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl From<i32> for QualityScore {
    fn from(value: i32) -> Self {
        Self(value as f64)
    }
}

impl From<i64> for QualityScore {
    fn from(value: i64) -> Self {
        Self(value as f64)
    }
}

impl From<f64> for QualityScore {
    fn from(value: f64) -> Self {
        Self(value)
    }
}

impl std::fmt::Display for QualityScore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.2}q", self.0)
    }
}

pub fn compute_mapq_from_diff(
    best_score: f64,
    second_best_score: Option<f64>,
    num_seeds: usize,
    _scale: f64,  // Calibration parameter, typically 10-30
) -> u8 {
    let s1 = best_score.max(1.0);
    let s2 = second_best_score.unwrap_or(0.0).max(0.0);
    let r = s2 / s1;
    let m = num_seeds as f64 / 10.0;
    let res = 60.0 * (1.0 - r) * m * s1.log10();
    log::debug!(
        "Computing MAPQ: best_score={}, second_best_score={:?}, num_seeds={}, r={:.3}, m={:.3}, raw_res={:.2}",
        best_score,
        second_best_score,
        num_seeds,
        r,
        m,
        res
    );
    res.max(0.0).min(100.0).round() as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "cannot be negative")]
    fn test_divergence_rejects_negative() {
        DivergenceScore::new(-5.0);
    }

    #[test]
    fn test_divergence_accepts_valid() {
        assert_eq!(DivergenceScore::new(10.0).value(), 10.0);
        assert_eq!(DivergenceScore::new(0.0).value(), 0.0);
    }

    #[test]
    fn test_divergence_ordering() {
        let low = DivergenceScore::new(1.0);
        let high = DivergenceScore::new(10.0);
        assert!(low < high); // Lower divergence is better, but sorts first
    }

    #[test]
    fn test_quality_ordering() {
        let low = QualityScore::new(1.0);
        let high = QualityScore::new(10.0);
        assert!(high > low); // Higher quality is better
        assert!(high.is_better_than(low));
        assert!(low.is_worse_than(high));
    }

    #[test]
    fn test_conversions() {
        let d: DivergenceScore = 42u32.into();
        assert_eq!(d.value(), 42.0);

        let q: QualityScore = (-10i32).into();
        assert_eq!(q.value(), -10.0);
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", DivergenceScore::new(3.5)), "3.50d");
        assert_eq!(format!("{}", QualityScore::new(-2.5)), "-2.50q");
    }
}
