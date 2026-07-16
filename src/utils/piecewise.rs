
/// Piecewise-linear function with `N` slope breakpoints.
///
/// Evaluates as:
///   `f(x) = base * x + sum_i slope_increments[i] * max(0, x - breakpoints[i])`
///
/// `base` is the initial slope (for x from 0 up to `breakpoints[0]`).
/// Each entry in `slope_increments` is the *change* in slope at the
/// corresponding breakpoint — positive to steepen, negative to flatten.
/// The function is always 0 at x ≤ 0.
///
/// # Example
///
/// ```
/// use parallax::utils::piecewise::Piecewise;
///
/// // Three-slope read-gap cost: 0.02/bp up to 15, 0.05/bp up to 60, 0.20/bp beyond.
/// let f = Piecewise {
///     base: 0.02,
///     breakpoints:      [15.0, 60.0],
///     slope_increments: [0.03, 0.15],
/// };
/// assert!((f.eval(0.0)  - 0.0 ).abs() < 1e-9);
/// assert!((f.eval(15.0) - 0.30).abs() < 1e-9);
/// assert!((f.eval(60.0) - 2.55).abs() < 1e-9);
/// assert!((f.eval(70.0) - 4.55).abs() < 1e-9);
/// ```
pub struct Piecewise<const N: usize> {
    pub base: f64,
    pub breakpoints: [f64; N],
    pub slope_increments: [f64; N],
}

impl<const N: usize> Piecewise<N> {
    /// Evaluate the function at `x`. Returns 0 for non-positive `x`.
    pub fn eval(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let mut v = self.base * x;
        for i in 0..N {
            v += self.slope_increments[i] * (x - self.breakpoints[i]).max(0.0);
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::Piecewise;

    const EPS: f64 = 1e-9;

    fn assert_close(got: f64, want: f64) {
        assert!((got - want).abs() < EPS, "got {got} want {want}");
    }

    // Zero breakpoints — just a linear function.
    #[test]
    fn zero_breakpoints_is_linear() {
        let f: Piecewise<0> = Piecewise { base: 0.5, breakpoints: [], slope_increments: [] };
        assert_close(f.eval(0.0), 0.0);
        assert_close(f.eval(10.0), 5.0);
        assert_close(f.eval(100.0), 50.0);
    }

    // Non-positive inputs always return 0.
    #[test]
    fn non_positive_returns_zero() {
        let f: Piecewise<2> = Piecewise {
            base: 0.02,
            breakpoints:      [15.0, 60.0],
            slope_increments: [0.03, 0.15],
        };
        assert_close(f.eval(0.0), 0.0);
        assert_close(f.eval(-1.0), 0.0);
        assert_close(f.eval(-100.0), 0.0);
    }

    // Two-breakpoint read-gap cost (the actual values used in DPConfig::default).
    // Slopes: 0.02 up to 15, 0.05 from 15–60, 0.20 beyond 60.
    #[test]
    fn read_gap_cost_two_breakpoints() {
        let f: Piecewise<2> = Piecewise {
            base: 0.02,
            breakpoints:      [15.0, 60.0],
            slope_increments: [0.03, 0.15],
        };
        // Before first breakpoint: pure base slope.
        assert_close(f.eval(10.0), 0.02 * 10.0);
        // At first breakpoint exactly.
        assert_close(f.eval(15.0), 0.02 * 15.0);
        // Between breakpoints: base + first increment.
        assert_close(f.eval(20.0), 0.02 * 20.0 + 0.03 * 5.0);
        // At second breakpoint exactly.
        assert_close(f.eval(60.0), 0.02 * 60.0 + 0.03 * 45.0);
        // Beyond second breakpoint.
        assert_close(f.eval(70.0), 0.02 * 70.0 + 0.03 * 55.0 + 0.15 * 10.0);
    }

    // One-breakpoint ref-dev cost (DPConfig::default): slope 0.01 up to 50,
    // then 0.001 beyond (slope drops by 0.009).
    #[test]
    fn ref_dev_cost_one_breakpoint_decreasing_slope() {
        let f: Piecewise<1> = Piecewise {
            base: 0.01,
            breakpoints:      [50.0],
            slope_increments: [-0.009],
        };
        // Before breakpoint.
        assert_close(f.eval(30.0), 0.01 * 30.0);
        // At breakpoint.
        assert_close(f.eval(50.0), 0.01 * 50.0);
        // Beyond breakpoint: slope drops to 0.001.
        assert_close(f.eval(100.0), 0.01 * 100.0 + (-0.009) * 50.0);
        assert_close(f.eval(100.0), 0.50 + 0.001 * 50.0);
    }

    // Function is continuous at breakpoints (no jump discontinuities).
    #[test]
    fn continuous_at_breakpoints() {
        let f: Piecewise<2> = Piecewise {
            base: 0.02,
            breakpoints:      [15.0, 60.0],
            slope_increments: [0.03, 0.15],
        };
        let delta = 1e-6;
        let bp1 = 15.0_f64;
        let bp2 = 60.0_f64;
        assert!((f.eval(bp1 - delta) - f.eval(bp1 + delta)).abs() < 1e-4);
        assert!((f.eval(bp2 - delta) - f.eval(bp2 + delta)).abs() < 1e-4);
    }
}
