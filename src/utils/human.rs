pub struct Human<T>(pub T);

macro_rules! impl_human_unsigned {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Human<$t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let mut n = self.0;
                let mut j = 0;
                let mut parts: [$t; 7] = [0; 7];
                while n >= 1000 {
                    parts[j] = n % 1000;
                    n /= 1000;
                    j += 1;
                }
                parts[j] = n;
                for i in (0..=j).rev() {
                    if i < j {
                        write!(f, ",{:03}", parts[i])?;
                    } else {
                        write!(f, "{}", parts[i])?;
                    }
                }
                Ok(())
            }
        }
    )+ };
}

macro_rules! impl_human_signed {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Human<$t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let mut n = self.0;
                if n < 0 {
                    write!(f, "-")?;
                }
                let mut j = 0;
                let mut parts: [$t; 7] = [0; 7];
                while n >= 1000 || n <= -1000 {
                    // n % 1000 can be negative; .abs() is safe since |remainder| <= 999
                    parts[j] = (n % 1000).abs();
                    n /= 1000;
                    j += 1;
                }
                // Final quotient has |n| < 1000, so .abs() won't overflow
                parts[j] = n.abs();
                for i in (0..=j).rev() {
                    if i < j {
                        write!(f, ",{:03}", parts[i])?;
                    } else {
                        write!(f, "{}", parts[i])?;
                    }
                }
                Ok(())
            }
        }
    )+ };
}

macro_rules! impl_human_float {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Human<$t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let n = self.0;
                let abs = if n < 0.0 { -n } else { n };
                let sign = if n < 0.0 { "-" } else { "" };
                if abs >= 1e9 as $t {
                    write!(f, "{sign}{:.2}G", abs / 1e9 as $t)
                } else if abs >= 1e6 as $t {
                    write!(f, "{sign}{:.2}M", abs / 1e6 as $t)
                } else if abs >= 1e3 as $t {
                    write!(f, "{sign}{:.2}K", abs / 1e3 as $t)
                } else {
                    write!(f, "{}", n)
                }
            }
        }
    )+ };
}

impl_human_unsigned!(u16, u32, u64, u128, usize);
impl_human_signed!(i16, i32, i64, i128, isize);
impl_human_float!(f32, f64);

/// Convenience trait for `.human()` syntax.
pub trait HumanReadable: Copy {
    fn human(self) -> Human<Self>;
}

macro_rules! impl_human_readable {
    ($($t:ty),+) => { $(
        impl HumanReadable for $t {
            #[inline]
            fn human(self) -> Human<$t> {
                Human(self)
            }
        }
    )+ };
}

impl_human_readable!(
    u16, u32, u64, u128, usize,
    i16, i32, i64, i128, isize,
    f32, f64
);

#[cfg(test)]
mod tests {
    use super::*;

    // ── Unsigned integers ───────────────────────────────────────────
    #[test]
    fn test_zero() {
        assert_eq!(format!("{}", Human(0u64)), "0");
    }

    #[test]
    fn test_small() {
        assert_eq!(format!("{}", Human(42u32)), "42");
        assert_eq!(format!("{}", Human(999u32)), "999");
    }

    #[test]
    fn test_thousands() {
        assert_eq!(format!("{}", Human(1_000u64)), "1,000");
        assert_eq!(format!("{}", Human(1_234u64)), "1,234");
        assert_eq!(format!("{}", Human(999_999u64)), "999,999");
    }

    #[test]
    fn test_millions() {
        assert_eq!(format!("{}", Human(1_000_000u64)), "1,000,000");
        assert_eq!(format!("{}", Human(1_234_567u64)), "1,234,567");
    }

    #[test]
    fn test_large() {
        assert_eq!(
            format!("{}", Human(1_234_567_890_123u64)),
            "1,234,567,890,123"
        );
    }

    #[test]
    fn test_usize() {
        assert_eq!(format!("{}", Human(12_345usize)), "12,345");
    }

    // ── Signed integers ─────────────────────────────────────────────
    #[test]
    fn test_negative_small() {
        assert_eq!(format!("{}", Human(-42i32)), "-42");
        assert_eq!(format!("{}", Human(-999i32)), "-999");
    }

    #[test]
    fn test_negative_thousands() {
        assert_eq!(format!("{}", Human(-1_234i32)), "-1,234");
        assert_eq!(format!("{}", Human(-999_999i64)), "-999,999");
    }

    #[test]
    fn test_negative_millions() {
        assert_eq!(format!("{}", Human(-1_234_567i64)), "-1,234,567");
    }

    #[test]
    fn test_i32_extremes() {
        assert_eq!(format!("{}", Human(i32::MAX)), "2,147,483,647");
        assert_eq!(format!("{}", Human(i32::MIN)), "-2,147,483,648");
    }

    #[test]
    fn test_i64_extremes() {
        assert_eq!(
            format!("{}", Human(i64::MAX)),
            "9,223,372,036,854,775,807"
        );
        assert_eq!(
            format!("{}", Human(i64::MIN)),
            "-9,223,372,036,854,775,808"
        );
    }

    #[test]
    fn test_signed_zero() {
        assert_eq!(format!("{}", Human(0i32)), "0");
    }

    // ── Floats ──────────────────────────────────────────────────────
    #[test]
    fn test_float_small() {
        assert_eq!(format!("{}", Human(42.5f64)), "42.5");
    }

    #[test]
    fn test_float_thousands() {
        assert_eq!(format!("{}", Human(1_500.0f64)), "1.50K");
        assert_eq!(format!("{}", Human(12_345.0f64)), "12.35K");
    }

    #[test]
    fn test_float_millions() {
        assert_eq!(format!("{}", Human(2_500_000.0f64)), "2.50M");
    }

    #[test]
    fn test_float_billions() {
        assert_eq!(format!("{}", Human(3_000_000_000.0f64)), "3.00G");
    }

    #[test]
    fn test_float_negative() {
        assert_eq!(format!("{}", Human(-1_500.0f64)), "-1.50K");
        assert_eq!(format!("{}", Human(-2_500_000.0f64)), "-2.50M");
        assert_eq!(format!("{}", Human(-42.5f64)), "-42.5");
    }

    #[test]
    fn test_f32() {
        assert_eq!(format!("{}", Human(1_500.0f32)), "1.50K");
    }

    #[test]
    fn test_human_trait() {
        assert_eq!(format!("{}", 1_234_567u64.human()), "1,234,567");
        assert_eq!(format!("{}", (-42i32).human()), "-42");
        assert_eq!(format!("{}", 2_500_000.0f64.human()), "2.50M");
    }
}