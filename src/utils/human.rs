#![allow(dead_code)] // This module is used in some configurations but not others

/// Human-readable display with K/M/G suffixes (floats, or integers promoted to f64).
pub struct Human<T>(pub T);

/// Comma-separated integer display (e.g. 1,234,567).
pub struct Commas<T>(pub T);

const SUFFIXES: [&str; 8] = ["", "K", "M", "G", "T", "P", "E", "Z"];

// ── Commas: unsigned ────────────────────────────────────────────────

macro_rules! impl_commas_unsigned {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Commas<$t> {
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

// ── Commas: signed ──────────────────────────────────────────────────

macro_rules! impl_commas_signed {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Commas<$t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let mut n = self.0;
                if n < 0 {
                    write!(f, "-")?;
                }
                let mut j = 0;
                let mut parts: [$t; 7] = [0; 7];
                while n >= 1000 || n <= -1000 {
                    parts[j] = (n % 1000).abs();
                    n /= 1000;
                    j += 1;
                }
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

impl_commas_unsigned!(u16, u32, u64, u128, usize);
impl_commas_signed!(i16, i32, i64, i128, isize);

// ── Human: float (K/M/G suffixes) ──────────────────────────────────

macro_rules! impl_human_float {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Human<$t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let n = self.0;
                if n.is_nan() || n.is_infinite() {
                    return write!(f, "{}", n);
                }

                let mut abs = if n < 0.0 { -n } else { n };
                if abs < 1000.0 {
                    return write!(f, "{}", n);
                }

                let sign = if n < 0.0 { "-" } else { "" };

                let mut i = 0;
                while abs >= 1000.0 && i < SUFFIXES.len() - 1 {
                    abs /= 1000.0;
                    i += 1;
                }
                write!(f, "{sign}{:.2}{}", abs, SUFFIXES[i])
            }
        }
    )+ };
}

impl_human_float!(f32, f64);

// ── Human: integers → promote to f64 and display with K/M/G ────────

macro_rules! impl_human_int {
    ($($t:ty),+) => { $(
        impl std::fmt::Display for Human<$t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                Human(self.0 as f64).fmt(f)
            }
        }
    )+ };
}

impl_human_int!(u16, u32, u64, u128, usize, i16, i32, i64, i128, isize);

// ── Convenience traits ──────────────────────────────────────────────

/// `.commas()` → `Commas<T>` (comma-separated integer display).
pub trait CommaReadable: Copy {
    fn commas(self) -> Commas<Self>;
}

macro_rules! impl_comma_readable {
    ($($t:ty),+) => { $(
        impl CommaReadable for $t {
            #[inline]
            fn commas(self) -> Commas<$t> {
                Commas(self)
            }
        }
    )+ };
}

impl_comma_readable!(u16, u32, u64, u128, usize, i16, i32, i64, i128, isize);

/// `.human()` → `Human<T>` (K/M/G for large values).
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
    u16, u32, u64, u128, usize, i16, i32, i64, i128, isize, f32, f64
);

#[cfg(test)]
mod tests {
    use super::*;

    // ── Commas: unsigned integers ───────────────────────────────────
    #[test]
    fn test_commas_zero() {
        assert_eq!(format!("{}", Commas(0u64)), "0");
    }

    #[test]
    fn test_commas_small() {
        assert_eq!(format!("{}", Commas(42u32)), "42");
        assert_eq!(format!("{}", Commas(999u32)), "999");
    }

    #[test]
    fn test_commas_thousands() {
        assert_eq!(format!("{}", Commas(1_000u64)), "1,000");
        assert_eq!(format!("{}", Commas(1_234u64)), "1,234");
        assert_eq!(format!("{}", Commas(999_999u64)), "999,999");
    }

    #[test]
    fn test_commas_millions() {
        assert_eq!(format!("{}", Commas(1_000_000u64)), "1,000,000");
        assert_eq!(format!("{}", Commas(1_234_567u64)), "1,234,567");
    }

    #[test]
    fn test_commas_large() {
        assert_eq!(
            format!("{}", Commas(1_234_567_890_123u64)),
            "1,234,567,890,123"
        );
    }

    #[test]
    fn test_commas_usize() {
        assert_eq!(format!("{}", Commas(12_345usize)), "12,345");
    }

    // ── Commas: signed integers ─────────────────────────────────────
    #[test]
    fn test_commas_negative_small() {
        assert_eq!(format!("{}", Commas(-42i32)), "-42");
        assert_eq!(format!("{}", Commas(-999i32)), "-999");
    }

    #[test]
    fn test_commas_negative_thousands() {
        assert_eq!(format!("{}", Commas(-1_234i32)), "-1,234");
        assert_eq!(format!("{}", Commas(-999_999i64)), "-999,999");
    }

    #[test]
    fn test_commas_negative_millions() {
        assert_eq!(format!("{}", Commas(-1_234_567i64)), "-1,234,567");
    }

    #[test]
    fn test_commas_i32_extremes() {
        assert_eq!(format!("{}", Commas(i32::MAX)), "2,147,483,647");
        assert_eq!(format!("{}", Commas(i32::MIN)), "-2,147,483,648");
    }

    #[test]
    fn test_commas_i64_extremes() {
        assert_eq!(format!("{}", Commas(i64::MAX)), "9,223,372,036,854,775,807");
        assert_eq!(
            format!("{}", Commas(i64::MIN)),
            "-9,223,372,036,854,775,808"
        );
    }

    #[test]
    fn test_commas_signed_zero() {
        assert_eq!(format!("{}", Commas(0i32)), "0");
    }

    // ── Human: floats ───────────────────────────────────────────────
    #[test]
    fn test_human_float_small() {
        assert_eq!(format!("{}", Human(42.5f64)), "42.5");
    }

    #[test]
    fn test_human_float_thousands() {
        assert_eq!(format!("{}", Human(1_500.0f64)), "1.50K");
        assert_eq!(format!("{}", Human(12_345.0f64)), "12.35K");
    }

    #[test]
    fn test_human_float_millions() {
        assert_eq!(format!("{}", Human(2_500_000.0f64)), "2.50M");
    }

    #[test]
    fn test_human_float_billions() {
        assert_eq!(format!("{}", Human(3_000_000_000.0f64)), "3.00G");
    }

    #[test]
    fn test_human_float_negative() {
        assert_eq!(format!("{}", Human(-1_500.0f64)), "-1.50K");
        assert_eq!(format!("{}", Human(-2_500_000.0f64)), "-2.50M");
        assert_eq!(format!("{}", Human(-42.5f64)), "-42.5");
    }

    #[test]
    fn test_human_f32() {
        assert_eq!(format!("{}", Human(1_500.0f32)), "1.50K");
    }

    // ── Human: integers (promoted to f64) ───────────────────────────
    #[test]
    fn test_human_int_small() {
        assert_eq!(format!("{}", Human(42u64)), "42");
    }

    #[test]
    fn test_human_int_thousands() {
        assert_eq!(format!("{}", Human(1_500u64)), "1.50K");
        assert_eq!(format!("{}", Human(12_345u64)), "12.35K");
    }

    #[test]
    fn test_human_int_millions() {
        assert_eq!(format!("{}", Human(2_500_000u64)), "2.50M");
    }

    #[test]
    fn test_human_int_negative() {
        assert_eq!(format!("{}", Human(-1_500i64)), "-1.50K");
    }

    // ── Trait methods ───────────────────────────────────────────────
    #[test]
    fn test_commas_trait() {
        assert_eq!(format!("{}", 1_234_567u64.commas()), "1,234,567");
        assert_eq!(format!("{}", (-42i32).commas()), "-42");
    }

    #[test]
    fn test_human_trait() {
        assert_eq!(format!("{}", 2_500_000u64.human()), "2.50M");
        assert_eq!(format!("{}", (-42i32).human()), "-42");
        assert_eq!(format!("{}", 2_500_000.0f64.human()), "2.50M");
    }
}
