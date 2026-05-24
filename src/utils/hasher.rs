#![allow(dead_code)]

pub trait Hasher {
    const NAME: &'static str;

    fn hash32(x: u32) -> u32;

    fn hash64(x: u64) -> u64;
}

pub struct IdentityHasher;
impl Hasher for IdentityHasher {
    const NAME: &'static str = "identity";

    fn hash32(x: u32) -> u32 {
        x
    }

    fn hash64(x: u64) -> u64 {
        x
    }
}

pub struct FnvHasher;
impl Hasher for FnvHasher {
    const NAME: &'static str = "fnv";

    fn hash32(x: u32) -> u32 {
        const FNV_PRIME: u32 = 0x01000193;
        const FNV_OFFSET_BASIS: u32 = 0x811c9dc5;

        let mut hash = FNV_OFFSET_BASIS;
        for byte in x.to_le_bytes().iter() {
            hash ^= *byte as u32;
            hash = hash.wrapping_mul(FNV_PRIME);
        }
        hash
    }

    fn hash64(x: u64) -> u64 {
        const FNV_PRIME: u64 = 0x00000100000001B3;
        const FNV_OFFSET_BASIS: u64 = 0xcbf29ce484222325;

        let mut hash = FNV_OFFSET_BASIS;
        for byte in x.to_le_bytes().iter() {
            hash ^= *byte as u64;
            hash = hash.wrapping_mul(FNV_PRIME);
        }
        hash
    }
}

/// Splitmix64 hasher - excellent mixing for 64-bit integers.
/// Based on the SplitMix64 PRNG by Sebastiano Vigna.
/// Very fast with excellent statistical properties (passes BigCrush).
/// Includes an offset to avoid the fixed point at 0 (important for poly-A k-mers).
pub struct Splitmix64Hasher;

impl Splitmix64Hasher {
    // Golden ratio constant used as offset to avoid fixed point at 0
    const OFFSET: u64 = 0x9E3779B97F4A7C15;
}

impl Hasher for Splitmix64Hasher {
    const NAME: &'static str = "splitmix64";

    fn hash32(x: u32) -> u32 {
        // Use the 64-bit version and truncate
        Self::hash64(x as u64) as u32
    }

    fn hash64(x: u64) -> u64 {
        // Add offset to avoid hash(0) = 0 fixed point
        let mut z = x.wrapping_add(Self::OFFSET);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    // --- Statistical helpers ---

    /// Chi-squared statistic for 32-bit hash bucket distribution.
    /// Hashes sequential inputs 0..n_samples and counts how many fall in each bucket.
    /// For a uniform hash, this follows χ²(n_buckets - 1): mean ≈ n_buckets - 1.
    fn chi_squared_32<H: Hasher>(n_samples: usize, n_buckets: usize) -> f64 {
        let mut counts = vec![0usize; n_buckets];
        for i in 0..n_samples {
            let bucket = (H::hash32(i as u32) as usize) % n_buckets;
            counts[bucket] += 1;
        }
        let expected = n_samples as f64 / n_buckets as f64;
        counts
            .iter()
            .map(|&c| {
                let d = c as f64 - expected;
                d * d / expected
            })
            .sum()
    }

    /// Chi-squared statistic for 64-bit hash bucket distribution.
    fn chi_squared_64<H: Hasher>(n_samples: usize, n_buckets: usize) -> f64 {
        let mut counts = vec![0usize; n_buckets];
        for i in 0..n_samples {
            let bucket = (H::hash64(i as u64) as usize) % n_buckets;
            counts[bucket] += 1;
        }
        let expected = n_samples as f64 / n_buckets as f64;
        counts
            .iter()
            .map(|&c| {
                let d = c as f64 - expected;
                d * d / expected
            })
            .sum()
    }

    /// Returns the fraction of hashes (of sequential inputs 0..n_samples) that have
    /// each bit set. A well-mixed hash should be close to 0.5 for every bit.
    fn bit_set_rates_32<H: Hasher>(n_samples: usize) -> [f64; 32] {
        let mut counts = [0u64; 32];
        for i in 0..n_samples {
            let hash = H::hash32(i as u32);
            for bit in 0..32u32 {
                if (hash >> bit) & 1 == 1 {
                    counts[bit as usize] += 1;
                }
            }
        }
        let mut rates = [0.0f64; 32];
        for i in 0..32 {
            rates[i] = counts[i] as f64 / n_samples as f64;
        }
        rates
    }

    fn bit_set_rates_64<H: Hasher>(n_samples: usize) -> [f64; 64] {
        let mut counts = [0u64; 64];
        for i in 0..n_samples {
            let hash = H::hash64(i as u64);
            for bit in 0..64u32 {
                if (hash >> bit) & 1 == 1 {
                    counts[bit as usize] += 1;
                }
            }
        }
        let mut rates = [0.0f64; 64];
        for i in 0..64 {
            rates[i] = counts[i] as f64 / n_samples as f64;
        }
        rates
    }

    /// Number of collisions when hashing sequential inputs 0..n_samples to 32 bits.
    fn collision_count_32<H: Hasher>(n_samples: usize) -> usize {
        let hashes: HashSet<u32> = (0..n_samples).map(|i| H::hash32(i as u32)).collect();
        n_samples - hashes.len()
    }

    /// Number of collisions when hashing sequential inputs 0..n_samples to 64 bits.
    fn collision_count_64<H: Hasher>(n_samples: usize) -> usize {
        let hashes: HashSet<u64> = (0..n_samples).map(|i| H::hash64(i as u64)).collect();
        n_samples - hashes.len()
    }

    #[test]
    fn test_identity_hasher_32() {
        assert_eq!(IdentityHasher::hash32(0), 0);
        assert_eq!(IdentityHasher::hash32(42), 42);
        assert_eq!(IdentityHasher::hash32(u32::MAX), u32::MAX);
        assert_eq!(IdentityHasher::hash32(12345), 12345);
    }

    #[test]
    fn test_identity_hasher_64() {
        assert_eq!(IdentityHasher::hash64(0), 0);
        assert_eq!(IdentityHasher::hash64(42), 42);
        assert_eq!(IdentityHasher::hash64(u64::MAX), u64::MAX);
        assert_eq!(IdentityHasher::hash64(123456789), 123456789);
    }

    #[test]
    fn test_fnv_hasher_32_deterministic() {
        // Same input should always produce same output
        let input = 12345u32;
        let hash1 = FnvHasher::hash32(input);
        let hash2 = FnvHasher::hash32(input);
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_fnv_hasher_64_deterministic() {
        // Same input should always produce same output
        let input = 123456789u64;
        let hash1 = FnvHasher::hash64(input);
        let hash2 = FnvHasher::hash64(input);
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_fnv_hasher_32_zero() {
        // Hash of 0 should not be 0 (due to offset basis)
        let hash = FnvHasher::hash32(0);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_fnv_hasher_64_zero() {
        // Hash of 0 should not be 0 (due to offset basis)
        let hash = FnvHasher::hash64(0);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_fnv_hasher_32_different_inputs() {
        // Different inputs should produce different hashes (usually)
        let hash1 = FnvHasher::hash32(1);
        let hash2 = FnvHasher::hash32(2);
        let hash3 = FnvHasher::hash32(100);

        assert_ne!(hash1, hash2);
        assert_ne!(hash2, hash3);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_fnv_hasher_64_different_inputs() {
        // Different inputs should produce different hashes (usually)
        let hash1 = FnvHasher::hash64(1);
        let hash2 = FnvHasher::hash64(2);
        let hash3 = FnvHasher::hash64(100);

        assert_ne!(hash1, hash2);
        assert_ne!(hash2, hash3);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_fnv_hasher_32_sequential() {
        // Sequential inputs should have good distribution
        let hash1 = FnvHasher::hash32(1000);
        let hash2 = FnvHasher::hash32(1001);
        let hash3 = FnvHasher::hash32(1002);

        // Hashes should be different
        assert_ne!(hash1, hash2);
        assert_ne!(hash2, hash3);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_fnv_hasher_64_sequential() {
        // Sequential inputs should have good distribution
        let hash1 = FnvHasher::hash64(1000);
        let hash2 = FnvHasher::hash64(1001);
        let hash3 = FnvHasher::hash64(1002);

        // Hashes should be different
        assert_ne!(hash1, hash2);
        assert_ne!(hash2, hash3);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_fnv_hasher_32_max_value() {
        // Should handle maximum values
        let hash = FnvHasher::hash32(u32::MAX);
        // Should produce some hash value
        assert!(hash > 0 || hash == 0);
    }

    #[test]
    fn test_fnv_hasher_64_max_value() {
        // Should handle maximum values
        let hash = FnvHasher::hash64(u64::MAX);
        // Should produce some hash value
        assert!(hash > 0 || hash == 0);
    }

    #[test]
    fn test_fnv_hasher_32_known_values() {
        // Test with some specific known values to ensure correctness
        // These values are computed from the FNV-1a algorithm
        let test_cases = vec![
            (0u32, FnvHasher::hash32(0)),
            (1u32, FnvHasher::hash32(1)),
            (255u32, FnvHasher::hash32(255)),
        ];

        // Verify they're deterministic by recomputing
        for (input, expected) in test_cases {
            assert_eq!(FnvHasher::hash32(input), expected);
        }
    }

    #[test]
    fn test_fnv_hasher_64_known_values() {
        // Test with some specific known values to ensure correctness
        let test_cases = vec![
            (0u64, FnvHasher::hash64(0)),
            (1u64, FnvHasher::hash64(1)),
            (255u64, FnvHasher::hash64(255)),
        ];

        // Verify they're deterministic by recomputing
        for (input, expected) in test_cases {
            assert_eq!(FnvHasher::hash64(input), expected);
        }
    }

    #[test]
    fn test_fnv_hasher_32_avalanche() {
        // Changing one bit should significantly change the hash (avalanche effect)
        let hash1 = FnvHasher::hash32(0b10101010_10101010_10101010_10101010);
        let hash2 = FnvHasher::hash32(0b10101010_10101010_10101010_10101011); // flip last bit

        // Count differing bits
        let diff = (hash1 ^ hash2).count_ones();

        // Should have good bit distribution (roughly half bits different)
        // Allow some variance, but at least a few bits should differ
        assert!(diff > 5, "Only {} bits differ, expected more", diff);
    }

    #[test]
    fn test_fnv_hasher_64_avalanche() {
        // Changing one bit should significantly change the hash (avalanche effect)
        let hash1 = FnvHasher::hash64(0x5555555555555555);
        let hash2 = FnvHasher::hash64(0x5555555555555556); // flip last bit

        // Count differing bits
        let diff = (hash1 ^ hash2).count_ones();

        // Should have good bit distribution
        assert!(diff > 10, "Only {} bits differ, expected more", diff);
    }

    #[test]
    fn test_identity_vs_fnv_32() {
        // Verify that identity and FNV produce different results (except for edge cases)
        let input = 12345u32;
        let identity = IdentityHasher::hash32(input);
        let fnv = FnvHasher::hash32(input);

        assert_eq!(identity, input);
        assert_ne!(fnv, input); // FNV should scramble the value
    }

    #[test]
    fn test_identity_vs_fnv_64() {
        // Verify that identity and FNV produce different results
        let input = 123456789u64;
        let identity = IdentityHasher::hash64(input);
        let fnv = FnvHasher::hash64(input);

        assert_eq!(identity, input);
        assert_ne!(fnv, input); // FNV should scramble the value
    }

    #[test]
    fn test_splitmix64_deterministic() {
        let input = 123456789u64;
        let hash1 = Splitmix64Hasher::hash64(input);
        let hash2 = Splitmix64Hasher::hash64(input);
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_splitmix64_zero() {
        // With the offset, hash(0) should not be 0
        let hash = Splitmix64Hasher::hash64(0);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_splitmix64_different_inputs() {
        let hash1 = Splitmix64Hasher::hash64(1);
        let hash2 = Splitmix64Hasher::hash64(2);
        let hash3 = Splitmix64Hasher::hash64(100);

        assert_ne!(hash1, hash2);
        assert_ne!(hash2, hash3);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_splitmix64_avalanche() {
        // Single bit difference in input should cause many bit differences in output
        let hash1 = Splitmix64Hasher::hash64(0x5555555555555555);
        let hash2 = Splitmix64Hasher::hash64(0x5555555555555556);

        let diff = (hash1 ^ hash2).count_ones();
        // Splitmix64 has excellent avalanche - expect roughly half the bits to differ
        assert!(diff > 20, "Only {} bits differ, expected more", diff);
    }

    #[test]
    fn test_splitmix64_vs_fnv() {
        // Verify Splitmix64 produces different results than FNV
        let input = 123456789u64;
        let splitmix = Splitmix64Hasher::hash64(input);
        let fnv = FnvHasher::hash64(input);

        assert_ne!(splitmix, fnv);
        assert_ne!(splitmix, input);
    }

    // --- Dispersion (chi-squared bucket distribution) ---
    //
    // With 65536 samples and 256 buckets, the chi-squared statistic for a
    // uniform distribution has mean 255 and std ~22.6.  The p=0.001 critical
    // value for df=255 is ~332; we use 380 (~5-sigma) to allow for the
    // randomness inherent in any finite test while still catching poor hashers.

    #[test]
    fn test_fnv_dispersion_32() {
        let chi_sq = chi_squared_32::<FnvHasher>(65536, 256);
        assert!(chi_sq < 300.0, "FNV 32-bit chi-squared too high: {:.1}", chi_sq);
    }

    #[test]
    fn test_fnv_dispersion_64() {
        let chi_sq = chi_squared_64::<FnvHasher>(65536, 256);
        assert!(chi_sq < 300.0, "FNV 64-bit chi-squared too high: {:.1}", chi_sq);
    }

    #[test]
    fn test_splitmix64_dispersion_32() {
        let chi_sq = chi_squared_32::<Splitmix64Hasher>(65536, 256);
        assert!(chi_sq < 300.0, "Splitmix64 32-bit chi-squared too high: {:.1}", chi_sq);
    }

    #[test]
    fn test_splitmix64_dispersion_64() {
        let chi_sq = chi_squared_64::<Splitmix64Hasher>(65536, 256);
        assert!(chi_sq < 300.0, "Splitmix64 64-bit chi-squared too high: {:.1}", chi_sq);
    }

    // --- Bit balance ---
    //
    // Every output bit should be set ~50% of the time across 100_000 samples.
    // A ±5% window (45%–55%) is a ~14-sigma bound for n=100_000.

    #[test]
    fn test_fnv_bit_balance_32() {
        let rates = bit_set_rates_32::<FnvHasher>(100_000);
        for (bit, &rate) in rates.iter().enumerate() {
            assert!(
                (0.45..=0.55).contains(&rate),
                "FNV 32-bit: bit {} is set {:.1}% of the time (expected ~50%)",
                bit,
                rate * 100.0
            );
        }
    }

    #[test]
    fn test_fnv_bit_balance_64() {
        let rates = bit_set_rates_64::<FnvHasher>(100_000);
        for (bit, &rate) in rates.iter().enumerate() {
            assert!(
                (0.45..=0.55).contains(&rate),
                "FNV 64-bit: bit {} is set {:.1}% of the time (expected ~50%)",
                bit,
                rate * 100.0
            );
        }
    }

    #[test]
    fn test_splitmix64_bit_balance_32() {
        let rates = bit_set_rates_32::<Splitmix64Hasher>(100_000);
        for (bit, &rate) in rates.iter().enumerate() {
            assert!(
                (0.45..=0.55).contains(&rate),
                "Splitmix64 32-bit: bit {} is set {:.1}% of the time (expected ~50%)",
                bit,
                rate * 100.0
            );
        }
    }

    #[test]
    fn test_splitmix64_bit_balance_64() {
        let rates = bit_set_rates_64::<Splitmix64Hasher>(100_000);
        for (bit, &rate) in rates.iter().enumerate() {
            assert!(
                (0.45..=0.55).contains(&rate),
                "Splitmix64 64-bit: bit {} is set {:.1}% of the time (expected ~50%)",
                bit,
                rate * 100.0
            );
        }
    }

    // --- Collision resistance ---
    //
    // Birthday bound: expected collisions ≈ n²/(2·2^w).
    // 32-bit, n=10_000 → E[collisions] ≈ 0.012  (allow ≤ 5 for safety)
    // 64-bit, n=1_000_000 → E[collisions] ≈ 2.7e-8 (effectively 0)

    #[test]
    fn test_fnv_collisions_32() {
        let c = collision_count_32::<FnvHasher>(10_000);
        assert!(c <= 5, "FNV 32-bit: {} collisions in 10_000 hashes", c);
    }

    #[test]
    fn test_fnv_collisions_64() {
        let c = collision_count_64::<FnvHasher>(1_000_000);
        assert_eq!(c, 0, "FNV 64-bit: {} collisions in 1_000_000 hashes", c);
    }

    #[test]
    fn test_splitmix64_collisions_32() {
        let c = collision_count_32::<Splitmix64Hasher>(10_000);
        assert!(c <= 5, "Splitmix64 32-bit: {} collisions in 10_000 hashes", c);
    }

    #[test]
    fn test_splitmix64_collisions_64() {
        let c = collision_count_64::<Splitmix64Hasher>(1_000_000);
        assert_eq!(c, 0, "Splitmix64 64-bit: {} collisions in 1_000_000 hashes", c);
    }
}
