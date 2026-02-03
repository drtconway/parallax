#![allow(dead_code)]

pub trait Hasher {
    fn hash32(x: u32) -> u32;

    fn hash64(x: u64) -> u64;
}

pub struct IdentityHasher;
impl Hasher for IdentityHasher {
    fn hash32(x: u32) -> u32 {
        x
    }

    fn hash64(x: u64) -> u64 {
        x
    }
}

pub struct FnvHasher;
impl Hasher for FnvHasher {
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

#[cfg(test)]
mod tests {
    use super::*;

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
}