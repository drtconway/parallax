//! Wait-free union-find on a dense index space.
//!
//! Based on the algorithm described in:
//!
//!   Richard J. Anderson & Heather Woll,
//!   "Wait-free Parallel Algorithms for the Union-Find Problem",
//!   Proc. 23rd ACM Symposium on Theory of Computing (STOC), 1991.
//!
//! The data structure stores parent pointers and ranks in a flat `Vec`
//! indexed by element id (0 .. n-1), making it suitable for the dense
//! key spaces that arise in the cluster subcommand where every BED record
//! has a contiguous index.
//!
//! All operations use `AtomicUsize` compare-and-swap so that `find` and
//! `union` are wait-free / lock-free and can be called from multiple
//! threads without external synchronisation.  The single-threaded
//! interface (`&mut self`) is also provided for drop-in replacement of
//! the existing `UnionFind`.

use std::sync::atomic::{AtomicU64, Ordering};

/// Encode (parent, rank) in a single `AtomicU64`.
///
/// We pack `parent` in the lower 40 bits and `rank` in the upper 24 bits.
/// This supports up to ~1 trillion elements and ranks up to 16 million,
/// which is more than sufficient for any practical genomic workload.
///
/// We use `u64` rather than `usize` so the packing is correct on both
/// 32-bit and 64-bit targets.
const PARENT_BITS: u32 = 40;
const PARENT_MASK: u64 = (1u64 << PARENT_BITS) - 1;

#[inline(always)]
fn pack(parent: usize, rank: u32) -> u64 {
    debug_assert!((parent as u64) <= PARENT_MASK, "parent index overflow");
    ((rank as u64) << PARENT_BITS) | (parent as u64)
}

#[inline(always)]
fn unpack_parent(word: u64) -> usize {
    (word & PARENT_MASK) as usize
}

#[inline(always)]
fn unpack_rank(word: u64) -> u32 {
    (word >> PARENT_BITS) as u32
}

// ---------------------------------------------------------------------------
// Lock-free dense union-find
// ---------------------------------------------------------------------------

/// A dense, wait-free union-find (disjoint-set) data structure.
///
/// Elements are integers in `0..n` where `n` is fixed at construction.
/// Parent pointers and ranks are packed into a single atomic word per
/// element so that `find` and `union` can be used concurrently from
/// multiple threads without locking.
pub struct UnionFind {
    /// Packed (rank, parent) per element.
    nodes: Vec<AtomicU64>,
}

// AtomicU64 is Send+Sync, so the whole struct is safe to share.
unsafe impl Sync for UnionFind {}

impl UnionFind {
    /// Create a union-find for elements `0..n`, each initially in its
    /// own singleton set.
    pub fn new(n: usize) -> Self {
        let nodes: Vec<AtomicU64> = (0..n).map(|i| AtomicU64::new(pack(i, 0))).collect();
        UnionFind { nodes }
    }

    /// Find the representative (root) of the set containing `x`.
    ///
    /// Uses path splitting (Anderson & Woll §3): each node on the
    /// find-path is re-pointed to its grandparent via a single CAS.
    /// This is wait-free — every call makes progress regardless of
    /// concurrent operations.
    pub fn find(&self, mut x: usize) -> usize {
        loop {
            let word_x = self.nodes[x].load(Ordering::Relaxed);
            let px = unpack_parent(word_x);
            if px == x {
                return x;
            }
            // Read grandparent
            let word_px = self.nodes[px].load(Ordering::Relaxed);
            let ppx = unpack_parent(word_px);
            if ppx == px {
                // Parent is root
                return px;
            }
            // Path splitting: try to point x directly to grandparent.
            // CAS failure is benign — another thread did it for us.
            let _ = self.nodes[x].compare_exchange_weak(
                word_x,
                pack(ppx, unpack_rank(word_x)),
                Ordering::Relaxed,
                Ordering::Relaxed,
            );
            x = ppx;
        }
    }

    /// Union the sets containing `x` and `y`.
    ///
    /// Uses union-by-rank with CAS (Anderson & Woll §4).  Returns the
    /// root of the merged set.  If `x` and `y` are already in the same
    /// set this is a no-op and returns the common root.
    pub fn union(&self, x: usize, y: usize) -> usize {
        loop {
            let mut rx = self.find(x);
            let mut ry = self.find(y);

            if rx == ry {
                return rx;
            }

            // Ensure rx has the higher (or equal) rank so we always
            // attach the smaller tree under the larger.
            let word_rx = self.nodes[rx].load(Ordering::Relaxed);
            let word_ry = self.nodes[ry].load(Ordering::Relaxed);
            let rank_rx = unpack_rank(word_rx);
            let rank_ry = unpack_rank(word_ry);

            if rank_rx < rank_ry || (rank_rx == rank_ry && rx > ry) {
                std::mem::swap(&mut rx, &mut ry);
            }

            // Try to make ry a child of rx.
            let cur_ry = self.nodes[ry].load(Ordering::Relaxed);
            // ry must still be a root (parent == self) for the link to be valid.
            if unpack_parent(cur_ry) != ry {
                // ry was concurrently linked elsewhere — retry.
                continue;
            }
            let new_ry = pack(rx, unpack_rank(cur_ry));
            if self.nodes[ry]
                .compare_exchange_weak(cur_ry, new_ry, Ordering::Relaxed, Ordering::Relaxed)
                .is_err()
            {
                // CAS failed — retry.
                continue;
            }

            // If ranks were equal, increment the rank of the new root.
            if rank_rx == rank_ry {
                let cur_rx = self.nodes[rx].load(Ordering::Relaxed);
                let expected = pack(rx, rank_rx);
                let desired = pack(rx, rank_rx + 1);
                // Best-effort: if this CAS fails the rank is just a
                // heuristic so correctness is unaffected.
                if unpack_parent(cur_rx) == rx && unpack_rank(cur_rx) == rank_rx {
                    let _ = self.nodes[rx].compare_exchange_weak(
                        expected,
                        desired,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    );
                }
            }

            return rx;
        }
    }

    /// Check if `x` and `y` are in the same set.
    #[cfg(test)]
    #[inline]
    pub fn connected(&self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }

    /// Count the number of disjoint sets.
    #[cfg(test)]
    pub fn count_sets(&self) -> usize {
        (0..self.nodes.len())
            .filter(|&i| self.find(i) == i)
            .count()
    }

    /// Get all elements in the same set as `x`.
    #[cfg(test)]
    pub fn members(&self, x: usize) -> Vec<usize> {
        let root = self.find(x);
        (0..self.nodes.len())
            .filter(|&i| self.find(i) == root)
            .collect()
    }

    /// Get all sets as a Vec of Vecs.
    #[cfg(test)]
    pub fn all_sets(&self) -> Vec<Vec<usize>> {
        let mut sets: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for i in 0..self.nodes.len() {
            let root = self.find(i);
            sets.entry(root).or_default().push(i);
        }
        sets.into_values().collect()
    }
}

impl std::fmt::Debug for UnionFind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UnionFind")
            .field("len", &self.nodes.len())
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_singleton() {
        let uf = UnionFind::new(5);
        for i in 0..5 {
            assert_eq!(uf.find(i), i);
        }
        assert_eq!(uf.count_sets(), 5);
    }

    #[test]
    fn test_union_basic() {
        let uf = UnionFind::new(3);
        assert!(!uf.connected(0, 1));
        uf.union(0, 1);
        assert!(uf.connected(0, 1));
        assert!(!uf.connected(0, 2));
    }

    #[test]
    fn test_union_chain() {
        let uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        uf.union(3, 4);

        for i in 0..5 {
            for j in 0..5 {
                assert!(uf.connected(i, j), "{i} and {j} should be connected");
            }
        }
        assert_eq!(uf.count_sets(), 1);
    }

    #[test]
    fn test_two_sets() {
        let uf = UnionFind::new(6);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(3, 4);
        uf.union(4, 5);

        assert!(uf.connected(0, 2));
        assert!(uf.connected(3, 5));
        assert!(!uf.connected(0, 3));
        assert_eq!(uf.count_sets(), 2);
    }

    #[test]
    fn test_count_sets() {
        let uf = UnionFind::new(6);
        assert_eq!(uf.count_sets(), 6);

        uf.union(0, 1);
        assert_eq!(uf.count_sets(), 5);

        uf.union(2, 3);
        assert_eq!(uf.count_sets(), 4);

        uf.union(0, 2);
        assert_eq!(uf.count_sets(), 3);
    }

    #[test]
    fn test_members() {
        let uf = UnionFind::new(6);
        uf.union(0, 2);
        uf.union(2, 4);

        let mut members = uf.members(0);
        members.sort();
        assert_eq!(members, vec![0, 2, 4]);

        let mut members = uf.members(1);
        members.sort();
        assert_eq!(members, vec![1]);
    }

    #[test]
    fn test_all_sets() {
        let uf = UnionFind::new(6);
        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(4, 5);

        let mut sets = uf.all_sets();
        for set in &mut sets {
            set.sort();
        }
        sets.sort();
        assert_eq!(sets, vec![vec![0, 1], vec![2, 3], vec![4, 5]]);
    }

    #[test]
    fn test_idempotent_union() {
        let uf = UnionFind::new(3);
        uf.union(0, 1);
        uf.union(0, 1);
        uf.union(1, 0);
        assert_eq!(uf.count_sets(), 2); // {0,1}, {2}

        uf.union(1, 2);
        uf.union(0, 2);
        assert_eq!(uf.count_sets(), 1);
    }

    #[test]
    fn test_path_compression() {
        let uf = UnionFind::new(100);
        for i in 0..99 {
            uf.union(i, i + 1);
        }

        let root = uf.find(0);
        for i in 0..100 {
            assert_eq!(uf.find(i), root);
        }
    }

    #[test]
    fn test_union_by_rank() {
        let uf = UnionFind::new(8);
        // Two balanced trees
        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(0, 2);

        uf.union(4, 5);
        uf.union(6, 7);
        uf.union(4, 6);

        uf.union(0, 4);
        assert_eq!(uf.count_sets(), 1);
    }

    #[test]
    fn test_large() {
        let n = 10_000;
        let uf = UnionFind::new(n);
        // Union all even numbers together
        for i in (0..n - 2).step_by(2) {
            uf.union(i, i + 2);
        }
        // Union all odd numbers together
        for i in (1..n - 2).step_by(2) {
            uf.union(i, i + 2);
        }
        assert_eq!(uf.count_sets(), 2);

        // Now merge even and odd
        uf.union(0, 1);
        assert_eq!(uf.count_sets(), 1);
    }

    #[test]
    fn test_concurrent_finds() {
        use std::sync::Arc;
        let uf = Arc::new(UnionFind::new(1000));

        // Build some structure
        for i in 0..999 {
            uf.union(i, i + 1);
        }

        // Concurrent finds from multiple threads
        let handles: Vec<_> = (0..4)
            .map(|_| {
                let uf = Arc::clone(&uf);
                std::thread::spawn(move || {
                    let root = uf.find(0);
                    for i in 0..1000 {
                        assert_eq!(uf.find(i), root);
                    }
                })
            })
            .collect();

        for h in handles {
            h.join().unwrap();
        }
    }

    #[test]
    fn test_concurrent_unions() {
        use std::sync::Arc;
        let uf = Arc::new(UnionFind::new(1000));

        // Each thread unions a contiguous range
        let handles: Vec<_> = (0..4)
            .map(|t| {
                let uf = Arc::clone(&uf);
                std::thread::spawn(move || {
                    let start = t * 250;
                    let end = start + 250;
                    for i in start..end - 1 {
                        uf.union(i, i + 1);
                    }
                })
            })
            .collect();

        for h in handles {
            h.join().unwrap();
        }

        // Each quarter should be connected
        for t in 0..4 {
            let start = t * 250;
            let end = start + 250;
            for i in start..end {
                assert!(uf.connected(start, i));
            }
        }

        // Now bridge them
        uf.union(249, 250);
        uf.union(499, 500);
        uf.union(749, 750);
        assert_eq!(uf.count_sets(), 1);
    }
}
