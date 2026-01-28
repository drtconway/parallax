#![allow(dead_code)]

use std::collections::HashMap;

/// A Union-Find (Disjoint Set Union) data structure.
///
/// Supports efficient union and find operations with:
/// - Path compression on find
/// - Union by rank
/// - Lazy initialization (elements are added on first use)
///
/// Time complexity: O(α(n)) amortized per operation, where α is the
/// inverse Ackermann function (effectively constant for all practical purposes).
#[derive(Debug, Clone, Default)]
pub struct UnionFind {
    /// Parent pointers. parent[i] == i means i is a root.
    /// Missing entries are implicitly self-parented.
    parent: HashMap<usize, usize>,
    /// Rank (upper bound on tree height) for union by rank.
    rank: HashMap<usize, usize>,
}

impl UnionFind {
    /// Create a new empty UnionFind.
    pub fn new() -> Self {
        UnionFind {
            parent: HashMap::new(),
            rank: HashMap::new(),
        }
    }

    /// Ensure an element exists in the structure.
    /// If it doesn't exist, it becomes its own set.
    fn ensure_exists(&mut self, x: usize) {
        if !self.parent.contains_key(&x) {
            self.parent.insert(x, x);
            self.rank.insert(x, 0);
        }
    }

    /// Find the representative (root) of the set containing `x`.
    ///
    /// Uses path compression: all nodes on the path to root are updated
    /// to point directly to the root.
    ///
    /// If `x` hasn't been seen before, it becomes its own set.
    pub fn find(&mut self, x: usize) -> usize {
        self.ensure_exists(x);
        
        let parent = self.parent[&x];
        if parent != x {
            let root = self.find(parent);
            self.parent.insert(x, root);
            root
        } else {
            x
        }
    }

    /// Union the sets containing `x` and `y`.
    ///
    /// Uses union by rank: the smaller tree is attached under the larger tree.
    /// Returns `true` if `x` and `y` were in different sets (union performed),
    /// `false` if they were already in the same set.
    pub fn union(&mut self, x: usize, y: usize) -> bool {
        let rx = self.find(x);
        let ry = self.find(y);

        if rx == ry {
            return false;
        }

        // Union by rank
        let rank_rx = self.rank[&rx];
        let rank_ry = self.rank[&ry];
        
        match rank_rx.cmp(&rank_ry) {
            std::cmp::Ordering::Less => {
                self.parent.insert(rx, ry);
            }
            std::cmp::Ordering::Greater => {
                self.parent.insert(ry, rx);
            }
            std::cmp::Ordering::Equal => {
                self.parent.insert(ry, rx);
                self.rank.insert(rx, rank_rx + 1);
            }
        }

        true
    }

    /// Check if `x` and `y` are in the same set.
    pub fn connected(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }

    /// Return the number of elements that have been added.
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    /// Check if the structure is empty.
    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// Count the number of disjoint sets.
    pub fn count_sets(&mut self) -> usize {
        let keys: Vec<_> = self.parent.keys().copied().collect();
        keys.iter()
            .filter(|&&i| self.find(i) == i)
            .count()
    }

    /// Get all elements in the same set as `x`.
    pub fn members(&mut self, x: usize) -> Vec<usize> {
        let root = self.find(x);
        let keys: Vec<_> = self.parent.keys().copied().collect();
        keys.into_iter()
            .filter(|&i| self.find(i) == root)
            .collect()
    }

    /// Get all sets as a Vec of Vecs.
    pub fn all_sets(&mut self) -> Vec<Vec<usize>> {
        let keys: Vec<_> = self.parent.keys().copied().collect();
        let mut sets: HashMap<usize, Vec<usize>> = HashMap::new();

        for i in keys {
            let root = self.find(i);
            sets.entry(root).or_default().push(i);
        }

        sets.into_values().collect()
    }
    
    /// Check if an element has been added to the structure.
    pub fn contains(&self, x: usize) -> bool {
        self.parent.contains_key(&x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let uf = UnionFind::new();
        assert_eq!(uf.len(), 0);
        assert!(uf.is_empty());
    }

    #[test]
    fn test_lazy_init() {
        let mut uf = UnionFind::new();
        assert!(!uf.contains(42));
        
        // Finding an element adds it
        assert_eq!(uf.find(42), 42);
        assert!(uf.contains(42));
        assert_eq!(uf.len(), 1);
    }

    #[test]
    fn test_find_self() {
        let mut uf = UnionFind::new();
        for i in [0, 5, 100, 1000] {
            assert_eq!(uf.find(i), i);
        }
        assert_eq!(uf.len(), 4);
    }

    #[test]
    fn test_union_basic() {
        let mut uf = UnionFind::new();

        // Initially all separate
        assert!(!uf.connected(0, 1));
        assert!(!uf.connected(1, 2));

        // Union 0 and 1
        assert!(uf.union(0, 1)); // Returns true, they were separate
        assert!(uf.connected(0, 1));
        assert!(!uf.connected(0, 2));

        // Union again should return false
        assert!(!uf.union(0, 1));
        assert!(!uf.union(1, 0));
    }

    #[test]
    fn test_union_chain() {
        let mut uf = UnionFind::new();

        // Create chain: 0-1-2-3-4
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        uf.union(3, 4);

        // All should be connected
        for i in 0..5 {
            for j in 0..5 {
                assert!(uf.connected(i, j), "{} and {} should be connected", i, j);
            }
        }
    }

    #[test]
    fn test_two_sets() {
        let mut uf = UnionFind::new();

        // Set 1: {0, 1, 2}
        uf.union(0, 1);
        uf.union(1, 2);

        // Set 2: {3, 4, 5}
        uf.union(3, 4);
        uf.union(4, 5);

        // Within sets
        assert!(uf.connected(0, 2));
        assert!(uf.connected(3, 5));

        // Between sets
        assert!(!uf.connected(0, 3));
        assert!(!uf.connected(2, 4));
    }

    #[test]
    fn test_count_sets() {
        let mut uf = UnionFind::new();
        assert_eq!(uf.count_sets(), 0);

        uf.find(0); // Add element 0
        assert_eq!(uf.count_sets(), 1);

        uf.find(1); // Add element 1
        assert_eq!(uf.count_sets(), 2);

        uf.union(0, 1);
        assert_eq!(uf.count_sets(), 1);

        uf.union(2, 3);
        assert_eq!(uf.count_sets(), 2);

        uf.union(0, 2); // Merge two sets
        assert_eq!(uf.count_sets(), 1);
    }

    #[test]
    fn test_members() {
        let mut uf = UnionFind::new();

        uf.union(0, 2);
        uf.union(2, 4);

        let mut members = uf.members(0);
        members.sort();
        assert_eq!(members, vec![0, 2, 4]);

        uf.find(1); // Add singleton
        let mut members = uf.members(1);
        members.sort();
        assert_eq!(members, vec![1]);
    }

    #[test]
    fn test_all_sets() {
        let mut uf = UnionFind::new();

        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(4, 5);

        let mut sets = uf.all_sets();
        // Sort each set and then sort the sets for deterministic comparison
        for set in &mut sets {
            set.sort();
        }
        sets.sort();

        assert_eq!(sets, vec![vec![0, 1], vec![2, 3], vec![4, 5]]);
    }

    #[test]
    fn test_sparse_elements() {
        let mut uf = UnionFind::new();

        // Use sparse, non-contiguous elements
        uf.union(100, 200);
        uf.union(200, 300);
        uf.union(1000, 2000);

        assert!(uf.connected(100, 300));
        assert!(uf.connected(1000, 2000));
        assert!(!uf.connected(100, 1000));

        assert_eq!(uf.len(), 5);
        assert_eq!(uf.count_sets(), 2);
    }

    #[test]
    fn test_path_compression() {
        let mut uf = UnionFind::new();

        // Create a long chain
        for i in 0..99 {
            uf.union(i, i + 1);
        }

        // After find with path compression, all should have same root
        let root = uf.find(0);
        for i in 0..100 {
            assert_eq!(uf.find(i), root);
        }

        // Verify path compression happened (parents should point directly to root)
        for i in 0..100 {
            assert_eq!(uf.parent[&i], root);
        }
    }

    #[test]
    fn test_union_by_rank() {
        let mut uf = UnionFind::new();

        // Build two balanced trees
        // Tree 1: union(0,1), union(2,3), union(0,2) -> root with rank 2
        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(0, 2);

        // Tree 2: union(4,5), union(6,7), union(4,6) -> root with rank 2
        uf.union(4, 5);
        uf.union(6, 7);
        uf.union(4, 6);

        // Now union the two trees
        uf.union(0, 4);

        // All should be in one set
        assert_eq!(uf.count_sets(), 1);
    }

    #[test]
    fn test_idempotent_union() {
        let mut uf = UnionFind::new();

        assert!(uf.union(0, 1));
        assert!(!uf.union(0, 1));
        assert!(!uf.union(1, 0));

        assert!(uf.union(1, 2));
        assert!(!uf.union(0, 2)); // 0 and 2 now connected through 1
    }
    
    #[test]
    fn test_default() {
        let uf: UnionFind = Default::default();
        assert!(uf.is_empty());
    }
}

