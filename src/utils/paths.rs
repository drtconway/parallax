//! Shortest path algorithms using Fibonacci heaps.
//!
//! This module provides an iterator-based interface for finding all shortest paths
//! in directed acyclic graphs (DAGs) with non-negative edge weights.
//!
//! # Example
//! ```ignore
//! use parallax::utils::paths::ShortestPaths;
//!
//! // Edges: (from, to, weight)
//! let edges = vec![
//!     (0, 1, 1.0),
//!     (0, 2, 4.0),
//!     (1, 2, 2.0),
//!     (1, 3, 5.0),
//!     (2, 3, 1.0),
//! ];
//!
//! for (target, distance, path) in ShortestPaths::new(4, 0, &edges) {
//!     println!("Path to {}: distance={}, path={:?}", target, distance, path);
//! }
//! ```

#![allow(dead_code)]

use ordered_float::OrderedFloat;

use super::fibomacci::FibHeap;
use super::heap::{HeapOrdering, Heapable};
use std::collections::VecDeque;

/// A weighted directed edge.
#[derive(Clone, Debug)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub weight: f64,
}

impl Edge {
    pub fn new(from: usize, to: usize, weight: f64) -> Self {
        Self { from, to, weight }
    }
}

impl From<(usize, usize, f64)> for Edge {
    fn from(tuple: (usize, usize, f64)) -> Self {
        Edge::new(tuple.0, tuple.1, tuple.2)
    }
}

/// State for Dijkstra's priority queue - min-heap by distance.
#[derive(Clone, Debug)]
struct DijkstraState {
    distance: f64,
    node: usize,
}

/// Min-heap configuration for Dijkstra's algorithm.
struct DijkstraHeap;

impl Heapable for DijkstraHeap {
    type Item = DijkstraState;
    const ORDERING: HeapOrdering = HeapOrdering::Min;

    fn cmp(lhs: &Self::Item, rhs: &Self::Item) -> std::cmp::Ordering {
        OrderedFloat(lhs.distance)
             .cmp(&OrderedFloat(rhs.distance))
    }
}

/// Result of shortest path computation for a single target node.
#[derive(Clone, Debug)]
pub struct ShortestPathResult {
    /// The target node.
    pub target: usize,
    /// The shortest distance from source to target.
    pub distance: f64,
    /// The path from source to target (inclusive).
    pub path: Vec<usize>,
}

/// Iterator over all shortest paths from a source node.
///
/// For DAGs with non-negative weights, this uses Dijkstra's algorithm with
/// a Fibonacci heap. When there are multiple shortest paths to a node
/// (i.e., ties), all such paths are yielded.
pub struct ShortestPaths {
    /// Shortest distance to each node (f64::INFINITY if unreachable).
    distances: Vec<f64>,
    /// All predecessors on shortest paths to each node.
    /// Multiple predecessors means multiple shortest paths.
    predecessors: Vec<Vec<usize>>,
    /// Source node.
    source: usize,
    /// Queue of (target, path_stack) for BFS over paths.
    /// path_stack is built backwards from target to source.
    path_queue: VecDeque<(usize, Vec<usize>)>,
    /// Targets we still need to enumerate paths for.
    pending_targets: VecDeque<usize>,
}

impl ShortestPaths {
    /// Create a new shortest paths iterator.
    ///
    /// # Arguments
    /// * `num_nodes` - Total number of nodes in the graph (0..num_nodes).
    /// * `source` - The source node to find paths from.
    /// * `edges` - Slice of (from, to, weight) tuples representing directed edges.
    ///
    /// # Panics
    /// Panics if any edge has a negative weight.
    pub fn new(num_nodes: usize, source: usize, edges: &[(usize, usize, f64)]) -> Self {
        // Build adjacency list
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); num_nodes];
        for &(from, to, weight) in edges {
            assert!(weight >= 0.0, "Edge weights must be non-negative");
            adj[from].push((to, weight));
        }

        Self::from_adjacency(num_nodes, source, &adj)
    }

    /// Create from an adjacency list representation.
    ///
    /// # Arguments
    /// * `num_nodes` - Total number of nodes.
    /// * `source` - The source node.
    /// * `adj` - Adjacency list where adj[u] contains (v, weight) for edge u->v.
    pub fn from_adjacency(num_nodes: usize, source: usize, adj: &[Vec<(usize, f64)>]) -> Self {
        let (distances, predecessors) = Self::dijkstra(num_nodes, source, adj);

        // Initialize pending targets (all reachable nodes except source)
        let pending_targets: VecDeque<usize> = (0..num_nodes)
            .filter(|&n| n != source && distances[n] < f64::INFINITY)
            .collect();

        ShortestPaths {
            distances,
            predecessors,
            source,
            path_queue: VecDeque::new(),
            pending_targets,
        }
    }

    /// Run Dijkstra's algorithm using Fibonacci heap.
    /// Returns (distances, predecessors) where predecessors[v] contains all
    /// nodes u such that the shortest path to v goes through u.
    fn dijkstra(
        num_nodes: usize,
        source: usize,
        adj: &[Vec<(usize, f64)>],
    ) -> (Vec<f64>, Vec<Vec<usize>>) {
        let mut distances = vec![f64::INFINITY; num_nodes];
        let mut predecessors: Vec<Vec<usize>> = vec![Vec::new(); num_nodes];
        let mut visited = vec![false; num_nodes];

        distances[source] = 0.0;

        let mut heap: FibHeap<DijkstraHeap> = FibHeap::new();
        heap.push(DijkstraState {
            distance: 0.0,
            node: source,
        });

        while let Some(state) = heap.pop() {
            let u = state.node;

            // Skip if already processed (we may have stale entries)
            if visited[u] {
                continue;
            }
            visited[u] = true;

            // Relax all outgoing edges
            for &(v, weight) in &adj[u] {
                let new_dist = distances[u] + weight;

                if new_dist < distances[v] {
                    // Found a strictly shorter path
                    distances[v] = new_dist;
                    predecessors[v].clear();
                    predecessors[v].push(u);
                    heap.push(DijkstraState {
                        distance: new_dist,
                        node: v,
                    });
                } else if (new_dist - distances[v]).abs() < 1e-10 && !predecessors[v].contains(&u) {
                    // Found an equally short path (tie)
                    predecessors[v].push(u);
                }
            }
        }

        (distances, predecessors)
    }

    /// Get the shortest distance to a node, or None if unreachable.
    pub fn distance_to(&self, target: usize) -> Option<f64> {
        if target < self.distances.len() && self.distances[target] < f64::INFINITY {
            Some(self.distances[target])
        } else {
            None
        }
    }

    /// Get all shortest distances (indexed by node).
    pub fn distances(&self) -> &[f64] {
        &self.distances
    }

    /// Check if a node is reachable from the source.
    pub fn is_reachable(&self, target: usize) -> bool {
        target < self.distances.len() && self.distances[target] < f64::INFINITY
    }

    /// Enumerate all paths to a specific target.
    /// Returns an iterator over all shortest paths to that target.
    pub fn paths_to(&self, target: usize) -> impl Iterator<Item = Vec<usize>> + '_ {
        PathsToTarget::new(self, target)
    }
}

impl Iterator for ShortestPaths {
    type Item = ShortestPathResult;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // First, try to complete a path we're currently building
            if let Some((node, mut path_so_far)) = self.path_queue.pop_front() {
                if node == self.source {
                    // We've reached the source - reverse to get source->target path
                    path_so_far.reverse();
                    let target = *path_so_far.last().unwrap();
                    return Some(ShortestPathResult {
                        target,
                        distance: self.distances[target],
                        path: path_so_far,
                    });
                } else {
                    // Expand this partial path with all predecessors
                    // Add all predecessor expansions to the queue
                    for (i, &pred) in self.predecessors[node].iter().enumerate() {
                        let mut new_path = if i == self.predecessors[node].len() - 1 {
                            // Reuse the vector for the last predecessor
                            path_so_far.clone()
                        } else {
                            path_so_far.clone()
                        };
                        new_path.push(pred);
                        self.path_queue.push_back((pred, new_path));
                    }
                    // Use the original path_so_far for the last one to avoid extra clone
                    if !self.predecessors[node].is_empty() {
                        continue; // Process the queue entries we just added
                    }
                }
            }

            // No partial paths - start a new target
            if let Some(target) = self.pending_targets.pop_front() {
                // Start building paths to this target
                self.path_queue.push_back((target, vec![target]));
                continue;
            }

            // No more targets
            return None;
        }
    }
}

/// Iterator over all shortest paths to a specific target.
struct PathsToTarget<'a> {
    paths: &'a ShortestPaths,
    target: usize,
    /// Stack of (current_node, path_so_far, predecessor_index)
    stack: Vec<(usize, Vec<usize>, usize)>,
}

impl<'a> PathsToTarget<'a> {
    fn new(paths: &'a ShortestPaths, target: usize) -> Self {
        let mut iter = PathsToTarget {
            paths,
            target,
            stack: Vec::new(),
        };

        // Initialize if target is reachable
        if paths.is_reachable(target) || target == paths.source {
            iter.stack.push((target, vec![target], 0));
        }

        iter
    }
}

impl<'a> Iterator for PathsToTarget<'a> {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((node, path, pred_idx)) = self.stack.pop() {
            if node == self.paths.source {
                // Complete path found
                let mut result = path;
                result.reverse();
                return Some(result);
            }

            let preds = &self.paths.predecessors[node];
            if pred_idx < preds.len() {
                // Push state for next predecessor at this node
                if pred_idx + 1 < preds.len() {
                    self.stack.push((node, path.clone(), pred_idx + 1));
                }

                // Explore this predecessor
                let pred = preds[pred_idx];
                let mut new_path = path;
                new_path.push(pred);
                self.stack.push((pred, new_path, 0));
            }
        }

        None
    }
}

/// Convenience function to find all shortest paths in a weighted DAG.
pub fn all_shortest_paths(
    num_nodes: usize,
    source: usize,
    edges: &[(usize, usize, f64)],
) -> ShortestPaths {
    ShortestPaths::new(num_nodes, source, edges)
}

/// Find the single shortest path to a target (or None if unreachable).
/// If multiple shortest paths exist, returns one of them.
pub fn shortest_path(
    num_nodes: usize,
    source: usize,
    target: usize,
    edges: &[(usize, usize, f64)],
) -> Option<(f64, Vec<usize>)> {
    let paths = ShortestPaths::new(num_nodes, source, edges);
    paths.paths_to(target).next().map(|p| {
        let dist = paths.distances[target];
        (dist, p)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_path() {
        // Linear graph: 0 -> 1 -> 2 -> 3
        let edges = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)];

        let paths = ShortestPaths::new(4, 0, &edges);
        assert_eq!(paths.distance_to(3), Some(3.0));

        let all: Vec<_> = paths.paths_to(3).collect();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0], vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_multiple_paths_same_length() {
        // Diamond graph with equal weights:
        //     1
        //    / \
        //   0   3
        //    \ /
        //     2
        let edges = vec![
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 3, 1.0),
            (2, 3, 1.0),
        ];

        let paths = ShortestPaths::new(4, 0, &edges);
        assert_eq!(paths.distance_to(3), Some(2.0));

        let mut all: Vec<_> = paths.paths_to(3).collect();
        all.sort();
        assert_eq!(all.len(), 2);
        assert!(all.contains(&vec![0, 1, 3]));
        assert!(all.contains(&vec![0, 2, 3]));
    }

    #[test]
    fn test_single_shorter_path() {
        // Diamond graph with unequal weights:
        //     1 (weight 1)
        //    / \
        //   0   3
        //    \ /
        //     2 (weight 5)
        let edges = vec![
            (0, 1, 1.0),
            (0, 2, 5.0),
            (1, 3, 1.0),
            (2, 3, 1.0),
        ];

        let paths = ShortestPaths::new(4, 0, &edges);
        assert_eq!(paths.distance_to(3), Some(2.0));

        let all: Vec<_> = paths.paths_to(3).collect();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0], vec![0, 1, 3]);
    }

    #[test]
    fn test_iterator_all_targets() {
        let edges = vec![
            (0, 1, 1.0),
            (0, 2, 2.0),
            (1, 3, 1.0),
            (2, 3, 1.0),
        ];

        let results: Vec<_> = ShortestPaths::new(4, 0, &edges).collect();

        // Should have paths to nodes 1, 2, 3 (not source 0)
        // 0->1 distance 1
        // 0->2 distance 2
        // 0->1->3 distance 2 (shorter than 0->2->3 = 3)

        // Node 1: one path
        assert!(results.iter().any(|r| r.target == 1 && r.path == vec![0, 1]));

        // Node 2: one path
        assert!(results.iter().any(|r| r.target == 2 && r.path == vec![0, 2]));

        // Node 3: one path (0->1->3 is shorter than 0->2->3)
        let paths_to_3: Vec<_> = results.iter().filter(|r| r.target == 3).collect();
        assert_eq!(paths_to_3.len(), 1);
        assert_eq!(paths_to_3[0].path, vec![0, 1, 3]);
        assert_eq!(paths_to_3[0].distance, 2.0);
    }

    #[test]
    fn test_unreachable_node() {
        let edges = vec![(0, 1, 1.0)];

        let paths = ShortestPaths::new(3, 0, &edges);
        assert_eq!(paths.distance_to(1), Some(1.0));
        assert_eq!(paths.distance_to(2), None);
        assert!(!paths.is_reachable(2));

        let all: Vec<_> = paths.paths_to(2).collect();
        assert!(all.is_empty());
    }

    #[test]
    fn test_source_to_source() {
        let edges = vec![(0, 1, 1.0)];

        let paths = ShortestPaths::new(2, 0, &edges);
        assert_eq!(paths.distance_to(0), Some(0.0));

        let all: Vec<_> = paths.paths_to(0).collect();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0], vec![0]);
    }

    #[test]
    fn test_complex_dag() {
        //        1
        //       /|\
        //      / | \
        //     0--2--4
        //      \ | /
        //       \|/
        //        3
        let edges = vec![
            (0, 1, 1.0),
            (0, 2, 1.0),
            (0, 3, 1.0),
            (1, 4, 1.0),
            (2, 4, 1.0),
            (3, 4, 1.0),
        ];

        let paths = ShortestPaths::new(5, 0, &edges);
        assert_eq!(paths.distance_to(4), Some(2.0));

        // Three equal-length paths to 4
        let mut all: Vec<_> = paths.paths_to(4).collect();
        all.sort();
        assert_eq!(all.len(), 3);
        assert!(all.contains(&vec![0, 1, 4]));
        assert!(all.contains(&vec![0, 2, 4]));
        assert!(all.contains(&vec![0, 3, 4]));
    }

    #[test]
    fn test_convenience_functions() {
        let edges = vec![(0, 1, 1.0), (1, 2, 2.0)];

        let result = shortest_path(3, 0, 2, &edges);
        assert!(result.is_some());
        let (dist, path) = result.unwrap();
        assert_eq!(dist, 3.0);
        assert_eq!(path, vec![0, 1, 2]);

        // Unreachable
        let result = shortest_path(3, 0, 2, &[]);
        assert!(result.is_none());
    }

    #[test]
    fn test_zero_weight_edges() {
        let edges = vec![(0, 1, 0.0), (1, 2, 0.0)];

        let paths = ShortestPaths::new(3, 0, &edges);
        assert_eq!(paths.distance_to(2), Some(0.0));

        let all: Vec<_> = paths.paths_to(2).collect();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0], vec![0, 1, 2]);
    }

    #[test]
    #[should_panic(expected = "non-negative")]
    fn test_negative_weight_panics() {
        let edges = vec![(0, 1, -1.0)];
        let _ = ShortestPaths::new(2, 0, &edges);
    }
}
