//! Multi-modal graph primitives.
//!
//! Nodes are typed by `Mode` (road junction, maritime waypoint / port, or
//! airport). Edges carry both a geographic length and a time cost; the final
//! LCP cost is a mode-specific function of both, plus any chokepoint premium.

use serde::{Deserialize, Serialize};

/// 0-indexed node handle. We cap at 4B which is well above the ~3M we expect.
pub type NodeId = u32;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u8)]
pub enum Mode {
    Road = 0,
    Maritime = 1,
    Air = 2,
    /// Intra-modal connector (e.g. port-to-road, airport-to-road).
    Transfer = 3,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Node {
    pub id: NodeId,
    pub lon: f32,
    pub lat: f32,
    pub mode: Mode,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Edge {
    pub from: NodeId,
    pub to: NodeId,
    pub length_km: f32,
    /// Time cost in hours, pre-chokepoint adjustment.
    pub time_hours: f32,
    pub mode: Mode,
}

/// Compressed-Sparse-Row adjacency list. Built once per (year, chokepoint
/// state) and held immutable during SSSP queries.
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph {
    pub nodes: Vec<Node>,
    /// Offsets into `edges`; node `u` has edges in `edges[offsets[u]..offsets[u+1]]`.
    pub offsets: Vec<u32>,
    pub edges: Vec<Edge>,
}

impl Graph {
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Iterate outgoing edges of a node.
    pub fn out_edges(&self, u: NodeId) -> &[Edge] {
        let start = self.offsets[u as usize] as usize;
        let end = self.offsets[u as usize + 1] as usize;
        &self.edges[start..end]
    }

    /// Build from an unsorted edge list. Edges must be pre-deduplicated.
    pub fn from_edges(nodes: Vec<Node>, mut edges: Vec<Edge>) -> Self {
        edges.sort_by_key(|e| e.from);
        let n = nodes.len();
        let mut offsets = vec![0u32; n + 1];
        for e in &edges {
            offsets[e.from as usize + 1] += 1;
        }
        for i in 1..=n {
            offsets[i] += offsets[i - 1];
        }
        Graph { nodes, offsets, edges }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiny_graph_roundtrip() {
        let nodes = vec![
            Node { id: 0, lon: 0.0, lat: 0.0, mode: Mode::Road },
            Node { id: 1, lon: 1.0, lat: 0.0, mode: Mode::Road },
            Node { id: 2, lon: 1.0, lat: 1.0, mode: Mode::Road },
        ];
        let edges = vec![
            Edge { from: 0, to: 1, length_km: 100.0, time_hours: 1.0, mode: Mode::Road },
            Edge { from: 1, to: 2, length_km: 100.0, time_hours: 1.0, mode: Mode::Road },
            Edge { from: 0, to: 2, length_km: 150.0, time_hours: 1.8, mode: Mode::Road },
        ];
        let g = Graph::from_edges(nodes, edges);
        assert_eq!(g.num_nodes(), 3);
        assert_eq!(g.num_edges(), 3);
        assert_eq!(g.out_edges(0).len(), 2);
        assert_eq!(g.out_edges(1).len(), 1);
        assert_eq!(g.out_edges(2).len(), 0);
    }
}
