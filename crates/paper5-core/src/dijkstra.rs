//! Dijkstra + many-to-many SSSP over the multi-modal graph.
//!
//! v0.1: plain binary-heap Dijkstra. Good enough for the 20-country prototype.
//! v0.2 (sprint day 5): bidirectional Dijkstra + node-contraction preprocessing
//! for the full planet graph.

use crate::chokepoints::ChokepointState;
use crate::graph::{Graph, Mode, NodeId};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Mode-specific freight rate ($/hour) — placeholder values calibrated from
/// Hummels-Schaur (2013) time-value-of-trade and standard freight indices.
/// Replace with empirically-sourced values in L3 before v0.9.
pub const COST_PER_HOUR: [f32; 4] = [
    30.0,   // Road (truck freight)
    8.0,    // Maritime (container)
    220.0,  // Air (bellyhold + dedicated freighter blend)
    5.0,    // Transfer (port/airport handling)
];

#[inline]
fn edge_cost(edge_mode: Mode, time_hours: f32, cp: &ChokepointState, flags: u8) -> f32 {
    let base = COST_PER_HOUR[edge_mode as usize] * time_hours;
    base * cp.multiplier(flags)
}

#[derive(Copy, Clone, PartialEq)]
struct HeapEntry {
    cost: f32,
    node: NodeId,
}

impl Eq for HeapEntry {}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Single-source shortest-path cost from `source` to every node reachable
/// under current chokepoint state.
pub fn sssp(graph: &Graph, source: NodeId, cp: &ChokepointState) -> Vec<f32> {
    let n = graph.num_nodes();
    let mut dist = vec![f32::INFINITY; n];
    dist[source as usize] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(HeapEntry { cost: 0.0, node: source });

    while let Some(HeapEntry { cost, node }) = heap.pop() {
        if cost > dist[node as usize] {
            continue;
        }
        for e in graph.out_edges(node) {
            // Chokepoint flags would come from edge metadata in v0.2;
            // for now we treat `0` as "no chokepoint" — wire up in L2.
            let w = edge_cost(e.mode, e.time_hours, cp, 0);
            let nd = cost + w;
            if nd < dist[e.to as usize] {
                dist[e.to as usize] = nd;
                heap.push(HeapEntry { cost: nd, node: e.to });
            }
        }
    }
    dist
}

/// Shortest path cost between a single (origin, destination) pair.
pub fn solve_pair(graph: &Graph, src: NodeId, dst: NodeId, cp: &ChokepointState) -> f32 {
    let dist = sssp(graph, src, cp);
    dist[dst as usize]
}

/// Many-to-many: return a |sources| × |targets| cost matrix in row-major layout.
/// Parallelized across sources via rayon.
pub fn solve_many_to_many(
    graph: &Graph,
    sources: &[NodeId],
    targets: &[NodeId],
    cp: &ChokepointState,
) -> Vec<f32> {
    use rayon::prelude::*;

    let ns = sources.len();
    let nt = targets.len();
    let mut out = vec![f32::INFINITY; ns * nt];

    let rows: Vec<Vec<f32>> = sources
        .par_iter()
        .map(|&s| {
            let dist = sssp(graph, s, cp);
            targets.iter().map(|&t| dist[t as usize]).collect()
        })
        .collect();

    for (i, row) in rows.into_iter().enumerate() {
        out[i * nt..(i + 1) * nt].copy_from_slice(&row);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Edge, Graph, Mode, Node};

    fn triangle() -> Graph {
        let nodes = vec![
            Node { id: 0, lon: 0.0, lat: 0.0, mode: Mode::Road },
            Node { id: 1, lon: 1.0, lat: 0.0, mode: Mode::Road },
            Node { id: 2, lon: 1.0, lat: 1.0, mode: Mode::Road },
        ];
        let edges = vec![
            Edge { from: 0, to: 1, length_km: 100.0, time_hours: 1.0, mode: Mode::Road },
            Edge { from: 1, to: 2, length_km: 100.0, time_hours: 1.0, mode: Mode::Road },
            Edge { from: 0, to: 2, length_km: 300.0, time_hours: 3.0, mode: Mode::Road },
        ];
        Graph::from_edges(nodes, edges)
    }

    #[test]
    fn takes_cheaper_two_hop_over_direct() {
        let g = triangle();
        let cp = ChokepointState::default();
        // 0->2 direct = 3h * 30 = 90;  0->1->2 = 2h * 30 = 60. Should pick 60.
        let c = solve_pair(&g, 0, 2, &cp);
        assert!((c - 60.0).abs() < 1e-4, "got {c}");
    }
}
