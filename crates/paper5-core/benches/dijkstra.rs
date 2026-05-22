//! Benchmark: SSSP on a synthetic grid graph of ~100K nodes, ~400K edges.
//! Run with `cargo bench -p paper5-core`.
//!
//! Target: < 20 ms per SSSP on an M-series MBP at v0.1 (plain Dijkstra).
//! If we beat 5 ms post-contraction, we're on pace for the full panel.

use criterion::{criterion_group, criterion_main, Criterion};
use paper5_core::{
    chokepoints::ChokepointState,
    dijkstra::sssp,
    graph::{Edge, Graph, Mode, Node},
};

fn synthetic_grid(side: u32) -> Graph {
    let n = side * side;
    let nodes: Vec<Node> = (0..n)
        .map(|i| Node {
            id: i,
            lon: (i % side) as f32,
            lat: (i / side) as f32,
            mode: Mode::Road,
        })
        .collect();

    let mut edges = Vec::with_capacity((n * 4) as usize);
    for y in 0..side {
        for x in 0..side {
            let u = y * side + x;
            if x + 1 < side {
                let v = y * side + (x + 1);
                edges.push(Edge { from: u, to: v, length_km: 1.0, time_hours: 0.02, mode: Mode::Road });
                edges.push(Edge { from: v, to: u, length_km: 1.0, time_hours: 0.02, mode: Mode::Road });
            }
            if y + 1 < side {
                let v = (y + 1) * side + x;
                edges.push(Edge { from: u, to: v, length_km: 1.0, time_hours: 0.02, mode: Mode::Road });
                edges.push(Edge { from: v, to: u, length_km: 1.0, time_hours: 0.02, mode: Mode::Road });
            }
        }
    }
    Graph::from_edges(nodes, edges)
}

fn bench_sssp(c: &mut Criterion) {
    let g = synthetic_grid(316); // ≈ 100K nodes
    let cp = ChokepointState::default();
    c.bench_function("sssp_100k", |b| b.iter(|| sssp(&g, 0, &cp)));
}

criterion_group!(benches, bench_sssp);
criterion_main!(benches);
