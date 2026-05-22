//! paper5-core — multi-modal LCP solver for the Effective Distance panel.
//!
//! See `ARCHITECTURE.md` in the repo root for the full design rationale.
//! This crate is the L2 layer: it takes a graph of road/maritime/air arcs
//! and a set of agglomeration centroids, and returns pair-wise least-cost-path
//! costs that the Python layer assembles into the eight distance variants.
//!
//! Scaffold as of v0.1 (sprint day 1).

pub mod graph;
pub mod dijkstra;
pub mod chokepoints;

pub use graph::{Graph, Edge, Mode, NodeId};
pub use dijkstra::{solve_pair, solve_many_to_many};
pub use chokepoints::ChokepointState;

/// Library version. Stamped into Parquet metadata so the reader can verify
/// which solver emitted a given output file.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
