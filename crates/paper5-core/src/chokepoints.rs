//! Chokepoint state machine for maritime arcs.
//!
//! Edges carry an 8-bit bitmask flagging which chokepoints they traverse.
//! At query time we look up the current state and apply the relevant
//! capacity factor / risk premium as a multiplicative cost adjustment.

use serde::{Deserialize, Serialize};

/// Bit flags on `Edge.chokepoint_flags` (wired in graph.rs in v0.2).
pub mod flag {
    pub const SUEZ: u8 = 1 << 0;
    pub const PANAMA: u8 = 1 << 1;
    pub const RED_SEA: u8 = 1 << 2;
    pub const BAB_EL_MANDEB: u8 = 1 << 3;
    pub const HORMUZ: u8 = 1 << 4;
    pub const MALACCA: u8 = 1 << 5;
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct ChokepointState {
    /// Ever-given 2021, Houthi diversions 2024, etc. `false` forces reroute
    /// via Cape of Good Hope (handled at graph-construction time, not here).
    pub suez_open: bool,
    /// Panama Canal drought-driven throughput reduction 2023–. Cost
    /// multiplier applied to Panama-flagged arcs.
    pub panama_capacity_factor: f32,
    /// Red Sea risk premium from Houthi attacks 2024. Multiplier applied
    /// to Red Sea-flagged arcs (insurance + rerouting overhead).
    pub red_sea_risk: f32,
    /// Generic residual (default 1.0).
    pub global_risk: f32,
}

impl Default for ChokepointState {
    fn default() -> Self {
        Self {
            suez_open: true,
            panama_capacity_factor: 1.0,
            red_sea_risk: 1.0,
            global_risk: 1.0,
        }
    }
}

impl ChokepointState {
    /// Resolve the multiplicative cost factor for an edge with the given
    /// chokepoint-flag bitmask.
    pub fn multiplier(&self, flags: u8) -> f32 {
        let mut m = self.global_risk;
        if flags & flag::SUEZ != 0 && !self.suez_open {
            // Closed canal — this should have been pruned at graph build;
            // defensive: make it prohibitively expensive.
            m *= 1e6;
        }
        if flags & flag::PANAMA != 0 {
            m *= self.panama_capacity_factor;
        }
        if flags & (flag::RED_SEA | flag::BAB_EL_MANDEB) != 0 {
            m *= self.red_sea_risk;
        }
        m
    }

    /// Canonical states for the three chokepoint shocks we study.
    pub fn suez_ever_given_2021() -> Self {
        Self { suez_open: false, ..Self::default() }
    }
    pub fn panama_drought_2023() -> Self {
        Self { panama_capacity_factor: 1.8, ..Self::default() }
    }
    pub fn houthi_red_sea_2024() -> Self {
        Self { red_sea_risk: 2.5, ..Self::default() }
    }
}
