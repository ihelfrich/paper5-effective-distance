# Paper 5 — Effective Distance: System Architecture

**Version:** 0.1 (2026-04-22, Sprint Day 1)
**Authors:** Ian T. Helfrich, Elizaveta Gonchar

This document specifies the computational architecture for the Paper 5 bilateral
distance panel. It is opinionated about language choice by layer: **Python is
the outer shell** (data I/O, econometrics, plotting, notebooks), **Rust is the
inner loop** (multi-modal graph construction, LCP solve, raster→agglomeration
aggregation). A thin PyO3 / `maturin` bridge stitches them together. DuckDB is
the storage substrate for every tabular artifact.

The target throughput for v0.9 (Apr 30): compute all eight distance variants
for a 100-country core × 2002–2022 (21 years) in **under 6 hours of wall time
on an M-series MacBook Pro**, so we can iterate on weighting parameters θ and
mode-cost calibrations without a full-day turnaround.

---

## 1. Layered Pipeline

```
 ┌─────────────────────────────────────────────────────────────────────┐
 │  L5  MANUSCRIPT / FIGURES          (LaTeX, matplotlib, plotly)      │
 │       manuscript/main.tex, figures/*.pdf                            │
 ├─────────────────────────────────────────────────────────────────────┤
 │  L4  GRAVITY ESTIMATION + WELFARE  (Python: pyfixest, Stata bridge) │
 │       src/paper5/gravity.py, src/paper5/welfare.py                  │
 ├─────────────────────────────────────────────────────────────────────┤
 │  L3  DISTANCE VARIANT ASSEMBLY     (Python: polars, duckdb)         │
 │       src/paper5/distance.py                                        │
 ├─────────────────────────────────────────────────────────────────────┤
 │  L2  LCP SOLVER + GRAPH BUILD     ★ RUST (paper5-core crate) ★      │
 │       crates/paper5-core/                                           │
 ├─────────────────────────────────────────────────────────────────────┤
 │  L1  RASTER → AGGLOMERATIONS       (Python: rioxarray, dask  or     │
 │                                      Rust: geo + ndarray)           │
 │       src/paper5/agglomerate.py  OR  crates/paper5-agglom/          │
 ├─────────────────────────────────────────────────────────────────────┤
 │  L0  DATA LOADERS                  (Python: pandas, rioxarray)      │
 │       src/paper5/data_loaders.py                                    │
 └─────────────────────────────────────────────────────────────────────┘
```

**Rule of thumb:** if the layer touches every pixel of a global raster, or
every arc of a continental graph, it lives in Rust. If it touches one row per
country-pair-year (~500K rows), it lives in Python.

---

## 2. Why Rust for L2 (and probably L1)

I (Ian) am going to be blunt about the math:

- **Graph size.** Multi-modal graph = road arcs from OSM (≈ 40M edges
  planet-wide, prune to 8–12M for arterials + highways + major secondary) +
  maritime arcs from GSHHG coastline graph + WPI port connectors (≈ 150K edges
  after port-to-port great-circle + chokepoint routing) + air arcs from
  OpenFlights routes (≈ 60K). Total ≈ 10M arcs, ≈ 3M nodes.
- **Queries.** 100 countries × 100 countries = 10K pairs. Each pair is
  actually a bipartite many-to-many query between origin-country agglomerations
  (≈ 20 per country average) and destination agglomerations (≈ 20). So ≈
  4M shortest-path queries per year × 21 years = **84M SSSP/A* calls**.
- **networkx.** Dijkstra on 10M edges = ~15 seconds per source in pure Python.
  84M calls / amortized ~20 sources per country-year = 4.2M effective
  source-runs. At 15s each = **700 days**. Hard no.
- **graph-tool / igraph (C++).** ~150ms per source. Gets us to ~7 days.
  Workable but tight, and the API is awkward for the many-to-many
  agglomeration-bundle pattern.
- **Rust with a domain-tuned contraction hierarchy or multi-level Dijkstra
  (RoutingKit-style).** ~1ms per source after preprocessing. Full panel in
  ~70 minutes. This is the target.

We do **not** need to reimplement RoutingKit. We have two credible options:

1. **Link against RoutingKit directly** (C++ library, MIT license) via a thin
   Rust wrapper. Proven on planet-scale OSM.
2. **Use the `petgraph` + custom bidirectional Dijkstra approach in Rust**,
   with node contraction we write ourselves. More work, but keeps the build
   hermetic and gives us full control over multi-modal edge-cost functions
   (time × mode-specific freight rate + chokepoint premiums).

**Decision for v0.9:** start with option 2 — a hand-rolled Rust solver using
`petgraph` for prototyping and a dense CSR adjacency for the hot loop. Option
1 is a v1.0 upgrade if we hit latency problems.

---

## 3. Crate Layout (Rust)

```
crates/
  paper5-core/            # LCP solver + graph types
    src/
      lib.rs              # pub API: build_graph, solve_pair, solve_all
      graph.rs            # CSR adjacency, multi-modal edge struct
      dijkstra.rs         # bidirectional Dijkstra + A* variants
      contraction.rs      # node contraction for CH-style preprocessing
      chokepoints.rs      # Suez/Panama/Red Sea state machine
      ffi.rs              # PyO3 bindings (feature-gated)
    Cargo.toml
  paper5-agglom/          # raster → agglomeration centroids (OPTIONAL, Phase 2)
    src/
      lib.rs
      worldpop.rs         # stream TIFF tiles, population-weighted clustering
      viirs.rs            # radiance aggregation per agglomeration footprint
    Cargo.toml
  paper5-py/              # PyO3 wheel: `import paper5_core`
    src/lib.rs
    Cargo.toml
    pyproject.toml        # maturin build
```

Key dependencies:
- `petgraph` — graph primitives (prototyping)
- `rayon` — parallel SSSP across source nodes
- `pyo3` + `maturin` — Python bindings
- `gdal` / `geo` — raster I/O (if we push L1 to Rust)
- `serde` + `bincode` — serialize preprocessed graph to disk

Build:
```bash
cd crates/paper5-py
maturin develop --release    # installs paper5_core into the active .venv
```

---

## 4. Data Contract Between Layers

All cross-layer data travels as **Parquet files** (Arrow in-memory where
possible). This keeps Python↔Rust boundary zero-copy and makes every
intermediate inspectable from DuckDB.

### L0 → L1: raster pointers
```
data/worldpop/ppp_{year}_100m_Aggregated.tif     # 40GB per year, already on HELFRICH-GD
data/ghs_pop/GHS_POP_E{year}_GLOBE_R2023A_54009_100.tif
data/viirs/VNL_v2_npp_{year}_global.tif
data/boundaries/gadm_410.gpkg                   # country polygons
```

### L1 → L2: agglomerations table
**Parquet schema** `data/derived/agglomerations_{year}.parquet`:
```
iso3: str               # ISO-3166 alpha-3
agglom_id: u32          # within-country index
lon: f32, lat: f32      # WGS84 centroid
pop: f64                # WorldPop-summed population
viirs_radiance: f32     # VIIRS sum over footprint
ghs_pop: f64            # GHS-POP validation
coverage_share: f32     # share of national pop represented
```
Target rows: ~5K per year (≈ 240 countries × ≈ 20 agglomerations).

### L2 → L3: pairwise LCP costs
**Parquet schema** `data/derived/lcp_{year}.parquet`:
```
iso_o: str, iso_d: str
agglom_o: u32, agglom_d: u32
mode: str               # "road" | "maritime" | "air" | "multi"
cost_hours: f32         # time-equivalent cost
cost_km: f32            # km-equivalent (for legacy comparability)
path_length_km: f32     # actual geographic path length
chokepoint_flags: u8    # bitmask: Suez=1, Panama=2, RedSea=4
```
Target rows: ~400M across 21 years (shard by year).

### L3 → L4: country-pair-year distance panel
**Parquet schema** `data/derived/distance_panel.parquet`:
```
iso_o: str, iso_d: str, year: u16
d_cepii_distw: f32
d_ovdl: f32
d_ovdl_t: f32
d_lcp_road: f32
d_lcp_multi: f32
d_lcp_multi_t: f32
d_lcp_act: f32
d_net_adj: f32
```
Target rows: 240² × 21 = 1.2M. Trivial.

### L4: gravity input
Join distance panel with BACI + CEPII covariates → wide DataFrame for
pyfixest. Already fits in memory.

---

## 5. Storage: DuckDB as the Control Plane

All Parquet files are registered as DuckDB views via a single
`data/catalog.duckdb` database file. This gives us:

- Ad-hoc SQL against any intermediate.
- Cross-year aggregations with no copy.
- A diff-able catalog file to check into git (the `.duckdb` itself stays out
  — only the `catalog.sql` schema-definition script is tracked).

```sql
-- data/catalog.sql
CREATE VIEW agglomerations AS
  SELECT * FROM read_parquet('data/derived/agglomerations_*.parquet',
                             filename=true);
CREATE VIEW lcp AS
  SELECT * FROM read_parquet('data/derived/lcp_*.parquet', filename=true);
CREATE VIEW distance_panel AS
  SELECT * FROM read_parquet('data/derived/distance_panel.parquet');
```

---

## 6. Reproducibility & Determinism

- **Seed.** All graph contractions and any stochastic subsampling seeded via
  `RNG_SEED=20260422` (today's date).
- **Hash-stamped outputs.** Every Parquet writer stamps
  `git_sha`, `crate_version`, `input_hash` into file metadata.
- **Snakemake DAG.** `workflow/Snakefile` wires L0→L4 with dependency-driven
  recomputation. `snakemake --forcerun compute_lcp_2024` re-does only 2024.
- **CI.** GitHub Actions runs a 5-country mini-panel (USA, MEX, CAN, DEU,
  CHN) × 2 years end-to-end on every PR, ≈ 3 minutes.

---

## 7. What Stays in Python (and why)

| Layer | Task | Why Python |
|------|------|-----------|
| L0 | BACI/Gravity CSV ingest | pandas/polars idiomatic; done once |
| L3 | variant assembly | joins + groupby, polars is fast enough |
| L4 | pyfixest gravity | the only credible PPML lib with 3-way FE and GPU demeaning |
| L4 | Stata bridge for WZ bias | `ppml_fe_bias` is Stata-only; shell out |
| L5 | plotting | matplotlib/plotly, no contest |

We will NOT rewrite pyfixest. We will NOT rewrite anything in L0/L3/L4/L5.

---

## 8. Risk: Rust Time-Budget

Risk: building a usable multi-modal LCP solver in Rust inside 8 days is
aggressive. Mitigation tiers:

1. **Day 3 prototype (Python/networkx, 20-country subset).** Prove the
   pipeline end-to-end with `d_lcp_road` only, 2022 only. This validates L0→L4
   data contracts without waiting on Rust.
2. **Day 4–5 Rust v1 (road only).** Get `paper5_core::solve_pair` working for
   the road graph. Benchmark vs networkx on the 20-country subset. Must beat
   networkx by ≥20× to justify; otherwise fall back to igraph/C++.
3. **Day 6 Rust v2 (multi-modal).** Add maritime + air layers.
4. **Day 7 full panel run + validation.**
5. **Day 8 estimation + manuscript.**

**Fallback.** If Rust is not ready by EoD Day 5, we ship v0.9 with
`d_lcp_multi_t` computed via `graph-tool` (Python bindings over Boost.Graph,
C++ under the hood). Slower preprocessing but same result; we lose the
"novel Rust crate on crates.io" Scientific Data hook.

---

## 9. Multi-Language Boundary Rules

To keep this from turning into a Tower of Babel:

- **No Rust that Python couldn't call.** Every Rust function we write has a
  PyO3 binding. No orphan Rust binaries.
- **No Python in the hot loop.** If a function is called > 10K times per
  variant-year, it lives in Rust.
- **DuckDB/Parquet is the lingua franca.** No pickled objects across language
  boundaries. No Arrow RecordBatches held in memory across a Rust↔Python call
  for more than one function invocation.
- **One test per FFI surface.** `tests/test_ffi_roundtrip.py` round-trips a
  10-row table through every exposed Rust function and checks byte-equality.

---

## 10. Open Architecture Questions

1. **GPU?** `pyfixest` supports CUDA demeaning. Worth wiring on the gravity
   step? Probably not for 1.2M rows; revisit if we add HS2 (100× rows).
2. **Distributed?** No. Single-node. The whole point of the Rust hot loop is
   that we don't need Spark/Ray.
3. **WASM?** Yes, eventually — compile `paper5-core` to WASM so the
   Scientific Data companion paper ships with an in-browser demo that routes a
   user-selected pair live. Post-v0.9 work.
4. **Does L1 need to be Rust?** Agglomeration clustering over a 40GB raster is
   borderline. Python+dask does this in ~10 min per year; Rust would be ~1
   min. We accept Python here for v0.9; revisit if we iterate on clustering
   parameters.

---

## 11. Immediate Next Steps (Sprint Day 2)

1. `cargo new --lib crates/paper5-core` — scaffold the Rust crate.
2. Define `Graph`, `Edge`, `Mode` types in `graph.rs`.
3. Port the networkx-free toy example (10 nodes, 20 edges) as an
   integration test.
4. Stand up `maturin develop` so `import paper5_core` works from the notebook.
5. Python-side: write `agglomerate.py` with WorldPop + Natural Breaks
   clustering per country.

---

_End of ARCHITECTURE.md v0.1. Updates tracked in git._
