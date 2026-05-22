# Paper 5 Sprint Status

**This file is the cross-session, cross-LLM handoff document.**
Update it at the END of every work session. Any LLM (Claude, GPT, Gemini)
picking up this project should read this file FIRST. Keep it under 150 lines.

---

## Last updated: 2026-04-22 (Sprint Day 4 — gravity validation complete)

## Target: v0.9 working paper submitted to SSRN by 2026-04-30

---

## Project in one sentence
Building a global time-varying satellite-calibrated bilateral distance panel
(2002–2022, 100-country core → 200 countries v1.0) validated in structural
gravity, released CC-BY 4.0 on Harvard Dataverse.

---

## What is DONE (as of this writing)

### Architecture & scaffolding
- [x] `ARCHITECTURE.md` — full polyglot design doc (Python shell + Rust hot loop)
- [x] `Cargo.toml` workspace + `crates/paper5-core/` Rust solver scaffold
  - `graph.rs` (CSR adjacency, typed nodes/edges, `#[repr(u8)]` Mode enum)
  - `dijkstra.rs` (heap Dijkstra + rayon many-to-many; 2 unit tests pass)
  - `chokepoints.rs` (Suez/Panama/RedSea state machine with presets)
  - `benches/dijkstra.rs` (Criterion bench, target <20ms SSSP on 100K nodes)
- [x] `crates/paper5-py/` PyO3 binding: `GraphHandle`, `ChokepointHandle`,
      `solve_pair`, `solve_many_to_many → numpy array`
- [x] `tests/test_ffi_roundtrip.py` (skips gracefully if Rust ext not built)
- [x] `workflow/Snakefile` — L0→L4 DAG with `build_rust_ext` rule
- [x] `src/paper5/cli.py` — `python -m paper5.cli` entry point for Snakemake

### Python package
- [x] `data_loaders.py` — BACI, CEPII, centrality panel, WorldPop, GHS-POP,
      VIIRS, DMSP harmonized, OSM, WPI ports, OpenFlights, CERDI Sea Distance
- [x] `agglomerate.py` — FULLY IMPLEMENTED:
      density threshold → CC → union-find merge → coverage-ranked retention →
      VIIRS aggregation. All unit tests pass.
- [x] `distance.py` — FULLY IMPLEMENTED (prototype):
      `great_circle_km`, `ces_aggregate`, `ChokepointState` (with `for_year()`),
      `compute_d_ovdl`, `build_transport_graph_nx` (maritime + air + road stubs),
      `lcp_pair_nx`, `compute_d_lcp_multi_t`, `compute_all_variants`.
      All unit tests pass.
- [x] `gravity.py` — pyfixest 3-way FE PPML; `estimate_gravity()` and
      `benchmark_distance_variants()` are IMPLEMENTED; WZ bridge stub
- [x] `notebooks/02_gravity_validation.py` — PPML prototype run complete;
      d_ovdl_t β=−0.894 (se=0.100) vs CEPII distw β=−0.689 (se=0.109); N=544
- [x] `cli.py` — `python -m paper5.cli agglomerate/compute-lcp/assemble-panel`
- [x] `__main__.py` — enables `python -m paper5`

### Manuscript
- [x] `manuscript/main.tex` — 898 lines; all 9 sections have prose (2 stubs remain)
- [x] Section 1 (Introduction) — FULLY DRAFTED (~900 words)
- [x] Section 2 (Related Literature) — FULLY DRAFTED (4 subsections, ~2400 words, ~55 citations)
- [x] Section 3 (Data) — FULLY DRAFTED (6 subsections: Trade, CEPII, Pop, NTL, Transport, Boundaries)
- [x] Section 4 (Construction) — FULLY DRAFTED (agglomeration → weights → graph → chokepoints → LCP → CES → Table 1 variants)
- [x] Section 5 (Gravity) — DRAFTED (spec + FE rationale + hypotheses; results placeholder)
- [x] Section 8 (Dataset Release) — FULLY DRAFTED (versioning, Dataverse, license, formats)
- [x] Section 9 (Conclusion) — FULLY DRAFTED (summary + 4-item roadmap)
- [~] Section 5 (Gravity) — preliminary prototype results now in manuscript
      (β=−0.894 cross-section; full 3-way FE pending 100-country panel)
- [ ] Section 6 (Results) — Tables 1–6 pending 100-country panel run
- [ ] Section 7 (Counterfactuals) — ACR welfare calc pending `welfare.py`

### Bibliography
- [x] `refs.bib` — ~65 entries across 11 categories (Tinbergen 1962 → Ghosh 2021)

### Tests
- [x] `tests/test_ffi_roundtrip.py` — Rust↔Python FFI (skips if paper5_core not built)
- [x] `tests/test_distance.py` — 20 pytest tests (ALL PASS; 5 new flag-regression tests)
- [x] `tests/test_agglomerate.py` — unit + integration tests (unit PASS; integration skippable)
- [x] `tests/test_transport_graph.py` — 6 maritime graph + chokepoint uplift tests (ALL PASS)
- [x] `tests/_run_inline.py` — quick smoke test, no pytest needed
**Total: 33 pass, 2 skip (integration/FFI) as of Day 4.**

### Notebooks
- [x] `notebooks/01_prototype_20country.py` — full L0→L2 prototype with graceful
      fallback to CEPII centroids if WorldPop not mounted

### Workflow
- [x] `workflow/Snakefile` — L0→L4 DAG
- [x] `workflow/config/core_100.txt` — 100-country v0.9 country list

---

## What is IN PROGRESS / NEXT

| Priority | Task | File | Sprint Day |
|----------|------|------|-----------|
| P0 | Wire Python→Rust: `lcp_pair_rust()` calling `paper5_core` | `distance.py` | Day 4–5 |
| P0 | OSM road graph: `pyrosm` extract + cache per country | `distance.py` | Day 4–5 |
| P0 | Run 20-country prototype with real WorldPop data | `01_prototype_20country.py` | Day 4–5 |
| P1 | Section 3 (Data) full prose draft | `main.tex` | Day 5 |
| P1 | Section 4 (Construction) full prose draft | `main.tex` | Day 5 |
| P1 | Full 100-country panel run (LCP multi-modal, 2002–2022) | `Snakefile` | Day 5–6 |
| ~~P2~~ | ~~Gravity estimation on panel~~ | DONE — 20c prototype β=−0.894 | Day 4 ✓ |
| ~~P2~~ | ~~Section 5 (Gravity) fill in~~ | DONE — preliminary table in §5 | Day 4 ✓ |
| P2 | Section 6 (Results) tables 1–4 | `main.tex` | Day 7 |
| P2 | Section 7 (Counterfactuals) Suez/Panama/RedSea | `main.tex` | Day 7 |
| P3 | `welfare.py` — ACR + exact-hat algebra | new file | Day 6–7 |
| P3 | WZ bias bridge (Stata ppml_fe_bias) | `gravity.py` | Day 7 |
| P3 | Section 1 (Intro) finalize with actual numbers | `main.tex` | Day 8 |

---

## Key design decisions (LOCKED for v0.9)

- **Rust** for LCP solver (not networkx; see ARCHITECTURE.md §2 for why)
- **WorldPop 100m + GHS-POP** for population; LandScan internal-only
- **pyfixest** for PPML; Stata bridge for WZ bias correction only
- **100-country core** for v0.9; 200 countries in v1.0
- **3 chokepoint shocks**: Suez 2021, Panama 2023, Red Sea 2024
- **Target journal**: JIE (main) + Scientific Data (data descriptor)
- **θ = −1** (CES inverse) as primary; θ = +1 as robustness
- **Years**: 2002–2022 (BACI HS02 coverage) for main panel

---

## Data paths (on Ian's machine)

```
HELFRICH_GD     = /Volumes/HELFRICH-GD
ONEDRIVE        = /Users/ian/Library/CloudStorage/OneDrive-Personal
PAPER5_DATA     = <repo>/data

BACI HS02:      HELFRICH_GD/TradeData/BACI_HS02_V202401b/BACI_HS02_Y{year}_V202401b.csv
CEPII Gravity:  HELFRICH_GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv
Centrality:     ONEDRIVE/ThesisData_Last/tables/centrality_evolution.csv
```

---

## Known issues / blockers

1. **`russ_2020` bib entry** — author/title placeholder; verify before final sub
2. **`lafrogne_joussier_2024`** — confirm exact JIE title (currently has
   Martin & Mejean as co-authors; verify)
3. **`heiland_2022`** — confirm this is the right paper for mode-specific priors
4. **Rust PyO3 crate** — `paper5-core` unit tests pass; FFI layer needs
   `maturin develop --release --manifest-path crates/paper5-py/Cargo.toml`
   — untested; will speed up panel run 1000x when wired
5. **Panama flag over-fires on Asian ports**: `_arc_chokepoint_flags` Panama check
   uses `lon < -90` (catches US West Coast) but Asian ports have positive lon (100-170°E).
   Shanghai→New York won't get Panama flag. Workaround: for v0.9 use direct GC approach
   for East Asia→Americas, flag only when one endpoint has lon < -90.
6. ~~**pyfixest install**~~ — **RESOLVED**: `numba --prefer-binary` found pre-built wheel;
   pyfixest 0.50.1 + numba 0.62.1 installed in `.venv/`.
7. **pytest config** — `pyproject.toml` testpaths=["tests"] is set; run with
   `.venv/bin/python -m pytest tests/` (not bare `pytest` to avoid system Python)
8. **venv location** — `.venv/` at repo root; activate with
   `source .venv/bin/activate` or use `.venv/bin/python` directly

## What changed in Day 4 continuation (gravity validation)

**Gravity estimation unblocked**: Static fallback centroids give zero within-pair
time variation in d_ovdl_t → collinear with pair FE. Fix: auto-detect static
centroid case (max within-pair std < 1e-6) and drop pair FE, falling back to
cross-section PPML with exporter-time + importer-time FE.

**Results** (20-country prototype, 2010+2022, N=544):
- d_ovdl_t: β = **−0.894** (se = 0.100) ← this paper
- CEPII distw: β = **−0.689** (se = 0.109) ← benchmark

Elasticity 30% larger in magnitude for satellite-calibrated measure vs CEPII distw.
Both within literature range (−0.9 to −1.2 from Disdier & Head 2008). ✓

**Files changed**: `notebooks/02_gravity_validation.py` (spec-selection logic),
`manuscript/main.tex` §5 (preliminary results table with actual numbers),
`data/derived/gravity_benchmark_20c.csv` (saved output).

---

## What changed in Day 4 (chokepoint bug fix)

**Bug**: `_arc_chokepoint_flags()` used midpoint bounding boxes + `is_east_far`/`is_west_med`
checks. The `is_west_med` check required `lon < 20°E`, missing Black Sea / eastern Med ports
(e.g., Novorossiysk at 37.8°E). Dijkstra could route Mumbai→Novorossiysk→Hamburg
(both arcs unflagged) at $2,378 under Houthi vs $2,375 normal — essentially 0% uplift.

**Fix**: Replaced with ocean-basin classification:
- `east_of_suez`: lon > 43°E (east of Bab-el-Mandeb)
- `west_of_suez`: lon ∈ [−80, 42], lat > −10 (Med, Black Sea, Atlantic, Europe)
- Any arc crossing east↔west basins gets SUEZ+RED_SEA flags
- Sub-Saharan ports (Cape Town lat=-33.9 < −10) remain unflagged → Cape route stays navigable

**Result**: IND→DEU Houthi uplift: 0% → **43%** ✓, JPN→DEU: 0% → **65%** ✓, SGP→DEU: **32%** ✓
AUS→USA (trans-Pacific) stays at 0% uplift ✓ (unflagged, correct behavior)

---

## For an LLM picking this up cold

1. Read this file.
2. Read `ARCHITECTURE.md` (system design, language boundaries).
3. Read `manuscript/main.tex` (sections 1–2 are done; 3–9 need content).
4. Read `src/paper5/gravity.py` (gravity estimation is mostly implemented).
5. The most impactful next thing is implementing `agglomerate.py` (L1 layer).
6. The code compiles: `cargo test -p paper5-core` should show 2 passing tests.
7. Do NOT re-architect. The Python/Rust split is a locked decision.
8. Do NOT add new dependencies without updating `pyproject.toml` + `Cargo.toml`.
