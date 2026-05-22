# Paper 5 — Effective Distance

**Satellite-Calibrated Bilateral Trade-Cost Panel (2000–2024)**

Elizaveta Gonchar · Ian T. Helfrich

---

Scale-up of Gonchar & Helfrich (2025), *"Trade in the Spotlight"* (SSRN 5202676) into a full multi-modal satellite-calibrated bilateral distance panel, validated in the Yotov three-way FE structural gravity framework, with welfare counterfactuals on the 2021–2024 chokepoint shocks (Suez, Panama, Red Sea). Released as public-good infrastructure.

## Repo layout (polyglot: Python outer shell + Rust inner loop)

```
Paper5_EffectiveDistance/
├── SCOPING_DOC.md        # Full project plan (living doc)
├── ARCHITECTURE.md       # System architecture + language-choice rationale
├── PROJECT_BRIEF.md      # Self-contained one-pager
├── README.md             # You are here
├── refs.bib              # Bibliography
├── manuscript/           # LaTeX source
│   └── main.tex
├── pyproject.toml        # Python package metadata
├── Cargo.toml            # Rust workspace root
├── crates/                       # Rust workspace
│   ├── paper5-core/              # Multi-modal LCP solver (the hot loop)
│   │   ├── src/{graph,dijkstra,chokepoints}.rs
│   │   └── benches/dijkstra.rs   # Criterion perf bench
│   └── paper5-py/                # PyO3 bindings → `import paper5_core`
│       └── src/lib.rs
├── src/paper5/                   # Python package
│   ├── data_loaders.py           # BACI, CEPII, WorldPop, VIIRS, OSM
│   ├── agglomerate.py            # raster → country agglomeration centroids
│   ├── distance.py               # eight distance-variant assembly
│   ├── gravity.py                # pyfixest three-way FE PPML
│   └── counterfactual.py         # ACR + Caliendo-Parro exact-hat (pending)
├── workflow/                     # Snakemake DAG (L0→L4 orchestration)
├── data/                         # Intermediate parquet/zarr (gitignored)
├── notebooks/                    # Exploration + replication demos
├── figures/                      # Output figures
├── tests/                        # pytest + Rust↔Python FFI round-trip tests
└── refs/                         # Reading notes, annotated PDFs
```

Rust handles the graph build + LCP solve (≈ 10M arcs, ≈ 84M SSSP calls across
the 21-year panel — infeasible in pure Python). Everything else — data
ingest, distance-variant assembly, gravity estimation, plotting — is Python.
DuckDB + Parquet is the lingua franca across language boundaries. See
`ARCHITECTURE.md` for the full rationale.

## Build

```bash
# Python side
pip install -e '.[dev]'

# Rust extension (requires rustc + maturin)
pip install maturin
maturin develop --release --manifest-path crates/paper5-py/Cargo.toml

# Rust perf benchmark
cargo bench -p paper5-core
```

## Target venues
- Main: *Journal of International Economics*
- Data descriptor: *Scientific Data*
- Backup ladder: RESTAT → JDE → EER

## Timeline
- **v0.9 WP: April 30, 2026** (sprint — 100-country core, Suez counterfactual)
- **v1.0 WP + JIE submission: ~June 1, 2026** (full 200 countries, all three counterfactuals, WZ bias corrections, *Scientific Data* descriptor)

## Data release (planned)
- Harvard Dataverse (citation of record) + Zenodo mirror
- CC-BY 4.0 data, MIT code
- Concept DOI + version DOI, calendar versioning (`v2026.0`)

See `SCOPING_DOC.md` for the full plan.
