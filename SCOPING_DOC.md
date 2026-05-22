# Paper 5 — Effective Distance: A Satellite-Calibrated Bilateral Trade Cost Panel

**Authors:** Elizaveta Gonchar, Ian T. Helfrich

**Target venues:** *Journal of International Economics* (main paper) + *Scientific Data* (data descriptor)

**Status:** Scoping, v0.1 — April 22, 2026

**Backup journal ladder:** RESTAT → JDE → EER

---

## 1. Contribution claim (vs. Gonchar-Helfrich 2025 SSRN baseline)

### What the SSRN paper (5202676) already establishes
- The **OVDL** distance measure (Origin-VIIRS-to-Destination-LandScan population-weighted centroid distance)
- Time-varying within-country activity weighting using annual VIIRS nightlights
- Cross-sectional US interstate application
- Global trade-flow panel prediction improvement vs. static CEPII distw
- The theoretical framing against Head-Mayer (2014) distance puzzle

### What Paper 5 adds as net-new
1. **Multi-modal transport-network distance layer** — OSM-derived least-cost-path over road/maritime/air with endogenous modal choice per pair-year. Subsumes pop-weighted, activity-weighted, and infrastructure-aware distances into a unified framework.
2. **Structural gravity validation in Yotov specification** — three-way FE (exporter-time, importer-time, pair), PPML with Weidner-Zylkin (2021) bias correction, intranational flows $X_{ii}$ included (Heid-Larch-Yotov 2021). The SSRN is prediction-focused; this is identification-focused.
3. **Welfare counterfactuals on named shocks** — Suez 2021 closure (7-day Ever Given), Panama Canal 2023 drought (draught restrictions), Houthi Red Sea 2024 rerouting. ACR formula + Caliendo-Parro exact-hat algebra. Shows the distance panel does consequential work.
4. **Public panel release** — 2000–2024, ~200 countries, multiple distance variants, CC-BY 4.0 on Harvard Dataverse + Zenodo, versioned (`v2026.0`). The SSRN paper does not ship reusable data.
5. **Companion software** — Python (`pyfixest`-based) + R (`fixest::fepois`) + Stata (`ppmlhdfe`) ports; reproducibility notebooks.
6. **Sector-heterogeneous robustness** — HS2 (97 sectors) with sector-specific θ_s from Fontagné-Guimbard-Orefice (2022).

**Paper 5 positioning sentence:** *"We scale the Gonchar-Helfrich (2025) nightlights-weighted distance concept into a full multi-modal satellite-calibrated bilateral trade-cost panel, validate it in the Yotov three-way FE structural gravity framework with Weidner-Zylkin bias corrections, and demonstrate consequential welfare divergence on the 2021–2024 chokepoint shocks. The panel is released as public-good infrastructure."*

---

## 2. Literature skeleton (15 anchor citations; full bib in `refs.bib`)

**Foundations**
- Anderson & van Wincoop (2003, AER) — structural gravity, multilateral resistance
- Head & Mayer (2002, 2014) — illusory border effects; handbook chapter
- Mayer & Zignago (2011, CEPII WP) — distw methodology we benchmark against
- Santos Silva & Tenreyro (2006, RESTAT) — PPML estimator

**Internal-geography / satellite-gravity frontier**
- Ramondo, Rodríguez-Clare & Saborío-Rodríguez (2016, AER) — internal trade frictions
- Donaldson & Hornbeck (2016, QJE) — market access from granular networks
- Allen & Arkolakis (2014, QJE; 2022, ECMA) — universal spatial gravity
- Fajgelbaum & Schaal (2020, ECMA) — optimal transport networks
- Henderson, Storeygard & Weil (2012, AER) — nightlights as economic proxy
- Agnosteva, Anderson & Yotov (2019, EER) — intra-national trade costs

**Estimation frontier**
- Yotov, Piermartini, Monteiro & Larch (2016, WTO/UNCTAD) — Advanced Guide
- Weidner & Zylkin (2021, JIE) — three-way PPML bias correction
- Correia, Guimarães & Zylkin (2020) — ppmlhdfe separation diagnostics
- Larch & Yotov (2024, World Economy) — evolution of structural gravity

**Welfare**
- Arkolakis, Costinot & Rodríguez-Clare (2012, AER) — ACR welfare formula
- Caliendo & Parro (2015, RESTud) — exact-hat tariff counterfactuals

**Competitors to benchmark against (not cite for novelty)**
- Gonchar & Helfrich (2025, SSRN 5202676) — **own prior; must cite as the direct predecessor**
- USITC ITPD-E-R03 / Dynamic Gravity (Borchert-Larch-Shikher-Yotov 2025)
- ESCAP-WB Trade Cost Database (Arvis-Duval-Shepherd-Utoktham, 2024 update)
- geopoliticaldistance.org (2025)
- IMF WP 25/93 "Nowcasting Global Trade from Space" (2025)
- Brancaccio, Kalouptsidi & Papageorgiou (2020, ECMA) — AIS maritime costs

---

## 3. Data stack

### Trade flows
- **BACI HS02 V202401b** (CEPII): 2002–2022, ~10.8M rows/year. Local at `/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/`.
- Aggregated to country-pair-year main spec; country-pair-year-HS2 robustness.
- Intranational flows: OECD TiVA ICIO for $X_{ii}$; UN National Accounts gross output for gap-filling.

### Gravity covariates
- **CEPII Gravity V202211** — 4.7M rows, 87 columns, 1948–2021. Local at `/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/`. Note: already contains custom `distw_harmonic_jh` / `distw_arithmetic_jh` variants.
- DESTA for PTA depth; Bailey-Strezhnev-Voeten UN voting; GSDB sanctions.

### Population (gridded)
- **WorldPop Constrained 100m**: primary, CC-BY 4.0, annual 2000–2024.
- **GHS-POP R2023A 100m**: cross-validation, JRC, CC-BY 4.0, epochs 1975/1980/.../2030.
- **HRSL 30m** (Meta/HDX): finer subnational where available, CC-BY 4.0.
- **LandScan 1km**: internal validation only (non-redistributable). Published as correlation statistics, not grids.

### Nighttime lights
- **VIIRS V2 annual composites** (NOAA/EOG): 2012–2024, 500m.
- **Li et al. (2020) harmonized DMSP-VIIRS NTL**: 2000–2018 (extended to 2021 by community). Used for 2000–2011 backcasting.
- Gas-flare masked; stray-light corrected.

### Transport networks (net-new vs SSRN)
- **OpenStreetMap planet** (ODbL): road/rail network. Processed via pyrosm/osmium. Geofabrik regional extracts for compute efficiency.
- **GSHHG coastlines** (LGPL): maritime shoreline.
- **World Port Index** (NGA, public domain): ~3,700 port nodes.
- **OpenFlights** (ODbL): airport graph.
- **CERDI Sea Distance** (CEPII, free): port-pair sea distances including Suez/Panama routing.
- **Kiel Trade Indicator** (IfW, free): container-flow validation.
- **GlobalFishingWatch AIS** (CC-BY-NC): vessel-track calibration.

### Boundaries
- **geoBoundaries ADM0/ADM1** (CC-BY 4.0) — for public redistribution. Avoids GADM non-commercial restriction.

### Network centrality (already built)
- **OneDrive `ThesisData_Last/tables/`** — 2002–2022 panel: Degree, Closeness, Betweenness, Eigenvector, PageRank. Plugs in as centrality-adjusted distance variant.

---

## 4. Methodology

### 4.1 Distance variants released

| Variant | Formula / construction | Data |
|---|---|---|
| `d_cepii_distw` | CEPII Head-Mayer (2002), θ=1 | Benchmark (for comparison) |
| `d_ovdl` | SSRN 2025 OVDL: pop-weighted + VIIRS-weighted centroid pair, great-circle | WorldPop + VIIRS |
| `d_ovdl_t` | Annualized OVDL with yearly WorldPop + VIIRS updates | WorldPop + VIIRS |
| `d_lcp_road` | Least-cost path on OSM road graph, pop-endpoint-weighted | OSM + WorldPop |
| `d_lcp_multi` | Multi-modal LCP (road + maritime + air), endogenous modal choice, pop-endpoint-weighted | OSM + GSHHG + WPI + OpenFlights |
| `d_lcp_multi_t` | Annualized multi-modal LCP (core time-varying output) | Full stack |
| `d_lcp_act` | Multi-modal LCP with VIIRS activity-weighted endpoints (instead of pop) | Full stack + VIIRS |
| `d_net_adj` | Centrality-adjusted `d_lcp_multi_t` (network-effective distance) | Full stack + centrality panel |

Primary specification: `d_lcp_multi_t`. All eight released in panel.

### 4.2 Distance formula (LCP multi-modal)

For country pair $(i,j)$ at year $t$:

$$d_{ij,t}^{\text{multi}} = \left[\sum_{k\in i}\sum_{\ell\in j}\, w_{k,t}^{(i)}\, w_{\ell,t}^{(j)}\, \left(\min_{m\in M}\, c_m^{\text{LCP}}(k,\ell; t)\right)^{\theta}\right]^{1/\theta}$$

where $k,\ell$ index agglomeration nodes within $i,j$; $w_{k,t}^{(i)} = \text{pop}_{k,t} / \sum_{k'\in i}\text{pop}_{k',t}$ (or VIIRS-weighted analog); $M = \{\text{road, maritime, air}\}$; $c_m^{\text{LCP}}$ is the time-to-cost transformed least-cost path on mode $m$'s graph (with Suez/Panama/Red-Sea arcs toggleable for counterfactuals); $\theta \in \{-1, +1\}$ reported both.

Agglomeration nodes: top-$N_i$ cities in country $i$ by population such that coverage $\geq 80\%$ national pop (cf. Mayer-Zignago 2011 20% threshold; we strengthen).

### 4.3 Structural gravity specification

Canonical Yotov / Larch-Yotov spec:

$$X_{ij,t} = \exp\left(\pi_{it} + \pi_{jt} + \mu_{ij} + \beta \cdot \ln d_{ij,t} + \gamma \cdot Z_{ij,t}\right) + \varepsilon_{ij,t}$$

Estimated by PPML (Santos Silva & Tenreyro 2006) with three-way FE (exporter-time, importer-time, pair) via `pyfixest::fepois` or `ppmlhdfe`. Pair FE absorb time-invariant friction (CEPII distw, language, colony). $d_{ij,t}$ identified from *time variation* — which is where our satellite-calibrated panel earns identification.

Include intranational observations $X_{ii,t}$ (Heid-Larch-Yotov 2021). Control $Z_{ij,t}$: time-varying PTA entry, sanctions on/off, contig × globalization-dummy-time interaction.

Bias correction: Weidner-Zylkin (2021) via `ppml_fe_bias` Stata; we compute in parallel and cross-report.

Sector-heterogeneous: repeat for each HS2, stack results, report pooled + sector-dispersed θ_s (Fontagné-Guimbard-Orefice 2022 values as priors).

### 4.4 Welfare counterfactuals

For shock $s$ (Suez 2021, Panama 2023, Red Sea 2024):
1. Re-compute $d_{ij,t}^{\text{multi}}$ with affected arc(s) removed/restricted.
2. Compute counterfactual trade flows via exact-hat algebra (Caliendo-Parro 2015).
3. ACR welfare change: $\hat{W}_i = \hat{\lambda}_{ii}^{-1/\theta}$.
4. Report country-level welfare deltas, GDP-weighted world average, and the set of winners/losers.
5. Validate against observed trade-flow anomalies in the shock quarter/year.

---

## 5. Identification strategy

The econometric hook: **CEPII distw is time-invariant, so pair FE absorb it. Our time-varying satellite-calibrated distance is identified in deviation from pair FE.** This sidesteps the primary critique of distance-in-gravity papers (endogeneity with pair-specific unobservables).

Further identification from **shock-based within-pair variation**: Suez/Panama/Red Sea closures generate plausibly exogenous arc-level distance variation for pairs whose least-cost path runs through the affected chokepoint. This supports a short-run IV / reduced-form check on the distance coefficient.

Pre-registered hypothesis (for publication credibility):
- **H1**: `d_lcp_multi_t` has significantly stronger gravity fit (AIC, pseudo-R², out-of-sample RMSE) than CEPII distw within pair-FE specification.
- **H2**: Chokepoint-affected pairs show distance-coefficient divergence in shock years consistent with arc-level distance changes.
- **H3**: ACR welfare counterfactuals using `d_lcp_multi_t` diverge from CEPII-based counterfactuals by ≥ X% on the three shocks, with divergence concentrated in exposed countries.

---

## 6. Software and compute

- **Primary language**: Python.
- **Gravity**: `pyfixest` (GPU demeaning via CuPy/LSMR, 2025). Fallback: `fixest::fepois` in R.
- **Stata cross-check**: `ppmlhdfe`, `ppml_fe_bias`, `ge_gravity2` for bias correction and GE counterfactual validation.
- **Network routing**: `networkx` + `osmnx` for moderate scale; `graph-tool` or `cugraph` (GPU) for full planet OSM.
- **Raster**: `rioxarray` + `xarray` + `dask` for WorldPop/GHS/VIIRS processing.
- **Storage**: Parquet + DuckDB for tabular; Zarr for raster time-series.
- **Reproducibility**: `uv` or `pixi` for Python env; `renv` for R; `docker` image for full reproducibility.
- **Compute budget**: ~200 GB intermediate storage; GPU optional but speeds demeaning 5-10×.

---

## 7. Dataset release plan

- **Hosting**: Harvard Dataverse (citation of record) + Zenodo (mirror, concept-DOI + version-DOI) + GitHub Pages landing.
- **License**: CC-BY 4.0 (data), MIT (code).
- **Versioning**: Calendar (`v2026.0`); annual minor updates.
- **Contents**: Panel CSV + Parquet (country-pair-year × 8 distance variants); codebook; `datapackage.json` (Frictionless); replication notebooks; CITATION.cff; SUCCESSION.md (name successor — Elizaveta).
- **Companion paper**: *Scientific Data* descriptor (~2000 words, DOI-anchored).
- **Launch**: NEP-INT + NEP-GEO announcements; VoxEU column pitch; direct outreach to Head, Mayer, Redding, Donaldson, Allen, Arkolakis, Yotov, Larch, Borchert.

---

## 8. Timeline

### Path A — Sprint WP by April 30, 2026 (v0.9) **[selected]**

| Day | Deliverable |
|---|---|
| Apr 22 (today) | Scoping doc finalized; bib seeded; folder structure set. |
| Apr 23 | Literature skeleton with 60 cites; methods section draft; decision on 100 vs 200 countries for v0.9. |
| Apr 24–25 | Transport-network pipeline: OSM extracts, GSHHG maritime, WPI ports, OpenFlights air. Core country coverage. |
| Apr 26 | LCP multi-modal distance computation for 100-country core × 2002–2022. Validation vs CEPII distw correlation. |
| Apr 27 | Modernize `estimate_gravity_models()` to `pyfixest` three-way FE; run on BACI HS02 aggregate. |
| Apr 28 | Weidner-Zylkin bias check via Stata bridge. Run Suez 2021 counterfactual. |
| Apr 29 | Draft intro, empirics, results; tables 1–4; figures 1–3. |
| Apr 30 | Abstract, polish, SSRN v0.9 submission; NEP-INT announcement. |

**v0.9 scope constraints:**
- 100-country core (OECD + top-50 BACI by trade value), not full 200
- Road + maritime + air; rail + pipeline deferred
- One counterfactual (Suez 2021); Panama + Red Sea in v1.0
- pyfixest primary; Weidner-Zylkin cross-check reported where available

### Path B — v1.0 submission-ready WP by June 1, 2026

| Week | Deliverable |
|---|---|
| May 1–7 | Scale to full 200 countries; add rail + pipeline modes. |
| May 8–14 | Full three counterfactuals (Suez, Panama, Red Sea); HS2 sector-heterogeneous θ_s robustness. |
| May 15–21 | Complete Weidner-Zylkin cross-check panel; sensitivity analyses. |
| May 22–28 | Scientific Data descriptor drafted; Dataverse + Zenodo prepared; CITATION.cff / SUCCESSION.md / datapackage.json. |
| May 29–Jun 1 | Final polish; v1.0 SSRN update; JIE submission; Dataverse public release with DOI; VoxEU pitch. |

### Post-v1.0 roadmap
- Paper 6: subnational Effective Distance (US, EU, India, China, Brazil) — reuses the pipeline
- Paper 7 (revival): FunProject TDA, now with real data via this pipeline
- Paper 8+: labor/skill redistribution layer (DOSE + HRSL); environmental-exposure-via-trade layer

---

## 9. Risks and mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| IMF or CEPII releases competing panel first | Medium | High | Ship v0.9 by Apr 30; claim the SSRN Gonchar-Helfrich 2025 seed as priority date. |
| OSM planet LCP compute blows up | Medium | Medium | Scope to 100-country core for v0.9; use Geofabrik extracts + DuckDB spatial; GPU cugraph if needed. |
| LandScan license issue flagged at submission | Low | Medium | Public panel uses WorldPop + GHS-POP; LandScan only reported as internal correlation statistic. |
| Weidner-Zylkin Stata bridge painful | Medium | Low | Run on aggregate first; HS2 robustness may not need WZ for v0.9. |
| Elizaveta's bandwidth / job market | Unknown | Medium | Ian confirm Apr 22; agree on work split; she leads spatial-econometrics layer, Ian leads transport + gravity pipeline. |
| Reviewer pushes for AIS / satellite ship-tracking | Medium | Low | GFW AIS aggregates cited; full AIS as explicit v2.0 extension. |
| Distance-elasticity estimate shifts implausibly | Low | High | Pre-registered H1–H3; sensitivity to θ, node threshold, transport-cost weights. |

---

## 10. Open questions for Elizaveta

1. Primary empirical responsibility split on Paper 5? Suggest: Ian on transport-network + gravity pipeline; Elizaveta on spatial-econometrics layer + US subnational extension.
2. Authorship order: Gonchar & Helfrich (SSRN order) or Helfrich & Gonchar (alphabetical sometimes, Ian-led on this extension)?
3. Target journal agreement: JIE or push higher to RESTud / AER? Recommend JIE for pace.
4. Any independent follow-up plans that might overlap with Paper 5 scope (from her research page)?
5. Endorsement for NBER submission (if either has affiliation)?

---

## 11. Immediate next actions (today/tomorrow)

- [ ] Ian: confirm Path A (v0.9 by Apr 30) with Elizaveta; get her sign-off on scope + authorship.
- [ ] Ian: share SSRN DOCX source and any existing figures/tables.
- [x] Claude: seed `refs.bib` with the 15 anchors + 45 supporting cites.
- [ ] Claude: draft literature review section (Section 2 of manuscript) — **Day 2**.
- [x] Claude: write data-loading pipeline stub (`data_loaders.py`).
- [x] Claude: write LCP pipeline architecture doc (`ARCHITECTURE.md`).
- [x] Claude: scaffold `paper5-core` Rust crate + PyO3 bindings (`crates/`). Compiles and tests pass.
- [ ] Claude: modernize `estimate_gravity_models()` — fork `TradeModels_Final.py` to `paper5.gravity` with pyfixest — **Day 2**.

### Architectural decision logged 2026-04-22

Stack is **polyglot by design**: Python outer shell (data I/O, econometrics,
plotting) + Rust inner loop (graph build, LCP solve). Scaffold landed today:

- `Cargo.toml` workspace, `crates/paper5-core` (solver) + `crates/paper5-py` (PyO3 bindings).
- Plain Dijkstra prototype + chokepoint state machine + criterion bench.
- `tests/test_ffi_roundtrip.py` validates the Rust↔Python seam.
- Rationale + fallback tiers documented in `ARCHITECTURE.md`.

Cost of this decision: ≈ 1 extra day of scaffolding + build-system complexity.
Payoff: full-panel LCP run in ~70 min instead of ~7 days, plus a publishable
Rust artifact for the Scientific Data companion descriptor.

---

**This doc is living. Edit freely; I'll update as we iterate.**
