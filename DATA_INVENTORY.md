# Data inventory — what we have access to, and where

**Status:** v0.1, surveyed 2026-05-20. Update whenever a new dataset is added.

This document is the single source of truth for "do we have it locally."
Before downloading, check here. Before running an analysis on a dataset, verify
its provenance is listed.

---

## 1. Bilateral trade flows

| Dataset | Path | Coverage | Notes |
|---|---|---|---|
| **BACI HS02 V202401b** | `/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/` | 2002–most-recent, bilateral 6-digit HS02 product, USD + qty | Primary country-pair trade data. CSV per year. |
| **BACI HS02 aggregate** | `~/1A_Helfrich_ThesisResearch_May2024/BACI_Trade_NetworkFiles/BACI_HS02_aggregate.csv` | 2002+ | Pre-aggregated to country-pair-year. Faster to load than building from raw. |
| **Master basic BACI** | `~/1A_Helfrich_ThesisResearch_May2024/BACI_Trade_NetworkFiles/master_basic_bacihs02.csv` | Multi-year | Country codes + names + flow totals. |
| **TRADHIST V4** | `/Volumes/HELFRICH-GD/TradeData/TRADHIST_v4.dta` | 1827–2014 long historical | For Paper 5's deep-history robustness if needed; not v0.9 critical. |
| **Tariff elasticities** | `/Volumes/HELFRICH-GD/TradeData/Tariff-Based Product-Level Trade Elasticities/` | Various | Product-level σ estimates; relevant for ACR welfare counterfactuals. |
| **FAF5 Highway Assignment** | `/Volumes/HELFRICH-GD/FAF5_Highway_Assignment_Results.zip` | 2017, 2022 | **Link-level**, not OD. Not directly usable for state-pair regressions. To get FAF5 state-pair OD, download FAF5 Regional Database from BTS separately. |
| **Liz's TradeFlows.dta** | `~/Library/CloudStorage/OneDrive-Personal/Research/GravityModelThesis/...` | Country-pair | Referenced in Liz's notebook. Exact location TBD — Liz can confirm. |

## 2. Bilateral gravity covariates

| Dataset | Path | Coverage | Notes |
|---|---|---|---|
| **CEPII Gravity V202211** | `/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv` | 1948–2021 bilateral | GDP, pop, RTA, language, contiguity, colony, legal-system, plus distance variants (`dist`, `distcap`, `distw_arithmetic`, `distw_harmonic`). This is the standard auxiliary dataset for gravity regressions. |
| **CEPII Countries V202211** | same dir / `Countries_V202211.csv` | ~250 countries | ISO3 mapping, capital coords, lat/lon, official languages, legal systems. |

## 3. Geospatial inputs — population

| Dataset | Path | Coverage | Notes |
|---|---|---|---|
| **WorldPop 1km UNadj (BRA, COL, IRL, PRT)** | `/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/worldpop/{bra,col,irl,prt}_ppp_YYYY_1km_UNadj.tif` | Annual 2001+ for these four countries | 80 files. Pulled for Katia's blended-finance project. |
| **WorldPop UK 2020 100m** | `/Volumes/HELFRICH-GD/UK_EconomicData/data/raw/WorldPop_UK_2020_100m.tif` | UK only, 2020 only, **100m resolution** | Higher-resolution sample; useful for UK sub-national distance robustness. |
| **WorldPop global (not local)** | (WorldPop FTP) | Annual, 1km, ~50 GB per year | Pull country-by-country as needed. |
| **GHS-POP 2020 R2023A 1km** | `/Volumes/HELFRICH-GD/TEG_data/inputs/ghsl/pop/GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif` | Global, 2020 only locally | JRC GHSL series. Full epochs 1975–2030 (5-yr) downloadable from JRC. ~2 GB per epoch at 1km global. |
| **GHS-SMOD 2020 R2023A 1km** | `/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/ghs_smod/` and `/Volumes/HELFRICH-GD/TEG_data/inputs/ghsl/GHS_SMOD_…2020….tif` | Global, 2020 | Settlement Model classification (urban centre, dense urban cluster, semi-dense, rural). Useful as urban/rural overlay if we want to test urban-vs-rural distance variants. |

## 4. Geospatial inputs — nightlights

| Dataset | Path | Coverage | Notes |
|---|---|---|---|
| **Harmonized DMSP/VIIRS (Li et al. 2020)** | `/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/nightlights/Harmonized_DN_NTL_*` | Annual global, **2010–2019** locally; Li 2020 product covers 1992–2018 originally with continuation extensions | The right product for cross-sensor consistency. `_calDMSP` files are calibrated DMSP for 2010-2013 overlap; `_simVIIRS` files are simulated-VIIRS-like for 2014-2019. |
| **VIIRS DNB raw** | (empty dir `/Volumes/HELFRICH-GD/viirs_data/`) | — | Placeholder; raw VIIRS DNB can be pulled from NOAA at month or annual composite level. |

## 5. Boundaries

| Dataset | Path | Coverage | Notes |
|---|---|---|---|
| Natural Earth admin0 | not yet local | Country polygons | Standard. Pull from `naturalearthdata.com`. |
| US Census state CB files | not yet local | US state polygons (50 + DC + territories) | Pull from `census.gov/geographies/mapping-files/`. |
| GADM | not yet local | Sub-national polygons globally | Optional for sub-national robustness. |

## 6. What's missing and needs to be pulled

In priority order for the v0.9 panel:

1. **WorldPop or GHS-POP for years other than 2020.** GHS-POP at 1km has 1975–2030 in 5-year epochs from JRC. For annual coverage, WorldPop pulled country-by-country from their FTP. Estimate: 1–2 days of pulling for the BACI country set.
2. **Harmonized nightlights for years not in 2010–2019.** The Li 2020 product has been updated through 2022 by Li et al. and others (Chen et al. 2021 extended; NOAA VIIRS-only 2020+). Estimate: half a day.
3. **CEPII GeoDist standalone.** Often more convenient than `Gravity_V202211` for distance-only use cases.
4. **FAF5 Regional Database** (state-pair OD) if Liz wants US-state-pair gravity. The CSV form is on BTS. Estimate: hour.
5. **Natural Earth admin0 polygons.** Trivial pull.

## 7. Software inventory

| Library | Status | Notes |
|---|---|---|
| `pyfixest` | ✓ in Paper 5 venv | PPML with hdfe FE absorption |
| `POT` (Python Optimal Transport) | ✓ just installed in Paper 5 venv (0.9.6) | EMD, Sinkhorn, sliced Wasserstein |
| `rasterio` | ✓ system Python | Raster I/O |
| `geopandas` | ✓ system Python | Vector boundaries |
| `scipy` | ✓ system Python 1.13 | rank stats, linalg |
| `statsmodels` | ✓ | OLS, GLM fallbacks |
| `linearmodels` | ✓ (per CLAUDE.md) | PanelOLS, IV |
| `numba` | ✓ in Paper 5 venv (req for pyfixest) | JIT for FE absorption |
| `pyfixest` (system Python) | ✗ not installed | We deliberately use the Paper 5 venv only |

## 8. Compute footprint estimates

| Pipeline | Spatial unit | Atoms per region | Method | Pairs | Wall time (rough) |
|---|---|---|---|---|---|
| Centroid-based (M1–M3, M5) | country (1km raster) | ~5e4–1e6 | spherical-mean | 200²/2 = 20,000 | < 1 hour |
| Exact W1 (M6, M7) | country (1km raster) | ~5e4 (subsampled) | `ot.emd2` | 20,000 | ~24 hours single-threaded; ~3 hours with rayon-style parallel |
| Sinkhorn (M8) | country (1km raster) | ~1e5 | `ot.sinkhorn2` ε=50 | 20,000 | ~6 hours |
| Sliced W (M9) | country (full atoms) | full atoms | random projections | 20,000 | < 1 hour |
| **US states only** | state (1km raster) | ~1e4 each | exact W1 | 2,550 | < 30 min |

The US state-pair panel is small enough that we can run all measures with exact OT in well under an hour total. The country-pair panel needs the Sinkhorn-or-sliced approximations to be tractable on a single machine.

---

## How to update this document

- When a new dataset is downloaded or located, add a row to the appropriate table.
- When a path moves, update — don't add. There should be exactly one canonical path per dataset.
- When something is downloaded that's not on this list, add it.
- Estimated compute footprints get refined once we have real benchmarks from running.
