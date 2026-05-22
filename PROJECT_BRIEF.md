# PROJECT BRIEF — Paper 5, Effective Distance

**One-pager, self-contained. For quick reference, coauthor handoff, and grant/job-talk framing.**

---

## Title
*Effective Distance: A Satellite-Calibrated Bilateral Trade-Cost Panel, 2000–2024*

## Authors
Elizaveta Gonchar · Ian T. Helfrich (Georgia Institute of Technology)

## One-sentence contribution
We build the first global, time-varying, satellite-calibrated, multi-modal bilateral effective-distance panel for 200 countries over 2000–2024, validate it in the structural gravity framework of Yotov et al. (2016) with Weidner-Zylkin (2021) bias correction, show it materially revises ACR welfare counterfactuals on the 2021 Suez closure / 2023 Panama drought / 2024 Houthi Red Sea crisis, and release it as CC-BY public-good infrastructure.

## Why now
1. CEPII's `distw` is static (2004 city populations); every gravity paper uses it despite documented time-variation in population and economic activity redistribution.
2. Three live chokepoint shocks (2021–2024) make the cost of ignoring time-varying transport-network geography empirically salient and press-legible.
3. Open satellite data (VIIRS, WorldPop, GHS-POP, HRSL) and open transport networks (OpenStreetMap) now make a global panel computationally feasible for a solo lab.
4. No competitor has released this panel (as of Apr 2026); IMF's WP 25/93 signals the space is warming. First-mover DOI matters.

## Precedent
Extends Gonchar & Helfrich (2025), *"Trade in the Spotlight: Enhancing Gravity Model Predictions with Nightlights and Population-Weighted Distance Measures,"* SSRN 5202676 — our April 2025 working paper introducing the OVDL (Origin-VIIRS-to-Destination-LandScan) measure. Paper 5 adds the transport-network layer, structural gravity validation, welfare counterfactuals, and public panel release.

## Data
BACI HS02 (CEPII) · Gravity V202211 (CEPII) · WorldPop Constrained 100m + GHS-POP R2023A + HRSL 30m · VIIRS V2 annual + Li et al. harmonized DMSP-VIIRS · OpenStreetMap planet · GSHHG · World Port Index · OpenFlights · CERDI Sea Distance · OECD TiVA · DESTA · GSDB · network centrality panel (Helfrich 2024 thesis, 2002–2022).

## Method (short)
Eight distance variants constructed, primary is `d_lcp_multi_t`:
$$d_{ij,t}^{\text{multi}} = \left[\sum_{k\in i,\ell\in j} w_{k,t}^{(i)} w_{\ell,t}^{(j)} \left(\min_{m} c_m^{\text{LCP}}(k,\ell;t)\right)^{\theta}\right]^{1/\theta}$$
with $k,\ell$ = agglomerations, $w$ = pop or VIIRS shares, $M$ = {road, maritime, air}, $c_m^{\text{LCP}}$ = least-cost path on OSM-derived transport graph with endogenous mode choice and arc-level shock toggling.

Gravity: three-way FE PPML (`pyfixest`), exporter-time × importer-time × pair, with intranational flows, Weidner-Zylkin bias correction, sector-heterogeneous θ_s robustness.

## Welfare counterfactuals
Suez 2021 (Ever Given, 7 days) · Panama Canal 2023 drought · Houthi Red Sea 2024 rerouting. ACR formula + Caliendo-Parro exact-hat algebra.

## Output
- Public panel: 200 countries × 25 years × 8 distance variants, ~40K rows/variant-year.
- Main paper: ~40 pages, target *JIE*.
- Data descriptor: ~2000 words, target *Scientific Data*.
- Software: Python (`paper5`), R mirror (`fixest`), Stata port (`ppmlhdfe`).

## Funding and affiliations
Georgia Institute of Technology. Seeking NBER affiliation for WP submission.

## Contact
Ian T. Helfrich — <ianthelfrich@gmail.com>
Elizaveta Gonchar — [Georgia Tech]
