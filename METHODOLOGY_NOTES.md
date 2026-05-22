# Methodology notes — known issues and design choices

*Last updated: 2026-05-20 evening (autonomous run)*

This document tracks methodological subtleties that surfaced while building
the empirical pipeline. Each entry: the issue, the current treatment, and
what would close the issue rigorously if time allowed.

## 1. Top-N-by-population atom subsampling biases d_eff(θ < 0) downward

**The issue.** Notebook 13 caps each country's atom cloud at `max_atoms = 5000`,
selected as the top-5000 highest-population cells. For large dense countries
(USA, China) this means we keep mostly metropolitan cells and discard the
spatially-extensive rural tail.

For θ ≥ 0 (arithmetic / geometric / sub-additive means) this is fine — the
discarded cells contribute small mass × moderate distance, so they barely
move the mean.

For θ < 0 (CES-consistent, gravity-relevant range), short distances dominate
the integrand. Keeping only metropolitan cells means we keep pairs of
near-neighbors (sub-10-km pairs within Houston, within NYC, within LA, etc.)
that pull d_eff down very fast. Intra-metro distances of 1–10 km, raised to
the θ = -5 power, swamp the inter-metro pairs.

**Empirical evidence.** First USA run with max_atoms=3000 returned
`d_eff(θ=-5) / d_HM = 0.11`, implying d_eff(θ=-5) ≈ 130 km. Raising to
max_atoms=5000 gave `ratio = 0.00` (rounded), implying d_eff(θ=-5) ≈ a few
km. Direction of bias is what theory predicts: more atoms in dense clusters
→ more short-distance pairs → smaller harmonic-style mean.

**Current treatment.** Use `max_atoms = 5000` as a placeholder.
**Document numerical sensitivity in notebook 16 (`sampling_bias_diagnostic.py`)**:
for the four small countries (NLD, BEL, CHE, DNK) where full-resolution is
tractable, compute d_eff(θ=-5) under three sampling rules — top-N, stratified
random by weight, and full-resolution. Report relative error.

**Resolution candidates** (in order of methodological preference):
  - **Block coarsening.** Down-sample the raster from 1 km to 5 km or 10 km
    by *summing* (not averaging) population over 5×5 or 10×10 blocks. Preserves
    total mass *and* spatial coverage. The downside: discards fine within-city
    structure that θ=-5 weights heavily. This is probably the right answer for
    the panel.
  - **Stratified weight-proportional random sampling.** Sample N cells with
    probability ∝ weight (not just the top-N). Preserves spatial distribution
    on average but introduces Monte Carlo variance.
  - **Hierarchical: full inside metros, coarsened outside.** Build a quadtree
    that adapts resolution to local density. Most accurate, most code.

The pre-analysis plan should commit to block coarsening at 5 km for the main
panel and report stratified-sampled robustness checks.

## 2. Geographic area vs populated area for Head-Mayer baseline

**The issue.** Head-Mayer's closed form `d_ii = 0.67 sqrt(area_i/π)` uses the
*country's total geographic area* — the conceptual analogue is a uniform-density
disk. For Russia, the geographic area is ~17M km², while the area of *populated*
GHS-POP cells is roughly 0.8M km² (most of Siberia is empty).

If we use populated area, d_HM falls dramatically and the ratio
`d_eff / d_HM` flatters our raster measure artificially.

**Current treatment.** Notebook 13 v3 computes geographic country area from the
Natural Earth polygon reprojected to Mollweide (`ESRI:54009`) and uses that for
d_HM. Both `area_km2` (geographic) and `pop_area_km2` (raster-populated)
are saved for transparency.

## 3. Production vs consumption asymmetry: the choice of nightlights proxy

**The issue.** Henderson, Storeygard & Weil (2012) show nightlights correlate
strongly with GDP at the country level but the *spatial distribution* of light
emission is dominated by industrial activity, road networks, gas flares
(Bluhm & Krause 2018), not by residential population per se. Using VIIRS
nightlights as the "production proxy" is defensible — it picks up factories,
ports, mining, refineries — but it also picks up commercial districts
(consumption activity).

**Current treatment.** Notebook 15 uses VIIRS as origin and GHS-POP as
destination. The asymmetry index ASI = log(d_ij / d_ji) is interpreted as a
*directional shift in effective distance* arising from the spatial divergence
of activity proxies, not strictly "production minus consumption."

**Resolution.** Sensitivity check using:
  - VIIRS × LandScan (alt consumption proxy)
  - VIIRS / pop ratio (a "specialization" index that strips out
    residential confound)
  - SRTM-corrected nightlights (Bluhm-Krause's "EconLight" — removes
    gas flares and offshore platforms)

## 4. Earth curvature and cosine-of-latitude correction

**The issue.** For countries spanning > 1000 km, planar (lat × 111.32 + lon ×
111.32 × cos) approximations break down. Russia spans ~9000 km east-west;
Canada spans ~5500 km; USA contiguous spans ~4500 km.

**Current treatment.** All pairwise distances use haversine great-circle in
`ground_cost_matrix(kind="gc_km")`. The cosine correction is applied for
equal-angle rasters (EPSG:4326), not for Mollweide (equal-area).

**Verified** in `tests/test_ot_distance.py::test_haversine_against_known_pairs`.

## 5. Intra-national trade flow proxy

**The issue.** We don't observe X_ii directly. The standard proxy is
`X_ii = max(0, GDP_i - sum_j X_ij)` — Wei (1996), Anderson-van Wincoop (2003).
This conflates true intra-national trade with services, government, and
non-tradeables.

**Current treatment.** Standard proxy. Notebook 14 acknowledges this as a
caveat; the comparative-statics (Spec A vs Spec B) holds the proxy fixed and
only swaps d_ii.

**Resolution.** A subsample using OECD ICIO inter-regional flows would give
true X_ii for ~40 OECD countries. Worth pursuing if Spec B shows a meaningful
home-bias shift.

## 6. Symmetry of the CES integral

**The issue.** The CES generalized mean of pairwise distances is symmetric
in (mu_origin, nu_dest) when both come from the same density (since
`∫∫ f(x)f(y) d(x,y)^θ` is symmetric under x↔y). For asymmetric origin-dest
densities (VIIRS_i × pop_j), it is NOT symmetric.

**Verified** in `tests/test_directional_and_velocity.py::test_asymmetric_when_densities_diverge`.

This is the empirical content of prediction P1: ASI_ij = log(d_ij / d_ji)
is non-zero in general.

## 7. Time-varying centroids: the Wasserstein-velocity prediction

**The issue.** P2 says year-over-year centroid shifts predict bilateral
trade growth. For this to identify, we need:
  - Nightlights or population time-series with at least 10–15 years
  - Plausibly exogenous centroid drift (industrial relocations,
    interior development) not driven directly by trade with a single partner

**Current treatment.** Use harmonized DMSP/VIIRS 1992–2024 from Li et al.
(2020). Identification via origin × year + destination × year + pair FE
(Yotov gravity), which absorbs any country-specific shock or pair-level
time-invariant unobservable. The remaining variation is *how the within-
country mass shift projects onto the i↔j direction*, year by year.

**Concern.** A trade-induced agglomeration response (Krugman-Venables)
could move mass toward the partner, producing reverse causation. We instrument
with weather / climate shocks driving structural shifts (severe drought
moving Brazilian Northeast population to the Southeast), or with
geological shocks (mining-belt openings).

---

## Status of empirical pipeline

| Notebook | Purpose | Status |
|---|---|---|
| 12 | Synthetic θ-divergence simulation | ✓ done |
| 13 | Real-country d_ii^eff at 5 thetas (27 countries) | RUNNING (v3, geo-area fix) |
| 14 | Border-effect pilot (Spec A vs Spec B) | Wired, waiting on 13 |
| 15 | Directional asymmetry pilot (VIIRS × pop) | Written, untested on real data |
| 16 | Sampling-bias diagnostic | Written, untested on real data |

| Module | Tests | Status |
|---|---|---|
| `ot_distance.ces_effective_distance` | 4 unit | ✓ pass |
| `region_raster.mask_raster_to_country` | (integration) | ✓ runs |
| `directional_asymmetry.*` | 4 unit | ✓ pass |
| `wasserstein_velocity.*` | 3 unit | ✓ pass |
| `gravity_panel.build_gravity_panel` | (integration) | ✓ runs in nb14 |
