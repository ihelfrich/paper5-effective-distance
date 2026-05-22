# Pre-Analysis Plan — Effective-Distance Horse Race

**Authors:** Elizaveta Gonchar, Ian T. Helfrich
**Drafted:** 2026-05-20
**Revised:** 2026-05-20 (afternoon) — corrected the measurement framework after verification + literature triangulation
**Status:** v0.2 — to be reviewed with EG before any analysis is run

This document specifies the experimental design *before* the analysis is run.
Its purpose is discipline: we commit to the metrics and decision rules in
advance so that the eventual write-up cannot be a post-hoc rationalization of
whatever results we happen to get.

---

## 1. Research question

> Among bilateral distance measures grounded in a spatially-disaggregated CES
> gravity framework, **which formulation captures the most economically
> meaningful variation in trade flows**, and **does an optimal-transport
> formulation produce empirically measurable improvements over centroid-based
> approximations**?

This question subdivides into three falsifiable predictions, restated formally
from the commented-out theory section of the SSRN draft.

**P1 (Asymmetry).** A directional, sub-nationally weighted distance measure
will explain residual asymmetries in bilateral trade flows that a symmetric
measure cannot. Operationally: in a regression of $\log |X_{ij} - X_{ji}|$ on
$|d_{ij} - d_{ji}|$, the coefficient on the directional measure is positive
and statistically significant, conditional on standard gravity controls.

**P2 (Internal distance).** A measure that handles intra-national distances by
the same mechanism as inter-national distances reduces the estimated *border
effect* relative to a measure that uses ad-hoc internal-distance proxies.

**P3 (Distance elasticity).** The estimated distance elasticity is larger in
magnitude for measures that incorporate sub-national heterogeneity than for
symmetric capital-to-capital or unweighted-centroid measures.

The OT-vs-centroid comparison is operationalized as a fourth prediction.

**P4 (OT improvement).** A Wasserstein-distance measure with economically
meaningful ground cost (great-circle distance raised to the trade-cost
elasticity power) explains additional variation in trade flows beyond a
centroid-based directional measure. Operationally: the within-$R^2$ from a
PPML gravity regression using the OT measure exceeds the centroid-measure
within-$R^2$ by at least 0.005, and the elasticity difference is significant
at $p < 0.05$ in a Davidson–MacKinnon $J$-test.

---

## 2. Datasets

### 2.1 Bilateral flow data

We use two independent flow datasets to test the same hypotheses on different
spatial scales. Each is a separate test, not a pooled sample.

**Country-pair panel.**
- **Source.** BACI HS02 V202401b (CEPII), aggregated from 6-digit product to
  bilateral country-year totals (USD).
- **Coverage.** 2002 onward; for v0.9 of the paper we use 2002–2022 (21 years).
- **Units.** ISO3 country codes.
- **Local path.** `/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/`.

**US state-pair panel** *(if available from Liz's working data).*
- **Source.** Liz's `TradeFlows.dta` augmented with state-level flows, OR
  FAF5 Regional Database (BTS) for FAF-zone-pair commodity flows.
- **Spatial units.** 50 states + DC (51 spatial units).
- **Time.** Whatever Liz has staged.
- **Note.** The FAF5 highway-assignment zip on `HELFRICH-GD` is link-level,
  not OD. The OD table is the FAF5 Regional Database — separate download.

### 2.2 Gravity covariates

- **CEPII Gravity V202211** (`Gravity_csv_V202211/Gravity_V202211.csv`):
  GDP, population, RTA, common language, contiguity, colonial ties, common
  legal system, distance variants (`dist`, `distcap`, `distw_arithmetic`,
  `distw_harmonic`).

### 2.3 Geospatial inputs for distance computation

**Population rasters.**
- **WorldPop UN-adjusted, 1km annual** (already on `HELFRICH-GD` for selected
  countries; will need to be pulled from WorldPop FTP for full global coverage).
- **GHS-POP 2020 R2023A 1km** (`GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif`)
  as a consistent global alternative. GHS-POP exists every 5 years 1975–2030
  in the JRC release; we will pull additional epochs for time-variation
  robustness.

**Nightlights rasters.**
- **Harmonized DMSP/VIIRS, 1km annual, Li et al. 2020** (already on
  `HELFRICH-GD/KatiaBlendedFinance/raster_cache/nightlights/`):
  2010–2019 in the cached set; will extend.

**Boundaries.**
- Natural Earth Admin0 for countries.
- US Census Cartographic Boundary Files for states.
- GADM for sub-national robustness (where applicable).

---

## 3. Distance measures to be tested

**REVISION 2026-05-20:** After verification by an independent methodology agent
(see SESSION_NOTES_2026-05-20.md), the correct framing of our contribution is
NOT optimal-transport. It is the **CES-aggregated effective distance** of
Head & Mayer (2002, 2014), extended from city-level to raster resolution.

The structurally-correct formula derived from CES iceberg-cost gravity is the
**generalized mean of pairwise distances** of order $\theta$:

$$d_{ij}^{\text{eff}}(\theta) = \left( \int\int f_i(x) g_j(y) d(x,y)^{\theta}\, dx\, dy \right)^{1/\theta}$$

where $\theta = \delta(1-\sigma)$ — for the iceberg distance-elasticity $\delta \approx 1$
and trade-cost elasticity $\sigma$, this is $1-\sigma$. For typical
$\sigma \in [4, 8]$, $\theta \in [-7, -3]$ — concave in distance, dominated
by short within-country pairs. The arithmetic mean ($\theta = 1$) is structurally
inconsistent; CEPII's `distw_harmonic` uses $\theta = -1$ as a reduced-form
choice from typical empirical estimates.

We test the following measures (replacing the earlier OT-flavored list):

We test **seven** measures. Each is a function of two regions
$i, j$ and (optionally) time $t$. All distances are in great-circle kilometers
(or in $\text{km}^\theta$ for iceberg variants).

### 3.1 Baseline (existing literature)
| ID    | Measure                       | Symmetric? | Time-varying? | Reference                            |
|-------|-------------------------------|:----------:|:-------------:|--------------------------------------|
| `M1`  | Capital-to-capital            | yes        | no            | Aitken (1973)                        |
| `M2`  | Unweighted centroid           | yes        | no            | standard                             |
| `M3`  | Pop-weighted centroid (CEPII) | yes        | no            | Mayer–Zignago (2011) `distw_arith`   |
| `M4`  | Hinz (2017) random sampling   | yes        | yes           | Hinz (2017)                          |

### 3.2 Contribution (existing draft)
| ID    | Measure                                                      | Symmetric? | Time-varying? |
|-------|--------------------------------------------------------------|:----------:|:-------------:|
| `M5`  | $d_{ij}^{NL}$: VIIRS-origin centroid → LandScan-dest centroid |    no      |     yes       |

### 3.3 Contribution (new, CES-aggregated raster resolution)
| ID    | Measure                                                                                          | Symmetric? | Time-varying? |
|-------|--------------------------------------------------------------------------------------------------|:----------:|:-------------:|
| `M6`  | $d^{\text{eff}}_{ij}(\theta=-1)$: raster CES, harmonic-mean (matches CEPII convention but at raster res.) | yes (sym weights) | yes |
| `M7`  | $d^{\text{eff}}_{ij}(\theta=-3)$: raster CES at $\sigma=4$ (low elasticity)                           | yes (sym weights) | yes |
| `M8`  | $d^{\text{eff}}_{ij}(\theta=-5)$: raster CES at $\sigma=6$ (mid; preferred baseline)                  | yes (sym weights) | yes |
| `M9`  | $d^{\text{eff}}_{ij}(\theta=-5, f=\text{VIIRS}, g=\text{LandScan})$: directional (origin VIIRS × destination LandScan) | **no** | yes |

Measures `M5`–`M9` form the contribution surface. `M5` is the existing draft's
contribution. `M6`–`M9` are the new measures whose performance is to be tested.
`M9` is an efficiency variant included only to validate the compute pipeline
against `M6`/`M7` (it should produce similar rank correlations on a held-out
subset; if it does, future work can use `M9` for scale).

---

## 4. Estimation specifications

### 4.1 Primary specification (preferred for headline results)

Three-way fixed-effects PPML on the country-pair panel:

$$
X_{ijt} = \exp\!\left[\,\alpha_{it} + \beta_{jt} + \gamma_{ij}^{(c)} + \delta \ln d_{ijt}^{(M)} \,\right] + \varepsilon_{ijt}
$$

where $X_{ijt}$ is bilateral trade flow (USD), $\alpha_{it}$ is
origin-by-year FE absorbing exporter multilateral resistance, $\beta_{jt}$
is destination-by-year FE absorbing importer MR, and $\gamma_{ij}^{(c)}$ is a
country-pair FE absorbing all time-invariant bilateral characteristics. The
parameter of interest is $\delta$, the distance elasticity.

When a distance measure is time-invariant (`M1`–`M3`), the pair FE absorbs
its level and $\delta$ is not identified; in that case we drop the pair FE
and add the standard time-invariant gravity controls (`contig`, `comlang_off`,
`colony`, `comlang_ethno`, `comleg_pretrans`, `rta`).

### 4.2 Robustness specifications

- **OLS log-gravity** (Silva-Tenreyro 2006 nonwithstanding) as a sanity check
- **PPML without pair FE** but with full bilateral controls
- **Gamma GLM** as a heteroscedasticity alternative
- **OLS-PPML horse race** via Davidson–MacKinnon $J$-test
- **Split sample** by income group (high-high vs high-low vs low-low)
- **Split sample** by region pair

### 4.3 Inference

Standard errors clustered at the pair level. Bootstrap standard errors with
1000 replicates for the difference-in-elasticity tests across distance
measures. Multiple-comparisons correction: Holm–Bonferroni for the 9-measure
family of tests on each metric.

---

## 5. Metrics

For each measure $M_k$ in $\{M1, ..., M9\}$, we report:

| Metric                                  | Definition                                                  | Decision direction |
|----------------------------------------|-------------------------------------------------------------|--------------------|
| Distance elasticity $\hat{\delta}_k$    | PPML coefficient with cluster-robust SE                     | reported           |
| Within-$R^2$                            | After absorbing FEs                                          | higher is better   |
| AIC, BIC                                | Information criteria                                         | lower is better    |
| OOS RMSE                                | Leave-one-year-out CV, then leave-one-country-out CV         | lower is better    |
| Davidson–MacKinnon $J$-test             | Pairwise tests between measures                              | reported           |
| Rank correlation matrix                 | Spearman $\rho$ across the 9 measures on bilateral distances | exploratory        |
| Residual asymmetry $R^2_{\text{asym}}$  | $R^2$ from regressing $\log|X_{ij}/X_{ji}|$ on $\log(d_{ij}/d_{ji})$ | higher is better (P1) |
| ACR welfare-counterfactual spread       | Range of welfare gains across measures under symmetric 10% tariff shock | reported           |

### 5.1 Decision rules

**P1 (Asymmetry).** Confirmed if a directional measure (`M5`–`M9`) produces
$R^2_{\text{asym}} \geq 0.05$ and the coefficient is positive at $p < 0.05$.

**P2 (Border effect).** Confirmed if the estimated home-bias coefficient
$\hat{\zeta}$ from a specification including internal trade is at least 25%
smaller in magnitude under `M5`–`M9` than under `M3`.

**P3 (Elasticity magnitude).** Confirmed if $|\hat{\delta}|$ under `M5`–`M9`
exceeds $|\hat{\delta}|$ under `M3` by at least 0.10 elasticity-unit and the
difference is significant at $p < 0.05$.

**P4 (OT improvement).** Confirmed if at least one of `M6`–`M9` produces:
(a) within-$R^2$ at least 0.005 higher than `M5`, AND
(b) Davidson–MacKinnon $J$-test rejects `M5` in favor of the OT measure at
$p < 0.05$.

If `P4` is not confirmed, the paper's contribution is the directional centroid
measure `M5`, not the OT measure. If `P4` is confirmed, the OT measure becomes
the headline result and `M5` is reframed as a fast approximation.

---

## 6. Identification & threats to validity

### 6.1 Endogeneity
- The new distance measures are constructed from satellite data, not from
  trade flow data; we do not estimate them jointly with trade flows.
- We will compute reverse-causality robustness: are the centroids "moving
  toward" trading partners (a Tobler / market-access effect)? If yes, lag the
  centroid by $\geq 5$ years.

### 6.2 Data quality
- VIIRS and DMSP have known issues: blooming, top-coding (DMSP), small-fire
  contamination. We use the Li et al. (2020) harmonized series specifically
  because it addresses cross-sensor inconsistencies. We will run robustness
  with raw VIIRS-DNB to confirm conclusions are not driven by harmonization.
- WorldPop and LandScan disagree at the sub-national level in some countries.
  We will run the analysis with each independently and report robustness.

### 6.3 Multiple testing
- Holm–Bonferroni correction within each family of tests.
- We report results for all 9 measures even when only the headline measure
  is reported in the main text. Robustness tables show every measure.

### 6.4 Power
- With $\sim 200$ countries and 21 years of bilateral flows, $N \approx 200^2 \cdot 21 \approx 840{,}000$
  observations (before dropping zeros and missing). At this sample size, an
  elasticity difference of 0.10 is detectable at conventional power for any
  reasonable variance structure.
- For the US state panel: $N \approx 51^2 \cdot T$. With $T=10$, $N \approx 26{,}000$
  — still adequate for elasticity comparisons but smaller than the country panel.

---

## 7. Replication & openness

- All code in `Paper5_EffectiveDistance/` with `pyproject.toml` and pinned
  dependencies.
- All intermediate datasets in `data/derived/` with provenance hashes.
- Final bilateral distance panel released on Harvard Dataverse, CC-BY 4.0,
  with all 9 measures + metadata.
- Pre-analysis plan committed to the repository before the first horse-race
  regression is run.

---

## 8. Out of scope (deliberately deferred)

These extensions are interesting and may justify follow-up papers but are not
in the v0.9 ship target:

- The multidimensional Riemannian cost manifold with endogenous Ricci-driven
  evolution (commented in the SSRN draft; theoretically rich, scope-mismatched
  for a gravity paper).
- Endogenous productivity growth and labor mobility dynamics.
- Three-dimensional cost (e.g., topographic elevation, infrastructure quality)
  beyond the two-dimensional great-circle baseline.
- Stochastic-OT formulations.

---

## 9. Authorship & roles

- **Elizaveta Gonchar:** lead empiricist; US panel construction; gravity
  estimation lead; manuscript writing lead.
- **Ian T. Helfrich:** OT theoretical derivation; distance-measure
  implementation; computational infrastructure; co-writer on theory section.

---

## 10. Open questions for EG before execution

1. Which trade-flow panel are you actively working with — country-level
   updated `TradeFlows.dta`, US-state-augmented version, FAF5 OD, or other?
2. Do you have a preferred PPML implementation (Stata `ppmlhdfe`, `pyfixest`,
   R `fixest`)? Different defaults can produce different SEs.
3. Are you committed to harmonized DMSP/VIIRS (Li 2020) or do you want to
   use VIIRS-only for post-2012 to avoid harmonization concerns?
4. Time-coverage preference for v0.9: 2002–2022 (BACI full range) vs
   2010–2022 (where Li 2020 harmonized lights start)?
5. Any prior on the OT cost choice ($W_1$ vs $W_p^\theta$)?
