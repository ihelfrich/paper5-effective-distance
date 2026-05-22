# Autonomous block deliverables — 2026-05-21

Single-page index for everything produced during the multi-hour autonomous
work block. All artifacts are also synced to Google Drive at
`AA_IanHelfrich_Rsrch_GTPHD/Paper5_EffectiveDistance_Theory/`.

---

## Paper 5a: the SSRN preprint

| Version | Pages | What's new |
|---|---|---|
| v0.1 (~19:09) | 9 | Math complete, Figure 1, empirical placeholder |
| v0.2 (~19:30) | 12 | Full empirical numbers from nb18d, Figure 2, three tables |
| v0.3 (~19:39) | 13 | Grounded intro opening with HM 2002 line 502 quote |
| v0.4 (~19:54) | 13 | 34-entry verified bibliography, strengthened lit positioning vs. EK/CP/CDK |
| **v0.5 (~20:12, current)** | **13** | **PPML preliminary results subsection added; abstract polished with concrete numbers** |

**Primary file**: `manuscript/paper_5a/paper_5a.pdf`
**LaTeX source**: `manuscript/paper_5a/paper_5a.tex`
**Bibliography**: `manuscript/paper_5a/paper_5a_refs.bib` (34 verified entries)

## Companion documents (Paper 5a)

| File | Description |
|---|---|
| `manuscript/paper_5a/visual_abstract.pdf` | 3-page landscape executive summary with both figures and the takeaway box |
| `manuscript/paper_5a/referee_responses.md` | 8-critique pre-emptive defense document |
| `manuscript/paper_5a/scoop_protection_email_draft.md` | v1 (Mayer) + v2 (Yotov) email drafts |
| `manuscript/paper_5a/SSRN_METADATA.md` | SSRN posting metadata block (title, abstract, JEL codes, keywords) |
| `manuscript/paper_5a/plain_english_summary.md` | Non-technical 2-page summary for policy / journalism audiences |
| `manuscript/paper_5a/paper_5a_refs_audit.md` | Bibliography verification log |
| `manuscript/paper_5a/figures/fig1_regularization_curve.{pdf,png}` | Regularization curve on uniform disk |
| `manuscript/paper_5a/figures/fig2_home_bias_vs_dmin.{pdf,png}` | Home-bias coefficient vs d_min, by year and θ |

## Teaching guide

| File | Pages | What's in it |
|---|---|---|
| `learning_guide/master_guide.pdf` (Volume 1) | 35 | Trade-econ foundations, CES gravity, internal distance, generalized means and the divergence (Parts I-IV) |
| `learning_guide/master_guide_vol2.pdf` (Volume 2) | **50** | Numerical evaluation, real-data construction, regression specification, empirical results, forward path (Parts V-IX). Fully developed chapters: Ch 1 (discretization), Ch 2 (uniform disk sweep), Ch 5 (BACI), Ch 6 (CEPII + Wei proxy), Ch 7 (population rasters), Ch 9 (OLS gravity with FE), Ch 10 (the collinearity trap), Ch 11 (PPML). Remaining 12 chapters are skeleton stubs. |

## Empirical apparatus

| Notebook / module | Description |
|---|---|
| `notebooks/18d_spec_b_prime_fixed.py` | Bug-fixed gravity regression with the formula API; sweep across years × θ × d_min |
| `notebooks/18e_spec_b_prime_ppml.py` | PPML robustness check; year 2010 only |
| `notebooks/paper5a_figures.py` | Figure generators for fig1 and fig2 |
| `notebooks/19c_regularization_sensitivity.py` | Uniform-disk regularization sweep (from earlier block) |
| `data_derived/spec_b_prime_fixed.csv` | 96-row sweep output |
| `data_derived/spec_b_prime_ppml.csv` | PPML sweep output (partial as of journaling) |
| `data_derived/regularization_dmin_cutoff.csv` | Disk simulation curve |

## Literature intel

| File | Description |
|---|---|
| `literature_intel/novelty_assessment_2026-05-21.md` | Master novelty assessment with all sub-agent findings |
| `literature_intel/headmayer_2014_novelty_check.md` | Agent close-read of HM 2014 |
| `literature_intel/pole_novelty_triangulation.md` | Agent triangulation across 6 candidate literatures |
| `literature_intel/structural_estimation_theta_check.md` | Agent close-read of EK 2002 + CP 2015 + CDK 2012 — confirms the conflation is HM-specific |

## Pre-GitHub-Pages stub

| File | Description |
|---|---|
| `docs/index.html` | Single-page HTML version of the abstract with both figures, math via KaTeX, and download links |
| `docs/fig1_regularization_curve.png` | Figure 1 image |
| `docs/fig2_home_bias_vs_dmin.png` | Figure 2 image |
| `docs/paper_5a.pdf` | Copy of v0.5 PDF |
| `docs/visual_abstract.pdf` | Copy of visual abstract |
| `docs/master_guide_vol2.pdf` | Copy of Volume 2 guide |

These files are ready to publish via GitHub Pages once the
`ihelfrich/paper5-effective-distance` repo is initialized. Awaiting the
"push it" command.

## Tasks closed during the block

- #29 Reframe paper around "wrong θ" thesis (refined HM-specific framing)
- #30 Check EK/CP/CDK for θ-conflation (NOVEL: the conflation is HM-specific)
- #31 Spec B-prime gravity regression (full sweep complete)
- #35 Debug Spec B-prime baseline sign-flip (root cause: sm.OLS with collinear log-GDP regressors)
- #36 Close out novelty check (NOVEL verdict confirmed)
- #41 Polish Paper 5a intro (grounded opening from HM 2002 line 502)
- #42 Investigate n_intra anomaly (pragmatic close; documented in Section 5.3)
- #43 Expand bibliography to 25+ entries (34 verified entries)
- #44 Build Volume 2 of teaching guide (50 pages; 8 fully developed chapters)
- #45 Write referee response document (8 critiques)
- #46 Build visual one-page abstract (3-page executive summary)
- #47 PPML robustness check (preliminary results: convergence warnings, direction matches OLS)

## Tasks still pending

- #9 US states panel horse-race notebook (not started in this block)
- #23 Replicate Borchert-Larch-Shikher-Yotov 2022 on ITPD-E R3 (not started)
- #24 Replicate Besedeš-Kohl-Lake 2020 from Mendeley (not started)
- #38 Initialize GitHub repo (local scaffold ready; awaiting "push it")
- #39 Send scoop-protection email to Mayer (Gmail draft ready; awaiting Ian's verify+send)

## Headline empirical finding (year 2010, the cleanest sample)

Home-bias coefficient `γ̂` across structural θ and regularization grain `d_min`:

| d_min (km) | θ = -3 | θ = -5 | θ = -7 |
|---|---|---|---|
| 0.5 | -1.155 (0.32x) | -2.919 (0.05x) | -3.360 (0.03x) |
| 1.0 | -0.365 (0.69x) | -2.054 (0.13x) | -2.494 (0.08x) |
| 2.0 | +0.308 (1.36x) | -1.293 (0.27x) | -1.733 (0.18x) |
| 5.0 | +1.114 (3.05x) | -0.378 (0.68x) | -0.819 (0.44x) |
| 10 | +1.601 (4.96x) | +0.209 (1.23x) | -0.231 (0.79x) |
| 50 | +2.336 (10.3x) | +1.282 (3.60x) | +0.842 (2.32x) |
| 500 | +2.610 (13.6x) | +2.031 (7.62x) | +1.590 (4.91x) |

**CEPII baseline (θ = -1)**: γ̂ = +2.132, multiplier 8.43x — the canonical
McCallum-Anderson-van Wincoop "border puzzle" finding.

The sign-flip threshold shifts cleanly with θ: between d_min ∈ (1, 2) km at
θ = -3 and between d_min ∈ (10, 20) km at θ = -7. The McCallum-AvW
multiplier sits at the high-d_min end of the sweep and disappears (or
flips sign) at fine spatial resolution.

PPML robustness (preliminary, with convergence caveats): the qualitative
pattern is preserved. Both OLS and PPML show the home coefficient moving
more negative as d_min shrinks. Magnitudes differ (PPML's distance
elasticity is smaller, baseline home is more negative) but the direction
is the same.

---

*Last updated: 2026-05-21 evening. Sole-author block by Ian T. Helfrich,
assisted by Claude (autonomous mode).*
