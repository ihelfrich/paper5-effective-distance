# Replication map for Paper 5a

Every table and figure in `manuscript/paper_5a/paper_5a.tex` traces to a
script and a saved output. The five CSVs in `replication/data/` are the
canonical results files; they were produced on the dates embedded in git
history and copied here from the compute archive on 2026-06-10 so the
repository reproduces the paper without the external drive.

## Tables

| Paper object | Producing script | Saved output |
|---|---|---|
| Table 1 (pole values) | direct evaluation, see Appendix B | analytic; verified 2026-06-10 |
| Panel PPML baseline + sweep tables | `notebooks/18f_spec_b_prime_panel_ppml.py` | `replication/data/spec_b_prime_panel_ppml.csv` |
| Cross-section OLS table | `notebooks/18d_spec_b_prime_fixed.py` | `replication/data/spec_b_prime_fixed.csv` |
| Cross-section PPML robustness (in text) | `notebooks/18e_spec_b_prime_ppml.py` | `replication/data/spec_b_prime_ppml.csv` |
| USA worked example (Section 4.1) | `country_d_eff()` in `notebooks/18f_spec_b_prime_panel_ppml.py` | `replication/data/spec_b_prime_country_d_eff.csv` (theta = -5 columns; the theta = -3 and -7 columns regenerate in seconds from the same function) |

## Figures

| Figure | Generator | Input data |
|---|---|---|
| Fig 1 (regularization curve) | `notebooks/paper5a_figures.py` | `replication/data/regularization_dmin_cutoff.csv`, produced by `notebooks/19c_regularization_sensitivity.py` (3 seeds, 10,000 atoms) |
| Fig 2 (panel home bias vs d_min) | `notebooks/paper5a_figures_panel.py` | `replication/data/spec_b_prime_panel_ppml.csv` |

## Raw inputs (not redistributed; see DATA_LICENSE.md)

- BACI HS02 V202401b (CEPII)
- CEPII Gravity V202211
- GHS-POP R2023A 1 km rasters (JRC)

## Known reproducibility caveats

1. The standard-error columns in `spec_b_prime_panel_ppml.csv` are
   unscaled IRLS artifacts (order 1e-07) and should not be used; the
   manuscript quotes no SEs from this file. Re-estimation with correct
   covariance is on the revision list.
2. The Wei (1996) internal-flow proxy in the paper text (USA 2010:
   $13.1T) reflects the corrected GDP-minus-goods-exports computation;
   the archived `*_oydy_input_*.csv` panel input predates that fix and
   carries GDP itself as the USA intra flow. The 18f pipeline builds its
   panel in memory; a regenerated input dump is on the revision list.
3. No execution log for 18f exists in `logs/`; rerun and save one before
   journal submission.
