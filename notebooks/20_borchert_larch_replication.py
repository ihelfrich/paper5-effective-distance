"""Replicate Borchert-Larch-Shikher-Yotov (2022) "Disaggregated Gravity" RIE.

Step 1: Load ITPD-E Release 3 and inspect its structure.
Step 2: Reproduce the headline three-way-FE PPML gravity for a selected
        industry (start with one sector, e.g. agriculture, to keep compute
        tractable while developing).
Step 3: Verify the headline coefficient estimates match the published paper.
Step 4: Repeat for all 170 industries.
Step 5 (next paper): Substitute the raster CES effective distance for CEPII
        centroid distance and re-estimate. Report which industries shift.

ITPD-E variables (Release 3, June 2025):
  exporter_iso3, importer_iso3, year, industry_id, industry_descr,
  trade (USD millions), flag_zero, flag_mirror, broad_sector
  - 264 countries (including non-sovereign units)
  - 170 industries: Agriculture (1-26), Mining/Energy (27-33), Manufacturing
    (34-153), Services (154-170)
  - 1986-2022 for goods (agriculture from 1986, others 1988-);
    2000-2022 for services
  - Includes domestic flows (i = j) — essential for Yotov-style identification

The structural gravity equation in BLSY is:
  X_{ij,t,k} = exp(α + θ log d_{ij} + β'Z_{ij} + π_{i,t} + π_{j,t}) + ε
where π_{i,t} are exporter-year FE (absorb output Y_i and outward MR Π_i),
π_{j,t} are importer-year FE (absorb expenditure E_j and inward MR P_j),
and Z_{ij} is the vector of standard gravity covariates (contig, comlang,
RTA, colony). Estimator: PPML (Poisson with log link, robust SE).

The trick: in PPML with three-way FE (exporter-year, importer-year, pair),
estimating directly is infeasible at the 264-country × 37-year scale.
Standard workaround: estimate by industry, with pair FE absorbed via
Mundlak / within transformation, or use the PPMLHDFE estimator's iterative
demeaning. We do not have PPMLHDFE in Python — use statsmodels.GLM with
formula-style C() and start with a small number of pair FE; expand once
the basic spec verifies.

Run:
  .venv/bin/python -u notebooks/20_borchert_larch_replication.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ITPDE = Path("/Volumes/HELFRICH-GD/TradeData/ITPDE/ITPDE_R03.csv")
CEPII = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")


def inspect_itpde_schema() -> None:
    """Read a small chunk to verify column names and types."""
    print(f"Inspecting {ITPDE.name} ({ITPDE.stat().st_size / 1e9:.1f} GB)...")
    chunk = pd.read_csv(ITPDE, nrows=1000)
    print(f"\nColumns ({len(chunk.columns)}):")
    for col in chunk.columns:
        print(f"  {col!r:30s}  dtype={chunk[col].dtype}  "
              f"n_unique={chunk[col].nunique()}  "
              f"example={chunk[col].iloc[0]!r}")
    print(f"\nFirst 3 rows:")
    print(chunk.head(3).to_string())


def load_single_industry(industry_id: int, *, year_min: int = 2010, year_max: int = 2019) -> pd.DataFrame:
    """Load all rows for one industry, one year window, into memory.

    For industry 1 (Wheat) at 2010-2019, ITPD-E has ~264*264*10 = ~700k rows,
    most of which are flag_zero=1 (no trade). PPML keeps the zeros.
    """
    print(f"\nLoading ITPD-E for industry_id = {industry_id}, "
          f"years {year_min}-{year_max}...")
    t0 = time.time()
    chunks = []
    for chunk in pd.read_csv(ITPDE, chunksize=2_000_000):
        sub = chunk[(chunk["industry_id"] == industry_id) &
                    (chunk["year"] >= year_min) & (chunk["year"] <= year_max)].copy()
        if len(sub) > 0:
            chunks.append(sub)
    if not chunks:
        return pd.DataFrame()
    df = pd.concat(chunks, ignore_index=True)
    print(f"  Loaded {len(df):,} rows in {time.time()-t0:.1f}s")
    return df


def load_cepii_panel(year_min: int, year_max: int) -> pd.DataFrame:
    """Load CEPII Gravity covariates for our years."""
    print(f"\nLoading CEPII Gravity {year_min}-{year_max}...")
    df = pd.read_csv(CEPII, dtype={"iso3_o": str, "iso3_d": str})
    df = df[(df["year"] >= year_min) & (df["year"] <= year_max)].copy()
    keep = ["year", "iso3_o", "iso3_d",
            "dist", "distcap", "distw_arithmetic", "distw_harmonic",
            "contig", "comlang_off", "comlang_ethno", "col_dep_ever",
            "gatt_o", "gatt_d", "rta"]
    keep = [c for c in keep if c in df.columns]
    return df[keep]


def merge_itpde_with_cepii(itpde: pd.DataFrame, cepii: pd.DataFrame) -> pd.DataFrame:
    """Standard left-join from ITPD-E (trade flows) to CEPII (covariates)."""
    merged = itpde.merge(
        cepii,
        left_on=["exporter_iso3", "importer_iso3", "year"],
        right_on=["iso3_o", "iso3_d", "year"],
        how="left",
    )
    # Drop pairs where CEPII has no covariates (mostly non-sovereign units)
    merged = merged.dropna(subset=["distw_harmonic"])
    return merged


def estimate_ppml_three_way_fe(panel: pd.DataFrame, *,
                                trade_col: str = "trade",
                                dist_col: str = "distw_harmonic",
                                covariates: tuple[str, ...] = ("contig", "comlang_off", "rta")) -> dict:
    """Estimate PPML gravity with exporter-year, importer-year FE.

    Pair FE not yet added (too many groups for statsmodels.GLM at this scale).
    """
    import statsmodels.api as sm

    df = panel.copy()
    df = df.dropna(subset=[trade_col, dist_col] + list(covariates))
    df["log_dist"] = np.log(df[dist_col].clip(lower=1.0))
    df["o_year"] = df["exporter_iso3"] + "_" + df["year"].astype(str)
    df["d_year"] = df["importer_iso3"] + "_" + df["year"].astype(str)
    # Drop singletons
    df = df.groupby("o_year").filter(lambda g: len(g) > 1)
    df = df.groupby("d_year").filter(lambda g: len(g) > 1)

    print(f"  Estimating PPML on {len(df):,} obs ({df['o_year'].nunique()} exp-yr, "
          f"{df['d_year'].nunique()} imp-yr FE)...")

    # Build dummy matrix (memory-intensive but works for ~500 FE per side)
    o_dummies = pd.get_dummies(df["o_year"], drop_first=True, dtype=np.float64)
    d_dummies = pd.get_dummies(df["d_year"], drop_first=True, dtype=np.float64)

    X = pd.concat(
        [df[["log_dist"] + list(covariates)], o_dummies, d_dummies], axis=1
    )
    X = sm.add_constant(X)
    y = df[trade_col].values

    t0 = time.time()
    try:
        model = sm.GLM(y, X, family=sm.families.Poisson()).fit(
            maxiter=100, tol=1e-7
        )
    except Exception as e:
        print(f"  ESTIMATION FAILED: {e}")
        return {"error": str(e)}
    print(f"  Done in {time.time()-t0:.0f}s.")

    out = {
        "n_obs": len(df),
        "coef_log_dist": float(model.params.get("log_dist", np.nan)),
        "se_log_dist": float(model.bse.get("log_dist", np.nan)),
    }
    for c in covariates:
        out[f"coef_{c}"] = float(model.params.get(c, np.nan))
        out[f"se_{c}"] = float(model.bse.get(c, np.nan))
    return out


def main():
    if not ITPDE.exists() or ITPDE.stat().st_size < 1e9:
        print(f"ITPDE_R03.csv not fully extracted yet at {ITPDE}")
        return

    inspect_itpde_schema()

    # Start small: one industry, 10 years
    INDUSTRY_TEST = 1  # Wheat (verify after inspecting schema)
    itpde = load_single_industry(INDUSTRY_TEST, year_min=2010, year_max=2019)
    if len(itpde) == 0:
        print(f"No data for industry {INDUSTRY_TEST}. Adjust id.")
        return

    cepii = load_cepii_panel(2010, 2019)
    print(f"CEPII rows: {len(cepii):,}")

    merged = merge_itpde_with_cepii(itpde, cepii)
    print(f"\nMerged: {len(merged):,} rows ({(merged['trade'] > 0).sum():,} positive)")
    print(f"Industry: {merged['industry_descr'].iloc[0]!r}")

    print("\n=== Baseline PPML (CEPII distw_harmonic) ===")
    res_baseline = estimate_ppml_three_way_fe(merged, dist_col="distw_harmonic")
    print(f"  log-distance coef: {res_baseline.get('coef_log_dist', np.nan):+.4f}  "
          f"(SE {res_baseline.get('se_log_dist', np.nan):.4f})")
    print(f"  contig coef:        {res_baseline.get('coef_contig', np.nan):+.4f}")
    print(f"  rta coef:           {res_baseline.get('coef_rta', np.nan):+.4f}")

    # Save the result
    out = Path(__file__).resolve().parents[1] / "data" / "derived" / "borchert_larch_replication_pilot.csv"
    pd.DataFrame([{"industry_id": INDUSTRY_TEST,
                   "industry_descr": merged['industry_descr'].iloc[0],
                   "n_obs": res_baseline.get("n_obs"),
                   **res_baseline}]).to_csv(out, index=False)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
