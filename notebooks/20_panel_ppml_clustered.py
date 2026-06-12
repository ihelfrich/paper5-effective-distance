"""Panel PPML sweep with proper cluster-robust inference.

Re-runs the 18f specification (pooled 2010/2015/2020, origin-year and
destination-year FE via pyfixest.fepois) with covariance clustered by
country pair (CRV1), which 18f omitted: its saved SEs came from the
default IID covariance and are unusable for inference on trade flows.

Also saves the intra-national panel input rows, so the Wei-proxy values
that enter the regression are an inspectable artifact (the paper quotes
USA 2010 = $13.1T; the only previously saved input predated the proxy
fix).

Outputs:
  replication/data/spec_b_prime_panel_ppml_clustered.csv
  replication/data/panel_intra_inputs.csv
  logs/nb20_<timestamp>.log (via tee at the call site)

Run:
  .venv/bin/python -u notebooks/20_panel_ppml_clustered.py
"""

from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

# Reuse the verified 18f pipeline (panel assembly, disk-interpolated
# d_eff) rather than duplicating it.
_spec = importlib.util.spec_from_file_location(
    "nb18f", _repo / "notebooks" / "18f_spec_b_prime_panel_ppml.py"
)
nb18f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nb18f)

OUT_DIR = _repo / "replication" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def estimate_clustered(df, dist_col, label):
    """fepois with CRV1 clustering on country pair."""
    import pyfixest as pf

    needed = ["X_ij_usd", dist_col, "iso3_o", "iso3_d", "year", "home"]
    df = df.dropna(subset=needed).copy()
    df = df[df[dist_col] > 0].copy()
    df["log_d"] = np.log(df[dist_col])
    df["home"] = pd.to_numeric(df["home"], errors="coerce").fillna(0)
    df["year"] = df["year"].astype(int)
    df["pair_id"] = df["iso3_o"] + "_" + df["iso3_d"]

    base = {"label": label, "n_obs": len(df),
            "n_intra": int(df["home"].sum()),
            "n_zero": int((df["X_ij_usd"] == 0).sum())}
    if len(df) == 0:
        return {**base, "converged": False}

    t0 = time.time()
    try:
        model = pf.fepois(
            "X_ij_usd ~ log_d + home | iso3_o^year + iso3_d^year",
            data=df,
            vcov={"CRV1": "pair_id"},
            iwls_tol=1e-8,
            iwls_maxiter=50,
        )
    except Exception as e:
        print(f"  [{label}]  FAILED: {type(e).__name__}: {e}")
        return {**base, "converged": False}

    coefs = model.coef()
    ses = model.se()
    pvals = model.pvalue()
    ci = model.confint()
    out = {
        **base,
        "coef_log_d": float(coefs["log_d"]),
        "se_log_d": float(ses["log_d"]),
        "p_log_d": float(pvals["log_d"]),
        "coef_home": float(coefs["home"]),
        "se_home": float(ses["home"]),
        "p_home": float(pvals["home"]),
        "ci_home_lo": float(ci.loc["home"].iloc[0]),
        "ci_home_hi": float(ci.loc["home"].iloc[1]),
        "home_multiplier": float(np.exp(coefs["home"])),
        "converged": True,
    }
    print(f"  [{label}]  N={out['n_obs']:,} intra={out['n_intra']} "
          f"log_d={out['coef_log_d']:+.4f}({out['se_log_d']:.4f}) "
          f"home={out['coef_home']:+.4f}({out['se_home']:.4f}) "
          f"ci=[{out['ci_home_lo']:+.3f},{out['ci_home_hi']:+.3f}] "
          f"({time.time()-t0:.0f}s)")
    return out


def main():
    print("=" * 70)
    print("nb20: panel PPML sweep, CRV1 clustered by country pair")
    print("=" * 70)

    area_lookup = nb18f.country_areas()
    print(f"[0] {len(area_lookup)} country areas")

    intra_dfs, inter_dfs = [], []
    for year in (2010, 2015, 2020):
        result = nb18f.build_year_panel(year, area_lookup)
        if result is None:
            continue
        intra_y, inter_y = result
        intra_dfs.append(intra_y)
        inter_dfs.append(inter_y)
        print(f"[1.{year}] intra={len(intra_y):,} inter={len(inter_y):,}")

    intra_panel = pd.concat(intra_dfs, ignore_index=True)
    inter_panel = pd.concat(inter_dfs, ignore_index=True)

    # Save the intra inputs for the Wei-proxy reconciliation.
    intra_panel.to_csv(OUT_DIR / "panel_intra_inputs.csv", index=False)
    usa = intra_panel[(intra_panel.iso3_o == "USA") & (intra_panel.year == 2010)]
    if len(usa):
        print(f"[check] USA 2010 Wei proxy: ${usa['X_ij_usd'].iloc[0]/1e12:.2f}T")

    inter_p = inter_panel.copy()
    inter_p["distance"] = inter_p["distw_harmonic"]

    rows = []
    for spec, col in (("A1_panel_baseline_cepii", "d_HM_cepii"),
                      ("A2_panel_baseline_polygon", "d_HM_polygon")):
        intra_b = intra_panel.copy()
        intra_b["distance"] = intra_b[col]
        panel = pd.concat([inter_p, intra_b], ignore_index=True, sort=False)
        res = estimate_clustered(panel, "distance", spec)
        res.update({"theta": None, "d_min": None, "spec": spec})
        rows.append(res)
        pd.DataFrame(rows).to_csv(
            OUT_DIR / "spec_b_prime_panel_ppml_clustered.csv", index=False)

    for theta in (-3.0, -5.0, -7.0):
        for d_min in (0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0):
            intra_s = intra_panel.copy()
            intra_s["distance"] = intra_s["area_km2"].apply(
                lambda a: nb18f.country_d_eff(a, d_min, theta)
                if np.isfinite(a) else np.nan)
            panel = pd.concat([inter_p, intra_s], ignore_index=True, sort=False)
            res = estimate_clustered(
                panel, "distance", f"Bprime_theta{theta}_dmin{d_min}")
            res.update({"theta": theta, "d_min": d_min, "spec": "B_prime_panel"})
            rows.append(res)
            pd.DataFrame(rows).to_csv(
                OUT_DIR / "spec_b_prime_panel_ppml_clustered.csv", index=False)

    df_out = pd.DataFrame(rows)
    print("\nFinal pivot, home coefficient (SE):")
    b = df_out[df_out["spec"] == "B_prime_panel"]
    print(b.pivot_table(index="d_min", columns="theta",
                        values="coef_home").round(3).to_string())
    print(b.pivot_table(index="d_min", columns="theta",
                        values="se_home").round(3).to_string())
    print("\nDone.")


if __name__ == "__main__":
    main()
