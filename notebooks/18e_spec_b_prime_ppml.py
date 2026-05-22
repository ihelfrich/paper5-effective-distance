"""Spec B-prime, PPML robustness check.

The biggest expected referee critique of Paper 5a (referee_responses.md
Critique #4) is that the OLS specification will not survive PPML. This
notebook tests that.

Specification:
  Y_ij = exp(alpha + beta * log_d + gamma * home + eta_i + mu_j) * eps_ij
  estimated by Poisson PML with country origin and destination FE.
  Uses statsmodels GLM(family=Poisson, link='log').

Sample: 2010 only (the cleanest panel with n_intra=10). One sweep over
theta=-5 and d_min in {0.5, 1, 5, 10, 50, 100, 500} km.

Expected output (if Paper 5a's central claim survives):
  - Within-year, home coefficient gamma_hat depends monotonically on d_min
  - Sign-flip threshold may shift but should still exist
  - Magnitude of gamma_hat may shrink (PPML is more efficient for zeros)
  but the QUALITATIVE pattern should hold

If the pattern fails to survive, that's important honest evidence and
Section 5.3 needs to acknowledge it.

Run: .venv/bin/python -u notebooks/18e_spec_b_prime_ppml.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Disk simulation reused from nb18d
DISK_SIM_T5 = {
    0.01: 0.063, 0.1: 0.631, 1: 6.288, 5: 29.449, 10: 51.855,
    50: 149.742, 100: 229.350, 200: 352.875, 500: 640.413,
}
DISK_R = 1000.0


def country_d_eff(area_km2, d_min_km, theta=-5.0):
    if not np.isfinite(area_km2) or area_km2 <= 0:
        return np.nan
    R_country = float(np.sqrt(area_km2 / np.pi))
    d_min_disk_equiv = d_min_km * (DISK_R / R_country)
    keys = sorted(DISK_SIM_T5.keys())
    vals = [DISK_SIM_T5[k] for k in keys]
    d_min_disk_equiv = float(np.clip(d_min_disk_equiv, keys[0], keys[-1]))
    d_eff_disk = float(np.interp(d_min_disk_equiv, keys, vals))
    d_eff_country = d_eff_disk * (R_country / DISK_R)
    cap = 0.9 * R_country
    return float(min(d_eff_country, cap))


def chunked_baci(year):
    f = Path(f"/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/BACI_HS02_Y{year}_V202401b.csv")
    if not f.exists():
        return None
    pieces = []
    for chunk in pd.read_csv(f, usecols=["i", "j", "v"], chunksize=1_000_000):
        chunk["v"] = pd.to_numeric(chunk["v"], errors="coerce")
        chunk = chunk.dropna(subset=["v"])
        pieces.append(chunk.groupby(["i", "j"], as_index=False)["v"].sum())
    agg = pd.concat(pieces).groupby(["i", "j"], as_index=False)["v"].sum()
    agg["X_ij_usd"] = agg["v"] * 1000
    return agg


def chunked_cepii(year):
    f = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")
    pieces = []
    for chunk in pd.read_csv(f, dtype={"iso3_o": str, "iso3_d": str, "year": "Int64"},
                             chunksize=500_000, low_memory=False):
        sub = chunk[chunk["year"] == year]
        if len(sub) > 0:
            pieces.append(sub)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def country_areas():
    from paper5.region_raster import load_country_boundaries
    bnd = load_country_boundaries().to_crs("ESRI:54009")
    return dict(zip(bnd["ISO3"], (bnd.geometry.area / 1e6).astype(float).values))


def estimate_ppml(df, dist_col, label, verbose=True):
    """PPML gravity regression with C(iso_o) + C(iso_d) FE."""
    import statsmodels.formula.api as smf

    needed = ["X_ij_usd", dist_col, "iso3_o", "iso3_d", "home"]
    df = df.dropna(subset=needed).copy()
    df = df[df[dist_col] > 0].copy()
    df["log_d"] = np.log(df[dist_col])
    df["home"] = pd.to_numeric(df["home"], errors="coerce").fillna(0)

    if len(df) == 0:
        return {"label": label, "n_obs": 0, "n_intra": 0,
                "coef_log_d": np.nan, "coef_home": np.nan,
                "home_multiplier": np.nan, "n_zero": 0}

    t0 = time.time()
    formula = "X_ij_usd ~ log_d + home + C(iso3_o) + C(iso3_d)"
    try:
        model = smf.glm(formula, data=df,
                        family=sm.families.Poisson(link=sm.families.links.Log())).fit(
                            method="bfgs", maxiter=200, tol=1e-6)
    except Exception as e:
        if verbose:
            print(f"  [{label}]  GLM FAILED: {e}")
        return {"label": label, "n_obs": len(df), "n_intra": int(df["home"].sum()),
                "coef_log_d": np.nan, "coef_home": np.nan,
                "home_multiplier": np.nan, "n_zero": int((df["X_ij_usd"] == 0).sum())}

    home_coef = float(model.params.get("home", np.nan))
    log_d = float(model.params.get("log_d", np.nan))
    if verbose:
        print(f"  [{label}]  N={len(df):,}  intra={int(df['home'].sum())}  "
              f"zeros={int((df['X_ij_usd']==0).sum())}  "
              f"log_d={log_d:+.4f}  home={home_coef:+.4f}  "
              f"mult={np.exp(home_coef):.2f}x  ({time.time()-t0:.0f}s)")
    return {
        "label": label, "n_obs": len(df), "n_intra": int(df["home"].sum()),
        "coef_log_d": log_d,
        "coef_home": home_coef,
        "home_multiplier": float(np.exp(home_coef)),
        "n_zero": int((df["X_ij_usd"] == 0).sum()),
    }


def main():
    print("=" * 70)
    print("nb18e: Spec B-prime PPML robustness check")
    print("Year 2010 only, theta = -5, d_min sweep")
    print("=" * 70)

    YEAR = 2010

    print(f"\n[0] Loading country areas...")
    area_lookup = country_areas()
    print(f"    {len(area_lookup)} countries")

    print(f"\n[1] Loading BACI {YEAR}...")
    baci = chunked_baci(YEAR)
    print(f"    {len(baci):,} pairs")

    print(f"\n[2] Loading CEPII Gravity {YEAR}...")
    cepii = chunked_cepii(YEAR)
    print(f"    {len(cepii):,} CEPII rows")

    # Build panel
    map_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates(subset=["iso3num_o"])
    map_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates(subset=["iso3num_d"])
    num_to_iso = {int(n): s for n, s in map_o.values}
    num_to_iso.update({int(n): s for n, s in map_d.values})
    baci["iso3_o"] = baci["i"].map(num_to_iso)
    baci["iso3_d"] = baci["j"].map(num_to_iso)
    baci = baci.dropna(subset=["iso3_o", "iso3_d"])

    cepii_self = cepii[cepii["iso3_o"] == cepii["iso3_d"]][["iso3_o", "distw_harmonic"]].copy()
    cepii_self = cepii_self.rename(columns={"distw_harmonic": "d_HM_cepii"})

    # Intra rows (Wei proxy)
    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exp_by_o = baci.groupby("iso3_o", as_index=False)["X_ij_usd"].sum().set_index("iso3_o")["X_ij_usd"].to_dict()
    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        if not np.isfinite(gdp) or gdp <= 0:
            continue
        exp = exp_by_o.get(iso, 0.0)
        area = area_lookup.get(iso, np.nan)
        intra_rows.append({
            "iso3_o": iso, "iso3_d": iso,
            "X_ij_usd": max(0.0, gdp - exp), "home": 1,
            "area_km2": area,
            "d_HM_polygon": 0.67 * np.sqrt(area / np.pi) if np.isfinite(area) else np.nan,
        })
    intra_df = pd.DataFrame(intra_rows)
    intra_df = intra_df.merge(cepii_self, on="iso3_o", how="left")

    inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]][
        ["iso3_o", "iso3_d", "distw_harmonic"]
    ].copy()
    inter = inter.merge(baci[["iso3_o", "iso3_d", "X_ij_usd"]],
                        on=["iso3_o", "iso3_d"], how="left")
    inter["X_ij_usd"] = inter["X_ij_usd"].fillna(0.0)  # PPML keeps zeros
    inter["home"] = 0

    print(f"    intra={len(intra_df):,} inter={len(inter):,}")
    print(f"    inter rows with X=0: {int((inter['X_ij_usd']==0).sum())}")

    rows = []

    # Baseline (CEPII intra)
    print(f"\n[3] Baseline (CEPII intra, theta=-1 equivalent):")
    intra_a = intra_df.copy()
    intra_a["distance"] = intra_a["d_HM_cepii"]
    inter_w = inter.copy()
    inter_w["distance"] = inter_w["distw_harmonic"]
    panel = pd.concat([inter_w, intra_a], ignore_index=True, sort=False)
    res = estimate_ppml(panel, dist_col="distance", label="PPML_baseline_cepii")
    res["d_min"] = None
    rows.append(res)

    # Sweep
    print(f"\n[4] Sweep at theta=-5:")
    for d_min in [0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0]:
        intra_swap = intra_df.copy()
        intra_swap["distance"] = intra_swap["area_km2"].apply(
            lambda a: country_d_eff(a, d_min, -5.0) if np.isfinite(a) else np.nan
        )
        panel = pd.concat([inter_w, intra_swap], ignore_index=True, sort=False)
        res = estimate_ppml(panel, dist_col="distance", label=f"PPML_theta-5_dmin{d_min}km")
        res["d_min"] = d_min
        rows.append(res)

    df_out = pd.DataFrame(rows)
    out_csv = OUT_DIR / "spec_b_prime_ppml.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"\nSaved to {out_csv}")
    print()
    print("=== PPML summary ===")
    print(df_out[["label", "d_min", "n_obs", "n_intra", "n_zero",
                  "coef_log_d", "coef_home", "home_multiplier"]].to_string(index=False))


if __name__ == "__main__":
    main()
