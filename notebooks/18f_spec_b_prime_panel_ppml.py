"""Spec B-prime PANEL PPML — the pivot per Gemini's critique.

The cross-sectional OLS / GLM specifications (nb18d, nb18e) have small
n_intra (6-10 per year) against ~500 FE dummies. A referee can dismiss
the finding as a small-N identification artifact. The fix: stack 2010,
2015, 2020 into a single panel and use origin-year × destination-year
FE (the YPMM 2016 canonical specification), giving n_intra ≈ 24 across
years and ~750 FE per side absorbed via pyfixest's alternating-projection
algorithm.

This is the proper structural-gravity specification:
  X_ij,t = exp(alpha + beta*log_d + gamma*home + eta_{i,t} + mu_{j,t}) * eps

estimated by Poisson PML with origin-year and destination-year FE
absorbed. Uses pyfixest.fepois (high-dim FE Poisson PML, the Python
analog of Stata's ppmlhdfe).

Output: data_derived/spec_b_prime_panel_ppml.csv (one row per
(theta, d_min) cell + baselines).

Run:
  .venv/bin/python -u notebooks/18f_spec_b_prime_panel_ppml.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")
LOG_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/logs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Disk simulation reused from nb18d
DISK_SIM_T5 = {
    0.01: 0.063, 0.1: 0.631, 1: 6.288, 5: 29.449, 10: 51.855,
    50: 149.742, 100: 229.350, 200: 352.875, 500: 640.413,
}
THETA_SCALE_VS_NEG5 = {-3: 3.25, -5: 1.0, -7: 0.75}
DISK_R = 1000.0


def country_d_eff(area_km2, d_min_km, theta):
    if not np.isfinite(area_km2) or area_km2 <= 0:
        return np.nan
    R_country = float(np.sqrt(area_km2 / np.pi))
    d_min_disk_equiv = d_min_km * (DISK_R / R_country)
    keys = sorted(DISK_SIM_T5.keys())
    vals = [DISK_SIM_T5[k] for k in keys]
    d_min_disk_equiv = float(np.clip(d_min_disk_equiv, keys[0], keys[-1]))
    d_eff_disk_t5 = float(np.interp(d_min_disk_equiv, keys, vals))
    d_eff_disk = d_eff_disk_t5 * THETA_SCALE_VS_NEG5.get(theta, 1.0)
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


def build_year_panel(year, area_lookup):
    """Build inter+intra panel for one year."""
    baci = chunked_baci(year)
    if baci is None:
        return None
    cepii = chunked_cepii(year)
    if len(cepii) == 0:
        return None

    map_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates(subset=["iso3num_o"])
    map_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates(subset=["iso3num_d"])
    num_to_iso = {int(n): s for n, s in map_o.values}
    num_to_iso.update({int(n): s for n, s in map_d.values})
    baci["iso3_o"] = baci["i"].map(num_to_iso)
    baci["iso3_d"] = baci["j"].map(num_to_iso)
    baci = baci.dropna(subset=["iso3_o", "iso3_d"])

    cepii_self = cepii[cepii["iso3_o"] == cepii["iso3_d"]][["iso3_o", "distw_harmonic"]].copy()
    cepii_self = cepii_self.rename(columns={"distw_harmonic": "d_HM_cepii"})

    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exp_by_o = baci.groupby("iso3_o", as_index=False)["X_ij_usd"].sum().set_index("iso3_o")["X_ij_usd"].to_dict()
    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        if not np.isfinite(gdp) or gdp <= 0:
            continue
        exp = exp_by_o.get(iso, 0.0)
        area = area_lookup.get(iso, np.nan)
        intra_rows.append({
            "year": year, "iso3_o": iso, "iso3_d": iso,
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
    inter["X_ij_usd"] = inter["X_ij_usd"].fillna(0.0)
    inter["home"] = 0
    inter["year"] = year
    inter["d_HM_cepii"] = inter["distw_harmonic"]
    inter["d_HM_polygon"] = inter["distw_harmonic"]  # placeholder, not used for inter
    inter["area_km2"] = np.nan
    return intra_df, inter


def estimate_panel_ppml(df, dist_col, label, verbose=True):
    """Pooled PPML panel with origin-year and destination-year FE.

    Uses pyfixest.fepois, the Python equivalent of Stata's ppmlhdfe.
    The FE specification iso3_o^year + iso3_d^year creates interacted
    fixed effects (one dummy per (origin, year) pair, one per
    (destination, year) pair). This is the YPMM 2016 canonical
    specification for structural gravity.
    """
    import pyfixest as pf

    needed = ["X_ij_usd", dist_col, "iso3_o", "iso3_d", "year", "home"]
    df = df.dropna(subset=needed).copy()
    df = df[df[dist_col] > 0].copy()
    df["log_d"] = np.log(df[dist_col])
    df["home"] = pd.to_numeric(df["home"], errors="coerce").fillna(0)
    df["year"] = df["year"].astype(int)

    if len(df) == 0:
        return {"label": label, "n_obs": 0, "n_intra": 0,
                "coef_log_d": np.nan, "coef_home": np.nan,
                "home_multiplier": np.nan, "n_zero": 0}

    t0 = time.time()
    # PPML with FE absorbed via alternating projections (Correia-Guimaraes-Zylkin)
    try:
        model = pf.fepois(
            "X_ij_usd ~ log_d + home | iso3_o^year + iso3_d^year",
            data=df,
            iwls_tol=1e-8,
            iwls_maxiter=50,
        )
        home_coef = float(model.coef()["home"])
        log_d = float(model.coef()["log_d"])
        try:
            home_se = float(model.se()["home"])
            log_d_se = float(model.se()["log_d"])
        except Exception:
            home_se = np.nan
            log_d_se = np.nan
        converged = True
    except Exception as e:
        if verbose:
            print(f"  [{label}]  FE-PPML FAILED: {type(e).__name__}: {e}")
        return {"label": label, "n_obs": len(df), "n_intra": int(df["home"].sum()),
                "coef_log_d": np.nan, "coef_home": np.nan,
                "se_log_d": np.nan, "se_home": np.nan,
                "home_multiplier": np.nan,
                "n_zero": int((df["X_ij_usd"] == 0).sum()),
                "converged": False}

    if verbose:
        print(f"  [{label}]  N={len(df):,}  intra={int(df['home'].sum())}  "
              f"zeros={int((df['X_ij_usd']==0).sum())}  "
              f"log_d={log_d:+.4f}({log_d_se:.4f})  "
              f"home={home_coef:+.4f}({home_se:.4f})  "
              f"mult={np.exp(home_coef):.3f}x  ({time.time()-t0:.0f}s)")
    return {
        "label": label, "n_obs": len(df), "n_intra": int(df["home"].sum()),
        "coef_log_d": log_d, "se_log_d": log_d_se,
        "coef_home": home_coef, "se_home": home_se,
        "home_multiplier": float(np.exp(home_coef)),
        "n_zero": int((df["X_ij_usd"] == 0).sum()),
        "converged": True,
    }


def main():
    print("=" * 70)
    print("nb18f: Spec B-prime POOLED PANEL PPML")
    print("Origin-year + Destination-year FE (YPMM 2016 specification)")
    print("Years: 2010, 2015, 2020 stacked")
    print("=" * 70)

    print("\n[0] Loading country areas...")
    area_lookup = country_areas()
    print(f"    {len(area_lookup)} countries")

    YEARS = [2010, 2015, 2020]
    intra_dfs = []
    inter_dfs = []
    for year in YEARS:
        print(f"\n[1.{year}] Loading and building panel for {year}...")
        result = build_year_panel(year, area_lookup)
        if result is None:
            print(f"   skipping {year}")
            continue
        intra_y, inter_y = result
        intra_dfs.append(intra_y)
        inter_dfs.append(inter_y)
        print(f"   intra={len(intra_y):,} inter={len(inter_y):,}")

    intra_panel = pd.concat(intra_dfs, ignore_index=True)
    inter_panel = pd.concat(inter_dfs, ignore_index=True)
    print(f"\nStacked panel:")
    print(f"  inter (positive + zero): {len(inter_panel):,}")
    print(f"  intra rows total: {len(intra_panel):,}")
    print(f"    of which positive Wei proxy: {(intra_panel['X_ij_usd'] > 0).sum():,}")
    print(f"    of which CEPII intra distance available: {intra_panel['d_HM_cepii'].notna().sum():,}")
    valid_intra = (intra_panel['X_ij_usd'] > 0) & intra_panel['d_HM_cepii'].notna()
    print(f"    BOTH: {valid_intra.sum():,}")

    rows = []

    # ── Baseline A1: CEPII intra distance ──
    print(f"\n[2a] Baseline A1 PANEL PPML (CEPII distw_harmonic for intra):")
    intra_a1 = intra_panel.copy()
    intra_a1["distance"] = intra_a1["d_HM_cepii"]
    inter_p = inter_panel.copy()
    inter_p["distance"] = inter_p["distw_harmonic"]
    panel_A1 = pd.concat([inter_p, intra_a1], ignore_index=True, sort=False)
    res = estimate_panel_ppml(panel_A1, dist_col="distance", label="PANEL_A1_baseline_cepii")
    res.update({"theta": None, "d_min": None, "spec": "A1_panel_baseline_cepii"})
    rows.append(res)

    # ── Baseline A2: polygon-area HM ──
    print(f"\n[2b] Baseline A2 PANEL PPML (polygon-area HM for intra):")
    intra_a2 = intra_panel.copy()
    intra_a2["distance"] = intra_a2["d_HM_polygon"]
    panel_A2 = pd.concat([inter_p, intra_a2], ignore_index=True, sort=False)
    res = estimate_panel_ppml(panel_A2, dist_col="distance", label="PANEL_A2_baseline_polygon")
    res.update({"theta": None, "d_min": None, "spec": "A2_panel_baseline_polygon"})
    rows.append(res)

    # ── Spec B-prime sweep ──
    print(f"\n[3] Spec B-prime PANEL PPML sweep (theta × d_min):")
    THETAS = [-3.0, -5.0, -7.0]
    D_MINS = [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0]
    for theta in THETAS:
        for d_min in D_MINS:
            intra_swap = intra_panel.copy()
            intra_swap["distance"] = intra_swap["area_km2"].apply(
                lambda a: country_d_eff(a, d_min, theta) if np.isfinite(a) else np.nan
            )
            panel_b = pd.concat([inter_p, intra_swap], ignore_index=True, sort=False)
            label = f"PANEL_Bprime_theta{theta}_dmin{d_min}km"
            res = estimate_panel_ppml(panel_b, dist_col="distance", label=label)
            res.update({"theta": theta, "d_min": d_min, "spec": "B_prime_panel"})
            rows.append(res)

    df_out = pd.DataFrame(rows)
    out_csv = OUT_DIR / "spec_b_prime_panel_ppml.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"\nSaved to {out_csv}")

    print("\n" + "=" * 70)
    print("PANEL PPML SUMMARY (year-pooled, origin-year + dest-year FE)")
    print("=" * 70)
    print(df_out[["label", "theta", "d_min", "n_obs", "n_intra",
                  "coef_log_d", "coef_home", "home_multiplier"]].to_string(index=False))

    b = df_out[df_out["spec"] == "B_prime_panel"].copy()
    if len(b) > 0:
        print("\nPanel home coefficient pivot:")
        print(b.pivot_table(index="d_min", columns="theta",
                            values="coef_home").round(3).to_string())
        print("\nPanel home multiplier pivot:")
        print(b.pivot_table(index="d_min", columns="theta",
                            values="home_multiplier").round(3).to_string())


if __name__ == "__main__":
    main()
