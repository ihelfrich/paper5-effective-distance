"""Spec B-prime, FIXED regression specification.

Diagnosis of nb18c bug (Task #35):
  nb18c used sm.OLS with manually-constructed design matrix that included
  log_gdp_o, log_gdp_d, contig, comlang_off ALONGSIDE origin and destination
  dummies. In a single-year cross-section, log_gdp_o is exactly collinear with
  origin FE (one value of log_gdp_o per origin), and log_gdp_d with destination
  FE. statsmodels' sm.OLS does NOT auto-drop colinear columns; it falls back to
  the pseudo-inverse, which produces unstable coefficients that get redistributed
  across the colinear set. This corrupted the home coefficient.

  nb18 v4 (the known-working version) used the formula API:
    smf.ols('log_X ~ log_d + home + C(iso_o) + C(iso_d)', data=df)
  which handles dummy creation and collinearity correctly, and crucially omits
  log_gdp_o and log_gdp_d (since the FE absorb them).

Fix: this notebook (nb18d) uses the nb18 v4 regression formula exactly, keeps the
(theta, d_min) sweep machinery from nb18c, and runs three years (2010, 2015, 2020).

The expected baseline output for year 2010:
  log_d  = -1.90 (SE ~0.02)
  home   = +2.32 (SE ~0.67)  → multiplier ~10x  (consistent with McCallum-AvW)

If we reproduce this, the bug is confirmed isolated to the regression spec.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")
FIG_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Disk simulation tables (reused from nb18c) ──────────────────────────────
DISK_SIM_T5 = {
    0.01: 0.063, 0.1: 0.631, 1: 6.288, 5: 29.449, 10: 51.855,
    50: 149.742, 100: 229.350, 200: 352.875, 500: 640.413,
}
THETA_SCALE_VS_NEG5 = {-3: 3.25, -5: 1.0, -7: 0.75}
DISK_R = 1000.0


def country_d_eff(area_km2: float, d_min_km: float, theta: float) -> float:
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


# ── Data loaders ────────────────────────────────────────────────────────────

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
    agg["year"] = year
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


def country_areas() -> dict:
    from paper5.region_raster import load_country_boundaries
    bnd = load_country_boundaries().to_crs("ESRI:54009")
    return dict(zip(bnd["ISO3"], (bnd.geometry.area / 1e6).astype(float).values))


# ── FIXED gravity regression: matches nb18 v4 exactly ───────────────────────

def estimate_gravity_fixed(df, *, dist_col, label, verbose=True):
    """Estimate log X_ij = a + b * log_d + c * home + origin_FE + dest_FE.

    No log_gdp_o or log_gdp_d (collinear with FE in single-year cross-section).
    Uses smf.ols formula API, which auto-drops collinear columns and handles
    dummy creation correctly.
    """
    needed = ["X_ij_usd", dist_col, "iso3_o", "iso3_d", "home"]
    df = df.dropna(subset=needed).copy()
    df = df[(df["X_ij_usd"] > 0) & (df[dist_col] > 0)].copy()
    df["log_X"] = np.log(df["X_ij_usd"])
    df["log_d"] = np.log(df[dist_col])
    df["home"] = pd.to_numeric(df["home"], errors="coerce").fillna(0)

    if len(df) == 0:
        return {"label": label, "n_obs": 0, "n_intra": 0,
                "coef_log_d": np.nan, "se_log_d": np.nan,
                "coef_home": np.nan, "se_home": np.nan,
                "home_multiplier": np.nan, "r2": np.nan}

    t0 = time.time()
    formula = "log_X ~ log_d + home + C(iso3_o) + C(iso3_d)"
    model = smf.ols(formula, data=df).fit()
    home_coef = float(model.params.get("home", np.nan))
    home_se   = float(model.bse.get("home", np.nan))
    log_d     = float(model.params.get("log_d", np.nan))
    log_d_se  = float(model.bse.get("log_d", np.nan))
    if verbose:
        print(f"  [{label}]  N={len(df):,}  intra={int(df['home'].sum())}  "
              f"log_d={log_d:+.4f}({log_d_se:.4f})  "
              f"home={home_coef:+.4f}({home_se:.4f})  mult={np.exp(home_coef):.2f}x  "
              f"R2={model.rsquared:.3f}  ({time.time()-t0:.0f}s)")
    return {
        "label": label, "n_obs": len(df), "n_intra": int(df["home"].sum()),
        "coef_log_d": log_d, "se_log_d": log_d_se,
        "coef_home": home_coef, "se_home": home_se,
        "home_multiplier": float(np.exp(home_coef)),
        "r2": float(model.rsquared),
    }


# ── Build panel ─────────────────────────────────────────────────────────────

def build_panel_for_year(year, baci, cepii, area_lookup):
    map_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates(subset=["iso3num_o"])
    map_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates(subset=["iso3num_d"])
    num_to_iso = {int(n): s for n, s in map_o.values}
    num_to_iso.update({int(n): s for n, s in map_d.values})
    baci_m = baci.copy()
    baci_m["iso3_o"] = baci_m["i"].map(num_to_iso)
    baci_m["iso3_d"] = baci_m["j"].map(num_to_iso)
    baci_m = baci_m.dropna(subset=["iso3_o", "iso3_d"])

    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exp_by_o = baci_m.groupby("iso3_o", as_index=False)["X_ij_usd"].sum().set_index("iso3_o")["X_ij_usd"].to_dict()

    # Pull CEPII's stored intra distance (its distw_harmonic for i==j; this is
    # the Head-Mayer closed form that nb18 v4 used and that the baseline should match).
    cepii_self = cepii[cepii["iso3_o"] == cepii["iso3_d"]][["iso3_o", "distw_harmonic"]].copy()
    cepii_self = cepii_self.rename(columns={"distw_harmonic": "d_HM_cepii"})

    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        if not np.isfinite(gdp) or gdp <= 0:
            continue
        exp = exp_by_o.get(iso, 0.0)
        area = area_lookup.get(iso, np.nan)
        intra_rows.append({
            "year": year, "iso3_o": iso, "iso3_d": iso,
            "X_ij_usd": max(0.0, gdp - exp), "home": 1,
            "gdp_o": gdp, "gdp_d": gdp,
            "area_km2": area,
            "d_HM_polygon": 0.67 * np.sqrt(area / np.pi) if np.isfinite(area) else np.nan,
        })
    intra_df = pd.DataFrame(intra_rows)
    intra_df = intra_df.merge(cepii_self, on="iso3_o", how="left")

    inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]][
        ["year", "iso3_o", "iso3_d", "gdp_o", "gdp_d", "distw_harmonic"]
    ].copy()
    inter = inter.merge(baci_m[["iso3_o", "iso3_d", "X_ij_usd"]],
                        on=["iso3_o", "iso3_d"], how="left")
    inter["X_ij_usd"] = inter["X_ij_usd"].fillna(0.0)
    inter["home"] = 0
    return intra_df, inter


def sweep_year(year, area_lookup, theta_values, d_min_values, use_cepii_intra_for_baseline=True):
    print(f"\n========== YEAR {year} ==========")
    print(f"[1] Loading BACI {year}...")
    baci = chunked_baci(year)
    if baci is None:
        return []
    print(f"    {len(baci):,} pairs")

    print(f"[2] Loading CEPII Gravity {year}...")
    cepii = chunked_cepii(year)
    if len(cepii) == 0:
        return []
    print(f"    {len(cepii):,} CEPII rows")

    intra_df, inter = build_panel_for_year(year, baci, cepii, area_lookup)
    print(f"    intra={len(intra_df):,} inter={len(inter):,}")

    rows = []

    # Baseline A1: CEPII's distw_harmonic for intra (matches nb18 v4)
    print(f"\n[3a] Baseline A1 (CEPII distw_harmonic for intra — matches nb18 v4):")
    intra_a1 = intra_df.copy()
    intra_a1["distance"] = intra_a1["d_HM_cepii"]
    inter_w = inter.copy()
    inter_w["distance"] = inter_w["distw_harmonic"]
    panel_A1 = pd.concat([inter_w, intra_a1], ignore_index=True, sort=False)
    res = estimate_gravity_fixed(panel_A1, dist_col="distance",
                                  label=f"A1_baseline_cepii_y{year}")
    res.update({"year": year, "theta": None, "d_min": None,
                "spec": "A1_baseline_cepii_intra"})
    rows.append(res)

    # Baseline A2: polygon-area HM for intra (what nb18c used)
    print(f"\n[3b] Baseline A2 (polygon-area HM for intra — what nb18c used):")
    intra_a2 = intra_df.copy()
    intra_a2["distance"] = intra_a2["d_HM_polygon"]
    panel_A2 = pd.concat([inter_w, intra_a2], ignore_index=True, sort=False)
    res = estimate_gravity_fixed(panel_A2, dist_col="distance",
                                  label=f"A2_baseline_polygon_y{year}")
    res.update({"year": year, "theta": None, "d_min": None,
                "spec": "A2_baseline_polygon_intra"})
    rows.append(res)

    # Spec B-prime sweep
    print(f"\n[4] Spec B-prime sweep (theta × d_min):")
    for theta in theta_values:
        for d_min in d_min_values:
            intra_swap = intra_df.copy()
            intra_swap["distance"] = intra_swap.apply(
                lambda r: country_d_eff(r["area_km2"], d_min, theta)
                          if np.isfinite(r["area_km2"]) else r["d_HM_polygon"],
                axis=1
            )
            panel_b = pd.concat([inter_w, intra_swap], ignore_index=True, sort=False)
            label = f"Bprime_y{year}_theta{theta}_dmin{d_min}km"
            res = estimate_gravity_fixed(panel_b, dist_col="distance", label=label)
            res.update({"year": year, "theta": theta, "d_min": d_min,
                        "spec": "B_prime"})
            rows.append(res)

    return rows


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("nb18d: Spec B-prime with FIXED regression (matches nb18 v4)")
    print("=" * 70)

    print("\n[0] Loading country boundaries (for polygon-area d_HM fallback)...")
    area_lookup = country_areas()
    print(f"    {len(area_lookup)} countries")

    YEARS = [2010, 2015, 2020]
    THETAS = [-3.0, -5.0, -7.0]
    D_MINS = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]

    all_rows = []
    for year in YEARS:
        rows = sweep_year(year, area_lookup, THETAS, D_MINS)
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    out_csv = OUT_DIR / "spec_b_prime_fixed.csv"
    df.to_csv(out_csv, index=False)
    print(f"\n\nSaved to {out_csv}")

    # Comparison summary
    print("\n" + "=" * 70)
    print("BASELINE COMPARISON: A1 (CEPII intra) vs A2 (polygon intra)")
    print("Expected for 2010 (from nb18 v4): log_d=-1.90, home=+2.32, mult≈10x")
    print("=" * 70)
    base = df[df["spec"].isin(["A1_baseline_cepii_intra", "A2_baseline_polygon_intra"])]
    print(base[["year", "spec", "n_obs", "n_intra",
                "coef_log_d", "coef_home", "home_multiplier", "r2"]].to_string())

    # Summary of Spec B-prime
    print("\n" + "=" * 70)
    print("SPEC B-PRIME: home coefficient across (theta, d_min)")
    print("=" * 70)
    b = df[df["spec"] == "B_prime"].copy()
    if len(b) > 0:
        print("Mean home coefficient across years:")
        print(b.pivot_table(index="d_min", columns="theta",
                            values="coef_home", aggfunc="mean").round(3).to_string())
        print("\nMean home multiplier (exp(coef)) across years:")
        print(b.pivot_table(index="d_min", columns="theta",
                            values="home_multiplier", aggfunc="mean").round(3).to_string())


if __name__ == "__main__":
    main()
