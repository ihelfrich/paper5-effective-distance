"""Spec B-prime: home-bias multiplier as a function of intra-national d_ii regularization.

The Head-Mayer 2002 closed-form intra-national distance formula has a pole at θ = -2.
For the structurally consistent θ = δ(1-σ) with σ ∈ [4, 8] (i.e., θ ∈ [-7, -3]), the
integral diverges. Whatever value gets assigned to d_ii is a regularization choice.

This notebook computes the gravity-regression home-bias multiplier as a function of
the regularization grain d_min. The d_eff per country is estimated from the
nb19c uniform-disk simulation, scaled to each country's R via the dimensionless
ratio d_min / R:

    R_country = sqrt(area_country / pi)
    ratio    = d_min / R_country
    d_eff    = R_country * disk_lookup(ratio)

where disk_lookup is interpolated from the disk simulation table at R = 1000 km.

For each d_min in a defensible range, we re-estimate the gravity model with the
corrected intra-national distance and report the home-bias multiplier.

Run:
  .venv/bin/python -u notebooks/18b_spec_b_prime.py
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

YEAR = 2010
OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Disk simulation table from nb19c (R = 1000 km, d_HM = 670 km)
# d_min (km) -> d_eff(theta=-5) (km)
DISK_SIM = {
    0.01: 0.063,
    0.1: 0.631,
    1: 6.288,
    5: 29.449,
    10: 51.855,
    50: 149.742,
    100: 229.350,
    200: 352.875,
    500: 640.413,
}
DISK_R = 1000.0  # km, the radius the simulation was run at


def country_d_eff(area_km2: float, d_min_km: float) -> float:
    """Estimate d_ii^eff(theta=-5) for a country of given area at given d_min.

    Uses the dimensionless ratio d_min / R_country to interpolate from the
    disk simulation, then scales back by R_country.
    """
    R_country = float(np.sqrt(area_km2 / np.pi))
    if R_country <= 0 or not np.isfinite(R_country):
        return np.nan
    # dimensionless: what would the same d_min look like in disk-relative units
    # (we scale our d_min, in km, by R_disk/R_country to get the equivalent
    # d_min the disk simulation was tested at)
    d_min_disk_equiv = d_min_km * (DISK_R / R_country)
    keys = sorted(DISK_SIM.keys())
    vals = [DISK_SIM[k] for k in keys]
    # Clip to interpolation range
    d_min_disk_equiv = float(np.clip(d_min_disk_equiv, keys[0], keys[-1]))
    d_eff_disk = float(np.interp(d_min_disk_equiv, keys, vals))
    # Scale back to country
    d_eff_country = d_eff_disk * (R_country / DISK_R)
    # Cap at Head-Mayer closed form (the regularization can't produce a value
    # larger than what the closed form would give at theta = 1)
    d_HM = 0.67 * R_country
    return float(min(d_eff_country, d_HM))


def chunked_baci(year):
    f = Path(f"/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/BACI_HS02_Y{year}_V202401b.csv")
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


def country_areas_from_polygons() -> dict:
    """Compute country area in km^2 from Natural Earth polygons (Mollweide)."""
    from paper5.region_raster import load_country_boundaries
    bnd = load_country_boundaries().to_crs("ESRI:54009")
    return dict(zip(bnd["ISO3"], (bnd.geometry.area / 1e6).astype(float).values))


def estimate_gravity(df, *, dist_col, label):
    df = df.copy()
    df = df[(df["X_ij_usd"] > 0) & (df[dist_col] > 0) &
            (df["gdp_o"] > 0) & (df["gdp_d"] > 0)].copy()
    df["log_X"] = np.log(df["X_ij_usd"])
    df["log_d"] = np.log(df[dist_col])
    df["log_gdp_o"] = np.log(df["gdp_o"])
    df["log_gdp_d"] = np.log(df["gdp_d"])

    o_dum = pd.get_dummies(df["iso3_o"], prefix="o", drop_first=True, dtype=np.float64)
    d_dum = pd.get_dummies(df["iso3_d"], prefix="d", drop_first=True, dtype=np.float64)

    X = pd.concat([df[["log_d", "home", "log_gdp_o", "log_gdp_d", "contig", "comlang_off"]],
                   o_dum, d_dum], axis=1)
    X = sm.add_constant(X)
    y = df["log_X"].values

    t0 = time.time()
    model = sm.OLS(y, X).fit()
    home_coef = float(model.params.get("home", np.nan))
    home_se   = float(model.bse.get("home", np.nan))
    log_d     = float(model.params.get("log_d", np.nan))
    log_d_se  = float(model.bse.get("log_d", np.nan))
    print(f"  [{label}]  N={len(df):,}  log_d={log_d:+.4f}({log_d_se:.4f})  "
          f"home={home_coef:+.4f}({home_se:.4f})  mult={np.exp(home_coef):.2f}x  "
          f"({time.time()-t0:.1f}s)")
    return {
        "label": label, "n_obs": len(df),
        "coef_log_d": log_d, "se_log_d": log_d_se,
        "coef_home": home_coef, "se_home": home_se,
        "home_multiplier": float(np.exp(home_coef)),
        "r2": float(model.rsquared),
    }


def main():
    print(f"=== Spec B-prime: home-bias multiplier vs intra-national regularization (year {YEAR}) ===\n")

    print(f"[1] Loading BACI {YEAR}...")
    baci = chunked_baci(YEAR)
    print(f"    {len(baci):,} pairs")

    print(f"[2] Loading CEPII Gravity {YEAR}...")
    cepii = chunked_cepii(YEAR)
    print(f"    {len(cepii):,} CEPII rows")

    # ISO num → ISO3 mapping
    map_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates(subset=["iso3num_o"])
    map_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates(subset=["iso3num_d"])
    num_to_iso = {int(n): s for n, s in map_o.values}
    num_to_iso.update({int(n): s for n, s in map_d.values})
    baci["iso3_o"] = baci["i"].map(num_to_iso)
    baci["iso3_d"] = baci["j"].map(num_to_iso)
    baci = baci.dropna(subset=["iso3_o", "iso3_d"])

    print(f"[3] Country areas from Natural Earth (Mollweide)...")
    area_lookup = country_areas_from_polygons()
    print(f"    {len(area_lookup)} countries")

    # Intra-national rows
    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exports_by_o = baci.groupby("iso3_o", as_index=False)["X_ij_usd"].sum().set_index("iso3_o")["X_ij_usd"].to_dict()

    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        if not np.isfinite(gdp) or gdp <= 0: continue
        exp = exports_by_o.get(iso, 0.0)
        x_ii = max(0.0, gdp - exp)
        area = area_lookup.get(iso, np.nan)
        d_HM_iso = 0.67 * np.sqrt(area / np.pi) if np.isfinite(area) else np.nan
        intra_rows.append({
            "year": YEAR, "iso3_o": iso, "iso3_d": iso,
            "X_ij_usd": x_ii, "home": 1, "gdp_o": gdp, "gdp_d": gdp,
            "contig": 0, "comlang_off": 1,
            "area_km2": area,
            "d_HM": d_HM_iso,
        })
    intra_df = pd.DataFrame(intra_rows)
    print(f"    {len(intra_df):,} intra-national rows")

    # Inter-national rows
    inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]][
        ["year", "iso3_o", "iso3_d", "gdp_o", "gdp_d",
         "distw_harmonic", "contig", "comlang_off"]
    ].copy()
    inter = inter.merge(baci[["iso3_o", "iso3_d", "X_ij_usd"]],
                        on=["iso3_o", "iso3_d"], how="left")
    inter["X_ij_usd"] = inter["X_ij_usd"].fillna(0.0)
    inter["home"] = 0
    inter["d_HM"] = inter["distw_harmonic"]
    print(f"    {len(inter):,} inter-national rows")

    # ─ Baseline (Spec A): CEPII Head-Mayer intra-national distance ─
    print("\n[4] Spec A baseline (CEPII Head-Mayer intra-national distance):")
    intra_df["distance"] = intra_df["d_HM"]
    inter["distance"] = inter["distw_harmonic"]
    panel = pd.concat([inter, intra_df], ignore_index=True, sort=False)
    res_A = estimate_gravity(panel, dist_col="distance", label="A_baseline")

    # ─ Spec B-prime sweeps: substitute d_eff(d_min) for intra-national ─
    print(f"\n[5] Spec B-prime: home multiplier vs d_min ∈ [1, 5, 10, 50, 100, 200] km:")
    results = [res_A]
    for d_min in [1, 5, 10, 50, 100, 200]:
        intra_swap = intra_df.copy()
        intra_swap["distance"] = intra_swap.apply(
            lambda r: country_d_eff(r["area_km2"], d_min) if np.isfinite(r["area_km2"]) else r["d_HM"],
            axis=1
        )
        panel_b = pd.concat([inter, intra_swap], ignore_index=True, sort=False)
        res = estimate_gravity(panel_b, dist_col="distance", label=f"B_prime_dmin_{d_min}km")
        results.append(res)

    df_results = pd.DataFrame(results)
    out = OUT_DIR / "spec_b_prime_regularization_sweep.csv"
    df_results.to_csv(out, index=False)
    print(f"\n[6] Saved to {out}")

    # ─ Print the table ─
    print(f"\n=== Home-bias multiplier as a function of intra-national regularization ===")
    print(f"  {'spec':<22}  {'home coef':>10}  {'home SE':>10}  {'multiplier':>12}  {'Δ vs A':>10}")
    base = res_A["coef_home"]
    for r in results:
        delta = r["coef_home"] - base
        print(f"  {r['label']:<22}  {r['coef_home']:>+10.4f}  {r['se_home']:>10.4f}  "
              f"{r['home_multiplier']:>9.2f}×  {delta:>+10.4f}")

    print(f"\n=== Interpretation ===")
    largest_drop = min(r["coef_home"] for r in results) - base
    if largest_drop < -0.5:
        print(f"  Home-bias multiplier shifts meaningfully under structural regularization.")
        print(f"  Largest log-coef change: {largest_drop:+.3f}.")
        print(f"  This supports the 'wrong θ' framing as quantitatively consequential.")
    elif largest_drop < -0.1:
        print(f"  Home-bias multiplier shifts modestly under regularization.")
        print(f"  Largest log-coef change: {largest_drop:+.3f}.")
        print(f"  The 'wrong θ' framing has empirical bite but smaller than a JIE editor would want.")
    else:
        print(f"  Home-bias multiplier is robust to regularization choice.")
        print(f"  Largest log-coef change: {largest_drop:+.3f}.")
        print(f"  The 'wrong θ' framing is theoretically real but empirically narrow.")


if __name__ == "__main__":
    main()
