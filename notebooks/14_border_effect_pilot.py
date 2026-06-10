"""Border-effect pilot: does raster-resolution d_ii^eff shrink the home-bias coefficient?

The cleanest possible test of the hypothesis. Keeping the inter-national
distance fixed at CEPII's distw_harmonic, we swap ONLY the intra-national
distance between two specifications:

  Spec A (baseline):   d_ii = 0.67 * sqrt(area_i/π)  [Head-Mayer closed form, as used by CEPII]
  Spec B (treatment):  d_ii = d_ii^eff(θ = -5)        [our raster-resolution CES measure]

Then estimate the gravity equation with home-bias indicator and compare δ
(the home-bias coefficient) across specifications.

Required inputs:
  1. BACI HS02 bilateral trade flows (2010 cross-section)
  2. CEPII Gravity covariates (GDP, distw, etc.)
  3. Intra-national flow proxy: Y_i - exports_i
  4. Per-country d_ii^eff from notebook 13

Run:
  .venv/bin/python notebooks/14_border_effect_pilot.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

# Paths
BACI_DIR = Path("/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b")
CEPII_GRAVITY = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")
DII_EFF_CSV = _repo / "data" / "derived" / "real_country_internal_distance.csv"

PILOT_YEAR = 2010  # GHS-POP epoch + BACI coverage + nightlights cached


def load_baci_year(year: int) -> pd.DataFrame:
    """Load BACI bilateral aggregated by country-pair-year (USD)."""
    f = BACI_DIR / f"BACI_HS02_Y{year}_V202401b.csv"
    if not f.exists():
        raise FileNotFoundError(f"BACI file not found: {f}")
    print(f"  Loading {f.name} (this may take ~30 seconds)...")
    # BACI columns: t (year), i (exporter ISO numeric), j (importer ISO numeric),
    #              k (HS6 product), v (value $1k), q (quantity tonnes)
    df = pd.read_csv(f, dtype={"t": int, "i": int, "j": int, "k": str},
                     usecols=["t", "i", "j", "v"])
    df = df.dropna(subset=["v"])
    df["v"] = pd.to_numeric(df["v"], errors="coerce")
    df = df.dropna(subset=["v"])
    # Aggregate to country-pair total
    agg = df.groupby(["i", "j"], as_index=False).agg(X_ij_usd=("v", "sum"))
    # BACI v is in thousands of USD
    agg["X_ij_usd"] *= 1000
    return agg


def load_cepii_gravity(year: int) -> pd.DataFrame:
    """Load CEPII Gravity covariates for a single year."""
    print(f"  Loading CEPII Gravity for {year}...")
    df = pd.read_csv(CEPII_GRAVITY, dtype={"iso3_o": str, "iso3_d": str, "year": int})
    df = df[df["year"] == year].copy()
    keep_cols = ["year", "iso3_o", "iso3_d", "iso3num_o", "iso3num_d",
                 "gdp_o", "gdp_d", "pop_o", "pop_d",
                 "dist", "distcap", "distw_arithmetic", "distw_harmonic",
                 "contig", "comlang_off", "col_dep_ever", "gatt_o", "gatt_d"]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df[keep_cols].copy()
    # CEPII Gravity reports GDP in thousands of current USD (USA 2010 is
    # 1.4964e10 there vs ~1.4964e13 actual dollars). BACI flows get converted
    # to USD in load_baci_year, so convert GDP to USD too; otherwise the Wei
    # proxy GDP - exports goes negative for every country with exports above
    # 0.1% of GDP and the intra-national sample collapses.
    for c in ("gdp_o", "gdp_d"):
        if c in df.columns:
            df[c] = df[c] * 1000.0
    return df


def country_area_km2_from_natural_earth(iso3_codes: list[str]) -> dict[str, float]:
    """Compute country area in km² from Natural Earth polygons (Mollweide)."""
    from paper5.region_raster import load_country_boundaries
    bnd = load_country_boundaries()
    # Reproject to Mollweide (equal-area)
    bnd_eq = bnd.to_crs("ESRI:54009")
    bnd_eq["area_km2"] = bnd_eq.geometry.area / 1e6
    return dict(zip(bnd_eq["ISO3"], bnd_eq["area_km2"]))


def assemble_panel(year: int) -> pd.DataFrame:
    """Build the gravity panel:
      - All bilateral flows (i != j) from BACI for the year
      - Intra-national flows (i == j) = max(0, GDP_i - sum_j BACI[i,j])
      - Bilateral distance = CEPII distw_harmonic for i != j, Head-Mayer for i = j (baseline)
      - GDP, contig, comlang, colony, gatt covariates from CEPII
    """
    print(f"\n[1/4] Loading BACI {year}...")
    baci = load_baci_year(year)
    print(f"      {len(baci):,} bilateral observations")

    print(f"\n[2/4] Loading CEPII Gravity {year}...")
    cepii = load_cepii_gravity(year)
    print(f"      {len(cepii):,} country-pair rows")

    # Map BACI's iso3num codes to ISO3 letters using CEPII's mapping.
    # CEPII has NaN iso3num for some territories; NaN keys make the mapper's
    # index non-unique, so drop them and cast to int to match BACI's codes.
    _m_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates()
    _m_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates()
    iso3num_to_iso3 = {int(k): v for k, v in _m_o.values}
    iso3num_to_iso3.update({int(k): v for k, v in _m_d.values})
    baci["iso3_o"] = baci["i"].map(iso3num_to_iso3)
    baci["iso3_d"] = baci["j"].map(iso3num_to_iso3)
    baci = baci.dropna(subset=["iso3_o", "iso3_d"])
    print(f"      {len(baci):,} BACI rows after ISO3 mapping")

    # Merge BACI flows onto CEPII covariates (inner join — both must have the pair)
    cepii_inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]].copy()
    panel_inter = cepii_inter.merge(
        baci[["iso3_o", "iso3_d", "X_ij_usd"]],
        on=["iso3_o", "iso3_d"], how="left"
    )
    # Many country-pair entries in CEPII won't have BACI flows (e.g., zeros).
    # Fill missing with 0.
    panel_inter["X_ij_usd"] = panel_inter["X_ij_usd"].fillna(0.0)
    panel_inter["home"] = 0
    print(f"\n[3/4] Inter-national panel: {len(panel_inter):,} obs, "
          f"{(panel_inter['X_ij_usd'] > 0).sum():,} with positive flows")

    # Intra-national: one row per country, X_ii = max(0, GDP_i - sum_j BACI[i,j])
    print(f"\n[4/4] Building intra-national flows (Y_i - exports)...")
    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exports_by_origin = baci.groupby("iso3_o", as_index=False)["X_ij_usd"].sum()
    exports_by_origin = exports_by_origin.rename(columns={"X_ij_usd": "exports_usd"})
    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        exp = exports_by_origin[exports_by_origin["iso3_o"] == iso]["exports_usd"].sum()
        x_ii = max(0.0, gdp - exp)
        intra_rows.append({"iso3_o": iso, "iso3_d": iso, "X_ij_usd": x_ii,
                           "gdp_o": gdp, "gdp_d": gdp, "home": 1})
    panel_intra = pd.DataFrame(intra_rows)
    print(f"      {len(panel_intra):,} intra-national obs")

    # Compute country area for Head-Mayer baseline d_ii
    print(f"      Computing country areas...")
    iso_set = set(panel_intra["iso3_o"].tolist())
    areas = country_area_km2_from_natural_earth(list(iso_set))
    panel_intra["area_km2"] = panel_intra["iso3_o"].map(areas)
    # Head-Mayer baseline: d_ii_HM = 0.67 * sqrt(area/π)
    panel_intra["d_ii_HM"] = 0.67 * np.sqrt(panel_intra["area_km2"] / np.pi)

    # Try to load our d_ii^eff measure from notebook 13
    if DII_EFF_CSV.exists():
        eff = pd.read_csv(DII_EFF_CSV)
        # Rename ISO column to align
        eff = eff.rename(columns={"iso3": "iso3_o"})
        if "d_theta_-5.0" in eff.columns:
            eff = eff[["iso3_o", "d_theta_-5.0"]].rename(columns={"d_theta_-5.0": "d_ii_eff_neg5"})
            panel_intra = panel_intra.merge(eff, on="iso3_o", how="left")
            n_with_eff = panel_intra["d_ii_eff_neg5"].notna().sum()
            print(f"      Merged d_ii^eff(θ=-5) for {n_with_eff}/{len(panel_intra)} countries")
        else:
            print(f"      WARNING: d_theta_-5.0 column not in {DII_EFF_CSV}")
            panel_intra["d_ii_eff_neg5"] = np.nan
    else:
        print(f"      d_ii^eff CSV not yet produced (waiting on notebook 13). Using HM only.")
        panel_intra["d_ii_eff_neg5"] = np.nan

    return panel_inter, panel_intra


def estimate_gravity(panel_inter, panel_intra, *, dist_intra_col: str, label: str):
    """Estimate log-OLS gravity with FE and a home-bias indicator."""
    import statsmodels.formula.api as smf

    # Bilateral distance: CEPII distw_harmonic for inter, dist_intra_col for intra
    inter = panel_inter.copy()
    inter["d_ij_km"] = inter["distw_harmonic"]
    intra = panel_intra.copy()
    intra["d_ij_km"] = intra[dist_intra_col]
    intra["distw_harmonic"] = intra[dist_intra_col]  # fill for consistency

    panel = pd.concat([inter, intra], ignore_index=True, sort=False)
    panel = panel[(panel["X_ij_usd"] > 0) &
                  (panel["d_ij_km"] > 0) &
                  (panel["gdp_o"] > 0) &
                  (panel["gdp_d"] > 0)].copy()
    panel = panel.dropna(subset=["d_ij_km"])
    if len(panel) == 0:
        print(f"  {label}: NO USABLE OBSERVATIONS")
        return None

    panel["log_X"] = np.log(panel["X_ij_usd"])
    panel["log_d"] = np.log(panel["d_ij_km"])
    panel["log_gdp_o"] = np.log(panel["gdp_o"])
    panel["log_gdp_d"] = np.log(panel["gdp_d"])

    print(f"\n  {label}")
    print(f"    Obs: {len(panel):,} ({panel['home'].sum():,} intra-national)")

    # Simple log-OLS with exporter and importer FE
    # (use 1.0 instead of large dummies for cleanliness in this pilot)
    formula = "log_X ~ log_d + log_gdp_o + log_gdp_d + home + C(iso3_o) + C(iso3_d)"
    try:
        model = smf.ols(formula, data=panel).fit(
            cov_type="cluster", cov_kwds={"groups": panel["iso3_o"] + "_" + panel["iso3_d"]}
        )
    except Exception as e:
        print(f"    ESTIMATION FAILED: {e}")
        return None

    # Extract key coefficients
    coef_dist = model.params.get("log_d", np.nan)
    se_dist = model.bse.get("log_d", np.nan)
    coef_home = model.params.get("home", np.nan)
    se_home = model.bse.get("home", np.nan)
    r2 = model.rsquared

    print(f"    log_d         : {coef_dist:+.4f}  (se {se_dist:.4f})")
    print(f"    home          : {coef_home:+.4f}  (se {se_home:.4f})")
    print(f"    Home-bias mult: exp(home) = {np.exp(coef_home):.2f}x")
    print(f"    R²            : {r2:.4f}")

    return {
        "label": label,
        "n_obs": len(panel),
        "n_intra": int(panel["home"].sum()),
        "coef_dist": coef_dist, "se_dist": se_dist,
        "coef_home": coef_home, "se_home": se_home,
        "home_bias_mult": float(np.exp(coef_home)),
        "r2": r2,
    }


def main():
    print(f"=== Border-effect pilot — year {PILOT_YEAR} ===")
    panel_inter, panel_intra = assemble_panel(PILOT_YEAR)

    # Run spec A: Head-Mayer closed-form for intra-national
    print("\n" + "="*70)
    print("SPEC A (BASELINE): d_ii = 0.67 * sqrt(area/π)  [Head-Mayer closed form]")
    res_A = estimate_gravity(panel_inter, panel_intra,
                              dist_intra_col="d_ii_HM",
                              label="Spec A — HM closed form")

    # Run spec B: our d_ii^eff(θ=-5) — only if available
    if panel_intra["d_ii_eff_neg5"].notna().sum() < 10:
        print("\n  d_ii^eff(θ=-5) not yet available for enough countries — Spec B skipped.")
        print("  Run notebook 13 first to populate real_country_internal_distance.csv")
        return

    print("\n" + "="*70)
    print("SPEC B (TREATMENT): d_ii = d_ii^eff(θ=-5)  [raster-resolution CES measure]")
    # Use d_ii^eff where available; fall back to HM where missing
    panel_intra["d_ii_combined"] = panel_intra["d_ii_eff_neg5"].fillna(panel_intra["d_ii_HM"])
    res_B = estimate_gravity(panel_inter, panel_intra,
                              dist_intra_col="d_ii_combined",
                              label="Spec B — raster d_ii^eff(θ=-5)")

    # Compare
    print("\n" + "="*70)
    print("COMPARISON")
    if res_A and res_B:
        delta_home = res_B["coef_home"] - res_A["coef_home"]
        delta_mult = res_B["home_bias_mult"] / res_A["home_bias_mult"]
        print(f"  Home-bias log-coef: A={res_A['coef_home']:+.4f}  B={res_B['coef_home']:+.4f}")
        print(f"  Δ(home-bias log-coef): {delta_home:+.4f}")
        print(f"  Home-bias multiplier: A={res_A['home_bias_mult']:.2f}x  B={res_B['home_bias_mult']:.2f}x")
        print(f"  Ratio of multipliers: {delta_mult:.3f}")
        if abs(delta_home) > 0.10:
            direction = "shrunk" if delta_home < 0 else "grew"
            print(f"\n  → The home-bias coefficient {direction} meaningfully under raster d_ii.")
            print(f"     {abs(delta_home):.2f} log-points reduction = {(1-delta_mult)*100:+.0f}% in multiplier terms.")
        else:
            print(f"\n  → No meaningful change in home-bias coefficient.")
            print(f"     Implication: sub-national heterogeneity in d_ii alone doesn't fix the border puzzle.")

    # Restricted comparison: identical intra sample for both specs (only the
    # countries that actually have a raster measure), so the home-coefficient
    # difference is driven purely by the d_ii swap, not by sample composition.
    print("\n" + "="*70)
    print("RESTRICTED COMPARISON (intra sample = countries with d_ii^eff only)")
    intra_eff = panel_intra[panel_intra["d_ii_eff_neg5"].notna()].copy()
    res_A_r = estimate_gravity(panel_inter, intra_eff,
                               dist_intra_col="d_ii_HM",
                               label="Spec A restricted — HM closed form")
    res_B_r = estimate_gravity(panel_inter, intra_eff,
                               dist_intra_col="d_ii_eff_neg5",
                               label="Spec B restricted — raster d_ii^eff(θ=-5)")
    if res_A_r and res_B_r:
        d_home = res_B_r["coef_home"] - res_A_r["coef_home"]
        print(f"\n  Restricted Δ(home-bias log-coef): {d_home:+.4f} "
              f"on {res_A_r['n_intra']} shared intra obs")


if __name__ == "__main__":
    main()
