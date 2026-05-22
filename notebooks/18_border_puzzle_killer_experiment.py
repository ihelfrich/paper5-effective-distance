"""Killer experiment: how much does the border-puzzle home-bias coefficient
shrink when we swap CEPII centroid distance for the structurally-consistent
CES effective distance for contiguous-country pairs?

The hypothesis from the harmonic-limit / contact-set degeneracy paper:
CEPII's distwces uses θ = -1, which is the reduced-form gravity slope, not
the structural δ(1-σ) ∈ [-7, -3]. For contiguous pairs (those with land
borders), the structurally-correct distance is two-to-three orders of
magnitude smaller than CEPII. Substituting this in should cause the
home-bias coefficient to fall, because intra-national distance is no longer
systematically smaller than inter-national distance for adjacent partners.

Pre-analysis specification:
  1. Build the standard 2010 gravity panel: BACI bilateral flows + CEPII
     covariates + intra-national rows with X_ii = max(0, GDP - exports).
  2. Estimate the baseline gravity: log X_ij = const + β log distw_harmonic
     + γ controls + home + exporter-year FE + importer-year FE.
  3. Replace distw_harmonic with d_eff(θ=-5) on the 34 contiguous pairs we
     have from nb13 v3. Re-estimate. Compare the home-bias coefficient.
  4. Robustness: also drop contiguous pairs entirely from the sample
     (extreme treatment).

Run:
  .venv/bin/python -u notebooks/18_border_puzzle_killer_experiment.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR = 2010
PILOT_COUNTRIES = [
    "USA","RUS","CAN","BRA","AUS","CHN","IND","IDN","JPN","FRA","DEU",
    "NLD","BEL","KOR","CHE","GBR","ITA","ESP","MEX","ARG","ZAF","EGY",
    "CHL","NOR","SWE","FIN","SAU",
]


def load_data():
    """Load BACI 2010, CEPII Gravity 2010, and our d_eff(θ=-5) bilateral matrix."""
    print("[1/5] Loading BACI 2010 (chunked, aggregating to country-pair)...")
    baci_path = "/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/BACI_HS02_Y2010_V202401b.csv"
    agg_partial = []
    for chunk in pd.read_csv(baci_path, usecols=["i", "j", "v"],
                             dtype={"i": int, "j": int},
                             chunksize=1_000_000):
        chunk["v"] = pd.to_numeric(chunk["v"], errors="coerce")
        chunk = chunk.dropna(subset=["v"])
        agg_chunk = chunk.groupby(["i", "j"], as_index=False)["v"].sum()
        agg_partial.append(agg_chunk)
    agg_all = pd.concat(agg_partial, ignore_index=True)
    agg = agg_all.groupby(["i", "j"], as_index=False)["v"].sum()
    agg["X_ij_kusd"] = agg["v"]
    agg["X_ij_usd"] = agg["X_ij_kusd"] * 1000
    print(f"  {len(agg):,} country-pair flows")

    print("[2/5] Loading CEPII Gravity 2010 (chunked, year-filtered)...")
    cepii_path = "/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv"
    chunks = []
    for chunk in pd.read_csv(cepii_path, dtype={"iso3_o": str, "iso3_d": str, "year": "Int64"},
                             chunksize=500_000, low_memory=False):
        sub = chunk[chunk["year"] == YEAR]
        if len(sub) > 0:
            chunks.append(sub)
    cepii = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    keep = ["iso3_o", "iso3_d", "iso3num_o", "iso3num_d",
            "gdp_o", "gdp_d", "dist", "distcap",
            "distw_arithmetic", "distw_harmonic",
            "contig", "comlang_off", "col_dep_ever", "rta"]
    keep = [c for c in keep if c in cepii.columns]
    cepii = cepii[keep]
    print(f"  {len(cepii):,} CEPII pair-year rows")

    # ISO num to ISO3 mapping — drop NaN keys/values
    map_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates(subset=["iso3num_o"])
    map_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates(subset=["iso3num_d"])
    num_to_iso = {int(n): s for n, s in map_o.values}
    num_to_iso.update({int(n): s for n, s in map_d.values})
    agg["iso_o"] = agg["i"].map(num_to_iso)
    agg["iso_d"] = agg["j"].map(num_to_iso)
    agg = agg.dropna(subset=["iso_o", "iso_d"])

    print("[3/5] Merging BACI flows onto CEPII covariates...")
    inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]].rename(
        columns={"iso3_o": "iso_o", "iso3_d": "iso_d"}
    )
    panel = inter.merge(agg[["iso_o", "iso_d", "X_ij_usd"]],
                        on=["iso_o", "iso_d"], how="left")
    panel["X_ij_usd"] = panel["X_ij_usd"].fillna(0.0)
    panel["home"] = 0
    print(f"  Inter-national panel: {len(panel):,} rows ({(panel['X_ij_usd'] > 0).sum():,} positive)")

    print("[4/5] Adding intra-national flows X_ii = max(0, GDP - exports)...")
    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exports_by_origin = agg.groupby("iso_o", as_index=False)["X_ij_usd"].sum()
    exports_by_origin = exports_by_origin.rename(columns={"X_ij_usd": "exports"})
    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        if not np.isfinite(gdp) or gdp <= 0:
            continue
        exp = exports_by_origin.loc[exports_by_origin["iso_o"] == iso, "exports"]
        exp_v = float(exp.iloc[0]) if len(exp) > 0 else 0.0
        x_ii = max(0.0, gdp - exp_v)
        intra_rows.append({
            "iso_o": iso, "iso_d": iso, "X_ij_usd": x_ii,
            "gdp_o": gdp, "gdp_d": gdp,
            "contig": 0, "comlang_off": 1, "col_dep_ever": 0, "rta": 0,
            "home": 1,
            # Head-Mayer closed form using geographic area
            # For now use CEPII's distw_harmonic which CEPII populates with the
            # closed-form for i=j (verified by spot check)
        })
    intra = pd.DataFrame(intra_rows)
    # Pull distw_harmonic for i=j from CEPII (it stores Head-Mayer there)
    cepii_self = cepii[cepii["iso3_o"] == cepii["iso3_d"]].rename(
        columns={"iso3_o": "iso_o"})[["iso_o", "distw_harmonic", "distw_arithmetic"]]
    intra = intra.merge(cepii_self, on="iso_o", how="left")
    intra["iso_d"] = intra["iso_o"]
    print(f"  Intra-national rows: {len(intra)}, with distw: {intra['distw_harmonic'].notna().sum()}")

    print("[5/5] Loading d_eff(θ=-5) bilateral from nb13 v3...")
    bilat = pd.read_csv(_repo / "data" / "derived" / "real_country_bilateral_distance.csv")
    bilat = bilat.rename(columns={"d_eff": "d_eff_neg5"})
    print(f"  {len(bilat)} ordered pairs in d_eff matrix (27-country pilot)")

    panel_full = pd.concat([panel, intra], ignore_index=True, sort=False)
    print(f"\nFinal panel: {len(panel_full)} rows ({panel_full['home'].sum()} intra-national)")
    return panel_full, bilat


def estimate_specs(panel: pd.DataFrame, bilat: pd.DataFrame):
    """Estimate four gravity specs. All use OLS on log-positive flows.

    Spec A (BASELINE): all distances from CEPII's distw_harmonic.
    Spec B (D_EFF FOR PILOT BILATERAL): substitute d_eff(θ=-5) on the 27-pilot pairs
       where we have it; keep CEPII for all other pairs.
    Spec C (D_EFF FOR CONTIG ONLY): substitute d_eff(θ=-5) ONLY on contiguous pairs
       in the pilot; keep CEPII for everything else (the cleanest test of the
       harmonic-limit prediction).
    Spec D (DROP CONTIG): drop contiguous pairs entirely from the sample.
    """
    import statsmodels.formula.api as smf

    # Attach d_eff
    bilat_keep = bilat.rename(columns={"iso_o": "iso_o", "iso_d": "iso_d"})
    panel = panel.merge(bilat_keep[["iso_o", "iso_d", "d_eff_neg5"]],
                        on=["iso_o", "iso_d"], how="left")

    print(f"\nPanel rows with d_eff_neg5: {panel['d_eff_neg5'].notna().sum()}")
    print(f"  intra-national rows with d_eff: {((panel['d_eff_neg5'].notna()) & (panel['home']==1)).sum()}")
    print(f"  inter-national rows with d_eff: {((panel['d_eff_neg5'].notna()) & (panel['home']==0)).sum()}")

    results = {}
    for spec_name in ["A_baseline", "B_d_eff_all_pilot", "C_d_eff_contig_only", "D_drop_contig"]:
        df = panel.copy()
        df["distance"] = df["distw_harmonic"]  # baseline default

        if spec_name == "A_baseline":
            pass
        elif spec_name == "B_d_eff_all_pilot":
            # Use d_eff wherever we have it (the pilot subset)
            df.loc[df["d_eff_neg5"].notna(), "distance"] = df["d_eff_neg5"]
        elif spec_name == "C_d_eff_contig_only":
            # Use d_eff ONLY on contiguous pairs
            contig_mask = (df["d_eff_neg5"].notna()) & (df["contig"] == 1)
            df.loc[contig_mask, "distance"] = df["d_eff_neg5"]
        elif spec_name == "D_drop_contig":
            df = df[df["contig"] != 1].copy()

        # Filter to estimation sample
        df = df[(df["X_ij_usd"] > 0) &
                (df["distance"] > 0) &
                (df["gdp_o"] > 0) & (df["gdp_d"] > 0)].copy()
        df["log_X"] = np.log(df["X_ij_usd"])
        df["log_d"] = np.log(df["distance"])
        df["log_gdp_o"] = np.log(df["gdp_o"])
        df["log_gdp_d"] = np.log(df["gdp_d"])

        # Estimate: log X ~ log_d + log_gdp_o + log_gdp_d + home + C(iso_o) + C(iso_d)
        # The exporter and importer fixed effects absorb multilateral resistance.
        try:
            formula = "log_X ~ log_d + home + C(iso_o) + C(iso_d)"
            model = smf.ols(formula, data=df).fit()
            coef_d = float(model.params.get("log_d", np.nan))
            se_d = float(model.bse.get("log_d", np.nan))
            coef_home = float(model.params.get("home", np.nan))
            se_home = float(model.bse.get("home", np.nan))
            home_mult = float(np.exp(coef_home))
            r2 = float(model.rsquared)
            n = int(model.nobs)
            n_intra = int(df["home"].sum())
            print(f"\n[{spec_name}]")
            print(f"  N = {n:,} (intra = {n_intra}), R² = {r2:.3f}")
            print(f"  log_d  = {coef_d:+.4f}  (SE {se_d:.4f})")
            print(f"  home   = {coef_home:+.4f}  (SE {se_home:.4f})   → multiplier {home_mult:.2f}x")
            results[spec_name] = {
                "n": n, "n_intra": n_intra, "r2": r2,
                "coef_d": coef_d, "se_d": se_d,
                "coef_home": coef_home, "se_home": se_home,
                "home_mult": home_mult,
            }
        except Exception as e:
            print(f"\n[{spec_name}]  ESTIMATION FAILED: {e}")
            results[spec_name] = {"error": str(e)}

    return results


def main():
    panel, bilat = load_data()
    results = estimate_specs(panel, bilat)

    out = pd.DataFrame(results).T
    out_path = OUT_DIR / "border_puzzle_killer_experiment.csv"
    out.to_csv(out_path)
    print(f"\nSaved to {out_path}")

    if "A_baseline" in results and "C_d_eff_contig_only" in results:
        a = results["A_baseline"]
        c = results["C_d_eff_contig_only"]
        if "coef_home" in a and "coef_home" in c:
            delta = c["coef_home"] - a["coef_home"]
            print(f"\n=== Border-puzzle shrinkage from contiguous-pair distance correction ===")
            print(f"  Spec A (CEPII baseline)        home = {a['coef_home']:+.4f}")
            print(f"  Spec C (d_eff on contig only)  home = {c['coef_home']:+.4f}")
            print(f"  Δ home coefficient             = {delta:+.4f}")
            print(f"  Multiplier change              = {a['home_mult']:.2f}x  →  {c['home_mult']:.2f}x")
            ratio = c['home_mult'] / a['home_mult']
            print(f"  Home-bias multiplier ratio     = {ratio:.4f}")
            if delta < -0.10:
                print(f"  → Home bias SHRUNK by {-delta*100:.0f} log points.")
            elif delta > 0.10:
                print(f"  → Home bias GREW by {delta*100:.0f} log points (unexpected).")
            else:
                print(f"  → Home bias roughly unchanged.")


if __name__ == "__main__":
    main()
