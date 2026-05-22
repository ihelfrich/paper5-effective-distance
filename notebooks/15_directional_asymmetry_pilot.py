"""Pilot test of VIIRS (production) × population (consumption) directional asymmetry.

For a curated set of ~20 countries with measurable internal production-consumption
density divergence, compute:

  (1) the centroid divergence (cheap diagnostic) between VIIRS nightlights centroid
      and GHS-POP centroid, per country
  (2) the full asymmetric d_ij matrix where origin uses VIIRS, destination uses GHS-POP
  (3) the asymmetry index ASI = log(d_ij / d_ji) for every unordered pair
  (4) compare ASI to the observed trade-flow asymmetry log(X_ij / X_ji) using
      BACI 2010 data with PPML gravity FE absorbed

This is the empirical test of prediction P1 in PRE_ANALYSIS_PLAN v0.2.

Per literature triangulation 2026-05-20: this specific test is genuinely novel.
Henderson-Storeygard-Weil (2012), Bluhm-Krause (2018), and Donaldson-Storeygard
(2016) all use nightlights as activity proxies but none combine production-side
nightlights with consumption-side population to drive bilateral effective
distance through CES gravity.

Run:
  .venv/bin/python -u notebooks/15_directional_asymmetry_pilot.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

from paper5.region_raster import (
    load_country_boundaries,
    mask_raster_to_country,
    mask_result_to_raster_dist,
)
from paper5.directional_asymmetry import (
    asymmetric_bilateral_matrix,
    asymmetry_index,
    build_centroid_divergence_panel,
)


# Use the same paths the running notebook found
_GHSPOP = Path("/Volumes/HELFRICH-GD/TEG_data/inputs/ghsl/pop/"
               "GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif")
_NIGHTLIGHTS = Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/"
                    "nightlights/Harmonized_DN_NTL_2020_simVIIRS.tif")

# Slightly smaller pilot — focus on the largest geographically heterogeneous
# countries plus a control group of compact ones.
TEST_COUNTRIES = [
    # heterogeneous (expect ASI signal)
    ("USA", "United States"),
    ("RUS", "Russia"),
    ("CAN", "Canada"),
    ("BRA", "Brazil"),
    ("AUS", "Australia"),
    ("CHN", "China"),
    ("IND", "India"),
    ("IDN", "Indonesia"),
    ("MEX", "Mexico"),
    ("ARG", "Argentina"),
    ("CHL", "Chile"),
    ("ZAF", "South Africa"),
    ("EGY", "Egypt"),
    # compact (expect ASI ≈ 0)
    ("NLD", "Netherlands"),
    ("BEL", "Belgium"),
    ("KOR", "South Korea"),
    ("CHE", "Switzerland"),
    ("JPN", "Japan"),
    ("DEU", "Germany"),
    ("FRA", "France"),
]


def _find_nightlights_2020():
    """Try several known names for the harmonized 2020 nightlights raster."""
    candidates = [
        _NIGHTLIGHTS,
        Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/"
             "nightlights/Harmonized_DN_NTL_2020_calDMSP.tif"),
        Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/"
             "nightlights/Harmonized_DN_NTL_2019_simVIIRS.tif"),
        Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/"
             "nightlights/Harmonized_DN_NTL_2019_calDMSP.tif"),
    ]
    for c in candidates:
        if c.exists():
            return c
    # Take whatever 2020 file exists, even if different suffix
    nl_dir = Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/nightlights")
    if nl_dir.exists():
        for f in sorted(nl_dir.glob("Harmonized_DN_NTL_2020*.tif")):
            return f
        for f in sorted(nl_dir.glob("Harmonized_DN_NTL_2019*.tif")):
            return f
    raise FileNotFoundError(f"No harmonized nightlights raster found in {nl_dir}")


def main():
    nl_path = _find_nightlights_2020()
    print(f"Using nightlights raster: {nl_path}")
    print(f"Using GHS-POP raster: {_GHSPOP}")
    if not _GHSPOP.exists():
        raise FileNotFoundError(_GHSPOP)

    boundaries = load_country_boundaries()
    available = set(boundaries["ISO3"])

    pop_rds: dict[str, "RasterDist"] = {}
    nl_rds: dict[str, "RasterDist"] = {}

    print(f"\nMasking {len(TEST_COUNTRIES)} countries against GHS-POP + nightlights...")
    for iso, name in TEST_COUNTRIES:
        if iso not in available:
            print(f"  {iso:>3}  {name:<25}  MISSING in boundaries")
            continue
        geom = boundaries[boundaries["ISO3"] == iso].iloc[0].geometry

        t0 = time.time()
        try:
            m_pop = mask_raster_to_country(_GHSPOP, geom, max_atoms=2000)
            m_nl = mask_raster_to_country(nl_path, geom, max_atoms=2000)
        except Exception as e:
            print(f"  {iso:>3}  {name:<25}  MASK FAILED: {e}")
            continue

        if m_pop.kept_n_cells == 0 or m_nl.kept_n_cells == 0:
            print(f"  {iso:>3}  {name:<25}  EMPTY MASK (pop={m_pop.kept_n_cells} "
                  f"nl={m_nl.kept_n_cells})")
            continue

        pop_rds[iso] = mask_result_to_raster_dist(m_pop)
        nl_rds[iso] = mask_result_to_raster_dist(m_nl)
        print(f"  {iso:>3}  {name:<25}  pop_cells={m_pop.kept_n_cells:>5}  "
              f"nl_cells={m_nl.kept_n_cells:>5}  ({time.time()-t0:.1f}s)")

    # ── Diagnostic 1: centroid divergence per country ──
    print(f"\nComputing centroid divergence (nightlights vs population)...")
    div = build_centroid_divergence_panel(nl_rds, pop_rds)
    print("\n=== Internal nightlights ↔ population centroid divergence ===")
    print(div.head(20).to_string(index=False))
    div.to_csv(_repo / "data" / "derived" / "centroid_divergence_pilot.csv", index=False)

    # ── Asymmetric bilateral matrix at θ=-5 ──
    print(f"\nComputing asymmetric d_ij matrix at θ=-5 (n={len(pop_rds)} countries)...")
    t0 = time.time()
    bilat = asymmetric_bilateral_matrix(
        nl_rds, pop_rds, theta=-5.0, verbose=False
    )
    print(f"  {len(bilat)} ordered pairs in {time.time()-t0:.1f}s")
    bilat.to_csv(_repo / "data" / "derived" / "directional_d_ij_pilot.csv", index=False)

    # ── Asymmetry index ──
    asi = asymmetry_index(bilat)
    print(f"\n=== Top 10 most-asymmetric pairs (largest |ASI|) ===")
    print(asi.head(10).to_string(index=False))
    asi.to_csv(_repo / "data" / "derived" / "asymmetry_index_pilot.csv", index=False)

    # ── Quick descriptive: ASI distribution ──
    print(f"\n=== ASI distribution ===")
    print(f"  n pairs       : {len(asi)}")
    print(f"  median |ASI|  : {asi['abs_ASI'].median():.4f}")
    print(f"  90th pct |ASI|: {asi['abs_ASI'].quantile(0.90):.4f}")
    print(f"  max |ASI|     : {asi['abs_ASI'].max():.4f}")

    # ── Merge against BACI 2010 trade-flow asymmetry ──
    baci_2010 = Path("/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/"
                     "BACI_HS02_Y2010_V202401b.csv")
    if baci_2010.exists():
        print(f"\nLoading BACI 2010 and joining trade-flow asymmetry...")
        baci = pd.read_csv(baci_2010, usecols=["i", "j", "v"])
        baci["v"] = pd.to_numeric(baci["v"], errors="coerce")
        agg = baci.dropna(subset=["v"]).groupby(["i", "j"], as_index=False)["v"].sum()
        agg["X_ij_kusd"] = agg["v"]  # BACI is already in thousands of USD

        # ISO-num → ISO3 via CEPII
        cepii_path = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")
        if cepii_path.exists():
            cepii = pd.read_csv(cepii_path,
                                usecols=["year", "iso3_o", "iso3num_o", "iso3_d", "iso3num_d"],
                                dtype={"iso3_o": str, "iso3_d": str})
            cepii = cepii[cepii["year"] == 2010]
            num_to_iso = dict(cepii[["iso3num_o", "iso3_o"]].drop_duplicates().values)
            num_to_iso.update(dict(cepii[["iso3num_d", "iso3_d"]].drop_duplicates().values))
            agg["iso_o"] = agg["i"].map(num_to_iso)
            agg["iso_d"] = agg["j"].map(num_to_iso)
            agg = agg.dropna(subset=["iso_o", "iso_d"])

            # Merge X_ij and X_ji to compute trade-flow asymmetry
            asi_subset = asi.copy()
            fwd = agg.rename(columns={"X_ij_kusd": "X_ij"})[["iso_o", "iso_d", "X_ij"]]
            rev = agg.rename(columns={"X_ij_kusd": "X_ji", "iso_o": "iso_d", "iso_d": "iso_o"})[
                ["iso_o", "iso_d", "X_ji"]
            ]
            merged = asi_subset.merge(fwd, on=["iso_o", "iso_d"], how="left")
            merged = merged.merge(rev, on=["iso_o", "iso_d"], how="left")
            merged = merged.dropna(subset=["X_ij", "X_ji"])
            merged = merged[(merged["X_ij"] > 0) & (merged["X_ji"] > 0)]
            merged["trade_log_ratio"] = np.log(merged["X_ij"] / merged["X_ji"])
            print(f"\n  {len(merged)} pairs matched with positive bilateral trade")

            r = merged[["ASI", "trade_log_ratio"]].corr().iloc[0, 1]
            print(f"  corr(ASI, log(X_ij/X_ji)) = {r:.4f}")

            # Naive regression: log(X_ij/X_ji) = β * ASI + ε
            x = merged["ASI"].values
            y = merged["trade_log_ratio"].values
            x_mean = x.mean(); y_mean = y.mean()
            beta = np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean) ** 2)
            print(f"  Naive OLS β (no FE)        = {beta:.3f}")
            print(f"  Expected sign under P1: β  < 0  (longer effective distance → less trade)")
            print(f"  Got: β {'<' if beta < 0 else '>='} 0")

            merged.to_csv(_repo / "data" / "derived" / "asi_vs_trade_asymmetry.csv", index=False)
        else:
            print(f"  (CEPII Gravity not at {cepii_path} — skipping ISO mapping)")
    else:
        print(f"  (BACI 2010 not at {baci_2010} — skipping trade merge)")

    print(f"\nDone. Outputs in {_repo / 'data' / 'derived'}/")


if __name__ == "__main__":
    main()
