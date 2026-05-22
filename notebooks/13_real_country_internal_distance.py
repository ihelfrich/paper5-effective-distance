"""First real-data test of d_ii^eff vs. Head-Mayer closed-form.

The synthetic simulations (notebooks 10, 11, 12) used shape primitives. This
notebook uses ACTUAL country polygons (Natural Earth admin0) and ACTUAL
population rasters (GHS-POP 2020) for a curated set of heterogeneous countries.

The hypothesis: for countries with strongly non-uniform internal geography
(US, Russia, Brazil, Australia, Canada, China), d_ii^eff at the
gravity-consistent θ=1-σ differs *meaningfully* from the Head-Mayer
0.67*sqrt(area/π) closed-form. For roughly-uniform / small countries
(Netherlands, Belgium, Korea), they should agree closely.

This is the first empirical confirmation that the raster-resolution
extension of Head-Mayer 2014 is worth running on the full panel.

Run:
  .venv/bin/python notebooks/13_real_country_internal_distance.py
"""

from __future__ import annotations

import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

from paper5.ot_distance import ces_effective_distance, head_mayer_closed_form
from paper5.region_raster import (
    load_country_boundaries,
    mask_raster_to_country,
    mask_result_to_raster_dist,
)


_GHSPOP_DIR = Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/ghs_pop/")
# Pick the most-recent extracted year automatically
_candidates = [_GHSPOP_DIR / f"GHS_POP_E{y}_GLOBE_R2023A_54009_1000_V1_0.tif"
               for y in (2020, 2015, 2010, 2005, 2000)]
GHS_POP_2020 = next((p for p in _candidates if p.exists()), _candidates[0])

# Curated test set: a mix of geographically heterogeneous and homogeneous countries
TEST_COUNTRIES = [
    # ISO3, name, expectation
    ("USA", "United States", "very heterogeneous (coastal concentrations)"),
    ("RUS", "Russia", "extremely heterogeneous (Moscow-west, Vladivostok-east)"),
    ("CAN", "Canada", "very heterogeneous (coastal/southern strip)"),
    ("BRA", "Brazil", "coastal concentration"),
    ("AUS", "Australia", "extreme coastal concentration"),
    ("CHN", "China", "east-coast concentration"),
    ("IND", "India", "moderately heterogeneous"),
    ("IDN", "Indonesia", "archipelagic"),
    ("JPN", "Japan", "moderately heterogeneous"),
    ("FRA", "France", "fairly compact"),
    ("DEU", "Germany", "fairly compact"),
    ("NLD", "Netherlands", "small, compact"),
    ("BEL", "Belgium", "small, compact"),
    ("KOR", "South Korea", "compact"),
    ("CHE", "Switzerland", "small, mountainous"),
    ("GBR", "United Kingdom", "moderate"),
    ("ITA", "Italy", "elongated peninsula"),
    ("ESP", "Spain", "moderate"),
    ("MEX", "Mexico", "elongated"),
    ("ARG", "Argentina", "elongated, southern population"),
    ("ZAF", "South Africa", "moderate"),
    ("EGY", "Egypt", "Nile-corridor concentration"),
    ("CHL", "Chile", "extremely elongated"),
    ("NOR", "Norway", "elongated, coastal"),
    ("SWE", "Sweden", "elongated, southern"),
    ("FIN", "Finland", "elongated, southern"),
    ("SAU", "Saudi Arabia", "interior emptiness"),
]


def compute_bilateral_matrix(rds: dict[str, "RasterDist"], theta: float) -> pd.DataFrame:
    """Compute pairwise d_ij^eff for every pair in the given dict of RasterDist."""
    isos = sorted(rds.keys())
    rows = []
    for i, iso_i in enumerate(isos):
        for iso_j in isos[i:]:  # exploit symmetry of the measure under symmetric weighting
            d = ces_effective_distance(rds[iso_i], rds[iso_j], theta=theta)
            rows.append({"iso_o": iso_i, "iso_d": iso_j, "d_eff": d})
            if iso_i != iso_j:
                rows.append({"iso_o": iso_j, "iso_d": iso_i, "d_eff": d})
    return pd.DataFrame(rows)


def main():
    print(f"Loading country boundaries...")
    boundaries = load_country_boundaries()
    print(f"  {len(boundaries)} countries")

    available_isos = set(boundaries["ISO3"].tolist())
    missing = [iso for iso, _, _ in TEST_COUNTRIES if iso not in available_isos]
    if missing:
        print(f"  WARNING: missing ISO3 codes: {missing}")

    print(f"\nProcessing {len(TEST_COUNTRIES)} countries against GHS-POP 2020...")
    if not GHS_POP_2020.exists():
        raise FileNotFoundError(f"GHS-POP raster not at {GHS_POP_2020}")

    THETAS = [1.0, -1.0, -3.0, -5.0, -7.0]
    results = []

    # Precompute geographic country areas (Mollweide equal-area)
    print(f"Precomputing geographic country areas (Mollweide projection)...")
    boundaries_eq = boundaries.to_crs("ESRI:54009")
    geo_area_lookup = dict(zip(
        boundaries_eq["ISO3"],
        (boundaries_eq.geometry.area / 1e6).astype(float).values
    ))

    for iso, name, expectation in TEST_COUNTRIES:
        if iso not in available_isos:
            continue
        row = boundaries[boundaries["ISO3"] == iso].iloc[0]
        geom = row.geometry

        t0 = time.time()
        # Block coarsening at 5km preserves total population and spatial
        # coverage; reduces atom count ~25x. Then cap at max_atoms=10000 as
        # safety. This is the methodology fix for the harmonic-dominance
        # bias documented in METHODOLOGY_NOTES.md §1.
        try:
            mask = mask_raster_to_country(
                GHS_POP_2020, geom, coarsen_factor=5, max_atoms=10000
            )
        except Exception as e:
            print(f"  {iso:>3}  {name:<25}  MASK FAILED: {e}")
            continue
        t_mask = time.time() - t0

        if mask.kept_n_cells == 0:
            print(f"  {iso:>3}  {name:<25}  NO POPULATION DATA (mask returned 0 cells)")
            continue

        rd = mask_result_to_raster_dist(mask)

        # Country area in km² — use GEOGRAPHIC country area (Natural Earth polygon
        # reprojected to Mollweide). The populated raster area (mask.full_area_km2)
        # is smaller because deserts/mountains have zero population; Head-Mayer's
        # closed form uses total country area, not populated area.
        geo_area_km2 = geo_area_lookup.get(iso, np.nan)
        pop_area_km2 = float(getattr(mask, "full_area_km2", 0.0) or mask.cell_areas_km2.sum())
        area_km2 = geo_area_km2 if np.isfinite(geo_area_km2) else pop_area_km2
        d_HM = head_mayer_closed_form(area_km2)

        t1 = time.time()
        d_thetas = {}
        for theta in THETAS:
            d_thetas[theta] = ces_effective_distance(rd, rd, theta=theta)
        t_ces = time.time() - t1

        result = {
            "iso3": iso, "name": name, "expectation": expectation,
            "area_km2": area_km2,  # geographic (Mollweide)
            "pop_area_km2": pop_area_km2,
            "n_cells": mask.kept_n_cells,
            "d_HM": d_HM,
            **{f"d_theta_{t}": d_thetas[t] for t in THETAS},
            "t_mask_s": t_mask, "t_ces_s": t_ces,
        }
        results.append(result)
        # Print row
        ratios = " ".join(f"{d_thetas[t]/d_HM:>5.2f}" for t in THETAS)
        print(f"  {iso:>3}  {name:<25}  cells={mask.kept_n_cells:>5}  "
              f"area={area_km2/1e6:>5.1f}Mkm²  HM={d_HM:>6.0f}km  "
              f"ratios(θ={THETAS}): {ratios}  ({t_mask:.1f}+{t_ces:.1f}s)")

    # Save intra-national results
    df = pd.DataFrame(results)
    out_csv = _repo / "data" / "derived" / "real_country_internal_distance.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nSaved {len(df)} rows to {out_csv}")

    # ── Cache RasterDist per country for bilateral computation ──
    # We'll re-mask cheaply since we already have everything in memory
    print(f"\nBuilding RasterDist cache for bilateral d_ij^eff matrix...")
    rd_cache = {}
    for iso, name, _ in TEST_COUNTRIES:
        if iso not in available_isos:
            continue
        if not any(r["iso3"] == iso for r in results):
            continue
        row = boundaries[boundaries["ISO3"] == iso].iloc[0]
        try:
            m = mask_raster_to_country(
                GHS_POP_2020, row.geometry, coarsen_factor=5, max_atoms=10000
            )
            rd_cache[iso] = mask_result_to_raster_dist(m)
        except Exception as e:
            print(f"  WARN: couldn't cache {iso}: {e}")
    print(f"  Cached {len(rd_cache)} country RasterDists")

    # Bilateral d_ij^eff at θ=-5 (gravity-consistent for σ=6)
    print(f"\nComputing bilateral d_ij^eff matrix at θ=-5...")
    t0 = time.time()
    bilat = compute_bilateral_matrix(rd_cache, theta=-5.0)
    print(f"  {len(bilat)} ordered pairs, {time.time()-t0:.1f}s")

    out_bilat = _repo / "data" / "derived" / "real_country_bilateral_distance.csv"
    bilat.to_csv(out_bilat, index=False)
    print(f"  Saved to {out_bilat}")

    # Quick sanity: compare to CEPII distw_harmonic for the same pairs
    try:
        cepii_path = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")
        if cepii_path.exists():
            print(f"\nComparing to CEPII distw_harmonic for the test set...")
            cepii = pd.read_csv(cepii_path, usecols=["year", "iso3_o", "iso3_d",
                                                     "distw_harmonic", "distw_arithmetic"])
            cepii = cepii[(cepii["year"] == 2010) &
                          (cepii["iso3_o"].isin(rd_cache.keys())) &
                          (cepii["iso3_d"].isin(rd_cache.keys())) &
                          (cepii["iso3_o"] != cepii["iso3_d"])]
            cepii = cepii.rename(columns={"iso3_o": "iso_o", "iso3_d": "iso_d"})
            merged = bilat[bilat["iso_o"] != bilat["iso_d"]].merge(
                cepii[["iso_o", "iso_d", "distw_harmonic", "distw_arithmetic"]],
                on=["iso_o", "iso_d"], how="inner"
            )
            merged = merged.dropna(subset=["distw_harmonic"])
            merged["ratio_eff_to_harmonic"] = merged["d_eff"] / merged["distw_harmonic"]
            merged["ratio_eff_to_arith"] = merged["d_eff"] / merged["distw_arithmetic"]
            print(f"  {len(merged)} pairs matched against CEPII")
            print(f"\n  Ratio d_eff(θ=-5) / CEPII distw_harmonic:")
            r = merged["ratio_eff_to_harmonic"].dropna()
            print(f"    median = {r.median():.3f}  mean = {r.mean():.3f}")
            print(f"    10%-90% = [{r.quantile(0.10):.3f}, {r.quantile(0.90):.3f}]")
            print(f"\n  Ratio d_eff(θ=-5) / CEPII distw_arithmetic:")
            r2 = merged["ratio_eff_to_arith"].dropna()
            print(f"    median = {r2.median():.3f}  mean = {r2.mean():.3f}")

            # Most-different pairs
            print(f"\n  Most-different pairs (largest |ratio - 1|):")
            merged["abs_dev"] = (merged["ratio_eff_to_harmonic"] - 1).abs()
            most_diff = merged.sort_values("abs_dev", ascending=False).head(8)
            for _, m in most_diff.iterrows():
                print(f"    {m['iso_o']} → {m['iso_d']}:  d_eff = {m['d_eff']:.0f} km,  "
                      f"CEPII = {m['distw_harmonic']:.0f} km,  ratio = {m['ratio_eff_to_harmonic']:.3f}")

            merged.to_csv(_repo / "data" / "derived" / "real_country_eff_vs_cepii.csv", index=False)
    except Exception as e:
        print(f"  Couldn't compare to CEPII: {e}")

    # Summary stats
    print("\n=== Summary statistics: d_theta / d_HM ratios ===")
    print(f"  {'theta':>6}  {'median':>8}  {'mean':>8}  {'min':>8}  {'max':>8}  "
          f"{'n diverge >20%':>14}")
    for theta in THETAS:
        col = f"d_theta_{theta}"
        if col not in df.columns:
            continue
        ratios = df[col] / df["d_HM"]
        n_diverge = ((ratios < 0.8) | (ratios > 1.2)).sum()
        print(f"  {theta:>6.1f}  {ratios.median():>8.3f}  {ratios.mean():>8.3f}  "
              f"{ratios.min():>8.3f}  {ratios.max():>8.3f}  {n_diverge:>14}/{len(df)}")

    # Concrete observation: most-different countries at θ=-5
    if "d_theta_-5.0" in df.columns:
        df["ratio_neg5"] = df["d_theta_-5.0"] / df["d_HM"]
        df_sorted = df.sort_values("ratio_neg5").reset_index(drop=True)
        print("\n=== Most-different countries at θ=-5 (where signal lives) ===")
        print(f"  {'iso':>3}  {'name':<25}  {'HM (km)':>10}  {'CES (km)':>10}  {'ratio':>6}  expectation")
        for _, r in df_sorted.head(8).iterrows():
            print(f"  {r['iso3']:>3}  {r['name']:<25}  {r['d_HM']:>10.0f}  "
                  f"{r['d_theta_-5.0']:>10.0f}  {r['ratio_neg5']:>6.3f}  {r['expectation']}")
        print("  ...")
        for _, r in df_sorted.tail(5).iterrows():
            print(f"  {r['iso3']:>3}  {r['name']:<25}  {r['d_HM']:>10.0f}  "
                  f"{r['d_theta_-5.0']:>10.0f}  {r['ratio_neg5']:>6.3f}  {r['expectation']}")


if __name__ == "__main__":
    main()
