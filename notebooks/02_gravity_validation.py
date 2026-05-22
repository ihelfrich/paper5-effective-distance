"""02_gravity_validation.py — gravity benchmark on 20-country prototype panel.

Run after 01_prototype_20country.py has saved agglomerations_prototype.parquet.
Requires BACI on HELFRICH-GD.

Spec selection (automatic):
  - If distance has within-pair time variation (real WorldPop centroids):
      PPML with iso_o^year + iso_d^year + iso_o^iso_d  (3-way FE, panel identification)
  - If centroids are static (fallback mode):
      PPML with iso_o^year + iso_d^year  (cross-section PPML; pair FE dropped to avoid
      perfect collinearity with time-invariant distance)

Outputs:
  - data/derived/gravity_benchmark_20c.csv   (beta, se, n_obs, spec per variant)
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from paper5.data_loaders import (
    HELFRICH_GD, PAPER5_DATA,
    load_baci_with_iso, load_cepii_gravity,
)
from paper5.distance import (
    Agglomeration, ChokepointState,
    compute_d_ovdl,
)

YEARS = [2010, 2022]
THETA = -1.0

CORE_20 = [
    "USA", "MEX", "CAN", "JAM",
    "DEU", "FRA", "POL",
    "NGA", "ETH", "ZAF",
    "SAU", "IND", "PAK",
    "CHN", "JPN", "VNM",
    "BRA", "COL",
    "AUS", "NZL",
]

FALLBACK_CENTROIDS = {
    "USA": (-95.71, 37.09, 331.0), "MEX": (-102.55, 23.63, 128.0),
    "CAN": (-96.80, 56.13, 38.0),  "JAM": (-77.30, 18.11, 3.0),
    "DEU": (10.45, 51.17, 83.0),   "FRA": (2.21, 46.23, 67.0),
    "POL": (19.15, 51.92, 38.0),   "NGA": (8.67, 9.08, 206.0),
    "ETH": (40.49, 9.15, 115.0),   "ZAF": (25.08, -29.00, 59.0),
    "SAU": (45.08, 23.89, 35.0),   "IND": (78.96, 20.59, 1380.0),
    "PAK": (69.35, 30.38, 220.0),  "CHN": (104.20, 35.86, 1411.0),
    "JPN": (138.25, 36.20, 126.0), "VNM": (108.28, 14.06, 97.0),
    "BRA": (-51.93, -14.24, 212.0), "COL": (-74.30, 4.57, 51.0),
    "AUS": (133.78, -25.27, 26.0), "NZL": (172.50, -41.50, 5.0),
}

# ── Step 1: Build d_ovdl_t panel ──────────────────────────────────────────────
print("Building d_ovdl_t panel for 20-country core (fallback centroids)...")

ovdl_rows = []
for year in YEARS:
    country_aggloms: dict[str, list[Agglomeration]] = {}
    for iso3, (lon, lat, pop) in FALLBACK_CENTROIDS.items():
        country_aggloms[iso3] = [
            Agglomeration(iso3=iso3, agglom_id=0, lon=lon, lat=lat,
                          pop=pop * 1e6, viirs=0.0, year=year)
        ]
    for iso_o in CORE_20:
        if iso_o not in country_aggloms:
            continue
        for iso_d in CORE_20:
            if iso_d == iso_o or iso_d not in country_aggloms:
                continue
            d = compute_d_ovdl(country_aggloms[iso_o], country_aggloms[iso_d],
                               theta=THETA, use_viirs=False)
            ovdl_rows.append({"iso_o": iso_o, "iso_d": iso_d, "year": year,
                              "d_ovdl_t": d})

ovdl_df = pd.DataFrame(ovdl_rows)
print(f"d_ovdl_t panel: {len(ovdl_df)} rows")

# ── Step 2: Load CEPII distw ───────────────────────────────────────────────────
print("Loading CEPII distw...")
try:
    cepii = load_cepii_gravity(filtered=False)
    cepii_dist = (cepii[cepii["iso3_o"].isin(CORE_20) & cepii["iso3_d"].isin(CORE_20)]
                  [["iso3_o", "iso3_d", "distw"]]
                  .rename(columns={"iso3_o":"iso_o","iso3_d":"iso_d","distw":"d_cepii_distw"})
                  .drop_duplicates(["iso_o","iso_d"]))
    print(f"CEPII pairs: {len(cepii_dist)}")
except Exception as e:
    print(f"CEPII load failed: {e}. Skipping.")
    cepii_dist = None

# ── Step 3: Load BACI trade flows ──────────────────────────────────────────────
print("Loading BACI trade flows...")
try:
    baci_frames = []
    for year in YEARS:
        try:
            df = load_baci_with_iso(year=year, aggregation="country-pair-year")
            df = df[df["iso_o"].isin(CORE_20) & df["iso_d"].isin(CORE_20)]
            baci_frames.append(df)
            print(f"  BACI {year}: {len(df)} pairs")
        except Exception as e:
            print(f"  BACI {year} failed: {e}")
    baci = pd.concat(baci_frames, ignore_index=True) if baci_frames else None
    if baci is not None:
        print(f"Total BACI: {len(baci)} rows, {baci['trade_value'].sum()/1e9:.1f}B USD")
except Exception as e:
    print(f"BACI load failed: {e}")
    baci = None

if baci is None:
    print("ERROR: BACI not available. Check HELFRICH-GD mount.")
    sys.exit(1)

# ── Step 4: Merge distance + trade ─────────────────────────────────────────────
print("Merging distance panel with trade flows...")

panel = baci.merge(ovdl_df, on=["iso_o","iso_d","year"], how="left")
if cepii_dist is not None:
    panel = panel.merge(cepii_dist, on=["iso_o","iso_d"], how="left")

# Require positive trade value and positive distance
panel = panel[panel["trade_value"] > 0].copy()
panel = panel[panel["d_ovdl_t"] > 0].copy()
print(f"Panel rows after filtering: {len(panel)}")
print(panel[["iso_o","iso_d","year","trade_value","d_ovdl_t"]].head(5))

# ── Step 5: Gravity estimation ─────────────────────────────────────────────────
print("\nRunning PPML gravity estimation...")

try:
    import pyfixest as pf

    results = []
    for dist_var in ["d_ovdl_t"] + (["d_cepii_distw"] if cepii_dist is not None else []):
        sub = panel.dropna(subset=[dist_var, "trade_value"]).copy()
        sub = sub[sub[dist_var] > 0].copy()

        if len(sub) < 50:
            print(f"  {dist_var}: too few observations ({len(sub)}), skipping")
            continue

        # If distance has no within-pair time variation (static fallback centroids),
        # drop pair FE to avoid perfect collinearity; still valid cross-section PPML.
        max_within_std = sub.groupby(["iso_o", "iso_d"])[dist_var].std().max()
        if max_within_std < 1e-6:
            formula = (f"trade_value ~ np.log({dist_var}) "
                       f"| iso_o^year + iso_d^year")
            spec_note = "cross-section (no pair FE — static centroids)"
        else:
            formula = (f"trade_value ~ np.log({dist_var}) "
                       f"| iso_o^year + iso_d^year + iso_o^iso_d")
            spec_note = "panel (3-way FE)"

        try:
            sub["pair_id"] = sub["iso_o"] + "_" + sub["iso_d"]
            fit = pf.fepois(fml=formula, data=sub,
                            vcov={"CRV1": "pair_id"})
            key = f"np.log({dist_var})"
            coefs = fit.coef()
            ses   = fit.se()
            b = float(coefs[key]) if key in coefs.index else float("nan")
            s = float(ses[key]) if key in ses.index else float("nan")
            n = int(fit._N)
            print(f"  {dist_var:20s}: β={b:.3f}  se={s:.3f}  N={n}  [{spec_note}]")
            results.append({"variant": dist_var, "beta": b, "se": s, "n_obs": n,
                            "spec": spec_note})
        except Exception as exc:
            print(f"  {dist_var}: estimation failed: {exc!r}")

    if results:
        res_df = pd.DataFrame(results)
        out_path = Path(__file__).parents[1] / "data" / "derived" / "gravity_benchmark_20c.csv"
        res_df.to_csv(out_path, index=False)
        print(f"\nSaved to {out_path}")
        print(res_df)

except ImportError:
    print("pyfixest not installed — install with: .venv/bin/pip install pyfixest")
