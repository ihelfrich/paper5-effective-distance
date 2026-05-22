"""Get a clean intra-national d_ii^eff(θ) using block-coarsened GHS-POP.

The killer experiment (nb18 v4) showed that the home-bias resolution lives
in intra-national distance correction, not bilateral contiguous-pair
correction. To quantify the resolution, we need an UNBIASED intra-national
d_ii^eff at θ=-5.

Block coarsening preserves total mass and spatial coverage while reducing
the atom count, avoiding the top-N concentration bias. We use 10 km blocks
on the 1 km GHS-POP Mollweide raster (100× atom reduction). For most
countries this gives ~5,000-50,000 atoms — feasible memory.

For each country, compute:
  - geographic area (from Natural Earth Mollweide polygon)
  - d_HM = 0.67 * sqrt(area/π) (CEPII intra-national closed form)
  - d_ii^eff at θ ∈ {1, -1, -3, -5, -7}

Compare. Save to /Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/.

Run:
  .venv/bin/python -u notebooks/19_intra_distance_block_coarsen.py
"""

from __future__ import annotations

import sys
import time
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

GHS_POP = Path("/Volumes/HELFRICH-GD/TEG_data/inputs/ghsl/pop/"
               "GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif")
OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived/")

TEST_COUNTRIES = [
    # Smaller subset to keep total runtime under ~30 min
    "USA", "RUS", "CAN", "BRA", "AUS", "CHN", "IND", "JPN", "FRA", "DEU",
    "NLD", "BEL", "KOR", "GBR", "ITA", "ESP", "MEX", "ARG", "CHL", "EGY",
]
THETAS = [1.0, -1.0, -3.0, -5.0, -7.0]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading boundaries...")
    bnd = load_country_boundaries()
    bnd_eq = bnd.to_crs("ESRI:54009")
    area_lookup = dict(zip(bnd_eq["ISO3"], (bnd_eq.geometry.area / 1e6).astype(float).values))
    print(f"  {len(bnd)} countries")

    rows = []
    t_start = time.time()
    for iso in TEST_COUNTRIES:
        if iso not in bnd["ISO3"].values:
            print(f"  {iso}: MISSING")
            continue
        geom = bnd[bnd["ISO3"] == iso].iloc[0].geometry
        area_geo = area_lookup.get(iso, np.nan)
        d_HM = head_mayer_closed_form(area_geo) if np.isfinite(area_geo) else np.nan

        t0 = time.time()
        try:
            # Block-coarsen at 10 km, then stratified-sample to keep under 12,000 atoms
            mask = mask_raster_to_country(
                GHS_POP, geom,
                coarsen_factor=10,           # 1km -> 10km blocks (100x atom reduction)
                max_atoms=12000,
                sample_strategy="stratified", # weight-proportional sampling preserves spatial coverage
            )
        except Exception as e:
            print(f"  {iso}: MASK FAILED: {e}")
            continue

        if mask.kept_n_cells == 0:
            print(f"  {iso}: empty mask")
            continue

        rd = mask_result_to_raster_dist(mask)
        t_mask = time.time() - t0

        d_eff = {}
        t0 = time.time()
        for theta in THETAS:
            d_eff[theta] = ces_effective_distance(rd, rd, theta=theta)
        t_ces = time.time() - t0

        rows.append({
            "iso3": iso,
            "area_km2": area_geo,
            "n_atoms": mask.kept_n_cells,
            "full_area_km2": mask.full_area_km2,
            "d_HM": d_HM,
            **{f"d_theta_{t}": d_eff[t] for t in THETAS},
            "ratio_neg5_to_HM": d_eff[-5.0] / d_HM if d_HM > 0 else np.nan,
            "ratio_neg3_to_HM": d_eff[-3.0] / d_HM if d_HM > 0 else np.nan,
            "t_mask_s": t_mask,
            "t_ces_s": t_ces,
        })
        print(f"  {iso:>3}  area={area_geo/1e6:>5.2f}Mkm² atoms={mask.kept_n_cells:>5}  "
              f"d_HM={d_HM:>6.0f}  d_-5={d_eff[-5.0]:>7.1f}  "
              f"ratio(-5/HM)={d_eff[-5.0]/d_HM if d_HM>0 else float('nan'):>6.3f}  "
              f"({t_mask:.1f}+{t_ces:.1f}s)")

    df = pd.DataFrame(rows)
    out = OUT_DIR / "intra_distance_block_coarsened.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved {len(df)} rows to {out}")
    print(f"Total wall-clock: {(time.time()-t_start)/60:.1f} min")

    print("\n=== Summary: d_ii^eff / d_HM ratios (block-coarsened, stratified) ===")
    for col, theta in [(c, t) for c, t in zip(
            [f"d_theta_{t}" for t in THETAS], THETAS)]:
        ratios = df[col] / df["d_HM"]
        print(f"  θ={theta:>5g}: median={ratios.median():.3f}  10th={ratios.quantile(0.1):.3f}  "
              f"90th={ratios.quantile(0.9):.3f}")


if __name__ == "__main__":
    main()
