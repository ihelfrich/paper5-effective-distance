"""Diagnostic: how does max_atoms top-N sampling bias d_eff(θ<0)?

Pick one small country (Netherlands), compute d_ii^eff at θ=-5 under several
sampling strategies:

  (a) top-N by population (N = 500, 1000, 2000, 5000, all ~40k cells)
  (b) stratified random (N = 500, 1000, 2000, 5000) — random subsample
      proportional to weight, preserving spatial coverage
  (c) full resolution (~40k atoms, ~1.6B pairs — only feasible for small
      countries)

If (a) at N=5000 differs from (c) by more than 5%, we have a real sampling
problem and need to switch to (b) or to a block-coarsening strategy for
the panel.

Run:
  .venv/bin/python -u notebooks/16_sampling_bias_diagnostic.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

from paper5.ot_distance import RasterDist, ces_effective_distance, head_mayer_closed_form
from paper5.region_raster import (
    load_country_boundaries,
    mask_raster_to_country,
)


# Pick small countries where full-resolution is tractable
SMALL_COUNTRIES = [
    ("NLD", "Netherlands"),
    ("BEL", "Belgium"),
    ("CHE", "Switzerland"),
    ("DNK", "Denmark"),
]

GHS_POP = Path("/Volumes/HELFRICH-GD/TEG_data/inputs/ghsl/pop/"
               "GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif")

SAMPLE_SIZES = [500, 1000, 2000, 5000, 10000, None]  # None = all
RNG = np.random.default_rng(42)


def stratified_subsample(coords, weights, n_target):
    """Random subsample proportional to weight. Preserves spatial distribution.

    Returns (coords_sub, weights_sub) where weights_sub renormalizes to 1.
    """
    if n_target is None or n_target >= len(weights):
        return coords, weights / weights.sum()
    idx = RNG.choice(len(weights), size=n_target, replace=False,
                     p=weights / weights.sum())
    return coords[idx], weights[idx] / weights[idx].sum()


def top_n_subsample(coords, weights, n_target):
    """Keep top-N cells by raster value. Concentrates in dense areas."""
    if n_target is None or n_target >= len(weights):
        return coords, weights / weights.sum()
    top = np.argsort(weights)[-n_target:]
    return coords[top], weights[top] / weights[top].sum()


def main():
    boundaries = load_country_boundaries()
    print(f"Loaded {len(boundaries)} country boundaries\n")

    results = []
    for iso, name in SMALL_COUNTRIES:
        if iso not in boundaries["ISO3"].values:
            print(f"  {iso}: MISSING in boundaries"); continue

        geom = boundaries[boundaries["ISO3"] == iso].iloc[0].geometry
        print(f"\n=== {iso} {name} ===")

        # Get the FULL mask (no max_atoms)
        t0 = time.time()
        m = mask_raster_to_country(GHS_POP, geom, max_atoms=None)
        print(f"  Full mask: {m.kept_n_cells:,} cells, "
              f"area = {m.full_area_km2:.0f} km², "
              f"loaded in {time.time()-t0:.1f}s")

        if m.kept_n_cells == 0:
            print(f"  EMPTY — skip")
            continue

        d_HM = head_mayer_closed_form(m.full_area_km2)
        print(f"  d_HM (closed-form): {d_HM:.0f} km")

        coords_full = np.column_stack([m.cell_lats, m.cell_lons])
        weights_full = m.raster_values * m.cell_areas_km2

        # Full-resolution baseline (capped at 12000 to keep pairwise feasible)
        N_FULL = min(12000, m.kept_n_cells)
        # Even for the "baseline" we have to subsample if N is large.
        # Use stratified random as the closest to truth.
        coords_base, w_base = stratified_subsample(coords_full, weights_full, N_FULL)
        rd_base = RasterDist(coords=coords_base, weights=w_base)
        t0 = time.time()
        d_baseline_neg5 = ces_effective_distance(rd_base, rd_base, theta=-5.0)
        t_base = time.time() - t0
        print(f"  Stratified-{N_FULL} baseline:  d_eff(θ=-5) = {d_baseline_neg5:.1f} km  "
              f"({t_base:.1f}s)")

        for n in SAMPLE_SIZES:
            if n is None or n >= m.kept_n_cells:
                continue
            # Top-N sampling
            c_top, w_top = top_n_subsample(coords_full, weights_full, n)
            rd_top = RasterDist(coords=c_top, weights=w_top)
            d_top = ces_effective_distance(rd_top, rd_top, theta=-5.0)

            # Stratified
            c_str, w_str = stratified_subsample(coords_full, weights_full, n)
            rd_str = RasterDist(coords=c_str, weights=w_str)
            d_str = ces_effective_distance(rd_str, rd_str, theta=-5.0)

            err_top = (d_top - d_baseline_neg5) / d_baseline_neg5
            err_str = (d_str - d_baseline_neg5) / d_baseline_neg5

            results.append({
                "iso": iso, "name": name, "n_atoms": n,
                "d_HM": d_HM,
                "d_top_neg5": d_top, "d_strat_neg5": d_str,
                "d_baseline_neg5": d_baseline_neg5,
                "rel_err_top": err_top, "rel_err_strat": err_str,
            })

            print(f"  N = {n:>5}: top-N = {d_top:>6.1f} km ({err_top:+.1%}),  "
                  f"strat = {d_str:>6.1f} km ({err_str:+.1%})")

    df = pd.DataFrame(results)
    out = _repo / "data" / "derived" / "sampling_bias_diagnostic.csv"
    df.to_csv(out, index=False)
    print(f"\n=== Summary ===")
    if len(df) > 0:
        for col in ["rel_err_top", "rel_err_strat"]:
            print(f"  {col:>18}: median = {df[col].median():+.2%}  "
                  f"abs median = {df[col].abs().median():.2%}  "
                  f"max abs = {df[col].abs().max():.2%}")
    print(f"\nSaved {len(df)} rows to {out}")


if __name__ == "__main__":
    main()
