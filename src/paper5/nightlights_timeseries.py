"""Time-series raster pipeline for nightlights 1992-2024.

Sources the Li et al. (2020) harmonized DMSP/VIIRS series and produces
per-country-per-year RasterDist objects. This is the input to:

  - Prediction P2 (Wasserstein velocity)
  - Question 1 (Disdier-Head missing-globalization test): does d_ij^eff
    change over time as internal density redistributes?

The harmonized series uses two naming conventions in the cache:
  Harmonized_DN_NTL_{year}_calDMSP.tif     (1992-2013, DMSP-calibrated)
  Harmonized_DN_NTL_{year}_simVIIRS.tif    (2014-2024, VIIRS-simulated to DMSP scale)

The simVIIRS files are produced by Li et al. via simulation; they are on the
DMSP digital-number scale (0-63) so the time series is comparable. Use the
calDMSP version for 1992-2013 and the simVIIRS version for 2014+.

Caveats:
  - The DMSP top-coding at DN=63 is severe in dense metros. Bluhm-Krause
    (2022) propose a Pareto correction. Apply post-mask if needed.
  - The series has known intercalibration issues at the DMSP-VIIRS transition.
    Robustness should restrict to either 1992-2013 or 2014-2024 alone.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from paper5.ot_distance import RasterDist
from paper5.region_raster import mask_raster_to_country, mask_result_to_raster_dist


NL_DIR = Path("/Volumes/HELFRICH-GD/KatiaBlendedFinance/raster_cache/nightlights/")


def nightlights_path(year: int) -> Path:
    """Path to the harmonized nightlights raster for the given year.

    Uses calDMSP for ≤ 2013, simVIIRS for ≥ 2014. Returns the file path,
    raising FileNotFoundError if absent.
    """
    suffix = "calDMSP" if year <= 2013 else "simVIIRS"
    p = NL_DIR / f"Harmonized_DN_NTL_{year}_{suffix}.tif"
    if not p.exists():
        # Fallback: try the other suffix
        other = "simVIIRS" if suffix == "calDMSP" else "calDMSP"
        p_alt = NL_DIR / f"Harmonized_DN_NTL_{year}_{other}.tif"
        if p_alt.exists():
            return p_alt
        raise FileNotFoundError(
            f"No nightlights raster for {year} at {p} or {p_alt}"
        )
    return p


def country_nightlights_year(
    iso3: str,
    year: int,
    boundaries: "gpd.GeoDataFrame",
    *,
    coarsen_factor: int = 5,
    max_atoms: int = 10000,
) -> Optional[RasterDist]:
    """Build a RasterDist for one country in one year from nightlights.

    Returns None if the country is missing from boundaries or the raster.
    """
    if iso3 not in boundaries["ISO3"].values:
        return None
    geom = boundaries[boundaries["ISO3"] == iso3].iloc[0].geometry
    try:
        path = nightlights_path(year)
    except FileNotFoundError:
        return None
    try:
        mask = mask_raster_to_country(
            path, geom, coarsen_factor=coarsen_factor, max_atoms=max_atoms,
        )
    except Exception as e:  # noqa
        return None
    if mask.kept_n_cells == 0:
        return None
    return mask_result_to_raster_dist(mask)


def build_country_timeseries(
    iso3: str,
    years: list[int],
    boundaries: "gpd.GeoDataFrame",
    *,
    coarsen_factor: int = 5,
    max_atoms: int = 10000,
    verbose: bool = False,
) -> dict[int, RasterDist]:
    """Build the year → RasterDist mapping for one country."""
    out = {}
    for y in years:
        t0 = time.time()
        rd = country_nightlights_year(
            iso3, y, boundaries,
            coarsen_factor=coarsen_factor, max_atoms=max_atoms,
        )
        if rd is None:
            continue
        out[y] = rd
        if verbose:
            print(f"    {iso3} {y}: {rd.coords.shape[0]} atoms ({time.time()-t0:.1f}s)")
    return out


def build_centroid_panel(
    iso3_codes: list[str],
    years: list[int],
    boundaries: "gpd.GeoDataFrame",
    *,
    coarsen_factor: int = 5,
    max_atoms: int = 10000,
    verbose: bool = True,
) -> pd.DataFrame:
    """Long-format: one row per (iso3, year) with the activity centroid.

    Output columns: iso3, year, centroid_lat, centroid_lon, total_mass, n_atoms.
    """
    rows = []
    for iso in iso3_codes:
        if verbose:
            print(f"\n{iso} ...", end="", flush=True)
        rds = build_country_timeseries(
            iso, years, boundaries,
            coarsen_factor=coarsen_factor, max_atoms=max_atoms,
        )
        for y, rd in rds.items():
            lat, lon = rd.centroid()
            rows.append({
                "iso3": iso, "year": y,
                "centroid_lat": lat, "centroid_lon": lon,
                "total_mass": float(rd.weights.sum() * rd.total_mass),
                "n_atoms": int(rd.coords.shape[0]),
            })
        if verbose:
            print(f" ({len(rds)} years)")
    return pd.DataFrame(rows)
