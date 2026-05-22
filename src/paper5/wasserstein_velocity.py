"""Wasserstein-velocity diagnostic for activity-density shifts.

The novel prediction: countries whose activity-density centroid is *moving toward*
each other year-over-year should experience faster bilateral trade growth than
otherwise-equivalent country pairs whose centroids are moving away from each other.

The intuition is geometric. If China's economic mass shifts inland (away from
the coast) over time, then in the Head-Mayer (2002) sense the *effective*
distance from China to, say, Vietnam grows — even if their geographic
centroids haven't moved. Conversely, if Brazil's nightlight mass is migrating
toward its southern industrial belt and Argentina's mass sits in the Pampas
to the south, both shifts shrink Brazil-Argentina effective distance.

For each country i and year t, define the activity centroid

    c_i(t) = ∑_x w_i(x, t) · x  / ∑_x w_i(x, t)

where w_i(x, t) is VIIRS nightlights in year t (alternatively, GHS-POP / WorldPop
counts). The annual velocity is

    v_i(t) = c_i(t) - c_i(t-1)        (units: degrees lat/lon, or km after projection)

For bilateral pair (i, j), the *velocity projection* onto the i→j direction is

    π_ij(t) = (v_i(t) - v_j(t)) · ((c_j(t-1) - c_i(t-1)) / |c_j(t-1) - c_i(t-1)|)

A negative π_ij means i, j centroids are moving toward each other (effective
distance shrinking). Positive means they're drifting apart.

Empirical claim (paper 5 prediction P2, per literature triangulation 2026-05-20):
in a panel gravity with origin-year, destination-year, and pair fixed effects,
the coefficient on π_ij(t) is significantly negative for trade-flow growth
d log X_ij / dt, conditional on standard time-varying controls (GDP-pair growth,
sanction / RTA / contig × year fixed effects).

This is not preempted because the literature (Allen-Arkolakis 2014;
Donaldson-Hornbeck 2016) treats distance as time-invariant geography. The
Krugman-Venables agglomeration literature lets economic mass move, but does
not relate inter-country mass-shift coherence to bilateral trade growth.

Run on at least 10–15 years of harmonized DMSP/VIIRS nightlights to identify.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from paper5.ot_distance import RasterDist, haversine_km


@dataclass
class CountryVelocity:
    iso3: str
    year: int
    centroid_lat: float
    centroid_lon: float
    total_mass: float  # sum of weights for normalization


def country_velocity_series(
    rds_by_year: dict[int, RasterDist],
    iso3: str,
) -> pd.DataFrame:
    """Compute (year, centroid_lat, centroid_lon, total_mass) for one country."""
    rows = []
    for year in sorted(rds_by_year.keys()):
        rd = rds_by_year[year]
        if rd.coords.shape[0] == 0:
            rows.append({"iso3": iso3, "year": year,
                         "centroid_lat": np.nan, "centroid_lon": np.nan,
                         "total_mass": 0.0})
            continue
        lat, lon = rd.centroid()
        rows.append({"iso3": iso3, "year": year,
                     "centroid_lat": lat, "centroid_lon": lon,
                     "total_mass": float(rd.weights.sum())})
    df = pd.DataFrame(rows)
    df = df.sort_values("year").reset_index(drop=True)
    # First differences
    df["v_lat_deg"] = df["centroid_lat"].diff()
    df["v_lon_deg"] = df["centroid_lon"].diff()
    return df


def bilateral_velocity_projection(
    centroids: pd.DataFrame,
) -> pd.DataFrame:
    """For each (i, j, t), compute v_proj = (v_i - v_j) · û_ij(t-1).

    Input
    -----
    centroids : DataFrame with columns iso3, year, centroid_lat, centroid_lon,
                v_lat_deg, v_lon_deg (output of country_velocity_series).

    Output
    ------
    DataFrame: iso_o, iso_d, year, v_proj_deg, v_proj_km, gc_dist_km_lag1.
    """
    # Long → year-keyed lookup
    centroids = centroids.dropna(subset=["centroid_lat", "centroid_lon"]).copy()
    isos = centroids["iso3"].unique()
    rows = []

    for year in sorted(centroids["year"].unique()):
        yr = centroids[centroids["year"] == year].set_index("iso3")
        yr_lag = centroids[centroids["year"] == year - 1].set_index("iso3")
        # Need both this year (for velocity) and last year (for direction)
        common = sorted(set(yr.index) & set(yr_lag.index))
        for iso_o in common:
            for iso_d in common:
                if iso_o == iso_d:
                    continue
                v_lat_diff = yr.loc[iso_o, "v_lat_deg"] - yr.loc[iso_d, "v_lat_deg"]
                v_lon_diff = yr.loc[iso_o, "v_lon_deg"] - yr.loc[iso_d, "v_lon_deg"]
                if not np.isfinite(v_lat_diff) or not np.isfinite(v_lon_diff):
                    continue
                # Direction unit vector from i to j at t-1
                dlat = yr_lag.loc[iso_d, "centroid_lat"] - yr_lag.loc[iso_o, "centroid_lat"]
                dlon = yr_lag.loc[iso_d, "centroid_lon"] - yr_lag.loc[iso_o, "centroid_lon"]
                norm = math.hypot(dlat, dlon)
                if norm < 1e-9:
                    continue
                u_lat, u_lon = dlat / norm, dlon / norm
                v_proj_deg = v_lat_diff * u_lat + v_lon_diff * u_lon
                # In km: approximate at midpoint latitude
                mid_lat = 0.5 * (yr_lag.loc[iso_o, "centroid_lat"]
                                 + yr_lag.loc[iso_d, "centroid_lat"])
                v_proj_lat_km = v_lat_diff * 111.32
                v_proj_lon_km = v_lon_diff * 111.32 * math.cos(math.radians(mid_lat))
                v_proj_km = v_proj_lat_km * u_lat + v_proj_lon_km * u_lon

                gc_km = float(np.asarray(haversine_km(
                    np.array([yr_lag.loc[iso_o, "centroid_lat"]]),
                    np.array([yr_lag.loc[iso_o, "centroid_lon"]]),
                    np.array([yr_lag.loc[iso_d, "centroid_lat"]]),
                    np.array([yr_lag.loc[iso_d, "centroid_lon"]]),
                )).reshape(-1)[0])

                rows.append({
                    "iso_o": iso_o, "iso_d": iso_d, "year": year,
                    "v_proj_deg": v_proj_deg, "v_proj_km": v_proj_km,
                    "gc_dist_km_lag1": gc_km,
                })

    return pd.DataFrame(rows)


def attach_to_gravity_panel(
    gravity_panel: pd.DataFrame,
    velocity_panel: pd.DataFrame,
) -> pd.DataFrame:
    """Merge v_proj_km into a gravity panel on (iso_o, iso_d, year)."""
    return gravity_panel.merge(
        velocity_panel[["iso_o", "iso_d", "year", "v_proj_km"]],
        on=["iso_o", "iso_d", "year"],
        how="left",
    )
