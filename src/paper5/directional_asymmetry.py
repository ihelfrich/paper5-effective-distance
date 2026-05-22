"""Production-consumption directional asymmetry in effective distance.

The standard Head-Mayer (2002, 2014) effective distance uses one density (typically
population) symmetrically: f_i = g_j = pop. That makes d_ij = d_ji by construction.

The geospatial novelty here: use VIIRS nightlights as the *origin* (production
proxy) and LandScan / GHS-POP as the *destination* (consumption proxy). Because
production density and consumption density disagree spatially (factories on
coasts, consumers in interior; light-emitting industry vs residential population),
the resulting CES effective distance is genuinely asymmetric:

    d_ij^eff (i exports to j)  =  (∑∑ viirs_i(x) · pop_j(y) · d(x,y)^θ )^(1/θ)
    d_ji^eff (j exports to i)  =  (∑∑ viirs_j(y) · pop_i(x) · d(x,y)^θ )^(1/θ)

These are NOT equal whenever the centroid of nightlights differs from the
centroid of population for at least one of i, j. Per Henderson-Storeygard-Weil
(2012) and Bluhm-Krause (2018), they routinely do.

Empirical claim (paper 5 prediction P1, per literature triangulation 2026-05-20):
the cross-sectional asymmetry index

    ASI_ij  =  log( d_ij^eff / d_ji^eff )

predicts the trade-flow asymmetry log(X_ij / X_ji) conditional on standard
gravity controls (GDP_i, GDP_j, contig, language, RTA, year FE) and even
conditional on origin-year × destination-year fixed effects in a Yotov-style
gravity. This has not been done in the literature.

This module does not estimate the regression — that lives in notebook 15
(`directional_asymmetry_pilot.py`). It builds the d_ij and d_ji effective
distance matrices given VIIRS + population rasters for a set of countries.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from paper5.ot_distance import RasterDist, ces_effective_distance


# ── Asymmetric bilateral distance ─────────────────────────────────────────────

def asymmetric_d_ij(
    origin_production: RasterDist,
    dest_consumption: RasterDist,
    *,
    theta: float = -5.0,
) -> float:
    """Compute one directional effective distance.

    The convention: `origin_production` is the activity-density at the origin
    (VIIRS nightlights, industrial-activity proxy) and `dest_consumption` is
    the activity-density at the destination (population, consumption proxy).
    """
    return ces_effective_distance(origin_production, dest_consumption, theta=theta)


def asymmetric_bilateral_matrix(
    production: dict[str, RasterDist],
    consumption: dict[str, RasterDist],
    *,
    theta: float = -5.0,
    countries: Optional[Sequence[str]] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Build the full (n × n) asymmetric distance matrix.

    Parameters
    ----------
    production : dict[ISO3 → RasterDist]
        Origin-side activity densities (typically VIIRS nightlights, per country).
    consumption : dict[ISO3 → RasterDist]
        Destination-side activity densities (typically GHS-POP / WorldPop, per country).
        Must have the same ISO3 keys as `production`.
    theta : float
        CES generalized-mean order. Default -5 (σ=6 gravity-consistent).
    countries : optional iterable
        Restrict to these ISO3 codes. Default: intersection of production & consumption.
    verbose : bool
        Print per-pair timing.

    Returns
    -------
    DataFrame with columns: iso_o, iso_d, d_ij_eff, theta.
        d_ij_eff is i → j (origin uses production_i, destination uses consumption_j).
        For internal flows (i==j), both come from country i but from different rasters,
        which produces an interesting d_ii^eff that captures internal
        production-consumption mismatch.
    """
    if countries is None:
        countries = sorted(set(production.keys()) & set(consumption.keys()))
    else:
        countries = [c for c in countries if c in production and c in consumption]

    rows = []
    n_pairs = len(countries) ** 2
    done = 0
    t0 = time.time()
    for iso_o in countries:
        for iso_d in countries:
            d = asymmetric_d_ij(
                production[iso_o], consumption[iso_d], theta=theta
            )
            rows.append({
                "iso_o": iso_o,
                "iso_d": iso_d,
                "d_ij_eff": d,
                "theta": theta,
            })
            done += 1
            if verbose and done % 50 == 0:
                rate = done / max(1e-6, time.time() - t0)
                eta = (n_pairs - done) / rate
                print(f"  {done}/{n_pairs}  rate={rate:.1f}/s  ETA {eta:.0f}s")
    return pd.DataFrame(rows)


def asymmetry_index(panel: pd.DataFrame) -> pd.DataFrame:
    """Add the log-asymmetry index ASI_ij = log(d_ij / d_ji).

    Returns one row per unordered pair {i, j} with i<j, and the index value.
    """
    p = panel.copy()
    # Wide: pivot to get d_ij and d_ji side by side
    fwd = p.rename(columns={"d_ij_eff": "d_ij"})[["iso_o", "iso_d", "d_ij"]]
    rev = p.rename(columns={"d_ij_eff": "d_ji", "iso_o": "iso_d", "iso_d": "iso_o"})[
        ["iso_o", "iso_d", "d_ji"]
    ]
    merged = fwd.merge(rev, on=["iso_o", "iso_d"], how="inner")
    merged = merged[merged["iso_o"] < merged["iso_d"]].copy()  # one direction per pair
    merged["ASI"] = np.log(merged["d_ij"] / merged["d_ji"])
    merged["abs_ASI"] = merged["ASI"].abs()
    return merged.sort_values("abs_ASI", ascending=False).reset_index(drop=True)


# ── Cross-raster centroid divergence (a simpler diagnostic) ──────────────────

def centroid_divergence(
    production: RasterDist, consumption: RasterDist
) -> tuple[float, float, float]:
    """Return (lat_diff_deg, lon_diff_deg, great_circle_km).

    This is the simpler M5 diagnostic — how far apart are the population
    centroid and the nightlights centroid within one country? Useful as a
    cheap sanity check before computing the full asymmetric matrix.
    """
    from paper5.ot_distance import haversine_km

    lat_p, lon_p = production.centroid()
    lat_c, lon_c = consumption.centroid()
    gc_km = float(
        np.asarray(haversine_km(
            np.array([lat_p]), np.array([lon_p]),
            np.array([lat_c]), np.array([lon_c]),
        )).reshape(-1)[0]
    )
    return (lat_p - lat_c, lon_p - lon_c, gc_km)


# ── Convenience: build a panel of intra-national centroid divergences ────────

def build_centroid_divergence_panel(
    production: dict[str, RasterDist],
    consumption: dict[str, RasterDist],
) -> pd.DataFrame:
    """One row per country with the internal P-C centroid divergence."""
    rows = []
    for iso in sorted(set(production.keys()) & set(consumption.keys())):
        try:
            lat_d, lon_d, gc_km = centroid_divergence(production[iso], consumption[iso])
            rows.append({
                "iso3": iso,
                "lat_diff_deg": lat_d,
                "lon_diff_deg": lon_d,
                "gc_km": gc_km,
            })
        except Exception as e:  # noqa
            rows.append({"iso3": iso, "lat_diff_deg": np.nan,
                         "lon_diff_deg": np.nan, "gc_km": np.nan})
    return pd.DataFrame(rows).sort_values("gc_km", ascending=False).reset_index(drop=True)
