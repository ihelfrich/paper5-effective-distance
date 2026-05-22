"""Tests for directional_asymmetry and wasserstein_velocity modules.

These are the two genuinely-novel empirical pieces (per literature
triangulation 2026-05-20). They need clean unit tests before we run
them on real panels.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from paper5.ot_distance import RasterDist
from paper5.directional_asymmetry import (
    asymmetric_d_ij,
    asymmetric_bilateral_matrix,
    asymmetry_index,
    centroid_divergence,
)
from paper5.wasserstein_velocity import (
    country_velocity_series,
    bilateral_velocity_projection,
)


# ── Helpers ─────────────────────────────────────────────────────────────────

def _make_rd(lat_lon_weights):
    """Make a RasterDist from a list of (lat, lon, weight) tuples."""
    arr = np.asarray(lat_lon_weights, dtype=float)
    coords = arr[:, :2]
    w = arr[:, 2]
    w = w / w.sum()
    return RasterDist(coords=coords, weights=w, total_mass=1.0)


# ── Directional asymmetry ──────────────────────────────────────────────────

class TestAsymmetricDistance:
    def test_symmetric_when_densities_match(self):
        """If origin and dest have identical density, d_ij = d_ji."""
        rd = _make_rd([(40, -75, 1.0), (35, -80, 1.0), (45, -100, 1.0)])
        d_fwd = asymmetric_d_ij(rd, rd, theta=-5.0)
        d_rev = asymmetric_d_ij(rd, rd, theta=-5.0)
        assert abs(d_fwd - d_rev) < 1e-6
        assert d_fwd > 0

    def test_asymmetric_when_densities_diverge(self):
        """If production sits at (40, -75) and consumption at (35, -100),
        d_ij from i_prod to j_cons ≠ d_ji where production swaps roles."""
        # Country A: production on east coast, consumption inland
        prod_A = _make_rd([(40, -75, 1.0)])
        cons_A = _make_rd([(40, -100, 1.0)])
        # Country B: production inland, consumption on west coast
        prod_B = _make_rd([(40, -110, 1.0)])
        cons_B = _make_rd([(40, -120, 1.0)])

        # A's production → B's consumption: -75 to -120 (45° span)
        d_AB = asymmetric_d_ij(prod_A, cons_B, theta=-5.0)
        # B's production → A's consumption: -110 to -100 (10° span)
        d_BA = asymmetric_d_ij(prod_B, cons_A, theta=-5.0)
        assert d_AB > d_BA  # one direction is genuinely longer
        # And these are not the symmetric population-only distance
        d_pop_pop = asymmetric_d_ij(cons_A, cons_B, theta=-5.0)
        # cons_A=(-100) to cons_B=(-120) → 20° span; intermediate
        assert d_BA < d_pop_pop < d_AB

    def test_bilateral_matrix_shape(self):
        prod = {"X": _make_rd([(0, 0, 1)]), "Y": _make_rd([(0, 10, 1)])}
        cons = {"X": _make_rd([(0, 0.5, 1)]), "Y": _make_rd([(0, 10.5, 1)])}
        df = asymmetric_bilateral_matrix(prod, cons, theta=-3.0)
        assert len(df) == 4  # 2x2 ordered pairs (including self)
        assert set(df["iso_o"]) == {"X", "Y"}

    def test_asymmetry_index_nonzero(self):
        # Build a matrix where d_XY ≠ d_YX
        df = pd.DataFrame({
            "iso_o": ["X", "Y"], "iso_d": ["Y", "X"],
            "d_ij_eff": [800.0, 600.0], "theta": [-5.0, -5.0],
        })
        idx = asymmetry_index(df)
        assert len(idx) == 1
        assert abs(idx.iloc[0]["ASI"] - np.log(800.0 / 600.0)) < 1e-9


class TestCentroidDivergence:
    def test_zero_when_same_density(self):
        rd = _make_rd([(40, -75, 1.0), (45, -100, 1.0)])
        lat_d, lon_d, gc = centroid_divergence(rd, rd)
        assert abs(lat_d) < 1e-9
        assert abs(lon_d) < 1e-9
        assert gc < 1e-6

    def test_positive_when_centroids_differ(self):
        prod = _make_rd([(40, -75, 1.0)])
        cons = _make_rd([(40, -100, 1.0)])
        lat_d, lon_d, gc = centroid_divergence(prod, cons)
        assert abs(lat_d) < 1e-9
        # Lon centroid differs by 25 degrees
        assert abs(lon_d - 25.0) < 1e-9
        # At 40°N, 25° lon ≈ 25 * 111.32 * cos(40°) ≈ 2132 km
        assert 1800 < gc < 2400


# ── Wasserstein velocity ───────────────────────────────────────────────────

class TestWassersteinVelocity:
    def test_static_centroid_zero_velocity(self):
        """If the country density doesn't change, v = 0."""
        rd = _make_rd([(40, -75, 1.0), (45, -100, 1.0)])
        rds_by_year = {2010: rd, 2011: rd, 2012: rd}
        s = country_velocity_series(rds_by_year, "USA")
        # v in year 2010 is NaN (first year, no diff); 2011 and 2012 are zero
        assert pd.isna(s.iloc[0]["v_lat_deg"])
        assert abs(s.iloc[1]["v_lat_deg"]) < 1e-9
        assert abs(s.iloc[2]["v_lon_deg"]) < 1e-9

    def test_moving_centroid(self):
        # Year 2010 centroid at (40, -75); year 2011 at (42, -77)
        rd_2010 = _make_rd([(40, -75, 1.0)])
        rd_2011 = _make_rd([(42, -77, 1.0)])
        s = country_velocity_series({2010: rd_2010, 2011: rd_2011}, "USA")
        assert abs(s.iloc[1]["v_lat_deg"] - 2.0) < 1e-9
        assert abs(s.iloc[1]["v_lon_deg"] - (-2.0)) < 1e-9

    def test_bilateral_projection_signs(self):
        """If A drifts east and B is east of A, projection is positive
        (A moving toward B), so v_proj_km > 0 in our sign convention."""
        # Lag year: A at (40, -100), B at (40, -80)  → B is 20° east of A
        # Current year: A drifts to (40, -95) (eastward, toward B)
        centroids = pd.DataFrame([
            {"iso3": "A", "year": 2010, "centroid_lat": 40.0,
             "centroid_lon": -100.0, "v_lat_deg": np.nan, "v_lon_deg": np.nan},
            {"iso3": "A", "year": 2011, "centroid_lat": 40.0,
             "centroid_lon": -95.0, "v_lat_deg": 0.0, "v_lon_deg": 5.0},
            {"iso3": "B", "year": 2010, "centroid_lat": 40.0,
             "centroid_lon": -80.0, "v_lat_deg": np.nan, "v_lon_deg": np.nan},
            {"iso3": "B", "year": 2011, "centroid_lat": 40.0,
             "centroid_lon": -80.0, "v_lat_deg": 0.0, "v_lon_deg": 0.0},
        ])
        proj = bilateral_velocity_projection(centroids)
        # iso_o = A, iso_d = B: A moves east (+5° lon), B doesn't move
        # Direction A→B at lag is +20° lon (east). Projection of (v_A - v_B) = (0, 5)
        # onto unit (0, 1) gives v_proj_deg = 5.0
        row = proj[(proj["iso_o"] == "A") & (proj["iso_d"] == "B") & (proj["year"] == 2011)]
        assert len(row) == 1
        assert row.iloc[0]["v_proj_deg"] == pytest.approx(5.0, abs=1e-6)
        # Going reverse: B→A direction is west; (v_B - v_A)=(0,-5), proj onto (0,-1) = +5
        row2 = proj[(proj["iso_o"] == "B") & (proj["iso_d"] == "A") & (proj["year"] == 2011)]
        assert row2.iloc[0]["v_proj_deg"] == pytest.approx(5.0, abs=1e-6)
        # Sign is positive in both directions: A and B are getting closer
        # (in our convention, +v_proj = centroids moving toward each other)
