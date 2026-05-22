"""Validation suite for ot_distance.py.

Every claim that one of the OT measures is correct must be checkable
against either (a) a closed-form result on a synthetic distribution, or
(b) a known limit of the centroid measures.

This file is the discipline that lets us trust the numbers downstream.
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from paper5.ot_distance import (
    EARTH_RADIUS_KM,
    CES_DISTANCE_FLOOR_KM,
    haversine_km,
    ground_cost_matrix,
    RasterDist,
    ces_effective_distance,
    measure_M1_capital,
    measure_M2_unweighted_centroid,
    measure_M3_pop_weighted_centroid,
    measure_M5_directional_NL_centroid,
    measure_M6_wasserstein_1,
    measure_M7_wasserstein_iceberg,
    measure_M9_sliced_wasserstein,
    all_measures,
    HAS_POT,
)


# ── Haversine sanity ─────────────────────────────────────────────────────────

class TestHaversine:
    def test_zero_distance_at_same_point(self):
        d = haversine_km(np.array([45.0]), np.array([10.0]),
                         np.array([45.0]), np.array([10.0]))
        assert abs(d[0]) < 1e-9

    def test_equator_quarter_circle(self):
        """Equator to (lat=0, lon=90) — quarter great circle = π/2 * R."""
        d = haversine_km(np.array([0.0]), np.array([0.0]),
                         np.array([0.0]), np.array([90.0]))
        expected = math.pi / 2 * EARTH_RADIUS_KM
        assert abs(d[0] - expected) < 1.0  # within 1 km

    def test_pole_to_equator(self):
        """North pole to equator = π/2 * R."""
        d = haversine_km(np.array([90.0]), np.array([0.0]),
                         np.array([0.0]), np.array([0.0]))
        expected = math.pi / 2 * EARTH_RADIUS_KM
        assert abs(d[0] - expected) < 1.0

    def test_atlanta_to_london(self):
        """Atlanta (33.749, -84.388) → London (51.507, -0.128). Known ≈ 6,770 km
        on a spherical Earth (radius 6371 km, Haversine). The slightly-larger
        ~6,860 km figure assumes WGS-84 ellipsoid + Vincenty, which we don't use.
        """
        d = haversine_km(np.array([33.749]), np.array([-84.388]),
                         np.array([51.507]), np.array([-0.128]))
        assert 6720 < d[0] < 6820


# ── Centroid sanity ──────────────────────────────────────────────────────────

class TestCentroid:
    def test_centroid_of_single_atom_is_itself(self):
        d = RasterDist(coords=np.array([[20.0, 30.0]]), weights=np.array([1.0]))
        lat, lon = d.centroid()
        assert abs(lat - 20.0) < 1e-9
        assert abs(lon - 30.0) < 1e-9

    def test_centroid_of_symmetric_two_points(self):
        d = RasterDist(
            coords=np.array([[0.0, -10.0], [0.0, 10.0]]),
            weights=np.array([0.5, 0.5])
        )
        lat, lon = d.centroid()
        assert abs(lat) < 1e-6
        assert abs(lon) < 1e-6

    def test_centroid_handles_dateline(self):
        """Two points on opposite sides of dateline; spherical mean handles it."""
        d = RasterDist(
            coords=np.array([[0.0, 179.0], [0.0, -179.0]]),
            weights=np.array([0.5, 0.5])
        )
        lat, lon = d.centroid()
        # Should land near (0, 180) or (0, -180), not at (0, 0).
        assert abs(lat) < 1e-6
        assert abs(abs(lon) - 180.0) < 1.0


# ── M1, M2, M3, M5 measure tests ─────────────────────────────────────────────

class TestCentroidMeasures:
    def test_M1_known_pair(self):
        """Atlanta to London ≈ 6,770 km (spherical Earth, Haversine)."""
        d = measure_M1_capital((33.749, -84.388), (51.507, -0.128))
        assert 6720 < d < 6820

    def test_M3_reduces_to_M2_when_uniform_weights(self):
        coords = np.array([[10.0, 20.0], [11.0, 21.0], [12.0, 22.0]])
        rd_a = RasterDist(coords=coords, weights=np.array([1/3, 1/3, 1/3]))
        rd_b = RasterDist(coords=coords + 5.0, weights=np.array([1/3, 1/3, 1/3]))
        d_M3 = measure_M3_pop_weighted_centroid(rd_a, rd_b)
        d_M2 = measure_M2_unweighted_centroid(rd_a.coords, rd_b.coords)
        assert abs(d_M3 - d_M2) < 1e-6

    def test_M5_is_asymmetric(self):
        """For asymmetric distributions, M5(i,j) != M5(j,i)."""
        # Origin VIIRS concentrated to the east
        viirs_o = RasterDist(coords=np.array([[40.0, 10.0]]), weights=np.array([1.0]))
        # Destination LandScan concentrated to the west
        landscan_d = RasterDist(coords=np.array([[40.0, -10.0]]), weights=np.array([1.0]))
        # Reverse:
        viirs_d = RasterDist(coords=np.array([[40.0, -8.0]]), weights=np.array([1.0]))  # different VIIRS centroid
        landscan_o = RasterDist(coords=np.array([[40.0, 8.0]]), weights=np.array([1.0]))  # different LandScan centroid

        d_ij = measure_M5_directional_NL_centroid(viirs_o, landscan_d)
        d_ji = measure_M5_directional_NL_centroid(viirs_d, landscan_o)
        assert d_ij != d_ji


# ── OT measure tests (POT required) ──────────────────────────────────────────

@pytest.mark.skipif(not HAS_POT, reason="POT not installed")
class TestWassersteinClosedForms:

    def test_W1_equals_distance_for_point_masses(self):
        """Two single-atom distributions: W1 = great-circle distance between them."""
        mu = RasterDist(coords=np.array([[0.0, 0.0]]), weights=np.array([1.0]))
        nu = RasterDist(coords=np.array([[0.0, 90.0]]), weights=np.array([1.0]))
        w1 = measure_M6_wasserstein_1(mu, nu)
        expected = math.pi / 2 * EARTH_RADIUS_KM
        assert abs(w1 - expected) < 1.0

    def test_W1_zero_for_identical_distributions(self):
        coords = np.array([[10.0, 20.0], [11.0, 22.0], [12.0, 24.0]])
        weights = np.array([0.5, 0.3, 0.2])
        mu = RasterDist(coords=coords, weights=weights)
        # POT.emd2 on identical distributions can return tiny numerical noise
        w1 = measure_M6_wasserstein_1(mu, mu)
        assert abs(w1) < 1e-6

    def test_W1_dominates_centroid_for_bimodal(self):
        """The key theoretical motivation: when distributions are bimodal with
        the same centroid, W1 > 0 but centroid distance = 0. This is exactly
        the case where OT captures what centroids miss.
        """
        # mu = uniform on two points at lat=0, lon=±10. Centroid = (0,0).
        mu = RasterDist(
            coords=np.array([[0.0, -10.0], [0.0, 10.0]]),
            weights=np.array([0.5, 0.5])
        )
        # nu = single atom at (0, 0). Centroid = (0,0).
        nu = RasterDist(coords=np.array([[0.0, 0.0]]), weights=np.array([1.0]))
        # Centroid-to-centroid distance is zero
        d_cent = measure_M3_pop_weighted_centroid(mu, nu)
        assert d_cent < 1.0  # essentially zero
        # But W1 is positive — each unit of mass in mu must travel 10° ≈ 1113 km
        w1 = measure_M6_wasserstein_1(mu, nu)
        expected = haversine_km(np.array([0.0]), np.array([10.0]),
                                np.array([0.0]), np.array([0.0]))[0]
        # Average of two equal-mass transports of the same distance
        assert abs(w1 - expected) < 1.0

    def test_W1_iceberg_equals_W1_at_theta_one(self):
        """M7 with theta=1, p=1 must equal M6."""
        rng = np.random.default_rng(0)
        coords_a = rng.normal(loc=[40, 10], scale=2.0, size=(20, 2))
        coords_b = rng.normal(loc=[45, 15], scale=2.0, size=(20, 2))
        mu = RasterDist(coords=coords_a, weights=np.ones(20) / 20)
        nu = RasterDist(coords=coords_b, weights=np.ones(20) / 20)
        w6 = measure_M6_wasserstein_1(mu, nu)
        w7 = measure_M7_wasserstein_iceberg(mu, nu, p=1, theta=1.0)
        assert abs(w6 - w7) < 1e-6

    def test_sliced_W_correlates_with_exact_W1(self):
        """M9 sliced Wasserstein should give similar magnitude to M6 on
        comparable distributions (not exactly equal — it's an approximation).
        Tested with rank correlation on 20 country-like pair samples.
        """
        from scipy.stats import spearmanr
        rng = np.random.default_rng(0)
        sw_values = []
        w1_values = []
        for _ in range(20):
            ca = rng.normal(loc=[rng.uniform(20, 50), rng.uniform(-100, 100)], scale=3.0, size=(50, 2))
            cb = rng.normal(loc=[rng.uniform(20, 50), rng.uniform(-100, 100)], scale=3.0, size=(50, 2))
            mu = RasterDist(coords=ca, weights=np.ones(50) / 50)
            nu = RasterDist(coords=cb, weights=np.ones(50) / 50)
            sw_values.append(measure_M9_sliced_wasserstein(mu, nu, n_projections=200))
            w1_values.append(measure_M6_wasserstein_1(mu, nu))
        rho, _ = spearmanr(sw_values, w1_values)
        # On synthetic data the rank correlation should be high (>0.85);
        # if it drops below this the sliced approximation is broken.
        assert rho > 0.85, f"Sliced W rank-corr with exact W1 too low: {rho:.3f}"


# ── End-to-end horse-race smoke ──────────────────────────────────────────────

class TestCESEffectiveDistance:
    """Validation suite for the CES-weighted effective distance.

    This is the structurally-correct gravity-consistent measure. Per Head-Mayer
    (2002), it's a generalized mean of order θ = 1-σ over pairwise distances,
    weighted by activity densities at both ends.
    """

    def test_two_point_masses_equals_distance_when_theta_nonzero(self):
        """For singleton distributions, d_eff = d(point_i, point_j) for any θ."""
        mu = RasterDist(coords=np.array([[0.0, 0.0]]), weights=np.array([1.0]))
        nu = RasterDist(coords=np.array([[0.0, 90.0]]), weights=np.array([1.0]))
        expected = math.pi / 2 * EARTH_RADIUS_KM  # quarter great-circle
        for theta in (-7, -5, -3, -1, 1):
            d_eff = ces_effective_distance(mu, nu, theta=theta)
            assert abs(d_eff - expected) < 5.0, f"theta={theta}: d_eff={d_eff}, expected~{expected}"

    def test_theta_1_recovers_arithmetic_mean(self):
        """At θ=1, CES generalized mean is the arithmetic expected pairwise distance."""
        rng = np.random.default_rng(0)
        coords_a = rng.normal(loc=[40, 10], scale=2.0, size=(15, 2))
        coords_b = rng.normal(loc=[45, 15], scale=2.0, size=(15, 2))
        wa = np.ones(15) / 15; wb = np.ones(15) / 15
        mu = RasterDist(coords=coords_a, weights=wa)
        nu = RasterDist(coords=coords_b, weights=wb)
        d_eff = ces_effective_distance(mu, nu, theta=1.0)
        # Compute arithmetic mean directly
        D = ground_cost_matrix(coords_a, coords_b, kind="gc_km")
        arith_mean = float(np.sum(np.outer(wa, wb) * D))
        assert abs(d_eff - arith_mean) < 0.5, f"d_eff={d_eff}, arith_mean={arith_mean}"

    def test_theta_minus_one_recovers_harmonic_mean(self):
        """At θ=-1, CES generalized mean is the weighted harmonic mean."""
        coords_a = np.array([[0.0, 0.0], [0.0, 2.0]])
        coords_b = np.array([[0.0, 1.0], [1.0, 3.0]])
        wa = np.array([0.25, 0.75])
        wb = np.array([0.6, 0.4])
        mu = RasterDist(coords=coords_a, weights=wa)
        nu = RasterDist(coords=coords_b, weights=wb)

        D = np.maximum(ground_cost_matrix(coords_a, coords_b, kind="gc_km"), CES_DISTANCE_FLOOR_KM)
        W = np.outer(wa, wb)
        harmonic = float(1.0 / np.sum(W / D))
        assert ces_effective_distance(mu, nu, theta=-1.0) == pytest.approx(harmonic, rel=1e-12)

    def test_theta_zero_limit_recovers_geometric_mean(self):
        """At θ=0 and near θ=0, CES returns the weighted geometric mean."""
        coords_a = np.array([[0.0, 0.0], [0.0, 2.0], [2.0, 0.0]])
        coords_b = np.array([[0.0, 1.0], [3.0, 4.0]])
        wa = np.array([0.2, 0.3, 0.5])
        wb = np.array([0.7, 0.3])
        mu = RasterDist(coords=coords_a, weights=wa)
        nu = RasterDist(coords=coords_b, weights=wb)

        D = np.maximum(ground_cost_matrix(coords_a, coords_b, kind="gc_km"), CES_DISTANCE_FLOOR_KM)
        W = np.outer(wa, wb)
        geometric = float(np.exp(np.sum(W * np.log(D))))
        for theta in (-1e-7, 0.0, 1e-7):
            assert ces_effective_distance(mu, nu, theta=theta) == pytest.approx(geometric, rel=1e-6)

    def test_negative_theta_logsumexp_stable_empirical_range(self):
        """Negative θ stays stable for θ ∈ [-7, -3] and distances from 0.5 to ~20,000 km."""
        km_to_lon_deg = lambda km: math.degrees(km / EARTH_RADIUS_KM)
        coords_a = np.array([[0.0, 0.0]])
        coords_b = np.array([
            [0.0, km_to_lon_deg(0.5)],
            [0.0, km_to_lon_deg(50.0)],
            [0.0, 90.0],
            [0.0, km_to_lon_deg(20000.0)],
        ])
        wa = np.array([1.0])
        wb = np.array([0.02, 0.18, 0.30, 0.50])
        mu = RasterDist(coords=coords_a, weights=wa)
        nu = RasterDist(coords=coords_b, weights=wb)
        D = np.maximum(ground_cost_matrix(coords_a, coords_b, kind="gc_km"), CES_DISTANCE_FLOOR_KM)
        W = np.outer(wa, wb)

        assert D.min() == pytest.approx(CES_DISTANCE_FLOOR_KM, rel=1e-6)
        assert D.max() > 19000
        for theta in (-7.0, -5.0, -3.0):
            log_terms = np.log(W[W > 0]) + theta * np.log(D[W > 0])
            m = np.max(log_terms)
            expected = float(np.exp((m + np.log(np.sum(np.exp(log_terms - m)))) / theta))
            observed = ces_effective_distance(mu, nu, theta=theta)
            assert np.isfinite(observed)
            assert observed == pytest.approx(expected, rel=1e-12)

    def test_distance_floor_does_not_change_values_above_floor(self):
        """The 0.5 km clip has no effect when all pairwise distances are already above it."""
        km_to_lon_deg = lambda km: math.degrees(km / EARTH_RADIUS_KM)
        coords_a = np.array([[0.0, 0.0], [0.0, km_to_lon_deg(25.0)]])
        coords_b = np.array([[0.0, km_to_lon_deg(5.0)], [0.0, km_to_lon_deg(100.0)]])
        wa = np.array([0.35, 0.65])
        wb = np.array([0.4, 0.6])
        mu = RasterDist(coords=coords_a, weights=wa)
        nu = RasterDist(coords=coords_b, weights=wb)
        D = ground_cost_matrix(coords_a, coords_b, kind="gc_km")
        W = np.outer(wa, wb)

        assert D.min() > CES_DISTANCE_FLOOR_KM
        for theta in (-7.0, -5.0, -3.0):
            expected = float(np.sum(W * D**theta) ** (1.0 / theta))
            assert ces_effective_distance(mu, nu, theta=theta) == pytest.approx(expected, rel=1e-12)

    def test_theta_negative_infinity_converges_to_min_distance(self):
        """As θ → -∞, the CES mean converges to min(d); check θ = -10, -50, -100."""
        km_to_lon_deg = lambda km: math.degrees(km / EARTH_RADIUS_KM)
        mu = RasterDist(coords=np.array([[0.0, 0.0]]), weights=np.array([1.0]))
        nu = RasterDist(
            coords=np.array([
                [0.0, km_to_lon_deg(10.0)],
                [0.0, km_to_lon_deg(20.0)],
                [0.0, km_to_lon_deg(100.0)],
            ]),
            weights=np.ones(3) / 3,
        )
        D = np.maximum(ground_cost_matrix(mu.coords, nu.coords, kind="gc_km"), CES_DISTANCE_FLOOR_KM)
        min_d = float(D.min())

        values = [ces_effective_distance(mu, nu, theta=theta) for theta in (-10.0, -50.0, -100.0)]
        assert values[0] > values[1] > values[2] > min_d
        assert values[-1] == pytest.approx(min_d, rel=0.02)

    def test_overlapping_support_is_regularized_by_distance_floor(self):
        """Exact zero raw distances are finite after flooring and converge to the floor."""
        coords_a = np.array([[0.0, 0.0], [0.0, 1.0]])
        coords_b = np.array([[0.0, 0.0], [0.0, 2.0]])
        wa = np.array([0.2, 0.8])
        wb = np.array([0.5, 0.5])
        mu = RasterDist(coords=coords_a, weights=wa)
        nu = RasterDist(coords=coords_b, weights=wb)
        raw = ground_cost_matrix(coords_a, coords_b, kind="gc_km")
        D = np.maximum(raw, CES_DISTANCE_FLOOR_KM)
        W = np.outer(wa, wb)

        assert raw.min() == pytest.approx(0.0)
        assert ces_effective_distance(mu, nu, theta=-1.0) == pytest.approx(1.0 / np.sum(W / D), rel=1e-12)
        values = [ces_effective_distance(mu, nu, theta=theta) for theta in (-10.0, -50.0, -100.0)]
        assert all(np.isfinite(v) for v in values)
        assert values[0] > values[1] > values[2] > CES_DISTANCE_FLOOR_KM
        assert values[-1] == pytest.approx(CES_DISTANCE_FLOOR_KM, rel=0.05)

    def test_planar_fallback_wraps_longitude_and_scales_by_latitude(self):
        """Small high-latitude antimeridian crossings should match haversine scale."""
        mu = RasterDist(coords=np.array([[60.0, 179.9]]), weights=np.array([1.0]))
        nu = RasterDist(coords=np.array([[60.0, -179.9]]), weights=np.array([1.0]))
        d_haversine = ces_effective_distance(mu, nu, theta=1.0, haversine_great_circle=True)
        d_planar = ces_effective_distance(mu, nu, theta=1.0, haversine_great_circle=False)

        assert d_haversine > CES_DISTANCE_FLOOR_KM
        assert d_planar == pytest.approx(d_haversine, rel=0.02)

    def test_concave_theta_dominated_by_short_distances(self):
        """For two distributions with both short and long pairs, θ << 0 should
        yield a value much closer to the short-pair distance than to the
        arithmetic mean."""
        # Origin: two atoms 50:50 — one close to dest, one far
        mu_coords = np.array([[40.0, 10.0], [40.0, 100.0]])  # one near, one 90° away
        mu = RasterDist(coords=mu_coords, weights=np.array([0.5, 0.5]))
        nu = RasterDist(coords=np.array([[40.0, 12.0]]), weights=np.array([1.0]))
        # Pair distances: ~170 km and ~7,500 km
        # Arithmetic mean ≈ 3,840; harmonic mean (θ=-1) ≈ 330; θ=-5 ≈ much closer to 170
        d_arith = ces_effective_distance(mu, nu, theta=1.0)
        d_harm = ces_effective_distance(mu, nu, theta=-1.0)
        d_concave5 = ces_effective_distance(mu, nu, theta=-5.0)
        d_concave7 = ces_effective_distance(mu, nu, theta=-7.0)
        # Concave should be monotonically closer to the short distance as θ → -∞
        assert d_concave7 < d_concave5 < d_harm < d_arith
        assert d_concave7 < 500, f"theta=-7 should be near the short distance (~170 km), got {d_concave7}"

    def test_head_mayer_closed_form_constants(self):
        """Verify the canonical 0.67 sqrt(area/π) formula gives ~0.376 sqrt(area)."""
        from paper5.ot_distance import head_mayer_closed_form
        # Test on US area (~9.8M km²): expected ~ 0.67 * sqrt(9.8e6 / π) ≈ 1180 km
        d_us = head_mayer_closed_form(9.8e6)
        assert 1100 < d_us < 1250
        # Cross-check the relationship: HM(area) ≈ 0.376 * sqrt(area)
        assert abs(d_us - 0.376 * math.sqrt(9.8e6)) < 10.0


@pytest.mark.skipif(not HAS_POT, reason="POT not installed")
def test_all_measures_runs_and_returns_floats():
    """Smoke test: all_measures() returns a dict of finite floats."""
    rng = np.random.default_rng(42)
    coords_o = rng.normal(loc=[40, -100], scale=2.0, size=(30, 2))
    coords_d = rng.normal(loc=[50, 10], scale=2.0, size=(30, 2))
    viirs_o = RasterDist(coords=coords_o, weights=np.ones(30) / 30)
    landscan_d = RasterDist(coords=coords_d, weights=np.ones(30) / 30)
    out = all_measures(viirs_o, landscan_d,
                       origin_capital=(38.9, -77.0), dest_capital=(51.5, -0.1))
    assert "M1_capital" in out
    assert "M3_pop_centroid" in out
    assert "M5_directional_NL" in out
    assert "M6_W1" in out
    assert "M9_sliced_W" in out
    for k, v in out.items():
        assert np.isfinite(v), f"{k} returned non-finite: {v}"
