"""test_distance.py — unit tests for Python-side distance functions.

These tests use only synthetic data (no external files). They run in CI
without any mounted drives or raster data.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

from paper5.distance import (
    great_circle_km,
    ces_aggregate,
    compute_d_ovdl,
    compute_bilateral_distance,
    ChokepointState,
    Agglomeration,
)


# ---------------------------------------------------------------------------
# great_circle_km
# ---------------------------------------------------------------------------

def test_gc_zero_same_point():
    assert great_circle_km(0.0, 0.0, 0.0, 0.0) == pytest.approx(0.0, abs=1e-6)


def test_gc_equatorial_degree():
    """1 degree of longitude at equator ≈ 111.32 km."""
    d = great_circle_km(0.0, 0.0, 1.0, 0.0)
    assert d == pytest.approx(111.32, abs=0.2)


def test_gc_known_pair():
    """London → New York ≈ 5,570 km (reference: multiple sources)."""
    d = great_circle_km(-0.1278, 51.5074, -74.0060, 40.7128)
    assert 5_500 < d < 5_650


def test_gc_symmetry():
    d1 = great_circle_km(10.0, 50.0, 100.0, 20.0)
    d2 = great_circle_km(100.0, 20.0, 10.0, 50.0)
    assert d1 == pytest.approx(d2, rel=1e-9)


# ---------------------------------------------------------------------------
# ces_aggregate
# ---------------------------------------------------------------------------

def test_ces_uniform_weights_and_costs():
    """With uniform weights and identical costs, CES = that cost regardless of θ."""
    costs = np.ones((3, 3)) * 100.0
    w_o = np.ones(3) / 3
    w_d = np.ones(3) / 3
    for theta in (-1.0, 1.0, 2.0):
        result = ces_aggregate(costs, w_o, w_d, theta=theta)
        assert result == pytest.approx(100.0, rel=1e-4), f"Failed for θ={theta}"


def test_ces_theta_minus1_harmonic():
    """θ=-1 is harmonic mean; verify with 2×2 known values."""
    # Costs: [[100, 200], [300, 400]], equal weights
    costs = np.array([[100.0, 200.0], [300.0, 400.0]])
    w = np.array([0.5, 0.5])
    # Expected: [0.25*100^-1 + 0.25*200^-1 + 0.25*300^-1 + 0.25*400^-1]^(-1)
    expected = 1.0 / (0.25 / 100 + 0.25 / 200 + 0.25 / 300 + 0.25 / 400)
    result = ces_aggregate(costs, w, w, theta=-1.0)
    assert result == pytest.approx(expected, rel=1e-5)


def test_ces_theta1_arithmetic():
    costs = np.array([[100.0, 200.0], [300.0, 400.0]])
    w = np.array([0.5, 0.5])
    expected = 0.25 * (100 + 200 + 300 + 400)
    result = ces_aggregate(costs, w, w, theta=1.0)
    assert result == pytest.approx(expected, rel=1e-5)


def test_ces_raises_theta_zero():
    costs = np.ones((2, 2))
    with pytest.raises(ValueError, match="θ=0"):
        ces_aggregate(costs, np.array([0.5, 0.5]), np.array([0.5, 0.5]), theta=0.0)


# ---------------------------------------------------------------------------
# ChokepointState
# ---------------------------------------------------------------------------

def test_chokepoint_default_multiplier_one():
    cp = ChokepointState()
    assert cp.multiplier(0) == pytest.approx(1.0)


def test_chokepoint_suez_closed():
    cp = ChokepointState.suez_ever_given_2021()
    from paper5.distance import ChokepointState as CS
    # Import the Rust flag constants via the Python module
    # In Python-side we check suez_open=False raises cost
    assert not cp.suez_open


def test_chokepoint_panama_multiplier():
    cp = ChokepointState.panama_drought_2023()
    assert cp.panama_capacity_factor == pytest.approx(1.8)


def test_chokepoint_for_year():
    assert ChokepointState.for_year(2015).panama_capacity_factor == 1.0
    assert ChokepointState.for_year(2023).panama_capacity_factor == pytest.approx(1.8)
    assert ChokepointState.for_year(2024).red_sea_risk == pytest.approx(2.5)


# ---------------------------------------------------------------------------
# compute_d_ovdl
# ---------------------------------------------------------------------------

def _make_agglom(iso3: str, lon: float, lat: float, pop: float, viirs: float = 1.0) -> Agglomeration:
    return Agglomeration(iso3=iso3, agglom_id=0, lon=lon, lat=lat,
                         pop=pop, viirs=viirs, year=2022)


def test_ovdl_single_agglom_is_great_circle():
    """With one agglomeration per country, d_ovdl = great-circle distance."""
    ao = [_make_agglom("AAA", 0.0, 0.0, 1.0)]
    ad = [_make_agglom("BBB", 10.0, 0.0, 1.0)]
    expected = great_circle_km(0.0, 0.0, 10.0, 0.0)
    result = compute_d_ovdl(ao, ad, theta=-1.0, use_viirs=False)
    assert result == pytest.approx(expected, rel=1e-4)


def test_ovdl_two_aggloms_weighted():
    """Two equal-pop aggloms give CES-weighted mean of pair-distances."""
    ao = [_make_agglom("A", 0.0, 0.0, 1.0), _make_agglom("A", 1.0, 0.0, 1.0)]
    ad = [_make_agglom("B", 10.0, 0.0, 1.0)]
    # Each pair has equal weight (1/2) × 1
    d00 = great_circle_km(0.0, 0.0, 10.0, 0.0)
    d10 = great_circle_km(1.0, 0.0, 10.0, 0.0)
    expected = (0.5 * d00 ** -1 + 0.5 * d10 ** -1) ** -1
    result = compute_d_ovdl(ao, ad, theta=-1.0)
    assert result == pytest.approx(expected, rel=1e-4)


def test_ovdl_viirs_vs_pop_differ_when_unequal():
    """Activity-weighted and pop-weighted should differ when VIIRS ≠ pop distribution."""
    # agglom 0: small pop, high VIIRS; agglom 1: large pop, low VIIRS
    ao = [
        _make_agglom("A", 0.0, 0.0, pop=1.0, viirs=10.0),
        _make_agglom("A", 5.0, 0.0, pop=9.0, viirs=1.0),
    ]
    ad = [_make_agglom("B", 100.0, 0.0, pop=1.0, viirs=1.0)]
    d_pop = compute_d_ovdl(ao, ad, theta=-1.0, use_viirs=False)
    d_viirs = compute_d_ovdl(ao, ad, theta=-1.0, use_viirs=True)
    # Agglom-0 is closer (lon=0) and has high VIIRS, so viirs-weighted < pop-weighted
    assert d_viirs != pytest.approx(d_pop, rel=0.01)


# ---------------------------------------------------------------------------
# _arc_chokepoint_flags — ocean basin regression tests
# ---------------------------------------------------------------------------

def test_flags_mumbai_hamburg_red_sea():
    """Mumbai → Hamburg arc must carry SUEZ+RED_SEA flags (ocean basin crossing)."""
    from paper5.distance import _arc_chokepoint_flags
    flags = _arc_chokepoint_flags(72.83, 18.93, 9.87, 53.58)
    assert flags & 4, "RED_SEA flag should be set for Mumbai→Hamburg"
    assert flags & 1, "SUEZ flag should be set for Mumbai→Hamburg"


def test_flags_novorossiysk_not_a_bypass():
    """Mumbai → Novorossiysk must also carry SUEZ+RED_SEA flags.

    Novorossiysk (Black Sea, lon=37.8°E) is in the Med/Black Sea basin;
    any route from the Indian Ocean to Novorossiysk must transit the Red Sea.
    This is the regression test for the 'Novorossiysk bypass' bug where
    Dijkstra could route Mumbai→Novorossiysk→Hamburg at no chokepoint cost.
    """
    from paper5.distance import _arc_chokepoint_flags
    flags = _arc_chokepoint_flags(72.83, 18.93, 37.77, 44.72)  # Mumbai→Novorossiysk
    assert flags & 4, "RED_SEA flag must be set for Mumbai→Novorossiysk (Black Sea basin)"
    assert flags & 1, "SUEZ flag must be set for Mumbai→Novorossiysk"


def test_flags_panama_la_rotterdam():
    """LA → Rotterdam arc must carry PANAMA flag (Pacific ↔ Atlantic)."""
    from paper5.distance import _arc_chokepoint_flags
    flags = _arc_chokepoint_flags(-118.26, 33.74, 4.48, 51.92)
    assert flags & 2, "PANAMA flag should be set for LA→Rotterdam"


def test_flags_cape_route_unflagged():
    """Mumbai → Cape Town and Cape Town → Hamburg must have no Suez/Red Sea flags.

    These arcs form the Cape of Good Hope bypass route. If they were flagged,
    Dijkstra could not find a valid reroute under Suez disruption.
    """
    from paper5.distance import _arc_chokepoint_flags
    # Mumbai → Cape Town (southern Africa, neither east nor west of Suez basin)
    f1 = _arc_chokepoint_flags(72.83, 18.93, 18.42, -33.92)
    assert not (f1 & 4), "Cape Town should NOT have Red Sea flag (Cape bypass route)"
    # Cape Town → Hamburg
    f2 = _arc_chokepoint_flags(18.42, -33.92, 9.87, 53.58)
    assert not (f2 & 4), "Cape Town→Hamburg should NOT have Red Sea flag"


def test_flags_sydney_la_no_red_sea():
    """Sydney → LA is a trans-Pacific route; no Suez/Red Sea flags."""
    from paper5.distance import _arc_chokepoint_flags
    flags = _arc_chokepoint_flags(151.21, -33.87, -118.26, 33.74)
    assert not (flags & 4), "Sydney→LA (trans-Pacific) should not carry Red Sea flag"
    assert not (flags & 1), "Sydney→LA should not carry Suez flag"
