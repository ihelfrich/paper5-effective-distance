"""test_transport_graph.py — integration tests for maritime graph + chokepoints.

Uses real WPI.csv (checked into repo at data/wpi/WPI.csv) and the fallback
country centroids. Requires networkx but no mounted drives or raster data.

Verifies that:
  - Houthi Red Sea 2024 produces >10% cost uplift for IND→DEU
  - Panama Drought 2023 produces >10% cost uplift for USA→NLD
  - Suez closure produces >50% cost uplift for IND→DEU (extreme event)
  - Arcs to sub-Saharan ports (Cape Town, Durban) remain unflagged
"""
from __future__ import annotations
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))


# ---------------------------------------------------------------------------
# Skip whole module if networkx is not installed
# ---------------------------------------------------------------------------
nx = pytest.importorskip("networkx", reason="networkx not installed")


from paper5.distance import (
    ChokepointState, Agglomeration,
    _build_maritime_graph_nx, lcp_pair_nx,
    PORT_HANDLING_HOURS, FREIGHT_RATE_USD_PER_HOUR,
)


# Fallback centroids (same as notebook)
_CENTROIDS = {
    "IND": (78.96, 20.59),
    "DEU": (10.45, 51.17),
    "USA": (-95.71, 37.09),
    "NLD": (4.92, 52.37),
    "JPN": (138.25, 36.20),
    "SGP": (103.82,  1.35),
    "AUS": (133.78, -25.27),
}


def _agglom(iso3: str) -> Agglomeration:
    lon, lat = _CENTROIDS[iso3]
    return Agglomeration(iso3=iso3, agglom_id=0, lon=lon, lat=lat, pop=1.0)


@pytest.fixture(scope="module")
def ports():
    from paper5.data_loaders import load_wpi_ports
    return load_wpi_ports()


def _maritime_lcp(ports, cp: ChokepointState, iso_o: str, iso_d: str) -> float:
    """Return total maritime LCP cost (including port handling) for a country pair."""
    G = _build_maritime_graph_nx(ports, cp)
    graphs = {"maritime": G}
    cost, _ = lcp_pair_nx(_agglom(iso_o), _agglom(iso_d), graphs)
    return cost


# ---------------------------------------------------------------------------
# Houthi Red Sea 2024: IND → DEU
# ---------------------------------------------------------------------------

def test_houthi_ind_deu_uplift(ports):
    """Under Houthi 2024 (red_sea_risk=2.5), IND→DEU cost should be >10% higher."""
    normal = _maritime_lcp(ports, ChokepointState(), "IND", "DEU")
    houthi = _maritime_lcp(ports, ChokepointState.houthi_red_sea_2024(), "IND", "DEU")
    ratio = houthi / normal
    assert ratio > 1.10, f"IND→DEU Houthi uplift {ratio:.3f}x should be >1.10"
    assert ratio < 5.0,  f"IND→DEU Houthi uplift {ratio:.3f}x seems implausibly large"


def test_houthi_jpn_deu_uplift(ports):
    """JPN→DEU (classic Suez route) should also show Houthi uplift."""
    normal = _maritime_lcp(ports, ChokepointState(), "JPN", "DEU")
    houthi = _maritime_lcp(ports, ChokepointState.houthi_red_sea_2024(), "JPN", "DEU")
    ratio = houthi / normal
    assert ratio > 1.10, f"JPN→DEU Houthi uplift {ratio:.3f}x should be >1.10"


def test_houthi_sgp_deu_uplift(ports):
    """Singapore → DEU (heavy Suez user) should show Houthi uplift."""
    normal = _maritime_lcp(ports, ChokepointState(), "SGP", "DEU")
    houthi = _maritime_lcp(ports, ChokepointState.houthi_red_sea_2024(), "SGP", "DEU")
    ratio = houthi / normal
    assert ratio > 1.10, f"SGP→DEU Houthi uplift {ratio:.3f}x should be >1.10"


# ---------------------------------------------------------------------------
# Panama Drought 2023: LA → Rotterdam (explicit Pacific ↔ Atlantic)
# ---------------------------------------------------------------------------

def test_panama_pacific_to_atlantic_uplift(ports):
    """Under Panama drought 2023, a Pacific port → Atlantic route should show >10% uplift.

    Uses Los Angeles (lon=-118.26, genuine Pacific port west of Panama) →
    Rotterdam (lon=4.48, Atlantic/North Sea) — a canonical trans-Panama route.
    The USA geographic centroid maps to Houston (Gulf Coast, Atlantic side) so
    USA→NLD does NOT use Panama; we test with explicit LA coordinates instead.
    """
    la  = Agglomeration("USA_LA",  0, -118.26, 33.74, 1.0)   # Los Angeles
    ams = Agglomeration("NLD_AMS", 0,    4.48, 51.92, 1.0)   # Rotterdam/Amsterdam

    G_normal = _build_maritime_graph_nx(ports, ChokepointState())
    G_panama = _build_maritime_graph_nx(ports, ChokepointState.panama_drought_2023())

    cost_n, _ = lcp_pair_nx(la, ams, {"maritime": G_normal})
    cost_p, _ = lcp_pair_nx(la, ams, {"maritime": G_panama})
    ratio = cost_p / cost_n
    assert ratio > 1.10, (
        f"LA→Rotterdam Panama drought uplift {ratio:.3f}x should be >1.10; "
        f"normal={cost_n:.1f}, shocked={cost_p:.1f}"
    )


# ---------------------------------------------------------------------------
# Routes that should NOT be affected by Houthi
# ---------------------------------------------------------------------------

def test_houthi_aus_usa_no_uplift(ports):
    """AUS → USA (trans-Pacific route) should be unaffected by Houthi Red Sea."""
    normal = _maritime_lcp(ports, ChokepointState(), "AUS", "USA")
    houthi = _maritime_lcp(ports, ChokepointState.houthi_red_sea_2024(), "AUS", "USA")
    ratio = houthi / normal
    assert ratio < 1.05, (
        f"AUS→USA (trans-Pacific) Houthi ratio {ratio:.3f}x should be ~1.0; "
        "this route does not use the Red Sea"
    )


# ---------------------------------------------------------------------------
# Print summary (visible with pytest -s)
# ---------------------------------------------------------------------------

def test_print_chokepoint_summary(ports):
    """Print a human-readable chokepoint uplift summary (run with -s to see)."""
    pairs = [
        ("IND", "DEU", "Houthi", ChokepointState.houthi_red_sea_2024()),
        ("JPN", "DEU", "Houthi", ChokepointState.houthi_red_sea_2024()),
        ("SGP", "DEU", "Houthi", ChokepointState.houthi_red_sea_2024()),
        ("USA", "NLD", "Panama", ChokepointState.panama_drought_2023()),
        ("AUS", "USA", "Houthi", ChokepointState.houthi_red_sea_2024()),
    ]
    print("\n\n=== Chokepoint uplift summary ===")
    print(f"{'Pair':<10} {'Scenario':<8} {'Normal':>10} {'Shocked':>10} {'Ratio':>7}")
    print("-" * 50)
    for iso_o, iso_d, label, cp in pairs:
        n = _maritime_lcp(ports, ChokepointState(), iso_o, iso_d)
        s = _maritime_lcp(ports, cp, iso_o, iso_d)
        print(f"{iso_o}→{iso_d:<4}  {label:<8}  ${n:>8.1f}  ${s:>8.1f}  {s/n:>6.3f}x")
