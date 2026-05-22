"""test_agglomerate.py — unit + integration tests for agglomerate.py.

Unit tests run with synthetic rasters (no real data needed).
Integration tests are marked with @pytest.mark.integration and
require WorldPop + boundaries on disk; they skip automatically in CI.
"""
from __future__ import annotations

import numpy as np
import pytest

from paper5.agglomerate import (
    AgglomerationConfig,
    _great_circle_km,
    _merge_components,
    _cell_area_km2,
)


# ---------------------------------------------------------------------------
# _cell_area_km2
# ---------------------------------------------------------------------------

def test_cell_area_equator():
    """At equator, 100m × 100m cell ≈ 0.01 km²."""
    lats = np.array([0.0])
    res_deg = 100 / 111_000  # 100m in degrees
    areas = _cell_area_km2(lats, res_deg=res_deg)
    # Should be close to (100m)² = 0.01 km²
    assert areas[0] == pytest.approx(0.01, rel=0.05)


def test_cell_area_shrinks_at_pole():
    """Cell area at lat=60° should be roughly half the equatorial area."""
    res_deg = 100 / 111_000
    area_eq = _cell_area_km2(np.array([0.0]), res_deg=res_deg)[0]
    area_60 = _cell_area_km2(np.array([60.0]), res_deg=res_deg)[0]
    # cos(60°) = 0.5, so area should be ≈ half
    assert area_60 / area_eq == pytest.approx(0.5, rel=0.02)


# ---------------------------------------------------------------------------
# _great_circle_km
# ---------------------------------------------------------------------------

def test_merge_gc_distance():
    d = _great_circle_km(0.0, 0.0, 10.0, 0.0)
    assert d == pytest.approx(1113.2, rel=0.01)


# ---------------------------------------------------------------------------
# _merge_components — union-find merging
# ---------------------------------------------------------------------------

def test_merge_nearby():
    """Two centroids 10 km apart should merge when threshold=25 km."""
    # (lon, lat, pop) — 10 km apart along equator
    centroids = [(0.0, 0.0, 100.0), (0.09, 0.0, 200.0)]
    groups = _merge_components(centroids, merge_distance_km=25.0)
    assert len(groups) == 1
    assert sorted(groups[0]) == [0, 1]


def test_merge_distant():
    """Two centroids 500 km apart should NOT merge at 25 km threshold."""
    centroids = [(0.0, 0.0, 100.0), (5.0, 0.0, 200.0)]
    groups = _merge_components(centroids, merge_distance_km=25.0)
    assert len(groups) == 2


def test_merge_three_chain():
    """A-B-C where d(A,B)=20km, d(B,C)=20km, d(A,C)=40km.
    All three should merge (chain through B)."""
    # Place along equator: ~0.18 deg ≈ 20 km
    centroids = [(0.0, 0.0, 1.0), (0.18, 0.0, 1.0), (0.36, 0.0, 1.0)]
    groups = _merge_components(centroids, merge_distance_km=25.0)
    assert len(groups) == 1
    assert len(groups[0]) == 3


# ---------------------------------------------------------------------------
# AgglomerationConfig defaults
# ---------------------------------------------------------------------------

def test_config_defaults():
    cfg = AgglomerationConfig()
    assert cfg.density_threshold_pp_km2 == 300.0
    assert cfg.merge_distance_km == 25.0
    assert cfg.coverage_target == 0.80
    assert cfg.max_agglomerations_per_country == 50
    assert cfg.min_pop_threshold == 5_000.0


# ---------------------------------------------------------------------------
# Integration test — skipped unless WorldPop is present
# ---------------------------------------------------------------------------

@pytest.mark.integration
def test_build_country_agglomerations_usa():
    """Integration: build agglomerations for USA using WorldPop 2010."""
    from pathlib import Path
    worldpop_path = Path("data/worldpop/ppp_2010_100m_Aggregated.tif")
    boundaries_path = Path("data/boundaries/gadm_410.gpkg")

    if not worldpop_path.exists() or not boundaries_path.exists():
        pytest.skip("WorldPop + boundaries not available.")

    import geopandas as gpd
    from paper5.agglomerate import build_country_agglomerations

    gdf = gpd.read_file(boundaries_path)
    iso_col = next(c for c in gdf.columns if c.upper() in ("GID_0", "ISO_A3", "ISO3"))
    gdf = gdf.rename(columns={iso_col: "iso3"}).set_index("iso3")

    if "USA" not in gdf.index:
        pytest.skip("USA not found in boundary file.")

    cfg = AgglomerationConfig(max_agglomerations_per_country=10)
    df = build_country_agglomerations(
        iso3="USA",
        worldpop_path=worldpop_path,
        viirs_path=None,
        boundary_geom=gdf.loc["USA", "geometry"],
        config=cfg,
        year=2010,
    )

    assert len(df) > 0, "Should find at least one agglomeration for USA"
    assert len(df) <= 10
    assert set(df.columns).issuperset({"iso3","agglom_id","lon","lat","pop","coverage_share"})
    # All centroids should be within rough USA bounding box
    assert df["lon"].between(-180, -60).all()
    assert df["lat"].between(18, 72).all()
    # Population-weighted centroid should be in continental US range
    total_pop = df["pop"].sum()
    assert total_pop > 1e8, "USA should have >100M population in WorldPop"
    print(f"\nUSA agglomerations ({len(df)}, coverage={df['coverage_share'].iloc[-1]:.1%}):")
    print(df[["agglom_id","lon","lat","pop"]].to_string(index=False))
