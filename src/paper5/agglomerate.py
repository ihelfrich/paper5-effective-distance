"""agglomerate.py — WorldPop/GHS-POP → within-country agglomeration centroids.

Clustering strategy (v0.2 — implemented):
    1. Open WorldPop 100m raster; clip to country bounding box; mask to polygon.
    2. Compute per-cell population density in pp/km² (accounts for lat-varying
       cell area).
    3. Threshold: cells < density_threshold set to zero.
    4. Connected components via scipy.ndimage.label on the binary threshold mask.
    5. Merge proto-agglomerations whose centroids are within merge_distance_km of
       each other (union-find by nearest-centroid iteration).
    6. Rank merged agglomerations by population; accumulate until
       coverage_target (default 80%) of national total.
    7. For each retained agglomeration compute:
       - Population-weighted centroid (lon, lat)
       - Sum VIIRS annual radiance over the footprint (requires separate raster)
       - Coverage share of national population.

v1.0 upgrade path: replace the thresholding+CC step with direct use of
GHS-UCDB Urban Centres; keep WorldPop for cross-validation.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class AgglomerationConfig:
    """Tuning parameters for the clustering algorithm."""

    density_threshold_pp_km2: float = 300.0
    """Minimum population density to be included in an agglomeration seed."""

    merge_distance_km: float = 25.0
    """Proto-agglomerations whose centroids are within this distance are merged."""

    coverage_target: float = 0.80
    """Stop adding agglomerations once this share of national population is covered."""

    max_agglomerations_per_country: int = 50
    """Hard cap regardless of coverage."""

    min_pop_threshold: float = 5_000.0
    """Discard proto-agglomerations below this absolute population count."""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _cell_area_km2(lat_deg: np.ndarray, res_deg: float = 100 / 111_000) -> np.ndarray:
    """Approximate cell area (km²) at a given latitude.

    WorldPop 100m has nominal resolution of ~8.3e-4 degrees at the equator.
    Area scales as cos(lat) in the x-direction; y-direction is constant.
    """
    lat_rad = np.deg2rad(lat_deg)
    # Earth radius ≈ 6371 km
    dy_km = res_deg * 111.32  # km per degree of latitude (constant)
    dx_km = res_deg * 111.32 * np.cos(lat_rad)  # shrinks toward poles
    return dy_km * dx_km


def _great_circle_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """Haversine great-circle distance in km."""
    R = 6371.0
    dlat = np.deg2rad(lat2 - lat1)
    dlon = np.deg2rad(lon2 - lon1)
    a = (np.sin(dlat / 2) ** 2
         + np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.sin(dlon / 2) ** 2)
    return 2 * R * np.arcsin(np.sqrt(a))


def _merge_components(
    centroids: list[tuple[float, float, float]],  # (lon, lat, pop)
    merge_distance_km: float,
) -> list[list[int]]:
    """Union-find merge: return groups of proto-agglomeration indices.

    Uses a simple O(n²) nearest-centroid merge; n ≤ ~500 per country so fine.
    """
    n = len(centroids)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in range(n):
        for j in range(i + 1, n):
            d = _great_circle_km(
                centroids[i][0], centroids[i][1],
                centroids[j][0], centroids[j][1],
            )
            if d <= merge_distance_km:
                union(i, j)

    groups: dict[int, list[int]] = {}
    for idx in range(n):
        root = find(idx)
        groups.setdefault(root, []).append(idx)
    return list(groups.values())


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_country_agglomerations(
    iso3: str,
    worldpop_path: Path,
    viirs_path: Optional[Path],
    boundary_geom,  # shapely Polygon/MultiPolygon for this country
    config: AgglomerationConfig,
    year: int,
) -> pd.DataFrame:
    """Identify agglomerations for a single country.

    Parameters
    ----------
    iso3 : ISO-3166 alpha-3.
    worldpop_path : Path to WorldPop 100m GeoTIFF for this year.
    viirs_path : Path to VIIRS annual composite GeoTIFF (may be None).
    boundary_geom : Shapely geometry for this country.
    config : AgglomerationConfig.
    year : Used to stamp output.

    Returns
    -------
    DataFrame with columns:
        iso3, agglom_id, lon, lat, pop, viirs_radiance, coverage_share, year
    """
    try:
        import rioxarray  # noqa: F401
        import xarray as xr
        from rasterio.features import geometry_mask
        import rasterio
    except ImportError as e:
        raise ImportError(
            "agglomerate.py requires rioxarray, rasterio, and xarray. "
            "Install with: pip install rioxarray rasterio xarray"
        ) from e

    try:
        from scipy.ndimage import label as cc_label
    except ImportError as e:
        raise ImportError("agglomerate.py requires scipy. pip install scipy") from e

    import rioxarray as rxr

    # ------------------------------------------------------------------
    # 1. Load and clip WorldPop raster to bounding box + mask to polygon
    # ------------------------------------------------------------------
    pop_da = rxr.open_rasterio(worldpop_path, masked=True).squeeze("band", drop=True)

    # Clip to bounding box first (fast; avoids loading global raster fully)
    minx, miny, maxx, maxy = boundary_geom.bounds
    pad = 0.05  # degrees — tiny buffer to avoid edge clipping
    pop_clipped = pop_da.sel(
        x=slice(minx - pad, maxx + pad),
        y=slice(maxy + pad, miny - pad),  # y is descending in GeoTIFFs
    )

    if pop_clipped.size == 0:
        warnings.warn(f"[{iso3}] Empty raster clip — bounding box outside raster coverage.")
        return pd.DataFrame(columns=["iso3", "agglom_id", "lon", "lat", "pop",
                                     "viirs_radiance", "coverage_share", "year"])

    pop_arr = pop_clipped.values.astype(np.float32)
    pop_arr = np.where(np.isnan(pop_arr) | (pop_arr < 0), 0.0, pop_arr)

    # Pixel coordinates
    lons = pop_clipped.x.values  # 1D
    lats = pop_clipped.y.values  # 1D, descending

    # ------------------------------------------------------------------
    # 2. Build polygon mask (True = inside country boundary)
    # ------------------------------------------------------------------
    with rasterio.open(worldpop_path) as src:
        transform = src.transform
    # Re-derive transform for the clipped array
    from rasterio.transform import from_bounds
    clip_transform = from_bounds(
        lons[0] - abs(lons[1] - lons[0]) / 2,
        lats[-1] - abs(lats[1] - lats[0]) / 2,
        lons[-1] + abs(lons[1] - lons[0]) / 2,
        lats[0] + abs(lats[1] - lats[0]) / 2,
        len(lons), len(lats),
    )
    country_mask = ~geometry_mask(
        [boundary_geom],
        transform=clip_transform,
        invert=False,
        out_shape=(len(lats), len(lons)),
    )
    pop_arr = pop_arr * country_mask

    # ------------------------------------------------------------------
    # 3. Density threshold (pp/km²)
    # ------------------------------------------------------------------
    lat_grid = np.broadcast_to(lats[:, np.newaxis], pop_arr.shape)
    res_deg = abs(float(lons[1] - lons[0])) if len(lons) > 1 else 100 / 111_000
    cell_area = _cell_area_km2(lat_grid, res_deg=res_deg)
    density = np.where(cell_area > 0, pop_arr / cell_area, 0.0)

    threshold_mask = (density >= config.density_threshold_pp_km2)

    # ------------------------------------------------------------------
    # 4. Connected components on threshold mask
    # ------------------------------------------------------------------
    labeled_arr, n_components = cc_label(threshold_mask)

    if n_components == 0:
        # Country has no dense agglomerations at this threshold — fall back to
        # the single population-weighted centroid of the whole country.
        total_pop = float(pop_arr.sum())
        if total_pop == 0:
            return pd.DataFrame(columns=["iso3", "agglom_id", "lon", "lat", "pop",
                                         "viirs_radiance", "coverage_share", "year"])
        lon_grid = np.broadcast_to(lons[np.newaxis, :], pop_arr.shape)
        w_lon = float(np.sum(lon_grid * pop_arr)) / total_pop
        w_lat = float(np.sum(lat_grid * pop_arr)) / total_pop
        return pd.DataFrame([{
            "iso3": iso3, "agglom_id": 0, "lon": w_lon, "lat": w_lat,
            "pop": total_pop, "viirs_radiance": 0.0,
            "coverage_share": 1.0, "year": year,
        }])

    # ------------------------------------------------------------------
    # 5. Compute proto-agglomeration statistics (centroid + pop)
    # ------------------------------------------------------------------
    lon_grid = np.broadcast_to(lons[np.newaxis, :], pop_arr.shape)
    proto: list[tuple[float, float, float]] = []  # (lon, lat, pop)

    for comp_id in range(1, n_components + 1):
        mask_comp = labeled_arr == comp_id
        comp_pop = float(pop_arr[mask_comp].sum())
        if comp_pop < config.min_pop_threshold:
            continue
        # Population-weighted centroid
        w = pop_arr[mask_comp]
        c_lon = float(np.sum(lon_grid[mask_comp] * w)) / comp_pop
        c_lat = float(np.sum(lat_grid[mask_comp] * w)) / comp_pop
        proto.append((c_lon, c_lat, comp_pop))

    if not proto:
        return pd.DataFrame(columns=["iso3", "agglom_id", "lon", "lat", "pop",
                                     "viirs_radiance", "coverage_share", "year"])

    # ------------------------------------------------------------------
    # 6. Merge nearby proto-agglomerations
    # ------------------------------------------------------------------
    groups = _merge_components(proto, config.merge_distance_km)
    merged: list[dict] = []
    for grp in groups:
        grp_pop = sum(proto[i][2] for i in grp)
        # Pop-weighted merged centroid
        m_lon = sum(proto[i][0] * proto[i][2] for i in grp) / grp_pop
        m_lat = sum(proto[i][1] * proto[i][2] for i in grp) / grp_pop
        merged.append({"lon": m_lon, "lat": m_lat, "pop": grp_pop})

    merged.sort(key=lambda d: d["pop"], reverse=True)

    # ------------------------------------------------------------------
    # 7. Retain until coverage target
    # ------------------------------------------------------------------
    national_pop = float(pop_arr.sum())
    cumulative = 0.0
    retained: list[dict] = []
    for agglom in merged:
        if (cumulative / national_pop >= config.coverage_target
                and len(retained) >= 1):
            break
        if len(retained) >= config.max_agglomerations_per_country:
            break
        retained.append(agglom)
        cumulative += agglom["pop"]

    # ------------------------------------------------------------------
    # 8. VIIRS radiance aggregation (optional)
    # ------------------------------------------------------------------
    viirs_sums: list[float] = [0.0] * len(retained)
    if viirs_path is not None and viirs_path.exists():
        try:
            viirs_da = rxr.open_rasterio(viirs_path, masked=True).squeeze("band", drop=True)
            viirs_clipped = viirs_da.sel(
                x=slice(minx - pad, maxx + pad),
                y=slice(maxy + pad, miny - pad),
            )
            if viirs_clipped.size > 0:
                viirs_arr = viirs_clipped.values.astype(np.float32)
                viirs_arr = np.where(np.isnan(viirs_arr) | (viirs_arr < 0), 0.0, viirs_arr)
                # Interpolate VIIRS pixel coordinates to the pop raster grid
                # (simple nearest-neighbor since both are ~100m)
                viirs_lons = viirs_clipped.x.values
                viirs_lats = viirs_clipped.y.values
                for idx, agglom in enumerate(retained):
                    # Find bounding pixels within merge_distance_km of centroid
                    lon_c, lat_c = agglom["lon"], agglom["lat"]
                    deg_radius = config.merge_distance_km / 111.32
                    lon_mask = np.abs(viirs_lons - lon_c) <= deg_radius
                    lat_mask = np.abs(viirs_lats - lat_c) <= deg_radius
                    patch = viirs_arr[np.ix_(lat_mask, lon_mask)]
                    viirs_sums[idx] = float(patch.sum())
        except Exception as exc:
            warnings.warn(f"[{iso3}] VIIRS aggregation failed: {exc}")

    # ------------------------------------------------------------------
    # 9. Build output DataFrame
    # ------------------------------------------------------------------
    rows = []
    for idx, agglom in enumerate(retained):
        rows.append({
            "iso3": iso3,
            "agglom_id": idx,
            "lon": round(agglom["lon"], 5),
            "lat": round(agglom["lat"], 5),
            "pop": round(agglom["pop"], 0),
            "viirs_radiance": round(viirs_sums[idx], 4),
            "coverage_share": round(cumulative / national_pop if national_pop > 0 else 1.0, 6),
            "year": year,
        })
    return pd.DataFrame(rows)


def build_all_agglomerations(
    year: int,
    iso_list: list[str],
    worldpop_path: Path,
    viirs_path: Path,
    boundaries_path: Path,
    config: Optional[AgglomerationConfig] = None,
    out_path: Optional[Path] = None,
) -> pd.DataFrame:
    """Build agglomeration table for all `iso_list` countries in `year`.

    Reads country polygons from `boundaries_path` (GeoPackage or shapefile;
    must have an "ISO_A3" or "GID_0" column for country codes).

    Parameters
    ----------
    year : e.g. 2010.
    iso_list : List of ISO-3166 alpha-3 codes.
    worldpop_path : Full path to WorldPop GeoTIFF for this year.
    viirs_path : Full path to VIIRS composite GeoTIFF for this year.
    boundaries_path : Path to GeoPackage / shapefile with country polygons.
    config : AgglomerationConfig (defaults if None).
    out_path : If provided, write Parquet here.

    Returns
    -------
    DataFrame with one row per (iso3, agglom_id, year).
    """
    try:
        import geopandas as gpd
    except ImportError as e:
        raise ImportError("build_all_agglomerations requires geopandas.") from e

    if config is None:
        config = AgglomerationConfig()

    # Load boundary file
    gdf = gpd.read_file(boundaries_path)
    # Normalize country-code column name
    iso_col = next(
        (c for c in gdf.columns if c.upper() in ("ISO_A3", "GID_0", "ISO3", "ISO_3166_1_ALPHA_3")),
        None,
    )
    if iso_col is None:
        raise ValueError(
            f"Could not find ISO-3 column in {boundaries_path}. "
            f"Found columns: {list(gdf.columns)}"
        )
    gdf = gdf.rename(columns={iso_col: "iso3"}).set_index("iso3")

    frames: list[pd.DataFrame] = []
    missing: list[str] = []

    for iso3 in iso_list:
        if iso3 not in gdf.index:
            missing.append(iso3)
            continue
        geom = gdf.loc[iso3, "geometry"]
        if hasattr(geom, "__geo_interface__") is False and not hasattr(geom, "geoms"):
            missing.append(iso3)
            continue

        try:
            df_country = build_country_agglomerations(
                iso3=iso3,
                worldpop_path=worldpop_path,
                viirs_path=viirs_path,
                boundary_geom=geom,
                config=config,
                year=year,
            )
            frames.append(df_country)
        except Exception as exc:
            warnings.warn(f"[{iso3}] agglomeration failed: {exc!r}")
            missing.append(iso3)

    if missing:
        warnings.warn(f"Skipped {len(missing)} countries: {missing[:10]}{'...' if len(missing)>10 else ''}")

    if not frames:
        raise RuntimeError(f"No agglomerations built for year={year}. Check paths and ISO list.")

    out_df = pd.concat(frames, ignore_index=True)

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_parquet(out_path, index=False)

    return out_df
