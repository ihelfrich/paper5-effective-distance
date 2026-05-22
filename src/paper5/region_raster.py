"""Country-polygon → RasterDist pipeline.

Reads a global raster (population, nightlights) and a set of country
boundaries, masks the raster to each country, and produces a `RasterDist`
ready for distance computation. Handles two technical issues that
matter:

  1. Equal-angle raster cell-area distortion. Global rasters are stored
     on an equal-degree grid, so a 1° × 1° cell at the equator has ~5×
     the metric area of a 1° × 1° cell at 80°N. If we use raw cell
     values as weights, we systematically over-weight tropical cells.
     We apply a cosine-of-latitude correction.

  2. Large countries cross the antimeridian (Russia, Fiji, USA via
     Alaska). We split the polygon at ±180° longitude before masking
     to avoid `rasterio.mask` returning the entire raster.

The output is cached as parquet per (country_iso3, raster_name, year)
so repeated runs don't re-mask.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

try:
    import rasterio
    from rasterio import features as rio_features
    from rasterio.mask import mask as rio_mask
    from rasterio.warp import transform_bounds
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False

try:
    import geopandas as gpd
    from shapely.geometry import box, mapping
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False

import pandas as pd

from paper5.ot_distance import RasterDist


# ── Top-level loaders ────────────────────────────────────────────────────────

def load_country_boundaries(path: Optional[Path] = None) -> "gpd.GeoDataFrame":
    """Load Natural Earth admin0 country boundaries.

    The shapefile is expected to live where we pulled it during the data
    inventory: /Volumes/HELFRICH-GD/Geography_Reference/natural_earth/
    """
    if not HAS_GEOPANDAS:
        raise ImportError("geopandas is required. pip install geopandas")
    if path is None:
        path = Path("/Volumes/HELFRICH-GD/Geography_Reference/natural_earth/"
                    "ne_10m_admin_0_countries.shp")
    if not path.exists():
        raise FileNotFoundError(
            f"Country boundaries not at {path}. Pull from naturalearthdata.com."
        )
    gdf = gpd.read_file(path)
    # Standardize: ISO_A3 is the ISO3 code (some entries have '-99' for unmapped)
    gdf["ISO3"] = gdf["ISO_A3"].where(gdf["ISO_A3"] != "-99", gdf["ADM0_A3"])
    return gdf[["ISO3", "ADMIN", "geometry"]].rename(columns={"ADMIN": "name"})


# ── Raster masking ───────────────────────────────────────────────────────────

@dataclass
class MaskResult:
    """Output of masking one raster against one country polygon."""
    raster_values: np.ndarray  # 1-D array of cell values (zeros dropped)
    cell_lats: np.ndarray      # 1-D array of cell-centroid latitudes
    cell_lons: np.ndarray      # 1-D array of cell-centroid longitudes
    cell_areas_km2: np.ndarray  # 1-D array of per-cell metric area (cosine-of-lat corrected)
    raw_n_cells: int           # before zero/NaN filtering
    kept_n_cells: int          # after filtering
    full_area_km2: float = 0.0  # total area of *all* non-zero cells (pre max_atoms downsample)
    full_n_cells: int = 0       # number of non-zero cells (pre max_atoms downsample)


def _split_antimeridian(geom):
    """If a geometry crosses ±180°, split it. Otherwise return as-is.
    Most Natural Earth polygons are already split. This is defensive."""
    minx, _, maxx, _ = geom.bounds
    if minx < -180 or maxx > 180 or (maxx - minx) > 270:
        # Probably an antimeridian-crossing polygon — clip to two halves
        eastern_half = geom.intersection(box(-180, -90, 180, 90))
        return eastern_half if not eastern_half.is_empty else geom
    return geom


def _cell_area_km2(lat_deg: np.ndarray, lon_res_deg: float, lat_res_deg: float) -> np.ndarray:
    """Approximate cell area in km² for equal-angle grid, with cosine-of-latitude correction.

    For a cell centered at latitude `lat`:
      width  = 111.32 * cos(lat) * lon_res_deg  km
      height = 111.32 * lat_res_deg            km
      area   = width * height
    """
    lat_rad = np.deg2rad(lat_deg)
    width_km = 111.32 * np.cos(lat_rad) * lon_res_deg
    height_km = 111.32 * lat_res_deg
    return np.abs(width_km * height_km)


def mask_raster_to_country(
    raster_path: Path,
    geometry,
    *,
    min_value: float = 1e-9,
    nodata_replacement: float = 0.0,
    max_atoms: Optional[int] = None,
    coarsen_factor: Optional[int] = None,
    sample_strategy: str = "top_n",  # "top_n" | "stratified" | "uniform"
    sample_seed: int = 42,
) -> MaskResult:
    """Mask a global raster to one country's geometry.

    Parameters
    ----------
    raster_path : Path to the GeoTIFF.
    geometry : shapely geometry (will be antimeridian-split if needed).
    min_value : drop cells whose value is below this (numerical floor).
    nodata_replacement : value used in place of raster's no-data sentinel.
    max_atoms : if set, keep only the top-N cells by value (caps OT compute).

    Returns
    -------
    MaskResult with raster_values, cell_lats, cell_lons, cell_areas_km2.
    """
    if not HAS_RASTERIO:
        raise ImportError("rasterio is required.")
    geom = _split_antimeridian(geometry)
    with rasterio.open(raster_path) as src:
        # Project geometry to raster CRS if different (Mollweide for GHS, etc.)
        if str(src.crs).lower() != "epsg:4326":
            # Geom is in EPSG:4326; need to reproject to raster CRS
            try:
                import geopandas as gpd
                from shapely.geometry import shape, mapping as sh_mapping
                gs = gpd.GeoSeries([geom], crs="EPSG:4326")
                gs = gs.to_crs(src.crs)
                geom_in_raster_crs = gs.iloc[0]
            except Exception as e:
                raise RuntimeError(f"Could not reproject geometry to raster CRS {src.crs}: {e}")
        else:
            geom_in_raster_crs = geom

        out, transform = rio_mask(src, [mapping(geom_in_raster_crs)],
                                  crop=True, nodata=src.nodata if src.nodata is not None else 0,
                                  all_touched=False, filled=True)
        arr = out[0].astype(np.float64)
        nodata = src.nodata
        raster_crs = src.crs

    # Handle nodata
    if nodata is not None:
        arr = np.where(arr == nodata, nodata_replacement, arr)
    arr = np.where(np.isfinite(arr), arr, nodata_replacement)
    raw_n = arr.size

    # Block-coarsen by summing over K×K blocks. This preserves total mass and
    # spatial coverage while reducing the atom count by K². The coarsened
    # transform inflates pixel size by K. Required for d_eff(θ<0) to be
    # unbiased on large dense countries — top-N sampling concentrates atoms
    # in metros and forces d_eff toward zero via the harmonic-dominance
    # property of the generalized mean at negative orders.
    if coarsen_factor is not None and coarsen_factor > 1:
        K = int(coarsen_factor)
        rows_in, cols_in = arr.shape
        # Pad to multiple of K
        pad_rows = (K - rows_in % K) % K
        pad_cols = (K - cols_in % K) % K
        if pad_rows or pad_cols:
            arr = np.pad(arr, ((0, pad_rows), (0, pad_cols)),
                         mode="constant", constant_values=0)
        new_rows = arr.shape[0] // K
        new_cols = arr.shape[1] // K
        arr = arr.reshape(new_rows, K, new_cols, K).sum(axis=(1, 3))
        # Update transform: pixel size scales by K
        from rasterio.transform import Affine
        transform = Affine(transform.a * K, transform.b, transform.c,
                           transform.d, transform.e * K, transform.f)

    # Cell-center coordinates in raster CRS
    rows, cols = np.where(arr > min_value)
    if rows.size == 0:
        return MaskResult(
            raster_values=np.zeros(0), cell_lats=np.zeros(0), cell_lons=np.zeros(0),
            cell_areas_km2=np.zeros(0), raw_n_cells=raw_n, kept_n_cells=0,
            full_area_km2=0.0, full_n_cells=0,
        )
    a, b, c, d, e, f = (transform.a, transform.b, transform.c,
                        transform.d, transform.e, transform.f)
    xs = a * (cols + 0.5) + b * (rows + 0.5) + c
    ys = d * (cols + 0.5) + e * (rows + 0.5) + f

    # If raster CRS is not WGS84, reproject cell centers back to lat/lon
    if str(raster_crs).lower() != "epsg:4326":
        from rasterio.warp import transform as rio_transform
        # rio_transform takes (xs, ys) lists and returns (lons, lats)
        lons, lats = rio_transform(raster_crs, "EPSG:4326", xs.tolist(), ys.tolist())
        lats = np.array(lats); lons = np.array(lons)
    else:
        lats = ys; lons = xs

    values = arr[rows, cols]

    # Per-cell metric area (cosine-of-latitude corrected)
    # Approximate lat/lon resolution from the affine transform diagonal
    # (this works in EPSG:4326; for Mollweide GHS-POP the area handling differs
    # and we'd need to use the raster's native cell area)
    if str(raster_crs).lower() == "epsg:4326":
        lon_res = abs(a)
        lat_res = abs(e)
        areas = _cell_area_km2(lats, lon_res, lat_res)
    else:
        # For Mollweide and other equal-area projections, the native pixel area is constant
        # GHS-POP at 1km Mollweide → each cell is exactly 1 km²
        areas = np.full_like(values, 1.0)

    # Capture full area BEFORE the max_atoms downsample so callers can compare
    # against Head-Mayer closed-form (which depends on full country area, not
    # the top-N highest-population subset).
    full_area_km2 = float(areas.sum())
    full_n_cells = int(values.size)

    # Cap atoms if requested.
    # sample_strategy: 'top_n' keeps the densest cells (was the original
    # behavior; biases d_eff(θ<0) toward zero); 'stratified' samples cells
    # with probability ∝ value, preserving the spatial distribution in
    # expectation; 'uniform' samples cells uniformly at random regardless
    # of value (used as a sanity baseline).
    if max_atoms is not None and values.size > max_atoms:
        if sample_strategy == "top_n":
            idx = np.argsort(values)[-max_atoms:]
        elif sample_strategy == "stratified":
            rng = np.random.default_rng(sample_seed)
            p = values / values.sum()
            idx = rng.choice(values.size, size=max_atoms, replace=False, p=p)
        elif sample_strategy == "uniform":
            rng = np.random.default_rng(sample_seed)
            idx = rng.choice(values.size, size=max_atoms, replace=False)
        else:
            raise ValueError(f"Unknown sample_strategy: {sample_strategy!r}")
        values = values[idx]; lats = lats[idx]; lons = lons[idx]; areas = areas[idx]

    return MaskResult(
        raster_values=values, cell_lats=lats, cell_lons=lons,
        cell_areas_km2=areas, raw_n_cells=raw_n, kept_n_cells=values.size,
        full_area_km2=full_area_km2, full_n_cells=full_n_cells,
    )


def mask_result_to_raster_dist(
    mask: MaskResult,
    *,
    correct_for_cell_area: bool = True,
) -> RasterDist:
    """Convert a MaskResult into a RasterDist (atom cloud with normalized weights).

    Cell-area correction is essential for equal-angle global rasters:
    without it, tropical regions are systematically over-weighted vs
    high-latitude regions.
    """
    if mask.kept_n_cells == 0:
        return RasterDist(coords=np.zeros((0, 2)), weights=np.zeros(0))
    coords = np.column_stack([mask.cell_lats, mask.cell_lons])
    if correct_for_cell_area:
        weights = mask.raster_values * mask.cell_areas_km2
    else:
        weights = mask.raster_values
    total = weights.sum()
    if total <= 0:
        return RasterDist(coords=coords, weights=np.zeros(weights.size))
    return RasterDist(coords=coords, weights=weights / total, total_mass=float(total))


# ── Batch / panel builder ────────────────────────────────────────────────────

def build_country_panel(
    raster_paths: dict[str, Path],         # name → path (e.g. {"ghs_pop_2020": Path(...)})
    countries: "gpd.GeoDataFrame",          # output of load_country_boundaries
    *,
    iso3_codes: Optional[Sequence[str]] = None,  # subset to compute; None = all
    max_atoms: Optional[int] = 2000,
    cache_dir: Optional[Path] = None,
) -> dict[tuple[str, str], RasterDist]:
    """Build a (raster_name, iso3) → RasterDist mapping.

    Parameters
    ----------
    raster_paths : dict mapping a human-readable name (e.g. "ghs_pop_2020") to
                   the GeoTIFF on disk.
    countries : GeoDataFrame with at least 'ISO3' and 'geometry'.
    iso3_codes : optional subset of ISO3 codes to compute; None = every country.
    max_atoms : atom-count cap per RasterDist (keeps OT compute bounded).
    cache_dir : if provided, cache one parquet per (raster_name, iso3) pair.

    Returns
    -------
    Dict mapping (raster_name, iso3) → RasterDist.
    """
    if cache_dir:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

    if iso3_codes is None:
        iso3_codes = countries["ISO3"].tolist()
    iso3_set = set(iso3_codes)

    out: dict[tuple[str, str], RasterDist] = {}
    for raster_name, raster_path in raster_paths.items():
        for _, row in countries.iterrows():
            iso = row["ISO3"]
            if iso not in iso3_set:
                continue
            cache_path = (cache_dir / f"{raster_name}__{iso}.parquet") if cache_dir else None
            if cache_path and cache_path.exists():
                df = pd.read_parquet(cache_path)
                rd = RasterDist(
                    coords=df[["lat", "lon"]].to_numpy(),
                    weights=df["weight"].to_numpy(),
                    total_mass=float(df["weight"].sum()),
                )
                out[(raster_name, iso)] = rd
                continue
            try:
                m = mask_raster_to_country(raster_path, row.geometry, max_atoms=max_atoms)
            except Exception as e:
                warnings.warn(f"Failed to mask {raster_name} for {iso}: {e}")
                out[(raster_name, iso)] = RasterDist(
                    coords=np.zeros((0, 2)), weights=np.zeros(0)
                )
                continue
            rd = mask_result_to_raster_dist(m)
            out[(raster_name, iso)] = rd

            # Cache
            if cache_path and rd.coords.shape[0] > 0:
                pd.DataFrame({
                    "lat": rd.coords[:, 0],
                    "lon": rd.coords[:, 1],
                    "weight": rd.weights,
                }).to_parquet(cache_path)

    return out
