"""data_loaders.py — standardized loaders for Paper 5 data sources.

All loaders return DataFrames with consistent column naming:
    iso_o, iso_d, year, ...

and handle path resolution via paper5.config (HELFRICH_GD, ONEDRIVE, PAPER5_DATA).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


# Path roots — parameterize via env vars in production; hard-coded for sprint.
HELFRICH_GD = Path("/Volumes/HELFRICH-GD")
ONEDRIVE = Path("/Users/ian/Library/CloudStorage/OneDrive-Personal")
PAPER5_DATA = Path(__file__).resolve().parents[2] / "data"


# ==================== TRADE ====================

def load_baci(
    year: Optional[int] = None,
    aggregation: str = "country-pair-year",
) -> pd.DataFrame:
    """Load BACI HS02 V202401b.

    Parameters
    ----------
    year : if None, return full panel.
    aggregation : one of {"country-pair-year", "country-pair-year-hs2", "raw"}
    """
    src = HELFRICH_GD / "TradeData" / "BACI_HS02_V202401b"
    if year is not None:
        files = [src / f"BACI_HS02_Y{year}_V202401b.csv"]
    else:
        files = sorted(src.glob("BACI_HS02_Y*_V202401b.csv"))

    dfs = [pd.read_csv(f, dtype={"t": int, "i": int, "j": int, "k": str})
           for f in files]
    df = pd.concat(dfs, ignore_index=True)

    # Rename to Paper 5 canonical columns
    df = df.rename(columns={"t": "year", "i": "iso_o_num", "j": "iso_d_num",
                            "k": "hs2", "v": "trade_value", "q": "quantity"})
    # HS6 collapse to HS2 for HS2 aggregation
    if aggregation == "country-pair-year-hs2":
        df["hs2"] = df["hs2"].str[:2]
        df = df.groupby(["year", "iso_o_num", "iso_d_num", "hs2"],
                        as_index=False)[["trade_value", "quantity"]].sum()
    elif aggregation == "country-pair-year":
        df = df.groupby(["year", "iso_o_num", "iso_d_num"],
                        as_index=False)[["trade_value", "quantity"]].sum()
    elif aggregation == "raw":
        pass
    else:
        raise ValueError(aggregation)

    return df


def load_cepii_gravity(
    filtered: bool = False,
    distw_col: str = "distw_harmonic",
) -> pd.DataFrame:
    """Load CEPII Gravity V202211 (or Ian's pre-filtered subset).

    V202211 renamed 'distw' → 'distw_harmonic' and 'distw_arithmetic'.
    We add a 'distw' alias so downstream code using the old name still works.
    """
    if filtered:
        path = ONEDRIVE / "ThesisData_Last" / "Gravity_V202211_filtered.csv"
    else:
        path = HELFRICH_GD / "TradeData" / "Gravity_csv_V202211" / "Gravity_V202211.csv"
    df = pd.read_csv(path, low_memory=False)
    # Alias distw_harmonic → distw for backwards compatibility
    if "distw_harmonic" in df.columns and "distw" not in df.columns:
        df["distw"] = df["distw_harmonic"]
    return df


def load_baci_iso_crosswalk() -> pd.DataFrame:
    """Return a DataFrame mapping BACI/ISO numeric codes → ISO-3166 alpha-3.

    BACI uses ISO 3166-1 numeric codes in columns i, j.
    The CEPII Countries_V202211.csv provides the numeric→alpha3 mapping.

    Returns DataFrame with columns: iso3num (int), iso3 (str).
    """
    path = HELFRICH_GD / "TradeData" / "Gravity_csv_V202211" / "Countries_V202211.csv"
    df = pd.read_csv(path)
    df = df[["iso3", "iso3num"]].dropna(subset=["iso3num"])
    df["iso3num"] = df["iso3num"].astype(int)
    return df.drop_duplicates("iso3num").set_index("iso3num")["iso3"].to_frame()


def load_baci_with_iso(
    year: Optional[int] = None,
    aggregation: str = "country-pair-year",
) -> pd.DataFrame:
    """Like load_baci() but with ISO-3166 alpha-3 columns (iso_o, iso_d) added."""
    df = load_baci(year=year, aggregation=aggregation)
    xwalk = load_baci_iso_crosswalk()
    df = df.merge(xwalk.rename(columns={"iso3": "iso_o"}),
                  left_on="iso_o_num", right_index=True, how="left")
    df = df.merge(xwalk.rename(columns={"iso3": "iso_d"}),
                  left_on="iso_d_num", right_index=True, how="left")
    return df


def load_centrality_panel() -> pd.DataFrame:
    """Load pre-computed network centrality panel (2002–2022)."""
    path = ONEDRIVE / "ThesisData_Last" / "tables" / "centrality_evolution.csv"
    return pd.read_csv(path)


# ==================== POPULATION ====================

def load_worldpop(year: int, resolution: str = "100m") -> "xarray.DataArray":  # type: ignore
    """Lazy-load WorldPop Constrained raster for a given year."""
    import rioxarray  # lazy import
    path = PAPER5_DATA / "worldpop" / f"ppp_{year}_{resolution}_Aggregated.tif"
    return rioxarray.open_rasterio(path, chunks={"x": 4096, "y": 4096})


def load_ghs_pop(year: int) -> "xarray.DataArray":  # type: ignore
    """Lazy-load GHS-POP R2023A raster for epoch."""
    import rioxarray
    path = PAPER5_DATA / "ghs_pop" / f"GHS_POP_E{year}_GLOBE_R2023A_54009_100.tif"
    return rioxarray.open_rasterio(path, chunks={"x": 4096, "y": 4096})


# ==================== NIGHTLIGHTS ====================

def load_viirs_annual(year: int) -> "xarray.DataArray":  # type: ignore
    """VIIRS V2 annual composite (2012–)."""
    import rioxarray
    path = PAPER5_DATA / "viirs" / f"VNL_v2_npp_{year}_global_vcmslcfg_c202101211500.average_masked.tif"
    return rioxarray.open_rasterio(path, chunks={"x": 4096, "y": 4096})


def load_dmsp_harmonized(year: int) -> "xarray.DataArray":  # type: ignore
    """Li et al. (2020) DMSP-VIIRS harmonized series (1992–2018)."""
    import rioxarray
    path = PAPER5_DATA / "dmsp_harmonized" / f"Harmonized_DN_NTL_{year}_calDMSP.tif"
    return rioxarray.open_rasterio(path, chunks={"x": 4096, "y": 4096})


# ==================== TRANSPORT ====================

def load_osm_planet(region: Optional[str] = None) -> Path:
    """Return path to OSM PBF for full planet or a Geofabrik region."""
    if region is None:
        return PAPER5_DATA / "osm" / "planet-latest.osm.pbf"
    return PAPER5_DATA / "osm" / f"{region}-latest.osm.pbf"


def load_wpi_ports() -> pd.DataFrame:
    """World Port Index — ~3,700 global ports with coordinates and attributes."""
    path = PAPER5_DATA / "wpi" / "WPI.csv"
    return pd.read_csv(path)


def load_openflights_airports() -> pd.DataFrame:
    """OpenFlights airport database (~10K airports, coordinates)."""
    path = PAPER5_DATA / "openflights" / "airports.dat"
    cols = ["airport_id", "name", "city", "country", "iata", "icao",
            "lat", "lon", "altitude", "timezone", "dst", "tz", "type", "source"]
    return pd.read_csv(path, header=None, names=cols)


def load_cerdi_sea_distance() -> pd.DataFrame:
    """CERDI port-pair sea distances."""
    path = PAPER5_DATA / "cerdi_sea" / "SeaDistance.csv"
    return pd.read_csv(path)
