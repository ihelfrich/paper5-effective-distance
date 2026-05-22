"""Gravity-panel construction utilities.

Builds a long-format panel of bilateral + intra-national trade flows ready for
gravity estimation. Reusable across all the empirical experiments (notebook 14
border-effect pilot, notebook 15 Wasserstein velocity, future panels).

The panel has one row per (origin, destination, year) with:
  - X_ij_usd       : USD trade flow (BACI for i≠j; Y_i − exports for i=j)
  - log_X          : log of X_ij_usd (only for X > 0)
  - gdp_o, gdp_d   : origin/destination GDP from CEPII Gravity
  - distw_harmonic : CEPII bilateral distance (population-weighted harmonic)
  - distw_arithmetic : CEPII bilateral distance (population-weighted arithmetic)
  - contig, comlang_off, col_dep_ever : standard gravity controls
  - home           : indicator for intra-national obs (i==j)
  - year           : panel year

Caller is responsible for adding any additional distance measures (e.g.,
raster-derived d_ii^eff or d_ij^eff) by joining on iso_o, iso_d, year.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

BACI_DIR_DEFAULT = Path("/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b")
CEPII_GRAVITY_DEFAULT = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")


# ── Loaders ──────────────────────────────────────────────────────────────────

def load_baci_year(year: int, baci_dir: Path = BACI_DIR_DEFAULT) -> pd.DataFrame:
    """Aggregate BACI HS02 bilateral flows to country-pair-year totals (USD)."""
    f = baci_dir / f"BACI_HS02_Y{year}_V202401b.csv"
    if not f.exists():
        raise FileNotFoundError(f"BACI year file missing: {f}")
    df = pd.read_csv(
        f, dtype={"t": int, "i": int, "j": int, "k": str},
        usecols=["t", "i", "j", "v"],
    )
    df = df.dropna(subset=["v"])
    df["v"] = pd.to_numeric(df["v"], errors="coerce")
    df = df.dropna(subset=["v"])
    agg = df.groupby(["i", "j"], as_index=False).agg(X_ij_usd=("v", "sum"))
    agg["X_ij_usd"] *= 1000  # BACI 'v' is in thousands of USD
    agg["year"] = year
    return agg


def load_cepii_gravity_for_year(year: int, path: Path = CEPII_GRAVITY_DEFAULT) -> pd.DataFrame:
    """Filter CEPII Gravity to a single year."""
    df = pd.read_csv(
        path, dtype={"iso3_o": str, "iso3_d": str, "year": int},
    )
    keep = ["year", "iso3_o", "iso3_d", "iso3num_o", "iso3num_d",
            "gdp_o", "gdp_d", "pop_o", "pop_d",
            "dist", "distcap", "distw_arithmetic", "distw_harmonic",
            "contig", "comlang_off", "comlang_ethno", "col_dep_ever",
            "gatt_o", "gatt_d", "rta"]
    keep = [c for c in keep if c in df.columns]
    return df[df["year"] == year][keep].copy()


# ── Country areas ────────────────────────────────────────────────────────────

def compute_country_areas_km2(geometries_path: Optional[Path] = None) -> dict[str, float]:
    """Compute country areas using Natural Earth admin0 polygons reprojected to Mollweide."""
    try:
        import geopandas as gpd
    except ImportError:
        warnings.warn("geopandas not installed; can't compute areas.")
        return {}
    if geometries_path is None:
        geometries_path = Path("/Volumes/HELFRICH-GD/Geography_Reference/natural_earth/"
                               "ne_10m_admin_0_countries.shp")
    if not geometries_path.exists():
        warnings.warn(f"Natural Earth admin0 not at {geometries_path}.")
        return {}
    bnd = gpd.read_file(geometries_path)
    bnd["ISO3"] = bnd["ISO_A3"].where(bnd["ISO_A3"] != "-99", bnd["ADM0_A3"])
    bnd_eq = bnd.to_crs("ESRI:54009")
    bnd_eq["area_km2"] = bnd_eq.geometry.area / 1e6
    return dict(zip(bnd_eq["ISO3"], bnd_eq["area_km2"]))


# ── Panel construction ──────────────────────────────────────────────────────

def build_gravity_panel(
    year: int,
    *,
    include_intra: bool = True,
    intra_distance_default: str = "head_mayer",  # 'head_mayer' or 'cepii_self'
    baci_dir: Path = BACI_DIR_DEFAULT,
    cepii_path: Path = CEPII_GRAVITY_DEFAULT,
) -> pd.DataFrame:
    """Build a single-year gravity panel ready for regression.

    Parameters
    ----------
    year : int
        Panel year (BACI starts at 2002; CEPII Gravity goes to 2021).
    include_intra : bool
        If True, append intra-national rows (one per country) with
        X_ii = max(0, GDP_i - exports_i) as the standard intra-national-flow
        proxy from Wei (1996) / Anderson & van Wincoop (2003).
    intra_distance_default : str
        Default intra-national distance to use:
        - 'head_mayer': d_ii = 0.67 * sqrt(area/π)  (CEPII baseline)
        - 'cepii_self': use CEPII's distw_harmonic for i==j (which CEPII fills
          with the Head-Mayer closed form, so functionally identical).
    baci_dir, cepii_path: data paths.

    Returns
    -------
    DataFrame with one row per (iso_o, iso_d) for the year.
    """
    # Load components
    baci = load_baci_year(year, baci_dir=baci_dir)
    cepii = load_cepii_gravity_for_year(year, path=cepii_path)

    # Map BACI iso-numeric → iso3 letters via CEPII
    iso3num_to_iso3 = dict(cepii[["iso3num_o", "iso3_o"]].drop_duplicates().values)
    iso3num_to_iso3.update(dict(cepii[["iso3num_d", "iso3_d"]].drop_duplicates().values))
    baci["iso3_o"] = baci["i"].map(iso3num_to_iso3)
    baci["iso3_d"] = baci["j"].map(iso3num_to_iso3)
    baci = baci.dropna(subset=["iso3_o", "iso3_d"])

    # Inter-national panel = CEPII (i ≠ j) with BACI flows merged in
    cepii_inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]].copy()
    panel = cepii_inter.merge(
        baci[["iso3_o", "iso3_d", "X_ij_usd"]],
        on=["iso3_o", "iso3_d"], how="left",
    )
    panel["X_ij_usd"] = panel["X_ij_usd"].fillna(0.0)
    panel["home"] = 0

    if include_intra:
        # Compute total exports per country
        exports_by_origin = baci.groupby("iso3_o", as_index=False)["X_ij_usd"].sum()
        exports_by_origin.rename(columns={"X_ij_usd": "exports_usd"}, inplace=True)
        gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)

        intra_rows = []
        areas = compute_country_areas_km2()
        for iso, gdp in gdp_lookup.items():
            if not np.isfinite(gdp) or gdp <= 0:
                continue
            exp = exports_by_origin.loc[exports_by_origin["iso3_o"] == iso, "exports_usd"]
            exp_v = float(exp.iloc[0]) if len(exp) > 0 else 0.0
            x_ii = max(0.0, gdp - exp_v)
            area_km2 = areas.get(iso, np.nan)
            row = {
                "year": year, "iso3_o": iso, "iso3_d": iso,
                "X_ij_usd": x_ii, "home": 1,
                "gdp_o": gdp, "gdp_d": gdp,
                "pop_o": np.nan, "pop_d": np.nan,
                "dist": np.nan, "distcap": np.nan,
                "distw_arithmetic": np.nan, "distw_harmonic": np.nan,
                "contig": 0, "comlang_off": 1, "col_dep_ever": 0,
                "area_km2": area_km2,
            }
            if intra_distance_default == "head_mayer":
                row["distw_harmonic"] = 0.67 * np.sqrt(area_km2 / np.pi) if np.isfinite(area_km2) else np.nan
                row["distw_arithmetic"] = row["distw_harmonic"]
                row["dist"] = row["distw_harmonic"]
            intra_rows.append(row)
        intra_df = pd.DataFrame(intra_rows)
        panel = pd.concat([panel, intra_df], ignore_index=True, sort=False)

    return panel


def add_log_columns(panel: pd.DataFrame) -> pd.DataFrame:
    """Add log_X, log_d, log_gdp_o, log_gdp_d for OLS regression."""
    p = panel.copy()
    valid_X = (p["X_ij_usd"] > 0).fillna(False)
    valid_d = (p["distw_harmonic"] > 0).fillna(False)
    valid_g = (p["gdp_o"] > 0).fillna(False) & (p["gdp_d"] > 0).fillna(False)
    p["log_X"] = np.where(valid_X, np.log(p["X_ij_usd"].clip(lower=1.0)), np.nan)
    p["log_d"] = np.where(valid_d, np.log(p["distw_harmonic"].clip(lower=1.0)), np.nan)
    p["log_gdp_o"] = np.where(valid_g, np.log(p["gdp_o"].clip(lower=1.0)), np.nan)
    p["log_gdp_d"] = np.where(valid_g, np.log(p["gdp_d"].clip(lower=1.0)), np.nan)
    return p


def swap_intra_distance(panel: pd.DataFrame, intra_d: dict[str, float]) -> pd.DataFrame:
    """Swap the intra-national distance values for specific countries.

    Parameters
    ----------
    panel : DataFrame from build_gravity_panel.
    intra_d : dict mapping ISO3 → new d_ii in km (e.g., from raster CES).

    Returns
    -------
    Panel with `distw_harmonic` updated for matching home rows.
    """
    p = panel.copy()
    home_mask = p["home"] == 1
    p.loc[home_mask, "distw_harmonic"] = p.loc[home_mask, "iso3_o"].map(intra_d).fillna(
        p.loc[home_mask, "distw_harmonic"]
    )
    return p
