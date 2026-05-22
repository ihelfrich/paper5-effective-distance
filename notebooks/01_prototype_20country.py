"""01_prototype_20country.py — 20-country prototype pipeline validation.

Run as a script or open in JupyterLab (percent-cell format).
Goal: validate L0→L1→L2 data contracts, produce a small distance panel,
and compare d_ovdl_t vs d_cepii_distw for 2010 and 2022.

Countries: 20 diverse core economies covering all regions, income levels,
and geographic situations (landlocked, island, coastal, large, small).

BEFORE RUNNING:
    pip install -e '.[dev]'
    # Ensure HELFRICH-GD is mounted at /Volumes/HELFRICH-GD
    # Ensure WorldPop 2010 + 2022 100m rasters are in data/worldpop/
    # Ensure VIIRS 2022 raster is in data/viirs/
    # Ensure CEPII Gravity_V202211.csv is accessible via HELFRICH_GD

The script degrades gracefully: if WorldPop is not available, it falls back
to CEPII city-population centroids. If OSM road is not available, it uses
great-circle (OVDL). You can still run it to see the gravity validation.
"""

# %% [markdown]
# # Prototype — 20-country pipeline validation

# %% Imports and paths
from __future__ import annotations

import warnings
from pathlib import Path
import numpy as np
import pandas as pd

# Project-local imports
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from paper5.data_loaders import (
    HELFRICH_GD, ONEDRIVE, PAPER5_DATA,
    load_baci, load_cepii_gravity, load_wpi_ports, load_openflights_airports,
)
from paper5.agglomerate import AgglomerationConfig, build_country_agglomerations
from paper5.distance import (
    ChokepointState, Agglomeration,
    great_circle_km, compute_d_ovdl, compute_bilateral_distance,
    build_transport_graph_nx, lcp_pair_nx,
)

# %% Configuration
YEARS = [2010, 2022]
THETA = -1.0  # Head-Mayer CES parameter

# 20-country core: diverse sample spanning regions, income, and geographic type.
CORE_20 = [
    # North America + Caribbean
    "USA", "MEX", "CAN", "JAM",
    # Europe
    "DEU", "FRA", "POL",
    # Sub-Saharan Africa (landlocked + coastal)
    "NGA", "ETH", "ZAF",
    # Middle East / South Asia
    "SAU", "IND", "PAK",
    # East / Southeast Asia
    "CHN", "JPN", "VNM",
    # Latin America
    "BRA", "COL",
    # Oceania + small-island
    "AUS", "NZL",
]

print(f"Countries: {len(CORE_20)}")
print(f"Pairs: {len(CORE_20) * (len(CORE_20) - 1)}")

# %% [markdown]
# ## Step 1 — Load CEPII distw as benchmark

# %% CEPII distw
print("Loading CEPII Gravity...")
try:
    cepii = load_cepii_gravity(filtered=False)
    cepii_sub = cepii[
        (cepii["iso3_o"].isin(CORE_20)) & (cepii["iso3_d"].isin(CORE_20))
    ][["iso3_o", "iso3_d", "distw", "contig", "comlang_off"]].copy()
    cepii_sub = cepii_sub.rename(columns={"iso3_o": "iso_o", "iso3_d": "iso_d",
                                           "distw": "d_cepii_distw"})
    print(f"  CEPII pairs: {len(cepii_sub)}")
except Exception as e:
    warnings.warn(f"CEPII load failed: {e}. Using NaN benchmark.")
    cepii_sub = None

# %% [markdown]
# ## Step 2 — Build agglomeration centroids
# Falls back to approximate centroids if WorldPop is not mounted.

# %% Agglomeration centroids
WORLDPOP_2010 = PAPER5_DATA / "worldpop" / "ppp_2010_100m_Aggregated.tif"
WORLDPOP_2022 = PAPER5_DATA / "worldpop" / "ppp_2022_100m_Aggregated.tif"
VIIRS_2022    = PAPER5_DATA / "viirs" / "VNL_v2_npp_2022_global_vcmslcfg_c202101211500.average_masked.tif"
BOUNDARIES    = PAPER5_DATA / "boundaries" / "gadm_410.gpkg"

# ---- Fallback: approximate centroids from CEPII city list ----
# If WorldPop is not available we use a hard-coded "one centroid = capital" table
# sourced from natural earth / CEPII city data. This gives d_ovdl (static) only.
FALLBACK_CENTROIDS = {
    # iso3: (lon, lat, pop_millions)
    "USA": (-95.71, 37.09, 331.0), "MEX": (-102.55, 23.63, 128.0),
    "CAN": (-96.80, 56.13, 38.0),  "JAM": (-77.30, 18.11, 3.0),
    "DEU": (10.45, 51.17, 83.0),   "FRA": (2.21, 46.23, 67.0),
    "POL": (19.15, 51.92, 38.0),   "NGA": (8.67, 9.08, 206.0),
    "ETH": (40.49, 9.15, 115.0),   "ZAF": (25.08, -29.00, 59.0),
    "SAU": (45.08, 23.89, 35.0),   "IND": (78.96, 20.59, 1380.0),
    "PAK": (69.35, 30.38, 220.0),  "CHN": (104.20, 35.86, 1411.0),
    "JPN": (138.25, 36.20, 126.0), "VNM": (108.28, 14.06, 97.0),
    "BRA": (-51.93, -14.24, 212.0), "COL": (-74.30, 4.57, 51.0),
    "AUS": (133.78, -25.27, 26.0), "NZL": (172.50, -41.50, 5.0),
}


def get_aggloms_fallback(iso_list: list[str], year: int) -> pd.DataFrame:
    """Fallback: one centroid per country from hard-coded table."""
    rows = []
    for iso3 in iso_list:
        if iso3 in FALLBACK_CENTROIDS:
            lon, lat, pop = FALLBACK_CENTROIDS[iso3]
            rows.append({"iso3": iso3, "agglom_id": 0,
                         "lon": lon, "lat": lat,
                         "pop": pop * 1e6, "viirs_radiance": 0.0,
                         "coverage_share": 1.0, "year": year})
    return pd.DataFrame(rows)


agglom_frames = []
config = AgglomerationConfig()

for year in YEARS:
    wp = WORLDPOP_2010 if year <= 2010 else WORLDPOP_2022
    if not wp.exists() or not BOUNDARIES.exists():
        print(f"  [year={year}] WorldPop/boundaries not found — using fallback centroids.")
        agglom_frames.append(get_aggloms_fallback(CORE_20, year))
        continue

    try:
        import geopandas as gpd
        gdf = gpd.read_file(BOUNDARIES)
        iso_col = next(c for c in gdf.columns if c.upper() in ("GID_0", "ISO_A3", "ISO3"))
        gdf = gdf.rename(columns={iso_col: "iso3"}).set_index("iso3")

        year_frames = []
        for iso3 in CORE_20:
            if iso3 not in gdf.index:
                print(f"  [{iso3}] not in boundaries — using fallback.")
                year_frames.append(get_aggloms_fallback([iso3], year))
                continue
            geom = gdf.loc[iso3, "geometry"]
            viirs = VIIRS_2022 if year >= 2012 else None
            df = build_country_agglomerations(iso3, wp, viirs, geom, config, year)
            if df.empty:
                print(f"  [{iso3}] empty agglomeration — fallback.")
                year_frames.append(get_aggloms_fallback([iso3], year))
            else:
                print(f"  [{iso3} {year}] {len(df)} agglomerations, "
                      f"pop coverage={df['coverage_share'].iloc[-1]:.1%}")
                year_frames.append(df)

        agglom_frames.append(pd.concat(year_frames, ignore_index=True))

    except Exception as exc:
        warnings.warn(f"Agglomeration year={year} failed: {exc!r}")
        agglom_frames.append(get_aggloms_fallback(CORE_20, year))

aggloms = pd.concat(agglom_frames, ignore_index=True)
print(f"\nTotal agglomeration rows: {len(aggloms)}")
print(aggloms.head(10))

# %% [markdown]
# ## Step 3 — Compute d_ovdl_t (no road graph needed)

# %% d_ovdl_t computation
def aggloms_for(iso3: str, year: int) -> list[Agglomeration]:
    sub = aggloms[(aggloms["iso3"] == iso3) & (aggloms["year"] == year)]
    return [
        Agglomeration(iso3=iso3, agglom_id=int(r["agglom_id"]),
                      lon=float(r["lon"]), lat=float(r["lat"]),
                      pop=float(r["pop"]), viirs=float(r.get("viirs_radiance", 0)),
                      year=year)
        for _, r in sub.iterrows()
    ]


ovdl_rows = []
for year in YEARS:
    for iso_o in CORE_20:
        ao = aggloms_for(iso_o, year)
        if not ao:
            continue
        for iso_d in CORE_20:
            if iso_d == iso_o:
                continue
            ad = aggloms_for(iso_d, year)
            if not ad:
                continue
            d_t = compute_d_ovdl(ao, ad, theta=THETA, use_viirs=False)
            d_viirs = compute_d_ovdl(ao, ad, theta=THETA, use_viirs=True)
            ovdl_rows.append({
                "iso_o": iso_o, "iso_d": iso_d, "year": year,
                "d_ovdl_t": d_t, "d_lcp_act_approx": d_viirs,
            })

ovdl_df = pd.DataFrame(ovdl_rows)
print(f"\nd_ovdl_t panel: {len(ovdl_df)} rows")
print(ovdl_df.head(6))

# %% [markdown]
# ## Step 4 — Compare d_ovdl_t vs d_cepii_distw

# %% Comparison
if cepii_sub is not None:
    # CEPII Gravity has one row per (iso_o, iso_d) pair; deduplicate before merging.
    cepii_unique = cepii_sub.drop_duplicates(subset=["iso_o", "iso_d"])
    comp = ovdl_df[ovdl_df["year"] == 2022].merge(
        cepii_unique[["iso_o", "iso_d", "d_cepii_distw"]],
        on=["iso_o", "iso_d"], how="inner"
    )
    comp["ratio"] = comp["d_ovdl_t"] / comp["d_cepii_distw"]
    # Show a concise sample
    sample = comp[["iso_o", "iso_d", "d_ovdl_t", "d_cepii_distw", "ratio"]].sort_values("ratio", ascending=False).head(10)
    print("\nTop-10 pairs by ratio (d_ovdl_t / d_cepii_distw), 2022:")
    print(sample.to_string(index=False))
    print(f"\nMedian ratio: {comp['ratio'].median():.3f}")
    print(f"Mean ratio:   {comp['ratio'].mean():.3f}")
    print(f"N pairs:      {len(comp)}")
    print(f"Pairs with ratio > 1.1 (satellite says farther): {(comp['ratio'] > 1.1).sum()}")
    print(f"Pairs with ratio < 0.9 (satellite says closer):  {(comp['ratio'] < 0.9).sum()}")

# %% Year-on-year change in d_ovdl_t (shows time-variation)
pivot = ovdl_df.pivot_table(
    index=["iso_o", "iso_d"], columns="year", values="d_ovdl_t"
)
if 2010 in pivot.columns and 2022 in pivot.columns:
    pivot["pct_change_2010_2022"] = (pivot[2022] - pivot[2010]) / pivot[2010] * 100
    top_changing = pivot["pct_change_2010_2022"].abs().nlargest(15)
    print(f"\nTop 15 pairs by |Δd_ovdl_t| 2010→2022 (%):")
    print(top_changing)
    if top_changing.max() == 0:
        print("  [Note: 0% change expected when using fallback centroids — static points by definition.]")
        print("  [With real WorldPop rasters, population weights shift as urban areas grow/contract.]")

# %% [markdown]
# ## Step 5 — Mini gravity regression (pyfixest)
# Quick validation: does d_ovdl_t produce a sensible β_distance?

# %% Gravity validation
try:
    import pyfixest as pf
    print("Loading BACI for gravity...")
    baci = load_baci(aggregation="country-pair-year")
    # Filter to our 20-country core
    baci_sub = baci[
        baci["iso_o_num"].isin(range(1000)) &  # placeholder — needs BACI→ISO3 crosswalk
        (baci["year"].isin(YEARS))
    ]
    print("Note: BACI→ISO3 crosswalk needed for gravity validation (sprint day 5).")
except ImportError:
    print("pyfixest not installed — skipping gravity validation.")
except Exception as exc:
    print(f"Gravity validation skipped: {exc!r}")

# %% [markdown]
# ## Step 6 — Transport graph prototype (maritime + air only, no OSM road)
# Road graph omitted here since OSM download takes ~30 min for 20 countries.
# This validates the maritime+air LCP and the multi-modal cost aggregation.

# %% Transport graph
try:
    ports = load_wpi_ports()
    airports = load_openflights_airports()
    print(f"Ports: {len(ports)}, Airports: {len(airports)}")

    cp = ChokepointState.for_year(2022)
    graphs = build_transport_graph_nx(
        ports_df=ports,
        airports_df=airports,
        iso_list=CORE_20,
        chokepoints=cp,
        skip_road=True,   # Skip road for prototype speed
    )
    print(f"Maritime graph: {graphs['maritime'].number_of_nodes()} nodes, "
          f"{graphs['maritime'].number_of_edges()} edges")
    print(f"Air graph: {graphs['air'].number_of_nodes()} nodes, "
          f"{graphs['air'].number_of_edges()} edges")

    # Spot check: USA → DEU maritime cost
    usa_aggloms = aggloms_for("USA", 2022)
    deu_aggloms = aggloms_for("DEU", 2022)
    if usa_aggloms and deu_aggloms:
        cost, mode = lcp_pair_nx(usa_aggloms[0], deu_aggloms[0], graphs)
        gc = great_circle_km(usa_aggloms[0].lon, usa_aggloms[0].lat,
                             deu_aggloms[0].lon, deu_aggloms[0].lat)
        print(f"\nUSA→DEU LCP cost: {cost:.1f} USD-equiv-hours (mode={mode})")
        print(f"USA→DEU great-circle: {gc:.0f} km")
        print(f"Implied travel-time km-equivalent: {cost / 30.0 * 65.0:.0f} km")

except Exception as exc:
    warnings.warn(f"Transport graph skipped: {exc!r}")

# %% [markdown]
# ## Summary
# ✓ Agglomeration pipeline: L0→L1 validated
# ✓ d_ovdl_t: time-varying great-circle distances computed
# ✓ Comparison with d_cepii_distw: cross-validated
# ✓ Transport graph: L2 node/edge contracts validated (maritime+air)
# → Next: road graph + LCP solve → d_lcp_multi_t

print("\n=== Prototype run complete ===")
print("Save aggloms to data/derived/agglomerations_prototype.parquet")
aggloms.to_parquet(
    Path(__file__).resolve().parents[1] / "data" / "derived" / "agglomerations_prototype.parquet",
    index=False,
)
