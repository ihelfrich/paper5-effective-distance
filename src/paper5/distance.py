"""distance.py — Multi-modal least-cost-path bilateral distance construction.

Eight distance variants:

    d_cepii_distw    CEPII Head-Mayer (2002), θ=1                 [benchmark]
    d_ovdl           Gonchar-Helfrich (2025) OVDL, great-circle   [prior]
    d_ovdl_t         OVDL annualized with yearly WorldPop+VIIRS   [prior, extended]
    d_lcp_road       LCP on OSM road graph                        [net-new]
    d_lcp_multi      LCP multi-modal (road+maritime+air)          [net-new]
    d_lcp_multi_t    Annualized multi-modal LCP                   [primary]
    d_lcp_act        Multi-modal LCP, VIIRS-activity-weighted     [net-new]
    d_net_adj        Centrality-adjusted d_lcp_multi_t            [net-new]

Pipeline:
    1. agglomerate.py → per-country agglomeration centroids (pop + VIIRS weighted)
    2. build_transport_graph() → multi-modal NetworkX graph (prototype) or
       paper5_core Rust GraphHandle (production)
    3. lcp_pair() / solve_many_to_many() → pairwise costs
    4. ces_aggregate() → country-pair-year distance scalars
    5. compute_all_variants() → distance_panel.parquet

Architecture note: The Rust solver (paper5_core) is the production path.
NetworkX is used only in the 20-country prototype (Day 3) to validate data
contracts before the Rust layer is fully wired.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Literal, Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Mode-specific freight cost rates ($/hour equivalent)
# Calibrated from Hummels & Schaur (2013): 1 day ≈ 0.8% ad valorem tariff,
# backed out to $/hour using median FOB shipment values.
# Road ≈ $30/hr, Maritime ≈ $8/hr, Air ≈ $220/hr (blended bellyhold/freighter)
# ---------------------------------------------------------------------------
FREIGHT_RATE_USD_PER_HOUR: dict[str, float] = {
    "road": 30.0,
    "maritime": 8.0,
    "air": 220.0,
    "transfer": 5.0,  # port / airport handling
}

# Average speed by mode (km/h) used to convert travel distance → time
SPEED_KM_PER_HOUR: dict[str, float] = {
    "road": 65.0,       # weighted average truck speed on mixed road classes
    "maritime": 24.0,   # ~13 knots for container ships
    "air": 850.0,       # cruise speed; door-to-door blended ~600
    "transfer": 0.0,    # handled separately as a fixed time penalty
}

PORT_HANDLING_HOURS: float = 24.0    # median port dwell time (both ends)
AIRPORT_HANDLING_HOURS: float = 4.0  # check-in + customs at both ends


# ---------------------------------------------------------------------------
# Chokepoint state (Python-side; mirrors crates/paper5-core/src/chokepoints.rs)
# ---------------------------------------------------------------------------

@dataclass
class ChokepointState:
    """State of maritime chokepoints. Used to adjust arc costs in the graph."""

    suez_open: bool = True
    suez_capacity_factor: float = 1.0    # <1.0 would impose queuing premium
    panama_capacity_factor: float = 1.0  # >1.0 during drought restriction
    red_sea_risk: float = 1.0            # Houthi insurance + rerouting premium
    global_risk: float = 1.0

    @classmethod
    def suez_ever_given_2021(cls) -> "ChokepointState":
        return cls(suez_open=False)

    @classmethod
    def panama_drought_2023(cls) -> "ChokepointState":
        return cls(panama_capacity_factor=1.8)

    @classmethod
    def houthi_red_sea_2024(cls) -> "ChokepointState":
        return cls(red_sea_risk=2.5)

    def multiplier(self, flags: int) -> float:
        """Return cost multiplier for an arc with the given chokepoint flag bitmask."""
        m = self.global_risk
        if flags & 1 and not self.suez_open:  # SUEZ closed
            m *= 1e6
        if flags & 2:  # PANAMA
            m *= self.panama_capacity_factor
        if flags & (4 | 8):  # RED_SEA or BAB_EL_MANDEB
            m *= self.red_sea_risk
        return m

    @classmethod
    def for_year(cls, year: int) -> "ChokepointState":
        """Return the canonical chokepoint state for a given year.
        Used when building annual distance panels."""
        if year == 2021:
            # Ever Given was March 23 – March 29; subnational episode.
            # Annual panel: small premium rather than full closure.
            return cls(suez_capacity_factor=1.05)
        if year == 2023:
            return cls.panama_drought_2023()
        if year == 2024:
            return cls.houthi_red_sea_2024()
        return cls()


# ---------------------------------------------------------------------------
# Agglomeration dataclass
# ---------------------------------------------------------------------------

@dataclass
class Agglomeration:
    iso3: str
    agglom_id: int
    lon: float
    lat: float
    pop: float
    viirs: float = 0.0
    year: int = 0


# ---------------------------------------------------------------------------
# Great-circle distance
# ---------------------------------------------------------------------------

def great_circle_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


def gc_cost(lon1: float, lat1: float, lon2: float, lat2: float, mode: str) -> float:
    """Great-circle cost in USD-equivalent hours."""
    d_km = great_circle_km(lon1, lat1, lon2, lat2)
    spd = SPEED_KM_PER_HOUR[mode]
    if spd == 0:
        return FREIGHT_RATE_USD_PER_HOUR[mode]
    hours = d_km / spd
    return FREIGHT_RATE_USD_PER_HOUR[mode] * hours


# ---------------------------------------------------------------------------
# CES aggregator (Head & Mayer 2002)
# ---------------------------------------------------------------------------

def ces_aggregate(
    costs: np.ndarray,      # shape (n_o, n_d)
    weights_o: np.ndarray,  # shape (n_o,)  — population or VIIRS shares
    weights_d: np.ndarray,  # shape (n_d,)
    theta: float = -1.0,
) -> float:
    """CES-aggregate pairwise costs to a scalar bilateral distance.

    d_ij = [Σ_k Σ_ℓ w_k^(i) w_ℓ^(j) c(k,ℓ)^θ]^(1/θ)

    θ = -1 corresponds to the Head-Mayer (2002) harmonic mean interpretation.
    θ = +1 is the arithmetic mean (reported as robustness).
    """
    if theta == 0:
        raise ValueError("θ=0 undefined; use θ→0 geometric mean limit.")
    w = np.outer(weights_o, weights_d)           # (n_o, n_d)
    powered = np.power(costs, theta) * w          # weighted c^θ
    agg = float(powered.sum())
    return float(agg ** (1.0 / theta))


# ---------------------------------------------------------------------------
# Transport graph builder — NetworkX prototype
# ---------------------------------------------------------------------------

def _arc_chokepoint_flags(lon1: float, lat1: float, lon2: float, lat2: float) -> int:
    """Determine which chokepoints a maritime arc passes through.

    Uses ocean-basin classification rather than midpoint bounding boxes, which
    correctly handles long arcs whose midpoint does not fall within the corridor
    (e.g., Mumbai→Novorossiysk midpoint is north of the Red Sea box but the
    actual sea route must transit Bab-el-Mandeb + Suez).

    Returns an integer bitmask:
        bit 0 = SUEZ (1)
        bit 1 = PANAMA (2)
        bit 2 = RED_SEA (4)
        bit 3 = BAB_EL_MANDEB (8)

    Ocean basin definitions
    -----------------------
    east_of_suez : Indian Ocean / Pacific / East Asian basin — lon > 43°E.
        43°E is just east of the Bab-el-Mandeb strait (~43.3°E).
    west_of_suez : Mediterranean, Black Sea, Baltic, North Sea, and Atlantic
        basin accessible via Suez (lon ∈ [−80, 42], lat > −10).
        Excludes: US/Canadian Pacific coast and Pacific Americas (lon < −80)
        which access Asia via trans-Pacific routes, not via Suez.

    An arc from east_of_suez to west_of_suez (or vice versa) MUST transit
    the Red Sea / Bab-el-Mandeb. Dijkstra will naturally route around this
    via sub-Saharan waypoints (Cape Town, Durban — neither is east nor west)
    if the chokepoint multiplier makes the Red Sea leg too expensive.
    """
    mid_lon = (lon1 + lon2) / 2.0
    mid_lat = (lat1 + lat2) / 2.0
    flags = 0

    # ── Panama Canal ──────────────────────────────────────────────────────────
    # One endpoint in deep Pacific (lon < −90) and the other in the Atlantic /
    # Caribbean / Europe (lon > −75).
    if (lon1 < -90.0 and lon2 > -75.0) or (lon2 < -90.0 and lon1 > -75.0):
        flags |= 2  # PANAMA

    # ── Suez Canal / Red Sea ──────────────────────────────────────────────────
    # "east of Suez": Indian Ocean, Bay of Bengal, Pacific, East Asia
    east1 = lon1 > 43.0
    east2 = lon2 > 43.0

    # "west of Suez": Mediterranean, Black Sea, Baltic/North Sea, Atlantic
    # (everything from Gibraltar to the US East Coast, but NOT Pacific Americas)
    west1 = (-80.0 <= lon1 <= 42.0) and (lat1 > -10.0)
    west2 = (-80.0 <= lon2 <= 42.0) and (lat2 > -10.0)

    if (east1 and west2) or (east2 and west1):
        # Arc crosses from Indian Ocean basin into Med/Atlantic or vice versa.
        # Sea route MUST pass through Bab-el-Mandeb → Red Sea → Suez
        # (or go around the Cape of Good Hope, which Dijkstra can find via
        # unflagged waypoints: Cape Town / Durban are neither east nor west).
        flags |= (1 | 4)  # SUEZ + RED_SEA

    # ── Local Red Sea / Suez box (arcs wholly within the corridor) ────────────
    if 10.0 <= mid_lat <= 32.0 and 32.0 <= mid_lon <= 45.0:
        flags |= 4  # RED_SEA
    if 29.0 <= mid_lat <= 32.0 and 31.0 <= mid_lon <= 33.0:
        flags |= (1 | 4)  # SUEZ + RED_SEA

    return flags


def _build_maritime_graph_nx(
    ports_df: pd.DataFrame,
    chokepoints: ChokepointState,
) -> "nx.Graph":
    """Build a fully-connected port-to-port maritime graph (great-circle costs).

    All 102 prototype ports are connected to every other port (full graph).
    Edge costs include chokepoint multipliers applied at build time, so the
    graph embeds the current chokepoint state — rebuild to change state.
    """
    import networkx as nx

    G = nx.Graph()
    for idx, row in ports_df.iterrows():
        G.add_node(int(idx), lon=float(row["longitude"]),
                   lat=float(row["latitude"]), country=str(row["country"]))

    port_ids = list(ports_df.index)
    lons = ports_df["longitude"].values
    lats = ports_df["latitude"].values

    for i, pi in enumerate(port_ids):
        for j, pj in enumerate(port_ids[i + 1:], start=i + 1):
            lon1, lat1 = float(lons[i]), float(lats[i])
            lon2, lat2 = float(lons[j]), float(lats[j])
            d = great_circle_km(lon1, lat1, lon2, lat2)
            base_cost = d / SPEED_KM_PER_HOUR["maritime"] * FREIGHT_RATE_USD_PER_HOUR["maritime"]

            # Identify chokepoint flags and apply multiplier
            flags = _arc_chokepoint_flags(lon1, lat1, lon2, lat2)
            cp_mult = chokepoints.multiplier(flags)

            # If Suez is closed, arcs that need Suez get a prohibitive penalty
            # (the Cape route will naturally be cheaper via alternative paths)
            if not chokepoints.suez_open and (flags & 1):
                cp_mult = 1e6

            cost = base_cost * cp_mult * chokepoints.global_risk
            G.add_edge(pi, pj, weight=cost, length_km=d, mode="maritime", flags=flags)

    return G


def _build_road_graph_nx(
    iso_list: list[str],
    use_cache: bool = True,
    cache_dir: Optional[Path] = None,
) -> "nx.Graph":
    """Download OSM road networks for `iso_list` countries and stitch into one graph.

    Uses osmnx with the 'drive' network type. For countries with large road networks
    this can be slow; use cache_dir to persist downloaded graphs.

    Returns a weighted DiGraph; edge weights are travel time in hours × road_rate.
    """
    import networkx as nx
    try:
        import osmnx as ox
    except ImportError as e:
        raise ImportError("Road graph requires osmnx: pip install osmnx") from e

    combined = nx.DiGraph()

    for iso3 in iso_list:
        cache_path = None
        if cache_dir is not None:
            cache_path = Path(cache_dir) / f"road_{iso3}.graphml"
            if use_cache and cache_path.exists():
                G = ox.load_graphml(cache_path)
                combined = nx.compose(combined, nx.DiGraph(G))
                continue

        try:
            # osmnx geocodes country name from ISO3 via Nominatim
            G = ox.graph_from_place(iso3, network_type="drive", simplify=True)
            if cache_path is not None:
                ox.save_graphml(G, filepath=cache_path)
        except Exception as exc:
            warnings.warn(f"[{iso3}] OSM road download failed: {exc!r}")
            continue

        # Convert edge travel_time (seconds) → cost
        G_proj = ox.add_edge_speeds(G)
        G_proj = ox.add_edge_travel_times(G_proj)
        for u, v, data in G_proj.edges(data=True):
            t_hours = data.get("travel_time", 0) / 3600.0
            data["weight"] = t_hours * FREIGHT_RATE_USD_PER_HOUR["road"]
            data["mode"] = "road"

        combined = nx.compose(combined, G_proj)

    return combined


def _build_air_graph_nx(
    airports_df: pd.DataFrame,
    routes_df: Optional[pd.DataFrame] = None,
    max_airports: int = 250,
) -> "nx.Graph":
    """Build airport-to-airport graph from OpenFlights data.

    airports_df: from data_loaders.load_openflights_airports().
    routes_df: Optional routes table; if None, uses great-circle between all
               airports within 15,000 km as proxy for served routes.
    max_airports: cap for the O(n²) fallback; keeps build time tractable.
                  Major international hubs only. Default 250 covers all
                  significant freight/passenger hubs globally.

    Edge weight = air cost (USD-equivalent hours) including airport handling.
    """
    import networkx as nx

    # Filter to real airports with IATA codes (excludes heliports, small airstrips)
    df = airports_df.copy()
    df = df[
        df["type"].fillna("airport").str.lower().isin(["airport"]) &
        df["iata"].notna() &
        ~df["iata"].isin([r"\N", "", "\\N"]) &
        df["lat"].notna() & df["lon"].notna()
    ]
    # If we have more than max_airports, keep a geographically spread subset.
    # Heuristic: take every Nth airport sorted by country+iata so we get
    # global coverage, not just one region.
    if len(df) > max_airports:
        df = df.sort_values(["country", "iata"]).iloc[::max(1, len(df)//max_airports)].head(max_airports)

    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_node(int(row["airport_id"]),
                   lon=float(row["lon"]), lat=float(row["lat"]),
                   iata=str(row["iata"]), country=str(row.get("country", "")))

    if routes_df is not None:
        for _, route in routes_df.iterrows():
            src, dst = route.get("src_airport_id"), route.get("dst_airport_id")
            if src not in G.nodes or dst not in G.nodes:
                continue
            lon1, lat1 = G.nodes[src]["lon"], G.nodes[src]["lat"]
            lon2, lat2 = G.nodes[dst]["lon"], G.nodes[dst]["lat"]
            d = great_circle_km(lon1, lat1, lon2, lat2)
            cost = (d / SPEED_KM_PER_HOUR["air"] + AIRPORT_HANDLING_HOURS) * FREIGHT_RATE_USD_PER_HOUR["air"]
            G.add_edge(src, dst, weight=cost, length_km=d, mode="air")
    else:
        # Great-circle fallback: all pairs under 15,000 km.
        # With max_airports=250 this is 250²/2 = 31K iterations — fast.
        node_list = list(G.nodes)
        for i, a in enumerate(node_list):
            for b in node_list[i + 1:]:
                lon1, lat1 = G.nodes[a]["lon"], G.nodes[a]["lat"]
                lon2, lat2 = G.nodes[b]["lon"], G.nodes[b]["lat"]
                d = great_circle_km(lon1, lat1, lon2, lat2)
                if d > 15_000:
                    continue
                cost = (d / SPEED_KM_PER_HOUR["air"] + AIRPORT_HANDLING_HOURS) * FREIGHT_RATE_USD_PER_HOUR["air"]
                G.add_edge(a, b, weight=cost, length_km=d, mode="air")

    return G


def build_transport_graph_nx(
    ports_df: pd.DataFrame,
    airports_df: pd.DataFrame,
    iso_list: list[str],
    chokepoints: ChokepointState,
    road_cache_dir: Optional[Path] = None,
    skip_road: bool = False,
) -> dict[str, "nx.Graph"]:
    """Build the three-mode transport graph for a given chokepoint state.

    Returns dict with keys "road", "maritime", "air", each a NetworkX graph.
    Node coordinates are lon/lat (WGS84). Edge weights are USD-equivalent hours.

    This is the Python/NetworkX prototype. Once paper5_core is compiled,
    `build_transport_graph_rust()` replaces this for production.
    """
    graphs: dict[str, "nx.Graph"] = {}

    if not skip_road:
        graphs["road"] = _build_road_graph_nx(iso_list, cache_dir=road_cache_dir)

    graphs["maritime"] = _build_maritime_graph_nx(ports_df, chokepoints)
    graphs["air"] = _build_air_graph_nx(airports_df)
    return graphs


# ---------------------------------------------------------------------------
# LCP solve — NetworkX prototype
# ---------------------------------------------------------------------------

def _nearest_node_nx(G: "nx.Graph", lon: float, lat: float) -> int:
    """Find nearest graph node to (lon, lat) by Euclidean proximity in degrees.
    Production version uses spatial index; prototype uses brute-force.
    """
    import networkx as nx
    best_id, best_d = None, float("inf")
    for n, data in G.nodes(data=True):
        n_lon = data.get("x", data.get("lon", 0))
        n_lat = data.get("y", data.get("lat", 0))
        d = (n_lon - lon) ** 2 + (n_lat - lat) ** 2
        if d < best_d:
            best_d = d
            best_id = n
    return best_id


def lcp_pair_nx(
    origin: Agglomeration,
    dest: Agglomeration,
    graphs: dict[str, "nx.Graph"],
    mode: Optional[Literal["road", "maritime", "air"]] = None,
) -> tuple[float, str]:
    """Compute LCP cost (USD-equivalent hours) between two agglomerations.

    If mode is None, computes the minimum over all available modes.
    Returns (cost, winning_mode).

    This is the NetworkX prototype; identical interface to the Rust path.
    """
    import networkx as nx

    modes_to_try = [mode] if mode is not None else list(graphs.keys())
    best_cost, best_mode = float("inf"), "none"

    for m in modes_to_try:
        G = graphs.get(m)
        if G is None:
            continue
        try:
            src_node = _nearest_node_nx(G, origin.lon, origin.lat)
            dst_node = _nearest_node_nx(G, dest.lon, dest.lat)
            cost = nx.shortest_path_length(G, src_node, dst_node, weight="weight")
            # Add inter-modal transfer overhead for maritime / air
            if m == "maritime":
                cost += PORT_HANDLING_HOURS * FREIGHT_RATE_USD_PER_HOUR["maritime"]
            elif m == "air":
                cost += AIRPORT_HANDLING_HOURS * FREIGHT_RATE_USD_PER_HOUR["air"]
            if cost < best_cost:
                best_cost = cost
                best_mode = m
        except nx.NetworkXNoPath:
            pass
        except Exception as exc:
            warnings.warn(f"LCP {m} {origin.iso3}→{dest.iso3}: {exc!r}")

    if best_cost == float("inf"):
        # Fall back to great-circle with road rate
        gc = great_circle_km(origin.lon, origin.lat, dest.lon, dest.lat)
        best_cost = gc / SPEED_KM_PER_HOUR["road"] * FREIGHT_RATE_USD_PER_HOUR["road"]
        best_mode = "gc_fallback"

    return best_cost, best_mode


# ---------------------------------------------------------------------------
# Country-pair aggregation
# ---------------------------------------------------------------------------

def compute_bilateral_distance(
    aggloms_o: list[Agglomeration],
    aggloms_d: list[Agglomeration],
    cost_fn: Callable[[Agglomeration, Agglomeration], float],
    theta: float = -1.0,
    weight_attr: str = "pop",   # "pop" or "viirs"
) -> float:
    """CES-aggregate pairwise costs to a scalar bilateral distance.

    cost_fn: callable(origin_agglom, dest_agglom) → float cost
    weight_attr: "pop" for standard variant, "viirs" for activity-weighted.
    """
    n_o, n_d = len(aggloms_o), len(aggloms_d)
    costs = np.zeros((n_o, n_d))
    for i, ao in enumerate(aggloms_o):
        for j, ad in enumerate(aggloms_d):
            costs[i, j] = cost_fn(ao, ad)

    # Normalise weights
    pop_o = np.array([getattr(a, weight_attr) for a in aggloms_o], dtype=float)
    pop_d = np.array([getattr(a, weight_attr) for a in aggloms_d], dtype=float)
    sum_o, sum_d = pop_o.sum(), pop_d.sum()
    if sum_o == 0:
        pop_o = np.ones(n_o) / n_o
    else:
        pop_o /= sum_o
    if sum_d == 0:
        pop_d = np.ones(n_d) / n_d
    else:
        pop_d /= sum_d

    return ces_aggregate(costs, pop_o, pop_d, theta=theta)


# ---------------------------------------------------------------------------
# OVDL (Gonchar-Helfrich 2025) — great-circle with VIIRS/pop weights
# ---------------------------------------------------------------------------

def compute_d_ovdl(
    aggloms_o: list[Agglomeration],
    aggloms_d: list[Agglomeration],
    theta: float = -1.0,
    use_viirs: bool = False,
) -> float:
    """OVDL distance: great-circle, pop-weighted (or VIIRS-weighted)."""
    weight_attr = "viirs" if use_viirs else "pop"

    def gc_cost_fn(ao: Agglomeration, ad: Agglomeration) -> float:
        return great_circle_km(ao.lon, ao.lat, ad.lon, ad.lat)

    return compute_bilateral_distance(aggloms_o, aggloms_d, gc_cost_fn, theta, weight_attr)


# ---------------------------------------------------------------------------
# Main computation functions
# ---------------------------------------------------------------------------

def compute_d_lcp_multi_t(
    agglomerations: pd.DataFrame,
    graphs_by_year: dict[int, dict],
    theta: float = -1.0,
    weight_attr: str = "pop",
) -> pd.DataFrame:
    """Compute d_lcp_multi_t for all country pairs × years.

    Parameters
    ----------
    agglomerations : DataFrame from build_all_agglomerations(), stacked across years.
                     Must have columns: iso3, agglom_id, lon, lat, pop, viirs_radiance, year.
    graphs_by_year : dict mapping year → {"road": G, "maritime": G, "air": G}.
    theta : CES parameter (default -1.0 = Head-Mayer harmonic).
    weight_attr : "pop" or "viirs" (for activity-weighted variant).

    Returns
    -------
    DataFrame: iso_o, iso_d, year, d_lcp_multi_t.
    """
    rows = []
    isos = sorted(agglomerations["iso3"].unique())

    for year, graphs in graphs_by_year.items():
        year_df = agglomerations[agglomerations["year"] == year]
        # Build agglomeration objects per country
        country_aggloms: dict[str, list[Agglomeration]] = {}
        for iso3, grp in year_df.groupby("iso3"):
            country_aggloms[iso3] = [
                Agglomeration(
                    iso3=iso3,
                    agglom_id=int(r["agglom_id"]),
                    lon=float(r["lon"]),
                    lat=float(r["lat"]),
                    pop=float(r["pop"]),
                    viirs=float(r.get("viirs_radiance", 0)),
                    year=year,
                )
                for _, r in grp.iterrows()
            ]

        for iso_o in isos:
            if iso_o not in country_aggloms:
                continue
            aggloms_o = country_aggloms[iso_o]
            for iso_d in isos:
                if iso_d == iso_o or iso_d not in country_aggloms:
                    continue
                aggloms_d = country_aggloms[iso_d]

                def cost_fn(ao: Agglomeration, ad: Agglomeration) -> float:
                    c, _ = lcp_pair_nx(ao, ad, graphs)
                    return c

                d = compute_bilateral_distance(aggloms_o, aggloms_d, cost_fn, theta, weight_attr)
                rows.append({"iso_o": iso_o, "iso_d": iso_d, "year": year, "d_lcp_multi_t": d})

    return pd.DataFrame(rows)


def compute_all_variants(
    agglomerations: pd.DataFrame,
    graphs_by_year: dict[int, dict],
    cepii_df: Optional[pd.DataFrame] = None,
    centrality_panel: Optional[pd.DataFrame] = None,
    theta: float = -1.0,
) -> pd.DataFrame:
    """Compute all eight distance variants and return a wide panel.

    Columns: iso_o, iso_d, year, d_cepii_distw, d_ovdl, d_ovdl_t,
             d_lcp_road, d_lcp_multi, d_lcp_multi_t, d_lcp_act, d_net_adj.
    """
    # --- d_lcp_multi_t (primary) ---
    panel = compute_d_lcp_multi_t(agglomerations, graphs_by_year, theta=theta, weight_attr="pop")

    # --- d_ovdl_t (VIIRS great-circle, time-varying) ---
    ovdl_rows = []
    isos = sorted(agglomerations["iso3"].unique())
    for year in agglomerations["year"].unique():
        year_df = agglomerations[agglomerations["year"] == year]
        country_aggloms: dict[str, list[Agglomeration]] = {}
        for iso3, grp in year_df.groupby("iso3"):
            country_aggloms[iso3] = [
                Agglomeration(iso3=iso3, agglom_id=int(r["agglom_id"]),
                              lon=float(r["lon"]), lat=float(r["lat"]),
                              pop=float(r["pop"]), viirs=float(r.get("viirs_radiance", 0)),
                              year=year)
                for _, r in grp.iterrows()
            ]
        for iso_o in isos:
            if iso_o not in country_aggloms:
                continue
            for iso_d in isos:
                if iso_d == iso_o or iso_d not in country_aggloms:
                    continue
                d_ovdl = compute_d_ovdl(country_aggloms[iso_o], country_aggloms[iso_d],
                                        theta=theta, use_viirs=False)
                d_ovdl_viirs = compute_d_ovdl(country_aggloms[iso_o], country_aggloms[iso_d],
                                              theta=theta, use_viirs=True)
                ovdl_rows.append({"iso_o": iso_o, "iso_d": iso_d, "year": year,
                                  "d_ovdl_t": d_ovdl, "d_lcp_act_approx": d_ovdl_viirs})

    ovdl_df = pd.DataFrame(ovdl_rows)
    panel = panel.merge(ovdl_df, on=["iso_o", "iso_d", "year"], how="outer")

    # --- d_cepii_distw (static; merge from cepii_df if supplied) ---
    if cepii_df is not None:
        distw = (
            cepii_df[["iso3_o", "iso3_d", "distw"]]
            .rename(columns={"iso3_o": "iso_o", "iso3_d": "iso_d", "distw": "d_cepii_distw"})
            .drop_duplicates(subset=["iso_o", "iso_d"])
        )
        panel = panel.merge(distw, on=["iso_o", "iso_d"], how="left")
    else:
        panel["d_cepii_distw"] = np.nan

    # --- d_net_adj (centrality-adjusted d_lcp_multi_t) ---
    if centrality_panel is not None:
        # cent: iso3, year, centrality_norm (normalized to mean=1)
        cp = centrality_panel.rename(columns={"iso3": "iso_o"})
        panel = panel.merge(cp[["iso_o", "year", "centrality_norm"]],
                            on=["iso_o", "year"], how="left")
        panel["d_net_adj"] = panel["d_lcp_multi_t"] / panel["centrality_norm"].fillna(1.0)
        panel = panel.drop(columns=["centrality_norm"])
    else:
        panel["d_net_adj"] = np.nan

    # --- Fill remaining variant columns with NaN (road-only; full multi) ---
    for col in ("d_ovdl", "d_lcp_road", "d_lcp_multi", "d_lcp_act"):
        if col not in panel.columns:
            panel[col] = np.nan

    # d_ovdl: time-invariant OVDL from first year's agglomerations (proxy for static)
    first_year = int(agglomerations["year"].min())
    fy_df = agglomerations[agglomerations["year"] == first_year]
    country_aggloms_fy: dict[str, list[Agglomeration]] = {}
    for iso3, grp in fy_df.groupby("iso3"):
        country_aggloms_fy[iso3] = [
            Agglomeration(iso3=iso3, agglom_id=int(r["agglom_id"]),
                          lon=float(r["lon"]), lat=float(r["lat"]),
                          pop=float(r["pop"]), viirs=float(r.get("viirs_radiance", 0)),
                          year=first_year)
            for _, r in grp.iterrows()
        ]
    ovdl_static: dict[tuple[str, str], float] = {}
    for iso_o in isos:
        if iso_o not in country_aggloms_fy:
            continue
        for iso_d in isos:
            if iso_d == iso_o or iso_d not in country_aggloms_fy:
                continue
            ovdl_static[(iso_o, iso_d)] = compute_d_ovdl(
                country_aggloms_fy[iso_o], country_aggloms_fy[iso_d], theta=theta
            )
    panel["d_ovdl"] = panel.apply(
        lambda r: ovdl_static.get((r["iso_o"], r["iso_d"]), np.nan), axis=1
    )

    col_order = [
        "iso_o", "iso_d", "year",
        "d_cepii_distw", "d_ovdl", "d_ovdl_t",
        "d_lcp_road", "d_lcp_multi", "d_lcp_multi_t",
        "d_lcp_act", "d_net_adj",
    ]
    return panel[[c for c in col_order if c in panel.columns]]
