"""Optimal-transport distance measures for the gravity horse race.

This module implements all nine distance measures defined in PRE_ANALYSIS_PLAN
section 3, with a uniform API:

    distance(measure_id, raster_origin, raster_dest, polygon_origin, polygon_dest,
             year=None, ground_cost='gc_km', **kwargs) -> float (km, or km^theta)

The module is deliberately panel-shape-agnostic. The caller supplies two raster
arrays plus their geographic envelopes; the module does not care whether the
units are countries, US states, or arbitrary polygons. This is what lets the
same code run Liz's country-pair panel and a US-state-pair panel without
forking.

Ground-cost discipline
----------------------
Every measure ultimately depends on a *ground cost* between two raster cells
(or two centroids). The economically meaningful ground cost for trade is the
great-circle distance in kilometres, which we compute via the Haversine
formula on a spherical Earth (radius 6371 km). For iceberg-cost OT variants,
we raise this to a power $\theta$ — the iceberg-cost trade-elasticity parameter.

Numerical OT discipline
-----------------------
We use POT (Python Optimal Transport, https://pythonot.github.io) for the
Wasserstein and Sinkhorn computations. POT's `ot.emd2` gives the exact
EMD (Wasserstein-1 squared cost as a 1-D minimization), `ot.sinkhorn2` gives
the regularized Sinkhorn divergence, and `ot.sliced_wasserstein_distance`
gives the sliced approximation.

Sanity: every OT computation is validated against either (a) a closed-form
result on synthetic test distributions, or (b) the centroid distance as a
degenerate-case limit. See `tests/test_ot_distance.py`.

References
----------
Cuturi (2013): Sinkhorn distances. NeurIPS.
Peyré & Cuturi (2019): Computational Optimal Transport. F&T in ML 11.
Mayer & Zignago (2011): CEPII GeoDist database notes.
Hinz (2017): A View from Outer Space — bilateral distances from nightlights.
Li, X., et al. (2020): Harmonized DMSP/VIIRS nightlights, 1992-2018.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np

# POT is a hard dependency of this module. Imported lazily to keep import-time
# light; an ImportError here is a fatal install problem the caller should fix.
try:
    import ot  # type: ignore
    HAS_POT = True
except ImportError:
    HAS_POT = False


EARTH_RADIUS_KM = 6371.0
CES_DISTANCE_FLOOR_KM = 0.5


# ── Ground costs ─────────────────────────────────────────────────────────────

def haversine_km(lat1: np.ndarray, lon1: np.ndarray,
                 lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    """Great-circle distance in kilometres on a spherical Earth.

    Vectorized; all four arguments may be scalars or compatibly-shaped arrays.
    Returns an array of the broadcast shape.
    """
    lat1, lon1, lat2, lon2 = [np.deg2rad(a) for a in (lat1, lon1, lat2, lon2)]
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    c = 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    return EARTH_RADIUS_KM * c


def ground_cost_matrix(coords_a: np.ndarray, coords_b: np.ndarray,
                       *, theta: float = 1.0,
                       kind: Literal["gc_km", "iceberg", "squared"] = "gc_km") -> np.ndarray:
    """Pairwise ground cost between two sets of (lat, lon) points.

    Parameters
    ----------
    coords_a : (n_a, 2) array of (lat, lon) in degrees
    coords_b : (n_b, 2) array of (lat, lon) in degrees
    theta : iceberg power; ignored unless kind == 'iceberg'
    kind :
        'gc_km' : pure great-circle kilometres
        'iceberg' : (great-circle km)^theta — the economically meaningful trade-cost form
        'squared' : (great-circle km)^2 — the W2-conventional choice

    Returns
    -------
    cost : (n_a, n_b) array
    """
    lat_a = coords_a[:, 0:1]; lon_a = coords_a[:, 1:2]
    lat_b = coords_b[None, :, 0]; lon_b = coords_b[None, :, 1]
    # Broadcast: (n_a, 1) and (1, n_b) → (n_a, n_b)
    d = haversine_km(lat_a, lon_a, lat_b.squeeze(0)[None, :], lon_b.squeeze(0)[None, :])
    if kind == "gc_km":
        return d
    if kind == "iceberg":
        return d ** theta
    if kind == "squared":
        return d ** 2
    raise ValueError(f"Unknown ground-cost kind: {kind!r}")


# ── Raster → distribution conversion ─────────────────────────────────────────

@dataclass
class RasterDist:
    """A raster reinterpreted as a probability distribution on lat/lon cells.

    `weights` are normalized to sum to 1 (or to a caller-supplied total mass
    for unbalanced-OT variants). `coords` is (n, 2) of cell-centroid lat/lons.

    The conversion treats each non-zero raster cell as one atom. For typical
    1-km rasters covering a single country this is on the order of 1e4–1e6
    atoms, which is feasible for `ot.emd2` up to ~1e4 and for Sinkhorn up to
    ~1e6 with care.
    """
    coords: np.ndarray  # (n, 2) lat/lon
    weights: np.ndarray  # (n,) sums to 1 (or to total_mass)
    total_mass: float = 1.0

    @classmethod
    def from_raster(cls, raster: np.ndarray, transform, *,
                    mask: Optional[np.ndarray] = None,
                    min_weight: float = 1e-9,
                    max_atoms: Optional[int] = None) -> "RasterDist":
        """Build a RasterDist from a 2-D raster and its affine transform.

        Parameters
        ----------
        raster : (rows, cols) array of nonnegative weights
        transform : rasterio-style affine transform mapping (row, col) → (lon, lat)
        mask : optional bool array of same shape; only True cells are kept
        min_weight : drop cells below this weight (numerical floor)
        max_atoms : if set, keep only the top-N cells by weight (caps compute)
        """
        if mask is None:
            mask = np.ones_like(raster, dtype=bool)
        valid = mask & (raster > min_weight) & np.isfinite(raster)
        rows, cols = np.where(valid)
        if rows.size == 0:
            return cls(coords=np.zeros((0, 2)), weights=np.zeros(0), total_mass=0.0)

        # affine: lon = a*col + b*row + c ; lat = d*col + e*row + f
        # For rasterio-style: transform = (a, b, c, d, e, f)
        try:
            a, b, c, d, e, f = (transform.a, transform.b, transform.c,
                                transform.d, transform.e, transform.f)
        except AttributeError:
            a, b, c, d, e, f = transform[:6]
        lons = a * (cols + 0.5) + b * (rows + 0.5) + c
        lats = d * (cols + 0.5) + e * (rows + 0.5) + f
        coords = np.column_stack([lats, lons])
        weights = raster[rows, cols].astype(np.float64)
        total = weights.sum()

        if max_atoms is not None and weights.size > max_atoms:
            # Keep the top-N cells by weight — this is a coarsening that
            # preserves the high-mass tail and discards the long tail.
            top = np.argsort(weights)[-max_atoms:]
            coords = coords[top]; weights = weights[top]

        weights = weights / weights.sum()
        return cls(coords=coords, weights=weights, total_mass=float(total))


    def centroid(self) -> tuple[float, float]:
        """Mass-weighted centroid (lat, lon) in degrees.

        This is the *spherical-mean* centroid: average the unit-sphere vectors,
        re-project. For small regions it agrees with the planar centroid; for
        large or boundary-crossing regions it's the correct answer.
        """
        if self.coords.shape[0] == 0:
            return (np.nan, np.nan)
        lat = np.deg2rad(self.coords[:, 0])
        lon = np.deg2rad(self.coords[:, 1])
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        xm = float(np.sum(self.weights * x))
        ym = float(np.sum(self.weights * y))
        zm = float(np.sum(self.weights * z))
        lat_c = math.atan2(zm, math.hypot(xm, ym))
        lon_c = math.atan2(ym, xm)
        return (math.degrees(lat_c), math.degrees(lon_c))


# ── The nine measures ─────────────────────────────────────────────────────────

def measure_M1_capital(lat_lon_origin_capital: tuple[float, float],
                       lat_lon_dest_capital: tuple[float, float]) -> float:
    """M1: Capital-to-capital great-circle distance (Aitken 1973)."""
    return float(haversine_km(
        np.array([lat_lon_origin_capital[0]]),
        np.array([lat_lon_origin_capital[1]]),
        np.array([lat_lon_dest_capital[0]]),
        np.array([lat_lon_dest_capital[1]]),
    )[0])


def measure_M2_unweighted_centroid(coords_origin: np.ndarray, coords_dest: np.ndarray) -> float:
    """M2: Unweighted (geographic) centroid-to-centroid distance.

    coords_* are any sets of points belonging to the region (e.g., raster
    cell centroids or polygon vertices); we compute the spherical mean of
    each and return the great-circle distance between them.
    """
    rd_o = RasterDist(coords=coords_origin, weights=np.ones(coords_origin.shape[0]) / max(1, coords_origin.shape[0]))
    rd_d = RasterDist(coords=coords_dest, weights=np.ones(coords_dest.shape[0]) / max(1, coords_dest.shape[0]))
    lat_o, lon_o = rd_o.centroid()
    lat_d, lon_d = rd_d.centroid()
    return float(haversine_km(np.array([lat_o]), np.array([lon_o]),
                              np.array([lat_d]), np.array([lon_d]))[0])


def measure_M3_pop_weighted_centroid(pop_origin: RasterDist, pop_dest: RasterDist) -> float:
    """M3: Population-weighted centroid-to-centroid (CEPII distw_arith style)."""
    lat_o, lon_o = pop_origin.centroid()
    lat_d, lon_d = pop_dest.centroid()
    return float(haversine_km(np.array([lat_o]), np.array([lon_o]),
                              np.array([lat_d]), np.array([lon_d]))[0])


def measure_M5_directional_NL_centroid(viirs_origin: RasterDist,
                                       landscan_dest: RasterDist) -> float:
    """M5: VIIRS-origin centroid → LandScan-destination centroid.

    The asymmetric measure from the existing draft. Origin uses production
    proxy (nightlights), destination uses consumption proxy (population).
    By construction d_M5(i,j) ≠ d_M5(j,i) when VIIRS and LandScan centroids
    differ within a country (which is always, modulo edge cases).
    """
    lat_o, lon_o = viirs_origin.centroid()
    lat_d, lon_d = landscan_dest.centroid()
    return float(haversine_km(np.array([lat_o]), np.array([lon_o]),
                              np.array([lat_d]), np.array([lon_d]))[0])


def measure_M6_wasserstein_1(mu_origin: RasterDist, nu_dest: RasterDist) -> float:
    """M6: Exact Wasserstein-1 distance under great-circle ground cost.

    Compute the exact EMD between mu_i^VIIRS (origin production density) and
    nu_j^LandScan (destination consumption density). The result is in
    kilometres — the same units as M1-M5 — so coefficients are comparable.
    """
    if not HAS_POT:
        raise ImportError("POT (Python Optimal Transport) is required for M6. Install with `pip install POT`.")
    if mu_origin.coords.shape[0] == 0 or nu_dest.coords.shape[0] == 0:
        return float("nan")
    C = ground_cost_matrix(mu_origin.coords, nu_dest.coords, kind="gc_km")
    return float(ot.emd2(mu_origin.weights, nu_dest.weights, C))


def measure_M7_wasserstein_iceberg(mu_origin: RasterDist, nu_dest: RasterDist,
                                   *, p: int = 1, theta: float = 1.0) -> float:
    """M7: Wasserstein-p distance under iceberg ground cost d^theta.

    The economically meaningful OT cost for gravity. Returns a number with
    units of (km^theta)^(1/p). For p=1, theta=1 this equals M6.
    For p=1, theta in {0.5, 1, 1.5, 2} we trace the iceberg-elasticity
    sensitivity directly.
    """
    if not HAS_POT:
        raise ImportError("POT is required for M7.")
    if mu_origin.coords.shape[0] == 0 or nu_dest.coords.shape[0] == 0:
        return float("nan")
    C = ground_cost_matrix(mu_origin.coords, nu_dest.coords, kind="iceberg", theta=theta)
    if p == 1:
        return float(ot.emd2(mu_origin.weights, nu_dest.weights, C))
    # For p > 1: compute C^p, then take p-th root of EMD
    Cp = C ** p
    emd = ot.emd2(mu_origin.weights, nu_dest.weights, Cp)
    return float(emd ** (1.0 / p))


def measure_M8_sinkhorn(mu_origin: RasterDist, nu_dest: RasterDist,
                        *, eps: float = 50.0, theta: float = 1.0,
                        max_iter: int = 1000) -> float:
    """M8: Sinkhorn divergence with entropic regularization eps.

    `eps` is in the same units as the ground cost (here, km^theta). A small
    eps approaches exact OT; a large eps approaches the mean cost. Default
    eps=50 is empirically a balanced choice for country-scale rasters.

    The Sinkhorn divergence S(mu,nu) = OT_eps(mu,nu) - 0.5*OT_eps(mu,mu) -
    0.5*OT_eps(nu,nu) is computed; this de-biases the regularized cost and
    is non-negative and zero iff mu == nu.
    """
    if not HAS_POT:
        raise ImportError("POT is required for M8.")
    if mu_origin.coords.shape[0] == 0 or nu_dest.coords.shape[0] == 0:
        return float("nan")
    C_mn = ground_cost_matrix(mu_origin.coords, nu_dest.coords, kind="iceberg", theta=theta)
    C_mm = ground_cost_matrix(mu_origin.coords, mu_origin.coords, kind="iceberg", theta=theta)
    C_nn = ground_cost_matrix(nu_dest.coords, nu_dest.coords, kind="iceberg", theta=theta)
    s_mn = ot.sinkhorn2(mu_origin.weights, nu_dest.weights, C_mn, eps, numItermax=max_iter)
    s_mm = ot.sinkhorn2(mu_origin.weights, mu_origin.weights, C_mm, eps, numItermax=max_iter)
    s_nn = ot.sinkhorn2(nu_dest.weights, nu_dest.weights, C_nn, eps, numItermax=max_iter)
    return float(s_mn - 0.5 * s_mm - 0.5 * s_nn)


def ces_effective_distance(mu_origin: RasterDist, nu_dest: RasterDist,
                            *, theta: float = -4.0,
                            haversine_great_circle: bool = True,
                            d_min_km: float = CES_DISTANCE_FLOOR_KM) -> float:
    """CES-weighted effective bilateral distance — the gravity-consistent measure.

    This is the structurally-correct internal- and external-distance object
    that enters the iceberg-cost CES gravity equation. Per Head & Mayer (2002)
    and the verification analysis: the aggregate trade cost τ_ij^(1-σ) equals
    the generalized mean of pairwise location costs:

        τ_ij^(1-σ) = ∫∫ f_i(x) g_j(y) d(x,y)^θ dx dy

    where θ = δ(1-σ) for distance elasticity δ and substitution elasticity σ.
    The implied effective distance with units of km is

        d_ij^eff = ( ∫∫ f_i(x) g_j(y) d(x,y)^θ dx dy )^(1/θ).

    For typical trade-elasticity estimates σ ∈ [4, 8] and δ ≈ 1, this is
    θ ∈ [-7, -3] — concave, harmonic-dominated, short distances weighted
    heavily. CEPII's `distw_ces` uses θ = -1 as a reduced-form choice; their
    arithmetic `distw_arith` is θ = 1.

    Parameters
    ----------
    mu_origin : RasterDist
        Origin density (typically production proxy: VIIRS nightlights).
    nu_dest : RasterDist
        Destination density (typically consumption proxy: LandScan population).
        For internal distance, pass the same RasterDist twice — but the
        weights need not be symmetric; you can use VIIRS as origin and
        LandScan as destination within the same country to capture
        production-to-consumption asymmetry.
    theta : float
        Generalized-mean order. CES-consistent value is 1 - σ. Common choices:
        -1 (CEPII default), -3 (σ=4 low elasticity), -5 (σ=6 mid), -7 (σ=8 high).
        theta = 1 reproduces arithmetic mean (M3 / M5 in this module).
        theta near 0 approaches the geometric mean (logarithmic).
    haversine_great_circle : bool
        If True, ground cost is haversine great-circle km. If False, planar
        approximation (faster, OK for small countries; wrong for Russia/Canada).
    d_min_km : float
        Positive raster-resolution floor applied before CES aggregation. The
        default 0.5 km is a sub-grid regularizer for zero/near-zero cell pairs.

    Returns
    -------
    Effective distance in km.

    Notes
    -----
    Numerical safety:
      - For theta < 0, cells at distance ~0 (same cell) produce d^theta = inf.
        We keep overlapping-support pairs but treat them as sub-grid pairs by
        flooring d at d_min_km. Thus the theta -> -infinity limit is the
        minimum floored pairwise distance, max(d_min_km, min(raw d)), not the
        literal raw minimum when supports overlap.
      - For theta < 0, use the log-sum-exp trick to avoid underflow.
      - For theta = 0 exactly, return the geometric mean (limit case).
    """
    if mu_origin.coords.shape[0] == 0 or nu_dest.coords.shape[0] == 0:
        return float("nan")
    if not np.isfinite(d_min_km) or d_min_km <= 0:
        raise ValueError("d_min_km must be finite and strictly positive.")
    # Pairwise ground cost (km)
    if haversine_great_circle:
        d = ground_cost_matrix(mu_origin.coords, nu_dest.coords, kind="gc_km")
    else:
        # Approximate planar — only for small regions
        lat_a = mu_origin.coords[:, 0:1]; lon_a = mu_origin.coords[:, 1:2]
        lat_b = nu_dest.coords[None, :, 0]; lon_b = nu_dest.coords[None, :, 1]
        lat_mid = 0.5 * (lat_a + lat_b)
        d_lat_km = (lat_a - lat_b) * 111.32
        dlon = (lon_a - lon_b + 180.0) % 360.0 - 180.0
        d_lon_km = dlon * 111.32 * np.cos(np.deg2rad(lat_mid))
        d = np.sqrt(d_lat_km**2 + d_lon_km**2)

    # Clip sub-grid distances to avoid blowup for theta < 0
    d = np.maximum(d, d_min_km)

    # Weight outer product
    W = np.outer(mu_origin.weights, nu_dest.weights)
    log_d = np.log(d)

    if abs(theta) < 1e-8:
        # Geometric-mean limit: exp(sum w * log(d))
        return float(np.exp(np.sum(W * log_d)))

    # Standard generalized-mean (Hölder mean) of order theta
    # GM = (sum w_xy * d_xy^theta) ^ (1/theta)
    if theta > 0:
        moment = np.sum(W * d**theta)
        return float(moment ** (1.0 / theta))
    else:
        # theta < 0: use log-sum-exp to maintain numerical stability
        positive_weight = W > 0
        if not np.any(positive_weight):
            return float("nan")
        log_terms = np.log(W[positive_weight]) + theta * log_d[positive_weight]
        m = np.max(log_terms)
        log_moment = m + np.log(np.sum(np.exp(log_terms - m)))
        # d_eff = exp(log_moment / theta) — note theta < 0 flips sign convention
        return float(np.exp(log_moment / theta))


def head_mayer_closed_form(area_km2: float) -> float:
    """CEPII / Head-Mayer closed-form internal distance: 0.67 * sqrt(area/π).

    This is the disk-uniform approximation. Equivalent to 0.376 * sqrt(area)
    (verified: 0.67 / sqrt(π) ≈ 0.378).

    Parameters
    ----------
    area_km2 : float
        Country area in square kilometres.

    Returns
    -------
    Closed-form internal distance in km.
    """
    import math
    return 0.67 * math.sqrt(area_km2 / math.pi)


def measure_M9_sliced_wasserstein(mu_origin: RasterDist, nu_dest: RasterDist,
                                  *, n_projections: int = 100,
                                  seed: int = 42) -> float:
    """M9: Sliced Wasserstein — fast approximation for the compute-budget regime.

    Project both 2-D distributions onto n_projections random 1-D lines, compute
    1-D Wasserstein on each, and average. O(n log n) per projection, vs O(n^3)
    for exact 2-D EMD. Used to validate that the rank-ordering of country pairs
    holds under cheap approximation — if so, a global panel becomes tractable.

    Note: this uses the (lat, lon) coordinate plane directly, NOT the
    great-circle metric. For small regions the planar approximation is fine;
    for very large regions (Russia, Canada, US) we should project to a local
    equal-area projection first. TODO in v2.
    """
    if not HAS_POT:
        raise ImportError("POT is required for M9.")
    if mu_origin.coords.shape[0] == 0 or nu_dest.coords.shape[0] == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    sw = ot.sliced_wasserstein_distance(
        mu_origin.coords, nu_dest.coords,
        a=mu_origin.weights, b=nu_dest.weights,
        n_projections=n_projections, seed=int(rng.integers(0, 2**31 - 1)),
    )
    # Convert degree-units to km via a local scale at the mean latitude.
    # Great-circle 1° latitude ≈ 111.32 km; 1° longitude ≈ 111.32 km * cos(lat).
    # We use the mean of the two distribution centroids' latitudes for the
    # longitude correction. This is a coarse km-equivalent — used only for
    # cross-measure scale alignment.
    lat_o, _ = mu_origin.centroid()
    lat_d, _ = nu_dest.centroid()
    lat_mean = 0.5 * (lat_o + lat_d) if np.isfinite(lat_o + lat_d) else 0.0
    km_per_degree = 111.32 * math.cos(math.radians(lat_mean))
    return float(sw * km_per_degree)


# ── Unified dispatch ─────────────────────────────────────────────────────────

def distance(measure_id: str, *, origin: RasterDist, dest: RasterDist,
             origin_capital: Optional[tuple[float, float]] = None,
             dest_capital: Optional[tuple[float, float]] = None,
             **kwargs) -> float:
    """Unified entry point. Dispatches to one of the measure functions above."""
    mid = measure_id.upper()
    if mid == "M1":
        if origin_capital is None or dest_capital is None:
            raise ValueError("M1 requires origin_capital and dest_capital (lat, lon) tuples.")
        return measure_M1_capital(origin_capital, dest_capital)
    if mid == "M2":
        return measure_M2_unweighted_centroid(origin.coords, dest.coords)
    if mid == "M3":
        return measure_M3_pop_weighted_centroid(origin, dest)
    if mid == "M5":
        return measure_M5_directional_NL_centroid(origin, dest)
    if mid == "M6":
        return measure_M6_wasserstein_1(origin, dest)
    if mid == "M7":
        return measure_M7_wasserstein_iceberg(origin, dest,
                                              p=kwargs.get("p", 1),
                                              theta=kwargs.get("theta", 1.0))
    if mid == "M8":
        return measure_M8_sinkhorn(origin, dest,
                                   eps=kwargs.get("eps", 50.0),
                                   theta=kwargs.get("theta", 1.0))
    if mid == "M9":
        return measure_M9_sliced_wasserstein(origin, dest,
                                             n_projections=kwargs.get("n_projections", 100),
                                             seed=kwargs.get("seed", 42))
    raise ValueError(f"Unknown measure_id: {measure_id!r}. Valid: M1, M2, M3, M5, M6, M7, M8, M9.")


# ── Convenience: full horse-race for a single pair ────────────────────────────

def all_measures(origin_viirs: RasterDist, dest_landscan: RasterDist,
                 *, origin_pop: Optional[RasterDist] = None,
                 dest_pop: Optional[RasterDist] = None,
                 origin_capital: Optional[tuple[float, float]] = None,
                 dest_capital: Optional[tuple[float, float]] = None,
                 thetas: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0)) -> dict[str, float]:
    """Run every measure for one origin-dest pair and return a dict.

    Production density is VIIRS at origin; consumption density is LandScan
    at destination. For symmetric measures we use population at origin AND
    destination (passed separately as origin_pop / dest_pop).
    """
    if origin_pop is None: origin_pop = dest_landscan  # fall back; caller should provide
    if dest_pop is None: dest_pop = dest_landscan
    out: dict[str, float] = {}
    if origin_capital and dest_capital:
        out["M1_capital"] = measure_M1_capital(origin_capital, dest_capital)
    out["M2_unweighted_centroid"] = measure_M2_unweighted_centroid(
        origin_viirs.coords, dest_landscan.coords)
    out["M3_pop_centroid"] = measure_M3_pop_weighted_centroid(origin_pop, dest_pop)
    out["M5_directional_NL"] = measure_M5_directional_NL_centroid(origin_viirs, dest_landscan)
    out["M6_W1"] = measure_M6_wasserstein_1(origin_viirs, dest_landscan)
    for th in thetas:
        out[f"M7_W1_theta{th:.1f}"] = measure_M7_wasserstein_iceberg(
            origin_viirs, dest_landscan, p=1, theta=th)
    out["M8_sinkhorn"] = measure_M8_sinkhorn(origin_viirs, dest_landscan, eps=50.0, theta=1.0)
    out["M9_sliced_W"] = measure_M9_sliced_wasserstein(origin_viirs, dest_landscan)
    return out
