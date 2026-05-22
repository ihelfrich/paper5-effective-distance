"""OT-vs-centroid divergence simulation.

The central theoretical claim of the proposed paper is that the
population-weighted centroid distance is a *first-moment approximation* of the
Wasserstein-1 transport cost between two activity densities. This simulation
asks the question that referee #1 will ask: how good is that approximation?

We answer empirically. Construct a corpus of synthetic "country pairs" whose
internal activity distributions span the realistic taxonomy of national
geographies (compact, bimodal, coastal-strip, archipelagic, donut). For each
pair compute both the centroid distance and the exact W1 distance. Then look
at:

  1. Scatter of W1 vs centroid distance. If the two are tightly linear, the
     centroid is a fine proxy and the OT machinery is not earning its compute.
  2. Distribution of the ratio W1 / centroid_distance. Concentrated near 1
     means good approximation; long-tailed or bimodal means the OT measure
     picks up real signal centroid misses.
  3. Spearman rank correlation. If the rank order is identical the OT measure
     adds no marginal information; if it differs the OT measure is genuinely
     a different number.
  4. The worst-offender cases — pairs where centroid and W1 disagree most —
     and what their geometry looks like. These become Figure 1 of the paper.

Run with the Paper 5 venv:
    .venv/bin/python notebooks/10_ot_vs_centroid_divergence.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

# Make the paper5 package importable when running from notebooks/
import sys
_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

from paper5.ot_distance import (
    RasterDist,
    measure_M3_pop_weighted_centroid,
    measure_M6_wasserstein_1,
    haversine_km,
)


# ── Synthetic country-shape generator ───────────────────────────────────────

@dataclass
class ShapeSpec:
    """A categorical descriptor of a synthetic country's internal geography."""
    kind: str  # 'compact' | 'bimodal' | 'coastal' | 'archipelagic' | 'donut' | 'elongated'
    n_atoms: int
    center_lat: float
    center_lon: float
    scale_km: float  # characteristic radius of the country
    seed: int


# Approximate scaling: 1 degree latitude = 111.32 km; longitude scales with cos(lat)
def km_to_deg_lat(km: float) -> float:
    return km / 111.32

def km_to_deg_lon(km: float, lat: float) -> float:
    return km / (111.32 * math.cos(math.radians(lat)))


def synth_country(spec: ShapeSpec) -> RasterDist:
    """Generate a synthetic country with given internal-density shape."""
    rng = np.random.default_rng(spec.seed)
    n = spec.n_atoms
    s_lat = km_to_deg_lat(spec.scale_km)
    s_lon = km_to_deg_lon(spec.scale_km, spec.center_lat)

    if spec.kind == "compact":
        # Single tight Gaussian — like Switzerland or Belgium
        lats = rng.normal(spec.center_lat, s_lat * 0.4, n)
        lons = rng.normal(spec.center_lon, s_lon * 0.4, n)
        w = rng.exponential(1.0, n)

    elif spec.kind == "bimodal":
        # Two distinct population centers — like a country with two megacities
        # (Brazil = Rio + São Paulo, USA = NYC + LA, Russia = Moscow + Vladivostok)
        n1 = n // 2
        n2 = n - n1
        off = s_lon * 1.2
        lats = np.concatenate([
            rng.normal(spec.center_lat, s_lat * 0.3, n1),
            rng.normal(spec.center_lat, s_lat * 0.3, n2),
        ])
        lons = np.concatenate([
            rng.normal(spec.center_lon - off, s_lon * 0.3, n1),
            rng.normal(spec.center_lon + off, s_lon * 0.3, n2),
        ])
        w = np.concatenate([rng.exponential(1.0, n1), rng.exponential(1.0, n2)])

    elif spec.kind == "coastal":
        # Activity along a strip — most population lives near one edge
        # (Chile, Australia, Peru)
        lats = rng.normal(spec.center_lat - s_lat * 0.5, s_lat * 0.2, n)
        lons = rng.uniform(spec.center_lon - s_lon, spec.center_lon + s_lon, n)
        w = rng.exponential(1.0, n)

    elif spec.kind == "archipelagic":
        # Several clusters across a region — Indonesia, Philippines, Japan
        n_clusters = 4
        clusters_per = n // n_clusters
        lats_list, lons_list = [], []
        for k in range(n_clusters):
            cl_lat = spec.center_lat + rng.uniform(-s_lat, s_lat)
            cl_lon = spec.center_lon + rng.uniform(-s_lon, s_lon)
            lats_list.append(rng.normal(cl_lat, s_lat * 0.15, clusters_per))
            lons_list.append(rng.normal(cl_lon, s_lon * 0.15, clusters_per))
        lats = np.concatenate(lats_list)
        lons = np.concatenate(lons_list)
        w = rng.exponential(1.0, lats.size)

    elif spec.kind == "donut":
        # Empty interior, ring of population — unusual but happens for
        # high-altitude or desert-centered countries (Iran, Algeria)
        theta = rng.uniform(0, 2 * math.pi, n)
        r = rng.normal(s_lat * 0.7, s_lat * 0.15, n)
        lats = spec.center_lat + r * np.sin(theta)
        lons = spec.center_lon + r * np.cos(theta) * (s_lon / s_lat)
        w = rng.exponential(1.0, n)

    elif spec.kind == "elongated":
        # Long thin country — Chile-like ratio
        lats = rng.normal(spec.center_lat, s_lat * 1.5, n)
        lons = rng.normal(spec.center_lon, s_lon * 0.2, n)
        w = rng.exponential(1.0, n)

    else:
        raise ValueError(f"Unknown shape kind: {spec.kind}")

    coords = np.column_stack([lats, lons])
    w = w / w.sum()
    return RasterDist(coords=coords, weights=w, total_mass=1.0)


# ── Experimental design ─────────────────────────────────────────────────────

SHAPES = ["compact", "bimodal", "coastal", "archipelagic", "donut", "elongated"]
N_PAIRS_PER_COMBO = 25       # how many random pairs of each (kind_i, kind_j)
ATOMS = 80                   # atoms per synthetic country (keeps OT tractable)


def run_simulation(seed: int = 0) -> dict:
    """Generate the corpus of synthetic country pairs and compute both
    measures for every pair. Returns a dict of arrays for analysis.
    """
    rng = np.random.default_rng(seed)

    rows = []  # (kind_i, kind_j, lat_i, lon_i, lat_j, lon_j, d_centroid, d_W1)

    # Realistic country centers spread across the globe.
    def random_center():
        return rng.uniform(-50, 60), rng.uniform(-150, 150)

    pair_idx = 0
    for kind_i in SHAPES:
        for kind_j in SHAPES:
            for _ in range(N_PAIRS_PER_COMBO):
                lat_i, lon_i = random_center()
                lat_j, lon_j = random_center()
                # Force pairs to be plausibly different countries
                if haversine_km(np.array([lat_i]), np.array([lon_i]),
                                np.array([lat_j]), np.array([lon_j]))[0] < 500:
                    continue
                spec_i = ShapeSpec(kind=kind_i, n_atoms=ATOMS,
                                   center_lat=lat_i, center_lon=lon_i,
                                   scale_km=rng.uniform(300, 1500),
                                   seed=int(rng.integers(0, 2**31 - 1)))
                spec_j = ShapeSpec(kind=kind_j, n_atoms=ATOMS,
                                   center_lat=lat_j, center_lon=lon_j,
                                   scale_km=rng.uniform(300, 1500),
                                   seed=int(rng.integers(0, 2**31 - 1)))
                rd_i = synth_country(spec_i)
                rd_j = synth_country(spec_j)
                d_cent = measure_M3_pop_weighted_centroid(rd_i, rd_j)
                d_w1 = measure_M6_wasserstein_1(rd_i, rd_j)
                rows.append((kind_i, kind_j, lat_i, lon_i, lat_j, lon_j,
                             d_cent, d_w1, spec_i.scale_km, spec_j.scale_km))
                pair_idx += 1

    arr = np.array([r[6:] for r in rows], dtype=float)  # (n, 4)
    kinds = [(r[0], r[1]) for r in rows]
    centers = np.array([r[2:6] for r in rows], dtype=float)

    return {
        "rows": rows,
        "kinds": kinds,
        "centers": centers,
        "d_cent": arr[:, 0],
        "d_w1": arr[:, 1],
        "scale_i": arr[:, 2],
        "scale_j": arr[:, 3],
        "n_pairs": len(rows),
    }


# ── Diagnostics ──────────────────────────────────────────────────────────────

def analyze(res: dict) -> dict:
    d_cent = res["d_cent"]
    d_w1 = res["d_w1"]
    # Pearson and Spearman correlations
    pearson = float(np.corrcoef(d_cent, d_w1)[0, 1])
    spearman, _ = spearmanr(d_cent, d_w1)
    # Ratio analysis
    ratio = d_w1 / d_cent
    # Per-shape-pair table
    per_kind = {}
    for kp in set(res["kinds"]):
        idx = [i for i, k in enumerate(res["kinds"]) if k == kp]
        per_kind[kp] = {
            "n": len(idx),
            "mean_ratio": float(np.mean(ratio[idx])),
            "median_ratio": float(np.median(ratio[idx])),
            "rank_corr": float(spearmanr(d_cent[idx], d_w1[idx])[0]) if len(idx) > 2 else float("nan"),
        }
    return {
        "pearson": pearson,
        "spearman": float(spearman),
        "ratio_mean": float(np.mean(ratio)),
        "ratio_median": float(np.median(ratio)),
        "ratio_p10": float(np.percentile(ratio, 10)),
        "ratio_p90": float(np.percentile(ratio, 90)),
        "per_kind": per_kind,
    }


# ── Figure ───────────────────────────────────────────────────────────────────

def make_figure(res: dict, summary: dict, output: Path) -> None:
    """Two-panel figure for the paper: scatter + ratio density."""
    plt.rcParams.update({
        "font.family": "Helvetica Neue, sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": "#1B1B1B",
        "xtick.color": "#5A5A5A",
        "ytick.color": "#5A5A5A",
    })
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: scatter
    ax = axes[0]
    d_cent = res["d_cent"]
    d_w1 = res["d_w1"]
    ax.scatter(d_cent, d_w1, s=8, alpha=0.32, color="#143969", edgecolors="none")
    lims = [0, max(d_cent.max(), d_w1.max()) * 1.02]
    ax.plot(lims, lims, color="#A52A2A", linewidth=1, linestyle="--", alpha=0.6,
            label="45° (centroid = W1)")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.set_xlabel("Centroid distance (km)", fontsize=11)
    ax.set_ylabel("Wasserstein-1 distance (km)", fontsize=11)
    ax.set_title(
        f"A.  Cross-pair comparison ({res['n_pairs']:,} synthetic country pairs)",
        loc="left", fontsize=12, fontweight="bold"
    )
    ax.text(0.04, 0.95,
            f"Pearson r = {summary['pearson']:.3f}\nSpearman ρ = {summary['spearman']:.3f}",
            transform=ax.transAxes, fontsize=10, color="#3A3A3A",
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#D8D8D8", linewidth=0.5))
    ax.legend(loc="lower right", frameon=False, fontsize=9)

    # Panel B: ratio density
    ax = axes[1]
    ratio = d_w1 / d_cent
    # Truncate extreme outliers for visualization
    clip = np.clip(ratio, 0, 3)
    ax.hist(clip, bins=60, color="#143969", alpha=0.65, edgecolor="white", linewidth=0.4)
    ax.axvline(1.0, color="#A52A2A", linewidth=1.5, linestyle="--", alpha=0.8,
               label="Centroid = W1 (ratio = 1)")
    ax.axvline(summary["ratio_median"], color="#1B1B1B", linewidth=1.5,
               label=f"Median = {summary['ratio_median']:.2f}")
    ax.set_xlim(0, 3)
    ax.set_xlabel("Ratio: Wasserstein-1 / centroid distance", fontsize=11)
    ax.set_ylabel("Number of country pairs", fontsize=11)
    ax.set_title("B.  Distribution of the ratio", loc="left",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", frameon=False, fontsize=9)

    # Footer
    fig.text(0.04, 0.96,
             "OT-vs-centroid divergence — when does the Wasserstein measure depart from the centroid proxy?",
             fontsize=14, fontweight="bold", color="#1B1B1B")
    fig.text(0.04, 0.92,
             "If centroid distance were a tight proxy for W1, panel A would be on the 45° line and panel B's mass would concentrate at ratio = 1.",
             fontsize=10.5, color="#555555", style="italic")
    fig.text(0.04, 0.03,
             f"Synthetic distributions over six geographic shapes (compact, bimodal, coastal, archipelagic, donut, elongated).  "
             f"Atoms per region: {ATOMS}.  Pairs per shape-combo: {N_PAIRS_PER_COMBO}.  Total pairs: {res['n_pairs']:,}.",
             fontsize=8, color="#8A8A8A", style="italic")

    fig.subplots_adjust(top=0.86, bottom=0.10, left=0.06, right=0.97, wspace=0.22)
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ── Worst-offender illustration figure ──────────────────────────────────────

def make_worst_offender_figure(res: dict, output: Path, n_show: int = 6) -> None:
    """Show the country-pair geometries where centroid and W1 disagree most.
    These become the supporting figure that justifies the OT framing.
    """
    d_cent = res["d_cent"]
    d_w1 = res["d_w1"]
    # Worst offenders: largest absolute |W1 - centroid| / centroid
    rel_diff = np.abs(d_w1 - d_cent) / d_cent
    worst_idx = np.argsort(rel_diff)[-n_show:][::-1]

    n_cols = 3
    n_rows = (n_show + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4.5 * n_rows))
    axes = np.atleast_1d(axes).flat

    for k, idx in enumerate(worst_idx):
        ax = axes[k]
        kind_i, kind_j = res["kinds"][idx]
        lat_i, lon_i, lat_j, lon_j = res["centers"][idx]
        # Regenerate the distributions for plotting (same seed via the row)
        row = res["rows"][idx]
        scale_i_km, scale_j_km = res["scale_i"][idx], res["scale_j"][idx]
        spec_i = ShapeSpec(kind=kind_i, n_atoms=ATOMS,
                           center_lat=lat_i, center_lon=lon_i,
                           scale_km=scale_i_km, seed=0)  # seed already burned; reproducibility loose here
        spec_j = ShapeSpec(kind=kind_j, n_atoms=ATOMS,
                           center_lat=lat_j, center_lon=lon_j,
                           scale_km=scale_j_km, seed=1)
        rd_i = synth_country(spec_i)
        rd_j = synth_country(spec_j)
        ax.scatter(rd_i.coords[:, 1], rd_i.coords[:, 0],
                   s=rd_i.weights * 600 * ATOMS, alpha=0.5, color="#143969",
                   edgecolors="none", label=f"Origin: {kind_i}")
        ax.scatter(rd_j.coords[:, 1], rd_j.coords[:, 0],
                   s=rd_j.weights * 600 * ATOMS, alpha=0.5, color="#A52A2A",
                   edgecolors="none", label=f"Dest: {kind_j}")
        # Centroids
        ci_lat, ci_lon = rd_i.centroid()
        cj_lat, cj_lon = rd_j.centroid()
        ax.scatter([ci_lon, cj_lon], [ci_lat, cj_lat],
                   s=150, color="black", marker="x", linewidths=2,
                   label="Centroids")
        ax.plot([ci_lon, cj_lon], [ci_lat, cj_lat],
                color="black", linewidth=1, linestyle=":")
        ax.set_title(
            f"{kind_i.title()} → {kind_j.title()}\n"
            f"centroid {d_cent[idx]:.0f} km   |   W1 {d_w1[idx]:.0f} km   "
            f"|   ratio {d_w1[idx]/d_cent[idx]:.2f}",
            fontsize=10
        )
        ax.legend(fontsize=7, loc="best", frameon=False)
        ax.set_xlabel("Longitude", fontsize=8)
        ax.set_ylabel("Latitude", fontsize=8)
        ax.tick_params(labelsize=7)

    fig.suptitle("Worst-offender pairs — where the OT measure most departs from the centroid",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    out_dir = _repo / "figures"
    out_dir.mkdir(exist_ok=True)
    derived = _repo / "data" / "derived"
    derived.mkdir(parents=True, exist_ok=True)

    print(f"Running OT-vs-centroid simulation with {len(SHAPES)**2 * N_PAIRS_PER_COMBO:,} planned pairs...")
    res = run_simulation(seed=2026)
    print(f"Completed {res['n_pairs']:,} pairs (some dropped for being too close).")

    summary = analyze(res)
    print()
    print("=== Aggregate diagnostics ===")
    print(f"  Pearson  r      = {summary['pearson']:.4f}")
    print(f"  Spearman ρ      = {summary['spearman']:.4f}")
    print(f"  Median W1/cent  = {summary['ratio_median']:.3f}")
    print(f"  Mean   W1/cent  = {summary['ratio_mean']:.3f}")
    print(f"  10-90 pct range = [{summary['ratio_p10']:.3f}, {summary['ratio_p90']:.3f}]")
    print()
    print("=== By shape-pair (showing |median ratio - 1| > 0.05) ===")
    per_kind = summary["per_kind"]
    flagged = [(kp, d) for kp, d in per_kind.items() if abs(d["median_ratio"] - 1.0) > 0.05]
    flagged.sort(key=lambda x: abs(x[1]["median_ratio"] - 1.0), reverse=True)
    for kp, d in flagged[:12]:
        print(f"  {kp[0]:>14} → {kp[1]:<14}  n={d['n']:>3}  "
              f"median ratio={d['median_ratio']:.3f}  rank corr ρ={d['rank_corr']:.3f}")

    # Save data and figures
    np.savez(derived / "ot_vs_centroid_simulation.npz",
             d_cent=res["d_cent"], d_w1=res["d_w1"])
    print(f"\nSaved data to {derived/'ot_vs_centroid_simulation.npz'}")

    fig_main = out_dir / "fig_ot_vs_centroid_divergence.png"
    make_figure(res, summary, fig_main)
    print(f"Saved main figure to {fig_main}")

    fig_worst = out_dir / "fig_worst_offenders.png"
    make_worst_offender_figure(res, fig_worst, n_show=6)
    print(f"Saved worst-offender figure to {fig_worst}")

    # ── Decision rule ────────────────────────────────────────────────────────
    print()
    print("=== Decision rule (per PRE_ANALYSIS_PLAN) ===")
    if summary["spearman"] > 0.98 and abs(summary["ratio_median"] - 1.0) < 0.05:
        print("  CENTROID IS A GOOD PROXY for W1 on synthetic data.")
        print("  Implication: OT compute may not be worth its cost on the global panel.")
        print("  Recommend: lean on the centroid measure (M5); use OT only for")
        print("  robustness / a small sample of country pairs as confirmation.")
    elif summary["spearman"] > 0.85:
        print("  CENTROID IS A MODERATE PROXY but the rank order is preserved.")
        print("  Implication: OT and centroid will give similar elasticities but")
        print("  may diverge in welfare counterfactuals. OT compute is justifiable.")
        print("  Recommend: run the full horse race; expect W1 to win in residual asymmetry.")
    else:
        print("  CENTROID IS A POOR PROXY for W1 on synthetic data.")
        print("  Implication: OT genuinely captures information centroid misses.")
        print("  Recommend: prioritize the OT measure as the paper's headline.")


if __name__ == "__main__":
    main()
