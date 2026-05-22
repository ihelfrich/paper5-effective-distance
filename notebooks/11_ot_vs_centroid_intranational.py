"""OT-vs-centroid divergence: the intra-national & close-neighbor regime.

Notebook 10 showed that for out-of-country pairs at typical international
separations (500+ km), W1 and centroid distance agree to within 1%. This
notebook tests the regime that simulation deliberately excluded:

  1. Same-country pairs (i = j) — what's the W1 of a country to itself?
     Centroid distance is zero, but W1 captures the internal-distance
     scaffolding that drives the border-effect prediction.
  2. Close neighbors with overlapping spatial support (centroids < 500 km).
  3. Atom-count sensitivity (50, 200, 500) — does resolution matter?
  4. Extreme bimodal / donut configurations.

If centroid still tracks W1 in these regimes, OT is dispensable everywhere
and the paper's contribution is the directional centroid measure, full stop.
If centroid diverges from W1 in these regimes, the OT framing earns its place
*specifically* in intra-national and border-effect contexts, which is exactly
the part of the gravity literature that's hardest to do right.

Run with the Paper 5 venv:
    .venv/bin/python notebooks/11_ot_vs_centroid_intranational.py
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

import sys
_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

from paper5.ot_distance import (
    RasterDist,
    measure_M3_pop_weighted_centroid,
    measure_M6_wasserstein_1,
    haversine_km,
)

# Reuse the shape generator from notebook 10
sys.path.insert(0, str(_repo / "notebooks"))
from importlib import import_module
nb10 = import_module("10_ot_vs_centroid_divergence")
synth_country = nb10.synth_country
ShapeSpec = nb10.ShapeSpec
SHAPES = nb10.SHAPES


# ── Experiment 1: intra-national (i = j) ────────────────────────────────────

def experiment_intranational(n_countries: int = 200, atoms: int = 100, seed: int = 2026):
    """For each synthetic country, compute W1(mu_i, mu_i) where the second
    distribution is mu_i resampled with a different seed (same shape, same
    parameters, but different stochastic realization). This is the internal
    'self-flow' distance the border-effect literature cares about.

    The centroid distance for i = i (same draw) is trivially zero. For two
    independent draws from the same distribution, the centroid distance is
    near zero (Monte-Carlo noise). The W1 is positive and meaningful.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n_countries):
        kind = SHAPES[rng.integers(0, len(SHAPES))]
        lat = rng.uniform(-50, 60); lon = rng.uniform(-150, 150)
        scale = rng.uniform(300, 1500)
        spec1 = ShapeSpec(kind=kind, n_atoms=atoms, center_lat=lat, center_lon=lon,
                          scale_km=scale, seed=int(rng.integers(0, 2**31 - 1)))
        spec2 = ShapeSpec(kind=kind, n_atoms=atoms, center_lat=lat, center_lon=lon,
                          scale_km=scale, seed=int(rng.integers(0, 2**31 - 1)))
        rd1 = synth_country(spec1)
        rd2 = synth_country(spec2)
        d_cent = measure_M3_pop_weighted_centroid(rd1, rd2)
        d_w1 = measure_M6_wasserstein_1(rd1, rd2)
        rows.append((kind, scale, d_cent, d_w1))
    return rows


# ── Experiment 2: close neighbors (overlapping spatial support) ─────────────

def experiment_close_neighbors(n_pairs: int = 250, atoms: int = 100, seed: int = 2027):
    """Pairs of countries whose centroids are between 50 and 500 km apart and
    whose spatial support overlaps. This is the border-effect / neighbor regime.
    """
    rng = np.random.default_rng(seed)
    rows = []
    while len(rows) < n_pairs:
        lat_i = rng.uniform(-50, 60); lon_i = rng.uniform(-150, 150)
        # Pick the second country center within 50-500 km in a random direction
        bearing = rng.uniform(0, 2 * math.pi)
        sep_km = rng.uniform(50, 500)
        # Approximate offset on the sphere
        d_lat_deg = (sep_km / 111.32) * math.cos(bearing)
        d_lon_deg = (sep_km / (111.32 * math.cos(math.radians(lat_i)))) * math.sin(bearing)
        lat_j = lat_i + d_lat_deg
        lon_j = lon_i + d_lon_deg
        if not (-85 < lat_j < 85):
            continue
        kind_i = SHAPES[rng.integers(0, len(SHAPES))]
        kind_j = SHAPES[rng.integers(0, len(SHAPES))]
        # Make countries of comparable scale to ensure spatial overlap
        scale = rng.uniform(200, 600)
        spec_i = ShapeSpec(kind=kind_i, n_atoms=atoms, center_lat=lat_i, center_lon=lon_i,
                           scale_km=scale, seed=int(rng.integers(0, 2**31 - 1)))
        spec_j = ShapeSpec(kind=kind_j, n_atoms=atoms, center_lat=lat_j, center_lon=lon_j,
                           scale_km=scale, seed=int(rng.integers(0, 2**31 - 1)))
        rd_i = synth_country(spec_i); rd_j = synth_country(spec_j)
        d_cent = measure_M3_pop_weighted_centroid(rd_i, rd_j)
        d_w1 = measure_M6_wasserstein_1(rd_i, rd_j)
        rows.append((kind_i, kind_j, sep_km, d_cent, d_w1))
    return rows


# ── Experiment 3: extreme-bimodal stress test ───────────────────────────────

def experiment_extreme_bimodal(n_pairs: int = 100, atoms: int = 100, seed: int = 2028):
    """Force the worst-case-for-centroid: two countries each with strong bimodal
    populations, where the modes themselves are far apart. Centroid is at the
    midpoint of the modes (essentially empty space), W1 must transport mass to
    the actual mode locations.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_pairs):
        lat_i = rng.uniform(-30, 50); lon_i = rng.uniform(-100, 100)
        lat_j = rng.uniform(-30, 50); lon_j = rng.uniform(-100, 100)
        # Force the separation to be moderate
        if haversine_km(np.array([lat_i]), np.array([lon_i]),
                        np.array([lat_j]), np.array([lon_j]))[0] < 500:
            continue
        # Extreme bimodal: lobes far from center
        s_lat = 5.0  # large scale to force strong separation between lobes
        n = atoms // 2
        def make_extreme_bimodal(cl, cn, seed_):
            r = np.random.default_rng(seed_)
            lats = np.concatenate([
                r.normal(cl, 0.5, n),
                r.normal(cl, 0.5, atoms - n),
            ])
            lons = np.concatenate([
                r.normal(cn - 8.0, 0.5, n),  # far left lobe
                r.normal(cn + 8.0, 0.5, atoms - n),  # far right lobe
            ])
            w = np.ones(atoms) / atoms
            return RasterDist(coords=np.column_stack([lats, lons]), weights=w)

        rd_i = make_extreme_bimodal(lat_i, lon_i, int(rng.integers(0, 2**31 - 1)))
        rd_j = make_extreme_bimodal(lat_j, lon_j, int(rng.integers(0, 2**31 - 1)))
        d_cent = measure_M3_pop_weighted_centroid(rd_i, rd_j)
        d_w1 = measure_M6_wasserstein_1(rd_i, rd_j)
        rows.append((d_cent, d_w1))
    return rows


# ── Experiment 4: atom-count sensitivity ────────────────────────────────────

def experiment_atom_sensitivity(seed: int = 2029):
    """Same pair, varying atom count. Does the OT-vs-centroid ratio change
    when we use higher-resolution distributions?
    """
    rng = np.random.default_rng(seed)
    atom_counts = [25, 50, 100, 200, 500, 1000]
    rows = []
    for atoms in atom_counts:
        # Repeat 20 random pairs at this resolution
        ratios = []
        for _ in range(20):
            lat_i, lon_i = rng.uniform(-50, 60), rng.uniform(-150, 150)
            lat_j, lon_j = rng.uniform(-50, 60), rng.uniform(-150, 150)
            if haversine_km(np.array([lat_i]), np.array([lon_i]),
                            np.array([lat_j]), np.array([lon_j]))[0] < 500:
                continue
            kind_i = SHAPES[rng.integers(0, len(SHAPES))]
            kind_j = SHAPES[rng.integers(0, len(SHAPES))]
            spec_i = ShapeSpec(kind=kind_i, n_atoms=atoms, center_lat=lat_i, center_lon=lon_i,
                               scale_km=rng.uniform(300, 1500),
                               seed=int(rng.integers(0, 2**31 - 1)))
            spec_j = ShapeSpec(kind=kind_j, n_atoms=atoms, center_lat=lat_j, center_lon=lon_j,
                               scale_km=rng.uniform(300, 1500),
                               seed=int(rng.integers(0, 2**31 - 1)))
            rd_i = synth_country(spec_i); rd_j = synth_country(spec_j)
            d_cent = measure_M3_pop_weighted_centroid(rd_i, rd_j)
            d_w1 = measure_M6_wasserstein_1(rd_i, rd_j)
            ratios.append(d_w1 / d_cent)
        rows.append((atoms, ratios))
    return rows


# ── Diagnostics and figure ──────────────────────────────────────────────────

def summarize(rows: list, fields: int) -> dict:
    arr = np.array([r[-2:] for r in rows], dtype=float)  # last two cols: cent, w1
    d_cent, d_w1 = arr[:, 0], arr[:, 1]
    # Guard against centroid = 0 (intra-national same-draw)
    finite_ratio_mask = d_cent > 1e-6
    if finite_ratio_mask.sum() < 2:
        return {"n": len(rows), "note": "centroid distance trivially zero or near-zero",
                "median_cent": float(np.median(d_cent)),
                "median_w1": float(np.median(d_w1)),
                "w1_minus_cent_median": float(np.median(d_w1 - d_cent))}
    ratio = d_w1[finite_ratio_mask] / d_cent[finite_ratio_mask]
    return {
        "n": len(rows),
        "pearson": float(np.corrcoef(d_cent, d_w1)[0, 1]) if d_cent.std() > 1e-9 else float("nan"),
        "spearman": float(spearmanr(d_cent, d_w1)[0]) if d_cent.std() > 1e-9 else float("nan"),
        "median_ratio": float(np.median(ratio)),
        "p10_ratio": float(np.percentile(ratio, 10)),
        "p90_ratio": float(np.percentile(ratio, 90)),
        "median_cent": float(np.median(d_cent)),
        "median_w1": float(np.median(d_w1)),
        "median_abs_diff": float(np.median(np.abs(d_w1 - d_cent))),
    }


def make_figure(intra, neighbors, extreme, atom_sens, output: Path):
    plt.rcParams.update({
        "font.family": "Helvetica Neue, sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": "#1B1B1B",
        "xtick.color": "#5A5A5A",
        "ytick.color": "#5A5A5A",
    })
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # A: intra-national — scatter cent vs w1 (cent is near-zero)
    ax = axes[0, 0]
    arr = np.array([r[-2:] for r in intra])
    ax.scatter(arr[:, 0], arr[:, 1], s=14, alpha=0.5, color="#143969", edgecolors="none")
    lim = max(arr.max(), 10)
    ax.plot([0, lim], [0, lim], color="#A52A2A", linestyle="--", linewidth=1, alpha=0.6, label="45°")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim * 1.1)
    ax.set_xlabel("Centroid distance (km)", fontsize=10)
    ax.set_ylabel("Wasserstein-1 distance (km)", fontsize=10)
    ax.set_title("A.  Intra-national (i = j, independent re-samples)\n"
                 "Centroid distance is near-zero by construction; W1 measures internal scaffolding",
                 loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    # B: close neighbors
    ax = axes[0, 1]
    arr = np.array([r[-2:] for r in neighbors])
    ax.scatter(arr[:, 0], arr[:, 1], s=14, alpha=0.5, color="#143969", edgecolors="none")
    lim = max(arr.max(), 100)
    ax.plot([0, lim], [0, lim], color="#A52A2A", linestyle="--", linewidth=1, alpha=0.6, label="45°")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim * 1.1)
    ax.set_xlabel("Centroid distance (km)", fontsize=10)
    ax.set_ylabel("Wasserstein-1 distance (km)", fontsize=10)
    ax.set_title("B.  Close neighbors (50-500 km, overlapping support)",
                 loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    # C: extreme bimodal stress test
    ax = axes[1, 0]
    arr = np.array([r[-2:] for r in extreme])
    ax.scatter(arr[:, 0], arr[:, 1], s=14, alpha=0.5, color="#143969", edgecolors="none")
    lim = max(arr.max(), 100)
    ax.plot([0, lim], [0, lim], color="#A52A2A", linestyle="--", linewidth=1, alpha=0.6, label="45°")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim * 1.1)
    ax.set_xlabel("Centroid distance (km)", fontsize=10)
    ax.set_ylabel("Wasserstein-1 distance (km)", fontsize=10)
    ax.set_title("C.  Extreme bimodal stress test\n"
                 "Both countries strongly bimodal — worst-case for centroid",
                 loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    # D: atom-count sensitivity
    ax = axes[1, 1]
    for atoms, ratios in atom_sens:
        if ratios:
            ax.scatter([atoms] * len(ratios), ratios, alpha=0.5, color="#143969", s=14, edgecolors="none")
    ax.axhline(1.0, color="#A52A2A", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Atoms per region", fontsize=10)
    ax.set_ylabel("Ratio: W1 / centroid", fontsize=10)
    ax.set_title("D.  Atom-count sensitivity\n"
                 "Does the W1-to-centroid ratio change with raster resolution?",
                 loc="left", fontsize=11, fontweight="bold")

    fig.text(0.04, 0.97,
             "OT-vs-centroid in the regime that matters: intra-national + close neighbors + extremes",
             fontsize=14, fontweight="bold", color="#1B1B1B")
    fig.subplots_adjust(top=0.92, bottom=0.06, left=0.05, right=0.97, hspace=0.30, wspace=0.20)
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    out_dir = _repo / "figures"
    out_dir.mkdir(exist_ok=True)

    print("Experiment 1: intra-national (same shape, independent re-samples)")
    intra = experiment_intranational(n_countries=120, atoms=100)
    s1 = summarize(intra, fields=4)
    print(f"  n={s1['n']}, median W1={s1['median_w1']:.1f} km, median centroid={s1['median_cent']:.1f} km")
    print(f"  median |W1 - cent|={s1['median_abs_diff']:.1f} km")

    print("\nExperiment 2: close neighbors (50-500 km, overlapping)")
    neighbors = experiment_close_neighbors(n_pairs=180, atoms=100)
    s2 = summarize(neighbors, fields=5)
    print(f"  n={s2['n']}, Pearson r={s2['pearson']:.4f}, Spearman ρ={s2['spearman']:.4f}")
    print(f"  median ratio W1/cent={s2['median_ratio']:.3f}, 10-90 pct=[{s2['p10_ratio']:.3f}, {s2['p90_ratio']:.3f}]")

    print("\nExperiment 3: extreme bimodal stress test")
    extreme = experiment_extreme_bimodal(n_pairs=80, atoms=100)
    s3 = summarize(extreme, fields=2)
    print(f"  n={s3['n']}, Pearson r={s3['pearson']:.4f}, Spearman ρ={s3['spearman']:.4f}")
    print(f"  median ratio W1/cent={s3['median_ratio']:.3f}, 10-90 pct=[{s3['p10_ratio']:.3f}, {s3['p90_ratio']:.3f}]")

    print("\nExperiment 4: atom-count sensitivity")
    atom_sens = experiment_atom_sensitivity()
    for atoms, ratios in atom_sens:
        if ratios:
            print(f"  atoms={atoms:>4}: n={len(ratios)} pairs, median ratio={np.median(ratios):.4f}, "
                  f"range=[{min(ratios):.4f}, {max(ratios):.4f}]")

    fig_path = out_dir / "fig_ot_vs_centroid_intranational.png"
    make_figure(intra, neighbors, extreme, atom_sens, fig_path)
    print(f"\nSaved figure to {fig_path}")

    print("\n=== Updated decision rule ===")
    intra_signal = s1["median_w1"] > 100  # W1 nontrivially positive when centroid ~ 0
    neighbors_diverge = abs(s2["median_ratio"] - 1.0) > 0.05
    extreme_diverge = abs(s3["median_ratio"] - 1.0) > 0.05
    if intra_signal:
        print("  Intra-national: W1 captures internal-distance scaffolding the centroid")
        print("  measure cannot. THIS IS WHERE THE OT FRAMING EARNS ITS PLACE.")
    else:
        print("  Intra-national: W1 also approaches zero. OT does not add information.")
    if neighbors_diverge or extreme_diverge:
        print("  Close neighbors / extreme bimodal: meaningful divergence detected.")
    else:
        print("  Close neighbors / extreme bimodal: centroid still tracks W1 tightly.")
    print()
    if intra_signal and not (neighbors_diverge or extreme_diverge):
        print("  ==> RECOMMENDED FRAMING: OT for INTRA-national distances and border effects.")
        print("  Centroid for INTER-national. This is a clean theoretical-empirical split.")
    elif intra_signal and (neighbors_diverge or extreme_diverge):
        print("  ==> OT genuinely captures distinct information in multiple regimes.")
        print("  Recommend full OT pipeline.")
    else:
        print("  ==> Centroid measure is sufficient everywhere. OT framing should be ")
        print("  dropped or relegated to a theoretical-motivation footnote.")


if __name__ == "__main__":
    main()
