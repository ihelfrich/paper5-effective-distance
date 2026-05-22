"""CES-weighted effective distance vs centroid measures — the corrected test.

Background:
  Notebook 10 found that the ARITHMETIC expected pairwise distance (θ=1) is
  essentially equal to the centroid distance at long inter-country range
  (Spearman ρ = 0.9998). I concluded that OT didn't add information.

  The verification agent caught a methodology error: the gravity-consistent
  measure is NOT θ=1 (arithmetic) — it's θ = 1-σ < 0, the CES-aggregated
  harmonic-mean-like generalization. For σ ∈ [4,8], that's θ ∈ [-7,-3].

  This notebook re-runs the divergence test with the structurally-consistent
  θ values. The hypothesis is that concave costs (θ < 0) DO pick up
  sub-national density structure that the arithmetic mean misses, because
  the short-distance pair weight explodes when θ is negative.

  Key empirical question: does d_ij^eff(θ=-5) differ meaningfully from
  - the centroid-based measures (M3, M5)?
  - the CEPII baseline (θ=-1)?
  - the arithmetic mean (θ=1)?

Run:
  .venv/bin/python notebooks/12_ces_theta_divergence.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))
sys.path.insert(0, str(_repo / "notebooks"))

from paper5.ot_distance import (
    ces_effective_distance,
    measure_M3_pop_weighted_centroid,
    haversine_km,
    RasterDist,
)

from importlib import import_module
nb10 = import_module("10_ot_vs_centroid_divergence")
synth_country = nb10.synth_country
ShapeSpec = nb10.ShapeSpec
SHAPES = nb10.SHAPES

THETAS = [1.0, 0.5, -0.5, -1.0, -3.0, -5.0, -7.0]
N_PAIRS_PER_COMBO = 20
ATOMS = 80


def run_simulation(seed: int = 2026):
    rng = np.random.default_rng(seed)
    rows = []
    for kind_i in SHAPES:
        for kind_j in SHAPES:
            for _ in range(N_PAIRS_PER_COMBO):
                lat_i = rng.uniform(-50, 60); lon_i = rng.uniform(-150, 150)
                lat_j = rng.uniform(-50, 60); lon_j = rng.uniform(-150, 150)
                if haversine_km(np.array([lat_i]), np.array([lon_i]),
                                np.array([lat_j]), np.array([lon_j]))[0] < 500:
                    continue
                spec_i = ShapeSpec(kind=kind_i, n_atoms=ATOMS, center_lat=lat_i, center_lon=lon_i,
                                   scale_km=rng.uniform(300, 1500),
                                   seed=int(rng.integers(0, 2**31 - 1)))
                spec_j = ShapeSpec(kind=kind_j, n_atoms=ATOMS, center_lat=lat_j, center_lon=lon_j,
                                   scale_km=rng.uniform(300, 1500),
                                   seed=int(rng.integers(0, 2**31 - 1)))
                rd_i = synth_country(spec_i); rd_j = synth_country(spec_j)
                d_cent = measure_M3_pop_weighted_centroid(rd_i, rd_j)
                row = {"kind_i": kind_i, "kind_j": kind_j, "d_cent": d_cent}
                for theta in THETAS:
                    row[f"d_theta_{theta}"] = ces_effective_distance(rd_i, rd_j, theta=theta)
                rows.append(row)
    return rows


def intra_national_simulation(seed: int = 2027, n_countries: int = 100):
    """Intra-national (i=j) test: how does d_ii^eff behave across θ values?
    Especially important: does it differ from the Head-Mayer 0.67*sqrt(area/π)
    closed-form for bimodal / coastal / heterogeneous countries?
    """
    from paper5.ot_distance import head_mayer_closed_form
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_countries):
        kind = SHAPES[rng.integers(0, len(SHAPES))]
        lat = rng.uniform(-30, 50); lon = rng.uniform(-150, 150)
        scale_km = rng.uniform(300, 1500)
        spec = ShapeSpec(kind=kind, n_atoms=ATOMS, center_lat=lat, center_lon=lon,
                         scale_km=scale_km,
                         seed=int(rng.integers(0, 2**31 - 1)))
        rd = synth_country(spec)
        # For intra-national, use the same RasterDist for origin and dest.
        # In a real run we'd use VIIRS for origin and LandScan for dest;
        # for the synthetic test we treat them as the same distribution.
        row = {"kind": kind, "scale_km": scale_km}
        # The Head-Mayer reference: 0.67 * sqrt(area/π) where area ≈ π * scale_km^2 (disk approx)
        approx_area = math.pi * scale_km**2
        row["d_HM"] = head_mayer_closed_form(approx_area)
        for theta in THETAS:
            row[f"d_theta_{theta}"] = ces_effective_distance(rd, rd, theta=theta)
        rows.append(row)
    return rows


def analyze(rows, title):
    print(f"\n=== {title} ===")
    if not rows:
        print("  (no rows)")
        return
    # Correlations of each theta with centroid distance
    d_cent = np.array([r["d_cent"] for r in rows]) if "d_cent" in rows[0] else None
    print(f"  n = {len(rows)} pairs")
    if d_cent is not None:
        print(f"  {'theta':>7}  {'pearson':>10}  {'spearman':>10}  {'median d':>10}  {'median ratio θ/centroid':>24}")
        for theta in THETAS:
            d_theta = np.array([r[f"d_theta_{theta}"] for r in rows])
            valid = np.isfinite(d_theta) & np.isfinite(d_cent) & (d_cent > 0)
            if valid.sum() < 10:
                continue
            r_pearson = pearsonr(d_cent[valid], d_theta[valid])[0]
            r_spearman = spearmanr(d_cent[valid], d_theta[valid])[0]
            ratio = d_theta[valid] / d_cent[valid]
            print(f"  {theta:>7.1f}  {r_pearson:>10.4f}  {r_spearman:>10.4f}  "
                  f"{np.median(d_theta[valid]):>10.0f}  {np.median(ratio):>24.3f}")
    else:
        # Intra-national: compare against Head-Mayer closed form
        d_HM = np.array([r["d_HM"] for r in rows])
        print(f"  {'theta':>7}  {'median ratio θ/HM':>20}  {'pearson':>10}  "
              f"{'>1.20':>8}  {'<0.80':>8}")
        for theta in THETAS:
            d_theta = np.array([r[f"d_theta_{theta}"] for r in rows])
            valid = np.isfinite(d_theta) & np.isfinite(d_HM) & (d_HM > 0)
            if valid.sum() < 10:
                continue
            ratio = d_theta[valid] / d_HM[valid]
            r_pearson = pearsonr(d_HM[valid], d_theta[valid])[0]
            print(f"  {theta:>7.1f}  {np.median(ratio):>20.3f}  {r_pearson:>10.4f}  "
                  f"{(ratio > 1.20).sum():>8d}  {(ratio < 0.80).sum():>8d}")


def make_figure(inter_rows, intra_rows, output: Path):
    plt.rcParams.update({
        "font.family": "Helvetica Neue, sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A: inter-national d_theta vs d_centroid
    ax = axes[0]
    d_cent = np.array([r["d_cent"] for r in inter_rows])
    colors = plt.cm.viridis(np.linspace(0, 0.95, len(THETAS)))
    for theta, color in zip(THETAS, colors):
        d_theta = np.array([r[f"d_theta_{theta}"] for r in inter_rows])
        ax.scatter(d_cent, d_theta, s=6, alpha=0.5, color=color, edgecolors="none",
                   label=f"θ = {theta:+.1f}")
    lims = [0, max(d_cent.max(), 20000)]
    ax.plot(lims, lims, color="#1B1B1B", linestyle="--", linewidth=0.7, alpha=0.4)
    ax.set_xlim(lims); ax.set_ylim(0, lims[1] * 1.05); ax.set_aspect("equal")
    ax.set_xlabel("Centroid distance (km)"); ax.set_ylabel("CES effective distance (km)")
    ax.set_title("A. Inter-national pairs: d_eff vs centroid by θ", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7, frameon=False)

    # Panel B: histograms of d_theta / d_centroid by θ
    ax = axes[1]
    for theta, color in zip(THETAS, colors):
        d_theta = np.array([r[f"d_theta_{theta}"] for r in inter_rows])
        ratio = d_theta / d_cent
        ratio = ratio[np.isfinite(ratio)]
        ax.hist(np.clip(ratio, 0, 2), bins=40, alpha=0.4, color=color,
                label=f"θ = {theta:+.1f}", histtype="step", linewidth=1.5)
    ax.axvline(1.0, color="#1B1B1B", linewidth=0.7, alpha=0.6)
    ax.set_xlabel("d_eff(θ) / d_centroid"); ax.set_ylabel("count")
    ax.set_title("B. Ratio distributions by θ — inter-national", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=7, frameon=False)

    # Panel C: intra-national d_ii^eff vs Head-Mayer closed form
    ax = axes[2]
    d_HM = np.array([r["d_HM"] for r in intra_rows])
    for theta, color in zip(THETAS, colors):
        d_theta = np.array([r[f"d_theta_{theta}"] for r in intra_rows])
        ax.scatter(d_HM, d_theta, s=10, alpha=0.5, color=color, edgecolors="none",
                   label=f"θ = {theta:+.1f}")
    lim = max(d_HM.max(), max(np.nanmax(np.array([r[f"d_theta_{t}"] for r in intra_rows]))
                              for t in THETAS), 100)
    ax.plot([0, lim], [0, lim], color="#1B1B1B", linestyle="--", linewidth=0.7, alpha=0.4)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim * 1.5); ax.set_aspect("equal")
    ax.set_xlabel("Head-Mayer 0.67·√(area/π) (km)")
    ax.set_ylabel("d_ii^eff (km)")
    ax.set_title("C. Intra-national: d_ii^eff vs Head-Mayer closed-form",
                 loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7, frameon=False)

    fig.suptitle("CES-weighted effective distance vs centroid and Head-Mayer baselines, by θ",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.subplots_adjust(top=0.93, bottom=0.10, left=0.05, right=0.97, wspace=0.25)
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    out_fig = _repo / "figures" / "fig_ces_theta_divergence.png"

    print(f"Running CES-weighted inter-national simulation with θ ∈ {THETAS}...")
    inter_rows = run_simulation()
    analyze(inter_rows, f"Inter-national pairs (n={len(inter_rows)})")

    print(f"\nRunning CES-weighted intra-national simulation...")
    intra_rows = intra_national_simulation()
    analyze(intra_rows, f"Intra-national / self-distance (n={len(intra_rows)})")

    make_figure(inter_rows, intra_rows, out_fig)
    print(f"\nSaved figure to {out_fig}")

    # ── Decision rule ───────────────────────────────────────────────────────
    print("\n=== Decision rule ===")
    # Inter-national: does theta=-5 differ from centroid?
    d_cent = np.array([r["d_cent"] for r in inter_rows])
    d_th_neg5 = np.array([r["d_theta_-5.0"] for r in inter_rows])
    ratio_neg5 = d_th_neg5 / d_cent
    p10, median, p90 = np.percentile(ratio_neg5[np.isfinite(ratio_neg5)], [10, 50, 90])
    print(f"  Inter-national θ=-5 / centroid: 10pct={p10:.3f}  median={median:.3f}  90pct={p90:.3f}")

    if abs(median - 1.0) < 0.05 and p90 - p10 < 0.10:
        print("  → θ=-5 and centroid agree at inter-national distance. No new bilateral signal.")
        print("  → Border-effect angle may still work via intra-national (Panel C).")
    else:
        print("  → θ=-5 DIFFERS materially from centroid at inter-national distance.")
        print("  → The OT machinery is doing real work even in standard bilateral pairs.")

    # Intra-national: does d_ii^eff(theta=-5) differ from Head-Mayer?
    d_HM = np.array([r["d_HM"] for r in intra_rows])
    d_ii_neg5 = np.array([r["d_theta_-5.0"] for r in intra_rows])
    ratio_ii = d_ii_neg5 / d_HM
    p10, median, p90 = np.percentile(ratio_ii[np.isfinite(ratio_ii)], [10, 50, 90])
    print(f"  Intra-national θ=-5 / HM-closed: 10pct={p10:.3f}  median={median:.3f}  90pct={p90:.3f}")
    if abs(median - 1.0) < 0.05 and p90 - p10 < 0.10:
        print("  → d_ii^eff(θ=-5) ≈ Head-Mayer closed-form. Sub-national heterogeneity doesn't matter.")
        print("  → Border-effect pilot likely null.")
    else:
        print("  → d_ii^eff(θ=-5) DIFFERS materially from Head-Mayer closed-form.")
        print("  → Sub-national heterogeneity matters. Border-effect pilot has signal.")


if __name__ == "__main__":
    main()
