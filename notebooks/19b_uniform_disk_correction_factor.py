"""Analytical correction factor: d_eff^CES(θ) / d_eff^HM for a uniform disk.

Strategy. Rather than recomputing d_ii^eff on every real country (which has
been hitting memory/time limits on large countries), we exploit the fact
that Head-Mayer 2002 defines d_HM as the closed-form expected distance from
the CENTER of a uniform disk to a uniformly-distributed point in the same
disk (i.e., 2R/3). Equivalently d_HM = 0.67 * sqrt(A/π).

For each θ, we can compute analytically (or by simulation) the
CES generalized-mean of pairwise distances over a uniform disk, then take
the ratio:

    r(θ) = d_eff_disk(θ) / d_HM_disk

This is a single number per θ (independent of disk radius — distances scale
linearly with R, so the ratio is scale-free for symmetric measures). Apply
r(θ) to any country's d_HM to get the structurally-consistent CES distance
under the uniform-disk approximation.

This understates real cross-country variation (real countries are not disks)
but it's the right first-order correction and it lets us run the gravity
regression IMMEDIATELY without rerunning nb19.

Run:
  .venv/bin/python -u notebooks/19b_uniform_disk_correction_factor.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

from paper5.ot_distance import RasterDist, ces_effective_distance, head_mayer_closed_form

OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")


def uniform_disk_atoms(n: int, radius_km: float = 1000.0, seed: int = 42) -> RasterDist:
    """Generate n uniform-on-disk atoms in (lat, lon) at the equator."""
    rng = np.random.default_rng(seed)
    # Uniform on a disk via rejection or √r transform
    r = radius_km * np.sqrt(rng.uniform(0, 1, n))
    theta = rng.uniform(0, 2*np.pi, n)
    # Convert km to lat/lon degrees at equator (cos(lat)=1)
    coords = np.column_stack([
        r * np.cos(theta) / 111.32,  # lat
        r * np.sin(theta) / 111.32,  # lon
    ])
    w = np.ones(n) / n
    return RasterDist(coords=coords, weights=w)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Uniform-disk correction factors ===")
    print("Computing d_eff(θ) / d_HM for a unit uniform disk")
    print("Each estimate uses 5000 atoms, averaged over 5 random seeds\n")

    R = 1000.0  # km — chosen scale; ratio is scale-invariant
    n_atoms = 5000
    n_seeds = 5
    thetas = [1.0, 0.5, 0.0, -0.5, -1.0, -2.0, -3.0, -5.0, -7.0, -10.0]

    d_HM_R = 0.67 * R  # since R = sqrt(A/π), d_HM = 0.67 R

    print(f"  R = {R} km, theoretical d_HM = 0.67*R = {d_HM_R:.1f} km")
    print(f"  Theoretical arith-mean (Borel 1925): 128R/(45π) = {128*R/(45*np.pi):.1f} km")
    print()

    rows = []
    for theta in thetas:
        estimates = []
        for seed in range(n_seeds):
            rd = uniform_disk_atoms(n_atoms, radius_km=R, seed=seed)
            d_eff = ces_effective_distance(rd, rd, theta=theta)
            estimates.append(d_eff)
        d_eff_mean = float(np.mean(estimates))
        d_eff_se = float(np.std(estimates) / np.sqrt(n_seeds))
        ratio = d_eff_mean / d_HM_R
        rows.append({
            "theta": theta,
            "d_eff_km": d_eff_mean,
            "d_eff_se": d_eff_se,
            "d_HM_km": d_HM_R,
            "ratio_eff_HM": ratio,
        })
        print(f"  θ = {theta:>6g}:  d_eff = {d_eff_mean:>7.1f} km (SE {d_eff_se:>5.2f})  "
              f"ratio = {ratio:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "uniform_disk_correction_factors.csv", index=False)
    print(f"\nSaved to {OUT_DIR / 'uniform_disk_correction_factors.csv'}")

    # Key takeaway numbers
    r_neg5 = df.loc[df["theta"] == -5.0, "ratio_eff_HM"].iloc[0]
    r_neg3 = df.loc[df["theta"] == -3.0, "ratio_eff_HM"].iloc[0]
    print(f"\n=== Key correction factors ===")
    print(f"  At θ = -3 (σ = 4):  d_eff / d_HM = {r_neg3:.3f}")
    print(f"  At θ = -5 (σ = 6):  d_eff / d_HM = {r_neg5:.3f}")
    print(f"\nApplying r(θ=-5) = {r_neg5:.3f} to any country's d_HM gives the")
    print(f"first-order structurally-consistent intra-national distance.")
    print()
    print("Caveat: this is the uniform-disk approximation. Real countries with")
    print("strong density heterogeneity (USA, Russia, China) will have d_eff")
    print("further from d_HM than the disk suggests. But the direction (smaller)")
    print("and order-of-magnitude is right.")


if __name__ == "__main__":
    main()
