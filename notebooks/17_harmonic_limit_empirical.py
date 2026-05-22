"""Empirical illustration of the harmonic-limit theorem.

For 702 ordered country pairs we have d_eff(θ=-5) from nb13 v3, plus CEPII's
distw_harmonic and distw_arithmetic. We show two regimes:

  Regime 1 (separated supports): d_eff(θ=-5) ≈ 0.98 × CEPII for non-contiguous pairs.
  Regime 2 (touching supports): d_eff(θ=-5) collapses to ~20-30 km for contiguous pairs.

This produces the central figure for the harmonic-limit theoretical note.

The "contiguous" classifier comes from CEPII's contig variable.

Run:
  .venv/bin/python notebooks/17_harmonic_limit_empirical.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))


def main():
    eff_cepii_csv = _repo / "data" / "derived" / "real_country_eff_vs_cepii.csv"
    if not eff_cepii_csv.exists():
        print(f"Need {eff_cepii_csv}")
        return

    df = pd.read_csv(eff_cepii_csv)
    print(f"Loaded {len(df)} pairs from nb13 v3 bilateral × CEPII merge")
    print(df.columns.tolist())

    # Join in contig from CEPII Gravity 2010
    cepii_path = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")
    cepii = pd.read_csv(
        cepii_path, dtype={"iso3_o": str, "iso3_d": str},
        usecols=["year", "iso3_o", "iso3_d", "contig"],
    )
    cepii = cepii[cepii["year"] == 2010].drop_duplicates(["iso3_o", "iso3_d"])
    cepii = cepii.rename(columns={"iso3_o": "iso_o", "iso3_d": "iso_d"})
    merged = df.merge(cepii[["iso_o", "iso_d", "contig"]], on=["iso_o", "iso_d"], how="left")
    print(f"  merged with contig: {merged['contig'].notna().sum()} matched")

    merged["regime"] = np.where(merged["contig"] == 1, "touching", "separated")

    print("\n=== Two regimes of CES effective distance at θ = -5 ===")
    for regime in ("separated", "touching"):
        sub = merged[merged["regime"] == regime].dropna(subset=["d_eff", "distw_harmonic"])
        ratio = sub["d_eff"] / sub["distw_harmonic"]
        print(f"\n  {regime.upper()} (n = {len(sub)}):")
        print(f"    d_eff(θ=-5)        : median {sub['d_eff'].median():>7.0f} km, "
              f"mean {sub['d_eff'].mean():>7.0f} km")
        print(f"    CEPII distw_harm   : median {sub['distw_harmonic'].median():>7.0f} km, "
              f"mean {sub['distw_harmonic'].mean():>7.0f} km")
        print(f"    ratio d_eff/CEPII  : median {ratio.median():.4f}, "
              f"10th pct {ratio.quantile(0.10):.4f}, "
              f"90th pct {ratio.quantile(0.90):.4f}")

    # Save the regime-tagged dataset
    out = _repo / "data" / "derived" / "harmonic_limit_two_regimes.csv"
    merged.to_csv(out, index=False)
    print(f"\nSaved to {out}")

    # Make the central figure: log d_eff vs log CEPII, colored by regime
    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        mpl.rcParams["font.family"] = "Charter"
        mpl.rcParams["font.size"] = 11
        mpl.rcParams["axes.spines.top"] = False
        mpl.rcParams["axes.spines.right"] = False

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        for regime, color in [("separated", "#4C7066"), ("touching", "#B5482A")]:
            sub = merged[merged["regime"] == regime].dropna(subset=["d_eff", "distw_harmonic"])
            ax1.scatter(sub["distw_harmonic"], sub["d_eff"],
                       c=color, alpha=0.5, s=14, label=f"{regime} (n={len(sub)})")
            ratio = sub["d_eff"] / sub["distw_harmonic"]
            ax2.hist(np.log10(ratio.replace([np.inf, -np.inf], np.nan).dropna()),
                     bins=40, alpha=0.55, color=color, label=f"{regime}")

        x = np.logspace(1, 4.5, 100)
        ax1.plot(x, x, "k--", linewidth=0.8, alpha=0.4, label="45° line")
        ax1.set_xscale("log"); ax1.set_yscale("log")
        ax1.set_xlabel("CEPII distw_harmonic (km)")
        ax1.set_ylabel(r"$d^{\mathrm{eff}}_{ij}(\theta = -5)$ (km)")
        ax1.set_title("Two regimes of CES effective distance")
        ax1.legend(loc="lower right", frameon=False)
        ax1.grid(True, alpha=0.2)

        ax2.set_xlabel(r"$\log_{10}(d^{\mathrm{eff}}/d^{\mathrm{CEPII}})$")
        ax2.set_ylabel("Density")
        ax2.set_title("Distribution of distance ratio by regime")
        ax2.legend(loc="upper left", frameon=False)
        ax2.axvline(0, color="black", linewidth=0.8, alpha=0.3)
        ax2.grid(True, alpha=0.2)

        fig.suptitle("CES effective distance vs. CEPII harmonic distance, "
                    "by border contiguity",
                    fontsize=13)
        fig.tight_layout()

        outfig = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/figures/harmonic_limit_two_regimes.pdf")
        outfig.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outfig, bbox_inches="tight", dpi=150)
        fig.savefig(outfig.with_suffix(".png"), bbox_inches="tight", dpi=200)
        print(f"\nFigure saved to {outfig}")
    except Exception as e:
        print(f"  Figure generation failed: {e}")


if __name__ == "__main__":
    main()
