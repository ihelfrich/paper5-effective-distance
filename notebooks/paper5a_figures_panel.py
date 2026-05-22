"""Regenerate Figure 2 (home-bias vs d_min) from the PANEL PPML output.

To run after nb18f completes.

Inputs:
  /Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived/spec_b_prime_panel_ppml.csv

Outputs:
  /Users/ian/1A_Helfrich_ThesisResearch_May2024/Paper5_EffectiveDistance/manuscript/paper_5a/figures/
    fig2_panel_home_bias_vs_dmin.{pdf,png}
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "Charter",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

COLORS = {
    "INK": "#1d1d1d",
    "RUST": "#B5482A",
    "SAGE": "#4C7066",
    "GOLD": "#C28B25",
    "VIOLET": "#5E4F7A",
    "DIM": "#7a7a7a",
}

REPO = Path(__file__).resolve().parents[1]
DRV = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")
FIG_OUT = REPO / "manuscript" / "paper_5a" / "figures"
FIG_OUT.mkdir(parents=True, exist_ok=True)


def fig2_panel_home_bias():
    path = DRV / "spec_b_prime_panel_ppml.csv"
    if not path.exists():
        print(f"[skip] {path} does not exist yet")
        return

    df = pd.read_csv(path)
    b = df[df["spec"] == "B_prime_panel"].copy()
    baseline = df[df["spec"] == "A1_panel_baseline_cepii"].copy()
    if len(b) == 0:
        print("[skip] no B_prime_panel rows yet")
        return

    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    theta_colors = {-3.0: COLORS["GOLD"], -5.0: COLORS["RUST"], -7.0: COLORS["VIOLET"]}

    for theta in sorted(b["theta"].unique()):
        sub = b[b["theta"] == theta].sort_values("d_min")
        if len(sub) == 0:
            continue
        ax.plot(sub["d_min"], sub["coef_home"],
                marker="o", color=theta_colors[theta],
                label=fr"$\theta = {int(theta)}$ (structurally consistent)")

    if len(baseline) > 0:
        base_home = float(baseline["coef_home"].iloc[0])
        ax.axhline(base_home, linestyle="--", color=COLORS["INK"],
                   linewidth=0.9,
                   label=fr"CEPII baseline ($\theta = -1$): $\hat\gamma = {base_home:+.2f}$")

    ax.axhline(0, linestyle=":", color=COLORS["DIM"], linewidth=0.5)

    ax.set_xscale("log")
    ax.set_xlabel(r"Regularization grain $d_{\min}$ (km)")
    ax.set_ylabel(r"Panel PPML home-bias coefficient $\hat{\gamma}$")
    ax.set_title(r"Pooled panel PPML home-bias coefficient $\hat{\gamma}$"
                 + "\n"
                 + r"vs $d_{\min}$, structurally consistent $\theta$, 2010-2015-2020 stacked",
                 fontsize=10, pad=8)
    ax.legend(loc="best", frameon=False, fontsize=8.5)

    out_pdf = FIG_OUT / "fig2_panel_home_bias_vs_dmin.pdf"
    out_png = FIG_OUT / "fig2_panel_home_bias_vs_dmin.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_pdf}")
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    fig2_panel_home_bias()
