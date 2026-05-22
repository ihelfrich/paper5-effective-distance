"""Figure generation for Paper 5a (the wrong-θ pole paper).

Two figures:
  Fig 1: Regularization sweep — d_eff(θ=-5) on uniform disk as a function of d_min.
         Shows that the implied intra-national distance is essentially linear
         in the regularization grain, with no d_min-independent fixed point.

  Fig 2: Gravity-regression home-bias coefficient as a function of d_min,
         faceted by θ and year. From nb18d outputs.

Inputs:
  /Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived/regularization_dmin_cutoff.csv
  /Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived/spec_b_prime_fixed.csv

Outputs (saved to manuscript/paper_5a/figures/):
  fig1_regularization_curve.pdf  +  .png
  fig2_home_bias_vs_dmin.pdf     +  .png
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

# Color palette: matches Ian's PoI palette (INK, RUST, SAGE, GOLD, VIOLET, DIM)
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


# ── Figure 1: Regularization curve ──────────────────────────────────────────

def fig1_regularization_curve():
    df = pd.read_csv(DRV / "regularization_dmin_cutoff.csv")

    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    # Main curve: d_eff vs d_min (both in km)
    ax.plot(df["d_min_km"], df["d_eff_km"],
            marker="o", color=COLORS["RUST"],
            label=r"$d_{ii}^{\rm eff}(\theta=-5)$, uniform disk $R=1000$ km")

    # Reference: Head-Mayer closed-form value at theta=-1 (analytical)
    d_HM_theta_neg1 = 500.0  # (1/2) * R for R=1000
    ax.axhline(d_HM_theta_neg1, linestyle="--", color=COLORS["SAGE"],
               linewidth=1.0,
               label=r"HM at $\theta=-1$: $d_{ii} = R/2 = 500$ km")

    # Reference: Head-Mayer at theta=+1 (arithmetic, 2R/3)
    d_HM_theta_pos1 = 2000.0 / 3.0
    ax.axhline(d_HM_theta_pos1, linestyle=":", color=COLORS["GOLD"],
               linewidth=1.0,
               label=r"HM at $\theta=+1$: $d_{ii} = 2R/3 \approx 667$ km")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Regularization grain $d_{\min}$ (km)")
    ax.set_ylabel(r"Implied $d_{ii}^{\rm eff}$ (km)")
    ax.set_title(r"Intra-national distance at gravity-consistent $\theta = -5$"
                 + "\n"
                 + r"on a uniform disk ($R = 1000$ km), vs grid resolution",
                 fontsize=10, pad=10)
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
    ax.set_xlim(0.005, 700)
    ax.set_ylim(0.05, 800)

    # Annotate the structural point
    ax.annotate(
        "No fixed value: any specific $d_{ii}$\nat $\\theta \\leq -2$ is conditional\non the discretization grain.",
        xy=(20, 100), xytext=(0.07, 200),
        fontsize=8.5, color=COLORS["INK"],
        ha="left",
        arrowprops=dict(arrowstyle="-", color=COLORS["DIM"], lw=0.6,
                       connectionstyle="arc3,rad=-0.2"),
    )

    out_pdf = FIG_OUT / "fig1_regularization_curve.pdf"
    out_png = FIG_OUT / "fig1_regularization_curve.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_pdf}")
    print(f"  wrote {out_png}")


# ── Figure 2: Home-bias coefficient vs d_min ────────────────────────────────

def fig2_home_bias_vs_dmin():
    path = DRV / "spec_b_prime_fixed.csv"
    if not path.exists():
        print(f"  [skip] {path} does not exist yet (nb18d not finished)")
        return

    df = pd.read_csv(path)
    b = df[df["spec"] == "B_prime"].copy()
    if len(b) == 0:
        print("  [skip] no Spec B-prime rows yet")
        return

    baselines = df[df["spec"] == "A1_baseline_cepii_intra"].copy()

    years = sorted(b["year"].unique())
    thetas = sorted(b["theta"].unique())

    fig, axes = plt.subplots(1, len(years), figsize=(4.5 * len(years), 4.2),
                             sharey=True, squeeze=False)
    axes = axes[0]

    theta_colors = {-3.0: COLORS["GOLD"], -5.0: COLORS["RUST"], -7.0: COLORS["VIOLET"]}

    for ax, year in zip(axes, years):
        b_y = b[b["year"] == year]
        base_y = baselines[baselines["year"] == year]
        if len(base_y) > 0:
            base_home = float(base_y["coef_home"].iloc[0])
            ax.axhline(base_home, linestyle="--", color=COLORS["INK"],
                       linewidth=0.9, label=f"baseline (CEPII): {base_home:+.2f}")

        for theta in thetas:
            sub = b_y[b_y["theta"] == theta].sort_values("d_min")
            if len(sub) == 0:
                continue
            ax.plot(sub["d_min"], sub["coef_home"],
                    marker="o", color=theta_colors[theta],
                    label=fr"$\theta = {int(theta)}$")

        ax.set_xscale("log")
        ax.set_xlabel(r"$d_{\min}$ (km)")
        ax.set_title(f"Year {year}", fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel(r"Home-bias coefficient $\hat{\gamma}$")
        ax.legend(loc="best", frameon=False, fontsize=8.5)

    fig.suptitle(r"Gravity-regression home-bias coefficient $\hat{\gamma}$"
                 + " across regularization grain and structural θ\n"
                 + r"(BACI-CEPII panel, $\log X_{ij} = \alpha + \beta \log d_{ij} + \gamma \mathbf{1}\{i=j\} + \eta_i + \mu_j$)",
                 fontsize=10, y=1.02)

    out_pdf = FIG_OUT / "fig2_home_bias_vs_dmin.pdf"
    out_png = FIG_OUT / "fig2_home_bias_vs_dmin.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_pdf}")
    print(f"  wrote {out_png}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("Figure 1: Regularization curve")
    fig1_regularization_curve()
    print("\nFigure 2: Home-bias vs d_min")
    fig2_home_bias_vs_dmin()


if __name__ == "__main__":
    main()
