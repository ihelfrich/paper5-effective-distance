"""Anderson-van Wincoop multilateral-resistance re-estimation.

Tests the conjecture in Section 6.4 of paper 5a: the AvW border-puzzle
resolution inherits the intra-national-distance regularization choice.

Design (conditional general equilibrium, YPMM ch. 2 style, endowments
held fixed):

  For each (theta, d_min) cell and the CEPII baseline:
    1. Estimate single-year (2010) PPML gravity
         X_ij = exp(beta log d_ij + gamma home_ij + eta_i + mu_j)
       with pyfixest.fepois, CRV1 clustered by pair.
    2. Form the friction composite phi_ij = exp(beta log d_ij +
       gamma home_ij). This is t_ij^(1-sigma); sigma never enters
       separately for trade-share predictions.
    3. Solve the AvW system for outward/inward resistances (Pi_i, P_j)
       by fixed-point iteration, using observed sales shares s_i and
       expenditure shares e_j.
    4. Re-solve with the home wedge removed (gamma -> 0 on intra cells).
    5. Report the conditional-GE home multiplier per country,
         M_i = X_ii(base) / X_ii(no-wedge),
       for USA, CAN, and the cross-country median.

  If the multiplier varies materially across (theta, d_min) cells, the
  AvW-style resolution is regularization-conditional, as conjectured.

Output: replication/data/avw_mr_reestimation.csv

Run:
  .venv/bin/python -u notebooks/21_avw_mr_reestimation.py
"""

from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

_spec = importlib.util.spec_from_file_location(
    "nb18f", _repo / "notebooks" / "18f_spec_b_prime_panel_ppml.py"
)
nb18f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nb18f)

OUT_DIR = _repo / "replication" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR = 2010
THETAS = (-3.0, -5.0, -7.0)
D_MINS = (0.5, 5.0, 50.0, 500.0)


def estimate_cell(panel, label):
    import pyfixest as pf

    df = panel.dropna(subset=["X_ij_usd", "distance"]).copy()
    df = df[df["distance"] > 0].copy()
    df["log_d"] = np.log(df["distance"])
    df["pair_id"] = df["iso3_o"] + "_" + df["iso3_d"]
    t0 = time.time()
    model = pf.fepois(
        "X_ij_usd ~ log_d + home | iso3_o + iso3_d",
        data=df, vcov={"CRV1": "pair_id"},
        iwls_tol=1e-8, iwls_maxiter=50,
    )
    beta = float(model.coef()["log_d"])
    gamma = float(model.coef()["home"])
    se_gamma = float(model.se()["home"])
    print(f"  [{label}] beta={beta:+.4f} gamma={gamma:+.4f}({se_gamma:.4f}) "
          f"N={len(df):,} ({time.time()-t0:.0f}s)")
    return beta, gamma, se_gamma, df


def solve_mr(phi, s, e, tol=1e-12, max_iter=5000, damp=0.5):
    """Fixed point for outward (Pi) and inward (P) resistance composites.

    phi[i, j] = t_ij^(1-sigma). Normalization P[0] = 1.
    System: Pi_i = sum_j phi_ij e_j / P_j ;  P_j = sum_i phi_ij s_i / Pi_i.
    """
    Pi = np.ones(len(s))
    P = np.ones(len(s))
    for it in range(max_iter):
        Pi_new = phi @ (e / P)
        P_new = phi.T @ (s / Pi_new)
        P_new = P_new / P_new[0]
        diff = max(np.max(np.abs(np.log(Pi_new / Pi))),
                   np.max(np.abs(np.log(P_new / P))))
        Pi = damp * Pi_new + (1 - damp) * Pi
        P = damp * P_new + (1 - damp) * P
        if diff < tol:
            return Pi, P, it
    return Pi, P, max_iter


def ge_home_multiplier(df, beta, gamma):
    """Conditional-GE intra-trade ratio, with vs without the home wedge."""
    intra_iso = sorted(df.loc[df["home"] == 1, "iso3_o"].unique())
    sub = df[df["iso3_o"].isin(intra_iso) & df["iso3_d"].isin(intra_iso)]
    mat_d = sub.pivot_table(index="iso3_o", columns="iso3_d",
                            values="distance", aggfunc="first")
    mat_d = mat_d.reindex(index=intra_iso, columns=intra_iso)
    # Require a complete distance matrix; iteratively drop the country
    # contributing the most missing cells.
    while mat_d.isna().any().any():
        worst = (mat_d.isna().sum(axis=1) + mat_d.isna().sum(axis=0)).idxmax()
        mat_d = mat_d.drop(index=worst, columns=worst)
    iso = list(mat_d.index)

    flows = sub.pivot_table(index="iso3_o", columns="iso3_d",
                            values="X_ij_usd", aggfunc="first")
    flows = flows.reindex(index=iso, columns=iso).fillna(0.0)
    X = flows.values
    keep = (X.sum(axis=1) > 0) & (X.sum(axis=0) > 0)
    iso = [c for c, k in zip(iso, keep) if k]
    X = X[np.ix_(keep, keep)]
    Y, E = X.sum(axis=1), X.sum(axis=0)
    D = mat_d.loc[iso, iso].values
    n = len(iso)
    Yw = Y.sum()
    s, e = Y / Yw, E / Yw

    home = np.eye(n)
    phi_base = np.exp(beta * np.log(D) + gamma * home)
    phi_cf = np.exp(beta * np.log(D))  # home wedge removed

    Pi_b, P_b, it_b = solve_mr(phi_base, s, e)
    Pi_c, P_c, it_c = solve_mr(phi_cf, s, e)

    diag_b = np.diag(phi_base) / (Pi_b * P_b)
    diag_c = np.diag(phi_cf) / (Pi_c * P_c)
    mult = diag_b / diag_c
    out = dict(zip(iso, mult))
    return {
        "n_countries": n,
        "iters": int(max(it_b, it_c)),
        "mult_median": float(np.median(mult)),
        "mult_USA": float(out.get("USA", np.nan)),
        "mult_CAN": float(out.get("CAN", np.nan)),
        "mult_p25": float(np.percentile(mult, 25)),
        "mult_p75": float(np.percentile(mult, 75)),
    }


def main():
    print("=" * 70)
    print("nb21: AvW multilateral-resistance re-estimation (2010, cond. GE)")
    print("=" * 70)
    area_lookup = nb18f.country_areas()
    intra_panel, inter_panel = nb18f.build_year_panel(YEAR, area_lookup)
    inter_p = inter_panel.copy()
    inter_p["distance"] = inter_p["distw_harmonic"]

    cells = [("CEPII_baseline", None, None, "d_HM_cepii")]
    for th in THETAS:
        for dm in D_MINS:
            cells.append((f"theta{th}_dmin{dm}", th, dm, None))

    rows = []
    for label, theta, d_min, base_col in cells:
        intra_c = intra_panel.copy()
        if base_col is not None:
            intra_c["distance"] = intra_c[base_col]
        else:
            intra_c["distance"] = intra_c["area_km2"].apply(
                lambda a: nb18f.country_d_eff(a, d_min, theta)
                if np.isfinite(a) else np.nan)
        panel = pd.concat([inter_p, intra_c], ignore_index=True, sort=False)
        try:
            beta, gamma, se_gamma, df = estimate_cell(panel, label)
            ge = ge_home_multiplier(df, beta, gamma)
        except Exception as exc:
            print(f"  [{label}] FAILED: {type(exc).__name__}: {exc}")
            rows.append({"label": label, "theta": theta, "d_min": d_min,
                         "converged": False})
            continue
        rows.append({"label": label, "theta": theta, "d_min": d_min,
                     "beta": beta, "gamma": gamma, "se_gamma": se_gamma,
                     "cond_mult_exp_gamma": float(np.exp(gamma)),
                     "converged": True, **ge})
        pd.DataFrame(rows).to_csv(OUT_DIR / "avw_mr_reestimation.csv",
                                  index=False)
        print(f"    GE home multiplier: median={ge['mult_median']:.4g} "
              f"USA={ge['mult_USA']:.4g} CAN={ge['mult_CAN']:.4g} "
              f"(n={ge['n_countries']}, iters={ge['iters']})")

    print("\nSaved replication/data/avw_mr_reestimation.csv")


if __name__ == "__main__":
    main()
