"""gravity.py — pyfixest three-way FE PPML with Weidner-Zylkin bias correction.

Replaces TradeModels_Final.estimate_gravity_models (statsmodels Poisson GLM)
with a spec that matches the Yotov-Piermartini-Monteiro-Larch (2016) gold
standard:

    X_{ij,t} = exp( π_{i,t} + π_{j,t} + μ_{ij} + β·ln d_{ij,t} + γ'·Z_{ij,t} ) + ε

- Three-way FE (exporter-time, importer-time, pair).
- Intranational flows X_{ii,t} included (Heid-Larch-Yotov 2021).
- PPML via pyfixest.fepois (Santos Silva-Tenreyro 2006).
- Weidner-Zylkin (2021) bias correction available via Stata bridge.

This module is scaffold-only as of v0.1; filled in during sprint days 5–7.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Optional

import pandas as pd

# pyfixest will be imported inside functions to keep module-import cheap
# during scoping/dry-run. See pyproject.toml for pinned version.


@dataclass
class GravitySpec:
    """Structural gravity specification."""

    distance_var: str  # e.g. "d_lcp_multi_t"
    controls: tuple[str, ...] = ()
    fe_exporter_time: bool = True
    fe_importer_time: bool = True
    fe_pair: bool = True
    include_intranational: bool = True
    theta_prior: Optional[float] = None  # sector-specific θ_s if set


def estimate_gravity(
    df: pd.DataFrame,
    spec: GravitySpec,
    cluster: Literal["pair", "exp-imp-time"] = "pair",
) -> "pyfixest.FixestResults":  # type: ignore[name-defined]
    """Estimate the three-way FE PPML gravity spec.

    Parameters
    ----------
    df : DataFrame with columns
        iso_o, iso_d, year, trade_value, {spec.distance_var}, *spec.controls
    spec : GravitySpec
    cluster : SE clustering.

    Returns
    -------
    pyfixest result object with .coef(), .se(), .tstat(), etc.
    """
    import pyfixest as pf

    rhs_terms: list[str] = [f"np.log({spec.distance_var})"]
    rhs_terms.extend(spec.controls)

    fe_terms: list[str] = []
    if spec.fe_exporter_time:
        fe_terms.append("iso_o^year")
    if spec.fe_importer_time:
        fe_terms.append("iso_d^year")
    if spec.fe_pair:
        fe_terms.append("iso_o^iso_d")

    formula = (
        f"trade_value ~ {' + '.join(rhs_terms)}"
        + (f" | {' + '.join(fe_terms)}" if fe_terms else "")
    )

    if cluster == "pair":
        cluster_col = "pair_id"
        df = df.assign(pair_id=df["iso_o"] + "_" + df["iso_d"])
    elif cluster == "exp-imp-time":
        # Two-way cluster via pyfixest's vcov='twoway' construct
        cluster_col = None
    else:
        raise ValueError(cluster)

    return pf.fepois(
        fml=formula,
        data=df,
        vcov={"CRV1": cluster_col} if cluster_col else "twoway",
    )


def benchmark_distance_variants(
    df: pd.DataFrame,
    distance_vars: Iterable[str],
    controls: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Re-estimate gravity across distance variants; return comparison table.

    Returns DataFrame with columns: variant, beta_dist, se, pseudo_r2, aic, n_obs.
    """
    rows = []
    for d in distance_vars:
        spec = GravitySpec(distance_var=d, controls=controls)
        res = estimate_gravity(df, spec)
        rows.append(
            {
                "variant": d,
                "beta_dist": float(res.coef()[f"np.log({d})"]),
                "se": float(res.se()[f"np.log({d})"]),
                "pseudo_r2": getattr(res, "pseudo_r2", float("nan")),
                "aic": getattr(res, "aic", float("nan")),
                "n_obs": int(res._N),
            }
        )
    return pd.DataFrame(rows)


def weidner_zylkin_bias_correct(
    stata_result_path: Path,
) -> pd.DataFrame:
    """Read Weidner-Zylkin bias-corrected estimates from Stata ppml_fe_bias output.

    We run ppml_fe_bias in Stata (see scripts/wz_bias.do) and read the
    exported bias-corrected coefficients here.
    """
    # Implemented in sprint day 6.
    raise NotImplementedError("WZ bridge pending; implement after pyfixest main spec.")
