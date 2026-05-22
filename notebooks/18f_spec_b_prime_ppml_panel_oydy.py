"""Pooled PPML Spec B-prime with origin-year and destination-year FE.

This is a non-overwriting follow-up to nb18e. The nb18e artifact on disk is a
single-year 2010 PPML check. This script builds the three-year panel used in
nb18d (2010, 2015, 2020), keeps zero trade flows, and estimates:

    X_ij,t = exp(beta log d_ij + gamma home_ij + origin-year FE
                 + destination-year FE) * eps_ij,t

using fixest::fepois. It writes timestamped input and result CSVs so reruns do
not clobber previous artifacts.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo / "src"))

OUT_DIR = Path("/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/data_derived")
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEARS = [2010, 2015, 2020]
D_MINS = [0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0]

DISK_SIM_T5 = {
    0.01: 0.063, 0.1: 0.631, 1: 6.288, 5: 29.449, 10: 51.855,
    50: 149.742, 100: 229.350, 200: 352.875, 500: 640.413,
}
DISK_R = 1000.0


def dist_col_name(d_min: float) -> str:
    return "dist_dmin_" + str(d_min).replace(".", "p")


def country_d_eff(area_km2: float, d_min_km: float) -> float:
    if not np.isfinite(area_km2) or area_km2 <= 0:
        return np.nan
    r_country = float(np.sqrt(area_km2 / np.pi))
    d_min_disk_equiv = d_min_km * (DISK_R / r_country)
    keys = sorted(DISK_SIM_T5.keys())
    vals = [DISK_SIM_T5[k] for k in keys]
    d_min_disk_equiv = float(np.clip(d_min_disk_equiv, keys[0], keys[-1]))
    d_eff_disk = float(np.interp(d_min_disk_equiv, keys, vals))
    d_eff_country = d_eff_disk * (r_country / DISK_R)
    return float(min(d_eff_country, 0.9 * r_country))


def country_areas() -> dict[str, float]:
    from paper5.region_raster import load_country_boundaries

    bnd = load_country_boundaries().to_crs("ESRI:54009")
    return dict(zip(bnd["ISO3"], (bnd.geometry.area / 1e6).astype(float).values))


def chunked_baci(year: int) -> pd.DataFrame | None:
    f = Path(f"/Volumes/HELFRICH-GD/TradeData/BACI_HS02_V202401b/BACI_HS02_Y{year}_V202401b.csv")
    if not f.exists():
        return None
    pieces = []
    for chunk in pd.read_csv(f, usecols=["i", "j", "v"], chunksize=1_000_000):
        chunk["v"] = pd.to_numeric(chunk["v"], errors="coerce")
        chunk = chunk.dropna(subset=["v"])
        pieces.append(chunk.groupby(["i", "j"], as_index=False)["v"].sum())
    agg = pd.concat(pieces).groupby(["i", "j"], as_index=False)["v"].sum()
    agg["X_ij_usd"] = agg["v"] * 1000
    return agg


def load_cepii_years(years: list[int]) -> dict[int, pd.DataFrame]:
    f = Path("/Volumes/HELFRICH-GD/TradeData/Gravity_csv_V202211/Gravity_V202211.csv")
    usecols = [
        "year", "iso3num_o", "iso3num_d", "iso3_o", "iso3_d",
        "gdp_o", "distw_harmonic",
    ]
    pieces = []
    for chunk in pd.read_csv(
        f,
        usecols=usecols,
        dtype={"iso3_o": str, "iso3_d": str, "year": "Int64"},
        chunksize=500_000,
        low_memory=False,
    ):
        sub = chunk[chunk["year"].isin(years)]
        if len(sub) > 0:
            pieces.append(sub)
    if not pieces:
        return {}
    all_rows = pd.concat(pieces, ignore_index=True)
    return {year: all_rows[all_rows["year"] == year].copy() for year in years}


def build_year_panel(year: int, cepii: pd.DataFrame, area_lookup: dict[str, float]) -> pd.DataFrame:
    print(f"[{year}] Loading BACI...")
    baci = chunked_baci(year)
    if baci is None:
        raise FileNotFoundError(f"Missing BACI file for {year}")

    map_o = cepii[["iso3num_o", "iso3_o"]].dropna().drop_duplicates(subset=["iso3num_o"])
    map_d = cepii[["iso3num_d", "iso3_d"]].dropna().drop_duplicates(subset=["iso3num_d"])
    num_to_iso = {int(n): s for n, s in map_o.values}
    num_to_iso.update({int(n): s for n, s in map_d.values})

    baci["iso3_o"] = baci["i"].map(num_to_iso)
    baci["iso3_d"] = baci["j"].map(num_to_iso)
    baci = baci.dropna(subset=["iso3_o", "iso3_d"])

    cepii_self = cepii[cepii["iso3_o"] == cepii["iso3_d"]][
        ["iso3_o", "distw_harmonic"]
    ].copy()
    cepii_self = cepii_self.rename(columns={"distw_harmonic": "dist_baseline_cepii"})

    gdp_lookup = dict(cepii[["iso3_o", "gdp_o"]].drop_duplicates().dropna().values)
    exp_by_o = (
        baci.groupby("iso3_o", as_index=False)["X_ij_usd"].sum()
        .set_index("iso3_o")["X_ij_usd"]
        .to_dict()
    )

    intra_rows = []
    for iso, gdp in gdp_lookup.items():
        if not np.isfinite(gdp) or gdp <= 0:
            continue
        exp = exp_by_o.get(iso, 0.0)
        area = area_lookup.get(iso, np.nan)
        row = {
            "year": year,
            "iso3_o": iso,
            "iso3_d": iso,
            "X_ij_usd": max(0.0, float(gdp) - float(exp)),
            "home": 1,
            "area_km2": area,
        }
        for d_min in D_MINS:
            row[dist_col_name(d_min)] = country_d_eff(area, d_min)
        intra_rows.append(row)

    intra = pd.DataFrame(intra_rows).merge(cepii_self, on="iso3_o", how="left")

    inter = cepii[cepii["iso3_o"] != cepii["iso3_d"]][
        ["year", "iso3_o", "iso3_d", "distw_harmonic"]
    ].copy()
    inter = inter.merge(
        baci[["iso3_o", "iso3_d", "X_ij_usd"]],
        on=["iso3_o", "iso3_d"],
        how="left",
    )
    inter["X_ij_usd"] = inter["X_ij_usd"].fillna(0.0)
    inter["home"] = 0
    inter["area_km2"] = np.nan
    inter["dist_baseline_cepii"] = inter["distw_harmonic"]
    for d_min in D_MINS:
        inter[dist_col_name(d_min)] = inter["distw_harmonic"]

    panel = pd.concat([inter, intra], ignore_index=True, sort=False)
    panel["oy"] = panel["iso3_o"].astype(str) + "_" + panel["year"].astype(str)
    panel["dy"] = panel["iso3_d"].astype(str) + "_" + panel["year"].astype(str)
    print(
        f"[{year}] rows={len(panel):,} intra={int(panel['home'].sum()):,} "
        f"zeros={int((panel['X_ij_usd'] == 0).sum()):,}"
    )
    return panel


def run_fixest(panel_csv: Path, out_csv: Path) -> None:
    specs = [
        {"label": "PPML_panel_baseline_cepii", "d_min": None, "dist_col": "dist_baseline_cepii"}
    ] + [
        {"label": f"PPML_panel_theta-5_dmin{d_min}km", "d_min": d_min, "dist_col": dist_col_name(d_min)}
        for d_min in D_MINS
    ]

    r_code = f"""
    suppressPackageStartupMessages(library(fixest))
    panel <- read.csv({json.dumps(str(panel_csv))}, stringsAsFactors = FALSE)
    specs <- data.frame(
      label = c({",".join(json.dumps(s["label"]) for s in specs)}),
      d_min = c({",".join("NA" if s["d_min"] is None else str(s["d_min"]) for s in specs)}),
      dist_col = c({",".join(json.dumps(s["dist_col"]) for s in specs)}),
      stringsAsFactors = FALSE
    )
    rows <- list()
    for (i in seq_len(nrow(specs))) {{
      spec <- specs[i,]
      dist_col <- spec$dist_col
      dat <- panel[is.finite(panel[[dist_col]]) & panel[[dist_col]] > 0 &
                   !is.na(panel$X_ij_usd) & is.finite(panel$X_ij_usd), ]
      dat$log_d <- log(dat[[dist_col]])
      fit <- tryCatch(
        fepois(X_ij_usd ~ log_d + home | oy + dy, data = dat, notes = FALSE),
        error = function(e) e
      )
      if (inherits(fit, "error")) {{
        rows[[i]] <- data.frame(
          label = spec$label, d_min = spec$d_min, n_obs = nrow(dat),
          n_obs_fit = NA_integer_, n_intra = sum(dat$home == 1),
          n_zero = sum(dat$X_ij_usd == 0), coef_log_d = NA_real_,
          se_log_d = NA_real_, coef_home = NA_real_, se_home = NA_real_,
          home_multiplier = NA_real_, error = fit$message
        )
      }} else {{
        cf <- coef(fit)
        se <- se(fit)
        gamma <- unname(cf["home"])
        beta <- unname(cf["log_d"])
        rows[[i]] <- data.frame(
          label = spec$label, d_min = spec$d_min, n_obs = nrow(dat),
          n_obs_fit = nobs(fit), n_intra = sum(dat$home == 1),
          n_zero = sum(dat$X_ij_usd == 0), coef_log_d = beta,
          se_log_d = unname(se["log_d"]), coef_home = gamma,
          se_home = unname(se["home"]), home_multiplier = exp(gamma),
          error = NA_character_
        )
      }}
    }}
    out <- do.call(rbind, rows)
    write.csv(out, {json.dumps(str(out_csv))}, row.names = FALSE)
    print(out[, c("label", "d_min", "n_obs", "n_obs_fit", "n_intra",
                  "n_zero", "coef_log_d", "coef_home", "home_multiplier")],
          row.names = FALSE)
    """
    subprocess.run(["Rscript", "-e", r_code], check=True)


def main() -> None:
    run_id = time.strftime("%H%M%S")
    panel_csv = OUT_DIR / f"spec_b_prime_ppml_panel_oydy_input_{run_id}.csv"
    out_csv = OUT_DIR / f"spec_b_prime_ppml_panel_oydy_{run_id}.csv"
    if panel_csv.exists() or out_csv.exists():
        raise FileExistsError("Refusing to overwrite an existing output file")

    print("=" * 70)
    print("18f: pooled PPML with origin-year and destination-year FE")
    print(f"Years: {YEARS}")
    print("=" * 70)

    print("[0] Loading country areas...")
    area_lookup = country_areas()
    print(f"    {len(area_lookup)} countries")

    print("[1] Loading CEPII Gravity rows...")
    cepii_by_year = load_cepii_years(YEARS)

    panels = [build_year_panel(year, cepii_by_year[year], area_lookup) for year in YEARS]
    panel = pd.concat(panels, ignore_index=True, sort=False)
    panel.to_csv(panel_csv, index=False)
    print(f"[2] Saved pooled panel input: {panel_csv}")
    print(
        f"[2] pooled rows={len(panel):,} intra={int(panel['home'].sum()):,} "
        f"zeros={int((panel['X_ij_usd'] == 0).sum()):,}"
    )

    print("[3] Estimating PPML models with fixest::fepois...")
    run_fixest(panel_csv, out_csv)
    print(f"[4] Saved estimates: {out_csv}")


if __name__ == "__main__":
    main()
