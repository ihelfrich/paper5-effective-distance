"""cli.py — Command-line entry points for the Paper 5 pipeline.

Called by the Snakemake DAG (workflow/Snakefile) via:
    python -m paper5.cli <subcommand> [options]

Subcommands:
    agglomerate       L1: raster → agglomeration centroids parquet
    compute-lcp       L2: multi-modal LCP via paper5_core Rust extension
    assemble-panel    L3: year-shards → distance_panel.parquet

Each subcommand is a thin wrapper that imports the relevant module and
delegates; the heavy lifting lives in agglomerate.py, distance.py.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def cmd_agglomerate(args: argparse.Namespace) -> None:
    from paper5.agglomerate import build_all_agglomerations, AgglomerationConfig
    df = build_all_agglomerations(
        year=args.year,
        iso_list=Path(args.iso_list).read_text().split(),
        worldpop_path=Path(args.worldpop),
        viirs_path=Path(args.viirs),
        boundaries_path=Path(args.boundaries),
        config=AgglomerationConfig(),
        out_path=Path(args.out),
    )
    print(f"[agglomerate] wrote {len(df)} rows → {args.out}", file=sys.stderr)


def cmd_compute_lcp(args: argparse.Namespace) -> None:
    from paper5.distance import compute_d_lcp_multi_t
    import polars as pl
    agglom = pl.read_parquet(args.agglom)
    # Stub: delegates to Rust extension via distance.py (sprint day 4).
    raise NotImplementedError(
        "compute-lcp delegates to paper5_core Rust extension; "
        "implement in distance.py sprint day 4."
    )


def cmd_assemble_panel(args: argparse.Namespace) -> None:
    from paper5.distance import compute_all_variants
    raise NotImplementedError("Sprint day 5.")


def main() -> None:
    p = argparse.ArgumentParser(prog="python -m paper5.cli")
    sub = p.add_subparsers(dest="cmd", required=True)

    # ---- agglomerate ----
    a = sub.add_parser("agglomerate")
    a.add_argument("--year", type=int, required=True)
    a.add_argument("--worldpop", required=True)
    a.add_argument("--viirs", required=True)
    a.add_argument("--boundaries", default="data/boundaries/gadm_410.gpkg")
    a.add_argument("--iso-list", required=True)
    a.add_argument("--out", required=True)

    # ---- compute-lcp ----
    b = sub.add_parser("compute-lcp")
    b.add_argument("--year", type=int, required=True)
    b.add_argument("--agglom", required=True)
    b.add_argument("--out", required=True)

    # ---- assemble-panel ----
    c = sub.add_parser("assemble-panel")
    c.add_argument("--out", required=True)

    args = p.parse_args()
    dispatch = {
        "agglomerate": cmd_agglomerate,
        "compute-lcp": cmd_compute_lcp,
        "assemble-panel": cmd_assemble_panel,
    }
    dispatch[args.cmd](args)


if __name__ == "__main__":
    main()
