"""paper5 — Effective Distance: a satellite-calibrated bilateral trade-cost panel.

Companion code for Gonchar & Helfrich (2026), "Effective Distance: A
Satellite-Calibrated Bilateral Trade-Cost Panel, 2000–2024."

Modules
-------
data_loaders : BACI, CEPII Gravity, WorldPop, GHS-POP, VIIRS, OSM, WPI, OpenFlights.
distance     : Multi-modal least-cost-path distance construction; 8 variants.
gravity      : pyfixest three-way FE PPML with Weidner-Zylkin bias correction.
counterfactual : ACR + Caliendo-Parro exact-hat algebra for chokepoint shocks.

Conventions
-----------
- All bilateral objects indexed as (i, j, t) with i = origin ISO3, j = destination ISO3, t = year.
- All distances in km unless suffixed _hours or _usd.
- Monetary values in current USD unless suffixed _real.
- Coordinate reference system: EPSG:4326 for raw rasters/shapes; EPSG:8857 (Equal Earth)
  for area and distance summaries.

Data paths (see config.py)
--------------------------
HELFRICH_GD = /Volumes/HELFRICH-GD
ONEDRIVE    = /Users/ian/Library/CloudStorage/OneDrive-Personal
PAPER5_DATA = ./data  (intermediate parquet/zarr; gitignored)
"""

__version__ = "0.1.0.dev"
__authors__ = ("Elizaveta Gonchar", "Ian T. Helfrich")
