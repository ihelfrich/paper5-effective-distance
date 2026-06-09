# Data licensing

Three different licenses apply to the contents of this repository.

## Code

All source code (Python in `src/` and `notebooks/`, Rust in `crates/`) is
released under the MIT License. See `LICENSE` in the repository root.

## Derived data

Derived data tables produced by this project (the contents of
`data/derived/`, including the effective-distance tables
`real_country_internal_distance.csv` and `real_country_eff_vs_cepii.csv`)
are released under the Creative Commons Attribution 4.0 International
license (CC-BY-4.0): https://creativecommons.org/licenses/by/4.0/

If you use the derived tables, cite the working paper listed in
`CITATION.cff`.

## Raw inputs

Raw input datasets are NOT redistributed in this repository and remain
under their original providers' terms:

- BACI (CEPII): http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37
- CEPII Gravity and GeoDist: http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele.asp
- WorldPop: https://www.worldpop.org/ (CC-BY-4.0)
- GHS-POP and GHS-SMOD (JRC GHSL): https://human-settlement.emergency.copernicus.eu/

The reproduction pipeline expects these to be downloaded separately; see
`DATA_INVENTORY.md` for paths and versions.
