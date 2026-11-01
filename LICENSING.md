# Licensing

This project uses two licenses, decided by the project owner on 2026-07-31:

- **Original code**: `src/`, `validate_artifacts.py`, and, in the Zenodo archive, the
  scripts under `figures/` (`*.R`, `*.py`, `*.sh`): [MIT License](LICENSE).
- **Documentation, figure sources, and derived public results**: `results/` and, in the
  Zenodo archive, `docs/` and the TikZ sources and panel data under `figures/`:
  [Creative Commons Attribution 4.0 International](LICENSE-CONTENT.md)
  (CC BY 4.0, <https://creativecommons.org/licenses/by/4.0/>).

- **Article materials**: `article_materials/` follows the same split; its own
  [Licensing](article_materials/LICENSING.md) note lists which parts are code, which are results, and
  the vendored model code it carries.

Third-party components are not relicensed by this project:

- Vendored model code keeps its upstream license: scFoundation (`src/v2/scf_vendor/`,
  Apache License 2.0, from `biomap-research/scFoundation`) and UCE (`src/v2/uce_vendor/`,
  MIT, from `snap-stanford/UCE`). Each directory carries the upstream license text.
- Foundation-model weights and upstream datasets (GSE174367, 10x PBMC multiome,
  GSE206767) are not redistributed; they remain governed by their own terms and
  are referenced by accession only.
