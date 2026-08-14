# Paper Review Targets

This directory is the active PeerJ Computer Science project for scReg-Eval. It is separate from `research/sc-fm-benchmark`, which is a different paper with its own PeerJ and IEEE venue variants.

## Human review target

Canonical PeerJ CS PDF (source of truth for the compiled report):

- `submission_peerj/flat_upload/manuscript.pdf` — flat-upload manuscript with
  `Figure1.pdf` (study design) through `Figure13.pdf` (coverage QC) plus
  `FigureA1.pdf`–`FigureA3.pdf`. Checksums: `SHA256SUMS.txt`.

The current editable source is:

- `manuscript.tex` — canonical editable manuscript source. Do not leave a
  competing `manuscript.pdf` in this directory.

## Package structure

- `submission_peerj/source/` — self-contained TeX source mirror. Its PDF is generated build output and is not the review target.
- `submission_peerj/flat_upload/` — upload-ready package and the only current PeerJ manuscript PDF for visual review.
- `submission_peerj/internal/figure_build/` — temporary standalone figure wrappers and PDFs; generated only.
- `submission_peerj/*.zip` — historical submission bundles; do not use for current review.
- `../release_candidate/` — historical audit capsules; do not use for current review.

## Excluded from current review

Do not inspect deleted/generated previews or historical archives as if they were current evidence. Rebuild from `manuscript.tex` with `src/v2/build_peerj_package.py` before creating a new package.
