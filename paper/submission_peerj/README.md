# scReg-Eval PeerJ submission package

Generated from the canonical manuscript by `src/v2/build_peerj_package.py`.
The package contains 13 pooled FM/readout rows, thirteen main-text figures,
three appendix figures, and eight tables. The builder first
regenerates `paper/panel_data.json` with `make_panel_data.py` and all canonical figure
and table fragments with `Rscript make_figs.R`, then rebuilds both source and flat mirrors.
The canonical human-review PDF is `flat_upload/manuscript.pdf`. Do not leave a competing `paper/manuscript.pdf`.

## Layout

- `source/`: self-contained editable article source, including the thirteen main-text
  figure fragments, three appendix fragments under `figs_extension/`, and eight tables.
- `flat_upload/`: upload-ready package with `manuscript.pdf`, `Figure1.pdf` through
  `Figure13.pdf`, `FigureA1.pdf`–`FigureA3.pdf`, table fragments, references, and checksums.
  Printed Figure 1 is study design (`Figure1.pdf`); printed Figure 13 is coverage QC
  (`Figure13.pdf`).
- `internal/figure_build/`: temporary standalone wrappers used to export the figure PDFs; generated only.
- `supplemental/`: PeerJ CS AI-in-code Supplemental Data S1–S3 (not members of `upload.zip`). Rebuild with `python3 src/v2/build_ai_disclosure_zips.py`.
- `HUMAN_GATES.md`: remaining author decisions and external actions (DONE/OPEN status).
- `PEERJ_FORM_TEMPLATE.md`: scientific fields.
- `PEERJ_SUBMIT_PASTE.md`: copy-ready blocks for the PeerJ CS online form.
- `SHA256SUMS.txt`: hashes for all files in `flat_upload/`.
- `supplemental/SHA256SUMS_AI.txt`: hashes for Supplemental Data S1–S3.

The flat manuscript is the current PeerJ upload variant. Convert and verify it against the current official PeerJ Computer Science template before submission (OPEN gate). Do not infer author, affiliation, funding, competing-interest, archive DOI, or cover-letter content from another paper.

## Rebuild

```bash
python3 src/v2/build_peerj_package.py
```

The builder fails on any `undefined citation`, `undefined reference`, or `Label(s) may have changed` in the LaTeX logs.

## Local checks

```bash
cd flat_upload
latexmk -pdf -interaction=nonstopmode manuscript.tex
sha256sum -c ../SHA256SUMS.txt
```

For current visual review, inspect only `flat_upload/manuscript.pdf` and its separately uploaded Figure PDFs. Historical ZIPs and release capsules are retained as provenance and are not current review targets.
