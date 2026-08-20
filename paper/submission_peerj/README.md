# scReg-Eval figure package

Figure PDFs, TikZ fragments, and table fragments generated from the canonical
`paper/figs/` tree. Full manuscripts are local-only and are not on the public
repository HEAD.

## Layout

- `source/`: self-contained editable TikZ source, including the thirteen main-text
  figure fragments, three appendix fragments under `figs_extension/`, and eight tables.
- `flat_upload/`: figure PDFs (`Figure1.pdf` through `Figure13.pdf`,
  `FigureA1.pdf`–`FigureA3.pdf`), table fragments, references, and checksums.
  Printed Figure 1 is study design (`Figure1.pdf`); printed Figure 13 is coverage QC
  (`Figure13.pdf`).
- `supplemental/`: AI-in-code Supplemental Data S1–S3. Rebuild with
  `python3 src/v2/build_ai_disclosure_zips.py`.
- `SHA256SUMS.txt`: hashes for files in `flat_upload/`.
- `supplemental/SHA256SUMS_AI.txt`: hashes for Supplemental Data S1–S3.

## Rebuild figures

```bash
cd paper
python3 make_panel_data.py   # if panel_data.json inputs changed
Rscript make_figs.R
```

See `paper/CANONICAL_BUILD.md`.
