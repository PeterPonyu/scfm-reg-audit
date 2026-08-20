# Canonical figure build (do not use stale mirrors)

Manuscript TeX/PDF sources are local-only and are not on the public repository HEAD.
The public tree keeps figure-generation code, TikZ fragments, and the bibliography.

## Source of truth on the public tree (edit these)

| Asset | Path |
| --- | --- |
| Figure R code | `paper/make_figs.R` |
| TikZ fragments (generated) | `paper/figs/fig*.tex` |
| Bibliography | `paper/references.bib` |
| Class | `paper/wlpeerj.cls` |

## Rebuild figures

```bash
cd paper
python3 make_panel_data.py   # if panel_data.json inputs changed
Rscript make_figs.R          # regenerates paper/figs/*.tex from R — never hand-edit those
```

Do not run `latexmk` in `paper/` on a public clone: full manuscripts are not shipped.

## Generated figure mirrors

| Mirror | Role |
| --- | --- |
| `paper/submission_peerj/source/` | Editable TikZ source package |
| `paper/submission_peerj/flat_upload/` | Figure PDFs + table fragments + checksums |

If `FigureN.pdf` is older than `paper/figs/` or `make_figs.R`, it is **stale**. Upload `Figure1.pdf` is study design; `Figure13.pdf` is coverage QC. Appendix figures are `FigureA1.pdf`–`FigureA3.pdf`. AI-in-code supplementals are three zip archives in `paper/submission_peerj/supplemental/` (Supplemental Data S1–S3), not inside `flat_upload/`.

## Not figure SoT (safe to ignore)

| Path | Notes |
| --- | --- |
| `paper/figs/*preview*.png` | Ad-hoc visual checks |
| `paper/Rplots.pdf`, `Rplots.pdf` | R accidental device output |
| `.omx/artifacts/` | Local visual-QA scratch |

## Panel D legend

Legend placement for Fig 8 panel D is controlled only in `make_figs.R` (`f7d`, `legend.position = "top"`). After editing R, re-run `Rscript make_figs.R` so every TikZ fragment updates.
