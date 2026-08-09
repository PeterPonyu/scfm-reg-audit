# Canonical manuscript build (do not use stale mirrors)

## Source of truth (edit these only)

| Asset | Path |
| --- | --- |
| Figure R code | `paper/make_figs.R` |
| TikZ fragments (generated) | `paper/figs/fig*.tex` |
| Manuscript | `paper/manuscript.tex` |
| Bibliography | `paper/references.bib` |
| Class | `paper/wlpeerj.cls` |
| **Human review PDF** | `paper/manuscript.pdf` |

## Rebuild (always this order)

```bash
cd paper
python3 make_panel_data.py   # if panel_data.json inputs changed
Rscript make_figs.R          # regenerates paper/figs/*.tex from R — never hand-edit those
cd ..
python3 src/v2/build_peerj_package.py   # syncs submission_peerj/ + refreshes manuscript.pdf
```

Or, for paper-only after R:

```bash
cd paper && latexmk -pdf -interaction=nonstopmode manuscript.tex
```

## Generated mirrors (do not edit; may lag until package rebuild)

| Mirror | Role |
| --- | --- |
| `paper/submission_peerj/source/` | Editable TikZ source package |
| `paper/submission_peerj/flat_upload/` | Upload: `manuscript.pdf` + `Figure1.pdf`…`Figure11.pdf` |
| `paper/submission_peerj/internal/figure_build/` | Standalone figure wrappers |
| `paper/submission_peerj/upload.zip` | Zip of flat_upload checksum set |
| `paper/submission_peerj/SHA256SUMS.txt` | Hashes for flat_upload |

If `flat_upload/manuscript.pdf` or `FigureN.pdf` is older than `paper/figs/` or `make_figs.R`, it is **stale**. Rebuild with `build_peerj_package.py`.

## Not manuscript SoT (safe to ignore / delete)

| Path | Notes |
| --- | --- |
| `paper/figs/*preview*.png` | Ad-hoc visual checks |
| `paper/Rplots.pdf`, `Rplots.pdf` | R accidental device output |
| `.omx/artifacts/` | Local visual-QA scratch |
| `paper/submission_peerj/_stale_archive/` | Parked outdated docx etc. |

## Panel D legend

Legend placement for Fig 8 panel D is controlled only in `make_figs.R` (`f7d`, `legend.position = "top"`). After editing R, you must re-run `Rscript make_figs.R` and the PeerJ package builder so **every** PDF copy updates.
