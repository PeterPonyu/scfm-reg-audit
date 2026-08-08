# scReg-Eval PeerJ submission package

Generated from the canonical manuscript by `src/v2/build_peerj_package.py`.
The package contains 13 pooled FM/readout rows and eleven figures. The builder first
regenerates `paper/panel_data.json` with `make_panel_data.py` and all canonical figure
and table fragments with `Rscript make_figs.R`, then rebuilds both source and flat mirrors.
The current human-review PDF is `flat_upload/manuscript.pdf`; it is rebuilt with figures
pre-rendered as standalone PDFs.

## Layout

- `source/`: self-contained editable article source, including the eleven figure fragments and two tables. Its PDF and LaTeX intermediates are generated build outputs, not the human-review target.
- `flat_upload/`: upload-ready package with `manuscript.pdf`, `Figure1.pdf` through `Figure11.pdf`, two table fragments, references, and checksums.
- `internal/figure_build/`: temporary standalone wrappers used to export the figure PDFs; generated only.
- `HUMAN_GATES.md`: remaining author decisions and external actions (DONE/OPEN status).
- `PEERJ_FORM_TEMPLATE.md`: copy-ready scientific fields; repository DOI filled.
- `SHA256SUMS.txt`: hashes for all files in `flat_upload/`.

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
