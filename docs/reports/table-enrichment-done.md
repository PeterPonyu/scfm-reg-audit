# Table enrichment — completed (manual after workflow fail)

Workflow `enrich-tables` paused at Verify because the Implement agent reported success without writing files. Work was extracted from the locked SPEC (`table_enrichment_spec.json`) and applied directly.

## Delivered

### Table 1 — `figs/table1_primary_fixed_panel.tex`
- Added **Support**: `both` / `both (neg)` / `neither` / `--` (baseline)
- Six `both` positives + brain Geneformer attention `both (neg)` match dual $q_M,q_D<0.05$

### Table 2 — `figs/table2_cross_tissue_observed.tex`
- Columns: Obs. ρ | Add. frac. | Resid. ρ | φ | Edge Jac. | Mean row ρ | Shared TFs
- Numbers match decomp + pair JSONs (e.g. Brain–PBMC 0.5249 / 0.7543 / 0.5073 / 0.5070 / 0.416 / 0.4001 / 446)

### Table 3 — `figs/table3_pertype_ranges.tex` (new)
- Brain/PBMC × full/non-degree: n, ρ min/max/median from audit `descriptive_summary`

### Manuscript
- Captions for `tab:cross`, `tab:primary`, `tab:pertype_ranges`
- In-text pointers in proxy and per-type subsections

### Build
- `Rscript make_figs.R` OK
- `manuscript.pdf` rebuild OK (19 pages)
- `build_peerj_package.py` after this note

## Generator
`paper/make_figs.R` table block (end of file) — not hand-edited table `.tex`.
