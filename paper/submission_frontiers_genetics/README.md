# Two report paths

This repository produces **two independent PDFs**. Do not mix their class files or build commands.

## 1. Original PhD / PeerJ CS report (scReg-Eval)

Source of truth: `paper/manuscript.tex` (PeerJ class). **Do not edit that file for the Frontiers draft.**

```bash
cd paper
latexmk -pdf -interaction=nonstopmode manuscript.tex
```

Output: `paper/manuscript.pdf`

Figures 1–12 and Tables 1–8 live under `paper/figs/` and are `\input` with a `figs/` prefix. Rebuild TikZ with `Rscript paper/make_figs.R` only when panel JSON or `make_figs.R` changed. Do **not** run `src/v2/build_peerj_package.py` for the Frontiers draft.

## 2. Frontiers in Genetics report (full conversion + appendix)

This folder is a **complete conversion** of the original paper into a default-LaTeX `article` draft for *Frontiers in Genetics* (Original Research / Computational Genomics), **plus** three construct-lane ATAC overlays as an **Appendix** of the same unified `.tex`.

It is **not** a three-figure-only paper. Main text keeps the original science (Figures 1–13, Tables 1–8, Methods, Results, Discussion), with a study-design schematic as Figure 1 and coverage QC as Figure 13. Appendix Figures A1–A3 are additional spleen / BMMC / Treg overlays.

```bash
cd paper/submission_frontiers_genetics
latexmk -pdf -interaction=nonstopmode manuscript.tex
```

Output: `paper/submission_frontiers_genetics/manuscript.pdf`

| Asset | Path |
| --- | --- |
| Unified manuscript | `manuscript.tex` (this folder) |
| Review PDF | `manuscript.pdf` |
| Bibliography (copy) | `references.bib` |
| Main TikZ (copies of `paper/figs/`) | `fig_study_design.tex`, `fig1_truth_construct.tex` … `fig12_protocol_pass_matrix.tex` |
| Tables (copies) | `table1_*.tex` … `table8_*.tex` |
| Appendix TikZ (copies of `paper/figs_extension/`) | `fig_ext1_construct_mantel.tex`, `fig_ext2_baselines_collectri.tex`, `fig_ext3_honesty_policy.tex` |
| Original PeerJ source (do not edit here) | `paper/manuscript.tex` |

This folder is **flat on purpose** (portal compilers cannot follow nested trees). Do not add `figures/` or `supplementary/` subdirectories. All `\input{...}` paths are basename-only.

Class: `\documentclass[11pt,a4paper]{article}` — not a journal class file. Single spacing, page numbers, and line numbers are on for review.

## What authors should edit

All author metadata is at the top of `manuscript.tex` in the **EDITABLE AUTHOR BLOCK**:

- `\AuthorName`
- `\AuthorAffiliation`
- `\AuthorEmail`

Those macros feed the title-page author line, the correspondence email, and the Author Contributions sentence. Add `\AuthorNameB` / `\AuthorAffiliationB` (and thread them into `\author{...}`) if co-authors join.

Also replace the Funding placeholder (`[Placeholder: replace with grant numbers if applicable.]`) with grant numbers when available.

## Regenerating TikZ copies

Main-text figures/tables (only if `paper/make_figs.R` or panel JSON changed):

```bash
# from repo root
Rscript paper/make_figs.R
cp -f paper/figs/fig_study_design.tex \
      paper/figs/fig{1,2,3,4,5,6,7,8,9,10,11,12}_*.tex \
      paper/figs/table{1,2,3,4,5,6,7,8}_*.tex \
      paper/submission_frontiers_genetics/
```

Appendix figures (only if claim-pack JSON or `paper/make_figs_extension.R` changed):

```bash
SCFM_BASE=. Rscript paper/make_figs_extension.R
cp -f paper/figs_extension/fig_ext1_construct_mantel.tex \
      paper/figs_extension/fig_ext2_baselines_collectri.tex \
      paper/figs_extension/fig_ext3_honesty_policy.tex \
      paper/submission_frontiers_genetics/
```

Then recompile this folder with `latexmk` as above.

## What not to do

- Do not rewrite `paper/manuscript.tex` for this draft.
- Do not compile separate supplementary `.tex` files.
- Do not `\input{../figs/...}` or `\input{../figs_extension/...}`; copies are local.
- Do not nest `figures/` or `supplementary/`.
- Do not run `src/v2/build_peerj_package.py` for this draft.
- Do not treat Appendix A1–A3 as the whole paper; they are extra construct-lane material.
