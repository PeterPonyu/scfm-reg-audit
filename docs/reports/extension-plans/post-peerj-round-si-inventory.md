# Post-PeerJ extension round — SI inventory + figure map

**Date:** 2026-08-12  
**Stats path:** R → tikzDevice → `paper/figs_extension/` (not PeerJ FIGURE_MAP)  
**Network:** 0 · **Support rows touched:** false  

## Statistical figures (R + TikZ) — delivered

| ID | File | Content |
|----|------|---------|
| **E1** | `paper/figs_extension/fig_ext1_construct_mantel.tex` | A ρ heatmap · B additive frac · C forest · D relevant peaks |
| **E2** | `paper/figs_extension/fig_ext2_baselines_collectri.tex` | A CollecTRI coverage · B baseline status · C raw vs panel edges · D freeze counts |
| **E3** | `paper/figs_extension/fig_ext3_honesty_policy.tex` | A HTAN blocked · B orphan lake · C BMMC coverage 0.8967/0.7646 · D freeze badge |

**Build:** `SCFM_BASE=. Rscript paper/make_figs_extension.R`  
**Preview:** `paper/figs_ext_preview.tex` (latexmk; optional)  
**Data SoT:** `construct_si_mantel.json` (9 pairs) + `honesty_policy.json`  

## Construct Mantel numbers (9 pairs)

| tissue | vs brain | vs PBMC | vs fibro |
|--------|--------:|--------:|---------:|
| DESCARTES spleen | 0.434 | 0.391 | 0.493 |
| BMMC | 0.529 | **0.890** | 0.473 |
| Treg pilot (filename meta) | 0.464 | 0.506 | 0.490 |

## Physical-object note

Statistical SI panels use **R+TikZ only**. Python may be used later only for true physical/world schematics (e.g. peak–genome geometry), not for ρ/coverage bars.

## Out of scope this delivery

- PeerJ fig1–12 rewrite / package zip  
- HTAN fragments×external peaks compute  
- Multi-type orphan join  
- Python statistical plots  
