# Visual verdict — extension SI figures (machine patrol)

**Date:** 2026-08-12  
**Scope:** New SI figures only (`paper/figs_extension/`).  
**Existing PeerJ fig1–12:** **no panel additions / no edits this wave** (git: only `figs_extension/*` added).  
**Method:** standalone pdflatex → PNG @150dpi → multimodal visual QA (`/visual-verdict` contract).  
**Style reference cousins:** PeerJ fig2 / fig11 grammar (4-panel patchwork, 11pt, BLUE/AQUA/YELLOW).

## Scope matrix

| Layer | Objects | Status |
|-------|---------|--------|
| L1 New figures | E1, E2, E3 (`fig_ext1–3`) | **Reviewed** |
| L2 Panels added onto old figures | — | **N/A** (none) |
| L3 Old figures unchanged | fig1–12 under `paper/figs/` | **Not in wave** (isolation check pass) |

## Per-figure verdicts (JSON contract)

### E1 `fig_ext1_construct_mantel`

```json
{
  "score": 92,
  "verdict": "pass",
  "category_match": true,
  "differences": [
    "Panel D title was truncated in first render; shortened to 'Linked peaks (PeerJ ref ~7--12k)'",
    "Heatmaps use discrete identity fills (no colorbar) vs PeerJ continuous legends"
  ],
  "suggestions": [
    "Optional: add a small discrete legend key for rho bands if human prefers",
    "Human final check of heatmap cell contrast on print PDF"
  ],
  "reasoning": "4-panel layout readable; no raster path watermarks; numbers match claim-pack; BMMC-PBMC 0.89 clear."
}
```

### E2 `fig_ext2_baselines_collectri`

```json
{
  "score": 93,
  "verdict": "pass",
  "category_match": true,
  "differences": [
    "Freeze panel mixes small integers (0/3/13) — intentional badge, not multi-scale science plot",
    "Status bars all 'ready' (limited visual variance)"
  ],
  "suggestions": [
    "Human may prefer freeze as table not bars — optional post-final polish"
  ],
  "reasoning": "Titles no longer collide; CollecTRI coverage and edge counts legible; freeze message readable."
}
```

### E3 `fig_ext3_honesty_policy`

```json
{
  "score": 92,
  "verdict": "pass",
  "category_match": true,
  "differences": [
    "C gene label sits near 0.90 gate line (readable but tight)",
    "C/D titles close at center gutter (no hard clip after shorten)"
  ],
  "suggestions": [
    "Human may want gate annotation only in title (already primarily title-driven)"
  ],
  "reasoning": "HTAN blocked / orphan empty-obs / BMMC coverage fail / Support=13 freeze all visible and correct."
}
```

## Aggregate machine gate

| Metric | Result |
|--------|--------|
| Min score | **92** (≥90 threshold) |
| Fail count | **0** |
| Machine patrol | **PASS** |
| Human final review | **READY** |

## Artifacts for human终审

| PNG | Path |
|-----|------|
| E1 | `paper/figs_extension/visual_qa/fig_ext1_construct_mantel-1.png` |
| E2 | `paper/figs_extension/visual_qa/fig_ext2_baselines_collectri-1.png` |
| E3 | `paper/figs_extension/visual_qa/fig_ext3_honesty_policy-1.png` |
| TeX | `paper/figs_extension/fig_ext{1,2,3}_*.tex` |
| Preview | `paper/figs_ext_preview.tex` |

## Isolation re-check

- `paper/figs/fig_ext*`: absent  
- `manuscript.tex` inputs of `figs_extension`: none  
- PeerJ package rebuild: not run  

**Handoff:** Machine visual-verdict **PASS**. Please do human visual final review on the three PNGs above.
