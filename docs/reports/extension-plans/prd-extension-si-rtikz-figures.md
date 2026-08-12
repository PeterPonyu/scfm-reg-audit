# RALPLAN: Extension SI figures via R + TikZ (post-PeerJ)

**Plan ID:** `prd-extension-si-rtikz-figures`  
**Status:** `executed` — completed in commit d8288bf (2026-08-12T11:24:24Z); E1–E3 R+TikZ delivered  
**Date:** 2026-08-12  
**Mode:** RALPLAN-DR SHORT · Architect ITERATE → plan revised · Critic re-check against this file  
**Hard UI:** **No Python plotting.** Path = **R (ggplot2/dplyr/patchwork) → tikzDevice → TeX → LaTeX** (`figure_helpers.R` / `make_figs.R` contract).  
**Hard science:** **No network download**; PeerJ Support 13 / primary `G_ATAC_*` / MANIFEST freeze.

---

## Architect + Critic closures (this revision)

| Item | Resolution |
|------|------------|
| P0 BMMC coverage SoT | **In-repo** `docs/reports/extension-claim-pack/honesty_policy.json` (memo constants). **Never** `results/enhancement_v1/gse194122_preflight.json` (absent). |
| P0 output isolation | Emit to **`paper/figs_extension/`** only — not `paper/figs/`. Custom emit or `figs_dir` override; do not use stock emit that hardcodes `paper/figs`. |
| P1 no MVP intermediate panel-data | **Banned for MVP.** No `make_panel_data_extension.R`, no `extension_panel_data.json`. |
| P1 E3 R-emitted only | ggplot/tikzDevice **or** R `writeLines` TeX — no hand-maintained schematic SoT. |
| Package freeze | Do **not** run `build_peerj_package.py`; FIGURE_MAP stays 12; no `\input{figs_extension` in `manuscript.tex` this wave. |

---

## 1. Requirements

Deep plan for SI visualization expansion after local-run extension wave:

1. Build SI figures from current extension data  
2. Data interfaces (schemas, paths, validation)  
3. Reasonableness + dock to frozen PeerJ figs/tables  
4. Gaps that need CPU/GPU before drawing is honest  

**Non-goals:** Python plots; D4/D5 fetch; rewrite PeerJ figs 1–12 / Support; default HTAN MACS.

---

## 2. RALPLAN-DR

### Principles (5)

1. **R+TikZ only** for figures — same 11pt/newtx/6.8in as PeerJ.  
2. **Extension SI additive** — `fig_ext*` under `figs_extension/`; not PeerJ Support.  
3. **Claim-pack JSON = statistical SoT** for Mantel numbers; honesty JSON for policy numbers.  
4. **Zero network**; heavy compute only if user later approves optional E4/G-CPU.  
5. **Disclosure over polish** (orphan filename meta, HTAN blocked, BMMC gate fail).

### Drivers (3)

1. Typography/pipeline continuity with PeerJ TikZ  
2. Traceability panel → JSON path + freeze flags  
3. Minimize compute while not inventing statistics  

### Options

| Opt | Description | Verdict |
|-----|-------------|---------|
| **A** | `paper/make_figs_extension.R` → **`paper/figs_extension/fig_ext*.tex`** + `figs_ext_preview.tex` | **Recommended** |
| **B** | Fold fig13+ into `make_figs.R` | Reject for freeze blur |
| **C** | Tables-only writeLines | Reject as sole path (OK for E3 subpanels) |

---

## 3. Data already on disk (R can draw MVP)

| Asset | Path | Use |
|-------|------|-----|
| Construct SI Mantel | `docs/reports/extension-claim-pack/construct_si_mantel.json` | E1: **9 rows** (3×3), ρ + additive + residual, n_tf=446 |
| Twin Mantel | `results/v2/extension/construct/{TAG}/mantel_vs_locked.json` | G-VAL-1 optional |
| G_ATAC meta | `.../G_ATAC_v2_{TAG}_meta.json` | E1D relevant_peaks |
| CollecTRI | `results/v2/extension/baselines/collectri_prior/summary.json` | E2A |
| Other baselines | `.../baselines/*/summary.json` | E2B status |
| HTAN blocked | `.../HTAN_GBM_C3N01334/htan_prepare_status.json` | E3A |
| Orphan pilot status | `.../GSE211155_treg/orphan_prepare_status.json` | E3B pilot fields |
| PeerJ docking RO | `paper/panel_data.json$third_tissue`, `results/v2/cross_tissue_additive_decomp_v2.json` | reference ticks only |
| Honesty (to create at execute) | `docs/reports/extension-claim-pack/honesty_policy.json` | E3 BMMC + lake counts + freeze |

### Mantel numbers (claim-pack)

| tissue | vs brain | vs PBMC | vs fibro | add. frac |
|--------|--------:|--------:|---------:|-----------|
| DESCARTES_spleen | 0.434 | 0.391 | 0.493 | 0.82–0.85 |
| BMMC GSE194122 | 0.529 | **0.890** | 0.473 | 0.46–0.69 |
| GSE211155_treg | 0.464 | 0.506 | 0.490 | 0.74–0.81 |

**Reasonableness:** same panel 446 TFs; BMMC–PBMC high but immune multiome-plausible; orphan single-type mid-ρ; HTAN no ρ row.

### Honesty JSON schema (execute step 1 — non-plot)

```json
{
  "schema_version": 1,
  "peerj_support_rows_touched": false,
  "bmmc": {
    "panel_policy": "frozen_brain_pinned_446x1200_disclosed_coverage",
    "gene_overlap": 1076, "gene_panel": 1200, "gene_coverage": 0.8967,
    "tf_overlap": 341, "tf_panel": 446, "tf_coverage": 0.7646,
    "gate_threshold": 0.9, "gate_status": "blocked",
    "source": "docs/reports/bmmc-panel-policy-memo.md"
  },
  "orphan_lake": {
    "n_h5ad": 115, "n_with_celltype_obs": 0,
    "pilot_gsm": "GSM6449881", "pilot_tag": "GSE211155_treg",
    "meta_source": "filename_sorted_population",
    "source": "docs/reports/extension-plans/construct/orphan_atac_inventory.md"
  },
  "htan": {
    "status": "blocked",
    "block_reason": "fragments_only_no_peak_matrix",
    "accepted_wave_completion": "E4.3b",
    "source": "results/v2/extension/construct/HTAN_GBM_C3N01334/htan_prepare_status.json"
  },
  "freeze": { "peerj_support_rows": 13, "primary_g_atac_mutable": false }
}
```

---

## 4. Figure set (R → TikZ)

Colors/theme from `make_figs.R` (BLUE/AQUA/YELLOW/VIOLET/…). Width 6.8in, 11pt.

### E1 `fig_ext1_construct_mantel.tex` (core)

- **A** Heatmap ρ (3 tissues × 3 locked proxies)  
- **B** Heatmap additive fraction  
- **C** Forest/dumbbell ρ; highlight BMMC–PBMC  
- **D** relevant_peaks bars (extension meta) + PeerJ third_tissue **reference** ticks (caption: not same total_peaks)

### E2 `fig_ext2_baselines_collectri.tex`

- **A** CollecTRI TF/gene coverage + edge keep rate from summary.json  
- **B** Baseline status strip (motif/degree/encode/collectri)  
- **C** *post-MVP only* (G-CPU-1 NPZ stats)

### E3 `fig_ext3_honesty_policy.tex` (required)

- **A** HTAN dual-path blocked (status JSON)  
- **B** Orphan lake + pilot filename meta (honesty + orphan_prepare)  
- **C** BMMC coverage vs 0.9 gate (honesty JSON only)  
- **D** Freeze badge Support=13 / peerj_support_rows_touched=false  

All R-emitted.

### E4 optional later

Multi-type orphan / HTAN external-peak Mantel — needs G-CPU-2/3.

---

## 5. R data interface (MVP)

```r
base <- Sys.getenv("SCFM_BASE", "..")
source(file.path(base, "src", "v2", "figure_helpers.R"))
setup_tikz_options()
J <- function(rel) jsonlite::fromJSON(file.path(base, rel), simplifyVector = FALSE)

si <- J("docs/reports/extension-claim-pack/construct_si_mantel.json")
stopifnot(isFALSE(isTRUE(si$peerj_support_rows_touched)))
stopifnot(length(si$rows) == 9L, as.integer(si$n_pair_rows_ok) >= 9L)
rows <- dplyr::bind_rows(lapply(si$rows, as.data.frame, stringsAsFactors = FALSE))

hon <- J("docs/reports/extension-claim-pack/honesty_policy.json")
stopifnot(isFALSE(isTRUE(hon$peerj_support_rows_touched)))

htan <- J("results/v2/extension/construct/HTAN_GBM_C3N01334/htan_prepare_status.json")
stopifnot(identical(htan$status, "blocked"))

# PeerJ RO docking only
panel <- J("paper/panel_data.json")
```

**Write:**

```
Rscript paper/make_figs_extension.R
→ paper/figs_extension/fig_ext{1,2,3}_*.tex
→ paper/figs_ext_preview.tex
# NEVER: paper/figs/fig_ext*
# NEVER: build_peerj_package.py in same chain
```

Emit helper must target `figs_extension/` (fork `emit` if `emit_figure()` hardcodes `paper/figs`).

---

## 6. Gaps / compute

| ID | Need | MVP? | Cost |
|----|------|------|------|
| **G-HON-1** | Write honesty_policy.json | **Required** | seconds, 0 net |
| G-CPU-1 | CollecTRI NPZ SI stats | Optional E2C | CPU minutes |
| G-CPU-2 | Multi-type orphan | E4 | CPU 0.5–2h |
| G-CPU-3 | HTAN fragments×external peaks | E4 | CPU hours |
| G-VAL-1 | Twin mantel assert | Optional | seconds |
| G-PKG-1 | Keep out of PeerJ zip | Process | 0 |

**MVP: no GPU, no heavy CPU.**

---

## 7. External peaks note (context only — not MVP)

Local peak universes exist (brain/fibro/BMMC/DESCARTES/cCRE). HTAN tar is fragments-only. Happy path = local count (0 traffic, CPU hours) or new download (forbidden). MVP uses **blocked honesty panel**.

---

## 8. Implementation steps (after user 批准执行 only)

1. Write `honesty_policy.json`  
2. Scaffold `make_figs_extension.R` + `figs_extension/` + custom emit  
3. E1 → E2 A/B → E3  
4. Preview latexmk only  
5. Verify isolation + freeze  
6. Update SI inventory map  

Staff: `$team` (honesty+E3 ∥ E1 ∥ E2) or `$ralph`.

---

## 9. Acceptance criteria

1. SI figures **only** R→tikzDevice and/or R writeLines TeX  
2. No Python plotting modules for this wave  
3. E1 ρ match claim-pack to 4 decimals; 9 rows  
4. E3C = honesty 0.8967 / 0.7646 with memo source  
5. Captions: orphan filename meta; HTAN blocked; BMMC P3  
6. R asserts freeze flags; HTAN blocked  
7. Zero network; Support/primary G_ATAC untouched  
8. `paper/figs_extension/fig_ext*.tex` exist; **`paper/figs/fig_ext*` must not**  
9. FIGURE_MAP length 12; no manuscript input of fig_ext; no package rebuild in SI chain  

---

## 10. Verification (after execute)

```bash
test -f docs/reports/extension-claim-pack/honesty_policy.json
cd paper && Rscript make_figs_extension.R
test -f figs_extension/fig_ext1_construct_mantel.tex
! ls figs/fig_ext*.tex 2>/dev/null
grep -E 'blocked|0\.8967|filename|Support' figs_extension/fig_ext3_honesty_policy.tex | head
cd .. && python -m pytest src/v2/tests/test_extension.py -q
python validate_artifacts.py
git status --short results/v2/G_ATAC* results/v2/*.public.json paper/submission_peerj/
# must NOT: python3 src/v2/build_peerj_package.py
```

---

## 11. Risks

| Risk | Mitigation |
|------|------------|
| tikz clip/metrics | figure_helpers 11pt newtx |
| SI = Support confusion | figs_extension + E3 badge |
| Invented coverage | honesty_policy.json only |
| Over-claim orphan | single-type caption |
| ρ=0.89 over-read | construct transfer language |
| Package mutation | forbid builder in SI chain |

---

## 12. ADR

- **Decision:** Option A — isolated R+TikZ SI module under `figs_extension/`; E1–E3 from claim-pack + honesty_policy + status JSON; no Python plots; no heavy compute for MVP.  
- **Drivers:** freeze; visual continuity; 9 Mantel pairs already on disk.  
- **Rejected:** B fig13+; C tables-only; missing preflight path.  
- **Consequences:** SI preview independent of PeerJ zip; honesty JSON is policy SoT.  
- **Follow-ups:** Critic APPROVE → user 批准执行 → implement.

---

## 13. Out of scope

Python plots · D4/D5 · MACS HTAN · FM Support inflation · hand-edit TikZ · package FIGURE_MAP expansion · MVP intermediate panel-data R script.

---

## 14. Status gate

**`pending approval`**

Do **not** implement until:

1. Critic **APPROVE** on **this** plan file, and  
2. User **「批准执行」** / approve via team or ralph.

**RALPLAN consensus state:** Architect ITERATE closed by this revision; Critic must re-read **this** `plan.md` (not stale docs copy).

---

## 15. Consensus record

| Role | Verdict | Notes |
|------|---------|-------|
| Planner | Draft Option A | R+TikZ SI, claim-pack SoT |
| Architect | ITERATE then closed | honesty JSON; figs_extension; no MVP panel-data |
| Critic | **APPROVE** | B1–B4 closed on this plan.md |

**RALPLAN gate:** complete → plan **`pending approval`** for **user execution approval only**.  
**Do not implement** until user: **「批准执行」** / approve via team|ralph.


---

## 16. Execution record

| Field | Value |
|-------|-------|
| Executed | 2026-08-12T11:24:24Z |
| Commit | `d8288bf` |
| Driver | `paper/make_figs_extension.R` |
| Outputs | `paper/figs_extension/fig_ext{1,2,3}_*.tex` |
| Verify | pytest 37/37; validate_artifacts PASS; no paper/figs leak |
| Ultragoal | G016–G018 complete (SI R+TikZ wave) |
