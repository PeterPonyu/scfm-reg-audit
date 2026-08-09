# Table enrichment plan (workflow: `enrich-tables`)

**Status:** Workflow launched. Pauses after Spec for human approval of columns, then Implement → Verify.

**Scope:** Enrich manuscript **tables** from already-computed `results/v2` stats (figures already polished). No new experiments. No retired bootstrap/MDE/power.

## Current tables

| Table | File | Columns now |
|-------|------|-------------|
| Table 1 primary | `figs/table1_primary_fixed_panel.tex` | Tissue, Readout, Partial ρ, p_M, q_M, p_D, q_D |
| Table 2 cross-tissue | `figs/table2_cross_tissue_observed.tex` | Tissue pair, Observed Spearman ρ, Shared TFs |

Generator: end of `paper/make_figs.R` (comment was “tables (unchanged)”).

## Locked enrichment targets

### Table 2 — cross-tissue (descriptive)

| Column | Source |
|--------|--------|
| Tissue pair | display labels |
| Observed Spearman ρ | decomp / pairwise |
| Additive fraction explained | `cross_tissue_additive_decomp_v2.json` |
| Residual Spearman | decomp |
| Binary support φ | decomp |
| Edge Jaccard | `cross_tissue_{brain_vs_pbmc,atac,pbmc_vs_fibroblast}.json` |
| Mean per-TF-row Spearman | same pairwise files |
| Shared TFs | 446 |

Example magnitudes (for verify): Brain–PBMC obs 0.5249, add_frac 0.75, φ 0.51, jacc 0.42, row_ρ 0.40.

### Table 1 — primary full-confound

Keep existing p/q columns; **add Support**:

- `both` if q_M < 0.05 and q_D < 0.05 (and note neg if ρ < 0)
- `neither` otherwise
- baseline: outside FM families (e.g. `baseline` or em-dash + caption)

### Table 3 (new, preferred) — per-type ranges

From `fixed_panel_audit_v2.json` → `per_cell_type.*.descriptive_summary`:

- Tissue × confound_spec (full / non_degree)
- n_rows, ρ min / max / median

Descriptive only; no p/q.

## Forbidden

- `power_analysis_v2`, `cross_tissue_bootstrap_v2` CIs, `stats_enhanced_v2`
- Per-type bootstrap CIs from `pertype_stats_enhanced_v2` (retired)

## Workflow phases

1. **Inventory** (read-only) — map gaps  
2. **Spec** (read-only) — exact numbers → `await_user` gate  
3. **Implement** (R + tex) — edit `make_figs.R` / `manuscript.tex`, `Rscript make_figs.R`, latexmk  
4. **Verify** — cell-level match to JSON + PDF build  

## How to operate

- Dashboard: `/workflows` (display name `enrich-tables`)
- After reviewing Spec: `/workflow resume enrich-tables`
- Definition: `.grok/workflows/enrich-tables.rhai`
