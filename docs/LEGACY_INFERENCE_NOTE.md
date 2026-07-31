# Legacy inference retirement note (scfm-reg-audit)

This note is the human-readable companion to `results/v2/inference_status_v2.json`.
It supersedes any claim in `paper/`, `docs/`, or earlier drafts that relied on the
bootstrap CI, MDE, superpopulation framing, or per-type global BH for the TF-block
confounded-controlled partial-Spearman statistic on the preregistered 446 TF × 1200
target panel.

## Status of legacy JSON outputs (kept verbatim, NOT authoritative)

| File | sha256 | Status |
|---|---|---|
| `results/v2/stats_enhanced_v2.json` | `bd84b5af0d81e74739495231ac8a5774f96197253d0efed69e374e50c948b39a` | retired, retained for audit |
| `results/v2/power_analysis_v2.json` | `ab6dd2e384ebc3244531cbedfec4a3b5074934881a624a354b3b42b7a52c0e9f` | retired, retained for audit |
| `results/v2/cross_tissue_bootstrap_v2.json` | `7f0ffa5a49196df2e843c3080a326b1edec54fbf63c4b5fea9bbc12e4dfbf750` | retired, retained for audit |
| `results/v2/pertype_stats_enhanced_v2.json` | `e50f4a552449c47c3a9e14e5787ac49ed981114c6d7f79f05c58b23db475ac86` | retired, retained for audit |

These hashes are pinned in `results/v2/inference_status_v2.json` and verified bit-for-bit
before any new authoritative output is generated. The files are intentionally not
modified or deleted; their byte-level integrity is required for the audit trail.

## Why they are retired

1. **Bootstrap CI was mis-calibrated.** The red-team review and an independent read
   of `src/v2/stats_enhanced_v2.py` confirm that the existing TF-block bootstrap draws
   one scalar index per sampled TF (`tf_pos[t] = np.where(tf_rows == t)[0]`, which
   returns a single-element array for each unique TF), then concatenates those scalar
   indices. The resulting `pos` array has length `n_tf = 446` per replicate and feeds
   `fm_v[pos]`, which has ~534K edges. The procedure effectively does NOT select whole
   TF clusters of edges, and the with-replacement draw produces duplicate scalar
   indices (~63% unique), which in turn produces tied ranks in `rankdata`. A percentile
   CI from a tied-rank bootstrap has no coverage guarantee, and the legacy outputs
   visibly violate it (e.g., observed 0.9055 vs CI [0.8505, 0.9023] in
   `power_analysis_v2.json` for brain α=1.0; observed 0.5245 vs CI [0.3084, 0.4869] in
   `cross_tissue_bootstrap_v2.json`).

2. **MDE / implied-alpha is a superpopulation claim.** The bootstrap CI is the
   load-bearing inference for the MDE curve and the FM-on-curve "implied alpha" mapping
   in `power_analysis_v2.json`. With the bootstrap broken, those numbers are also
   unsupported.

3. **Seeds are non-deterministic.** `stats_enhanced_v2.py` used
   `hash(tag + label) % 2**31` to seed the bootstrap RNG. Python's built-in `hash()`
   is salted per process (`PYTHONHASHSEED`), so two invocations of the same script on
   the same inputs produce different CI numerics.

4. **"Four models" was overclaimed.** The legacy outputs report 4 FM families in the
   pooled rows, but `results/v2/` does not contain per-cell-type UCE graphs, does not
   contain any PBMC-pooled scGPT graph, and does not contain any per-cell-type scGPT
   graphs. The new `model_scope_decision_v2.json` records the matrix explicitly.

## Per-type is DESCRIPTIVE EXPLORATORY ROBUSTNESS ONLY

`pertype_stats_enhanced_v2.json` previously carried per-cell-type global BH-adjusted
q-values across the full per-type family. **The new audit does NOT do that.** Per-type
rows in `fixed_panel_audit_v2.json` are DESCRIPTIVE: they report the observed confound-
controlled partial-Spearman rho per (tissue, cell_type, model_family, readout, spec)
along with a `descriptive_summary` block (range / median / mean / std / sign counts).
There is NO `p_mc`, NO `q`, NO `bh_q_family`, and NO BH adjustment across per-type rows.
This is a deliberate scope reduction to keep per-type wall time in the seconds range
and to drop the over-extended randomization framework that the red-team review
flagged.

## The 8 pooled BH family IDs

Pooled rows are organized into exactly EIGHT independent BH families, indexed by
`tissue ∈ {brain, pbmc}` × `confound_spec ∈ {full, non_degree}` × `null_type ∈
{mantel, degree}`. Each pooled FM row carries two `family_id` fields (`family_id_mantel`
and `family_id_degree`) and the matching `family_id` field inside the Mantel and
degree-preserving summaries.

| tissue | spec | null type | family_id |
|---|---|---|---|
| brain | full | mantel | `brain_pooled_full_confound_mantel` |
| brain | full | degree | `brain_pooled_full_confound_degree` |
| brain | non_degree | mantel | `brain_pooled_non_degree_confound_mantel` |
| brain | non_degree | degree | `brain_pooled_non_degree_confound_degree` |
| pbmc | full | mantel | `pbmc_pooled_full_confound_mantel` |
| pbmc | full | degree | `pbmc_pooled_full_confound_degree` |
| pbmc | non_degree | mantel | `pbmc_pooled_non_degree_confound_mantel` |
| pbmc | non_degree | degree | `pbmc_pooled_non_degree_confound_degree` |

BH-FDR is computed per family (independent) and stored as `bh_q_family` in each
summary block. **No BH is computed across these eight families** — each is its own
preregistered primary-or-sensitivity family.

## What replaces them

| Audience need | Replaced by |
|---|---|
| Pooled point effect + finite Monte-Carlo p_mc + BH | `results/v2/fixed_panel_audit_v2.json` `pooled.<tissue>.<primary_family>` (full) and `pooled.<tissue>.<sensitivity_family>` (non-degree). Eight family IDs above. |
| Per-cell-type exploratory robustness (effect sizes only) | `results/v2/fixed_panel_audit_v2.json` `per_cell_type.<tissue>.<full|non_degree_confound>` rows + `descriptive_summary` block. NO p_mc / q / BH. |
| Cross-tissue construct reproducibility | `results/v2/fixed_panel_audit_v2.json` `cross_tissue_construct_reproducibility` (observed Spearman only; legacy Mantel p_mc retained in `cross_tissue_bootstrap_v2.json` for cross-reference). |
| Pipeline sensitivity under injected signal | `results/v2/fixed_panel_signal_injection_v2.json` (axis_aligned_pipeline_sensitivity; per-replicate point effect only; no CI/MDE/power/exclusion). |
| Audit / status / hash pin | `results/v2/inference_status_v2.json` + `results/v2/model_scope_decision_v2.json`. |

## Forbidden vocabulary in new authoritative outputs

`bootstrap_ci95`, `clears_zero`, `implied_alpha`, `mde_*`, `coverage`, `power`,
`exclusion`, and any other superpopulation-inference language. New outputs use `p_mc`
(plus-one Monte-Carlo randomization p-value) with explicit seed, `N_perm`, and
resolution `1/(N_perm+1)`. They use `regulatory_potential_proxy` instead of "truth"
when referring to the motif/accessibility proxy. No `p < resolution` claim is made.

## Realism gate (production N=999)

The N=99 pooled benchmark (`benchmark_n99.py`) on the real cached graphs
(Ng=1200, ~534K edges) measures ~27.7s for brain pooled + ~22.7s for PBMC pooled
with shared-batch optimization. Scaling to N_PERM_POOLED_MANTEL=999 (resolution 0.001)
and N_PERM_POOLED_DEG=999 (resolution 0.001) projects ~280s (brain) + ~230s (pbmc)
for pooled; per-type descriptive adds ~14s (constant, no randomization). Total
production wall ≈ 8.7 min, well under the 2 h budget. Throughput is reported as
`wall / (N_perm × N_rows_in_shared_batch)` and shared-batch wall as `wall / N_perm`.

## Authoritative hashes at retirement

The hashes above are the byte-level state of the legacy files at the moment of this
retirement. They will be reverified before any future claim that depends on them.