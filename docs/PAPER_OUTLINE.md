# Paper claim register — fixed-panel audit

Updated 2026-07-31 after verification of the enlarged audit and TF-disjoint probe.

## Current status (2026-07-31)

The current manuscript reports 13 pooled model/readout rows, six figures, and the supervised
TF-disjoint probe. The probe paired sign-flip q-values are UCE vs co-expression .0325, random-init
floor .0300, scGPT .255, Geneformer embedding .5125, and Geneformer attention .65. The UCE contrast
is not model-specific because the random-init floor is also supported. Do not state that every FM is
indistinguishable from co-expression, and do not claim that the probe establishes regulatory recovery.

## Historical snapshot (2026-07-24)

The following claim register was the pre-expansion snapshot. It remains here for audit history; the
current manuscript and JSON coverage take precedence where counts or status differ.

## Authority

Only these files support inferential manuscript claims (hashes pinned 2026-07-31 after the UCE
normalization repair and deterministic-seed reruns):

- `results/v2/fixed_panel_audit_v2.json`
  - SHA-256: `442e0d877f5665ed6345106f46e9397e2a8fd448d85c5307e2b98ec5db57774a`
- `results/v2/fixed_panel_signal_injection_v2.json`
  - SHA-256: `96129465691d5e34f310fed4ac04fb06c15e0f92144ace3ba9a6ef5de980253e`
- `results/v2/inference_status_v2.json`
  - SHA-256: `1c1c1ad77157460fa6aafb18b152eb473d46724ce84eb4b2d1fbc4fa55300405`
- `results/v2/brain_coexp_baseline_null_v2.json` (explicit integer seeds)
  - SHA-256: `5122f75ffe96f7bdc03973a302ffd94e96750ac880539e6644adfad4ef5231f2`
- `results/v2/tf_probe_pair_stats_v2.json` (schema 2, deterministic family seeds)
  - SHA-256: `94690d2a9aa0ac7368307fcd0e00ad31101a9c9ade8e320af747100d3eef95ba`
- `results/v2/pbmc_uce_pooled_v2.json` (CP10k-log1p normalized co-expression control)
  - SHA-256: `80e5f4aa4352b1c8b95afc587c8f21da34a8f125ef1674e70d2f93651f838a5f`

The four legacy files pinned in `LEGACY_INFERENCE_NOTE.md` remain unchanged for audit, but their
bootstrap CI, MDE, coverage, exclusion, and TF-superpopulation interpretations are retired.

## Primary design

- Fixed panel: 446 unique TFs by 1,200 genes.
- Reference: accessibility/motif **regulatory-potential proxy**, not causal regulatory truth.
- Statistic: partial Spearman correlation after rank transformation.
- Full controls: co-expression, target peak count, gene length, detection rate, GC, TF out-degree,
  target in-degree.
- Nulls: gene-label permutation and within-TF-row target shuffle.
- Monte Carlo: `N=999`, two-sided plus-one correction, resolution 0.001.
- BH: eight families = tissue × confound specification × null type.
- Primary specification: full confounds. Non-degree is sensitivity only.

## Primary pooled claims

| Tissue | Readout | partial rho | Mantel p/q | row-shuffle p/q | Interpretation |
|---|---|---:|---:|---:|---|
| Brain | Geneformer embedding | 0.004420 | .003/.006 | .001/.002 | supported by both nulls |
| Brain | scFoundation encoder | 0.000619 | .633/.633 | .643/.643 | unsupported |
| Brain | UCE encoder | 0.004612 | .001/.003 | .001/.002 | supported by both nulls |
| Brain | scGPT encoder | 0.012407 | .001/.003 | .001/.002 | supported by both nulls |
| Brain | Geneformer KO | -0.000933 | .407/.4884 | .366/.4392 | unsupported |
| Brain | Geneformer KO control | -0.001887 | .100/.150 | .091/.1365 | unsupported |
| PBMC | Geneformer embedding | 0.004254 | .001/.003 | .002/.006 | supported by both nulls |
| PBMC | Geneformer attention | -0.000176 | .866/.866 | .839/.839 | unsupported |
| PBMC | scFoundation encoder | 0.001115 | .393/.5895 | .359/.5385 | unsupported |

Allowed summary: four weak, readout-specific primary alignments and five unsupported primary rows.
Do not summarize this as either “scFMs encode regulation” or “no scFM encodes regulation.”

## Secondary and descriptive claims

- Non-degree sensitivity: no row is supported by both nulls; several are supported only by the
  row-shuffle null. Report null dependence, not discoveries.
- Per-cell-type: 58 FM rows across two specifications; descriptive only, no p/q/BH or population CI.
- Cross-tissue proxy: observed Spearman 0.524858, 0.524496, and 0.462325; descriptive only.
- Injection: 11 alpha levels × 30 replicates × 2 tissues = 660 unique-seed runs. This is an
  axis-aligned pipeline-sensitivity diagnostic, not power, MDE, CI, coverage, or exclusion analysis.

## Forbidden claim patterns

- “regulatory truth,” “causal truth,” or a causal ceiling for the proxy;
- global positive or global negative claims about scFMs;
- TF-block bootstrap CI as current uncertainty;
- Mantel or row-shuffle null described as anti-conservative in general;
- MDE, implied alpha, power, coverage, or excluded effect range;
- per-cell-type significance or cross-tissue uncertainty;
- population generalization beyond the fixed panel and stated randomizations.

## Current manuscript assets

- Canonical source: `paper/manuscript.tex`
- Data-driven generator: `paper/make_figs.R`
- Current visuals: six figures and two tables from authoritative JSON.
- Build: `Rscript make_figs.R`, then `latexmk -g -pdf manuscript.tex` from `paper/`.
- Author identity, affiliation, funding, and CRediT declarations are now populated in the manuscript;
  the repository DOI remains gated on public release.
