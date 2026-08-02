# Full-pipeline run log — 2026-08-02

One complete external rerun of the scReg-Eval statistics and probe pipeline
(layer 2 of `FULL_RERUN.md`: statistics + probe, with the FM graphs treated as
pinned inputs). Closes the round-2 review SHOULD item "publish one external
full-pipeline run log (environment + input/output hashes)".

- Worktree: clean `git worktree` at commit `9864a6a` (v0.2.9-era code; the
  statistics/probe stages are untouched by the later v0.3.0 figure commit)
- Wall clock: 11:58:09 → 12:41:52 EDT (43m 43s), all 9 stages exit 0
- Untracked runtime inputs symlinked in: `data/annotation`, `data/multiome`,
  `data/uce`, `data/genome`, `data/motifs`, `data/scfoundation`, and the 31
  pinned `results/v2/*.npz` caches

## Environment

- Host: Linux 7.0.0-28-generic x86_64; Intel Core Ultra 9 275HX;
  NVIDIA RTX 5090 Laptop 24 GB (driver 595.84) — unused by these CPU stages
- Python 3.13.5; numpy 2.2.6; scipy 1.16.3; anndata 0.12.10; torch 2.12.0+cu130;
  scikit-learn 1.8.0; pandas 2.3.3; pyfaidx 0.9.0.4
- R 4.3.3; ggplot2 4.0.3; dplyr 1.2.1; patchwork 1.3.2; jsonlite 2.0.0;
  tikzDevice 0.12.6 (figures not regenerated in this run)

Matches the locked environment in `FULL_RERUN.md` §1 exactly.

## Input hashes (SHA-256, verified this run)

All 15 dataset/model files re-hashed on 2026-08-02 and match the tables in
`FULL_RERUN.md` §2–§3 exactly (GSE174367, AD RNA, PBMC RNA/ATAC/raw bundle,
GSE206767; Geneformer V2-104M safetensors + 3 dicts; scGPT best_model.pt +
vocab.json; scFoundation model.pt; UCE 4layer_model.torch + human_esm2.pt).
The 31 pinned `results/v2/*.npz` FM/proxy graph caches were hashed as run
inputs; their values are archived alongside this log in the repository run
artifacts (contact the authors if byte-level cache verification is needed —
the caches are regenerable per `FULL_RERUN.md` §4).

## Stage log

| stage | start (EDT) | exit | duration |
|---|---|---|---|
| run_fixed_panel_audit | 11:58:09 | 0 | 26m 29s |
| brain_coexp_baseline_null | 12:24:38 | 0 | 4m 58s |
| pbmc_coexp_baseline_null | 12:29:36 | 0 | 4m 26s |
| subdivide_injection | 12:34:02 | 0 | 0m 59s |
| tf_disjoint_split | 12:35:01 | 0 | 1s |
| build_pair_features | 12:35:02 | 0 | 2s |
| run_pair_probe | 12:35:04 | 0 | 12s |
| pair_probe_stats | 12:35:16 | 0 | 6m 30s |
| pair_probe_sensitivity | 12:41:46 | 0 | 6s |

## Output comparison (rerun vs released `results/v2/*.json`)

Field-level comparison over every shared key (floats compared exactly):

| file | result |
|---|---|
| fixed_panel_audit_v2.json | 4,567 shared keys; 11 diffs, ALL embedded absolute paths in provenance fields (worktree vs repo location); zero numeric diffs |
| inference_status_v2.json | 182 shared keys; 2 diffs: one embedded path inside a verification command string, one cascading file hash caused by the path-embedded JSON above; zero numeric diffs |
| brain_coexp_baseline_null_v2.json | identical (13 keys) |
| pbmc_coexp_baseline_null_v2.json | identical (13 keys) |
| injection_subdivided_v2.json | 1,145 shared keys; 1 numeric diff — see finding below |
| tf_probe_pair_eval_v2.json | identical (248 keys) |
| tf_probe_pair_stats_v2.json | identical (91 keys) |
| tf_probe_pair_sensitivity_v2.json | identical (68 keys) |

**Finding (fixed):** the released `injection_subdivided_v2.json` carried one
stale copied field — `observed_effects[12]` (PBMC UCE encoder) read 0.00588,
the pre-renormalization value, while the authoritative
`fixed_panel_audit_v2.json` (and the manuscript) report 0.005927 after the
cp10k/log1p co-expression-control repair. The rerun regenerated 0.005927,
confirming the released file was simply never refreshed after that repair.
The field is provenance-only: no figure, table, or manuscript number reads it
(figure ladder curves use the `rows[*].replicate_runs` subdivision data, which
matched exactly). Corrected to 0.005927 on 2026-08-02.

## Conclusion

The pipeline reproduces every published statistic from the pinned inputs.
The only divergences are embedded absolute paths (an artifact of running in a
temporary worktree) and the single stale provenance field documented above.
