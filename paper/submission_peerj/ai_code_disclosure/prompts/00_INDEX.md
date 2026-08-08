# Prompt log — AI-assisted code changes (session-log reconstruction)

This index maps each AI-assisted change in the disclosure package to the user instruction that
drove it. Prompts are reconstructed faithfully from the session records of the 2026-07 revision
cycle; wording is condensed, not verbatim.

## 01 — PBMC UCE co-expression normalization repair

**Files:** `pbmc_uce_eval_v2.py` (before: commit 81f2e23; after: current release)

**User instruction (2026-07-31):** "An independent review found the PBMC UCE co-expression
control skips the library-size normalization used by every other row. Fix the control to the
canonical CP10k→log1p contract, keep the verified UCE embedding, version the cache, add tests,
and rerun only what is affected."

**AI-assisted change:** added `normalized_log_counts()` (library-size CP10k then log1p), a
`co_normalization_version` cache field with provenance rejection of stale caches, a legacy-cache
migration that reuses the verified UCE embedding and rebuilds only the co-expression graph, and
schema-2 result metadata. Follow-up reruns changed PBMC UCE full-spec ρ from 0.005880 to
0.005927 with all other pooled rows bit-identical.

## 02 — Deterministic probe seeds

**Files:** `pair_probe_stats.py` (before: commit 6ac7252; after: current release)

**User instruction (2026-07-31):** "The probe's Mantel null seeds use Python's salted hash(), so
q-values are not reproducible across processes. Replace them with explicit deterministic seeds
and record the seed contract in the output."

**AI-assisted change:** replaced `abs(hash(fam))` seed derivation with
`seed_root * 1000 + sorted_family_index`, recorded `mantel_seed` per family and a
`seed_contract` string, bumped the stats schema to v2, and preserved the existing deterministic
sign-flip stream (`seed_root + 1`). Point estimates and sign-flip q-values were unchanged; Mantel
q-values shifted only within Monte Carlo tolerance.

## 03 — Layout-independent project root

**Files:** `fixed_panel_audit.py` (before: commit 2c2876b; after: current release)

**User instruction (2026-08-01):** "The release capsule's bundled tests fail because the module
computes the project root as two directories up from src/v2, which is wrong inside the capsule
layout. Make root resolution work from both layouts without hardcoding."

**AI-assisted change:** added `_project_root()`, which honors `SCREG_PROJECT_ROOT` and otherwise
walks upward from the module until it finds `data/manifest/shared_genes.v2.json` (dev layout
`src/v2/` and capsule layout `src/`). The release builder now runs the bundled test suite after
every rebuild and fails on error.

---

Full release and validator: <https://github.com/PeterPonyu/scfm-reg-audit> ·
<https://doi.org/10.5281/zenodo.21724336>
