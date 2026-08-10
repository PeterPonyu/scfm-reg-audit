# AI-in-Code Disclosure Package — scReg-Eval

This package accompanies the manuscript's Declaration of generative AI use and answers the
PeerJ Computer Science form item *Use of Artificial Intelligence (AI) in Computer Code*.

## Scope

This is a **representative, curated disclosure**, not a full development history. It covers the
parts of the analysis code where AI tools (Anthropic Claude models accessed through an
agentic coding assistant) made material engineering contributions during the 2026-07
revision cycle:

1. `pbmc_uce_eval_v2.py` — PBMC UCE graph generation. AI-assisted repair of the co-expression
   normalization contract (library-size CP10k followed by log1p), cache provenance/versioning,
   and a legacy-cache migration path.
2. `pair_probe_stats.py` — TF-disjoint probe statistics. AI-assisted replacement of salted
   `hash()` seeds with an explicit, documented per-family seed contract (schema v2).
3. `fixed_panel_audit.py` — core fixed-panel statistics module. AI-assisted project-root
   resolution (works from both the development and the release-capsule layouts) and hardening of
   null-semantics documentation.

No AI system generated scientific data, randomization outputs, or any reported number. All
results are recomputed by the deterministic pipelines and verified by the test suite
and the released artifact validator.

## Contents

- `before/` — the actual pre-revision versions of the three scripts, taken from the project's
  own version-control history (commit hashes recorded in `MANIFEST.json`). Personal machine
  paths are replaced by `${SCREG_DATA_ROOT}` / `${SCREG_PROJECT_ROOT}` placeholders.
- `after/` — tip copies of the same three scripts from `src/v2/` at package rebuild time
  (`allow_pickle=False` on NPZ loads). **Do not treat `after/` as an independent SoT**; the
  production loaders live in the repository at `src/v2/*.py`.
- `prompts/` — an ordered prompt log reconstructed from the session records, mapping each
  AI-assisted change to the files it touched.
- `MANIFEST.json` — file inventory, roles, source commits, and SHA-256.

## How to read the delta

Each `before/after` pair is a real, inspectable git diff of one file. The prompt log in
`prompts/00_INDEX.md` lists, for each pair, the user instruction and the resulting change at a
high level. The full current codebase is public at <https://github.com/PeterPonyu/scfm-reg-audit>
and archived at <https://doi.org/10.5281/zenodo.21724336>.
