# AI-in-code disclosure — scReg-Eval

This package answers the PeerJ Computer Science form item on use of AI in
computer code. It is a **limited, representative disclosure**, not a full
development history.

## What is disclosed

Three analysis scripts received limited engineering edits with a large language
model during the 2026-07 revision cycle:

1. `pbmc_uce_eval_v2.py` — co-expression normalization contract (CP10k then log1p).
2. `pair_probe_stats.py` — deterministic per-family seed contract (schema v2).
3. `fixed_panel_audit.py` — layout-independent project-root resolution.

No AI system generated scientific data, randomization outputs, or reported
numbers. Results are recomputed by the deterministic pipelines.

## Three supplemental archives (upload these)

PeerJ CS asks for original code, the prompts used, and the resulting code as
separate compressed files. Build them with:

```bash
python3 src/v2/build_ai_disclosure_zips.py
```

The archives are written to `paper/submission_peerj/supplemental/`,
**not** inside `flat_upload/` or `upload.zip`:

| PeerJ SI name | Archive |
| --- | --- |
| Supplemental Data S1 | `Supplemental_Data_S1_original_code.zip` |
| Supplemental Data S2 | `Supplemental_Data_S2_prompts.zip` |
| Supplemental Data S3 | `Supplemental_Data_S3_resulting_code.zip` |

Do not unpack these into GitHub Pages. The site assembler copies only
`paper/submission_peerj/flat_upload/manuscript.pdf` and the Frontiers PDF.

## Working tree

- `before/` — pre-revision versions from git history (commits in `MANIFEST.json`).
- `after/` — the same three scripts after the disclosed edits. Production code
  remains `src/v2/*.py`.
- `prompts/00_INDEX.md` — condensed user instructions mapped to those edits.
  This is a reconstructed log, not a dump of every assistant session.

Full current code: https://github.com/PeterPonyu/scfm-reg-audit
Archive: https://doi.org/10.5281/zenodo.21724336
