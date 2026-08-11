# BMMC panel-policy memo (no FM audit)

**Date:** 2026-08-11  
**Status:** planning / extension infra only  
**Constraint:** Do **not** start BMMC full FM audit into PeerJ Support (PRD Option C gated).

## Facts

- BMMC NeurIPS multiome (`GSE194122`) is already local (~5.7 GB h5ad).
- Frozen panel is brain-pinned **446 TF × 1200 genes**.
- Preflight (enhancement_v1): gene coverage ≈ **0.897**, TF coverage ≈ **0.765** → blocked for naive remap.
- PeerJ freeze: 13 full-spec FM rows, 8 BH families; adding BMMC FM implies BH **8 → 12** and N≠13 (re-SAP).

## Policy options (choose before any FM suite)

| Option | Action | Estimand identity | PeerJ impact |
|--------|--------|-------------------|--------------|
| **P1 — Filter to frozen panel** | Keep 446×1200; disclose coverage; drop uncovered edges | Same estimand; coverage incomplete | Extension paper / post-submit; still re-SAP if Support grows |
| **P2 — Re-freeze panel** | Rebuild gene/TF panel for multi-tissue coverage | New estimand identity | Not for current PeerJ package |
| **P3 — Defer** | Keep BMMC as `lane=construct` dry-run only | Unchanged | **Default now** |

## Recommendation

Stay on **P3** until PeerJ `0x30`–`0x32` closeout. Infra registers `bmmc` as construct dry-run under `configs/extension/tissues.yaml`. Fallback tissue if remapping intractable later: **GSE200046** BM/CD34 multiome (GEO fetch).

## Explicit non-goals

- No Support-row inflation in this revision.
- No merge of Cancer/Dev RNA lakes into `G_ATAC`.
- No MANIFEST primary lock edits for BMMC.
