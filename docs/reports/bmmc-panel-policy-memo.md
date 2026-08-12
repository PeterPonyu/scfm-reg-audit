# BMMC panel-policy memo (no FM audit)

**Date:** 2026-08-12 (decision refresh; draft 2026-08-11)  
**Status:** **written decision — P3 stay**  
**Constraint:** Do **not** start BMMC full FM audit into PeerJ Support (PRD Option C gated).

## Decision

**Stay P3:** BMMC (`GSE194122`) remains `role=construct_candidate`, `lane=construct` only.

| Field | Value |
|-------|--------|
| **panel_policy** | `frozen_brain_pinned_446x1200_disclosed_coverage` |
| **Panel** | Frozen brain-pinned **446 TF × 1200 genes** (no re-freeze) |
| **Allowed** | Extension construct SI only: prepare peaks → overlay `G_ATAC` under `results/v2/extension/construct/GSE194122/` → Mantel/decomp vs locked proxies |
| **Forbidden** | FM Support rows; PeerJ public JSON edits; re-SAP; Option C full FM audit |

**Rationale:** Preflight gene/TF coverage on the frozen panel fails the ≥90% replication gate (see Coverage). Re-freezing (P2) would change estimand identity and is out of scope for the current PeerJ package. Promoting BMMC into Support (Option C / full FM) would expand BH families 8→12 and N≠13, forcing re-SAP. Construct-lane SI on the frozen panel with **disclosed incomplete coverage** preserves estimand identity for transfer diagnostics without inflating Support.

**Written decision (one line):** P3 stay — construct SI only on frozen 446×1200 with disclosed coverage; no FM Support / no re-SAP / no Option C.

### Explicit non-goals

- **Option C** full FM audit into PeerJ Support — **gated** indefinitely under this decision (requires separate re-SAP + panel-policy revision).
- No Support-row inflation in this revision.
- No merge of Cancer/Dev RNA lakes into `G_ATAC`.
- No MANIFEST primary lock edits for BMMC.
- No PeerJ public JSON edits for BMMC.

## Coverage (disclosed)

Source: `enhancement_v1` postflight / monorepo  
`results/enhancement_v1/gse194122_preflight.json`  
(`post_download_inspection.frozen_panel_coverage`; also cited in tissues notes).

| Axis | Overlap / panel | Coverage | ≥90% gate |
|------|-----------------|----------|-----------|
| Genes | **1076 / 1200** | **0.8967** (≈0.897) | fail |
| TFs | **341 / 446** | **0.7646** (≈0.765) | fail |

Gate status: `blocked` (`frozen_gene_coverage_below_90_percent`, `frozen_tf_coverage_below_90_percent`).

Overlay build (G006) used the full frozen panel dimensions: `n_tf=446`, `n_genes=1200`  
(`results/v2/extension/construct/GSE194122/G_ATAC_v2_GSE194122_meta.json`).  
Coverage above is the **GEX/feature overlap preflight** disclosure, not a resize of the panel contract.

## Construct progress (G006)

Executed under extension construct (not PeerJ Support):

| Artifact | Path |
|----------|------|
| Overlay G_ATAC | `results/v2/extension/construct/GSE194122/G_ATAC_v2_GSE194122.npz` |
| Meta | `…/G_ATAC_v2_GSE194122_meta.json` |
| Mantel vs locked | `…/mantel_vs_locked.json` |
| Additive decomp row | `…/additive_decomp_row.json` |
| Summary | `…/construct_summary.json` |

Mantel pairs (all `status=ok`, `n_tf_common=446`):

| Pair | observed Spearman | frac. explained by additive marginals |
|------|-------------------|----------------------------------------|
| GSE194122 vs GSE174367 (brain) | 0.529 | 0.608 |
| GSE194122 vs PBMC10k | 0.890 | 0.460 |
| GSE194122 vs GSE206767 (fibroblast) | 0.473 | 0.689 |

`construct_summary.json`: `peerj_support_rows_unchanged=13`, `peerj_decomp_rows_unchanged=3`.

## Facts

- BMMC NeurIPS multiome (`GSE194122`) is already local (~5.7 GB h5ad).
- Frozen panel is brain-pinned **446 TF × 1200 genes**.
- PeerJ freeze: 13 full-spec FM rows, 8 BH families; adding BMMC FM implies BH **8 → 12** and N≠13 (re-SAP).

## Policy options (historical menu)

| Option | Action | Estimand identity | PeerJ impact |
|--------|--------|-------------------|--------------|
| **P1 — Filter to frozen panel** | Keep 446×1200; disclose coverage; drop uncovered edges | Same estimand; coverage incomplete | Extension paper / post-submit; still re-SAP if Support grows |
| **P2 — Re-freeze panel** | Rebuild gene/TF panel for multi-tissue coverage | New estimand identity | Not for current PeerJ package |
| **P3 — Defer FM / construct SI** | Keep BMMC as `lane=construct` only; frozen panel + disclosed coverage | Unchanged | **Chosen** — no Support / no re-SAP |

P3 operationalizes the frozen-panel + disclosed-coverage *panel* choice without P1’s Support path: SI under `results/v2/extension/` only.

## Registry

- `src/v2/extension/configs/tissues.json` / `tissues.yaml`: `bmmc.panel_policy = frozen_brain_pinned_446x1200_disclosed_coverage`
- Fallback tissue if remapping ever required later: **GSE200046** BM/CD34 multiome (GEO fetch) — not activated by this decision.
