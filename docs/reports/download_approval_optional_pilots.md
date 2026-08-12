# Download approval document — optional cancer/dev / construct pilots

**Date:** 2026-08-11 (updated 2026-08-12: D4/D5 → `pending_large`)  
**Status:** Approval / policy + **fail-closed download infrastructure**  
**Hard stop default:** Do not curl/wget/aria2 D4/D5 until verbal approve.  
**Related:** `optional_cancer_dev_download_costs.md`, `bmmc-panel-policy-memo.md`,  
`${DESKTOP_DATA}/datasets/extension_pilots/manifests/FETCH_LOG.md`

**Phase-2 seam (fail-closed stub):** `src/v2/extension/download_gate.py` +
`cli.py download --plan-id …` refuse without `SCREG_DOWNLOAD_APPROVED=1` and a
matching `SCREG_DOWNLOAD_PLAN_ID`. Even when those env gates pass, the CLI **does
not fetch** — it only prints a manual recipe / writes a dry-run plan under
`results/v2/extension/download-plans/`.

**Separate execute path (still gated):**  
`src/v2/extension/scripts/fetch_approved_plan.py` (wrapper:
`scripts/fetch_plan.sh`) can perform HTTP fetch **only** with env match **and**
`--execute`. Do **not** run `--execute` for D4/D5 until a human says so.
`fetch_optional_pilots.sh` remains a demoted pointer (exit **2**).

**Ceremony vs capability:** `SCREG_DOWNLOAD_*` is an operator checklist ceremony
(filled approval row + matching plan id), **not** a security/auth boundary. Exit
**0** from the CLI means recipe/dry-run emitted only — never “download succeeded.”

**Data-root fallback:** Extension path resolution (`src/v2/extension/paths.py`)
uses `DESKTOP_DATA` → `SCREG_DATA_ROOT` → `~/Desktop/data`. Core v2 tooling may
still fall back to repo `data/` in some paths — do not assume the same default.

---

## Why D4/D5 were previously called “cannot download”

This was **project policy**, not a technical impossibility:

| Plan | Prior wording | Actual reason blocked | Technical readiness |
|------|---------------|-----------------------|---------------------|
| **D4** | Rejected / refuse lake | Whole-lake HTAN/DESCARTES (e.g. `GSE149683_RAW` ~29 GB, File_S6 ~4.3 GB) exceeds pilot hard gates; PeerJ freeze + disk risk | URLs + dest paths now in `PLAN_REGISTRY`; fetch script ready |
| **D5** | Rejected / refuse | Cancer28/Dev27 RNA lakes (~90 GB class) already largely local; **forbidden** for Support / `G_ATAC` (RNA∩ATAC empty); re-fetch pointless for FM Support | Local inventory recipe ready; still identity-forbidden for Support |

**Now:** D4/D5 are **`pending_large` / approval-gated pending infrastructure** (same
ceremony as D1/D3), not “impossible forever.”

---

## Dual gates (policy + technical)

Every large plan must clear **both**:

1. **Policy risk gate** — PeerJ freeze / disk / identity (esp. D5 never → Support/`G_ATAC`).  
2. **Technical readiness gate** — exact URLs (or Synapse shortlist), Desktop destinations, size estimates, checksum notes, copy-paste commands.

Verbal one-liner required before any D4/D5 `--execute`, e.g. **「批准 D4」** / **「批准 D5」** / **「批准 D4+D5」**.

---

## Purpose

Record **what may be fetched later**, under what size/policy gates, and who must approve.
Executable extension code prefers **local assets**. Large lakes wait for explicit
human go-ahead after infrastructure is ready (this document + registry + scripts).

---

## Approval matrix

| ID | Asset | Network | Disk | May feed `G_ATAC` / Support? | Approver gate | Decision now |
|----|-------|--------:|-----:|------------------------------|---------------|--------------|
| D0 | Nothing (local-only) | 0 | 0 | N/A — use existing fibro/brain/PBMC/`G_ATAC` | — | **Default** |
| D1 | DESCARTES spleen RDS (tiny construct SI) | ~0.1 GB | ~0.2 GB | Construct SI under `results/v2/extension/` only | Explicit user go-ahead + URL | **Pending** (infra ready; may already be local) |
| D2 | BMMC multiome | 0 (already local) | ~6 GB | Construct OK; **full FM Support blocked** until panel policy | Panel memo (P1/P2/P3) | **No new download** |
| D3 | HTAN open single-sample pilot | ~0.5–5 GB | ~1–8 GB | Extension-only; never inflate PeerJ 13-row SAP | Explicit go-ahead + Synapse/open check | **Pending** |
| D4 | Whole HTAN / DESCARTES RAW / File_S6 lakes | ~33 GB GEO (+ TB HTAN) | same | **Forbidden** for Support / locked `G_ATAC`; lake_blocked until subset | Verbal 「批准 D4」 + env + disk check | **`pending_large`** |
| D5 | Cancer 28 / Dev 27 RNA lakes | 0 (local inventory) | ~90 GB | **Forbidden** for Support / `G_ATAC` even after inventory | Verbal 「批准 D5」 + env | **`pending_large`** (inventory; no Support wiring) |

---

## Policy rules

1. **No multi-GB automated fetches** in agent sessions without a filled approval row above **and** a verbal approve for `pending_large`.  
2. **27/28 RNA lakes** must not enter Support or `G_ATAC` construction — **even after D5 approve**.  
3. **BMMC full FM audit** stays gated (`bmmc-panel-policy-memo.md`); construct lane + infra OK.  
4. PeerJ freeze: do not rewrite `results/*.public.json` contracts, MANIFEST primary locks, or 13-row SAP for new tissues.  
5. After any future approved fetch, raw data lives under Desktop `data/`; derived graphs only under `results/v2/extension/`.  
6. D4 HTAN **cohort-wide** remains Synapse/TB — not an auto-URL; shortlist open samples only.

---

## Cost summary (order-of-magnitude)

See full tables in `optional_cancer_dev_download_costs.md`.

| Scenario | Network | 5090 hours (construct) | Remote GPU |
|----------|--------:|------------------------:|------------|
| D1 tiny RDS | ~0.1 GB | ~0.5–2 h | N/A |
| D2 BMMC local | 0 | ~1–4 h construct; FM tens of hours (gated) | N/A |
| D3 HTAN sample | ~0.5–5 GB | ~2–10 h | Usually N/A |
| D4 RAW+S6 | ~33 GB | days+ if naively FM’d; prefer subset construct | Transfer-heavy |
| D5 RNA lakes | 0 (local) | N/A for `G_ATAC`; coexp-side only | N/A |

---

## Human checklist before any fetch

- [ ] Approval ID selected (D1/D3/D4/D5)  
- [ ] Exact URL / Synapse ID recorded (D4 GEO URLs are in `PLAN_REGISTRY`)  
- [ ] Expected compressed size < gate (D1: 2 GB; D3: &lt;10 GB; D4: expect ~33 GB GEO)  
- [ ] Free disk checked (D4: ≥ ~40 GB recommended)  
- [ ] Confirm write path is Desktop data, **not** PeerJ submission package  
- [ ] Confirm construct code will consume local path via env (`ATAC_FILE` / `SCREG_*`)  
- [ ] **D5:** confirm no Support / `G_ATAC` wiring after inventory  

---

## Copy-paste commands (after verbal approve)

### Dry-run recipe (no network) — any plan

```bash
export SCREG_DOWNLOAD_APPROVED=1
export SCREG_DOWNLOAD_PLAN_ID=D4   # or D5 / D1 / D3
python src/v2/extension/cli.py download --plan-id "$SCREG_DOWNLOAD_PLAN_ID"
# writes results/v2/extension/download-plans/<ID>.dry_run.json
```

### Real fetch / inventory (still requires env + `--execute`)

```bash
# D4 — GEO RAW (~29 GB) + File_S6 (~4.3 GB). HTAN cohort skipped (Synapse manual).
export SCREG_DOWNLOAD_APPROVED=1 SCREG_DOWNLOAD_PLAN_ID=D4
python src/v2/extension/scripts/fetch_approved_plan.py --plan-id D4            # dry-run
python src/v2/extension/scripts/fetch_approved_plan.py --plan-id D4 --execute  # ONLY after 「批准 D4」

# Optional single asset:
python src/v2/extension/scripts/fetch_approved_plan.py --plan-id D4 --execute --asset-id GSE149683_File_S6

# D5 — local inventory only (no GEO re-fetch by default; still forbidden for Support/G_ATAC)
export SCREG_DOWNLOAD_APPROVED=1 SCREG_DOWNLOAD_PLAN_ID=D5
python src/v2/extension/scripts/fetch_approved_plan.py --plan-id D5 --execute
```

Destinations (D4):

- `${DESKTOP_DATA}/datasets/extension_pilots/descartes_lake/GSE149683_RAW.tar`
- `${DESKTOP_DATA}/datasets/extension_pilots/descartes_lake/GSE149683_File_S6.Cicero_gene_activity_scores_by_cell_type.csv.gz`

D5 paths (already local):

- `${DESKTOP_DATA}/datasets/CancerDatasets{,2}/`
- `${DESKTOP_DATA}/datasets/DevelopmentDatasets{,2}/`

---

## Post-approve local steps (no CLI fetch)

### D1 — DESCARTES spleen (after human go-ahead + manual fetch)

1. Place RDS (and/or converted h5ad) under  
   `${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/`.  
   Preferred peaks file: `descartes_spleen_peaks.h5ad` with `var_names` as `chrom:start-end`.
2. Validate bridge (fail-closed if absent; **never downloads**):  
   `python src/v2/extension/cli.py descartes-bridge`
3. Build overlay `G_ATAC` (honors `SCREG_EXTENSION_OUT`; PeerJ lock on):  
   `TAG=DESCARTES_spleen ATAC_FILE=<h5ad> SCREG_EXTENSION_OUT=results/v2/extension/construct/DESCARTES_spleen SCREG_PEERJ_SUPPORT_LOCK=1 python src/v2/build_atac_graph_v2.py`
4. Mantel/decomp:  
   `python src/v2/extension/cli.py construct --tissue descartes_spleen --execute --write`

### D2 — BMMC (already local; P3 = no FM Support)

1. Confirm h5ad under `${DESKTOP_DATA}/external/scfm-reg-audit/gse194122/` (or `SCREG_BMMC_H5AD`).
2. Extract/rename ATAC peaks into the extension overlay:  
   `python src/v2/extension/cli.py prepare-bmmc --execute --write`
3. Run the `build_command` printed by prepare-bmmc / construct dry-run  
   (`SCREG_EXTENSION_OUT=results/v2/extension/construct/GSE194122`, `SCREG_PEERJ_SUPPORT_LOCK=1`).  
   Do **not** overwrite locked `results/v2/G_ATAC_v2_*.npz`.
4. Mantel/decomp:  
   `python src/v2/extension/cli.py construct --tissue bmmc --execute --write`  
   Full FM Support remains gated (`bmmc-panel-policy-memo.md` **P3**).

### D4 / D5 — after 「批准 D*」 only

- Use `fetch_approved_plan.py` as above.  
- D4: keep under `extension_pilots/descartes_lake/`; registry role `lake_blocked` until a future subset construct job.  
- D5: inventory only; **never** feed Support / `G_ATAC`.

---

## 中文摘要（审批用）

- **过去说 D4/D5「不能下」= 政策闸**（整湖体量 / PeerJ 冻结 / RNA 不能进 Support），**不是技术做不到**。  
- **现在：** D4/D5 = `pending_large`，与 D1/D3 一样要审批 + env；CLI 仍只出配方；真下载走 `fetch_approved_plan.py --execute`。  
- **默认可做：** 本地 fibro / brain / PBMC 的 extension construct / baseline / claim-pack。  
- **BMMC：** 已在硬盘，不要再下；完整 FM 进 Support 仍被面板政策卡住。  
- **可后续批准的小下载：** D1 脾脏 RDS；D3 单个 HTAN 开放样本。  
- **可后续批准的大湖（需口头批准）：** D4 GEO RAW+S6；D5 仅盘点本地 RNA 湖，且**永远不得**进 Support/`G_ATAC`。  

---

## Sign-off

| Role | Name | Date | Decision |
|------|------|------|----------|
| Owner | | | ☐ D0 only / ☐ Approve D1 / ☐ Approve D3 / ☐ Approve D4 / ☐ Approve D5 / ☐ Approve D4+D5 / ☐ Other |
| Notes | | | |
