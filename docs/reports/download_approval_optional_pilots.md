# Download approval document — optional cancer/dev / construct pilots

**Date:** 2026-08-11  
**Status:** Approval / policy only — **no automated download constructor**  
**Hard stop (2026-08-11):** User/leader revoked fetch execution — **D0 local-only**. Do not curl/wget/aria2.  
**Related:** `optional_cancer_dev_download_costs.md`, `bmmc-panel-policy-memo.md`

---

## Purpose

Record **what may be fetched later**, under what size/policy gates, and who must approve.
This wave builds **executable extension code on local assets**. Download tooling is
intentionally **not** a deliverable (`fetch_optional_pilots.sh` is demoted to a pointer).

---

## Approval matrix

| ID | Asset | Network | Disk | May feed `G_ATAC` / Support? | Approver gate | Decision now |
|----|-------|--------:|-----:|------------------------------|---------------|--------------|
| D0 | Nothing (local-only) | 0 | 0 | N/A — use existing fibro/brain/PBMC/`G_ATAC` | — | **Default** |
| D1 | DESCARTES spleen RDS (tiny construct SI) | ~0.1 GB | ~0.2 GB | Construct SI under `results/v2/extension/` only | Explicit user go-ahead + URL | **Pending** |
| D2 | BMMC multiome | 0 (already local) | ~6 GB | Construct OK; **full FM Support blocked** until panel policy | Panel memo (P1/P2/P3) | **No new download** |
| D3 | HTAN open single-sample pilot | ~0.5–5 GB | ~1–8 GB | Extension-only; never inflate PeerJ 13-row SAP | Explicit go-ahead + Synapse/open check | **Pending / not this wave** |
| D4 | Whole HTAN / DESCARTES RAW / File_S6 lakes | 10s–100s GB | same | **Forbidden** for Support / `G_ATAC` | Must refuse | **Rejected** |
| D5 | Cancer 28 / Dev 27 RNA lakes | 0 (local) | ~90 GB | **Forbidden** for Support / `G_ATAC` (RNA∩ATAC empty) | Must refuse | **Rejected** |

---

## Policy rules

1. **No multi-GB automated fetches** in agent sessions without a filled approval row above.  
2. **27/28 RNA lakes** must not enter Support or `G_ATAC` construction.  
3. **BMMC full FM audit** stays gated (`bmmc-panel-policy-memo.md`); construct lane + infra OK.  
4. PeerJ freeze: do not rewrite `results/*.public.json` contracts, MANIFEST primary locks, or 13-row SAP for new tissues.  
5. After any future approved fetch, raw data lives under Desktop `data/`; derived graphs only under `results/v2/extension/`.

---

## Cost summary (order-of-magnitude)

See full tables in `optional_cancer_dev_download_costs.md`.

| Scenario | Network | 5090 hours (construct) | Remote GPU |
|----------|--------:|------------------------:|------------|
| D1 tiny RDS | ~0.1 GB | ~0.5–2 h | N/A |
| D2 BMMC local | 0 | ~1–4 h construct; FM tens of hours (gated) | N/A |
| D3 HTAN sample | ~0.5–5 GB | ~2–10 h | Usually N/A |
| D4/D5 lakes | refuse | — | — |

---

## Human checklist before any fetch

- [ ] Approval ID (D1/D3) selected  
- [ ] Exact URL / Synapse ID recorded  
- [ ] Expected compressed size < gate (D1: 2 GB hard stop; D3: confirm <10 GB)  
- [ ] Confirm write path is Desktop data, **not** PeerJ submission package  
- [ ] Confirm construct code will consume local path via env (`ATAC_FILE` / `SCREG_*`)  

---

## 中文摘要（审批用）

- **本阶段不需要“下载构建代码”**；只要这份审批文档 + 成本估计。  
- **默认可做：** 用已有本地 fibro / brain / PBMC 的 `G_ATAC` 跑 extension construct / baseline / claim-pack。  
- **BMMC：** 已在硬盘，**不要再下**；完整 FM 进 Support 仍被面板政策卡住。  
- **可后续批准的小下载：** DESCARTES 脾脏 RDS（约 0.1 GB）或单个 HTAN 开放样本（约 0.5–5 GB）。  
- **明确拒绝：** 整湖 HTAN/DESCARTES RAW，以及 27/28 癌/发育 RNA lake 进 `G_ATAC`/Support。  

---

## Sign-off

| Role | Name | Date | Decision |
|------|------|------|----------|
| Owner | | | ☐ D0 only / ☐ Approve D1 / ☐ Approve D3 / ☐ Other |
| Notes | | | |
