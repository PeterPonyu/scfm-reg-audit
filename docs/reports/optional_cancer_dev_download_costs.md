# Optional cancer/dev download cost model（估计 / 不下载体）

**Date:** 2026-08-11  
**Scope:** Order-of-magnitude only. **No downloads in this session.**  
**Hardware assumption:** local RTX **5090 24 GB**; remote GPU usually **N/A** if local is used.  
**Chinese-friendly summary at bottom.**

## Scenarios

### S1 — Tiny construct SI（DESCARTES spleen RDS）

| Resource | Estimate | Notes |
|----------|----------|-------|
| Network | **~0.08–0.15 GB** | ~100 MB compressed RDS if not already local |
| Disk | **~0.2 GB** | RDS + derived peaks/cache |
| 5090 hours | **~0.5–2 h** | motif scan + `G_ATAC` + Mantel vs locked proxies (CPU-heavy motifs; GPU lightly used) |
| Remote GPU | **N/A** | Local sufficient |
| PeerJ risk | Low if SI-only under `results/v2/extension/` | |

### S2 — One multiome tissue（BMMC；已在本地）

| Resource | Estimate | Notes |
|----------|----------|-------|
| Network | **0 GB** | Already on disk (~5.7 GB) |
| Disk | **~6–12 GB** | h5ad + extension NPZ/graphs |
| 5090 hours — construct only | **~1–4 h** | proxy/`G_ATAC` + Mantel/decomp; no FM Support |
| 5090 hours — full FM audit | **~20–80 h** | rough; 5–8 readouts × nulls × BH; **gated / not this wave** |
| Remote GPU | **N/A** (local) or **same wall-time** if offloaded | Bandwidth to move 6 GB+ caches may dominate |
| PeerJ risk | **High** if Support/BH mutated | Stay construct-lane until re-SAP |

### S3a — Modest cancer multiome pilot（HTAN 单样本级）

| Resource | Estimate | Notes |
|----------|----------|-------|
| Network | **~0.5–5 GB** | One open Synapse/L3–L4 sample; not cohort dump |
| Disk | **~1–8 GB** | peaks + RNA partner if paired |
| 5090 hours | **~2–10 h** construct; **tens of hours** if FM | Synapse auth / dbGaP may block |
| Remote GPU | Usually **N/A** | |
| PeerJ risk | Extension-only | |

### S3b — Whole cancer/dev lake（对照：不要做）

| Resource | Estimate | Notes |
|----------|----------|-------|
| Network | **10s–100s GB** | DESCARTES RAW / File_S6 / HTAN cohort |
| Disk | **same order** | Desktop lakes already ~90 GB RNA-only |
| 5090 hours | **weeks** if naively FM’d | Estimand mismatch for RNA lakes |
| Remote GPU | Costly + transfer | |
| Policy | **Forbidden** for Support / `G_ATAC` | RNA∩ATAC empty |

## Rough compute decomposition（5090 24 GB）

| Step | VRAM | Time scale (one tissue) |
|------|------|-------------------------|
| Motif scan / peak×TF | low–mid (often CPU) | 0.5–3 h |
| `G_ATAC` assemble + Mantel | low | 0.2–1 h |
| Coexp / simple baselines | low | <1 h |
| Single FM readout embed | ~8–22 GB | 1–8 h depending on model |
| Full dual-null FM suite | peaks near 24 GB | many hours–days |

## Decision heuristic

1. **Need SI construct numbers soon?** → S1 (~0.1 GB network) after PeerJ freeze.  
2. **Need 3rd-tissue same-estimand story?** → S2 compute-only + panel policy; still not PeerJ Support.  
3. **Cancer narrative?** → S3a one sample max; never S3b into Support.

---

## 中文摘要（给决策用）

- **小下载（脾脏 DESCARTES RDS）**：流量约 **0.1 GB**，本地 5090 大约 **0.5–2 小时**做 construct；远程 GPU 一般不需要。  
- **BMMC**：已经在硬盘上，**再下载 = 0**；只算算力。Construct 约 **1–4 小时**；完整 FM 审计大约 **几十小时**量级，且会动 SAP/BH，**当前论文冻结期内不做**。  
- **HTAN 级小试点**：流量大约 **0.5–5 GB**；整湖/整队列是 **几十到上百 GB**，且 RNA lake **不能**进 Support/`G_ATAC`。  
- **远程 GPU**：本机 5090 够用时通常 **N/A**；若上云，传输大 h5ad/缓存往往比算力更贵。
