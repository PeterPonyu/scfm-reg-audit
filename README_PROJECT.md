# scfm-reg-audit — 单细胞 FM 的调控-vs-共表达 正交审计

**独立项目。与已投稿的 `research/sc-fm-benchmark`（IEEE JBHI/TCBB、PeerJ）完全分离。**

## 边界声明（避免与已投稿论文混淆）
- 本项目**不复用** `research/` 下任何结果、图、表、稿件；只共享 `~/Desktop/data`、`~/Desktop/labs` 里的**公开原始数据**。
- 已投稿论文的论点是"FM 相对基线的表现 + 评测惯例"。**本项目论点不同**：用 ATAC/motif 构成的 regulatory-potential proxy 审计"FM 学到的是调控信号还是共表达"。两者不得交叉引用彼此的未发表结论。
- 所有产物留在 `projects/scfm-reg-audit/` 内。

## 一句话
scReg-Eval：一个固定面板（446 TF × 1,200 基因）、混杂显式、null 语义明确的审计协议，回答"在受控条件下 FM 基因图与 accessibility/motif proxy 的一致性是否超出共表达基线"。当前交付是协议与审计结果（PeerJ CS 目标），不是因果真值判定或 recalibration 工具。

## 当前状态（2026-07-31）
- 13 个 pooled model/readout 行（brain 8 + PBMC 5），两个组织 × 两种 confound spec × 两种 null，全部入 `results/v2/fixed_panel_audit_v2.json`。
- 主稿 `paper/manuscript.tex`：6 图 2 表，13 页，含 TF-disjoint probe（Fig 6）。
- 模型范围：Geneformer（embed/attn/KO/floor）、scGPT、scFoundation、UCE。**CellPLM 已正式排除**：其 `OmicsEmbedder` 在 contextual 层前即聚合基因轴，没有可比的 per-gene contextual state（见 `docs/SCREG_EVAL_PROTOCOL.md` 架构适用性门槛）。
- PBMC UCE co-expression control 已修复为 CP10k→log1p（`co_normalization_version=cp10k_log1p_v1`）；brain baseline 与 probe stats 的种子已全部显式固定。
- 原始 DESIGN 中的 scReg-Fix、因果预训练臂、cancer/development case study 不在本篇范围内。

## 锁定的决策（2026-07-08，后经修订）
| 项 | 决定 |
|---|---|
| FM 阵容 | 4 个可比 family：Geneformer, scGPT, scFoundation, UCE（CellPLM 架构性排除） |
| 因果预训练臂 | 搁置（不在本篇范围） |
| 生物 case study | 搁置（不在本篇范围） |
| 算力 | 见 `compute/GPU_PLAN.md`；本篇全部结果 CPU 可复核，模型图缓存已固定 |

## 目录
- `docs/DESIGN.md` — 原始设计 spec（历史文档；当前范围见 SCREG_EVAL_PROTOCOL.md）
- `docs/SCREG_EVAL_PROTOCOL.md` — 当前协议定义与适用范围
- `docs/PAPER_OUTLINE.md` — 论文 claim register 与权威结果清单
- `compute/GPU_PLAN.md` — 历史算力配置记录
- `src/` `data/` `results/` — 代码 / 语料 manifest / 产物
- `LICENSE`（MIT，原创代码）、`LICENSE-CONTENT.md`（CC BY 4.0，稿件/图/派生结果）、`LICENSING.md`

## 数据来源（只读引用，不搬运）
- 公开数据：GSE174367（brain snATAC）、10x Genomics 10k PBMC multiome、GSE206767（fibroblast ATAC）
- scATAC 辅助：`~/Desktop/data/datasets/ATAC_data`（只读）

## Prior art（已核实真实，必须在 related work 引用并陈述 delta）
- Kendiukhov, arXiv:2602.17532 (2026-02) — attention captures co-expression not regulation
- Kendiukhov, arXiv:2603.02952 (2026-03) — SAEs, minimal regulatory logic
- 同作者 2603.01752 / 2602.22247 / 2603.11940；VCBench (bioRxiv 2026)
- **Delta**：它们用 CRISPR+TRRUST（RNA 邻近、带共表达混杂）；本项目用 ATAC/motif regulatory-potential proxy + 固定面板随机化推断 + 完整 readout 报告。
