# Orphan ATAC meta inventory (G015)

**Tissue registry key:** `orphan_atac_placeholder`  
**Local root:** `${DESKTOP_DATA}/datasets/ATAC_data/`  
**Scope:** inventory only — **do not** build `G_ATAC`; no network; no peak calling; no PeerJ mutation  
**Evidence date:** 2026-08-12  

Registry notes (`src/v2/extension/configs/tissues.json`):

> ~111 peak matrices local; obs often empty. Construct hooks only after cell-type metadata recovery. No Support rows.

This inventory updates counts and ranks **candidates for a future construct SI** after meta rescue — not an implement plan to build graphs now.

---

## 1. File counts

Root size: **~36–37 GB** under `${DESKTOP_DATA}/datasets/ATAC_data/`.

| Class | Count | Notes |
|-------|------:|-------|
| **`.h5ad` peak / feature matrices** | **115** | All top-level (no nested h5ad duplicates) |
| **`.h5` (Cell Ranger-style)** | **38** | Mostly under `*_RAW/` dirs; several series already mirrored as h5ad |
| **`.mtx.gz` peak matrices** | **2** | `GSE190162_RAW` only (MCF7 ± CPI1612) |
| **`peaks.bed(.gz)`** | **2** | Same GSE190162 pair |
| **Barcode TSV sidecars** | **2** | GSE190162 barcodes only |
| **Annotation / cell-type CSV·TSV** | **0** | No dedicated meta tables in tree |
| **Peak-matrix candidates (h5ad+h5+mtx)** | **155** | Union for “matrix-like” inventory |

Registry “~111” is slightly stale vs current **115** h5ad files (same lake; drift expected).

### 1.1 Related local paths (outside ATAC_data, for completeness)

| Path | Role |
|------|------|
| `${DESKTOP_DATA}/datasets/ATAC_data/` | `orphan_atac_placeholder` + brain / fibroblast locked matrices |
| `${DESKTOP_DATA}/datasets/extension_pilots/descartes_spleen/descartes_spleen_peaks.h5ad` | Already construct-wired (not orphan) |
| `${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/*.tar.gz` | HTAN fragments-only pilot (G013; not peak matrix) |
| `${DESKTOP_DATA}/external/scfm-reg-audit/gse194122/…` | BMMC multiome (construct candidate; not orphan lake) |

**Already claimed (not “orphan” for new construct SI):**

- `GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad` → primary brain audit / locked `G_ATAC_v2_GSE174367`
- `GSE206767_filtered_peak_bc_matrix.h5ad` → construct fibroblast / locked `G_ATAC_v2_GSE206767`

---

## 2. Cell-type-like `obs` columns

**Method:** `anndata.read_h5ad(..., backed="r")`; inspect `obs.columns` only (no full `X` load).  
**Coverage:** **all 115** top-level `.h5ad` files.

| Metric | Value |
|--------|------:|
| Files scanned | 115 |
| Files with **any** non-index `obs` columns | **0** |
| Files with cell-type-like columns (`cell_type`, `cluster`, `annotation`, …) | **0** |
| Sampled “largest 40” with CT-like columns | **0 / 40** |
| Remaining 75 with CT-like columns | **0 / 75** |

Typical structure (example: brain and many GSM peaks):

- `obs` keys: only `_index` (barcodes)
- `uns`: often `file_info` only
- `var_names`: mostly builder-friendly `chrN:start-end` (peak matrices)
- Exception class: `*_filtered_feature_bc_matrix.h5ad` (e.g. GSM8685152/5155/5177) → **gene-like** `var_names` (`Xkr4`, …), not ATAC peaks — **unsuitable** as peak G_ATAC input without a different pipeline

**Sidecar meta search:** recursive `*meta*`, `*annot*`, `*cell*type*`, csv/tsv under `ATAC_data` → only the two GSE190162 barcode lists. **No recoverable cell-type map is on disk.**

### 2.1 Implication

The placeholder note is confirmed and strengthened: **obs is empty lake-wide**, not “often empty.”  
**Construct hooks / G_ATAC builds for orphan tissues must wait on cell-type metadata recovery** (or an explicit type-agnostic construct waiver in a future plan). **Do not build G_ATAC from this inventory alone.**

---

## 3. Series clusters (for ranking)

| Cluster (filename / RAW) | ~n h5ad | Rough biology (from names only) |
|--------------------------|--------:|----------------------------------|
| GSM6449* (GSE211155_RAW) | 8 | Human Treg / PBMC / NSCLC TIL scATAC |
| GSM8248* (GSE266511_RAW) | 9 | BCG / immune LN time course |
| GSM7734* + related MGH | 7+ | MGH samples + Naive/D2/D8 timepoints |
| GSM7884* | 6 | Lobe / macular ATAC (retina-ish) |
| GSM8462* | 5 | AL\* ATAC multi-sample |
| GSM8900* | 6 | sample1–6 ATAC multi-sample |
| GSM7777* | 4 | CTRL / 2DG / DON / AOA metabolic |
| GSM7062* / GSM7064* | 5+ | D-series time course |
| GSM6052* | 3 | E13 / E16 / P0 — **likely mouse** development |
| GSM8685* | 3 | `feature_bc` — **gene features, not peaks** |
| GSE198730* | 2 | ESC / aPSM (likely mouse) |
| Singleton large | 1 | GSE174367 (claimed), GSE206767 (claimed) |

RAW-only residual (not fully converted to top-level h5ad): GSE192947 (11 h5), GSE199556 (2), GSE225803 (1), GSE275786 expression txt, GSE284492, GSE292195, GSE190162 mtx+bed.

---

## 4. Ranked top candidates for **future** construct SI

Ranking criteria (no G_ATAC built):

1. Not already locked as brain/fibroblast primary construct  
2. Peak-style `var_names` (`chr:start-end`), not gene feature matrix  
3. Prefer **human** (low first-peak genomic offsets / naming) over clear mouse-dev sets  
4. Multi-sample series (replicates / conditions) for SI narrative  
5. Larger `n_obs` / `n_vars` when meta eventually lands  
6. Still blocked on **meta** — rank is “best when annotations appear,” not “build now”

| Rank | Candidate set | Why | Blocker |
|-----:|---------------|-----|---------|
| **1** | **GSM6449\* Treg/PBMC/NSCLC TIL (×8)** | Human immune multi-sample; clear cell-program story; peak names look hg38-like; pairs with known PBMC audit narrative **without** reusing locked PBMC10k G_ATAC | Empty obs; need GEO/paper cell-type map |
| **2** | **GSM8248\* BCG / LN time course (×9)** | Largest multi-sample immune series; time/condition structure for SI | Empty obs |
| **3** | **GSM8900\* sample1–6 (×6)** | Multi-sample ATAC; large n_obs (up to ~20k); peak matrix | Biology less obvious from names; empty obs |
| **4** | **GSM8462\* AL\* ATAC (×5)** | Multi-donor/sample; large peak sets (~200k vars) | Empty obs; organism confirm needed |
| **5** | **GSM7884\* lobe/macular (×6)** | Tissue-distinct (retina/eye) construct diversity vs brain/fibro/immune | Empty obs |
| **6** | **GSM7777\* metabolic CTRL/2DG/DON/AOA (×4)** | Clean condition panel; large peak universes | Empty obs; may be non-human — verify before SI |
| **7** | **GSM7734\* MGH + Naive/D2/D8** | Clinical-ish multi-sample + timepoints | Empty obs; mixed stems |
| **8** | **GSE198730 ESC/aPSM (×2)** | Development angle | Likely **mouse**; empty obs; genome policy |
| **—** | **GSM8685\* feature_bc (×3)** | Large n_obs but **gene features** | Wrong modality for peak G_ATAC |
| **—** | **GSM6052\* E13/E16/P0** | Large matrices | **Mouse** development; not first human construct SI |
| **—** | GSE174367 / GSE206767 | Best matrices | **Already claimed** — not orphan |

### 4.1 Illustrative single-file stats (backed read; empty obs throughout)

| File | n_obs | n_vars | ~MB | Notes |
|------|------:|-------:|----:|-------|
| GSE174367_…h5ad | 143401 | 219070 | 1739 | Claimed brain |
| GSM8685152_…feature_bc… | 27166 | 31053 | 516 | Gene features — skip for peaks |
| GSM8900551_sample4_… | 20020 | 64152 | 339 | Orphan multi-sample rank 3 |
| GSM8397466_WT1_… | 18157 | 148834 | 550 | Large singleton WT |
| GSE206767_… | 15259 | 275448 | 700 | Claimed fibroblast |
| GSM7734290_Naive_… | 12883 | 103294 | 540 | Time-course family |
| GSM6449881_HD4_CD4CD25… | 3473 | 101785 | 382 | Treg series rank 1 |
| GSM8248670_8WKs_PosBCG… | 6718 | 127081 | 503 | BCG series rank 2 |

---

## 5. Go / no-go for orphan construct implement

| Question | Answer |
|----------|--------|
| Can we pick a matrix and build G_ATAC today? | **NO** — zero cell-type meta lake-wide |
| Is the orphan lake empty of peaks? | **NO** — 115 peak-like h5ads are present |
| Next useful work | **Meta rescue** (offline GEO companion tables, or paper Supplements already on disk elsewhere) then re-scan `obs` / join on barcodes |
| Build G_ATAC in this task? | **Explicitly no** |

**Recommendation one-liner:** Hold `orphan_atac_placeholder`; prioritize meta recovery for GSM6449 Treg/PBMC series (then GSM8248 BCG) before any construct SI or G_ATAC build.

---

## 6. Method notes / constraints honored

- Zero network  
- `backed='r'` only; no full matrix loads  
- No peak calling; no G_ATAC NPZ writes  
- No edits under `results/v2/G_ATAC_*` or PeerJ freezes  
- Counts from local filesystem 2026-08-12  

---

*G015 inventory-only deliverable.*
