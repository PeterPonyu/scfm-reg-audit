# HTAN GBM C3N-01334 — D3 peak-set follow-up design (G013)

**Plan:** D3 open single-sample pilot  
**Pinned IDs:** `tissue_id=htan_gbm_pilot`, `g_atac_tag=HTAN_GBM_C3N01334`, `role=construct_candidate`  
**Scope:** design only — no D4/D5 fetch, no peak calling, no primary `results/v2/G_ATAC_*` edits  
**Inputs:** `design-spike.md`, `src/v2/extension/htan_prepare.py`, live `htan_prepare_status.json`  
**Evidence date:** 2026-08-12 (local only; zero network)

---

## 1. Current-wave status (completed)

| Checkpoint | Result |
|------------|--------|
| Local tar | Present under `${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/` (~698 MB) |
| Tar layout | **fragments-only** (`outs/fragments.tsv.gz`; no peaks BED / peak×bc matrix / cell-type meta) |
| `prepare-htan` | Implemented; dry-run / `--write` → structured **`status=blocked`** |
| `peak_set_strategy` | **`blocked`** |
| `block_reason` | `fragments_only_no_peak_matrix` |
| Wave completion | **`accepted_wave_completion=E4.3b`** |
| D3 download gate | `decision=no_new_download`, `fetchable=False`, `assets=[]` |
| `build_command` | `null` (no builder-ready peak h5ad on disk) |
| PeerJ / FM Support | `fm_support_allowed=false`; construct-candidate only |

**Conclusion:** The current wave is **done**. E4.3b blocked is the intentional success mode for fragments-only inputs. Happy path remains dual-path code only (pre-existing peak h5ad under pilot dir or extension overlay); it is not active for this tar.

Reference status artifact:

- `results/v2/extension/construct/HTAN_GBM_C3N01334/htan_prepare_status.json`

---

## 2. Peak-set strategies — offline options for a **FUTURE** wave only

These options are **not** in-scope for the completed E4.3b wave. They are design notes if a later implement gate re-opens HTAN GBM after new local inputs appear.

### 2.1 `external_peaks` (preferred future path if inputs appear)

**Idea:** Obtain or reuse a peak universe BED / peak×cell matrix; count HTAN fragments into that universe (or supply a finished h5ad with builder-ready `chrom:start-end` `var_names`).

| Source class | Feasibility | Notes |
|--------------|-------------|-------|
| **Locked construct tissues peak universe reuse** | Medium | Brain (`GSE174367`, ~219k peaks, hg38, `chrN:start-end`) and fibroblast (`GSE206767`, ~275k peaks) are locked construct-valid proxies with builder-ready names. **Brain is tissue-closest to GBM** for a shared CRE universe; PBMC (`PBMC10k` G_ATAC exists as locked NPZ) is immune, not neural — weaker biological match. Reuse = *peak coordinates only*, not cell labels or G_ATAC edges. |
| **Paper / GEO companion peak matrix for GSE240822** | Unknown locally | Design spike: cell-type annotations may exist outside this tar; **not present** under `sample_pilot/`. Any fetch is a **new plan**, not D3 re-download. |
| **Drop-in prebuilt h5ad** | Already coded | `htan_prepare._find_local_peak_h5ad` looks for `HTAN_GBM_C3N01334_atac_peaks.h5ad` / `atac_peaks.h5ad` under pilot or extension overlay. If present + builder-ready → `ready_to_prepare` / `prepared_present` (E4.3a). |

**Cost / risk:**

- Cost: fragment→peak counting is non-trivial CPU/RAM on ~700 MB compressed fragments; still far cheaper than unsupervised MACS2 design debates.
- Risk: **cross-tissue peak universe** (brain AD snATAC → GBM snATAC) changes peak definition and confounds Mantel vs tissues that use native peaks; must be labeled as `peak_set_strategy=external_peaks` with explicit source tag.
- Risk: still **no cell-type meta** for HTAN barcodes → type-aware SAP remains weak or N/A unless external barcode→type map is supplied.

**Hard precondition for implement gate:** external peak BED **and/or** peak×cell h5ad with `chrom:start-end` names, plus preferably barcode-level cell-type CSV.

### 2.2 `fixed_bins`

**Idea:** Tile hg38 (e.g. 5 kb bins) and count fragments per bin×barcode.

| Aspect | Assessment |
|--------|------------|
| Feasibility | Mechanically possible from fragments alone |
| Fit to builder | **Poor** — `build_atac_graph_v2` peak contract expects CRE-like peaks (`normalize_peak_name` / `peak_name_is_builder_ready`); bins are non-standard |
| Cost | Full matrix build (large `n_vars` if genome-wide) |
| Risk | Method choice confounds construct transfer vs all peak-based G_ATAC; high interpretability cost |

**Verdict for future:** research-only / ablation; **not** recommended as default HTAN construct path.

### 2.3 `call_peaks` (MACS2 / Signac-style)

**Idea:** Call peaks from fragments, then build peak×bc matrix.

| Aspect | Assessment |
|--------|------------|
| Feasibility | Possible in principle from Cell Ranger fragments |
| Cost | High CPU; new dependency surface (MACS2/Signac); full extract of fragments |
| Risk | Peak definition becomes a free method parameter that confounds Mantel vs locked tissues; QC thresholds multiply sensitivity |
| Policy | **Must not** become default `prepare-htan --execute` behavior |

**Verdict for future:** only under an explicit, scoped implement plan with frozen MACS parameters and a documented peak BED artifact — never silent auto-call.

### 2.4 `blocked` (current and default)

Already shipped. Correct when:

- tar is `fragments_only`, and
- no external/prebuilt peak h5ad on disk.

---

## 3. Explicit non-goals (this and immediate next gate)

| Non-goal | Reason |
|----------|--------|
| **No D4 fetch** | `htan_whole_lake` / DESCARTES lake remain `lake_blocked` / `pending_large`; TB-scale Synapse/RAW — not a D3 fix |
| **No D5 fetch** | Cancer/Development RNA lakes are `rna_lake_only`; forbidden for G_ATAC |
| **No default MACS in `prepare-htan`** | `execute` flag is reserved for happy-path extract only; code comments and CLI help forbid auto peak calling |
| **No PeerJ SAP / FM Support inflation** | Construct-candidate only; `SCREG_PEERJ_SUPPORT_LOCK=1` on any future overlay build |
| **No silent peak invention** | Do not write fake peak h5ad from fragments without a named future-wave strategy |

---

## 4. Go / no-go for next implement gate

### Recommendation: **NO-GO** (default)

| Gate question | Answer |
|---------------|--------|
| Is E4.3b wave complete? | **Yes** |
| Are builder-ready peaks local for HTAN C3N-01334? | **No** |
| Is cell-type meta local for HTAN barcodes? | **No** |
| Should next wave implement fragment→peak counting or MACS? | **No** without new inputs |
| Should D4/D5 open to “fix” HTAN? | **No** |

**NO-GO without:**

1. **External peak BED or prebuilt peak×cell h5ad** with builder-ready `chrom:start-end` `var_names` staged under pilot_dir / extension overlay, **and**
2. Preferably **barcode-level cell-type metadata** compatible with HTAN barcodes (or an explicit waiver that construct SI is peak-only / type-agnostic).

### Conditional GO (future only)

| Path | Minimum artifacts | Notes |
|------|-------------------|-------|
| A. Drop-in h5ad | `HTAN_GBM_C3N01334_atac_peaks.h5ad` + optional `…_cell_meta.csv.gz` | Re-run `prepare-htan` → expect E4.3a; then locked-overlay `build_command` + `construct --tissue htan_gbm_pilot` |
| B. External universe + count | Peak BED (e.g. brain GSE174367 var→BED export) + counted matrix | New implement plan; strategy tag `external_peaks`; no MACS default |
| C. Native MACS | Documented MACS2 params + peak BED + matrix | Highest risk; separate approval |

Until A/B/C inputs exist offline: **keep status=blocked; do not open an implement gate for peak calling.**

---

## 5. Relation to locked peak universes (design hint only)

| Locked tag | Local matrix (hint) | Role if reused for HTAN |
|------------|---------------------|-------------------------|
| `GSE174367` (brain) | `${DESKTOP_DATA}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad` | Best **coordinate** universe for neural/GBM (hg38, colon peak names). Reuse peaks, not cells/G_ATAC. |
| `PBMC10k` | Locked NPZ only in `results/v2/` (RNA partner monorepo) | Immune peak set — weak match for GBM. |
| `GSE206767` (fibroblast) | `${DESKTOP_DATA}/datasets/ATAC_data/GSE206767_filtered_peak_bc_matrix.h5ad` | Construct tissue; stromal, not GBM. |

**Do not** treat reusing brain peak coordinates as “running brain G_ATAC again.” Overlay builds must stay under `results/v2/extension/construct/HTAN_GBM_C3N01334/` with `SCREG_PEERJ_SUPPORT_LOCK=1`.

---

## 6. Decision table (locked for G013)

| Question | Decision |
|----------|----------|
| Current wave | **Complete via E4.3b blocked (fragments-only)** |
| Next implement gate | **NO-GO** until external peak BED/h5ad (+ preferably cell-type meta) |
| Default `prepare-htan` peak strategy | **`blocked`** |
| Default MACS2 | **Forbidden** |
| D4/D5 | **Out of scope** |
| Future preferred strategy if unblocked | **`external_peaks`** (drop-in h5ad or count into declared universe) |
| `fixed_bins` / `call_peaks` | Future research only; not default |

---

*G013 design-only deliverable. No network; no peak calling; no PeerJ mutation; no G_ATAC build.*
