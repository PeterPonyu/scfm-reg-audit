# HTAN GBM C3N-01334 — peak-set design spike (G008 D-design)

**Plan:** D3 open single-sample pilot  
**Pinned IDs:** `tissue_id=htan_gbm_pilot`, `g_atac_tag=HTAN_GBM_C3N01334`, `role=construct_candidate`  
**Canonical unpack root:** `${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/`  
**Scope:** design + tar inventory only (no full `prepare-htan` impl; no `tissues.json` write; no D4/D5)

---

## 1. Tar inventory (local, zero network)

### Asset

| Field | Value |
|-------|-------|
| Path | `${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/GSM7710026_C3N-01334_CPT0125220004_snATAC_GBM.tar.gz` |
| Size (ls) | **698M** (`731888013` bytes) |
| FETCH sidecar | `…tar.gz.FETCH.json` (plan=D3, GEO GSM7710026 / GSE240822 open FTP) |
| sha256 (FETCH) | `6e5c0248b03f935061efdc4297178947a20ea2c13acd967148405ce227769b25` |

### Members (`tar -tzf` — **exactly 3 entries**)

```
C3N-01334_CPT0125220004_snATAC_GBM/
C3N-01334_CPT0125220004_snATAC_GBM/outs/
C3N-01334_CPT0125220004_snATAC_GBM/outs/fragments.tsv.gz   # ~698 MB compressed member
```

### Keyword scan (`fragment|peak|barcode|meta|h5|mtx|bed|tsv|csv`)

Only hit: **`outs/fragments.tsv.gz`**.

### Present vs absent

| Artifact class | Present? | Notes |
|----------------|----------|-------|
| `fragments.tsv.gz` | **YES** | Sole data payload; Cell Ranger ATAC 2.0.0 |
| Peaks BED / `peaks.bed` / `peak_annotation` | **NO** | Not in tar |
| Peak×barcode matrix (`filtered_peak_bc_matrix`, `.h5`, `.mtx`) | **NO** | Not in tar |
| Barcodes TSV / singlecell.csv | **NO** | Only barcodes embedded as col-4 of fragments |
| Cell-type metadata | **NO** | Not in tar (GEO GSE240822 may publish annotations *outside* this tar; not local here) |
| RNA / multiome h5ad | **NO** | snATAC fragments only |

### Fragment format (stream peek, no full extract)

- **Pipeline header:** `cellranger-atac-2.0.0`; ref `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` → **hg38**
- **51** comment lines (`# …`), then standard 10x columns:
  - `chrom`, `start`, `end`, `barcode`, `count` (tab-separated, 5 cols)
- Example data row: `chr1  10067  10198  AGCGTGCTCAGGGTTT-1  1`
- First **1e6** data rows → **~69.7k** unique barcodes (sample only; full count needs full scan)

### Extraction policy for this spike

- **Did not** fully extract `fragments.tsv.gz` into `unpacked/` (≈700 MB compressed; design evidence sufficient from stream).
- Empty dir may exist at  
  `${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/unpacked/`  
  — optional later partial extract only if D-build needs a small header dump; **prefer** streaming / `tar -xOzf` for inventory.

**Conclusion:** this pilot tar is **fragments-only**. There is no builder-ready peak matrix and no local cell-type meta.

---

## 2. Peak-set strategy recommendation

| Strategy | Feasibility for this tar | Verdict |
|----------|--------------------------|---------|
| **`call_peaks`** | Possible in principle (MACS2/Signac from fragments) | **Out of wave scope** — heavy compute, new dependency surface, peak definition becomes a method choice that confounds construct Mantel |
| **`external_peaks`** | Need an external BED + count-from-fragments | Requires choosing a reference peak universe (e.g. paper CREs / another tissue’s peaks) **and** still building counts; not in tar |
| **`fixed_bins`** | Tiling genome → count fragments | Non-standard for `build_atac_graph_v2` peak contract; still full matrix build cost |
| **`blocked`** | Always available when local inputs cannot yield builder-ready peaks | **Recommended for G008 wave completion (E4.3b)** |

### Recommended path: **blocker-first (E4.3b acceptable)**

Lean **`status=blocked`** as the primary success mode for this wave:

1. Local asset is confirmed (tar present, fragments-only layout proven).
2. D3 becomes **`no_new_download`** (asset already on disk under `sample_pilot/`).
3. `prepare-htan` emits structured blocked JSON (not a silent skip).
4. Happy-path h5ad remains a **future dual-path** if peaks/meta appear later — do not invent peaks in D-build.

**Do not** ship a half-built peak caller in D-build unless a later wave explicitly expands scope.

---

## 3. Dual-path success contract

Mirror `prepare-bmmc` / `descartes-bridge`: default dry-run status JSON; optional `--write` under extension overlay.

### Happy path (future / if inputs grow)

**Preconditions (all required):**

- Peak×cell matrix available as h5ad (or convertible) with **builder-ready** `var_names`:
  - `chrom:start-end` (or hyphen form fixable via `normalize_peak_name` → `peak_name_is_builder_ready`)
- Optional but preferred: cell-type meta → `Barcode,Cell.Type` CSV.gz (like BMMC)

**Outputs:**

- `results/v2/extension/construct/HTAN_GBM_C3N01334/HTAN_GBM_C3N01334_atac_peaks.h5ad`
- Optional: `…/HTAN_GBM_C3N01334_cell_meta.csv.gz`
- Status: `prepared` / `prepared_present` / `ready_to_prepare`
- `build_command` via `builder_env_command(...)` with:
  - `SCREG_EXTENSION_OUT=results/v2/extension/construct/HTAN_GBM_C3N01334` (or tag root under `results/v2/extension/`)
  - `SCREG_PEERJ_SUPPORT_LOCK=1`
  - `TAG=HTAN_GBM_C3N01334`
  - `ATAC_FILE=<absolute prepared h5ad>`
  - `META_FILE=<meta or none>`

### Blocker path (expected **now** — wave success)

**Trigger:** tar inventory classifies as `fragments_only` (or peaks absent + no external peak matrix on disk).

**Status JSON shape (normative sketch):**

```json
{
  "schema_version": 1,
  "plan_id": "D3",
  "tissue_id": "htan_gbm_pilot",
  "g_atac_tag": "HTAN_GBM_C3N01334",
  "role": "construct_candidate",
  "network_fetch_performed": false,
  "status": "blocked",
  "block_reason": "fragments_only_no_peak_matrix",
  "peak_set_strategy": "blocked",
  "peak_set_strategies_considered": [
    "call_peaks",
    "external_peaks",
    "fixed_bins",
    "blocked"
  ],
  "accepted_wave_completion": "E4.3b",
  "pilot_dir": "${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/",
  "tar_path": "${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/GSM7710026_C3N-01334_CPT0125220004_snATAC_GBM.tar.gz",
  "tar_present": true,
  "inventory": {
    "n_members": 3,
    "has_fragments": true,
    "has_peaks": false,
    "has_peak_bc_matrix": false,
    "has_cell_type_meta": false,
    "fragments_member": "C3N-01334_CPT0125220004_snATAC_GBM/outs/fragments.tsv.gz",
    "pipeline": "cellranger-atac-2.0.0",
    "genome": "hg38"
  },
  "fm_support_allowed": false,
  "panel_policy_gate": "construct_candidate_only",
  "message": "Local D3 tar is Cell Ranger fragments-only; no peaks/matrix/cell-type meta in archive. Blocker path is intentional wave completion (E4.3b).",
  "next_steps": [
    "Do not call peaks in default prepare-htan.",
    "Optional future: supply external peak h5ad with chrom:start-end var_names under pilot_dir or extension overlay.",
    "D3 download gate: decision=no_new_download (asset already local).",
    "Never inflate PeerJ 13-row SAP / FM Support for this tissue."
  ],
  "build_command": null
}
```

Exit codes (align with BMMC):

- `0` — status reported successfully, including intentional `blocked`
- `2` — hard failure (missing tar when expected, inspect crash, path confinement violation)

---

## 4. CLI shape: `prepare-htan`

Mirror `prepare-bmmc` / `descartes-bridge`.

### Module (D-build)

- `src/v2/extension/htan_prepare.py` (name flexible; keep next to `bmmc_prepare.py`)
- Wire in `src/v2/extension/cli.py`:

```text
python src/v2/extension/cli.py prepare-htan [--tar PATH] [--pilot-dir PATH] [--execute] [--write]
```

| Flag | Behavior |
|------|----------|
| default | Dry-run inventory + status JSON to stdout |
| `--write` | Write `htan_prepare_status.json` under `results/v2/extension/construct/HTAN_GBM_C3N01334/` (confined) |
| `--execute` | Only perform real peak extract/rename **if** happy-path inputs exist; **must not** auto-run MACS2/peak calling |
| `--tar` / `--pilot-dir` | Override defaults under `…/htan/sample_pilot/` |

### Constants (pinned)

```python
TISSUE_ID = "htan_gbm_pilot"
G_ATAC_TAG = "HTAN_GBM_C3N01334"
DEFAULT_PILOT_DIR = DESKTOP_DATA / "datasets/extension_pilots/htan/sample_pilot"
DEFAULT_TAR_GLOB = "GSM7710026_C3N-01334_*_snATAC_GBM.tar.gz"  # or exact name
```

### Inventory function (core)

```python
def inventory_htan_tar(tar_path: Path) -> dict: ...
# tar -tzf style listing via tarfile; classify:
#   fragments_only | peak_matrix_present | unknown
```

Prefer `tarfile` + optional stream of first N fragment lines; **do not** require full unpack for blocked path.

### Reuse from `paths.py`

- `normalize_peak_name` / `peak_name_is_builder_ready` — happy path only
- `builder_env_command` — happy path only
- `assert_confined_write_path` / `redact_path` / `DESKTOP_DATA` / `HEAVY_ARTIFACT_ROOT`
- Add `LOCAL_ATAC_HINTS["htan_gbm_pilot"]` + optional `HTAN_PILOT_DIR` (D-build)

---

## 5. D3 download_gate post-state (required)

**Today (pre-build):** `decision=pending`, `fetchable=True`, recipe points at obsolete `htan_pilot/`.

**After D-build lands design+prepare+gate update**, D3 entry must become:

```python
"D3": {
    "title": "HTAN open single-sample pilot",
    "decision": "no_new_download",
    "fetchable": False,
    "rejected": False,
    "assets": [],
    "manual_recipe": [
        "No new download — tar already under "
        "${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/ "
        "(GSM7710026_C3N-01334_*_snATAC_GBM.tar.gz).",
        "Inventory / prepare: python src/v2/extension/cli.py prepare-htan [--write]",
        "Expect status=blocked (fragments-only) unless a builder-ready peak h5ad is supplied.",
        "If happy-path peaks appear: run emitted build_command "
        "(SCREG_EXTENSION_OUT + SCREG_PEERJ_SUPPORT_LOCK=1).",
        "Then construct-only: python src/v2/extension/cli.py construct "
        "--tissue htan_gbm_pilot --execute --write",
        "Never inflate PeerJ 13-row SAP / FM Support.",
    ],
}
```

**Do not** re-introduce fetch URLs into D3 assets. FETCH.json next to the tar is provenance only.

Also update any docs that still say `htan_pilot/` → **`htan/sample_pilot/`**.

---

## 6. Draft `tissues.json` entry (do **not** apply in D-design)

G008-build will first-insert after design lands. Draft only:

```json
"htan_gbm_pilot": {
  "accession": "GSM7710026",
  "series": "GSE240822",
  "sample_piece": "C3N-01334",
  "cancer": "GBM",
  "modality": "snATAC",
  "lane": "construct",
  "role": "construct_candidate",
  "genome": "hg38",
  "peerj_freeze": false,
  "panel_policy": "frozen_brain_pinned_446x1200",
  "fetch_plan_id": "D3",
  "fetch_status": "local_fragments_only",
  "g_atac_tag": "HTAN_GBM_C3N01334",
  "local_path_hint": "${DESKTOP_DATA}/datasets/extension_pilots/htan/sample_pilot/",
  "local_tar": "GSM7710026_C3N-01334_CPT0125220004_snATAC_GBM.tar.gz",
  "notes": "D3 open GEO snATAC pilot (HTAN/CPTAC GBM). Tar is Cell Ranger fragments-only (no peaks/matrix/cell-type meta). prepare-htan defaults to status=blocked (E4.3b). Construct-candidate only; no FM Support / PeerJ SAP inflation. Optional future: external peak matrix with chrom:start-end var_names."
}
```

Mirror the same block in `tissues.yaml` when inserting.

Registry: ensure `htan_gbm_pilot` is **not** in forbidden-role deny lists for construct dry-run (role is `construct_candidate`, not `lake_blocked`). `htan_whole_lake` remains D4 / `lake_blocked`.

---

## 7. Risks

| Risk | Severity | Mitigation |
|------|----------|------------|
| **Fragments-only tar** | High (blocks happy path) | Accept E4.3b blocked status as wave success |
| **No cell-type meta in tar** | High for type-aware SAP | Blocker; do not invent labels. GEO annotations exist for GBM/ccRCC per HTAN README but are **not** in this archive — future optional fetch is a **new** plan, not D3 re-download of the same tar |
| **Peak calling cost / method sensitivity** | High if forced | Do not default to `call_peaks`; keep strategy enum explicit |
| **Wrong pilot path (`htan_pilot/` vs `sample_pilot/`)** | Medium | Pin `…/htan/sample_pilot/` everywhere (gate recipe, LOCAL_ATAC_HINTS, design) |
| **Accidental network via `download --plan-id D3`** | Medium | Flip gate to `no_new_download` / `fetchable=False` / `assets=[]` |
| **PeerJ Support inflation** | Policy | `fm_support_allowed=false`; construct-only role |
| **Full unpack disk** | Low–med | Inventory via `tarfile`/stream; no mandatory unpack for blocked path |
| **Barcode suffix / multiplet conventions** | Low for blocked | Defer until happy-path matrix build |

---

## 8. Patterns to copy (existing code)

| Pattern | Source |
|---------|--------|
| Status JSON + `--execute` / `--write` CLI | `src/v2/extension/bmmc_prepare.py` |
| Local readiness / fail-closed absent | `src/v2/extension/descartes_bridge.py` |
| D2-style `no_new_download` recipe | `src/v2/extension/download_gate.py` D2 |
| Peak name normalize + builder env | `src/v2/extension/paths.py` (`normalize_peak_name`, `builder_env_command`) |
| CLI subparser wiring | `src/v2/extension/cli.py` (`prepare-bmmc`, `descartes-bridge`) |

---

## 9. D-build implementation checklist (next agent)

1. **[ ]** Implement `src/v2/extension/htan_prepare.py`
   - Constants: `TISSUE_ID`, `G_ATAC_TAG`, pilot/tar defaults under `…/htan/sample_pilot/`
   - `inventory_htan_tar()` via `tarfile` (no full extract)
   - `prepare_status()` → dual path: default **`blocked`** for fragments-only; happy path only if peak h5ad found + builder-ready names
   - **No** MACS2 / peak-calling in default `--execute`
2. **[ ]** Wire `prepare-htan` in `cli.py` (flags mirror `prepare-bmmc`)
3. **[ ]** Update `download_gate.py` **D3** → `decision=no_new_download`, `fetchable=False`, `assets=[]`, recipe → `prepare-htan` + `…/htan/sample_pilot/`
4. **[ ]** `paths.py`: `LOCAL_ATAC_HINTS["htan_gbm_pilot"]`, optional `HTAN_PILOT_DIR`
5. **[ ]** **First-insert** `htan_gbm_pilot` into `configs/tissues.json` + `tissues.yaml` (from §6 draft)
6. **[ ]** Optional: `construct_hooks.py` branch for `htan_gbm_pilot` (emit prepare-htan next step when G_ATAC absent)
7. **[ ]** Tests (local, no network):
   - Inventory fixture or mock tar listing → `fragments_only`
   - Status JSON keys / `status=="blocked"`
   - D3 gate decision fields
   - Path confinement for `--write`
8. **[ ]** Smoke: `python src/v2/extension/cli.py prepare-htan` → exit 0, `status=blocked`
9. **[ ]** Smoke: `python src/v2/extension/cli.py download --plan-id D3` → recipe only, no fetch, `no_new_download`
10. **[ ]** Do **not** implement D4/D5; do **not** call peaks; do **not** write peak h5ad from fragments in this wave unless inputs magically include peaks

---

## 10. Design decision (locked for G008)

| Question | Decision |
|----------|----------|
| Recommended path | **Blocker lean (E4.3b)** |
| Peak-set strategy now | **`blocked`** |
| Happy path | Kept as dual-path contract for future external peak h5ad |
| Full `prepare-htan` in D-design? | **No** — design only |
| `tissues.json` write in D-design? | **No** — draft only (§6) |
| D3 re-download? | **No** — post-state `no_new_download` |

---

*Evidence date: 2026-08-12. Inventory commands run against local tar only; zero network.*
