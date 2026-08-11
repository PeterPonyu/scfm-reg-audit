# Extension infra wave: local compute overlay (zero-download)

**Branch:** `w8/extension-infra-zero-dl`  
**Plan:** `prd-deepen-data-roles-and-carryout` Option **B′**  
**Ultragoal planId:** `deepen-data-roles-carryout`

## Built / stabilized

| Path | Role |
|------|------|
| `src/v2/extension/configs/{tissues,methods}.json` | Registry (BMMC=`construct_candidate`) |
| `src/v2/extension/registry.py` | load + dry-run + RNA-lake deny |
| `src/v2/extension/paths.py` | local ATAC / locked `G_ATAC` / `SCREG_EXTENSION_OUT` resolve |
| `src/v2/build_atac_graph_v2.py` | honors `SCREG_EXTENSION_OUT` + PeerJ lock |
| `src/v2/extension/cli.py` | first-class CLI (+ descartes-bridge / prepare-bmmc) |
| `src/v2/extension/emit_claim_pack.py` | JSON claim pack from public JSON |
| `docs/reports/extension-claim-pack/` | SI tables (docs overlay) |
| `src/v2/extension/construct_hooks.py` | Mantel/decomp on **existing** local NPZ |
| `src/v2/extension/descartes_bridge.py` | D1 local RDS/h5ad → ATAC_FILE readiness |
| `src/v2/extension/bmmc_prepare.py` | BMMC multiome → peak ATAC under extension |
| `src/v2/extension/baseline_stubs.py` | runnable Tier A–C emitters |
| `docs/reports/download_approval_optional_pilots.md` | download approval + post-approve steps |
| `src/v2/extension/scripts/fetch_optional_pilots.sh` | demoted pointer (no curl) |
| `src/v2/tests/test_extension.py` | hard merge gate |

Heavy runtime target: `results/v2/extension/` (local overlay; already MANIFEST-safe via `results/v2/`).

## Commands (local-only; **no download**)

Do **not** run dataset fetches. `fetch_optional_pilots.sh` is an approval pointer
and exits without network I/O. Approval docs under `docs/reports/download_approval_*.md`
remain documentation only until explicit human sign-off.

```bash
# 1) registry + deny-lake policy
python src/v2/extension/cli.py registry --json

# 2) SI claim-pack from frozen public JSON (zero recompute)
python src/v2/extension/cli.py claim-pack

# 3) construct Mantel/decomp on locked local NPZ (fibro / brain / PBMC)
python src/v2/extension/cli.py construct --tissue fibroblast --execute --write

# 4) Tier A–C baseline emitters from local G_ATAC / ENCODE JSON
python src/v2/extension/cli.py baselines --all --execute --write --proxy-tag GSE174367

# 5) tests + PeerJ freeze gate
python -m pytest src/v2/tests/test_extension.py -q
python validate_artifacts.py

# approval pointer only (must print DISABLED; never curls)
./src/v2/extension/scripts/fetch_optional_pilots.sh
```

## 100% carry-out checklist (extension lane)

| Facility | Status |
|----------|--------|
| Registry + CLI | **complete** |
| Claim-pack emit | **complete** |
| Construct `--execute` on local NPZ / fibro | **complete** |
| Baseline emitters (Tier A–C) | **complete** (CollecTRI skips without local cache) |
| Extension tests + `validate_artifacts.py` | **complete** |
| Download constructor | **intentionally absent** (docs-only) |
| BMMC full FM / PeerJ 13-row rewrite | **out of scope / blocked** |

## Explicit non-claims

- Construct `--execute` ≠ `build_atac_graph_v2` success  
- BMMC ≠ PeerJ `primary_audit` until panel policy + re-SAP  
- No 27/28 → Support / `G_ATAC`  
- No multi-GB download constructor in this package  
- **No fetch execution** without a filled approval row + human go-ahead

## Non-goals

- PeerJ package / MANIFEST primary lock edits  
- BMMC full FM audit  
- Forced Fig 1–12 redraw (SI/extension only until human go-ahead)  
