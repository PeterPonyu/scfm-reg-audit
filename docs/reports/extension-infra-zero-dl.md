# Extension infra wave: local compute overlay (zero-download)

**Branch:** `w8/extension-infra-zero-dl`  
**Plan:** `prd-deepen-data-roles-and-carryout` Option **B′**  
**Ultragoal planId:** `deepen-data-roles-carryout`

## Built / stabilized

| Path | Role |
|------|------|
| `src/v2/extension/configs/{tissues,methods}.json` | Registry (BMMC=`construct_candidate`) |
| `src/v2/extension/registry.py` | load + dry-run + RNA-lake deny |
| `src/v2/extension/paths.py` | local ATAC / locked `G_ATAC` resolve |
| `src/v2/extension/cli.py` | first-class CLI |
| `src/v2/extension/emit_claim_pack.py` | JSON claim pack from public JSON |
| `docs/reports/extension-claim-pack/` | SI tables (docs overlay) |
| `src/v2/extension/construct_hooks.py` | Mantel/decomp on **existing** local NPZ |
| `src/v2/extension/baseline_stubs.py` | runnable Tier A–C emitters |
| `docs/reports/download_approval_optional_pilots.md` | download approval only |
| `src/v2/extension/scripts/fetch_optional_pilots.sh` | demoted pointer (no curl) |
| `src/v2/tests/test_extension.py` | hard merge gate |

Heavy runtime target: `results/v2/extension/` (local overlay; already MANIFEST-safe via `results/v2/`).

## Commands

```bash
python src/v2/extension/cli.py registry --json
python src/v2/extension/cli.py claim-pack
python src/v2/extension/construct_hooks.py --tissue fibroblast --execute --write
python src/v2/extension/baseline_stubs.py --all --execute --write
python -m unittest src.v2.tests.test_extension
python validate_artifacts.py
./src/v2/extension/scripts/fetch_optional_pilots.sh   # approval pointer only
```

## Explicit non-claims

- Construct `--execute` ≠ `build_atac_graph_v2` success  
- BMMC ≠ PeerJ `primary_audit` until panel policy + re-SAP  
- No 27/28 → Support / `G_ATAC`  
- No multi-GB download constructor in this package  

## Non-goals

- PeerJ package / MANIFEST primary lock edits  
- BMMC full FM audit  
- Forced Fig 1–12 redraw (SI/extension only until human go-ahead)  
