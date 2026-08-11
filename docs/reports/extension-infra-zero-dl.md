# Extension infra wave: zero-download scaffolding

**Branch:** `w8/extension-infra-zero-dl`  
**Plan:** `prd-infra-data-extension` Step 2 (+ local claim pack)  
**Ultragoal planId:** `local-infra-zero-dl`

## Built

| Path | Goal |
|------|------|
| `src/v2/extension/configs/{tissues,methods}.json` | G001 registry |
| `src/v2/extension/registry.py` | load + dry-run + RNA-lake deny |
| `src/v2/extension/emit_claim_pack.py` | G002 JSON claim pack |
| `docs/reports/extension-claim-pack/` | G002 emitted SI tables |
| `src/v2/extension/construct_hooks.py` | G003 fibro-style stubs |
| `src/v2/extension/baseline_stubs.py` | G003 Tier A–C stubs |
| `docs/reports/extension-plans/` | G003 dry-run plans |
| `src/v2/extension/scripts/fetch_optional_pilots.sh` | estimate-only fetch helper |
| `docs/reports/optional_cancer_dev_download_costs.md` | G004 cost model |
| `docs/reports/bmmc-panel-policy-memo.md` | panel policy (no FM) |

Heavy runtime NPZ target (not committed): `results/v2/extension/` (already a local-worktree overlay).

## Commands

```bash
python src/v2/extension/registry.py --json
python src/v2/extension/emit_claim_pack.py
python src/v2/extension/construct_hooks.py --tissue bmmc --write
python src/v2/extension/baseline_stubs.py --all --write
python validate_artifacts.py
./src/v2/extension/scripts/fetch_optional_pilots.sh   # estimates only
```

## Non-goals (this wave)

- No PeerJ package / MANIFEST primary lock edits  
- No multi-GB downloads  
- No 27/28 → Support  
- No BMMC full FM audit  
