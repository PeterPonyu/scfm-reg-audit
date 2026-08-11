# Extension registries (local compute overlay)

Config-driven tissue and method registries for **extension-lane** work that must
not mutate the PeerJ freeze (`brain` + `pbmc` FM audit, fibroblast construct,
13 full-spec rows, 8 BH families). JSON is authoritative; YAML is a mirror.

| File | Role |
|------|------|
| `tissues.json` | Tissue ids, modality, lane (`construct` \| `audit`), panel policy, paths |
| `methods.json` | Graph/method ids, BH membership, estimand notes |
| `*.yaml` | Optional human-edit mirrors (loader uses JSON) |

## Rules

1. PeerJ validators continue to pin `tissues ∈ {brain, pbmc}`, `full_rows == 13`.
2. Heavy extension artifacts write under `results/v2/extension/` (local overlay; MANIFEST-safe).
3. SI claim-pack tables write under `docs/reports/extension-claim-pack/`.
4. Role tag `out_of_scope` (Cancer/Dev RNA lakes 27/28) must never feed Support / `G_ATAC`.
5. BMMC is `role=construct_candidate` (not `primary_audit`) until panel policy + re-SAP.
6. Construct `--execute` = Mantel/decomp on existing local `G_ATAC` NPZ (not graph rebuild).

```bash
python src/v2/extension/cli.py registry --json
python src/v2/extension/cli.py claim-pack
```
