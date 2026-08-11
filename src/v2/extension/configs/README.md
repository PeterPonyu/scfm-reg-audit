# Extension registries (zero-download scaffolding)

Config-driven tissue and method registries for **extension-lane** work that must
not mutate the PeerJ freeze (`brain` + `pbmc` FM audit, fibroblast construct,
13 full-spec rows, 8 BH families).

| File | Role |
|------|------|
| `tissues.json` | Tissue ids, modality, lane (`construct` \| `audit`), panel policy, paths |
| `methods.json` | Graph/method ids, BH membership, estimand notes |
| `*.yaml` | Optional human-edit mirrors (loader uses JSON) |

## Rules

1. PeerJ validators continue to pin `tissues ∈ {brain, pbmc}`, `full_rows == 13`.
2. Heavy extension artifacts write under `results/v2/extension/` (gitignored local overlay).
3. SI claim-pack tables write under `docs/reports/extension-claim-pack/`.
4. Role tag `out_of_scope` (Cancer/Dev RNA lakes 27/28) must never feed Support / `G_ATAC`.
5. BMMC may be registered as `lane=construct` with explicit `panel_policy`; full FM audit into PeerJ Support is gated (PRD Option C).

Load via `python src/v2/extension/registry.py`.
