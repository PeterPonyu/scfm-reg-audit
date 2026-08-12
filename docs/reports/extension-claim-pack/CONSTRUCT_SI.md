# Construct SI — Mantel vs locked G_ATAC

Generated: 2026-08-12T09:30Z

Construct-lane Mantel Spearman and additive-fraction vs locked GSE174367 / PBMC10k / GSE206767. Extension SI only — does not add PeerJ Support FM rows or mutate cross-tissue decomp denominators.

- peerj_support_rows_touched: **False**
- Tissues packed: **3**
- OK pair rows: **9**
- Locked proxies: `['GSE174367', 'PBMC10k', 'GSE206767']`

## Per-tissue pairs

| tissue_id | TAG | locked_proxy | n_tf | observed Spearman | fraction additive | status |
|---|---|---|---:|---:|---:|---|
| descartes_spleen | DESCARTES_spleen | GSE174367 | 446 | 0.433944 | 0.8301 | ok |
| descartes_spleen | DESCARTES_spleen | PBMC10k | 446 | 0.391066 | 0.8462 | ok |
| descartes_spleen | DESCARTES_spleen | GSE206767 | 446 | 0.493373 | 0.8243 | ok |
| bmmc | GSE194122 | GSE174367 | 446 | 0.529094 | 0.6079 | ok |
| bmmc | GSE194122 | PBMC10k | 446 | 0.890248 | 0.4603 | ok |
| bmmc | GSE194122 | GSE206767 | 446 | 0.472997 | 0.689 | ok |
| orphan_treg_gse211155 | GSE211155_treg | GSE174367 | 446 | 0.464455 | 0.7419 | ok |
| orphan_treg_gse211155 | GSE211155_treg | PBMC10k | 446 | 0.506009 | 0.805 | ok |
| orphan_treg_gse211155 | GSE211155_treg | GSE206767 | 446 | 0.489878 | 0.7659 | ok |

## Sources

- `results/v2/extension/construct/DESCARTES_spleen/mantel_vs_locked.json`
- `results/v2/extension/construct/DESCARTES_spleen/additive_decomp_row.json`
- `results/v2/extension/construct/GSE194122/mantel_vs_locked.json`
- `results/v2/extension/construct/GSE194122/additive_decomp_row.json`
- `results/v2/extension/construct/GSE211155_treg/mantel_vs_locked.json`
- `results/v2/extension/construct/GSE211155_treg/additive_decomp_row.json`

