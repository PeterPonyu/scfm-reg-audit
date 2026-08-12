# Local claim pack (JSON-only, zero download)

Generated: 2026-08-12T09:30Z

Extension SI numbers from existing `results/*.public.json`. Does **not** change PeerJ Support denominators.

## Dual-null rarity (C1)

- Full-spec FM rows: **13**
- Dual-null full (both BH q < 0.05): **7**
- Dual-null with ρ > 0: **6**
- Observed counts by tissue: `{'brain': 4, 'pbmc': 3}`

Independence OC P(≥1 dual-null row):
- brain: **0.0002** (n_sims=20000, family_n=8)
- pbmc: **0.0003** (n_sims=20000, family_n=5)

_If independence OC rate of ≥1 dual-null row is far below the event of observing 7 dual-null rows across tissues, dual-null Support is not explained by naive independent double-dipping alone — still not causal recovery._

## FM − baseline (C5)

- 4/7 dual-null full rows have ρ_full > tissue co-expression baseline; 4/13 FM rows beat baseline ignoring dual-null.
- Beats baseline: **4/13**
- Dual ∧ beats baseline: **4**

## Protocol-pass gates (C2)

- Frozen PeerJ protocol-pass: **0/13**
- Partial conjunction (excl. Multi-RO): **0** (diagnostic only; not a SAP pass)
- Gate counts: `{'dual_full': 7, 'dual_nondeg': 0, 'concordance': 5, 'nd_same_sign': 7, 'rho_gt_baseline': 4}`

Frozen PeerJ claim remains 0/13 protocol-pass (SAP §4 / Table 4). Multi-RO sign is not reconstructed here (null). Partial conjunction excluding Multi-RO is diagnostic only.

## Construct SI Mantel (extension overlay)

- See `CONSTRUCT_SI.md` / `construct_si_mantel.json`.
- Tissues: **3** (DESCARTES_spleen, GSE194122, GSE211155_treg)
- OK pairs vs locked ['GSE174367', 'PBMC10k', 'GSE206767']: **9**
- peerj_support_rows_touched: **False**

| TAG | locked_proxy | observed Spearman | fraction additive |
|---|---|---:|---:|
| DESCARTES_spleen | GSE174367 | 0.433944 | 0.8301 |
| DESCARTES_spleen | PBMC10k | 0.391066 | 0.8462 |
| DESCARTES_spleen | GSE206767 | 0.493373 | 0.8243 |
| GSE194122 | GSE174367 | 0.529094 | 0.6079 |
| GSE194122 | PBMC10k | 0.890248 | 0.4603 |
| GSE194122 | GSE206767 | 0.472997 | 0.689 |
| GSE211155_treg | GSE174367 | 0.464455 | 0.7419 |
| GSE211155_treg | PBMC10k | 0.506009 | 0.805 |
| GSE211155_treg | GSE206767 | 0.489878 | 0.7659 |

