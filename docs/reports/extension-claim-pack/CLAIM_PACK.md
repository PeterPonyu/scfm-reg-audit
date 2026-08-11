# Local claim pack (JSON-only, zero download)

Generated: 2026-08-11T15:19Z

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

