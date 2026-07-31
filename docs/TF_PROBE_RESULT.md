# TF-disjoint probe — result, 2026-07-30

**Status:** Ran locally, complete. Remote GPU box is **not needed** for this arm.

**Question (from TF_DISJOINT_PROBE_SPEC.md).** Do FMs carry *recoverable* regulatory
signal under a supervised linear probe trained on disjoint TFs? This closes the
"you did not probe for it" rebuttal against the unsupervised negative.

**Answer.** No recoverable regulatory structure. Under confound adjustment, no FM family recovers
ATAC-proxy structure on held-out TFs above the null floor. Two family-level contrasts separate
from the co-expression baseline under paired sign-flip tests — UCE (q=0.0325) and the random-init
floor (q=0.0300) — but both sit at zero while the baseline sits below it, so neither is an FM
advantage; the floor result specifically rules out a model-specific reading of UCE. The spec's
"probe ties co-expression" branch applies in substance. Null seeds are explicit and deterministic
(`pair_probe_stats.py` schema 2).

## Design as executed

The spec assumed per-gene embeddings would be probed. Those are not on disk — the
pooled caches store gene x gene graphs only, one 1200x1200 matrix per family. Rather
than re-run inference, the probe was reformulated at pair level, which answers the
same question without new FM compute:

- **Sample** = ordered pair (TF, target gene). 311 train TFs x 1200 genes = 373,200
  train pairs; 135 test TFs x 1200 = 162,000 test pairs.
- **Features** (per family, from that family's graph only): `w` = S[tf, gene],
  `w_rev` = S[gene, tf], `tf_deg`, `g_deg`, `rank_g` = within-TF rank of `w`.
- **Target** = pooled PBMC ATAC-proxy G[tf, gene], averaged over 7 cell types.
- **Model** = RidgeCV (25 alphas, standardised features), fit on train-TF pairs only.
- **Readout** = per-test-TF Spearman between prediction and target, averaged.

Nothing in the feature construction is indexed by TF identity, so a probe fitted on
the 311 train TFs applies to the 135 held-out TFs unchanged. Train/test TF sets are
verified disjoint at load time, and every cache's gene order is checked against
`shared_genes.v2.json` before use.

### Two arms, because `all` lets the probe cheat

With all five features the probe puts almost no weight on the edge weights and
loads instead on `g_deg` (coef -0.012 vs `w` 0.001 for Geneformer). `g_deg` is a
gene-level quantity correlated with the gene-level confounds, so the probe reaches
rho = +0.11 marginally while learning nothing about edges. `edge_only`
(`w`, `w_rev`, `rank_g`) withholds that shortcut and is therefore the primary arm.

## Numbers

`results/v2/tf_probe_pair_eval_v2.json`, primary arm `edge_only`:

| family | marginal rho | adjusted rho | perm p | BH q |
|---|---|---|---|---|
| geneformer_attn | +0.0704 | -0.0102 | 0.001 | 0.002 |
| co_expression (baseline) | +0.0606 | -0.0120 | 0.001 | 0.002 |
| geneformer_embed | +0.0525 | -0.0163 | 0.001 | 0.002 |
| scGPT_encoder | +0.0281 | -0.0191 | 0.001 | 0.002 |
| UCE_encoder | +0.0100 | +0.0003 | 0.880 | 0.880 |
| random_floor | +0.0046 | +0.0023 | 0.285 | 0.342 |

Marginal rho is small for everything and *ordered with the baseline in the middle*:
Geneformer attention beats co-expression by +0.010 and scGPT loses to it by -0.033.
Paired sign-flip contrasts against co-expression: no FM family is significantly
better (geneformer_attn q=0.65, geneformer_embed q=0.51, scGPT q=0.26). The two
families that *do* separate significantly are UCE (q=0.033) and random_floor
(q=0.030) — both because they sit at zero while the baseline sits slightly below it,
which is not an FM advantage.

## The negative sign is over-correction, not reverse biology

Four families have adjusted rho significantly *below* the permutation null
(null width +-0.0022, so -0.012 is many null SDs out). This needed ruling out before
any claim, and `results/v2/tf_probe_pair_sensitivity_v2.json` does it — adjusted rho
across confound subsets:

| family | none | atac_construction | full | detv_only |
|---|---|---|---|---|
| co_expression | +0.0606 | -0.0026 | -0.0120 | +0.0635 |
| geneformer_embed | +0.0525 | -0.0120 | -0.0163 | +0.0281 |
| geneformer_attn | +0.0704 | -0.0054 | -0.0102 | +0.0907 |
| scGPT_encoder | +0.0281 | -0.0204 | -0.0191 | -0.0020 |
| UCE_encoder | +0.0100 | -0.0024 | +0.0003 | +0.0052 |
| random_floor | +0.0046 | +0.0021 | +0.0023 | +0.0006 |

`peakcount`, `genelen` and `gc` enter the ATAC-proxy's own construction, so
conditioning on them partly conditions on the target's definition. Conditioning on
those three alone (`atac_construction`) already collapses every family to within
+-0.02 of zero; adding `detv` moves nothing materially. The residual negative bias is
small relative to the unadjusted +0.06, and `random_floor` stays at zero throughout,
which is what a mild over-correction looks like — not a reversed biological signal.
Read the sign as "no signal, slightly over-corrected", not "anti-correlated".

## What this licenses

- **Claim supported:** a supervised linear probe trained on disjoint TFs does not extract
  confound-independent regulatory structure from FM gene-gene geometry above the null floor. Two
  family-level contrasts do separate from the co-expression baseline under paired sign-flip tests
  (UCE q=0.0325, random-init floor q=0.0300), but both sit at zero while the baseline sits below
  it — neither is an FM advantage, and the floor result rules out a model-specific reading of UCE.
  scGPT, Geneformer embedding, and Geneformer attention do not separate (q=0.255, 0.5125, 0.65).
- **Claim NOT supported:** "FMs contain no regulatory information." This probe is
  linear, operates on pairwise graph summaries rather than raw embeddings, and scores
  against an ATAC *proxy* that is itself a co-expression-family object (see
  `scfm-reg-audit-project` memory). A nonlinear probe on raw per-gene embeddings
  remains untested.
- **Proxy caveat carries over.** The 2026-07-13 peer review already established that
  `G_ATAC` is a co-expression-family proxy, not regulatory ground truth. That
  ceiling bounds this result exactly as it bounds the unsupervised one.

## Reproduce

```
python src/v2/tf_disjoint_split.py        # 311/135 split, seed 20260730
python src/v2/build_pair_features.py      # pair blocks + targets, ~30 s
python src/v2/run_pair_probe.py           # both arms, ~10 s
python src/v2/pair_probe_stats.py         # 999 perms x 6 families, ~6 min
python src/v2/pair_probe_sensitivity.py   # confound grid, ~5 s
```

Total under 10 minutes on CPU. The spec's Phase 4 remote-GPU staging is moot for
this arm; the remote box is still required only for the per-gene-embedding probe, if
that gets built.
