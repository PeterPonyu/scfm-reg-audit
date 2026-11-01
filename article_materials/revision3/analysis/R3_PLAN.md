# ARRAY-D-26-05613 — Revision 3 analysis plan (fixed before any new result)

Written 2026-10-31, before any of the analyses below was run. The SHA-256 of this file is
recorded in `analysis/R3_PLAN.sha256` at creation. Every new analysis is a supplement to the
registered corrective analysis and the post hoc 9,999-draw extension; neither is replaced.

Common conventions for all new statistics
- Statistic: partial Spearman correlation between a model graph and the corrective
  motif-accessibility proxy over directed non-self TF-target pairs, exactly as in the
  registered corrective code (`screg_validation.batched_null.batched_components`).
- Components: gene-label relabelling (M) and within-TF-row target shuffle (D), restricted
  (TF status x peak-count quintile x GC tertile) or unrestricted; own degree recomputed per draw.
- Draws: 9,999 per component (plus-one rule; floor 1e-4). Seed 20261001, a fresh generator
  per stream, as registered. Probabilities two-sided; pIUT = max(pM, pD).
- Selection: adjusted value strictly below 0.05, BY and Holm, within the named family.
- Every family below is fully measured. No slot enters a family at an administrative p = 1.

R3-1 Calibration of the restricted scheme (Reviewer 1, comment 1)
- Same synthetic construction and seeds as the published stress test (200 replicates per
  scenario, 999 transformations, raw p <= 0.05), contexts: historical mixed-source, historical
  PBMC, and the genuine-brain brain-selected panel.
- Scenarios: independent noise; shared heteroskedastic TF/target shocks scaled by the
  stratifying variables (the published scenario); the same shocks scaled instead by RNA
  detection, which is not a stratifier (off-stratifier scenario).
- Nulls: unrestricted and restricted on identical synthetic readouts. The unrestricted
  counts must reproduce the published counts exactly before restricted counts are read.
- Rule: a component is called nominal when its Wilson 95% interval contains 0.05 and its
  rate lies in [0.025, 0.075] (Bradley's liberal criterion). If the restricted joint test is
  not nominal under the published shock scenario in a historical context, the 5/13 count is
  reported only as a sensitivity result and supports no selection claim.

R3-2 Common-support analysis (comment 2)
- Per context, S = panel genes with positive RNA variance and finite gene-level scores for
  every readout; E = all non-self TF-target pairs within S with finite scores for all six
  readouts (Geneformer attention is pair-level). Degree, strata and quantile cut points are
  recomputed within S.
- Family: 3 contexts x 6 readouts x 2 degree modes x 2 transformation classes = 72 slots.
- Sensitivity: the five gene-level readouts on all pairs in S (attention excluded), 60 slots.

R3-3 Symmetric baseline (comment 5)
- Model readout and co-expression are each related to the proxy with the identical
  structural covariates (peak count, gene length, RNA detection, linked-window GC; plus own
  degree in the own-degree mode), with co-expression in neither nuisance set.
- Same pairs E, same M/D draws (9,999). Difference = rho_FM - rho_coexpression with a
  TF-cluster bootstrap 95% percentile interval (2,000 resamples of TF rows, seed 20261031).
  The interval is descriptive and conditional on the panel; it is not a donor interval.

R3-4 Ten-design degree grid on source-matched contexts (comment 4)
- The ten historical designs (none, linear own, out only, in only, log, rank, natural
  spline with five fixed knots, weighted degree, fibroblast-mix proxy degree, other-tissue
  proxy degree) on the three genuine contexts, pairs E, all six readouts, unrestricted
  components, 9,999 draws. External proxies use the corrective construction on the same panel.
- Families: within context x design (6 tests) and pooled across all 180 tests.

R3-5 Disjoint discovery and evaluation samples (comment 3)
- The 18 GSE174367 sample identifiers are assigned to two halves of nine, stratified by
  diagnosis and batch, seed 20261031; each half is discovery once and evaluation once
  (two folds). Cells are the already-registered 103 nuclei per sample.
- Discovery half: panel selection with the registered rule (all eligible motif TFs, fill to
  1,200 genes by detection with symbol tie-break, detection >= 1%). Evaluation half: model
  inference for the six readouts, co-expression, and the proxy from evaluation-half ATAC
  nuclei. Common support (R3-2 rule) within the evaluation half.
- Family: 2 folds x 6 readouts x 2 degree modes x 2 transformation classes = 48 slots.
- Pretraining exposure is not removed by this split and is reported by model.

R3-6 Five-donor bone-marrow analysis (comment 3)
- Gene set: genes with finite scores in all five evaluation donors for the five gene-level
  readouts; attention remains unavailable (pair coverage 73-93% per donor).
- Per donor x readout x degree (reference, own, none) x transformation class, 9,999 draws;
  family 5 x 5 x 3 x 2 = 150 slots. Donor summaries as registered: equal-donor bootstrap
  (2,000 resamples) of the reference-degree restricted effect and of the FM-minus-baseline
  difference, and the exact paired sign-flip test. With five donors the smallest two-sided
  sign-flip p is 0.0625, so donor-level significance at 0.05 is unattainable by design;
  the donor intervals are the inferential summary.

R3-7 scFoundation weights (comment 6)
- Local file SHA-256 against the Hugging Face LFS record of the mirror it came from, then
  tensor-by-tensor comparison with the official BioMap `models.ckpt` (`cell` key) if the
  official file is obtained. No equivalence is claimed beyond what is compared.

Interpretation fixed in advance
- Source-matched non-rejections are not equivalence, and anti-conservative components
  would make them more, not less, conservative; only positive selections depend on
  component validity.
- A selection in any new family is reported with its sign; negative alignment is not
  positive regulatory recovery.

## Addendum A (2026-10-31, written after R3-1 calibration and before any R3-2..R3-6 result existed)
Calibration result (R3-1): under the requested shared-shock control the restricted joint
test is nominal in all three contexts (8/200 mixed-source, 5/200 PBMC, 5/200 genuine brain)
and remains nominal with detection-scaled shocks (6/200, 8/200, 9/200); the unrestricted
joint test is not (18/200 mixed-source, 35/200 and 33/200 genuine brain; PBMC 13/200 is
borderline). The restricted target-shuffle component alone is anti-conservative (16-36/200).
Decision: for every new source-matched family (R3-2, R3-5, R3-6), the restricted
transformation class is the primary inferential basis; unrestricted results are reported
beside it as a sensitivity analysis. Families and slot counts in the plan are unchanged.
The historical 5/13 restricted count passes the pre-specified calibration rule; it remains
historical mixed-source/PBMC evidence with the brain interpretation withdrawn.

## Addendum B (2026-10-31, written after every R3-2..R3-6 family result; post hoc)
After the family results were known, an internal check noted that R3-1 calibrated only one of the
three common-support contexts and none of the bone-marrow donors. The R3-1 construction, seeds,
scenarios and rule were therefore applied, outside this plan, to the remaining common-support
contexts (new panel in PBMC, old panel in brain) and to the two bone-marrow donors with restricted
selections (13272 and 28483); the runs finished at 08:10, after the last family result (06:42).
These four calibrations are reported as post hoc. Their records carry the hash of the original
plan because the same code produced them. No family, slot count, transformation class or
selection rule was changed. The hash recorded for Addendum A covers the text up to the end of
Addendum A.
