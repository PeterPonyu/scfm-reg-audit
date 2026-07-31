# scReg-Eval — a regulation-aware evaluation protocol for scRNA foundation models

Status: v1, validated on 2 tissues (2026-07-13) · Deliverable #1 of DESIGN.md §8

## What it answers

*Does a frozen scRNA foundation model's gene-gene graph encode gene regulation beyond what
plain co-expression already gives you?* Not "does the FM cluster cells well" — specifically:
regulation, tested against a truth that cannot be satisfied by covariation alone.

## Why existing evaluations are insufficient

Most FM gene-network evaluations (including the Kendiukhov series this project's Sprint 0 pilot
initially replicated) score an FM graph against a truth built from the same or a highly correlated
data type — CRISPR perturbation transcriptomes, TRRUST, or co-accessibility — all of which share
co-expression's own statistical structure. A model that has merely learned "genes that go up
together" can look regulation-aware under those truths. scReg-Eval's truth is **sequence-derived**:
a TF motif match inside a chromatin-accessible peak is true or false independent of any RNA
covariation, in either the FM or the co-expression baseline.

## Protocol

### 1. Construct the truth graph `G_ATAC`

For a scATAC(-adjacent) dataset with a peak×cell count matrix on a known genome build:

1. **Peak → gene link**: assign each peak to a gene if its midpoint falls in
   `[TSS − 2kb, gene-end]` (strand-aware).
2. **Motif → peak**: scan each peak's genomic sequence (resized to a 500bp accessible core) for
   JASPAR CORE vertebrate PWM matches (MOODS, p < 1e-5, both strands). A TF "hits" a peak if any of
   its motifs match.
3. **Per-cell-type accessibility gate**: for each cell type T, weight peaks by mean log1p accessibility
   in T's cells.
4. **Edge score**: `G_ATAC_T[TF, target] = Σ_peaks-of-target  accessibility_T(peak) · motif-hit(TF, peak)`.

This is **directed** (TF rows only), **asymmetric**, and **orthogonal to co-expression by
construction** — nothing about it depends on cross-cell RNA covariation.

**Validity checks (run before trusting the truth graph):**
- Non-degeneracy: SVD top-1 energy fraction of the row-normalized TF×target matrix should be well
  below 1 (rank-1 collapse means "every TF hits every peak" — a motif-scan p-value/peak-width bug).
- TF distinctness: mean pairwise cosine between TF target-profiles should be well below 1.
- Spot-check a handful of canonical TFs against literature-known targets.
- Cross-tissue reproducibility (if ≥2 datasets available): the truth itself should replicate
  (Mantel test) — a low-reproducibility truth invalidates any downstream FM verdict.

### 2. Build the comparison graphs

On a **frozen, pre-registered gene manifest** (§3) — never select genes by ATAC variance or any
other statistic entangled with the outcome:

- `G_coexp` = `|Pearson|` co-expression from RNA (the confound baseline).
- `G_FM` — **three independent readouts**, because any one can be a degenerate proxy for the
  model's actual regulatory content:
  - **embedding**: cosine graph of per-gene contextual embeddings.
  - **attention**: last-layer gene-gene attention, mean over heads, symmetrized.
  - **in-silico perturbation**: token-deletion shift in every other gene's embedding, *with a
    random-gene-deletion control subtracted* (isolates the causal/perturbation signal from the
    token-removal/rank-shift artifact).

Preprocessing note (this is where a naive rebuild silently breaks): rank genes for Geneformer-style
models by **non-log** CP10k / gene-median — the tokenizer's rank encoding is defined on that scale,
not on log-normalized values. Verify this against the model's official tokenizer, not by inspection.

**Architecture applicability gate — check this before investing in a new FM's loader.** Not every
scFM architecture exposes a per-gene *contextual* embedding at all. Concretely: does the model process
a per-gene token sequence through its transformer/encoder (Geneformer, scGPT, scFoundation, UCE all
do — a per-gene hidden state is recoverable even when only cell-level output is officially exposed),
or does it pool gene expression into a single vector *before* any contextual processing (e.g. CellPLM's
`OmicsEmbedder`, which computes a single expression-weighted sum over a static gene-embedding table via
one `sparse.mm`, with no per-gene hidden state anywhere downstream)? For the latter class, the only
"gene embedding" available is a static, cell-independent lookup vector — testing that is not the same
claim as testing the other FMs' contextual representations, and would be misleading to report
side-by-side with them. Read the embedding/encoder forward pass before committing engineering time.

**Readout sanity check**: compute `Spearman(G_FM, G_coexp)`. If it's near zero or negative, the
readout is likely degenerate (e.g. embedding-cosine dominated by anisotropy) — do not base a claim
on a readout that doesn't even reproduce the well-established "attention ≈ co-expression" prior.

### 3. Frozen gene manifest

Universe = `hg38-coord genes ∩ FM-tokenizable ∩ ≥1 ATAC peak in window ∩ RNA-detected ≥1% of cells`,
sorted deterministically, sha256-hashed. Cap by a **pre-registered** rule (e.g. "keep all TFs, fill
remainder by RNA detection rate") — never by a statistic computed from the truth graph itself (the
original failure mode this protocol was built to avoid: an ATAC-variance-ranked gene cap silently
contaminates the eval set with the outcome).

### 4. Statistics

On the directed TF→target pair set `P = {(i,j) : i is a TF, i ≠ j}`:

- **Marginal**: `Spearman(G_coexp, G_ATAC)`, `Spearman(G_FM, G_ATAC)`.
- **Decisive (conditional)**: partial Spearman of `G_FM` vs `G_ATAC` controlling `G_coexp`
  (rank-residualize both on rank(coexp), correlate residuals). This is the number that answers the
  actual question — does the FM add anything *beyond* co-expression.
- **Null**: Mantel gene-label permutation (permute `G_ATAC`'s gene ordering, N≥1000) — graph-graph
  edges are non-independent, so a naive per-edge p-value is invalid.
- **Confound control** — mandatory before interpreting any positive partial: regress out target
  peak count, gene length, GC%, RNA detection rate, TF out-degree, target in-degree (all of these
  correlate with "how easy is this gene to find a spurious hit on," not regulation). Re-run the
  partial correlation on the residuals. **A positive partial that does not survive this step is a
  hubness artifact, not a regulatory finding** — this protocol's own pilot found exactly this in
  the in-silico-perturbation readout (z=2.9 raw, −0.002 after confound control) and it would have
  been a false positive without this step.
- Report per-cell-type as well as pooled, with a **confidence-gated** annotation (assign a cluster
  to a type only if `top_marker_zscore − second_marker_zscore ≥ threshold`; drop ungated clusters —
  do not force an argmax) and **exclude marker genes from the evaluation pair set** (they'd otherwise
  make the coexp/FM graphs trivially agree with the annotation-derived truth).

### 5. Decision rule

| pattern | reading |
|---|---|
| FM partial ≈ 0, survives confounds, stable across readouts + cell types | FM encodes no regulation beyond co-expression (confirmatory negative) |
| FM partial > 0, survives Mantel null AND confound control, in ≥2 independent readouts | genuine partial regulatory encoding — report which readout, effect size |
| FM partial > 0 but vanishes after confound control | readout artifact — report as a methods finding (what it correlates with instead), not a regulatory claim |
| truth graph itself doesn't replicate cross-tissue / fails non-degeneracy check | do not proceed to an FM verdict — fix the truth construct first |

## Reference implementation

`src/v2/` in this repo: `motif_utils.py` (JASPAR/MOODS), `build_atac_graph_v2.py` (truth
construction, parameterized by `ATAC_FILE`/`META_FILE`/`TAG`), `freeze_gene_manifest.py`,
`fm_readout.py` (corrected FM preprocessing, all 3 readouts), `crossmodal_v2.py` +
`confound_regression_v2.py` (test suite), `cross_tissue_atac_v2.py` (truth reproducibility),
`inspect_gatac_v2.py` (non-degeneracy checks), `ko_confound_check.py` (perturbation-readout
confound gate). `pbmc_eval_v2.py` is a from-scratch second-dataset application (10x PBMC
Multiome — paired RNA+ATAC on the same cells) demonstrating the protocol generalizes beyond the
original AD-brain pilot. `fm_readout_scf.py` + `crossmodal_scf_v2.py`/`pbmc_eval_scf_v2.py`/
`pertype_fm_scf_v2.py`/`scf_confound_check.py` extend the protocol to a 3rd, architecturally
distinct FM (scFoundation, encoder-only — its decoder needs an uninstallable dependency).

## Current status (2026-07-31)

The protocol is currently exercised on 13 pooled model/readout rows across brain and PBMC, with
six figures and an integrated TF-disjoint probe. Current coverage includes Geneformer embedding,
attention, knockout, and position-control readouts; scFoundation; UCE; scGPT; and a random-init
Geneformer floor. CellPLM remains out of pooled scope because its architecture does not expose a
contextual per-gene state. The probe's paired sign-flip q-values are UCE vs co-expression .0325,
random-init floor .0300, scGPT .255, Geneformer embedding .5125, and Geneformer attention .65.
These are mixed family-level results, not evidence that all FM readouts are null or that UCE has
model-specific regulatory recovery.

## Historical validation snapshot (2026-07-13)

The following three-FM summary records the earlier validation stage. It is retained for provenance;
it is not the current coverage statement.

**Historical 3-FM, 2-tissue validation summary:** Geneformer, scGPT, scFoundation × AD brain
(unpaired cross-study) + PBMC Multiome (paired single-dataset) × up to 3 readouts each (embedding,
attention, in-silico perturbation) — every positive partial correlation observed (scFoundation
pooled brain z=2.97 p=.002; perturbation-readout z=2.9 p=.004) collapsed to ≈0 under confound
control; every negative/null result stayed null. Co-expression's own regulatory recovery converges
to ~0.001–0.009 across both tissues once confounded. UCE was architecturally valid for this
protocol but had not yet been integrated at that stage; CellPLM was evaluated and found
architecturally incompatible (see above). Both UCE integration and the probe are now represented
in the current manuscript/results.


- Peak→gene linkage is a fixed genomic window, not a co-accessibility- or Hi-C-derived link
  (Cicero/LinkPeaks integration is a v2 upgrade, not yet implemented).
- Truth is unpaired-by-cell-type across independent studies unless the input is a multiome dataset;
  paired-cell calibration (§ methodology validation) is currently a single-dataset check, not a
  general guarantee.
- Effect sizes throughout are small (r² ≈ 0.001–0.01) — the protocol is well-powered to detect
  *whether* a signal exists and survives confounding, not to resolve fine-grained effect magnitudes.
