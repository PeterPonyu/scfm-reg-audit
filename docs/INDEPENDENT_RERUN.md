# Independent rerun instructions

This repository contains the fixed-panel motif-accessibility audit and its public numerical
summaries. Inputs are project-local or explicitly configured; no other research repository
is discovered implicitly. The archived release and the subsequent manuscript supplements
have separate contents. The instructions below distinguish calculations supported by the
released summaries from reconstruction steps that still require verified upstream inputs.

## Validated public-summary layer

Use the declared Python environment (`requirements.txt`; reference Python 3.13.7):

```bash
python -m pip install -r requirements.txt
python src/revision_05613/multiplicity.py --summary-only
python src/revision_05613/mde.py --summary-only
python -m unittest discover -s src/v2/tests -p test_independent_release.py -v
```

The validation performed for these commands used an existing local environment, not a newly
installed environment. No dependency was added. Installation is a setup recipe, not a claim
that a fresh environment installation was validated.

The summary commands write respectively:

- `output/revision_checks/multiplicity_robustness.json`
- `output/revision_checks/mde.json`

Use `--output /path/to/check.json` to choose a separate output location. The frozen JSON in
`results/revision_05613/` is not overwritten by default. Input selection is explicit configured
file first, then `results/v2/<name>.json`, then `results/<name>.public.json`. An invalid
explicit configuration fails instead of silently substituting another file. For a guaranteed
public-only check, even in a development tree containing private/cache outputs, use:

```bash
SCREG_AUDIT_JSON="$PWD/results/fixed_panel_audit_v2.public.json" \
  python src/revision_05613/multiplicity.py --summary-only
SCREG_AUDIT_JSON="$PWD/results/fixed_panel_audit_v2.public.json" \
SCREG_PROBE_JSON="$PWD/results/tf_probe_pair_stats_v2.public.json" \
  python src/revision_05613/mde.py --summary-only
```

The outputs hash the actual JSON bytes selected. Public redaction changes provenance paths
and consequently input hashes; archived pre-redaction input hashes need not equal public-file
hashes. The regression fixtures compare the numerical fields directly.

Multiplicity recomputes BH, BY, Holm and the already archived per-row IUT corrections from
published N=999 Monte Carlo p-values. Analytic MDE recomputes the stated normal-approximation
formula using published null standard deviations. These are arithmetic checks conditional on
the released summaries. They do not regenerate randomization samples, observed graph
statistics, empirical MDE, Westfall–Young, or new manuscript supplemental analyses, and are
not independent validation of the biological interpretation or population power.

## PBMC 10x H5 conversion

The PBMC input is the 10x 10k PBMC Multiome Chromium X reference processed with Cell Ranger
ARC 2.0.0 against GRCh38. The [original filtered matrix](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5)
and [official summary](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_summary.csv)
identify 10,970 cells and 111,743 peaks. The retained matrix is 166,323,468 bytes, with SHA-256
`3897be5c916a66def9273049f6c5418ae042111f1b64ce769198253f56bf356d`.
On 28 September 2026, its byte count and locally recomputed multipart checksum matched the
public object. This verifies the specific input object, not the downstream inference pipeline.

For an already obtained 10x Cell Ranger ARC filtered feature-barcode H5 containing both
`Gene Expression` and `Peaks` features:

```bash
python src/v2/split_pbmc_multiome.py \
  --input /path/to/pbmc10k_multiome_filtered.h5 \
  --out-dir data/multiome
```

This writes `pbmc10k_rna.h5ad` and `pbmc10k_atac.h5ad`, retaining the same barcode order.
RNA uses gene symbols with duplicate names made unique; ATAC uses the peak feature IDs.
Missing source files or a missing modality fail before output creation. A small 10x-format
fixture verifies counts, barcodes and feature names. Full public PBMC source acquisition and
model inference were not executed in the summary verification. Repeated successful conversion into the
same output directory replaces those two files; use a new directory to retain prior inputs.

## Configured inputs and fail-closed paths

- `SCREG_DATA_ROOT`: explicit brain dataset root; otherwise the project's `data/`.
- `SCFM_BRAIN_ATAC`: explicit brain ATAC H5AD for the audit/shared-null/revision drivers;
  otherwise `data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad`
  under the configured data root.
- `SCREG_PBMC_ATAC`: explicit PBMC ATAC H5AD for the audit/shared-null/revision drivers;
  otherwise `data/multiome/pbmc10k_atac.h5ad` under the project.
- `SCREG_GENE_COORDS`: explicit gene coordinates for the shared-null script; otherwise
  project `data/annotation/gene_coords_hg38.tsv`. This override is not a universal setting
  for all historical scripts; those continue to expect their documented project-local inputs.
- `META_FILE`: cell-type metadata CSV.gz with `Barcode` and `Cell.Type` columns for
  `build_atac_graph_v2.py`; otherwise project `data/annotation/atac_cell_meta.csv.gz`.
  Missing metadata fails. `META_FILE=none` is the explicit choice for a pooled proxy and is
  a different analysis; it is not a substitute for reproducing the archived cell-type proxy.
  Metadata with no matching labelled ATAC barcodes also fails.
- The proxy builder's ATAC selector is `ATAC_FILE`, with `TAG` controlling its output names.
  Set both explicitly for each tissue. Its default data root is `SCREG_DATA_ROOT` or project
  `data/`. The builder continues to need local coordinates, motif MEME and hg38 FASTA inputs.

For example, this selects supplied inputs for the brain proxy; it is not a validated full-run
recipe or a guarantee that the supplied files match the original frozen inputs:

```bash
ATAC_FILE=/path/to/brain_atac.h5ad \
META_FILE=/path/to/brain_cell_metadata.csv.gz \
TAG=GSE174367 python src/v2/build_atac_graph_v2.py
```

Missing PBMC or gene-coordinate inputs are reported explicitly; a neighbouring research
workspace is never used to supply them.

## Remaining requirements for a full rerun

The exact brain RNA accession, selected cells, donor mapping and deterministic preprocessing
are still unverified. A surviving processed file and matching detection values establish
file continuity, but do not supply those missing source records. Cell-type metadata, gene
coordinates, motif profiles and the reference genome also require documented sources and
matching digests before a raw-input reconstruction can be accepted.

Model inference requires checkpoint-specific environments, identities and input contracts.
For PBMC, modality conversion and the scFoundation inference stage must finish before their
graph consumers run. The public summaries and replicate-level null archive do not contain
every upstream graph or confound cache. Deep randomization analyses use unrounded observed
statistics from those inputs; substituting rounded public statistics can change an
exceedance count near a tail threshold.

Resolution, degree-grid and injection analyses need verification at their own execution
level. Successful summary calculations establish arithmetic consistency of the released
records. They do not reproduce model inference, independent biological validation or
population uncertainty. Historical figure sources can be rendered, while the current
manuscript supplements supply the revised data-driven figures and the later statistical
analyses separately from the immutable archive.

A complete rerun record should identify the code revision, source inputs and digests,
installed environment, precomputed inputs, executed stages and exit status. Until that
record is available, the supported result is public-summary recalculation with an incomplete
raw-to-model reconstruction.
