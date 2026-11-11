# Rechecking the historical numerical summaries

These instructions cover the original fixed-panel analysis. Its brain-labelled RNA was later
identified as ascending aorta, and its brain interpretation is withdrawn. The source-corrected
pipeline and Revision 3 results are described in [article_materials](../article_materials/README.md).

The commands below recalculate statistics from released summaries. A full reconstruction also
requires the original upstream inputs, as described below.

## Recalculate the published summaries

Use the declared Python environment (`requirements.txt`; reference Python 3.13.7):

```bash
python -m pip install -r requirements.txt
python src/revision_05613/multiplicity.py --summary-only
python src/revision_05613/mde.py --summary-only
python -m unittest discover -s src/v2/tests -p test_independent_release.py -v
```

These commands were checked in an existing environment. Installation in a fresh environment
was not tested.

The summary commands write respectively:

- `output/revision_checks/multiplicity_robustness.json`
- `output/revision_checks/mde.json`

Use `--output /path/to/check.json` to choose a separate output location. The frozen JSON in
`results/revision_05613/` is not overwritten by default. Input selection is explicit configured
file first, then `results/v2/<name>.json`, then `results/<name>.public.json`. An invalid
explicit configuration fails instead of silently substituting another file. To select the published
inputs even when other cached outputs are present, use:

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

## Input paths

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

This example selects inputs for the brain ATAC proxy. Check their identities against the
recorded inputs before using it to reproduce the analysis:

```bash
ATAC_FILE=/path/to/brain_atac.h5ad \
META_FILE=/path/to/brain_cell_metadata.csv.gz \
TAG=GSE174367 python src/v2/build_atac_graph_v2.py
```

Missing PBMC or gene-coordinate inputs produce an error. The scripts do not search other
research directories for replacements.

## Remaining requirements for a full rerun

The original RNA source was identified as ascending aorta. The historical summary checks do
not reconstruct the original cell selection, donor mapping or preprocessing. A full rerun also
requires documented cell-type metadata, gene coordinates, motifs and reference-genome files
with matching digests. The source-corrected pipeline has its own inputs and instructions in
`article_materials/`.

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
installed environment, precomputed inputs, executed stages and exit status. Without that
record, these commands establish summary-level numerical agreement only.

