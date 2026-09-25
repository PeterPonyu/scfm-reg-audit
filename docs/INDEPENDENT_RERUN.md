# Independent rerun instructions

This document describes the bounded independence patch for ARRAY-D-26-05613. It does not
change or expand the immutable v0.5.1 archive. No other paper repository is an implicit input.
External scientific inputs and their original citations/licenses remain necessary.
Historical registry keys such as `peerj_freeze` remain for schema compatibility; they
label the frozen panel and do not connect to another manuscript or repository.

## Validated public-summary layer

Use the declared Python environment (`requirements.txt`; reference Python 3.13.7):

```bash
python -m pip install -r requirements.txt
python src/revision_05613/multiplicity.py --summary-only
python src/revision_05613/mde.py --summary-only
python -m unittest discover -s src/v2/tests -p test_independent_release.py -v
```

The validation performed for this patch used an existing local environment, not a newly
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
model inference were not executed for this patch. Repeated successful conversion into the
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

The historical `SCREG_MONOREPO_DATA` sibling fallback has been removed. Merely having a
neighboring research workspace cannot satisfy missing PBMC or coordinates inputs.

## Remaining requirements for a full rerun

The archived `docs/FULL_RERUN.md` supplies historical context but does not close all of these
requirements. Do not infer a missing accession or a deterministic converter from a filename.

1. Supply the exact public brain RNA accession(s), subset and deterministic preprocessing
   that produced the `ad_hm_prepped` input. This patch does not invent the missing provenance.
2. Provide the cell-type metadata source/converter and digest; pin gene coordinates, JASPAR
   motif file and the exact hg38 FASTA build/digest. Keep attribution and upstream licenses.
3. Complete the ordered model-specific environment, checkpoint and input contracts. These
   include PBMC conversion and `pbmc_eval_scf_v2.py` before consumers of
   `G_scf_pbmc_pooled.npz`; optional model packages are not installed by this patch.
4. Reconstruct and verify the graph/confound caches in `results/v2/`. The public JSON and
   replicate-level null archive do not contain every input used by `common.load_groups()`.
   Deep multiplicity currently reconstructs unrounded observed statistics from these caches;
   the rounded public statistic is not an interchangeable substitute near a tail threshold.
5. Validate resolution, degree-grid, power and other graph-dependent analyses separately.
   Passing the public-summary fixture tests does not establish these stages ran.
6. Close the data-to-figure chain. The frozen TikZ sources can be rendered; the historical
   R figure driver still references a different path layout and is not a validated renderer
   of the newly revised manuscript figures. Figure 1/14/15 and new IUT supplements belong
   to the final accompanying manuscript package, not retrospectively to v0.5.1.

A future full-run record should identify its code revision, pinned inputs and hashes, installed
environment, use of precomputed caches, executed stages and exit status. Until that record
exists, the supported wording is **validated public-summary recalculation with incomplete
raw-to-model reconstruction**, not an independently reproduced full analysis.
