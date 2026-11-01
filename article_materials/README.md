# Fixed-panel audit of motif-accessibility alignment: analysis materials

This folder accompanies the article "A fixed-panel audit of motif-accessibility alignment in
single-cell RNA foundation-model gene graphs". It contains the code and saved numerical results
for the source-matched analyses reported in the article. The archive version of the folder adds
the larger per-draw arrays and retains the files of the earlier versions, including the historical
results shown in the appendix of the article.

## Contents

`corrective_pipeline/` prepares the genuine-brain (GSE174367) and PBMC contexts from the public
raw files, selects the gene panels, builds the proximal motif-accessibility proxy, runs the six
model readouts (Geneformer embedding and attention, random-initialized Geneformer, scGPT,
scFoundation and UCE; adapters and their licences are in `screg_validation/vendor`) and computes
the gene-label (M) and target-shuffle (D) randomization statistics. `corrective_plan.json` and
`plan.json` record the registered designs, `resources.example.json` lists the inputs the pipeline
expects, and `INPUT_CONTRACTS.md` describes the input formats and the local run record that the
download and inference steps read.

`registered_bone_marrow/` contains the registered bone-marrow analysis of evaluation donor 16710:
the panel, the preparation summary, the registered rows and, for each degree design and
transformation class, the observed rows and all 9,999 saved draws (`draw-progress.json`).

`resolution_9999/` contains the scripts and results of the post hoc 9,999-draw extension of the
registered genuine-brain analysis; each stream folder in `results/` holds its rows and all 9,999
saved draws (`draw-progress.json`).

`revision3/analysis/` contains the analyses added in the final revision: calibration of the
component tests under synthetic stress (`r3_calibration.py`, which imports the stress-test helpers
of the earlier version from `extend_injection_power.py`), the common-support analysis and the
symmetric co-expression comparison (`r3_common_support.py`), the ten degree designs
(`r3_grid.py`, with external proxies from `r3_external_proxies.py`), disjoint discovery and
evaluation samples (`r3_split.py`), the five-donor bone-marrow analysis (`r3_bmmc.py`), the exact
reproduction of the registered streams (`r3_reproduce_registered.py`), the weight check
(`r3_scf_verify.py`), aggregation and multiplicity adjustment (`r3_aggregate.py`), table
generation (`r3_tables.py`) and the figure code (`viz/`). `R3_PLAN.md` lists the decision rules,
which were fixed and checksummed before these analyses were run; addendum A records the decision
taken after the calibration result, and addendum B records the four calibration contexts added
after the family results. Input locations are set in `paths.local.json` (schema in
`paths.example.json`) or by environment variables, as described in `r3paths.py`.

`revision3/evidence/` contains every reported statistic of the new analyses: `R3_RESULTS.json`
(families, adjusted values and summaries), one file per randomization stream with its observed
statistics, component probabilities and null moments, the reproduction checks against the
registered results (repeated with the final code in `reproduction/final_code/`), and the weight
check. In the archive version the same folders also hold the
saved reference distributions of the calibration and common-support streams, the external proxy
graphs, and the proxy, co-expression and model graphs of the evaluation halves (`.npz`). The reference distributions of the
degree-design, sample-half and five-donor bone-marrow streams were not saved; rerunning these
streams with their fixed seeds regenerates them.

`revision3/inputs/` and `revision3/historical_evidence/` contain the inputs of the earlier analyses
that the table and figure code reads (CTCF labels, context definitions, exposure levels and the
plotted data of the historical comparisons), with only the fields the code uses.
`revision3/corrective_run/` holds the panel list here and, in the archive version, the gene-level
graphs of the corrective run (proxy, co-expression and the six model readouts for each context;
no cell-level data).

`revision3/figures/` contains the figures of the article.

## Reproducing the results

Python 3.12 or later with NumPy, SciPy, h5py, pandas, scikit-learn, PyTorch and MOODS is
required. The public inputs are GSE174367 (snRNA-seq counts and metadata, snATAC-seq peak
matrix and metadata), the 10x Genomics 10k PBMC multiome reference, GSE194122 (processed
bone-marrow multiome), GSE206767 (scATAC peak matrix), the JASPAR 2024 CORE vertebrate motifs,
the hg38 genome and the pretrained weights from the model developers. The corrective pipeline
turns these into the per-context panels, proxies and model scores that the new analyses read
(`corrective_run`); the calibration of the two historical contexts and the appendix figures read
the historical inputs of the earlier versions of the archive. The inputs can be obtained with
`corrective_pipeline/run.py`, which downloads with the Python standard library (optionally through
the proxy named in `SC_FETCH_PROXY`) after checking each file's length against its registered size and, where one is registered, its
checksum, or from the sources listed in `resources.example.json`.

Seeds are fixed in the code: 20261001 for the corrective, common-support, degree-design and
sample-half randomization streams; 20260930 for the bone-marrow streams and donor intervals;
20261031 for the sample-half assignment, the TF-cluster bootstrap and the synthetic readouts of
the five analysis contexts; and the published synthetic seeds for the two historical contexts.
Every analysis script writes the files found in `revision3/evidence`. With the archive version
unpacked, the seven new tables and the eight figures generated by the revision figure code can be
rebuilt from the saved inputs, apart from the Arial font files (set `font_dir` or `ARIAL_FONT_DIR`):

    python revision3/analysis/r3_tables.py
    cd revision3/analysis && python -m viz.build_r3

This rebuild was tested on a clean copy of the prepared archive: the seven generated tables matched
the article's source tables byte for byte, and the eight generated figures matched the article's
figure PDFs when rendered at 144 dpi. Three historical appendix figures are included as PDFs;
the revision figure command does not regenerate them. The aggregation script can also be run
separately from the saved inputs.

## Notes

Files are copied unchanged except those listed in `publication_edits.json`. In those, paths are
shortened to file names, machine names are replaced, local approval notes are removed, the run-record
messages are worded for general use and the download step no longer requires a local gateway; the
list records the SHA-256 of each original file, and `resolution_9999/scripts/verify_resolution.py`
accepts these run-time digests (it prints its report; `--out` writes the report files to a folder). No statistic, draw or numerical result is changed.

Pretrained weights are not redistributed. The evaluated scFoundation file is identified by its
SHA-256 in `revision3/evidence/scfoundation_verification.json`, which also records the published
checksum of the mirror from which it was obtained.
