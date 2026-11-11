# scReg-Eval

scReg-Eval tests whether gene graphs from single-cell RNA foundation models align with a chromatin-accessibility and sequence-motif reference after accounting for co-expression and graph structure. The reference measures proximal regulatory potential, not causal regulation.

The code and results for Revision 3 of "A fixed-panel audit of motif-accessibility alignment in single-cell RNA foundation-model gene graphs" are in [`article_materials/`](article_materials/README.md). This folder contains the source-corrected pipeline, analysis scripts and numerical summaries. The [Zenodo archive](https://doi.org/10.5281/zenodo.21724336) also contains the larger saved randomization arrays and gene graphs. The folder README lists which arrays were saved and which analyses require a rerun.

The RNA labelled as brain in the original analysis came from ascending aorta, GSE213740. We withdrew its brain interpretation. The corrected analyses use brain data from GSE174367, PBMC and bone marrow from GSE194122. The comparisons did not establish stronger reference alignment for any tested zero-shot readout than for co-expression from the same cells.

The other folders retain the original 1,200-gene, 446-transcription-factor analysis. Seven of its thirteen readouts passed the original dual-null screen, with six positive and one negative association; none met all five operational criteria. These are historical counts, and the brain interpretation is withdrawn. The [rerun instructions](docs/INDEPENDENT_RERUN.md) explain how to recalculate the historical multiple-testing and sensitivity summaries.

[CITATION.cff](CITATION.cff) identifies the software and archive. Public datasets and pretrained weights are available from their original providers.

Original software is licensed under MIT. Results and documentation use CC BY 4.0. Third-party code retains its upstream licences; see [LICENSING.md](LICENSING.md).
