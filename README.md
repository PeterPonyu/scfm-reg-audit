# scReg-Eval

scReg-Eval audits whether gene graphs extracted from single-cell RNA foundation models align with a reference built from chromatin accessibility and sequence motifs beyond what co-expression and graph structure already explain. The reference measures regulatory potential rather than causal regulation.

The materials for the article "A fixed-panel audit of motif-accessibility alignment in single-cell RNA foundation-model gene graphs" are in [`article_materials/`](article_materials/README.md). They contain the corrective pipeline from the public raw inputs through panel selection, proxy construction, model inference and randomization; the saved randomization statistics of every reported analysis, including all 9,999-draw reference distributions of the post hoc extension and of the registered bone-marrow analysis; and the code of the analyses added in the final revision.

Source verification showed that the RNA used as brain in the first version of this analysis came from ascending aorta (GSE213740). That brain interpretation is withdrawn. The source-matched analyses use genuine brain (GSE174367), PBMC and bone marrow (GSE194122); in them no tested zero-shot readout is demonstrably more aligned with the reference than co-expression from the same cells.

The rest of this repository records the original analysis, which used a fixed panel of 1,200 genes including 446 transcription factors. In that analysis seven of thirteen readouts passed the original dual-null screen, including six positive and one negative association, and none met all five operational criteria; its brain rows rest on the misidentified RNA. Its released numerical summaries support independent recalculation of the documented multiple-testing and sensitivity summaries.

The [numerical archive](https://doi.org/10.5281/zenodo.21724336) preserves the version history, and its newest version adds the article materials. [Independent rerun instructions](docs/INDEPENDENT_RERUN.md) describe the summary checks of the original analysis. The [citation record](CITATION.cff) identifies the software and archive. Public datasets and pretrained weights remain with their original providers.

Original software uses the MIT License. Results and documentation use CC BY 4.0; included third-party code retains its upstream licenses, as described in [Licensing](LICENSING.md).
