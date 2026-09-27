# scReg-Eval

scReg-Eval audits the alignment of gene graphs from single-cell RNA foundation models with a reference built from chromatin accessibility and sequence motifs. The reference measures regulatory potential rather than causal regulation. The demonstration uses a fixed panel of 1,200 genes, including 446 transcription factors, with brain and peripheral blood mononuclear cell data.

The analysis separates residual rank alignment, comparisons with co-expression, and an operational recovery criterion. It adjusts for measured expression and structural covariates and evaluates two randomizations. Seven of thirteen readouts pass the original dual-null screen, including six positive and one negative association; none meets all five recovery criteria. These results apply to the available frozen graphs and do not establish a general ranking of model capability.

The repository is independent of other research repositories. Inputs must be local to this project or explicitly supplied. Released numerical summaries support independent recalculation of the documented multiple-testing and sensitivity summaries. A complete reconstruction from raw data through model inference has not been verified: the exact brain RNA preprocessing and metadata provenance remain incomplete. File validation and successful figure rendering do not close that gap.

The [numerical archive](https://doi.org/10.5281/zenodo.21724336) preserves the version history. [Version 0.5.1](https://doi.org/10.5281/zenodo.22958412) includes the underlying review analyses, protocol, analysis plan, historical execution instructions, figure sources and file digests. Subsequent manuscript figure and statistical supplements are separate from that immutable archive. Public datasets and pretrained weights remain with their original providers.

[Independent rerun instructions](docs/INDEPENDENT_RERUN.md) describe the available summary checks, required inputs and unresolved reconstruction steps. The source, public reference results and pinned environment specification are included here. The [citation record](CITATION.cff) identifies the software and archive.

Original software uses the MIT License. Results and documentation use CC BY 4.0; included third-party code retains its upstream licenses, as described in [Licensing](LICENSING.md).
