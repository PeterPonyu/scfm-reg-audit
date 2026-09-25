# scReg-Eval

A reusable fixed-panel protocol and software capsule for auditing regulatory alignment in
single-cell RNA foundation-model gene graphs.

**[Numerical archive](https://doi.org/10.5281/zenodo.21724336)**

## What the protocol delivers

- **A frozen audit panel.** 1,200 genes, including the 446 transcription factors with a JASPAR
  motif that pass the panel filters, fixed before any sweep and shipped as a manifest so a rerun
  uses the same edge set (filter cascade: `results/revision_05613/panel_cascade.json`).
- **Expression-independent edge weights.** Edge weights come from a sequence- and
  accessibility-derived regulatory-potential proxy without co-expression weights. Panel
  selection and peak availability remain expression-adjacent design choices; the reference
  is not proven biologically independent of expression.
- **Confound control on that edge set.** Co-expression, target peak count, gene length, GC
  content, RNA detection rate, transcription-factor out-degree, and target in-degree.
- **Two structurally different randomizations.** A gene-label permutation and a within-row,
  TF-row-preserving target shuffle (target in-degree is not preserved), each a plus-one
  Monte Carlo at N=999. The original Dual-null Support screen intersects separate
  corrections within eight component families. Conjunction error control instead requires
  valid rowwise max-p values followed by a suitable correction within tissue × specification.
- **Separate reporting gates.** A Support census and a five-gate protocol-pass are reported
  as distinct outputs, so residual unusual alignment and a passable recovery claim never
  collapse into one verdict.
- **A prespecified degree sensitivity**, so a reader can tell how much of a Support result is
  conditional on degree covariates.
- **Two worked instances**, brain snATAC+RNA and paired PBMC multiome, with the constants
  frozen before the sweep and archived numerical outputs with analysis source code.

The panel, randomizations and gates are specified so that a new tissue, panel or model can be
scored without changing the protocol. So far the protocol has been run on one panel and two
tissues.

## Contents

- `src/` — audit implementations and automated tests; `src/revision_05613/` holds the analyses
  added in the first review round (see below).
- `data/manifest/` — the frozen shared-gene panel.
- `results/*.public.json` — public reference outputs used by validation and the examples;
  `results/revision_05613/` holds the revision outputs.
- `requirements.txt` — pinned Python versions of the reference environment.
- `ENVIRONMENT.example` — optional runtime-path configuration.
- `LICENSING.md` — code is MIT, results and documentation are CC BY 4.0, and vendored model
  code keeps its upstream license.

This tree carries the code and the public outputs. The Zenodo archive linked above adds the
protocol and analysis-plan documents and the full-rerun recipe with dataset and checkpoint
hashes (`docs/`; see the remaining input gaps below), the figure sources (`figures/`), and the SHA-256 digest of every file
(`MANIFEST.json`, `SHA256SUMS`). Datasets, model weights, and caches stay with their original
sources.

## Independent rerun scope

This repository runs independently of other paper repositories: inputs must be project-local
or explicitly configured. Upstream datasets, reference genomes and pretrained weights remain
scientific inputs with their original attribution and licenses.

- **Archive inspection:** public JSON outputs and the panel are in this tree. The immutable
  v0.5.1 Zenodo record additionally contains documents, frozen figure sources and a separate
  archive of replicate-level null statistics. A digest check validates bytes, not a rerun.
- **Validated summary checks:** the commands below recompute BH/BY/Holm and the existing
  per-row intersection–union corrections from released N=999 p-values, plus analytic MDE
  from released null standard deviations. Fixture tests compare these values with the
  corresponding archived results. No graph cache, raw NPZ or model inference is needed.
  These checks do not regenerate the p-values or establish the null models' validity.
- **Incomplete full rerun:** new null generation, graph-based Westfall–Young and model/proxy
  reconstruction need additional derived inputs, graph caches and upstream data/weights.
  Exact brain RNA preprocessing provenance and several input pins/conversion steps remain
  incomplete. No clean raw-to-model-to-statistics rerun is claimed.
- **Figures and new manuscript supplements:** rendering archived fixed TikZ sources is
  distinct from rebuilding figures from result JSON. The revised Figure 1, new Figures 14/15
  and additional manuscript intersection–union source package are not retroactively part of
  v0.5.1. Their final manifest and availability belong to the accompanying manuscript package.

```bash
python src/revision_05613/multiplicity.py --summary-only
python src/revision_05613/mde.py --summary-only
```

Both commands write to `output/revision_checks/`, leaving archived JSON unchanged. See
[Independent rerun instructions](docs/INDEPENDENT_RERUN.md) for input precedence, conversion
commands, validation scope and the remaining gaps. The archived `FULL_RERUN.md` is a historical
recipe; it is not evidence that those gaps are closed.

## Checks

```bash
pip install -r requirements.txt
python -m unittest discover -s src/tests
python -m unittest discover -s src/v2/tests
python validate_artifacts.py
```

`validate_artifacts.py` checks the public boundary: no manuscript or submission folders, no
private paths or keys, and only finite numbers in the JSON outputs. In the Zenodo archive,
`sha256sum -c SHA256SUMS` checks every file against its digest. Tests that need external data
run when their documented inputs are present.

## Revision analyses (ARRAY-D-26-05613, first review round)

The table identifies the archived analysis outputs, not a claim that each can be regenerated
from this checkout alone. The two summary commands above write separate check outputs.

| Script in `src/revision_05613/` | Output in `results/revision_05613/` | What it answers |
| --- | --- | --- |
| `resolution_check.py` | `resolution_n9999/` | 52 FM component, four original-baseline and eleven probe tests rerun at N = 9,999; their prefixes reproduce the 67 specified p-values; thirteen baseline-contrast tests remain at N = 999 |
| `multiplicity.py` | `multiplicity_robustness.json` | dual-null Support under BH, Benjamini–Yekutieli, Holm, intersection–union and Westfall–Young |
| `degree_grid.py` | `degree_grid_n999/` | the 13 rows under ten ways of controlling graph degree |
| `mde.py`, `power_sim.py` | `mde.json`, `power_sim_n999_r20/` | marginal normal-shift sensitivity approximations and a planted-signal joint-rejection check; no effect-size exclusion bound |
| `panel_cascade.py` | `panel_cascade.json` | the filter cascade behind the 1,200-gene panel |
| `motif_expectation.py` | `motif_expectation.json` | the expected number of random motif hits per peak |
| `brain_baseline_fix.py` | `brain_baseline_fix.json` | brain co-expression baseline with GC from brain peaks |
| `attention_guard.py` | `attention_guard.json` | deterministic rerun of the attention-omission guard |
| `make_appendix_b.py` | `appendix_b_tables.tex` | the manuscript's Appendix B tables, written from the JSON files |

Every audit test's seed, N, null mean, null standard deviation and exceedance count are in
`results/fixed_panel_audit_v2.public.json`; the revision analyses draw from the
`SeedSequence([20260924, k])` streams recorded in their outputs. Except for the explicit summary checks above, these analyses read cached gene
graphs, proxy graphs and/or covariate inputs of the full pipeline (`results/v2/`), which are not
redistributed here; the per-replicate null arrays of the revision runs are in the separate Zenodo archive.
Null arrays alone do not supply the unrounded observed statistics and graph-derived confounds
required by the current graph-based Westfall–Young entry point.

## Citation

See `CITATION.cff`, or cite the archive directly: <https://doi.org/10.5281/zenodo.21724336>.

## License

Original software is released under the MIT License.
