# scReg-Eval

A reusable fixed-panel protocol and software capsule for auditing regulatory alignment in
single-cell RNA foundation-model gene graphs.

**[Numerical archive](https://doi.org/10.5281/zenodo.21724336)** ·
**[Companion capsule: tumor locked-proxy transfer](https://github.com/PeterPonyu/locked-proxy-transfer-tumor-chromatin)**

## What the protocol delivers

- **A frozen audit panel.** 1,200 genes, including the 446 transcription factors with a JASPAR
  motif that pass the panel filters, fixed before any sweep and shipped as a manifest so a rerun
  uses the same edge set (filter cascade: `results/revision_05613/panel_cascade.json`).
- **An expression-free reference.** Edge weights come from a sequence- and
  accessibility-derived regulatory-potential proxy, so the reference does not inherit
  co-expression structure from the models under audit.
- **Confound control on that edge set.** Co-expression, target peak count, gene length, GC
  content, RNA detection rate, transcription-factor out-degree, and target in-degree.
- **Two structurally different randomizations.** A gene-label permutation and a within-row,
  degree-preserving target shuffle, each a plus-one Monte Carlo at N=999, read together as
  Dual-null Support and corrected within eight prespecified families.
- **Separate reporting gates.** A Support census and a five-gate protocol-pass are reported
  as distinct outputs, so residual unusual alignment and a passable recovery claim never
  collapse into one verdict.
- **A prespecified degree sensitivity**, so a reader can tell how much of a Support result is
  conditional on degree covariates.
- **Two worked instances**, brain snATAC+RNA and paired PBMC multiome, with the constants
  frozen before the sweep and the scripts that render every reported value.

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
hashes (`docs/`), the figure sources (`figures/`), and the SHA-256 digest of every file
(`MANIFEST.json`, `SHA256SUMS`). Datasets, model weights, and caches stay with their original
sources.

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

| Script in `src/revision_05613/` | Output in `results/revision_05613/` | What it answers |
| --- | --- | --- |
| `resolution_check.py` | `resolution_n9999/` | every test rerun at N = 9,999 with its own seed; the first 999 replicates reproduce all 67 published p-values |
| `multiplicity.py` | `multiplicity_robustness.json` | dual-null Support under BH, Benjamini–Yekutieli, Holm, intersection–union and Westfall–Young |
| `degree_grid.py` | `degree_grid_n999/` | the 13 rows under ten ways of controlling graph degree |
| `mde.py`, `power_sim.py` | `mde.json`, `power_sim_n999_r20/` | minimum detectable effect of each test and a planted-signal check |
| `panel_cascade.py` | `panel_cascade.json` | the filter cascade behind the 1,200-gene panel |
| `motif_expectation.py` | `motif_expectation.json` | the expected number of random motif hits per peak |
| `brain_baseline_fix.py` | `brain_baseline_fix.json` | brain co-expression baseline with GC from brain peaks |
| `attention_guard.py` | `attention_guard.json` | deterministic rerun of the attention-omission guard |
| `make_appendix_b.py` | `appendix_b_tables.tex` | the manuscript's Appendix B tables, written from the JSON files |

Every audit test's seed, N, null mean, null standard deviation and exceedance count are in
`results/fixed_panel_audit_v2.public.json`; the revision analyses draw from the
`SeedSequence([20260924, k])` streams recorded in their outputs. The scripts read the cached gene
graphs, proxy graphs and covariate inputs of the full pipeline (`results/v2/`), which are not
redistributed here; the per-replicate null arrays of the revision runs are in the Zenodo archive.

## Citation

See `CITATION.cff`, or cite the archive directly: <https://doi.org/10.5281/zenodo.21724336>.

## License

Original software is released under the MIT License.
