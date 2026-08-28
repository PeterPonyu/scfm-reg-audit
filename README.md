# scReg-Eval

A reusable fixed-panel protocol and software capsule for auditing regulatory alignment in
single-cell RNA foundation-model gene graphs.

**[Numerical archive](https://doi.org/10.5281/zenodo.21724336)** ·
**[Companion capsule: tumor locked-proxy transfer](https://github.com/PeterPonyu/locked-proxy-transfer-tumor-chromatin)**

## What the protocol delivers

- **A frozen audit panel.** 1,200 genes, including 446 in-range transcription factors, fixed
  before any sweep and shipped as a manifest so a rerun uses the same edge set.
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

The same panel, randomizations, and gates attach to a new tissue, a new panel, or a new model
without changing the protocol.

## Contents

- `src/` — audit implementations and automated tests.
- `data/manifest/` — the frozen shared-gene panel.
- `results/*.public.json` — public reference outputs used by validation and the examples.
- `requirements.txt` — pinned Python versions of the reference environment.
- `ENVIRONMENT.example` — optional runtime-path configuration.

This tree carries the code and the public inputs. The full numerical record, including the
per-row outputs and content digests, lives in the Zenodo archive linked above; datasets, model
weights, and caches stay with their original sources.

## Checks

```bash
pip install -r requirements.txt
python -m unittest discover -s src/tests
python -m unittest discover -s src/v2/tests
python validate_artifacts.py
```

`validate_artifacts.py` revalidates the capsule against the MANIFEST digests. Tests that need
external data run when their documented inputs are present.

## Citation

See `CITATION.cff`, or cite the archive directly: <https://doi.org/10.5281/zenodo.21724336>.

## License

Original software is released under the MIT License.
