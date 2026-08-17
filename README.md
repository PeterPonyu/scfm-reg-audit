# scReg-Eval fixed-panel audit capsule v0.4.0

This is a sanitized **audit capsule** for the scReg-Eval manuscript. It validates the published
numerical artifacts and carries the audited statistical code; it is not a raw-data reproduction
environment.

## Included

- fixed-panel statistical implementation, PBMC UCE pipeline, baseline and probe scripts, and tests;
- frozen 446-TF / 1,200-gene manifest;
- path-scrubbed authoritative JSON derivatives (pooled audit, injection, inference status,
  co-expression baselines, TF-disjoint probe, cell-type invariance, subdivided calibration);
- manuscript source, active data-driven figure fragments, and the figure generator;
- licenses (MIT for original code, CC BY 4.0 for manuscript/figures/derived results), citation
  metadata, and this validator with a SHA-256 manifest.

## Excluded

Raw/processed H5AD/H5 data, genome FASTA, model weights, cached NPZ graphs, model-inference/vendor
code, pilot code, session state, build logs/caches, internal review records, and the four retired
numerical JSON files are excluded. See `docs/PAPER_OUTLINE.md` and `docs/LEGACY_INFERENCE_NOTE.md`.

## Validate

```bash
python validate_artifacts.py
```

## Test

```bash
pip install numpy scipy anndata pyfaidx  # test dependencies
python -m unittest src.tests.test_fixed_panel_audit
```

The bundled suite is the statistical-contract suite (`test_fixed_panel_audit.py`). Legacy-hash,
model-scope, and real-cache integration cases skip automatically when their retired or
machine-local inputs are absent; the development repository runs the full suite including
model/cache contract tests. The release builder runs this command after every rebuild and fails
on any error.

## Reproduction boundary

The capsule validates published numerical artifacts. Re-running production requires external public
datasets and cached model/readout graphs not redistributed here. See `docs/FULL_RERUN.md` in the
development repository for the full-rerun recipe. Set `SCFM_BRAIN_ATAC` to the brain
ATAC H5AD path. Provenance strings use `${SCFM_PROJECT_ROOT}`, `${SCFM_DATA_ROOT}`, and
`${SCFM_BRAIN_ATAC}` placeholders.
