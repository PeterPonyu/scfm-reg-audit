# scReg-Eval code capsule

This repository contains executable source code and small, path-scrubbed public inputs for fixed-panel regulatory graph audits. It is a software release surface only.

## Contents

- `src/`: audit implementations and automated tests.
- `data/manifest/`: the frozen shared-gene panel.
- `results/*.public.json`: small public reference outputs used by validation and examples.
- `requirements.txt`: pinned Python versions of the reference environment.
- `ENVIRONMENT.example`: optional runtime-path configuration.

Large datasets, model weights, caches, manuscript sources, figures, submission packages, internal review material, and local execution state are intentionally excluded.

## Checks

```bash
pip install -r requirements.txt
python -m unittest discover -s src/tests
python -m unittest discover -s src/v2/tests
python validate_artifacts.py
```

External-data tests skip when their documented inputs are unavailable.

## Publication boundary

GitHub Pages and GitHub Actions are not used for this capsule. The repository tree is the complete public surface.

## License

Original software is released under the MIT License.
