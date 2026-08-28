# scReg-Eval code capsule

This repository contains executable source code and small, path-scrubbed public inputs for fixed-panel regulatory graph audits.

## Contents

- `src/`: audit implementations and automated tests.
- `data/manifest/`: the frozen shared-gene panel.
- `results/*.public.json`: small public reference outputs used by validation and examples.
- `requirements.txt`: pinned Python versions of the reference environment.
- `ENVIRONMENT.example`: optional runtime-path configuration.

Large datasets, model weights, caches, and local execution state are excluded.

Companion code capsules:

- [Tumor locked-proxy transfer](https://github.com/PeterPonyu/locked-proxy-transfer-tumor-chromatin)

## Checks

```bash
pip install -r requirements.txt
python -m unittest discover -s src/tests
python -m unittest discover -s src/v2/tests
python validate_artifacts.py
```

External-data tests skip when their documented inputs are unavailable.

## License

Original software is released under the MIT License.
