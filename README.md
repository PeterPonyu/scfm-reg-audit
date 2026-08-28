# scReg-Eval code capsule

This repository contains executable source code and small, path-scrubbed public inputs for fixed-panel regulatory graph audits. It is the public implementation wall, not a manuscript host.

## Contents

- `src/`: audit implementations and automated tests.
- `data/manifest/`: the frozen shared-gene panel.
- `results/*.public.json`: small public reference outputs used by validation and examples.
- `requirements.txt`: pinned Python versions of the reference environment.
- `ENVIRONMENT.example`: optional runtime-path configuration.

Large datasets, model weights, caches, journal submission packages, internal review material, and local execution state are excluded.

Independent public full-text walls (compiled PDF plus source):

- [scReg-Eval fixed-panel manuscript](https://github.com/PeterPonyu/scfm-reg-paper-a-fixed-panel)
- [Tumor locked-proxy transfer](https://github.com/PeterPonyu/locked-proxy-transfer-tumor-chromatin)
- [Cross-organ chromatin maps](https://github.com/PeterPonyu/scfm-reg-paper-u-organ-transfer)

See `REPO_TOPOLOGY.md`.

## Checks

```bash
pip install -r requirements.txt
python -m unittest discover -s src/tests
python -m unittest discover -s src/v2/tests
python validate_artifacts.py
```

External-data tests skip when their documented inputs are unavailable.

## Publication boundary

GitHub Pages and GitHub Actions are not used for this capsule. Each upcoming manuscript has its own public repository; venue upload decks stay local.

## License

Original software is released under the MIT License.
