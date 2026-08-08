#!/usr/bin/env python
"""Confound-subset sensitivity for the pair-level TF-disjoint probe.

The primary readout residualises the probe's prediction and the ATAC-proxy target
against four gene-level confounds. Three of them (peakcount, genelen, gc) enter
the construction of the ATAC proxy itself, so the adjusted rho is sensitive to
which subset is conditioned on. This script reports the whole grid rather than a
single number, so the reported sign cannot be an artefact of one choice.
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import RidgeCV

ROOT = Path(__file__).resolve().parents[2]
PAIR_DIR = ROOT / "results/v2/tf_probe_pair"
EVAL = ROOT / "results/v2/tf_probe_pair_eval_v2.json"
OUT = ROOT / "results/v2/tf_probe_pair_sensitivity_v2.json"

ALPHAS = np.logspace(-3, 3, 25)
CONF_NAMES = ["peakcount", "genelen", "gc", "detv"]
SUBSETS = {
    "none": [],
    "atac_construction": [0, 1, 2],
    "full": [0, 1, 2, 3],
    "detv_only": [3],
}


def log(*values):
    print(f"[{time.strftime('%H:%M:%S')}]", *values, flush=True)


def adjusted_mean(y_test, pred, conf, cols):
    """Mean per-TF Spearman after residualising both sides on a confound subset."""
    n_tf, n_genes = y_test.shape
    if cols:
        c = conf[:, cols]
        cz = (c - c.mean(0)) / np.where(c.std(0) == 0, 1.0, c.std(0))
        design = np.column_stack([np.ones(n_genes), cz])
    else:
        design = np.ones((n_genes, 1))
    pinv = np.linalg.pinv(design)

    def resid(vec):
        r = rankdata(vec)
        return r - design @ (pinv @ r)

    out = np.full(n_tf, np.nan)
    for i in range(n_tf):
        a, b = resid(y_test[i]), resid(pred[i])
        if np.std(a) > 0 and np.std(b) > 0:
            out[i] = spearmanr(a, b).statistic
    return float(np.nanmean(out)), float(np.nanstd(out))


def main():
    """Write the family x confound-subset grid of adjusted rho."""
    ev = json.loads(EVAL.read_text())
    arm = ev["primary_arm"]
    families = sorted(ev["results"][arm])
    n_genes, n_tf = ev["n_genes"], ev["n_test_tfs"]

    prov = json.loads((PAIR_DIR / "provenance.json").read_text())
    cols = [prov["feature_names"].index(f) for f in ev["arms"][arm]]

    tg = np.load(PAIR_DIR / "pair_targets.npz", allow_pickle=False)
    y_train = tg["y_train"].astype(np.float64)
    y_test = tg["y_test"].astype(np.float64).reshape(n_tf, n_genes)
    conf = tg["conf_gene"].astype(np.float64)

    grid = {}
    for fam in families:
        d = np.load(PAIR_DIR / f"{fam}_pairs.npz", allow_pickle=False)
        xtr = d["X_train"].astype(np.float64)[:, cols]
        xte = d["X_test"].astype(np.float64)[:, cols]
        mu, sd = xtr.mean(0), xtr.std(0)
        sd = np.where(sd == 0, 1.0, sd)
        model = RidgeCV(alphas=ALPHAS, scoring="neg_mean_squared_error")
        model.fit((xtr - mu) / sd, y_train)
        pred = model.predict((xte - mu) / sd).reshape(n_tf, n_genes)

        row = {}
        for label, subset in SUBSETS.items():
            mean, std = adjusted_mean(y_test, pred, conf, subset)
            row[label] = {"adjusted_rho_mean": mean, "adjusted_rho_std": std}
        grid[fam] = row
        log(f"{fam:20s} " + " ".join(
            f"{lbl}={row[lbl]['adjusted_rho_mean']:+.4f}" for lbl in SUBSETS))

    out = {
        "schema_version": 1,
        "design": prov["design"],
        "arm": arm,
        "arm_features": ev["arms"][arm],
        "confound_names": CONF_NAMES,
        "subsets": {k: [CONF_NAMES[i] for i in v] for k, v in SUBSETS.items()},
        "note": (
            "peakcount/genelen/gc enter the ATAC-proxy construction, so "
            "`atac_construction` and `full` partly condition on the target's own "
            "definition; `none` is the unadjusted upper bound."
        ),
        "n_test_tfs": n_tf,
        "grid": grid,
    }
    OUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    log(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
