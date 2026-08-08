#!/usr/bin/env python
"""Train and evaluate the TF-disjoint pair-level probe.

One RidgeCV per family, fitted on train-TF pairs and applied to held-out test
TF pairs. Readouts per test TF: marginal Spearman, confound-adjusted Spearman,
and (for FM families) partial Spearman conditioning on the co-expression
probe's prediction for the same pair.
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import RidgeCV

ROOT = Path(__file__).resolve().parents[2]
PAIR_DIR = ROOT / "results/v2/tf_probe_pair"
OUT = ROOT / "results/v2/tf_probe_pair_eval_v2.json"

ALPHAS = np.logspace(-3, 3, 25)
BASELINE = "co_expression"

# Feature-set arms. `all` lets the probe use gene-level degree columns, which are
# themselves confound-correlated and let it shortcut past the edge weights;
# `edge_only` withholds them so the readout reflects the graph's own edges.
ARMS = {
    "all": [0, 1, 2, 3, 4],
    "edge_only": [0, 1, 4],
}
PRIMARY_ARM = "edge_only"


def log(*values):
    print(f"[{time.strftime('%H:%M:%S')}]", *values, flush=True)


def residualise(vec, design):
    """Rank-residualise vec against an intercept-augmented design matrix."""
    r = rankdata(vec)
    coef, *_ = np.linalg.lstsq(design, r, rcond=None)
    return r - design @ coef


def per_tf_rho(y_true, y_pred, n_tf, n_genes, design=None):
    """Spearman per test TF, optionally after residualising both sides."""
    yt = y_true.reshape(n_tf, n_genes)
    yp = y_pred.reshape(n_tf, n_genes)
    out = []
    for i in range(n_tf):
        a, b = yt[i], yp[i]
        if design is not None:
            a, b = residualise(a, design), residualise(b, design)
        if np.std(b) == 0 or np.std(a) == 0:
            out.append(np.nan)
            continue
        out.append(spearmanr(a, b).statistic)
    return np.asarray(out, dtype=float)


def partial_rho(y_true, y_pred, y_ctrl, n_tf, n_genes):
    """Partial Spearman of (y_true, y_pred) given the baseline prediction."""
    yt = y_true.reshape(n_tf, n_genes)
    yp = y_pred.reshape(n_tf, n_genes)
    yc = y_ctrl.reshape(n_tf, n_genes)
    out = []
    for i in range(n_tf):
        design = np.column_stack([np.ones(n_genes), rankdata(yc[i])])
        a = residualise(yt[i], design)
        b = residualise(yp[i], design)
        if np.std(a) == 0 or np.std(b) == 0:
            out.append(np.nan)
            continue
        out.append(spearmanr(a, b).statistic)
    return np.asarray(out, dtype=float)


def summarise(name, rhos):
    """Mean/std/median of a per-TF rho vector, ignoring degenerate TFs."""
    ok = rhos[np.isfinite(rhos)]
    return {
        f"{name}_mean": float(np.mean(ok)),
        f"{name}_std": float(np.std(ok)),
        f"{name}_median": float(np.median(ok)),
        f"{name}_n_ok": int(ok.size),
    }


def main():
    """Fit every family's probe and write the evaluation JSON."""
    prov = json.loads((PAIR_DIR / "provenance.json").read_text())
    families = prov["families"]
    n_genes = prov["n_genes"]
    n_test_tf = prov["n_test_tfs"]

    tg = np.load(PAIR_DIR / "pair_targets.npz", allow_pickle=False)
    y_train = tg["y_train"].astype(np.float64)
    y_test = tg["y_test"].astype(np.float64)
    conf_gene = tg["conf_gene"].astype(np.float64)

    # Confound design replicated across test TFs, standardised, with intercept.
    cz = (conf_gene - conf_gene.mean(0)) / np.where(conf_gene.std(0) == 0, 1.0, conf_gene.std(0))
    design_conf = np.column_stack([np.ones(n_genes), cz])
    log(f"train pairs {y_train.shape[0]}, test pairs {y_test.shape[0]}, test TFs {n_test_tf}")

    results = {arm: {} for arm in ARMS}
    for arm, cols in ARMS.items():
        log(f"--- arm: {arm} ({len(cols)} features) ---")
        preds = {}
        for fam in families:
            d = np.load(PAIR_DIR / f"{fam}_pairs.npz", allow_pickle=False)
            names = [str(v) for v in d["feature_names"]]
            xtr = d["X_train"].astype(np.float64)[:, cols]
            xte = d["X_test"].astype(np.float64)[:, cols]
            mu, sd = xtr.mean(0), xtr.std(0)
            sd = np.where(sd == 0, 1.0, sd)
            model = RidgeCV(alphas=ALPHAS, scoring="neg_mean_squared_error")
            model.fit((xtr - mu) / sd, y_train)
            preds[fam] = model.predict((xte - mu) / sd)

            marg = per_tf_rho(y_test, preds[fam], n_test_tf, n_genes)
            adj = per_tf_rho(y_test, preds[fam], n_test_tf, n_genes, design=design_conf)
            results[arm][fam] = {
                "alpha": float(model.alpha_),
                "coef": {names[c]: float(v) for c, v in zip(cols, model.coef_)},
                "train_r2": float(model.score((xtr - mu) / sd, y_train)),
                **summarise("marginal_rho", marg),
                **summarise("adjusted_rho", adj),
            }
            log(f"  {fam:20s} alpha={model.alpha_:8.3f} "
                f"marginal={results[arm][fam]['marginal_rho_mean']:+.4f} "
                f"adjusted={results[arm][fam]['adjusted_rho_mean']:+.4f}")

        for fam in families:
            if fam == BASELINE:
                continue
            part = partial_rho(y_test, preds[fam], preds[BASELINE], n_test_tf, n_genes)
            results[arm][fam].update(summarise("partial_rho_given_coexp", part))
            for key in ("marginal_rho", "adjusted_rho"):
                results[arm][fam][f"delta_{key}_vs_coexp"] = float(
                    results[arm][fam][f"{key}_mean"] - results[arm][BASELINE][f"{key}_mean"]
                )

    evaluation = {
        "schema_version": 1,
        "design": prov["design"],
        "baseline": BASELINE,
        "primary_arm": PRIMARY_ARM,
        "arms": {k: [prov["feature_names"][c] for c in v] for k, v in ARMS.items()},
        "n_genes": n_genes,
        "n_test_tfs": n_test_tf,
        "n_train_tfs": prov["n_train_tfs"],
        "split_seed": prov["split_seed"],
        "confounds": ["peakcount", "genelen", "gc", "detv"],
        "results": results,
    }
    OUT.write_text(json.dumps(evaluation, indent=2, sort_keys=True) + "\n")
    log(f"wrote {OUT.relative_to(ROOT)}")

    log(f"primary arm ({PRIMARY_ARM}) ranked by confound-adjusted rho")
    for fam, r in sorted(results[PRIMARY_ARM].items(), key=lambda kv: -kv[1]["adjusted_rho_mean"]):
        log(f"  {fam:20s} adjusted={r['adjusted_rho_mean']:+.4f} "
            f"marginal={r['marginal_rho_mean']:+.4f}")


if __name__ == "__main__":
    main()
