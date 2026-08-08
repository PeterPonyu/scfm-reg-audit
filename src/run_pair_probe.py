#!/usr/bin/env python
"""Train and evaluate the TF-disjoint pair-level probe.

One RidgeCV per family, fitted on train-TF pairs and applied to held-out test
TF pairs. Readouts per test TF: marginal Spearman, confound-adjusted Spearman,
and (for FM families) partial Spearman conditioning on the co-expression
probe's prediction for the same pair.
"""
import json
import os
import sys

import numpy as np
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pair_probe_common as ppc  # noqa: E402
from pair_probe_common import BASELINE, PAIR_DIR, ROOT, blocks, log  # noqa: E402

OUT = ROOT / "results/v2/tf_probe_pair_eval_v2.json"

# Feature-set arms. `all` lets the probe use gene-level degree columns, which are
# themselves confound-correlated and let it shortcut past the edge weights;
# `edge_only` withholds them so the readout reflects the graph's own edges.
ARMS = {
    "all": [0, 1, 2, 3, 4],
    "edge_only": [0, 1, 4],
}
PRIMARY_ARM = "edge_only"


def partial_rho(y_true, y_pred, y_ctrl):
    """Partial Spearman of (y_true, y_pred) given the baseline prediction."""
    out = np.full(y_true.shape[0], np.nan)
    n_genes = y_true.shape[1]
    for i in range(y_true.shape[0]):
        design = np.column_stack([np.ones(n_genes), rankdata(y_ctrl[i])])
        residualise = ppc.rank_residualiser(design)
        a, b = residualise(y_true[i]), residualise(y_pred[i])
        if np.std(a) > 0 and np.std(b) > 0:
            out[i] = spearmanr(a, b).statistic
    return out


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
    y_test = blocks(tg["y_test"], n_test_tf, n_genes)

    # Confound design replicated across test TFs, standardised, with intercept.
    adjust = ppc.rank_residualiser(ppc.confound_design(tg["conf_gene"], n_genes))
    log(f"train pairs {y_train.shape[0]}, test pairs {y_test.size}, test TFs {n_test_tf}")

    results = {arm: {} for arm in ARMS}
    for arm, cols in ARMS.items():
        log(f"--- arm: {arm} ({len(cols)} features) ---")
        preds = {}
        for fam in families:
            fit = ppc.fit_family_probe(fam, cols, y_train)
            preds[fam] = blocks(fit.prediction, n_test_tf, n_genes)

            marg = ppc.per_tf_rho(y_test, preds[fam])
            adj = ppc.per_tf_rho(y_test, preds[fam], residualise=adjust)
            results[arm][fam] = {
                "alpha": fit.alpha,
                "coef": {fit.feature_names[c]: float(v) for c, v in zip(cols, fit.coef)},
                "train_r2": fit.train_r2,
                **summarise("marginal_rho", marg),
                **summarise("adjusted_rho", adj),
            }
            log(f"  {fam:20s} alpha={fit.alpha:8.3f} "
                f"marginal={results[arm][fam]['marginal_rho_mean']:+.4f} "
                f"adjusted={results[arm][fam]['adjusted_rho_mean']:+.4f}")

        for fam in families:
            if fam == BASELINE:
                continue
            part = partial_rho(y_test, preds[fam], preds[BASELINE])
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
