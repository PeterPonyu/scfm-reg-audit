#!/usr/bin/env python
"""Permutation nulls and BH correction for the pair-level TF-disjoint probe.

Null: within each held-out TF, permute the gene labels of the probe's prediction
and recompute the confound-adjusted per-TF Spearman. This preserves each TF's
target degree and the confound design while destroying the pairing the probe
claims to have learnt. The test statistic is the across-TF mean, so the null is
built on the same aggregate that is reported.

Also reports a paired FM-vs-co-expression contrast: the per-TF difference in
adjusted rho, tested by sign-flipping the paired differences.
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
OUT = ROOT / "results/v2/tf_probe_pair_stats_v2.json"

N_PERM = 999
SEED_ROOT = 20260730
ALPHAS = np.logspace(-3, 3, 25)
BASELINE = "co_expression"


def log(*values):
    print(f"[{time.strftime('%H:%M:%S')}]", *values, flush=True)


def bh(pvals):
    """Benjamini-Hochberg q-values, order preserved."""
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0, 1)
    return q


def make_residualiser(design):
    """Return f(vec) -> rank residuals against a fixed design matrix."""
    pinv = np.linalg.pinv(design)

    def f(vec):
        r = rankdata(vec)
        return r - design @ (pinv @ r)

    return f


def adjusted_per_tf(y_true, y_pred, resid):
    """Confound-adjusted Spearman per TF for (n_tf, n_gene) blocks."""
    out = np.full(y_true.shape[0], np.nan)
    for i in range(y_true.shape[0]):
        a, b = resid(y_true[i]), resid(y_pred[i])
        if np.std(a) > 0 and np.std(b) > 0:
            out[i] = spearmanr(a, b).statistic
    return out


def main():
    """Run the null for every family in the primary arm and write stats JSON."""
    ev = json.loads(EVAL.read_text())
    arm = ev["primary_arm"]
    families = sorted(ev["results"][arm])
    n_genes = ev["n_genes"]
    n_tf = ev["n_test_tfs"]

    prov = json.loads((PAIR_DIR / "provenance.json").read_text())
    cols = [prov["feature_names"].index(f) for f in ev["arms"][arm]]

    tg = np.load(PAIR_DIR / "pair_targets.npz", allow_pickle=False)
    y_train = tg["y_train"].astype(np.float64)
    y_test = tg["y_test"].astype(np.float64).reshape(n_tf, n_genes)
    cg = tg["conf_gene"].astype(np.float64)
    cz = (cg - cg.mean(0)) / np.where(cg.std(0) == 0, 1.0, cg.std(0))
    resid = make_residualiser(np.column_stack([np.ones(n_genes), cz]))
    log(f"arm={arm} features={ev['arms'][arm]} test TFs={n_tf} perms={N_PERM}")

    obs, per_tf = {}, {}
    for fam_index, fam in enumerate(families):
        d = np.load(PAIR_DIR / f"{fam}_pairs.npz", allow_pickle=False)
        xtr = d["X_train"].astype(np.float64)[:, cols]
        xte = d["X_test"].astype(np.float64)[:, cols]
        mu, sd = xtr.mean(0), xtr.std(0)
        sd = np.where(sd == 0, 1.0, sd)
        model = RidgeCV(alphas=ALPHAS, scoring="neg_mean_squared_error")
        model.fit((xtr - mu) / sd, y_train)
        pred = model.predict((xte - mu) / sd).reshape(n_tf, n_genes)
        per_tf[fam] = adjusted_per_tf(y_test, pred, resid)
        obs[fam] = float(np.nanmean(per_tf[fam]))

        fam_seed = SEED_ROOT * 1000 + fam_index
        rng = np.random.default_rng(fam_seed)
        null = np.empty(N_PERM)
        for b in range(N_PERM):
            shuffled = np.stack([rng.permutation(pred[i]) for i in range(n_tf)])
            null[b] = np.nanmean(adjusted_per_tf(y_test, shuffled, resid))
        # Two-sided: the claim is "differs from chance", in either direction.
        p = (1 + np.sum(np.abs(null) >= abs(obs[fam]))) / (1 + N_PERM)
        obs[fam] = {
            "adjusted_rho_mean": obs[fam],
            "null_mean": float(null.mean()),
            "null_std": float(null.std()),
            "null_abs_p95": float(np.percentile(np.abs(null), 95)),
            "mantel_p": float(p),
            "mantel_seed": int(fam_seed),
        }
        log(f"  {fam:20s} obs={obs[fam]['adjusted_rho_mean']:+.4f} "
            f"null={null.mean():+.4f}+-{null.std():.4f} p={p:.4f}")

    qs = bh([obs[f]["mantel_p"] for f in families])
    for fam, q in zip(families, qs):
        obs[fam]["mantel_q"] = float(q)
        obs[fam]["significant_q05"] = bool(q < 0.05)

    contrasts = {}
    for fam in families:
        if fam == BASELINE:
            continue
        diff = per_tf[fam] - per_tf[BASELINE]
        diff = diff[np.isfinite(diff)]
        rng = np.random.default_rng(SEED_ROOT + 1)
        signs = rng.choice([-1.0, 1.0], size=(N_PERM, diff.size))
        null = (signs * diff).mean(axis=1)
        p = (1 + np.sum(np.abs(null) >= abs(diff.mean()))) / (1 + N_PERM)
        contrasts[fam] = {
            "paired_delta_mean": float(diff.mean()),
            "paired_delta_std": float(diff.std()),
            "n_tf": int(diff.size),
            "signflip_p": float(p),
        }
    cq = bh([contrasts[f]["signflip_p"] for f in contrasts])
    for fam, q in zip(contrasts, cq):
        contrasts[fam]["signflip_q"] = float(q)
        contrasts[fam]["significant_q05"] = bool(q < 0.05)
        log(f"  {fam:20s} vs {BASELINE}: delta={contrasts[fam]['paired_delta_mean']:+.4f} "
            f"q={contrasts[fam]['signflip_q']:.4f}")

    stats = {
        "schema_version": 2,
        "design": prov["design"],
        "arm": arm,
        "arm_features": ev["arms"][arm],
        "readout": "confound-adjusted per-TF Spearman, averaged over held-out TFs",
        "null": "within-TF gene-label permutation of probe predictions",
        "n_perm": N_PERM,
        "seed_root": SEED_ROOT,
        "seed_contract": "per_family_seed = seed_root*1000 + sorted_family_index; signflip_seed = seed_root+1",
        "n_test_tfs": n_tf,
        "n_train_tfs": ev["n_train_tfs"],
        "families": obs,
        "contrasts_vs_baseline": contrasts,
    }
    OUT.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    log(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
