#!/usr/bin/env python
"""Shared pieces of the TF-disjoint pair-level probe (fit + readout + nulls).

`run_pair_probe.py` (evaluation) and `pair_probe_stats.py` (permutation nulls)
both fit one RidgeCV per family on the same cached pair features and score
held-out TFs with rank-residualised Spearman correlations; everything they have
in common lives here.
"""
import re
from pathlib import Path
from typing import List, NamedTuple

import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import RidgeCV

import audit_utils as au

ROOT = Path(__file__).resolve().parents[2]
PAIR_DIR = ROOT / "results/v2/tf_probe_pair"
EVAL = ROOT / "results/v2/tf_probe_pair_eval_v2.json"
STATS = ROOT / "results/v2/tf_probe_pair_stats_v2.json"

ALPHAS = np.logspace(-3, 3, 25)
BASELINE = "co_expression"
FAMILY_RE = re.compile(r"\A[A-Za-z0-9_]+\Z")

log = au.log
bh = au.bh_qvalues


def family_pairs_path(family):
    """Resolve a family's pair file, rejecting names that escape PAIR_DIR."""
    if not FAMILY_RE.match(family):
        raise ValueError(f"invalid family name: {family!r}")
    return PAIR_DIR / f"{family}_pairs.npz"


def rank_residualiser(design):
    """Return f(vec) -> rank residuals against a fixed intercept-augmented design."""
    pinv = np.linalg.pinv(design)

    def residualise(vec):
        r = rankdata(vec)
        return r - design @ (pinv @ r)

    return residualise


def confound_design(conf_gene, n_genes):
    """Intercept plus standardised gene-level confounds, replicated per test TF."""
    conf = np.asarray(conf_gene, dtype=np.float64)
    sd = conf.std(0)
    z = (conf - conf.mean(0)) / np.where(sd == 0, 1.0, sd)
    return np.column_stack([np.ones(n_genes), z])


class ProbeFit(NamedTuple):
    """One family's fitted probe and its held-out predictions."""

    alpha: float
    coef: np.ndarray
    train_r2: float
    prediction: np.ndarray
    feature_names: List[str]


def fit_family_probe(family, cols, y_train) -> ProbeFit:
    """Fit one family's RidgeCV on standardised train features; predict test pairs."""
    d = np.load(family_pairs_path(family), allow_pickle=False)
    xtr = d["X_train"].astype(np.float64)[:, cols]
    xte = d["X_test"].astype(np.float64)[:, cols]
    mu, sd = xtr.mean(0), xtr.std(0)
    sd = np.where(sd == 0, 1.0, sd)
    xtr_z = (xtr - mu) / sd
    model = RidgeCV(alphas=ALPHAS, scoring="neg_mean_squared_error")
    model.fit(xtr_z, y_train)
    return ProbeFit(
        alpha=float(model.alpha_),
        coef=model.coef_,
        train_r2=float(model.score(xtr_z, y_train)),
        prediction=model.predict((xte - mu) / sd),
        feature_names=[str(v) for v in d["feature_names"]],
    )


def per_tf_rho(y_true, y_pred, residualise=None):
    """Spearman per test TF for (n_tf, n_gene) blocks, optionally after residualising.

    Degenerate TFs (zero variance on either side) yield NaN.
    """
    out = np.full(y_true.shape[0], np.nan)
    for i in range(y_true.shape[0]):
        a, b = y_true[i], y_pred[i]
        if residualise is not None:
            a, b = residualise(a), residualise(b)
        if np.std(a) > 0 and np.std(b) > 0:
            out[i] = spearmanr(a, b).statistic
    return out


def blocks(vector, n_tf, n_genes):
    """Reshape a flat pair vector into (n_tf, n_genes) per-TF blocks."""
    return np.asarray(vector, dtype=np.float64).reshape(n_tf, n_genes)


def permutation_p(null, observed, n_perm):
    """Two-sided plus-one Monte-Carlo p-value: the claim is 'differs from chance'."""
    return au.mc_pvalue(null, observed, n_perm)[0]
