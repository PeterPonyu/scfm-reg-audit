#!/usr/bin/env python
"""Cross-tissue additive-marginal decomposition for the regulatory-potential proxy.

For each cross-tissue pair, decomposes the observed fixed-panel Spearman into:
- the part predictable from a TF-row + gene-column additive (marginal) structure
  fit on one tissue's proxy and applied to the other, and
- the residual agreement once each tissue's own additive structure is removed.

This quantifies how much of the "reproducibility" is sequence-fixed marginal
structure (JASPAR library + peak-gene proximity, tissue-independent by design)
versus tissue-specific edge structure. No randomization; descriptive only.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
import fixed_panel_audit as fpa  # noqa: E402

OUT = Path(fpa.OUT)
TAGS = ["GSE174367", "PBMC10k", "GSE206767"]
PAIRS = [("GSE174367", "PBMC10k"), ("GSE174367", "GSE206767"), ("PBMC10k", "GSE206767")]
N_ITER = 50


def load_consensus(tag):
    Z = np.load(OUT / f"G_ATAC_v2_{tag}.npz", allow_pickle=False)
    types = [str(t) for t in Z["types"]]
    G = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float64)
    return G, np.array(Z["tf_rows"])


def additive_fit(M, tf_rows, n_iter=N_ITER):
    """Fit M[i,j] ~= mu + r_i + c_j over the TF-row edge set by alternating centering."""
    sub = M[tf_rows, :].copy()
    mask = np.ones_like(sub, dtype=bool)
    for k, tf in enumerate(tf_rows):
        mask[k, tf] = False  # non-self edges only
    mu = sub[mask].mean()
    r = np.zeros(len(tf_rows))
    c = np.zeros(M.shape[1])
    for _ in range(n_iter):
        r_new = np.where(mask.sum(1) > 0, (sub - mu - c[None, :]).sum(1) / mask.sum(1), 0.0)
        c_new = np.where(mask.sum(0) > 0, (sub - mu - r_new[:, None]).sum(0) / mask.sum(0), 0.0)
        if np.allclose(r_new, r) and np.allclose(c_new, c):
            r, c = r_new, c_new
            break
        r, c = r_new, c_new
    return mu, r, c


def pair_stats(a, b, Ga, tf_a, Gb, tf_b):
    tf_common = np.intersect1d(tf_a, tf_b)
    Ng = Ga.shape[0]
    ii = np.repeat(tf_common, Ng)
    jj = np.tile(np.arange(Ng), len(tf_common))
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    x, y = Ga[ii, jj], Gb[ii, jj]

    observed = float(spearmanr(x, y).statistic)

    mu, r, c = additive_fit(Ga, tf_common)
    row_of = {tf: k for k, tf in enumerate(tf_common)}
    pred = mu + r[[row_of[t] for t in ii]] + c[jj]
    additive_pred_rho = float(spearmanr(pred, y).statistic)

    mu_b, r_b, c_b = additive_fit(Gb, tf_common)
    pred_b = mu_b + r_b[[row_of[t] for t in ii]] + c_b[jj]
    resid_rho = float(spearmanr(x - pred, y - pred_b).statistic)

    xb, yb = (x > 0).astype(float), (y > 0).astype(float)
    binary_phi = float(np.corrcoef(xb, yb)[0, 1])

    return {
        "pair": [a, b],
        "n_tf_common": int(len(tf_common)),
        "observed_spearman": round(observed, 6),
        "additive_pred_spearman": round(additive_pred_rho, 6),
        "residual_spearman_after_own_additive_fits": round(resid_rho, 6),
        "binary_support_phi": round(binary_phi, 6),
        "fraction_explained_by_additive_marginals": round(additive_pred_rho / observed, 4)
        if observed else None,
    }


def main():
    consensus = {}
    tfs = {}
    for tag in TAGS:
        consensus[tag], tfs[tag] = load_consensus(tag)
    rows = [pair_stats(a, b, consensus[a], tfs[a], consensus[b], tfs[b]) for a, b in PAIRS]
    doc = {
        "analysis": "cross_tissue_additive_decomp",
        "design": ("TF-row + gene-column additive fit by alternating centering over the fixed "
                   "non-self TF-target edge set; descriptive only, no randomization."),
        "n_iter_max": N_ITER,
        "rows": rows,
    }
    out = OUT / "cross_tissue_additive_decomp_v2.json"
    out.write_text(json.dumps(doc, indent=1) + "\n")
    for r in rows:
        print(r)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
