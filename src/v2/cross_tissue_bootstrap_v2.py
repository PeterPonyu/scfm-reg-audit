#!/usr/bin/env python
"""
scfm-reg-audit v2 — bootstrap CI on the cross-tissue TRUTH reproducibility (Fig 1B), for
consistency with the rest of the pipeline now that every other headline number has one. The
Mantel z-scores there were already enormous (z=40-51) so this isn't expected to overturn anything,
but the audit's own standard ("no bare point estimates") should apply uniformly, and it's cheap
given the graphs are cached. Block-bootstrap over TF nodes, same design as stats_enhanced_v2.py.
"""
import os, json, hashlib, time, numpy as np
from scipy.stats import spearmanr
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
NBOOT = int(os.environ.get("NBOOT", "2000"))

man = json.load(open(f"{ROOT}/data/manifest/shared_genes.v2.json")); genes = man["genes"]; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]


def consensus(tag):
    Z = np.load(f"{OUT}/G_ATAC_v2_{tag}.npz", allow_pickle=True)
    ts = [str(t) for t in Z["types"]]
    return np.mean([Z[f"G_{t}"] for t in ts], axis=0).astype(np.float32), np.array(Z["tf_rows"])


def bootstrap_pair(tagA, tagB, seed):
    GA, tfA = consensus(tagA); GB, tfB = consensus(tagB)
    tf_common = np.intersect1d(tfA, tfB)
    ii = np.repeat(tf_common, Ng); jj = np.tile(np.arange(Ng), len(tf_common)); m = ii != jj; ii, jj = ii[m], jj[m]
    a1, a2 = GA[ii, jj], GB[ii, jj]
    observed = float(spearmanr(a1, a2).statistic)

    rng = np.random.default_rng(seed)
    tf_pos = {t: np.where(tf_common == t)[0] for t in tf_common}
    vals = []
    for _ in range(NBOOT):
        samp = rng.choice(tf_common, size=len(tf_common), replace=True)
        pos = np.concatenate([tf_pos[t] for t in samp])
        v = spearmanr(a1[pos], a2[pos]).statistic
        if np.isfinite(v): vals.append(v)
    vals = np.array(vals)
    ci_lo, ci_hi = np.percentile(vals, [2.5, 97.5])
    log(f"  {tagA} vs {tagB}: observed={observed:.4f} CI95=[{ci_lo:.4f},{ci_hi:.4f}] (n_tf_common={len(tf_common)})")
    return dict(pair=[tagA, tagB], observed=round(observed, 4), n_tf_common=int(len(tf_common)),
               bootstrap_ci95=[round(float(ci_lo), 4), round(float(ci_hi), 4)], bootstrap_se=round(float(vals.std()), 4))


if __name__ == "__main__":
    res = [
        bootstrap_pair("GSE174367", "PBMC10k", 1),
        bootstrap_pair("GSE174367", "GSE206767", 2),
        bootstrap_pair("PBMC10k", "GSE206767", 3),
    ]
    json.dump(res, open(f"{OUT}/cross_tissue_bootstrap_v2.json", "w"), indent=2)
    log(f"SAVED {OUT}/cross_tissue_bootstrap_v2.json")
