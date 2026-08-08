#!/usr/bin/env python
"""
scfm-reg-audit v2 — cross-tissue reproducibility of the G_ATAC v2 regulatory truth (CPU-only).
Does the motif->accessible-peak->TF->target graph replicate across two independent tissues
(AD brain GSE174367 vs GSE206767)? High agreement => the truth (and thus the FM-null verdict)
generalizes; low => tissue-specific / noisy, a caveat. Same frozen manifest, same pair set P.
"""
import os, sys, json, hashlib, numpy as np
from scipy.stats import spearmanr, rankdata
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
man = json.load(open(f"{ROOT}/data/manifest/shared_genes.v2.json")); genes = man["genes"]; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]

def consensus(tag):
    Z = np.load(f"{OUT}/G_ATAC_v2_{tag}.npz", allow_pickle=True)
    ts = [str(t) for t in Z["types"]]
    return np.mean([Z[f"G_{t}"] for t in ts], axis=0).astype(np.float32), np.array(Z["tf_rows"])

TAG1 = sys.argv[1] if len(sys.argv) > 1 else "GSE174367"
TAG2 = sys.argv[2] if len(sys.argv) > 2 else "GSE206767"
OUTNAME = sys.argv[3] if len(sys.argv) > 3 else "cross_tissue_atac_v2.json"
G1, tf = consensus(TAG1); G2, _ = consensus(TAG2)
ii = np.repeat(tf, Ng); jj = np.tile(np.arange(Ng), len(tf)); m = ii != jj; ii, jj = ii[m], jj[m]
a1, a2 = G1[ii, jj], G2[ii, jj]

both_nz = ((a1 > 0) & (a2 > 0)).mean(); either_nz = ((a1 > 0) | (a2 > 0)).mean()
rho = float(spearmanr(a1, a2).statistic)
# Mantel null
rng = np.random.default_rng(20260713); null = []
for _ in range(1000):
    p = rng.permutation(Ng); null.append(spearmanr(a1, G2[p[ii], p[jj]]).statistic)
nd = np.array(null); z = float((rho - nd.mean()) / (nd.std() + 1e-9)); pperm = float((np.sum(np.abs(nd) >= abs(rho)) + 1) / 1001)
# per-TF row agreement
per_tf = []
for t in tf:
    r1, r2 = G1[t], G2[t]
    if (r1 > 0).sum() > 5 and (r2 > 0).sum() > 5: per_tf.append(spearmanr(r1, r2).statistic)
res = dict(tissues=[TAG1, TAG2], n_pairs=int(len(ii)),
           truth_spearman=round(rho, 4), mantel_z=round(z, 2), p_perm=round(pperm, 4),
           edge_jaccard=round(float(both_nz / (either_nz + 1e-9)), 3),
           mean_per_tf_row_spearman=round(float(np.mean(per_tf)), 4), n_tf_rows=len(per_tf))
json.dump(res, open(f"{OUT}/{OUTNAME}", "w"), indent=2)
for k, v in res.items(): print(f"  {k}: {v}")
print(f"SAVED {OUT}/{OUTNAME}")
