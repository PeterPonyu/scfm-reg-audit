#!/usr/bin/env python
"""Sanity-check G_ATAC v2: is it a differentiated regulatory graph, or rank-1 collapsed?
Guards against 'every TF hits every peak' -> every TF row proportional to target accessibility."""
import os, json, numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
Z = np.load(f"{ROOT}/results/v2/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
genes = [str(g) for g in Z["genes"]]; types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G = np.mean([Z[f"G_{t}"] for t in types], axis=0)
Ng = len(genes)

sub = G[tf_rows]                                   # (nTF x Ng) TF rows only
nz_per_tf = (sub > 0).sum(1)
print(f"genes {Ng} | TF rows {len(tf_rows)} | types {len(types)}")
print(f"edge density on TF rows: {(sub>0).mean():.4f}")
print(f"targets/TF: mean {nz_per_tf.mean():.0f}  median {np.median(nz_per_tf):.0f}  "
      f"min {nz_per_tf.min()}  max {nz_per_tf.max()}")

# rank-1 collapse test: SVD energy in top component
Sc = sub / (np.linalg.norm(sub, axis=1, keepdims=True) + 1e-9)   # row-normalized
s = np.linalg.svd(Sc, compute_uv=False)
print(f"SVD top-1 energy frac: {(s[0]**2/ (s**2).sum()):.3f}  (near 1 => rank-1 collapse/BAD)")
print(f"SVD top-5 energy frac: {((s[:5]**2).sum()/(s**2).sum()):.3f}")

# TF-profile distinctness: mean pairwise cosine between TF target profiles (lower = more distinct)
import numpy.random as npr
idx = npr.default_rng(0).choice(len(tf_rows), size=min(200, len(tf_rows)), replace=False)
C = Sc[idx] @ Sc[idx].T; iu = np.triu_indices(len(idx), 1)
print(f"mean pairwise TF-profile cosine (200 sampled): {C[iu].mean():.3f}  (near 1 => TFs indistinct)")

# spot-check a few TFs' top targets
tf_names = {int(r): genes[int(r)] for r in tf_rows}
for want in ["SPI1", "NRF1", "CTCF", "NEUROD2", "SOX10", "MEF2C"]:
    if want in genes and genes.index(want) in set(tf_rows.tolist()):
        r = genes.index(want); row = G[r]; top = np.argsort(-row)[:6]
        print(f"  {want} top targets:", [(genes[j], round(float(row[j]), 2)) for j in top if row[j] > 0])
