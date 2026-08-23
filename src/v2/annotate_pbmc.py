#!/usr/bin/env python
"""Confidence-gated PBMC cell-type annotation (canonical markers, cell-level z-score argmax,
same GAP-gating convention as the brain annotation). Writes a Barcode,Cell.Type CSV matching
the format build_atac_graph_v2.py / pilot's annotate_and_match.py already consume."""
import time, gzip, numpy as np, anndata as ad, scipy.sparse as sp
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
RNA = f"{ROOT}/data/multiome/pbmc10k_rna.h5ad"
GAP = 0.6

MARKERS = {
 "CD4T": ["IL7R", "CD3D", "CD3E", "CD4"], "CD8T": ["CD8A", "CD8B", "CD3D", "CD3E"],
 "NK": ["GNLY", "NKG7", "KLRD1", "NCAM1"], "B": ["MS4A1", "CD79A", "CD79B", "CD19"],
 "MonoCD14": ["CD14", "LYZ", "S100A8", "S100A9"], "MonoCD16": ["FCGR3A", "MS4A7"],
 "DC": ["FCER1A", "CST3", "CLEC9A"], "Platelet": ["PPBP", "PF4"],
}

A = ad.read_h5ad(RNA)
X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
rs = np.asarray(X.sum(1)).ravel(); rs[rs == 0] = 1
Xn = X.multiply(1e4 / rs[:, None]).tocsr(); Xn.data = np.log1p(Xn.data)
rsym = {str(s): i for i, s in enumerate(A.var_names)}
tnames = list(MARKERS)
score = np.zeros((A.shape[0], len(tnames)))
for ti, t in enumerate(tnames):
    gi = [rsym[g] for g in MARKERS[t] if g in rsym]
    log(f"  {t}: markers found {len(gi)}/{len(MARKERS[t])}")
    if not gi: continue
    score[:, ti] = np.asarray(Xn[:, gi].mean(1)).ravel()
zs = (score - score.mean(0)) / (score.std(0) + 1e-8)
order = np.argsort(-zs, axis=1)
top, second = order[:, 0], order[:, 1]
gap = zs[np.arange(len(zs)), top] - zs[np.arange(len(zs)), second]
assign = np.where((zs[np.arange(len(zs)), top] > 0) & (gap >= GAP), np.array(tnames)[top], "NA")
from collections import Counter
vc = Counter(assign.tolist())
log(f"assigned (GAP={GAP}): {dict(vc)}")

with gzip.open(f"{ROOT}/data/multiome/pbmc_cell_meta.csv.gz", "wt") as f:
    f.write("Barcode,Cell.Type\n")
    for bc, t in zip(A.obs_names, assign): f.write(f"{bc},{t}\n")
log(f"SAVED {ROOT}/data/multiome/pbmc_cell_meta.csv.gz")
