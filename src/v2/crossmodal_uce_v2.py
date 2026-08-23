#!/usr/bin/env python
"""
scfm-reg-audit v2 — UCE (4th FM) vs construct-valid G_ATAC v2, pooled AD brain.
Same test suite as crossmodal_scf_v2.py (marginal, partial|coexp, Mantel null); confound-controlled
partial computed via a follow-up script (uce_confound_check.py) mirroring scf_confound_check.py.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, scipy.sparse as sp
from scipy.stats import spearmanr, rankdata
import fm_readout as fr, fm_readout_uce as fuce
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
N_CELLS = int(os.environ.get("N_CELLS", "2000")); BATCH = int(os.environ.get("UCE_BATCH", "6"))

man = json.load(open(MANI)); genes = man["genes"]; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
F = np.load(f"{OUT}/fmgraphs_pooled_v2.npz"); G_co = F["co"]
log(f"manifest {Ng} genes, G_ATAC + G_coexp loaded")

rd = fuce.UCEReadout(genes, batch=BATCH)

A = ad.read_h5ad(RNA)
X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)   # RAW counts (UCE's official weighting needs raw, not CP10k)
rsym = {str(s): i for i, s in enumerate(A.var_names)}
uni_cols = np.array([rsym.get(g, -1) for g in rd.uni_syms])
present = uni_cols >= 0
Xraw = X[:, np.where(present, uni_cols, 0)].tocsr()
Xraw = Xraw.multiply(present[None, :]).tocsr()
log(f"universe genes present in RNA: {int(present.sum())}/{len(rd.uni_syms)}")

rng = np.random.default_rng(20260713)
cells = rng.choice(Xraw.shape[0], size=min(N_CELLS, Xraw.shape[0]), replace=False)
log(f"UCE embedding over {len(cells)} cells (batch={BATCH})")
t0 = time.time()
E = rd.gene_embed(Xraw, cells, uni_cols)
covered = int((np.abs(E).sum(1) > 0).sum())
log(f"embedded in {time.time()-t0:.0f}s | manifest genes covered: {covered}/{Ng}")
G_uce = fr.FMReadout.cos_graph(E)
np.savez(f"{OUT}/G_uce_pooled.npz", G=G_uce, covered=covered)

ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
a = G_atac[ii, jj]; co = G_co[ii, jj]; uce = G_uce[ii, jj]

def partial_spear(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def resid(v, c): c1 = np.c_[np.ones_like(c), c]; b = np.linalg.lstsq(c1, v, rcond=None)[0]; return v - c1 @ b
    return float(np.corrcoef(resid(rx, rz), resid(ry, rz))[0, 1])

obs = dict(
    uce_vs_atac=float(spearmanr(uce, a).statistic),
    uce_vs_coexp=float(spearmanr(uce, co).statistic),
    uce_partial_given_coexp=partial_spear(uce, a, co),
    coexp_vs_atac=float(spearmanr(co, a).statistic),
)
rng2 = np.random.default_rng(13); NPERM = int(os.environ.get("NPERM", "1000")); null = []
for _ in range(NPERM):
    p = rng2.permutation(Ng); null.append(partial_spear(uce, G_atac[p[ii], p[jj]], co))
nd = np.array(null); o = obs["uce_partial_given_coexp"]
mantel = dict(observed=round(o, 4), null_mean=round(float(nd.mean()), 4), null_sd=round(float(nd.std()), 4),
              z=round(float((o - nd.mean()) / (nd.std() + 1e-9)), 2),
              p_perm=round(float((np.sum(np.abs(nd) >= abs(o)) + 1) / (NPERM + 1)), 4))

res = dict(fm="UCE_4layer", n_cells=int(len(cells)), manifest_genes_covered=covered,
           n_pairs=int(len(ii)), observed={k: round(v, 4) for k, v in obs.items()}, mantel_partial=mantel)
json.dump(res, open(f"{OUT}/crossmodal_uce_v2.json", "w"), indent=2)
log("=== UCE CROSS-MODAL (pooled, AD brain) ===")
for k, v in obs.items(): log(f"  {k}: {round(v,4)}")
log(f"  MANTEL uce_partial: z={mantel['z']} p={mantel['p_perm']}")
log(f"SAVED {OUT}/crossmodal_uce_v2.json")
