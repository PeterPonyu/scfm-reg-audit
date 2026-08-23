#!/usr/bin/env python
"""
scfm-reg-audit v2 — scFoundation (3rd FM) vs construct-valid G_ATAC v2, pooled AD brain.

Same test suite as crossmodal_v2.py / confound_regression_v2.py (marginal, partial|coexp,
confound-controlled partial, Mantel null) applied to a 3rd, architecturally distinct FM
(scFoundation: autobin value encoding + read-depth tokens, not rank-based like Geneformer/
scGPT) via the encoder-only readout in fm_readout_scf.py. Reuses the already-built G_ATAC v2
and cached G_coexp (fmgraphs_pooled_v2.npz) — only the new FM graph is computed here.
CPU (GPU contended by a concurrent unrelated job); N_CELLS/CAP tuned for CPU wall-clock.
"""
import os, json, hashlib, time, numpy as np, torch
from scipy.stats import spearmanr, rankdata
import fm_readout as fr, fm_readout_scf as fscf
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
N_CELLS = int(os.environ.get("N_CELLS", "1500")); CAP = int(os.environ.get("SCF_CAP", "512"))
BATCH = int(os.environ.get("SCF_BATCH", "8")); NTHREAD = int(os.environ.get("NTHREAD", "20"))
torch.set_num_threads(NTHREAD)

man = json.load(open(MANI)); genes = man["genes"]; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
F = np.load(f"{OUT}/fmgraphs_pooled_v2.npz"); G_co = F["co"]           # reuse cached pooled co-expression
log(f"manifest {Ng} genes, G_ATAC + G_coexp loaded")

A, Xc, Xl, _ = fr.load_norm(RNA)
rsym = {str(s): k for k, s in enumerate(A.var_names)}
ri = np.array([rsym[g] for g in genes])
Xc_g = Xc[:, ri].tocsr()
rng = np.random.default_rng(20260713)
cells = rng.choice(Xc_g.shape[0], size=min(N_CELLS, Xc_g.shape[0]), replace=False)
log(f"scFoundation encoder embedding over {len(cells)} cells (CPU, cap={CAP}, batch={BATCH})")

rd = fscf.SCFReadout(genes, batch=BATCH, cap=CAP)
assert len(rd.manifest_ok) == Ng, f"manifest overlap {len(rd.manifest_ok)}/{Ng}"
t0 = time.time()
E = rd.gene_embed(Xc_g, cells)
covered = int((np.abs(E).sum(1) > 0).sum())
log(f"embedded in {time.time()-t0:.0f}s | manifest genes covered: {covered}/{Ng}")
G_scf = fr.FMReadout.cos_graph(E)
np.savez(f"{OUT}/G_scf_pooled.npz", G=G_scf, covered=covered)

ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
a = G_atac[ii, jj]; co = G_co[ii, jj]; scf = G_scf[ii, jj]

def partial_spear(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def resid(v, c): c1 = np.c_[np.ones_like(c), c]; b = np.linalg.lstsq(c1, v, rcond=None)[0]; return v - c1 @ b
    return float(np.corrcoef(resid(rx, rz), resid(ry, rz))[0, 1])

obs = dict(
    scf_vs_atac=float(spearmanr(scf, a).statistic),
    scf_vs_coexp=float(spearmanr(scf, co).statistic),
    scf_partial_given_coexp=partial_spear(scf, a, co),
    coexp_vs_atac=float(spearmanr(co, a).statistic),
)

rng2 = np.random.default_rng(11); NPERM = int(os.environ.get("NPERM", "1000")); null = []
for _ in range(NPERM):
    p = rng2.permutation(Ng); null.append(partial_spear(scf, G_atac[p[ii], p[jj]], co))
nd = np.array(null); o = obs["scf_partial_given_coexp"]
mantel = dict(observed=round(o, 4), null_mean=round(float(nd.mean()), 4), null_sd=round(float(nd.std()), 4),
              z=round(float((o - nd.mean()) / (nd.std() + 1e-9)), 2),
              p_perm=round(float((np.sum(np.abs(nd) >= abs(o)) + 1) / (NPERM + 1)), 4))

res = dict(fm="scFoundation_encoder_only", n_cells=int(len(cells)), cap=CAP, manifest_genes_covered=covered,
           n_pairs=int(len(ii)), observed={k: round(v, 4) for k, v in obs.items()}, mantel_partial=mantel)
json.dump(res, open(f"{OUT}/crossmodal_scf_v2.json", "w"), indent=2)
log("=== scFoundation CROSS-MODAL (pooled, AD brain) ===")
for k, v in obs.items(): log(f"  {k}: {round(v,4)}")
log(f"  MANTEL scf_partial: z={mantel['z']} p={mantel['p_perm']}")
log(f"SAVED {OUT}/crossmodal_scf_v2.json")
