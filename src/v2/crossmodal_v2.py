#!/usr/bin/env python
"""
scfm-reg-audit v2 — re-run the decisive test against the CONSTRUCT-VALID G_ATAC.

Pooled AD-brain, frozen manifest (sha-checked), corrected FM preprocessing (fm_readout).
Evaluated on the REGULATORY PAIR SET  P = {(i,j): i is a TF, i!=j}  (directed TF->target edges).
Reports, directly comparable to the pilot's crossmodal_pooled.json:
  marginal Spearman(coexp, ATAC) and Spearman(FM, ATAC), partial Spearman(FM, ATAC | coexp),
  plus Mantel gene-label permutation nulls (N).
Key question: does the pilot headline (coexp rho=0.097, FM partial~0) survive when the truth is
sequence-grounded regulation instead of co-accessibility covariation?
"""
import os, json, hashlib, time, numpy as np
from scipy.stats import spearmanr, rankdata
import fm_readout as fr
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"
GATAC = f"{ROOT}/results/v2/G_ATAC_v2_GSE174367.npz"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
OUT = f"{ROOT}/results/v2"; NPERM = int(os.environ.get("NPERM", "1000"))

man = json.load(open(MANI)); genes = man["genes"]
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"], "manifest hash mismatch"
Ng = len(genes); log(f"manifest {Ng} genes sha={man['sha256'][:12]} OK")

# ---- G_ATAC v2 consensus (mean over cell types) ----
Z = np.load(GATAC, allow_pickle=False)
zgenes = [str(g) for g in Z["genes"]]; assert zgenes == genes, "G_ATAC gene order != manifest"
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
log(f"G_ATAC v2 consensus {G_atac.shape} over {len(types)} types, {len(tf_rows)} TF rows")

# ---- regulatory pair set P (TF i -> target j, i != j) ----
ii, jj = [], []
tfset = set(tf_rows.tolist())
for i in tf_rows:
    for j in range(Ng):
        if j != i: ii.append(i); jj.append(j)
ii = np.array(ii); jj = np.array(jj)
log(f"pair set P: {len(ii)} directed TF->target pairs")
a = G_atac[ii, jj]

# ---- RNA -> manifest gene matrices (corrected normalization) ----
A, Xc, Xl, _ = fr.load_norm(RNA)
rsym = {str(s): k for k, s in enumerate(A.var_names)}
ri = np.array([rsym[g] for g in genes])
Xc_g = Xc[:, ri].tocsr(); Xl_g = Xl[:, ri].tocsr()
log(f"RNA cells {Xc_g.shape[0]} -> manifest genes")

# ---- graphs (cached) ----
fmc = f"{OUT}/fmgraphs_pooled_v2.npz"
if os.path.exists(fmc):
    F = np.load(fmc); G_co, G_gf, G_sg = F["co"], F["gf"], F["sg"]; log("loaded cached FM/coexp graphs")
else:
    G_co = fr.gene_coexp(Xl_g.toarray())
    rd = fr.FMReadout(genes, batch=8)
    cells = np.arange(Xc_g.shape[0])
    log("Geneformer embed…"); Egf = rd.geneformer(Xc_g, cells); G_gf = fr.FMReadout.cos_graph(Egf)
    log("scGPT embed…");      Esg = rd.scgpt(Xl_g, cells);      G_sg = fr.FMReadout.cos_graph(Esg)
    np.savez(fmc, co=G_co, gf=G_gf, sg=G_sg); log("cached FM/coexp graphs")

co = G_co[ii, jj]; gf = G_gf[ii, jj]; sg = G_sg[ii, jj]

def partial_spear(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def resid(v, c): c1 = np.c_[np.ones_like(c), c]; b = np.linalg.lstsq(c1, v, rcond=None)[0]; return v - c1 @ b
    return float(np.corrcoef(resid(rx, rz), resid(ry, rz))[0, 1])

obs = dict(
    coexp_vs_atac=float(spearmanr(co, a).statistic),
    geneformer_vs_atac=float(spearmanr(gf, a).statistic),
    scgpt_vs_atac=float(spearmanr(sg, a).statistic),
    geneformer_partial_given_coexp=partial_spear(gf, a, co),
    scgpt_partial_given_coexp=partial_spear(sg, a, co),
    geneformer_vs_coexp=float(spearmanr(gf, co).statistic),
    scgpt_vs_coexp=float(spearmanr(sg, co).statistic),
)

# ---- Mantel gene-label permutation null (permute G_ATAC gene ordering) ----
rng = np.random.default_rng(20260713)
keys = ["coexp_vs_atac", "geneformer_vs_atac", "geneformer_partial_given_coexp"]
null = {k: [] for k in keys}
for _ in range(NPERM):
    perm = rng.permutation(Ng)
    ap = G_atac[perm[ii], perm[jj]]
    null["coexp_vs_atac"].append(spearmanr(co, ap).statistic)
    null["geneformer_vs_atac"].append(spearmanr(gf, ap).statistic)
    null["geneformer_partial_given_coexp"].append(partial_spear(gf, ap, co))
stats = {}
for k in keys:
    nd = np.array(null[k]); o = obs[k]
    stats[k] = dict(observed=round(o, 4), null_mean=round(float(nd.mean()), 4), null_sd=round(float(nd.std()), 4),
                    z=round(float((o - nd.mean()) / (nd.std() + 1e-9)), 2),
                    p_perm=round(float((np.sum(np.abs(nd) >= abs(o)) + 1) / (NPERM + 1)), 4))

res = dict(construct="G_ATAC_v2 motif->accessible-peak->TF->target", pair_set="TF->target directed",
           n_pairs=int(len(ii)), n_tf=int(len(tf_rows)), n_genes=Ng, manifest_sha=man["sha256"],
           observed={k: round(v, 4) for k, v in obs.items()}, mantel=stats, n_perm=NPERM)
json.dump(res, open(f"{OUT}/crossmodal_v2.json", "w"), indent=2)
log("=== CROSS-MODAL v2 (construct-valid, pooled) ===")
for k, v in obs.items(): log(f"  {k}: {round(v,4)}")
for k, v in stats.items(): log(f"  MANTEL {k}: obs={v['observed']} z={v['z']} p={v['p_perm']}")
log(f"SAVED {OUT}/crossmodal_v2.json")
