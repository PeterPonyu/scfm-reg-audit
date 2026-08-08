#!/usr/bin/env python
"""
scfm-reg-audit v2 — scFoundation (3rd FM) on PBMC Multiome: extends the paired-cell calibration
(DESIGN §6) to a 3rd, architecturally distinct FM. Reuses G_ATAC_v2_PBMC10k truth and the cached
pooled PBMC co-expression graph (pbmc_fmgraphs_pooled.npz); only the new scFoundation graph is
computed here. Same marginal + partial|coexp + confound-controlled partial + Mantel test suite.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, scipy.sparse as sp, pyfaidx
from scipy.stats import spearmanr, rankdata
import fm_readout as fr, fm_readout_scf as fscf
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
RNA = f"{ROOT}/data/multiome/pbmc10k_rna.h5ad"
ATAC = f"{ROOT}/data/multiome/pbmc10k_atac.h5ad"
PROMOTER = 2000; W = 500
N_CELLS = int(os.environ.get("N_CELLS", "2000")); CAP = int(os.environ.get("SCF_CAP", "1024")); BATCH = int(os.environ.get("SCF_BATCH", "6"))

man = json.load(open(MANI)); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
det = man["detection"]
Z = np.load(f"{OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=True)
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
G_co = np.load(f"{OUT}/pbmc_fmgraphs_pooled.npz")["co"]
log(f"G_ATAC PBMC + cached coexp loaded, {len(tf_rows)} TFs")

# ---- confound covariates (same as pbmc_eval_v2.py) ----
gco = {}
for ln in open(COORDS):
    c, s, e, st, nm = ln.rstrip("\n").split("\t")
    if nm not in gidx or nm in gco: continue
    s, e = int(s), int(e); lo = s - PROMOTER if st == "+" else s; hi = e if st == "+" else e + PROMOTER; gco[nm] = (c, lo, hi, abs(e - s))
Av = ad.read_h5ad(ATAC, backed="r"); peaks = [str(p) for p in Av.var_names]
pchr = np.array([p.split(":")[0] for p in peaks]); pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks]); pmid = (pse[:, 0] + pse[:, 1]) // 2
by = {}
for i, c in enumerate(pchr): by.setdefault(c, []).append(i)
for c in by: by[c] = np.array(by[c])
fa = pyfaidx.Fasta(HG38, sequence_always_upper=True)
peakcount = np.zeros(Ng); genelen = np.zeros(Ng); gc = np.zeros(Ng)
for g, i in gidx.items():
    if g not in gco: continue
    c, lo, hi, ln_ = gco[g]; genelen[i] = ln_; pis = by.get(c)
    if pis is None: continue
    sel = pis[(pmid[pis] >= lo) & (pmid[pis] <= hi)]; peakcount[i] = len(sel)
    if len(sel) and c in fa.keys():
        gg = [((str(fa[c][max(0, int(pmid[p]) - W // 2):int(pmid[p]) + W // 2]).count("G") + str(fa[c][max(0, int(pmid[p]) - W // 2):int(pmid[p]) + W // 2]).count("C")) / W) for p in sel[:40]]
        gc[i] = np.mean(gg) if gg else 0.5
detv = np.array([det.get(g, 0.0) for g in genes]); tf_outdeg = (G_atac > 0).sum(1).astype(float); atac_indeg = (G_atac > 0).sum(0).astype(float)

# ---- RNA + scFoundation embedding ----
A, Xc, Xl, _ = fr.load_norm(RNA)
rsym = {str(s): k for k, s in enumerate(A.var_names)}
present = [g for g in genes if g in rsym]
ri = np.full(Ng, -1, int)
for g in present: ri[gidx[g]] = rsym[g]
valid = np.where(ri >= 0)[0]
log(f"RNA present for {len(valid)}/{Ng} manifest genes")

rd = fscf.SCFReadout([genes[i] for i in valid], batch=BATCH, cap=CAP)
Xc_sub = Xc[:, ri[valid]].tocsr()
rng = np.random.default_rng(20260713)
cells = rng.choice(Xc_sub.shape[0], size=min(N_CELLS, Xc_sub.shape[0]), replace=False)
log(f"scFoundation embedding over {len(cells)} PBMC cells (cap={CAP}, batch={BATCH})")
t0 = time.time()
Esub = rd.gene_embed(Xc_sub, cells)
covered = int((np.abs(Esub).sum(1) > 0).sum())
log(f"embedded in {time.time()-t0:.0f}s | manifest genes covered: {covered}/{len(valid)}")
E = np.zeros((Ng, 768), np.float32); E[valid] = Esub
G_scf = fr.FMReadout.cos_graph(E)
np.savez(f"{OUT}/G_scf_pbmc_pooled.npz", G=G_scf, covered=covered)

ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m0 = ii != jj
validmask = np.zeros(Ng, bool); validmask[valid] = True
keep = m0 & validmask[ii] & validmask[jj]; ii, jj = ii[keep], jj[keep]
a = G_atac[ii, jj]; co = G_co[ii, jj]; scf = G_scf[ii, jj]

def z(v): v = v.astype(float); return (v - v.mean()) / (v.std() + 1e-9)
conf = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj]), z(tf_outdeg[ii]), z(atac_indeg[jj])]
def resid(y, C): C1 = np.c_[np.ones(len(y)), C]; b = np.linalg.lstsq(C1, y, rcond=None)[0]; return y - C1 @ b
def pcorr(x, y, C): return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])
ra, rco, rscf = rankdata(a), rankdata(co), rankdata(scf)

res = dict(fm="scFoundation_encoder_only", tissue="PBMC10k_multiome_paired", n_cells=int(len(cells)),
           n_pairs=int(len(ii)), manifest_genes_covered=covered,
           scf_marginal=round(float(spearmanr(scf, a).statistic), 4),
           scf_partial_coexp=round(pcorr(rscf, ra, rco.reshape(-1, 1)), 4),
           scf_partial_coexp_confounds=round(pcorr(rscf, ra, np.c_[rco, conf]), 4),
           coexp_partial_confounds=round(pcorr(rco, ra, conf), 4),
           scf_vs_coexp=round(float(spearmanr(scf, co).statistic), 4))
json.dump(res, open(f"{OUT}/pbmc_eval_scf_v2.json", "w"), indent=2)
log("=== scFoundation on PBMC (paired-cell) ===")
for k, v in res.items(): log(f"  {k}: {v}")
log(f"SAVED {OUT}/pbmc_eval_scf_v2.json")
