#!/usr/bin/env python
"""Confound-control gate for the UCE cross-modal result (same rigor as scf_confound_check.py /
ko_confound_check.py — required before trusting any positive partial). CPU."""
import os, json, hashlib, numpy as np, anndata as ad, pyfaidx
from scipy.stats import rankdata, spearmanr
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"; COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
ATAC = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
PROM = 2000; W = 500
man = json.load(open(f"{ROOT}/data/manifest/shared_genes.v2.json")); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
det = man["detection"]
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False); types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
G_co = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")["co"]
G_uce = np.load(f"{OUT}/G_uce_pooled.npz")["G"]

gco = {}
for ln in open(COORDS):
    c, s, e, st, nm = ln.rstrip("\n").split("\t")
    if nm not in gidx or nm in gco: continue
    s, e = int(s), int(e); lo = s - PROM if st == "+" else s; hi = e if st == "+" else e + PROM; gco[nm] = (c, lo, hi, abs(e - s))
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
ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
a, co, uce = G_atac[ii, jj], G_co[ii, jj], G_uce[ii, jj]
def z(v): v = v.astype(float); return (v - v.mean()) / (v.std() + 1e-9)
conf = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj]), z(tf_outdeg[ii]), z(atac_indeg[jj])]
def resid(y, C): C1 = np.c_[np.ones(len(y)), C]; b = np.linalg.lstsq(C1, y, rcond=None)[0]; return y - C1 @ b
def pcorr(x, y, C): return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])
ra, rco, ruce = rankdata(a), rankdata(co), rankdata(uce)
out = dict(
    uce_marginal=round(float(spearmanr(uce, a).statistic), 4),
    uce_partial_coexp=round(pcorr(ruce, ra, rco.reshape(-1, 1)), 4),
    uce_partial_coexp_confounds=round(pcorr(ruce, ra, np.c_[rco, conf]), 4),
    coexp_partial_confounds=round(pcorr(rco, ra, conf), 4),
)
json.dump(out, open(f"{OUT}/uce_confound_check.json", "w"), indent=2)
for k, v in out.items(): print(k, v)
print("SAVED uce_confound_check.json")
