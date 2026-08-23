#!/usr/bin/env python
"""
scfm-reg-audit v2 — CONFOUND-CONTROL edge regression (review: 'required for any beyond-coexp claim').

Does the FM-vs-regulation verdict survive controlling co-expression AND edge confounds?
Confounds (per directed TF->target edge i->j):
  target peak count, target gene length, target RNA detection, target peak GC%, TF out-degree.
Partial Spearman via rank-residualization on the full covariate design. CPU-only.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, scipy.sparse as sp, pyfaidx
from scipy.stats import rankdata, spearmanr
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
ATAC = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
PROMOTER = 2000; WIDTH = 500

man = json.load(open(MANI)); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
det = man["detection"]
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
F = np.load(f"{OUT}/fmgraphs_pooled_v2.npz"); G_co, G_gf = F["co"], F["gf"]
log(f"loaded graphs, {len(tf_rows)} TFs")

# ---- per-gene covariates: length, peak count, peak GC ----
gco = {}
for ln in open(COORDS):
    c, s, e, st, nm = ln.rstrip("\n").split("\t")
    if nm not in gidx or nm in gco: continue
    s, e = int(s), int(e); lo = s - PROMOTER if st == "+" else s; hi = e if st == "+" else e + PROMOTER
    gco[nm] = (c, lo, hi, abs(e - s))
Av = ad.read_h5ad(ATAC, backed="r"); peaks = [str(p) for p in Av.var_names]
pchr = np.array([p.split(":")[0] for p in peaks]); pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks])
pmid = (pse[:, 0] + pse[:, 1]) // 2
by_chr = {}
for i, c in enumerate(pchr): by_chr.setdefault(c, []).append(i)
for c in by_chr: by_chr[c] = np.array(by_chr[c])
fa = pyfaidx.Fasta(HG38, sequence_always_upper=True)
peakcount = np.zeros(Ng); genelen = np.zeros(Ng); gc = np.zeros(Ng)
for g, i in gidx.items():
    if g not in gco: continue
    c, lo, hi, ln_ = gco[g]; genelen[i] = ln_
    pis = by_chr.get(c)
    if pis is None: continue
    sel = pis[(pmid[pis] >= lo) & (pmid[pis] <= hi)]; peakcount[i] = len(sel)
    if len(sel) and c in fa.keys():
        gcs = []
        for p in sel[:40]:
            m = int(pmid[p]); seq = str(fa[c][max(0, m - WIDTH // 2):m + WIDTH // 2])
            if seq: gcs.append((seq.count("G") + seq.count("C")) / len(seq))
        gc[i] = np.mean(gcs) if gcs else 0.5
detv = np.array([det.get(g, 0.0) for g in genes])
tf_outdeg = (G_atac > 0).sum(1).astype(float)
atac_indeg = (G_atac > 0).sum(0).astype(float)
log("covariates built")

# ---- pair set P and vectors ----
ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
a = G_atac[ii, jj]; co = G_co[ii, jj]; gf = G_gf[ii, jj]
def z(v):
    v = v.astype(float); s = v.std(); return (v - v.mean()) / (s + 1e-9)
# 6-covariate spec matches scf_confound_check.py / uce_confound_check.py / ko_confound_check.py
# (this script originally used a different spec -- log1p(genelen) instead of target in-degree --
# reconciled 2026-07-13 so all 4 FM confound checks are directly comparable)
conf = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj]), z(tf_outdeg[ii]), z(atac_indeg[jj])]

ra, rco, rgf = rankdata(a), rankdata(co), rankdata(gf)
def resid(y, C):
    C1 = np.c_[np.ones(len(y)), C]; b = np.linalg.lstsq(C1, y, rcond=None)[0]; return y - C1 @ b
def pcorr(x, y, C): return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])

res = dict(
    n_pairs=int(len(ii)), n_tf=int(len(tf_rows)),
    fm_vs_atac_marginal=round(float(spearmanr(gf, a).statistic), 4),
    fm_partial_coexp_only=round(pcorr(rgf, ra, rco.reshape(-1, 1)), 4),
    fm_partial_coexp_plus_confounds=round(pcorr(rgf, ra, np.c_[rco, conf]), 4),
    coexp_vs_atac_marginal=round(float(spearmanr(co, a).statistic), 4),
    coexp_partial_confounds_only=round(pcorr(rco, ra, conf), 4),
)
json.dump(res, open(f"{OUT}/confound_regression_v2.json", "w"), indent=2)
log("=== CONFOUND-CONTROL REGRESSION (pooled, TF->target) ===")
for k, v in res.items(): log(f"  {k}: {v}")
log(f"SAVED {OUT}/confound_regression_v2.json")
