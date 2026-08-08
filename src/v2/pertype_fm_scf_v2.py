#!/usr/bin/env python
"""
scfm-reg-audit v2 — PER-CELL-TYPE scFoundation readout (GPU) vs G_ATAC v2, gated annotation.
Extends pertype_fm_v2.py (Geneformer embedding+attention) to the 3rd FM. Same gated annotation,
marker exclusion, marginal+partial|coexp test, per brain type. GPU now free (~24GB) so full CAP.
"""
import os, json, hashlib, time, numpy as np
from scipy.stats import spearmanr, rankdata
import fm_readout as fr, fm_readout_scf as fscf
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
GAP = float(os.environ.get("GAP", "0.75")); MINCELL = 150
CAP = int(os.environ.get("SCF_CAP", "1024")); BATCH = int(os.environ.get("SCF_BATCH", "6"))
MARKERS = {
 "ODC": ["MOBP", "MOG", "PLP1", "MBP", "ST18", "CTNNA3"], "EX": ["SLC17A7", "SATB2", "RBFOX3", "NRGN", "SLC17A6"],
 "INH": ["GAD1", "GAD2", "SLC32A1", "DLX1"], "ASC": ["AQP4", "GFAP", "SLC1A2", "ALDH1L1"],
 "MG": ["CSF1R", "P2RY12", "CX3CR1", "C1QA", "C1QB"], "OPC": ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "VCAN"],
 "PER.END": ["CLDN5", "FLT1", "PECAM1", "PDGFRB", "RGS5"],
}
ALL_MARKERS = {g for v in MARKERS.values() for g in v}

man = json.load(open(MANI)); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
atac_types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"]); Gatac = {t: Z[f"G_{t}"] for t in atac_types}

A, Xc, Xl, _ = fr.load_norm(RNA, ctcol="cell_type")
rsym = {str(s): k for k, s in enumerate(A.var_names)}; ri = np.array([rsym[g] for g in genes])
Xc_g = Xc[:, ri].tocsr(); Xl_g = Xl[:, ri].tocsr()
clusters = A.obs["cell_type"].astype(str).values; cl_ids = sorted(set(clusters)); tnames = list(MARKERS)
score = np.zeros((len(cl_ids), len(tnames)))
for ti, t in enumerate(tnames):
    gi = [rsym[g] for g in MARKERS[t] if g in rsym]
    ms = np.asarray(Xl[:, gi].mean(1)).ravel()
    for ci, c in enumerate(cl_ids): score[ci, ti] = ms[clusters == c].mean()
zs = (score - score.mean(0)) / (score.std(0) + 1e-8)
cl2type = {}
for ci, c in enumerate(cl_ids):
    o = np.argsort(-zs[ci])
    cl2type[c] = tnames[o[0]] if (zs[ci, o[0]] > 0 and zs[ci, o[0]] - zs[ci, o[1]] >= GAP) else None
rna_type = np.array([cl2type[c] for c in clusters])
log(f"gated: {sum(1 for v in cl2type.values() if v)}/{len(cl_ids)} clusters assigned")

mask_gene = np.array([g not in ALL_MARKERS for g in genes])
ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows))
keep = (ii != jj) & mask_gene[ii] & mask_gene[jj]; ii, jj = ii[keep], jj[keep]

def partial(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def rsd(v, c): c1 = np.c_[np.ones_like(c), c]; b = np.linalg.lstsq(c1, v, rcond=None)[0]; return v - c1 @ b
    return float(np.corrcoef(rsd(rx, rz), rsd(ry, rz))[0, 1])

rd = fscf.SCFReadout(genes, batch=BATCH, cap=CAP); rows = []
for t in atac_types:
    cids = np.where(rna_type == t)[0]
    if len(cids) < MINCELL: log(f"  skip {t}: {len(cids)} cells"); continue
    a = Gatac[t][ii, jj]
    Gco = fr.gene_coexp(Xl_g[cids].toarray()); co = Gco[ii, jj]
    t0 = time.time()
    Escf = rd.gene_embed(Xc_g, cids)
    covered = int((np.abs(Escf).sum(1) > 0).sum())
    Gscf = fr.FMReadout.cos_graph(Escf); scf = Gscf[ii, jj]
    np.savez(f"{OUT}/brain_scfgraphs_{t}.npz", scf=Gscf)  # cache for downstream confound+bootstrap analysis
    rows.append(dict(cell_type=t, n=int(len(cids)), covered=covered, sec=round(time.time() - t0, 0),
                     coexp_vs_atac=round(float(spearmanr(co, a).statistic), 4),
                     scf_vs_atac=round(float(spearmanr(scf, a).statistic), 4),
                     scf_partial=round(partial(scf, a, co), 4),
                     scf_vs_coexp=round(float(spearmanr(scf, co).statistic), 4)))
    log(f"  {t}: n={len(cids)} covered={covered}/{Ng} coexp={rows[-1]['coexp_vs_atac']} "
        f"scf_partial={rows[-1]['scf_partial']} scf_vs_coexp={rows[-1]['scf_vs_coexp']} ({rows[-1]['sec']:.0f}s)")

def mean(k): return round(float(np.mean([r[k] for r in rows])), 4) if rows else None
summary = dict(gap=GAP, cap=CAP, per_type=rows, markers_excluded=True,
               mean=dict((k, mean(k)) for k in ["coexp_vs_atac", "scf_partial", "scf_vs_coexp"]))
json.dump(summary, open(f"{OUT}/pertype_fm_scf_v2.json", "w"), indent=2)
log(f"=== PER-TYPE scFoundation means: {summary['mean']} ===")
log(f"SAVED {OUT}/pertype_fm_scf_v2.json")
