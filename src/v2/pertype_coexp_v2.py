#!/usr/bin/env python
"""
scfm-reg-audit v2 — PER-CELL-TYPE coexp vs regulation, with a CONFIDENCE-GATED annotation
(fixes the pilot's forced marker->argmax HIGH). CPU-only; FM-per-type readout is GPU-deferred.

- Annotate ad_hm clusters -> brain type by marker z-score, assign ONLY if (z_top - z_2nd) >= GAP.
  Ungated clusters are dropped (not forced).
- Marker genes are EXCLUDED from the evaluation pair set (avoids markers-in-eval circularity).
- For each ATAC type with a confidently-matched RNA population: per-type |Pearson| coexp vs
  G_ATAC_v2[type] on the regulatory pair set P (TF->target). Reports per-type + mean.
"""
import os, json, hashlib, time, numpy as np
from scipy.stats import spearmanr
import fm_readout as fr
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
GAP = float(os.environ.get("GAP", "0.75"))

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
atac_types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
Gatac = {t: Z[f"G_{t}"] for t in atac_types}

A, Xc, Xl, _ = fr.load_norm(RNA, ctcol="cell_type")
rsym = {str(s): k for k, s in enumerate(A.var_names)}; ri = np.array([rsym[g] for g in genes])
Xl_g = Xl[:, ri].tocsr()
clusters = A.obs["cell_type"].astype(str).values
cl_ids = sorted(set(clusters)); tnames = list(MARKERS)

# marker z-score per (cluster, type)
score = np.zeros((len(cl_ids), len(tnames)))
for ti, t in enumerate(tnames):
    gi = [rsym[g] for g in MARKERS[t] if g in rsym]
    if not gi: continue
    ms = np.asarray(Xl[:, gi].mean(1)).ravel()
    for ci, c in enumerate(cl_ids): score[ci, ti] = ms[clusters == c].mean()
zs = (score - score.mean(0)) / (score.std(0) + 1e-8)
cl2type = {}
for ci, c in enumerate(cl_ids):
    order = np.argsort(-zs[ci]); top, second = order[0], order[1]
    if zs[ci, top] > 0 and (zs[ci, top] - zs[ci, second]) >= GAP:
        cl2type[c] = tnames[top]
    else:
        cl2type[c] = None                                  # ungated -> dropped
rna_type = np.array([cl2type[c] for c in clusters])
assigned = {c: t for c, t in cl2type.items() if t}
log(f"gated annotation (GAP={GAP}): {len(assigned)}/{len(cl_ids)} clusters assigned: {assigned}")

# evaluation pair set (exclude marker genes)
mask_gene = np.array([g not in ALL_MARKERS for g in genes])
ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); keep = (ii != jj) & mask_gene[ii] & mask_gene[jj]
ii, jj = ii[keep], jj[keep]
log(f"pair set P (markers excluded): {len(ii)} pairs")

rows = []
for t in atac_types:
    cids = np.where(rna_type == t)[0]
    if len(cids) < 150: log(f"  skip {t}: {len(cids)} confidently-typed RNA cells"); continue
    Gc = fr.gene_coexp(Xl_g[cids].toarray())
    a = Gatac[t][ii, jj]; co = Gc[ii, jj]
    rows.append(dict(cell_type=t, n_rna=int(len(cids)),
                     coexp_vs_atac=round(float(spearmanr(co, a).statistic), 4)))
    log(f"  {t}: n={len(cids)} coexp_vs_atac={rows[-1]['coexp_vs_atac']}")

summary = dict(gap=GAP, annotation={c: t for c, t in cl2type.items()},
               n_assigned_clusters=len(assigned), markers_excluded_from_eval=True,
               per_type=rows,
               mean_coexp_vs_atac=round(float(np.mean([r["coexp_vs_atac"] for r in rows])), 4) if rows else None,
               note="FM per-type readout GPU-deferred (embedding+attention)")
json.dump(summary, open(f"{OUT}/pertype_coexp_v2.json", "w"), indent=2)
log(f"=== PER-TYPE (gated) mean coexp_vs_atac: {summary['mean_coexp_vs_atac']} ===")
log(f"SAVED {OUT}/pertype_coexp_v2.json")
