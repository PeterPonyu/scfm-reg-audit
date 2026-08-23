#!/usr/bin/env python
"""
scfm-reg-audit v2 — confound ABLATION: which covariate(s) actually drive the raw->confound-
controlled drop? Answers a specific open question flagged in the statistical rigor audit: the
confound set includes the truth graph's OWN degree (TF out-degree, target in-degree, both computed
from G_ATAC) -- defensible, but does controlling for a graph's own degree over-correct and remove
real signal along with the artifact? Reports partial rho under nested covariate sets:
  (0) none (marginal)              (1) coexp only
  (2) coexp + non-degree confounds (peakcount, genelen, gc, detection -- gene-identity properties,
      independent of G_ATAC's edge values)
  (3) coexp + non-degree + degree confounds (tf_outdeg, atac_indeg -- the FULL set used everywhere
      else in this pipeline)
If (2)->(3) does most of the work, degree IS the dominant confound (worth flagging explicitly,
not hiding). If (1)->(2) already kills most of the raw signal, the non-degree confounds are doing
the heavy lifting and the "controlling for the graph's own degree" concern is less load-bearing.
CPU-only, reuses cached graphs.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, pyfaidx
from scipy.stats import rankdata
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
PROM = 2000; W = 500

man = json.load(open(f"{ROOT}/data/manifest/shared_genes.v2.json")); genes = man["genes"]
gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
det = man["detection"]


def build_confounds(atac_file):
    gco = {}
    for ln in open(COORDS):
        c, s, e, st, nm = ln.rstrip("\n").split("\t")
        if nm not in gidx or nm in gco: continue
        s, e = int(s), int(e); lo = s - PROM if st == "+" else s; hi = e if st == "+" else e + PROM
        gco[nm] = (c, lo, hi, abs(e - s))
    Av = ad.read_h5ad(atac_file, backed="r"); peaks = [str(p) for p in Av.var_names]
    pchr = np.array([p.split(":")[0] for p in peaks]); pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks])
    pmid = (pse[:, 0] + pse[:, 1]) // 2
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
    detv = np.array([det.get(g, 0.0) for g in genes])
    return peakcount, genelen, gc, detv


def z(v):
    v = v.astype(float); return (v - v.mean()) / (v.std() + 1e-9)


def resid(y, C):
    C1 = np.c_[np.ones(len(y)), C]; b = np.linalg.lstsq(C1, y, rcond=None)[0]; return y - C1 @ b


def pcorr(x, y, C=None):
    if C is None or C.shape[1] == 0: return float(np.corrcoef(x, y)[0, 1])
    return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])


def ablate(tf_rows, G_atac, G_co, G_fm, atac_file, label, control_coexp=True):
    """control_coexp=False when `fm` IS co-expression itself -- controlling coexp for coexp is
    degenerate (residual ~0 by construction), so that row skips straight to the structural
    confounds without a co-expression-partialling stage."""
    ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
    peakcount, genelen, gc, detv = build_confounds(atac_file)
    tf_outdeg = (G_atac > 0).sum(1).astype(float); atac_indeg = (G_atac > 0).sum(0).astype(float)
    fm, atac, co = G_fm[ii, jj], G_atac[ii, jj], G_co[ii, jj]
    rfm, ratac, rco = rankdata(fm), rankdata(atac), rankdata(co)
    nondeg = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj])]
    deg = np.c_[z(tf_outdeg[ii]), z(atac_indeg[jj])]

    if control_coexp:
        stages = dict(
            marginal=pcorr(rfm, ratac),
            coexp_only=pcorr(rfm, ratac, rco.reshape(-1, 1)),
            coexp_plus_nondegree=pcorr(rfm, ratac, np.c_[rco, nondeg]),
            coexp_plus_full=pcorr(rfm, ratac, np.c_[rco, nondeg, deg]),
        )
        log(f"  [{label}] marginal={stages['marginal']:+.4f} -> +coexp={stages['coexp_only']:+.4f} "
            f"-> +nondeg={stages['coexp_plus_nondegree']:+.4f} -> +degree(full)={stages['coexp_plus_full']:+.4f}")
    else:
        stages = dict(
            marginal=pcorr(rfm, ratac),
            coexp_only=None,   # N/A: fm IS coexp here, self-partialling is degenerate
            coexp_plus_nondegree=pcorr(rfm, ratac, nondeg),
            coexp_plus_full=pcorr(rfm, ratac, np.c_[nondeg, deg]),
        )
        log(f"  [{label}] marginal={stages['marginal']:+.4f} -> (coexp N/A, self-test) "
            f"-> +nondeg={stages['coexp_plus_nondegree']:+.4f} -> +degree(full)={stages['coexp_plus_full']:+.4f}")
    return dict(label=label, control_coexp=control_coexp,
               **{k: (round(v, 4) if v is not None else None) for k, v in stages.items()})


if __name__ == "__main__":
    results = []
    Zb = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
    types_b = [str(t) for t in Zb["types"]]; tf_b = np.array(Zb["tf_rows"])
    Gb = np.mean([Zb[f"G_{t}"] for t in types_b], axis=0).astype(np.float32)
    Fb = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")
    ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
    log("=== BRAIN confound ablation ===")
    results.append(ablate(tf_b, Gb, Fb["co"], Fb["co"], ATAC_B, "coexp_vs_atac", control_coexp=False))
    results.append(ablate(tf_b, Gb, Fb["co"], Fb["gf"], ATAC_B, "geneformer_embed"))
    Sb = np.load(f"{OUT}/G_scf_pooled.npz")["G"]; results.append(ablate(tf_b, Gb, Fb["co"], Sb, ATAC_B, "scfoundation"))
    Ub = np.load(f"{OUT}/G_uce_pooled.npz")["G"]; results.append(ablate(tf_b, Gb, Fb["co"], Ub, ATAC_B, "uce"))
    Kb = np.load(f"{OUT}/G_ko_v2.npz"); results.append(ablate(tf_b, Gb, Fb["co"], Kb["G_ko"], ATAC_B, "geneformer_ko_raw"))

    json.dump(results, open(f"{OUT}/confound_ablation_v2.json", "w"), indent=2)
    log(f"SAVED {OUT}/confound_ablation_v2.json")
