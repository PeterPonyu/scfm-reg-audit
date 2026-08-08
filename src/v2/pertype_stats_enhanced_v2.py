#!/usr/bin/env python
"""
scfm-reg-audit v2 — per-cell-type confound control + TF-block bootstrap CI (the Fig 5 numbers
never had confound control at all, let alone uncertainty quantification -- this closes both gaps
at once). Per-type degree confounds (TF out-degree, target in-degree) are derived from that TYPE's
own truth graph, matching how the pooled analysis derives them from the pooled consensus.

Scope: brain (5 types x coexp/geneformer-embed/geneformer-attn/scfoundation = 20 combos) + PBMC
(7 types x coexp/geneformer-embed/geneformer-attn = 21 combos; PBMC per-type scFoundation was never
computed, only pooled -- not in scope here). Bootstrap only (B=2000) -- full Mantel permutation at
this combo count (~41) would take hours; bootstrap CI is what was asked for and is the primary
uncertainty statement per the pooled analysis's own conclusion (Mantel null is anti-conservative
for this graph structure anyway, so skipping it here loses little).
"""
import os, json, hashlib, time, numpy as np, anndata as ad, pyfaidx
from scipy.stats import rankdata, spearmanr
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
PROM = 2000; W = 500
NBOOT = int(os.environ.get("NBOOT", "2000"))
MINCELL = 150

man = json.load(open(f"{ROOT}/data/manifest/shared_genes.v2.json")); genes = man["genes"]
gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
det = man["detection"]
MARKERS_GENES = {"MOBP", "MOG", "PLP1", "MBP", "ST18", "CTNNA3", "SLC17A7", "SATB2", "RBFOX3", "NRGN", "SLC17A6",
                 "GAD1", "GAD2", "SLC32A1", "DLX1", "AQP4", "GFAP", "SLC1A2", "ALDH1L1", "CSF1R", "P2RY12",
                 "CX3CR1", "C1QA", "C1QB", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "VCAN", "CLDN5", "FLT1",
                 "PECAM1", "PDGFRB", "RGS5"}


def build_confounds_cached(atac_file, cache={}):
    if atac_file in cache: return cache[atac_file]
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
    cache[atac_file] = (peakcount, genelen, gc, detv)
    return cache[atac_file]


def z(v):
    v = v.astype(float); return (v - v.mean()) / (v.std() + 1e-9)


def resid(y, C):
    C1 = np.c_[np.ones(len(y)), C]; b = np.linalg.lstsq(C1, y, rcond=None)[0]; return y - C1 @ b


def pcorr(x, y, C):
    return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])


def confound_stat(fm, atac, coexp, ii, jj, peakcount, genelen, gc, detv, tf_outdeg, atac_indeg, control_coexp=True):
    conf = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj]), z(tf_outdeg[ii]), z(atac_indeg[jj])]
    C = np.c_[rankdata(coexp), conf] if control_coexp else conf
    return pcorr(rankdata(fm), rankdata(atac), C)


def analyze_type(tissue, ctype, tf_rows, G_atac_type, G_co, G_fm, atac_file, readout_label, n_cells, mask_markers=False):
    control_coexp = readout_label != "coexp_vs_atac"   # fm IS coexp for that row -- self-partialling is degenerate
    ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj
    if mask_markers:
        mgi = np.array([g not in MARKERS_GENES for g in genes])
        m = m & mgi[ii] & mgi[jj]
    ii, jj = ii[m], jj[m]
    peakcount, genelen, gc, detv = build_confounds_cached(atac_file)
    tf_outdeg = (G_atac_type > 0).sum(1).astype(float); atac_indeg = (G_atac_type > 0).sum(0).astype(float)
    fm_v, atac_v, co_v = G_fm[ii, jj], G_atac_type[ii, jj], G_co[ii, jj]

    observed = confound_stat(fm_v, atac_v, co_v, ii, jj, peakcount, genelen, gc, detv, tf_outdeg, atac_indeg, control_coexp)
    marginal = float(spearmanr(fm_v, atac_v).statistic)

    tf_list = np.unique(tf_rows); rng = np.random.default_rng(hash(tissue + ctype + readout_label) % (2**31))
    tf_pos = {t: np.where(tf_rows == t)[0] for t in tf_list}
    nboot_vals = []
    for _ in range(NBOOT):
        samp_tf = rng.choice(tf_list, size=len(tf_list), replace=True)
        pos = np.concatenate([tf_pos[t] for t in samp_tf])
        try:
            v = confound_stat(fm_v[pos], atac_v[pos], co_v[pos], ii[pos], jj[pos], peakcount, genelen, gc, detv, tf_outdeg, atac_indeg, control_coexp)
            if np.isfinite(v): nboot_vals.append(v)
        except Exception: continue
    nboot_vals = np.array(nboot_vals)
    ci_lo, ci_hi = (np.percentile(nboot_vals, [2.5, 97.5]) if len(nboot_vals) > 10 else (np.nan, np.nan))
    row = dict(tissue=tissue, cell_type=ctype, readout=readout_label, n_cells=int(n_cells),
              marginal=round(marginal, 4), confound_controlled=round(float(observed), 4),
              bootstrap_ci95=[round(float(ci_lo), 4), round(float(ci_hi), 4)],
              bootstrap_se=round(float(nboot_vals.std()), 4) if len(nboot_vals) else None,
              crosses_zero=bool(ci_lo <= 0 <= ci_hi) if np.isfinite(ci_lo) else None)
    log(f"  [{tissue}/{ctype}/{readout_label}] n={n_cells} marginal={marginal:+.4f} confound_ctrl={observed:+.4f} "
        f"CI=[{ci_lo:+.4f},{ci_hi:+.4f}]")
    return row


if __name__ == "__main__":
    results = []

    # ---- BRAIN: 5 types x {coexp, geneformer_embed, geneformer_attn, scfoundation} ----
    Zb = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
    types_b = [str(t) for t in Zb["types"]]; tf_b = np.array(Zb["tf_rows"])
    ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
    ncells_b = json.load(open(f"{OUT}/pertype_fm_v2.json"))["per_type"]
    ncells_b = {r["cell_type"]: r["n"] for r in ncells_b}
    log("=== BRAIN per-type ===")
    for t in types_b:
        fpath = f"{OUT}/brain_fmgraphs_{t}.npz"
        if not os.path.exists(fpath): log(f"  skip {t}: no cached FM graph"); continue
        Ga = Zb[f"G_{t}"].astype(np.float32)
        F = np.load(fpath); n = ncells_b.get(t, 0)
        if n < MINCELL: log(f"  skip {t}: n={n} < {MINCELL}"); continue
        results.append(analyze_type("brain", t, tf_b, Ga, F["co"], F["co"], ATAC_B, "coexp_vs_atac", n, mask_markers=True))
        results.append(analyze_type("brain", t, tf_b, Ga, F["co"], F["gf"], ATAC_B, "geneformer_embed", n, mask_markers=True))
        results.append(analyze_type("brain", t, tf_b, Ga, F["co"], F["at"], ATAC_B, "geneformer_attn", n, mask_markers=True))
        sfpath = f"{OUT}/brain_scfgraphs_{t}.npz"
        if os.path.exists(sfpath):
            S = np.load(sfpath)
            results.append(analyze_type("brain", t, tf_b, Ga, F["co"], S["scf"], ATAC_B, "scfoundation", n, mask_markers=True))

    # ---- PBMC: 7 types x {coexp, geneformer_embed, geneformer_attn} ----
    Zp = np.load(f"{OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=True)
    types_p = [str(t) for t in Zp["types"]]; tf_p = np.array(Zp["tf_rows"])
    ATAC_P = f"{ROOT}/data/multiome/pbmc10k_atac.h5ad"
    ncells_p = json.load(open(f"{OUT}/pbmc_eval_v2.json"))["per_type_coexp"]
    ncells_p = {r["cell_type"]: r["n"] for r in ncells_p}
    log("=== PBMC per-type ===")
    for t in types_p:
        fpath = f"{OUT}/pbmc_fmgraphs_{t}.npz"
        if not os.path.exists(fpath): log(f"  skip {t}: no cached FM graph"); continue
        Ga = Zp[f"G_{t}"].astype(np.float32)
        F = np.load(fpath); n = ncells_p.get(t, 0)
        if n < MINCELL: log(f"  skip {t}: n={n} < {MINCELL}"); continue
        results.append(analyze_type("pbmc", t, tf_p, Ga, F["co"], F["co"], ATAC_P, "coexp_vs_atac", n))
        results.append(analyze_type("pbmc", t, tf_p, Ga, F["co"], F["gf"], ATAC_P, "geneformer_embed", n))
        results.append(analyze_type("pbmc", t, tf_p, Ga, F["co"], F["at"], ATAC_P, "geneformer_attn", n))

    json.dump(results, open(f"{OUT}/pertype_stats_enhanced_v2.json", "w"), indent=2)
    log(f"=== DONE: {len(results)} combos ===")
    n_crosses = sum(1 for r in results if r["crosses_zero"])
    log(f"CI crosses zero in {n_crosses}/{len(results)} combos")
    log(f"SAVED {OUT}/pertype_stats_enhanced_v2.json")
