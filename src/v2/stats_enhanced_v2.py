#!/usr/bin/env python
"""
scfm-reg-audit v2 — statistical rigor pass on the confound-CONTROLLED headline numbers.

Closes 3 gaps flagged in review of the pipeline so far:
  1. No CI/SE anywhere -> block-bootstrap CI (resample TF nodes w/ replacement, B=2000) on every
     confound-controlled partial rho.
  2. The Mantel permutation test was only ever run on RAW (pre-confound) values. The headline
     "confound-controlled partial rho ~ 0" claim was never itself significance-tested -> gene-label
     Mantel permutation applied DIRECTLY to the confound-controlled statistic (N=1000), with the
     degree-based confounds (TF out-degree, target in-degree) re-derived under each permutation
     (they're computed FROM G_ATAC, so a label permutation reassigns them consistently:
     outdeg_perm[i] = outdeg_orig[perm[i]]); the non-degree confounds (peak count, gene length, GC%,
     detection rate) are properties of gene identity independent of G_ATAC's edge structure and stay
     fixed under permutation.
  3. No multiple-testing correction across the ~10 tests run -> BH-FDR q-values.
  4. A DEGREE-PRESERVING null (stronger than gene-label permutation) as a secondary check on the
     headline brain results: independently permute each TF row's edge weights across targets. This
     preserves each TF's exact out-degree/value multiset by construction (nothing like a rank-1
     collapse could inflate it) and tests whether the SPECIFIC target assignment carries signal
     beyond aggregate magnitude, not just whether gene labels are exchangeable.

Reuses cached graphs (results/v2/*.npz) -- no GPU/CPU FM re-embedding needed.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, pyfaidx
from scipy.stats import rankdata
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
PROM = 2000; W = 500
NBOOT = int(os.environ.get("NBOOT", "2000"))
NPERM = int(os.environ.get("NPERM", "1000"))
NDEG = int(os.environ.get("NDEG", "500"))

man = json.load(open(f"{ROOT}/data/manifest/shared_genes.v2.json")); genes = man["genes"]
gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
det = man["detection"]


def build_confounds(atac_file):
    """peakcount, genelen, gc (gene identity, NOT G_atac-derived -> fixed under permutation)."""
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


def pcorr(x, y, C):
    return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])


def confound_stat(fm, atac, coexp, ii, jj, peakcount, genelen, gc, detv, tf_outdeg, atac_indeg, control_coexp=True):
    """control_coexp=True: partial FM-vs-truth controlling co-expression + structural confounds
    (the FM test). control_coexp=False: partial [fm]-vs-truth controlling structural confounds
    ONLY (used when `fm` IS co-expression itself -- controlling coexp for coexp is degenerate)."""
    conf = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj]), z(tf_outdeg[ii]), z(atac_indeg[jj])]
    C = np.c_[rankdata(coexp), conf] if control_coexp else conf
    return pcorr(rankdata(fm), rankdata(atac), C)


def analyze(tag, tf_rows, G_atac, G_co, G_fm, atac_file, label, control_coexp=True):
    ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
    peakcount, genelen, gc, detv = build_confounds(atac_file)
    tf_outdeg_full = (G_atac > 0).sum(1).astype(float); atac_indeg_full = (G_atac > 0).sum(0).astype(float)
    fm_v, atac_v, co_v = G_fm[ii, jj], G_atac[ii, jj], G_co[ii, jj]

    observed = confound_stat(fm_v, atac_v, co_v, ii, jj, peakcount, genelen, gc, detv, tf_outdeg_full, atac_indeg_full, control_coexp)
    log(f"  [{label}] observed confound-controlled partial rho = {observed:.4f}")

    # ---- 1. block-bootstrap CI (resample TF nodes with replacement) ----
    t0 = time.time()
    tf_list = np.unique(tf_rows); nboot_vals = []
    rng = np.random.default_rng(hash(tag + label) % (2**31))
    tf_pos = {t: np.where(tf_rows == t)[0] for t in tf_list}   # positions in ii/jj for each TF
    for _ in range(NBOOT):
        samp_tf = rng.choice(tf_list, size=len(tf_list), replace=True)
        pos = np.concatenate([tf_pos[t] for t in samp_tf])
        try:
            v = confound_stat(fm_v[pos], atac_v[pos], co_v[pos], ii[pos], jj[pos], peakcount, genelen, gc, detv, tf_outdeg_full, atac_indeg_full, control_coexp)
            if np.isfinite(v): nboot_vals.append(v)
        except Exception: continue
    nboot_vals = np.array(nboot_vals)
    ci_lo, ci_hi = np.percentile(nboot_vals, [2.5, 97.5])
    log(f"  [{label}] bootstrap CI ({len(nboot_vals)}/{NBOOT} ok, {time.time()-t0:.0f}s): [{ci_lo:.4f}, {ci_hi:.4f}]")

    # ---- 2. gene-label Mantel permutation test on the CONFOUND-CONTROLLED statistic ----
    t0 = time.time(); null = []
    rng2 = np.random.default_rng((hash(tag + label + "perm")) % (2**31))
    for _ in range(NPERM):
        perm = rng2.permutation(Ng)
        atac_perm = G_atac[perm[ii], perm[jj]]
        outdeg_perm = tf_outdeg_full[perm]; indeg_perm = atac_indeg_full[perm]
        v = confound_stat(fm_v, atac_perm, co_v, ii, jj, peakcount, genelen, gc, detv, outdeg_perm, indeg_perm, control_coexp)
        null.append(v)
    null = np.array(null)
    zscore = (observed - null.mean()) / (null.std() + 1e-9)
    p_perm = (np.sum(np.abs(null) >= abs(observed)) + 1) / (NPERM + 1)
    log(f"  [{label}] Mantel-on-confound-controlled ({time.time()-t0:.0f}s): z={zscore:.2f} p={p_perm:.4f}")

    return dict(tag=tag, label=label, observed=round(float(observed), 4),
               bootstrap_ci95=[round(float(ci_lo), 4), round(float(ci_hi), 4)],
               bootstrap_n_ok=int(len(nboot_vals)), bootstrap_se=round(float(nboot_vals.std()), 4),
               mantel_confound_controlled=dict(null_mean=round(float(null.mean()), 4), null_sd=round(float(null.std()), 4),
                                               z=round(float(zscore), 2), p_perm=round(float(p_perm), 4)))


def degree_preserving_null(tf_rows, G_atac, G_co, G_fm, atac_file, label, ndeg=NDEG, control_coexp=True):
    """Row-shuffle null: independently permute each TF row's edge weights across targets. Preserves
    each TF's exact out-degree/value multiset by construction; tests whether the SPECIFIC
    target assignment carries signal, a stronger/different null than gene-label permutation."""
    ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
    peakcount, genelen, gc, detv = build_confounds(atac_file)
    tf_outdeg_full = (G_atac > 0).sum(1).astype(float)
    fm_v, atac_v, co_v = G_fm[ii, jj], G_atac[ii, jj], G_co[ii, jj]
    atac_indeg_full = (G_atac > 0).sum(0).astype(float)
    observed = confound_stat(fm_v, atac_v, co_v, ii, jj, peakcount, genelen, gc, detv, tf_outdeg_full, atac_indeg_full, control_coexp)

    t0 = time.time(); rng = np.random.default_rng((hash(label + "deg")) % (2**31)); null = []
    for _ in range(ndeg):
        Gp = G_atac.copy()
        for t in tf_rows:
            row = Gp[t].copy(); nz = np.where(np.arange(Ng) != t)[0]
            perm_vals = rng.permutation(row[nz]); Gp[t, nz] = perm_vals
        atac_p = Gp[ii, jj]; indeg_p = (Gp > 0).sum(0).astype(float)  # outdeg per TF unchanged by row-shuffle
        v = confound_stat(fm_v, atac_p, co_v, ii, jj, peakcount, genelen, gc, detv, tf_outdeg_full, indeg_p, control_coexp)
        null.append(v)
    null = np.array(null)
    zscore = (observed - null.mean()) / (null.std() + 1e-9)
    p_perm = (np.sum(np.abs(null) >= abs(observed)) + 1) / (ndeg + 1)
    log(f"  [{label}] DEGREE-PRESERVING null ({ndeg} shuffles, {time.time()-t0:.0f}s): observed={observed:.4f} z={zscore:.2f} p={p_perm:.4f}")
    return dict(label=label, observed=round(float(observed), 4), n_shuffles=ndeg,
               null_mean=round(float(null.mean()), 4), null_sd=round(float(null.std()), 4),
               z=round(float(zscore), 2), p_perm=round(float(p_perm), 4))


def bh_fdr(pvals):
    pvals = np.array(pvals); n = len(pvals); order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


if __name__ == "__main__":
    results = []

    # ---- BRAIN ----
    Zb = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
    types_b = [str(t) for t in Zb["types"]]; tf_b = np.array(Zb["tf_rows"])
    Gb = np.mean([Zb[f"G_{t}"] for t in types_b], axis=0).astype(np.float32)
    Fb = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")
    ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
    log("=== BRAIN ===")
    results.append(analyze("brain", tf_b, Gb, Fb["co"], Fb["co"], ATAC_B, "coexp_vs_atac(confounds_only)", control_coexp=False))
    results.append(analyze("brain", tf_b, Gb, Fb["co"], Fb["gf"], ATAC_B, "geneformer_embed"))
    Sb = np.load(f"{OUT}/G_scf_pooled.npz")["G"]; results.append(analyze("brain", tf_b, Gb, Fb["co"], Sb, ATAC_B, "scfoundation"))
    Ub = np.load(f"{OUT}/G_uce_pooled.npz")["G"]; results.append(analyze("brain", tf_b, Gb, Fb["co"], Ub, ATAC_B, "uce"))
    Kb = np.load(f"{OUT}/G_ko_v2.npz"); results.append(analyze("brain", tf_b, Gb, Fb["co"], Kb["G_ko"], ATAC_B, "geneformer_ko_raw"))
    results.append(analyze("brain", tf_b, Gb, Fb["co"], Kb["G_ko_ctrl"], ATAC_B, "geneformer_ko_posctrl"))

    # ---- PBMC ----
    Zp = np.load(f"{OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=True)
    types_p = [str(t) for t in Zp["types"]]; tf_p = np.array(Zp["tf_rows"])
    Gp = np.mean([Zp[f"G_{t}"] for t in types_p], axis=0).astype(np.float32)
    Fp = np.load(f"{OUT}/pbmc_fmgraphs_pooled.npz")
    ATAC_P = f"{ROOT}/data/multiome/pbmc10k_atac.h5ad"
    log("=== PBMC ===")
    results.append(analyze("pbmc", tf_p, Gp, Fp["co"], Fp["gf"], ATAC_P, "geneformer_embed"))
    results.append(analyze("pbmc", tf_p, Gp, Fp["co"], Fp["at"], ATAC_P, "geneformer_attn"))
    Sp = np.load(f"{OUT}/G_scf_pbmc_pooled.npz")["G"]; results.append(analyze("pbmc", tf_p, Gp, Fp["co"], Sp, ATAC_P, "scfoundation"))

    # ---- BH-FDR across the full family ----
    pvals = [r["mantel_confound_controlled"]["p_perm"] for r in results]
    q = bh_fdr(pvals)
    for r, qi in zip(results, q): r["bh_q"] = round(float(qi), 4)

    json.dump(results, open(f"{OUT}/stats_enhanced_v2.json", "w"), indent=2)
    log("=== SUMMARY (confound-controlled, bootstrap CI + Mantel-on-adjusted-stat + BH-q) ===")
    for r in results:
        m = r["mantel_confound_controlled"]
        log(f"  {r['tag']:6s} {r['label']:32s} obs={r['observed']:+.4f} CI95=[{r['bootstrap_ci95'][0]:+.4f},{r['bootstrap_ci95'][1]:+.4f}] z={m['z']:+.2f} p={m['p_perm']:.4f} q={r['bh_q']:.4f}")

    # ---- degree-preserving null on the headline brain results ----
    log("=== DEGREE-PRESERVING NULL (secondary, stronger; headline brain results) ===")
    deg_results = []
    deg_results.append(degree_preserving_null(tf_b, Gb, Fb["co"], Fb["gf"], ATAC_B, "brain_geneformer_embed"))
    deg_results.append(degree_preserving_null(tf_b, Gb, Fb["co"], Sb, ATAC_B, "brain_scfoundation"))
    deg_results.append(degree_preserving_null(tf_b, Gb, Fb["co"], Ub, ATAC_B, "brain_uce"))
    deg_results.append(degree_preserving_null(tf_b, Gb, Fb["co"], Fb["co"], ATAC_B, "brain_coexp_self", control_coexp=False))
    json.dump(deg_results, open(f"{OUT}/degree_preserving_null_v2.json", "w"), indent=2)
    log(f"SAVED {OUT}/stats_enhanced_v2.json + {OUT}/degree_preserving_null_v2.json")
