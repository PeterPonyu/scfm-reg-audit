#!/usr/bin/env python
"""
scfm-reg-audit v2 -- POSITIVE CONTROL + minimum-detectable-effect (MDE) / power analysis.

Why this exists: the headline is a NULL (no FM recovers the regulatory truth beyond co-expression,
confound-controlled). A null is only informative if the test could have detected a real effect had
one existed -- and our own TF-block bootstrap finding (SE ~0.045) makes the reviewer's first
objection sharp: "your test can't resolve anything, so the null says nothing." This script answers
it head-on by *injecting* a controlled, known dose of genuine regulatory signal and showing the
exact same pipeline recovers it dose-dependently, then reads off the smallest dose whose 95%
TF-block bootstrap CI clears zero (the minimum detectable effect). Real FMs are then placed on that
calibration curve: if every FM sits below the MDE, the null is quantified, not just asserted.

Design of the injected signal (this is the crux of making alpha interpretable):
  atac_resid = residual of rank(truth) after regressing out [co-expression, structural confounds].
  This is EXACTLY the quantity an FM is tested for -- regulatory structure orthogonal to
  co-expression and to peak density / gene length / GC% / detection / network degree. We inject a
  fraction alpha of it into an otherwise-noise "model" graph:
      G_synth(alpha) = alpha * z(atac_resid) + (1-alpha) * z(gaussian noise)
  and measure the confound-controlled partial rho of G_synth vs the truth, controlling co-expression
  + confounds -- the identical statistic used on every real FM. So alpha is directly "how much
  genuine orthogonal regulatory signal the model encodes," on the same axis as the FM results.

Pure numpy on cached graphs (results/v2/*.npz) -- no GPU, no FM re-embedding.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, pyfaidx
from scipy.stats import rankdata
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
PROM = 2000; W = 500
NBOOT = int(os.environ.get("NBOOT", "1000"))
ALPHAS = [0.0, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.75, 1.0]

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

def pcorr(x, y, C):
    return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])


def confound_controlled(fm_v, atac_v, co_v, conf_mat, control_coexp=True):
    """Identical statistic to stats_enhanced_v2.confound_stat, but conf_mat is prebuilt (edge-indexed)
    so the injection sweep doesn't rebuild it every alpha."""
    C = np.c_[rankdata(co_v), conf_mat] if control_coexp else conf_mat
    return pcorr(rankdata(fm_v), rankdata(atac_v), C)


def run_tissue(tag, tf_rows, G_atac, G_co, atac_file, fm_observed):
    ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
    peakcount, genelen, gc, detv = build_confounds(atac_file)
    tf_outdeg = (G_atac > 0).sum(1).astype(float); atac_indeg = (G_atac > 0).sum(0).astype(float)
    conf_mat = np.c_[z(peakcount[jj]), z(genelen[jj]), z(detv[jj]), z(gc[jj]), z(tf_outdeg[ii]), z(atac_indeg[jj])]
    atac_v, co_v = G_atac[ii, jj], G_co[ii, jj]

    # the injectable signal: truth's rank, orthogonalized to coexp + confounds (what an FM is tested for)
    atac_resid = resid(rankdata(atac_v), np.c_[rankdata(co_v), conf_mat])
    atac_resid_z = z(atac_resid)

    rng = np.random.default_rng((hash(tag + "power") % (2**31)))
    noise = rng.standard_normal(len(atac_v))

    tf_list = np.unique(tf_rows); tf_pos = {t: np.where(tf_rows == t)[0] for t in tf_list}

    def boot_ci(synth):
        vals = []
        for _ in range(NBOOT):
            samp = rng.choice(tf_list, size=len(tf_list), replace=True)
            pos = np.concatenate([tf_pos[t] for t in samp])
            try:
                v = confound_controlled(synth[pos], atac_v[pos], co_v[pos], conf_mat[pos])
                if np.isfinite(v): vals.append(v)
            except Exception: continue
        vals = np.array(vals)
        return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), float(vals.std())

    curve = []
    log(f"=== {tag}: injection dose-response (B={NBOOT}) ===")
    for a in ALPHAS:
        synth = a * atac_resid_z + (1 - a) * z(noise)
        obs = confound_controlled(synth, atac_v, co_v, conf_mat)
        lo, hi, se = boot_ci(synth)
        curve.append(dict(alpha=a, observed=round(obs, 4), ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                          se=round(se, 4), clears_zero=bool(lo > 0)))
        log(f"  alpha={a:.2f}  rho={obs:+.4f}  CI=[{lo:+.4f},{hi:+.4f}]  clears0={lo>0}")

    # empirical MDE: smallest alpha whose bootstrap CI lower bound > 0
    mde_alpha = next((c["alpha"] for c in curve if c["clears_zero"]), None)
    # the recovered rho at that alpha (the minimum detectable confound-controlled effect on the rho axis)
    mde_rho = next((c["observed"] for c in curve if c["clears_zero"]), None)
    # analytic MDE cross-check: 80% power, two-sided 0.05 -> effect ~ (1.96+0.84)*SE at that scale
    se0 = curve[0]["se"]  # SE near the null, the relevant sampling noise
    mde_analytic = round(2.80 * se0, 4)

    # place real FMs on the curve: implied alpha by linear interpolation of observed(alpha)
    alphas_arr = np.array([c["alpha"] for c in curve]); rhos_arr = np.array([c["observed"] for c in curve])
    order = np.argsort(rhos_arr)
    fm_on_curve = {}
    for name, rho in fm_observed.items():
        implied_a = float(np.interp(rho, rhos_arr[order], alphas_arr[order]))
        fm_on_curve[name] = dict(observed_rho=rho, implied_alpha=round(implied_a, 4),
                                 above_mde=bool(mde_rho is not None and rho >= mde_rho))
        log(f"  FM {name:28s} rho={rho:+.4f} -> implied alpha={implied_a:.3f}  above_MDE={mde_rho is not None and rho >= mde_rho}")

    log(f"  MDE (empirical): alpha>={mde_alpha}  (recovered rho={mde_rho});  analytic 80%-power MDE~{mde_analytic}")
    return dict(tag=tag, curve=curve, mde_alpha_empirical=mde_alpha, mde_rho_empirical=mde_rho,
                se_at_null=se0, mde_rho_analytic_80pct=mde_analytic, fm_on_curve=fm_on_curve)


if __name__ == "__main__":
    enh = {r["label"]: r["observed"] for r in json.load(open(f"{OUT}/stats_enhanced_v2.json")) if r["tag"] == "brain"}
    enh_p = {r["label"]: r["observed"] for r in json.load(open(f"{OUT}/stats_enhanced_v2.json")) if r["tag"] == "pbmc"}

    out = {}
    # ---- BRAIN (headline figure) ----
    Zb = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
    types_b = [str(t) for t in Zb["types"]]; tf_b = np.array(Zb["tf_rows"])
    Gb = np.mean([Zb[f"G_{t}"] for t in types_b], axis=0).astype(np.float32)
    Fb = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")
    ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
    fm_b = {"Geneformer (embed)": enh.get("geneformer_embed"),
            "scFoundation": enh.get("scfoundation"),
            "UCE": enh.get("uce"),
            "co-expression": enh.get("coexp_vs_atac(confounds_only)")}
    fm_b = {k: v for k, v in fm_b.items() if v is not None}
    out["brain"] = run_tissue("brain", tf_b, Gb, Fb["co"], ATAC_B, fm_b)

    # ---- PBMC (robustness) ----
    Zp = np.load(f"{OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=True)
    types_p = [str(t) for t in Zp["types"]]; tf_p = np.array(Zp["tf_rows"])
    Gp = np.mean([Zp[f"G_{t}"] for t in types_p], axis=0).astype(np.float32)
    Fp = np.load(f"{OUT}/pbmc_fmgraphs_pooled.npz")
    ATAC_P = f"{ROOT}/data/multiome/pbmc10k_atac.h5ad"
    fm_p = {"Geneformer (embed)": enh_p.get("geneformer_embed"),
            "scFoundation": enh_p.get("scfoundation")}
    fm_p = {k: v for k, v in fm_p.items() if v is not None}
    out["pbmc"] = run_tissue("pbmc", tf_p, Gp, Fp["co"], ATAC_P, fm_p)

    json.dump(out, open(f"{OUT}/power_analysis_v2.json", "w"), indent=2)
    log(f"SAVED {OUT}/power_analysis_v2.json")
