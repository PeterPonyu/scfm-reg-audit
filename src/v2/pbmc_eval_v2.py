#!/usr/bin/env python
"""
scfm-reg-audit v2 — PBMC Multiome: 2nd-TISSUE replication + PAIRED-CELL CALIBRATION, one pass.

This dataset is a single 10x cellranger-arc combined matrix: RNA and ATAC come from the SAME
10,970 cells (barcode-matched), unlike the AD-brain analysis (2 independent studies matched only
by cell-type label). So this run does double duty:
  (a) 2nd tissue (blood, non-brain, non-AD) — replicates the FM-vs-regulation verdict.
  (b) Paired-cell calibration — the RNA/ATAC pairing removes cross-study/cross-annotation mismatch,
      testing whether the type-pooled 'unpaired' methodology used for GSE174367+ad_hm gives the
      same qualitative verdict when there is no cross-study confound to begin with (DESIGN §6).

Same frozen manifest, same test suite as crossmodal_v2.py + confound_regression_v2.py (marginal,
partial|coexp, confound-controlled partial, Mantel null), pooled + per-type.
"""
import os, json, hashlib, time, numpy as np, anndata as ad, scipy.sparse as sp
from scipy.stats import spearmanr, rankdata
import fm_readout as fr
import pbmc_cache
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"; HG38 = f"{ROOT}/data/genome/hg38.fa"
RNA = f"{ROOT}/data/multiome/pbmc10k_rna.h5ad"
ATAC = f"{ROOT}/data/multiome/pbmc10k_atac.h5ad"
PROMOTER = 2000; W = 500; MIN_CELLS = 150
POOL_CAP = int(os.environ.get("POOL_CAP", "4000"))
POOLED_ONLY = os.environ.get("PBMC_POOLED_ONLY", "0") == "1"
PREFLIGHT_ONLY = os.environ.get("PBMC_PREFLIGHT_ONLY", "0") == "1"
PREFLIGHT_REPORT = os.environ.get("PBMC_PREFLIGHT_REPORT")
BATCH = int(os.environ.get("FM_BATCH", "4"))

man = json.load(open(MANI)); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
manifest_sha = hashlib.sha256(("\n".join(genes)).encode()).hexdigest()
assert manifest_sha == man["sha256"]
det = man["detection"]
Z = np.load(f"{OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=True)
zgenes = [str(g) for g in Z["genes"]]; assert zgenes == genes
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac_pooled = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
log(f"G_ATAC PBMC {G_atac_pooled.shape} over {len(types)} types: {types}")

# ---- confound covariates (PBMC's own peak set) ----
confound_cache = os.environ.get("PBMC_CONFOUND_CACHE")
if confound_cache:
    peakcount, genelen, detv, gc = pbmc_cache.load_confound_cache(
        confound_cache, genes, manifest_sha)
    log(f"confound covariates loaded from {confound_cache}")
else:
    import pyfaidx
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
    detv = np.array([det.get(g, 0.0) for g in genes])
tf_outdeg = (G_atac_pooled > 0).sum(1).astype(float); atac_indeg = (G_atac_pooled > 0).sum(0).astype(float)
log("confound covariates ready")

# ---- pair set + RNA ----
ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
def z(v): v = v.astype(float); return (v - v.mean()) / (v.std() + 1e-9)
def resid(y, C): C1 = np.c_[np.ones(len(y)), C]; b = np.linalg.lstsq(C1, y, rcond=None)[0]; return y - C1 @ b
def pcorr(x, y, C): return float(np.corrcoef(resid(x, C), resid(y, C))[0, 1])
def partial(x, y, z_): return pcorr(rankdata(x), rankdata(y), rankdata(z_).reshape(-1, 1))

A, Xc, Xl, _ = fr.load_norm(RNA)
import gzip as _gz
_meta = {}
with _gz.open(f"{ROOT}/data/multiome/pbmc_cell_meta.csv.gz", "rt") as _f:
    _hdr = _f.readline().rstrip("\n").split(","); _bi = _hdr.index("Barcode"); _ci = _hdr.index("Cell.Type")
    for _ln in _f:
        _p = _ln.rstrip("\n").split(","); _meta[_p[_bi]] = _p[_ci]
labels = np.array([_meta.get(bc, "NA") for bc in A.obs_names])
log(f"labels loaded: {int((labels!='NA').sum())}/{len(labels)} cells typed")
rsym = {str(s): k for k, s in enumerate(A.var_names)}
present = [g for g in genes if g in rsym]
ri = np.full(Ng, -1, int)
for g in present: ri[gidx[g]] = rsym[g]
valid = ri >= 0
log(f"RNA present for {valid.sum()}/{Ng} manifest genes")

def build_graphs(cell_ids, cache_path):
    if os.path.exists(cache_path):
        F = np.load(cache_path); return F["co"], F["gf"], F["at"]
    sub = np.where(valid)[0]
    Xl_sub = Xl[cell_ids][:, ri[sub]].tocsr(); Xc_sub = Xc[cell_ids][:, ri[sub]].tocsr()
    Gco_s = fr.gene_coexp(Xl_sub.toarray())
    rd = fr.FMReadout([genes[s] for s in sub], batch=BATCH)
    cells = np.arange(len(cell_ids))
    Gemb_s = fr.FMReadout.cos_graph(rd.geneformer(Xc_sub, cells))
    Gatt_s = rd.geneformer_attention(Xc_sub, cells)
    Gco = np.zeros((Ng, Ng), np.float32); Gemb = np.zeros((Ng, Ng), np.float32); Gatt = np.zeros((Ng, Ng), np.float32)
    Gco[np.ix_(sub, sub)] = Gco_s; Gemb[np.ix_(sub, sub)] = Gemb_s; Gatt[np.ix_(sub, sub)] = Gatt_s
    np.savez(cache_path, co=Gco, gf=Gemb, at=Gatt)
    return Gco, Gemb, Gatt

def build_pooled_scgpt_graphs(cell_ids, cache_path):
    rna_sha256 = pbmc_cache.sha256_file(RNA)
    cached = pbmc_cache.load_scgpt_cache(
        cache_path, cell_ids, genes, manifest_sha, 20260713, POOL_CAP, rna_sha256)
    if cached is not None:
        return cached
    sub = np.where(valid)[0]
    Xl_sub = Xl[cell_ids][:, ri[sub]].tocsr()
    Gco_s = fr.gene_coexp(Xl_sub.toarray())
    rd = fr.FMReadout([genes[s] for s in sub], batch=BATCH)
    cells = np.arange(len(cell_ids))
    log("scGPT embed…")
    Gsg_s = fr.FMReadout.cos_graph(rd.scgpt(Xl_sub, cells))
    Gco = np.zeros((Ng, Ng), np.float32)
    Gsg = np.zeros((Ng, Ng), np.float32)
    Gco[np.ix_(sub, sub)] = Gco_s
    Gsg[np.ix_(sub, sub)] = Gsg_s
    pbmc_cache.write_scgpt_cache(
        cache_path, Gco, Gsg, cell_ids, genes, manifest_sha, 20260713, POOL_CAP,
        rna_sha256)
    return Gco, Gsg

def run_test(a_, co_, fm_, conf_, tag):
    obs = dict(coexp_vs_atac=float(spearmanr(co_, a_).statistic), fm_vs_atac=float(spearmanr(fm_, a_).statistic),
               fm_vs_coexp=float(spearmanr(fm_, co_).statistic), fm_partial_coexp=partial(fm_, a_, co_),
               coexp_partial_confounds=pcorr(rankdata(co_), rankdata(a_), conf_),
               fm_partial_coexp_confounds=pcorr(rankdata(fm_), rankdata(a_), np.c_[rankdata(co_), conf_]))
    return {f"{tag}__{k}": round(v, 4) for k, v in obs.items()}

# ---- POOLED (valid manifest genes only, same-cell paired) ----
allcells = pbmc_cache.select_pool_cell_ids(A.shape[0], POOL_CAP, 20260713)
log(f"pooled FM readout on {len(allcells)} cells (capped, GPU-memory-shared with a concurrent job)")
Gco, Gemb, Gatt = build_graphs(allcells, f"{OUT}/pbmc_fmgraphs_pooled.npz")
sub = np.where(valid)[0]; iiv = ii[valid[ii] & valid[jj]]; jjv = jj[valid[ii] & valid[jj]]
conf = np.c_[z(peakcount[jjv]), z(genelen[jjv]), z(detv[jjv]), z(gc[jjv]), z(tf_outdeg[iiv]), z(atac_indeg[jjv])]
a_p = G_atac_pooled[iiv, jjv]; co_p = Gco[iiv, jjv]
res = dict(n_pairs=int(len(iiv)), n_tf=int(len(tf_rows)), tissue="PBMC10k_multiome_paired", types=types)
res.update(run_test(a_p, co_p, Gemb[iiv, jjv], conf, "embed"))
res.update(run_test(a_p, co_p, Gatt[iiv, jjv], conf, "attn"))
reference_metrics = os.environ.get("PBMC_REFERENCE_METRICS")
if reference_metrics:
    pbmc_cache.verify_reference_metrics(res, reference_metrics)
    log(f"reference metric gate passed against {reference_metrics}")
if PREFLIGHT_ONLY:
    if not reference_metrics:
        raise RuntimeError("PBMC_PREFLIGHT_ONLY requires PBMC_REFERENCE_METRICS")
    if PREFLIGHT_REPORT:
        pbmc_cache.write_preflight_report(
            PREFLIGHT_REPORT, res, reference_metrics, manifest_sha)
        log(f"preflight report saved to {PREFLIGHT_REPORT}")
    log("CPU preflight complete; exiting before scGPT model load")
    raise SystemExit(0)
Gco_sg, Gsg = build_pooled_scgpt_graphs(allcells, f"{OUT}/pbmc_scgpt_pooled_v2.npz")
res.update(run_test(a_p, Gco_sg[iiv, jjv], Gsg[iiv, jjv], conf, "scgpt"))
log("=== POOLED PBMC (paired-cell) ===")
for k, v in res.items():
    if isinstance(v, float): log(f"  {k}: {v}")

# ---- PER-TYPE ----
pertype = []
for t in ([] if POOLED_ONLY else types):
    cid = np.where(np.array(labels) == t)[0] if labels is not None else np.array([])
    if len(cid) < MIN_CELLS: log(f"  skip {t}: {len(cid)} cells"); continue
    Gc, _, _ = build_graphs(cid, f"{OUT}/pbmc_fmgraphs_{t}.npz")
    Gat = Z[f"G_{t}"]
    a_t = Gat[iiv, jjv]; co_t = Gc[iiv, jjv]
    coexp_vs_atac = round(float(spearmanr(co_t, a_t).statistic), 4)
    row = dict(cell_type=t, n=int(len(cid)), coexp_vs_atac=coexp_vs_atac)
    pertype.append(row)
    log(f"  {t}: n={len(cid)} coexp_vs_atac={coexp_vs_atac}")

res["per_type_coexp"] = pertype
json.dump(res, open(f"{OUT}/pbmc_eval_v2.json", "w"), indent=2)
log(f"SAVED {OUT}/pbmc_eval_v2.json")
