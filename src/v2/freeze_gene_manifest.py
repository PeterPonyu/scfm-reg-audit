#!/usr/bin/env python
"""
scfm-reg-audit v2 — freeze the pre-registered shared gene manifest (fixes pilot `shared[:1000]` bias).

Universe (deterministic, order = sorted symbol):
  hg38-coord protein-coding  ∩  Geneformer-tokenizable  ∩  scGPT-vocab
  ∩  >=1 ATAC peak in [TSS-2kb, gene end]  ∩  RNA-detected in >= DETECT_FRAC of ad_hm cells.
Cap to N_MAX (fits scGPT MAXLEN) by a PRE-REGISTERED rule: keep all TFs (JASPAR symbols) that pass,
then fill remaining slots by RNA detection rate (modality-neutral, NOT the ATAC test signal).
Writes data/manifest/shared_genes.v2.json with a sha256 of the sorted gene list.
"""
import os, json, hashlib, pickle, time, numpy as np, anndata as ad, scipy.sparse as sp
import motif_utils as mu
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"
MEME = f"{ROOT}/data/motifs/JASPAR2024_CORE_vertebrates.meme"
ATAC = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
GF = f"{DATA_ROOT}/models/Geneformer/geneformer"; CK = f"{DATA_ROOT}/models/scGPT-human"
OUTDIR = f"{ROOT}/data/manifest"; os.makedirs(OUTDIR, exist_ok=True)
PROMOTER = 2000; DETECT_FRAC = float(os.environ.get("DETECT_FRAC", "0.01")); N_MAX = int(os.environ.get("N_MAX", "1200"))

# ---- tokenizable vocabularies ----
gtok = pickle.load(open(f"{GF}/token_dictionary_gc104M.pkl", "rb"))
gmed = pickle.load(open(f"{GF}/gene_median_dictionary_gc104M.pkl", "rb"))
gn2i = pickle.load(open(f"{GF}/gene_name_id_dict_gc104M.pkl", "rb"))
svoc = json.load(open(f"{CK}/vocab.json"))
def gf_ok(s): e = gn2i.get(s); return e is not None and e in gtok and e in gmed

# ---- TF symbols from JASPAR ----
mot = mu.parse_meme(MEME); TF = set()
for d in mot.values():
    for t in mu.tf_symbols(d["name"]): TF.add(t)
log("JASPAR TF symbols:", len(TF))

# ---- gene coords (first occurrence per symbol) ----
genes = {}
for ln in open(COORDS):
    chrom, s, e, strand, name = ln.rstrip("\n").split("\t")
    if name in genes: continue
    s, e = int(s), int(e)
    lo = s - PROMOTER if strand == "+" else s
    hi = e if strand == "+" else e + PROMOTER
    genes[name] = (chrom, lo, hi)
log("coord genes:", len(genes))

# ---- ATAC peaks -> which genes have >=1 peak ----
Av = ad.read_h5ad(ATAC, backed="r")
pchr = np.array([str(p).split(":")[0] for p in Av.var_names])
pse = np.array([[int(x) for x in str(p).split(":")[1].split("-")] for p in Av.var_names])
pmid = (pse[:, 0] + pse[:, 1]) // 2
by_chr = {}
for i, c in enumerate(pchr): by_chr.setdefault(c, []).append(i)
for c in by_chr: by_chr[c] = np.array(by_chr[c])
has_peak = set()
for name, (chrom, lo, hi) in genes.items():
    pis = by_chr.get(chrom)
    if pis is None: continue
    m = pmid[pis]
    if np.any((m >= lo) & (m <= hi)): has_peak.add(name)
log("genes with >=1 ATAC peak:", len(has_peak))

# ---- RNA detection rate (ad_hm) ----
R = ad.read_h5ad(RNA)
X = R.X.tocsc() if sp.issparse(R.X) else sp.csr_matrix(R.X).tocsc()
rsym = {str(s): j for j, s in enumerate(R.var_names)}
ncell = X.shape[0]
det = {}
for name, j in rsym.items():
    col = X[:, j]
    det[name] = float(col.getnnz()) / ncell
log("RNA genes:", len(rsym), "| ncells:", ncell)

# ---- assemble universe ----
uni = [g for g in genes
       if gf_ok(g) and g in svoc and g in has_peak and det.get(g, 0.0) >= DETECT_FRAC]
log("universe (all filters):", len(uni))

# ---- pre-registered cap: all TFs first, then fill by detection rate ----
uni_tf = sorted([g for g in uni if g in TF])
uni_non = sorted([g for g in uni if g not in TF], key=lambda g: -det[g])
if len(uni_tf) >= N_MAX:
    chosen = sorted(uni_tf, key=lambda g: -det[g])[:N_MAX]
else:
    chosen = uni_tf + uni_non[:max(0, N_MAX - len(uni_tf))]
chosen = sorted(set(chosen))                     # deterministic final order
n_tf = sum(g in TF for g in chosen)
log(f"manifest: {len(chosen)} genes ({n_tf} TFs)")

sha = hashlib.sha256(("\n".join(chosen)).encode()).hexdigest()
out = dict(
    n_genes=len(chosen), n_tf=n_tf, sha256=sha,
    params=dict(promoter_bp=PROMOTER, detect_frac=DETECT_FRAC, n_max=N_MAX,
                jaspar="JASPAR2024_CORE_vertebrates", coords="gencode.v44.basic",
                rna="ad_hm_prepped", atac="GSE174367"),
    cap_rule="keep all TFs passing filters; fill remainder by RNA detection rate (desc)",
    genes=chosen,
    tf_flags={g: (g in TF) for g in chosen},
    detection={g: round(det[g], 4) for g in chosen},
)
json.dump(out, open(f"{OUTDIR}/shared_genes.v2.json", "w"), indent=2)
log(f"SAVED {OUTDIR}/shared_genes.v2.json  sha256={sha[:16]}…")
