#!/usr/bin/env python
"""
scfm-reg-audit v2 — CONSTRUCT-VALID G_ATAC (fixes the load-bearing CRIT).

Replaces the pilot's |Pearson| gene-activity co-accessibility (a co-expression-family estimand)
with a SEQUENCE-GROUNDED, co-expression-orthogonal regulatory truth:

  G_ATAC_T[i, j] = Σ_p  a_T[p] · L[j, p] · HT[p, i]         (raw, log1p at readout — Spearman-invariant)

  i = TF gene (manifest gene with >=1 JASPAR motif), j = target gene, p = peak
  L[j,p] = 1 if peak p ∈ [TSS-2kb, gene-end] of gene j        (peak->gene link)
  HT[p,i] = 1 if any motif of TF i scores a MOODS hit (p<1e-4, either strand) in peak p's hg38 sequence
  a_T[p]  = mean peak count over cells of type T (log1p)       (cell-type accessibility gate)

Edges come from DNA motif presence in accessible peaks, NOT from cross-cell covariation ->
orthogonal to co-expression by construction. Asymmetric (TF rows only). Per brain cell type.
Env: SAMPLE_PEAKS (test), CELL_CAP (accessibility subsample), MOTIF_P.
"""
import os, sys, json, time, gzip, hashlib, numpy as np, anndata as ad, scipy.sparse as sp
import multiprocessing as mp
import pyfaidx, motif_utils as mu
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

PEAK_WIDTH = int(os.environ.get("PEAK_WIDTH", "500"))   # resize peaks to accessible core (motifmatchr convention)
NPROC = int(os.environ.get("NPROC", str(max(1, mp.cpu_count() - 4))))

# --- worker globals for multiprocessed motif scan ---
_W = {}
def _winit(meme, use_ids, p, hg38):
    mot = mu.parse_meme(meme); use = {m: mot[m] for m in use_ids}
    mats, thr, order = mu.build_scanner(use, p=p)
    _W.update(mats=mats, thr=thr, order=order, fa=pyfaidx.Fasta(hg38, sequence_always_upper=True))
def _wscan(arg):
    k, chrom, a, b = arg
    fa = _W["fa"]
    if chrom not in fa.keys(): return (k, ())
    seq = str(fa[chrom][max(0, a):b])
    if not seq: return (k, ())
    return (k, tuple(mu.scan_seq_hits(seq, _W["mats"], _W["thr"], _W["order"])))

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"
COORDS = f"{ROOT}/data/annotation/gene_coords_hg38.tsv"
MEME = f"{ROOT}/data/motifs/JASPAR2024_CORE_vertebrates.meme"
HG38 = f"{ROOT}/data/genome/hg38.fa"
ATAC = os.environ.get("ATAC_FILE", f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")
META = os.environ.get("META_FILE", os.path.join(ROOT, "../../research/sc-fm-benchmark/raw_pulls/scatac/atac_cell_meta.csv.gz"))
TAG = os.environ.get("TAG", "GSE174367")
OUT = f"{ROOT}/results/v2"; os.makedirs(OUT, exist_ok=True)
CACHE = f"{ROOT}/data/motifs"
PROMOTER = 2000; CELL_CAP = int(os.environ.get("CELL_CAP", "8000")); MIN_CELLS = 300
MOTIF_P = float(os.environ.get("MOTIF_P", "1e-5")); SAMPLE = int(os.environ.get("SAMPLE_PEAKS", "0"))

# ---- manifest ----
man = json.load(open(MANI)); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}
Ng = len(genes); tf_flags = man["tf_flags"]; TFgenes = [g for g in genes if tf_flags[g]]
log(f"manifest {Ng} genes ({len(TFgenes)} TFs) sha={man['sha256'][:12]}")

# ---- coords for manifest genes ----
gco = {}
for ln in open(COORDS):
    chrom, s, e, strand, name = ln.rstrip("\n").split("\t")
    if name not in gidx or name in gco: continue
    s, e = int(s), int(e)
    lo = s - PROMOTER if strand == "+" else s; hi = e if strand == "+" else e + PROMOTER
    gco[name] = (chrom, lo, hi)

# ---- ATAC peaks (coords only, backed) ----
Av = ad.read_h5ad(ATAC, backed="r"); peaks = [str(p) for p in Av.var_names]
pchr = np.array([p.split(":")[0] for p in peaks])
pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks])
pmid = (pse[:, 0] + pse[:, 1]) // 2
by_chr = {}
for i, c in enumerate(pchr): by_chr.setdefault(c, []).append(i)
for c in by_chr: by_chr[c] = np.array(by_chr[c])

# ---- peak -> gene links (manifest genes) ----
rows_g, cols_p = [], []
for g, (chrom, lo, hi) in gco.items():
    pis = by_chr.get(chrom)
    if pis is None: continue
    sel = pis[(pmid[pis] >= lo) & (pmid[pis] <= hi)]
    rows_g.extend([gidx[g]] * len(sel)); cols_p.extend(sel.tolist())
rel = np.array(sorted(set(cols_p))); reidx = {p: k for k, p in enumerate(rel)}
if SAMPLE: rel = rel[:SAMPLE]; reidx = {p: k for k, p in enumerate(rel)}
nRel = len(rel)
keep = [(g, p) for g, p in zip(rows_g, cols_p) if p in reidx]
L = sp.csr_matrix((np.ones(len(keep), np.float32), ([g for g, _ in keep], [reidx[p] for _, p in keep])),
                  shape=(Ng, nRel))
log(f"peak->gene links {len(keep)} | relevant peaks {nRel}")

# ---- motifs restricted to manifest TFs; TF -> motif_ids ----
mot = mu.parse_meme(MEME)
tf2mot = {}; use_mot = {}
for mid, d in mot.items():
    for tf in mu.tf_symbols(d["name"]):
        if tf in gidx and tf_flags[tf]:
            tf2mot.setdefault(tf, set()).add(mid); use_mot[mid] = d
log(f"manifest TFs with >=1 motif: {len(tf2mot)} | motifs used: {len(use_mot)}")

# ---- scan relevant peaks (resized to PEAK_WIDTH core) -> motif hits (cached) ----
mot_order = sorted(use_mot); mcol = {m: k for k, m in enumerate(mot_order)}
ck = hashlib.sha256(f"{man['sha256']}|{nRel}|{MOTIF_P}|{PEAK_WIDTH}|{mot_order}".encode()).hexdigest()[:16]
cpath = f"{CACHE}/peak_motif_{ck}.npz"
if os.path.exists(cpath):
    Z = np.load(cpath, allow_pickle=False); H = sp.csr_matrix((Z["Hdata"], Z["Hind"], Z["Hptr"]), shape=tuple(Z["Hshape"])); mot_order = list(Z["mot_order"])
    log(f"loaded cached peak×motif {H.shape} from {os.path.basename(cpath)}")
else:
    w2 = PEAK_WIDTH // 2
    args = [(k, pchr[p], int(pmid[p]) - w2, int(pmid[p]) + w2) for k, p in enumerate(rel)]
    hr, hc = [], []; t0 = time.time()
    with mp.Pool(NPROC, initializer=_winit, initargs=(MEME, mot_order, MOTIF_P, HG38)) as pool:
        for k, hits in pool.imap_unordered(_wscan, args, chunksize=64):
            for mid in hits: hr.append(k); hc.append(mcol[mid])
    H = sp.csr_matrix((np.ones(len(hr), np.float32), (hr, hc)), shape=(nRel, len(mot_order)))
    np.savez(cpath, Hdata=H.data, Hind=H.indices, Hptr=H.indptr, Hshape=H.shape, mot_order=np.array(mot_order))
    log(f"scanned {nRel} peaks @ {PEAK_WIDTH}bp p={MOTIF_P} with {NPROC} procs -> {H.nnz} hits "
        f"({H.nnz/max(1,nRel):.1f}/peak, {time.time()-t0:.0f}s), cached")

# ---- peak × TF (OR over a TF's motifs) ----
mcol = {m: k for k, m in enumerate(mot_order)}
htr, htc = [], []
for tf, mids in tf2mot.items():
    cols = [mcol[m] for m in mids if m in mcol]
    if not cols: continue
    hit_peaks = H[:, cols].getnnz(axis=1).nonzero()[0]
    htr.extend(hit_peaks.tolist()); htc.extend([gidx[tf]] * len(hit_peaks))
HT = sp.csr_matrix((np.ones(len(htr), np.float32), (htr, htc)), shape=(nRel, Ng))
log(f"peak×TF matrix {HT.shape} nnz={HT.nnz}")

# ---- per-cell-type (or pooled) accessibility over relevant peaks ----
A = ad.read_h5ad(ATAC)
bc = np.array([str(b) for b in A.obs_names])
if META not in ("none", "") and os.path.exists(META):
    meta = {}
    with gzip.open(META, "rt") as f:
        hdr = f.readline().rstrip("\n").split(","); bi = hdr.index("Barcode"); ci = hdr.index("Cell.Type")
        for ln in f:
            q = ln.rstrip("\n").split(","); meta[q[bi]] = q[ci]
    lab = np.array([meta.get(b, "NA") for b in bc])
else:
    lab = np.array(["ALL"] * len(bc)); log("no META -> POOLED single accessibility (type=ALL)")
X = (A.X.tocsc() if sp.issparse(A.X) else sp.csr_matrix(A.X).tocsc())[:, rel].tocsr().astype(np.float32)
from collections import Counter
rng = np.random.default_rng(20260713); vc = Counter(lab.tolist())
types = [t for t, c in vc.most_common() if c >= MIN_CELLS and t != "NA"]
log("cell types:", [(t, vc[t]) for t in types])

graphs = {}
for t in types:
    cid = np.where(lab == t)[0]
    if len(cid) > CELL_CAP: cid = rng.choice(cid, size=CELL_CAP, replace=False)
    aT = np.asarray(X[cid].mean(0)).ravel(); aT = np.log1p(aT).astype(np.float32)   # (nRel,)
    Aw = L.multiply(aT[None, :]).tocsr()                    # genes × peaks, accessibility-weighted
    G = (Aw @ HT).T.toarray().astype(np.float32)            # (Ng TF-rows × Ng targets)
    np.fill_diagonal(G, 0.0)
    graphs[t] = G
    log(f"  G_ATAC[{t}] n={len(cid)} nnz_rows={int((G.sum(1)>0).sum())}")

tf_rows = np.array([gidx[g] for g in TFgenes])
np.savez(f"{OUT}/G_ATAC_v2_{TAG}.npz",
         genes=np.array(genes), types=np.array(types), tf_rows=tf_rows,
         **{f"G_{t}": graphs[t] for t in types})
json.dump(dict(construct="motif->accessible-peak->TF->target", tag=TAG, atac=os.path.basename(ATAC),
               n_genes=Ng, n_tf=len(TFgenes), relevant_peaks=int(nRel), peak_motif_hits=int(H.nnz),
               motif_p=MOTIF_P, manifest_sha=man["sha256"], types={t: int(vc[t]) for t in types}),
          open(f"{OUT}/G_ATAC_v2_{TAG}_meta.json", "w"), indent=2)
log(f"SAVED {OUT}/G_ATAC_v2_{TAG}.npz ({len(types)} types, {Ng} genes, {len(TFgenes)} TF rows)")
