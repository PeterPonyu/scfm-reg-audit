#!/usr/bin/env python
"""
scfm-reg-audit v2 — IN-SILICO KO readout (the 3rd/last FM readout; could overturn the null).

Geneformer in-silico perturbation: delete a TF's token from each cell, re-embed, and measure how
much each other gene's contextual embedding shifts. influence(TF g -> target j) = mean over cells of
(1 - cos(e_j^baseline, e_j^{-g})). This builds a TF->target FM INFLUENCE graph — the causal/attention-
free readout — and tests it vs the construct-valid G_ATAC v2, marginal + partial | co-expression.

Position-shift control: also delete a RANDOM non-TF gene per cell and subtract its mean influence
(nets out the rank-shift artifact of token deletion). GPU. N_CELLS cells (default 300), all TFs.
"""
import os, json, hashlib, time, numpy as np, torch
from scipy.stats import spearmanr, rankdata
import fm_readout as fr
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
dev = "cuda" if torch.cuda.is_available() else "cpu"

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
GFM = f"{DATA_ROOT}/models/Geneformer/Geneformer-V2-104M"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
N_CELLS = int(os.environ.get("N_CELLS", "300")); BATCH = int(os.environ.get("BATCH", "24"))

man = json.load(open(MANI)); genes = man["genes"]; gidx = {g: i for i, g in enumerate(genes)}; Ng = len(genes)
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
G_co = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")["co"]

A, Xc, Xl, _ = fr.load_norm(RNA)
rsym = {str(s): k for k, s in enumerate(A.var_names)}; ri = np.array([rsym[g] for g in genes])
Xc_g = Xc[:, ri].tocsr()
rng = np.random.default_rng(20260713)
cells = rng.choice(Xc_g.shape[0], size=min(N_CELLS, Xc_g.shape[0]), replace=False)
rd = fr.FMReadout(genes, batch=BATCH); gf_tokens, gf_median = rd.gf_tokens, rd.gf_median
tok2gene = {int(gf_tokens[i]): i for i in range(Ng)}
tf_set = set(tf_rows.tolist())

from transformers.models.bert.modeling_bert import BertModel
mdl = BertModel.from_pretrained(GFM, output_hidden_states=True, add_pooling_layer=False).to(dev).half().eval()
CLS, EOS, PAD = rd.gtok["<cls>"], rd.gtok["<eos>"], rd.gtok["<pad>"]

def embed_seqs(seqs):
    """seqs: list of token lists (already CLS..EOS). Returns list of per-position hidden vectors."""
    L = max(len(s) for s in seqs); ic = np.full((len(seqs), L), PAD, np.int64); am = np.zeros((len(seqs), L), np.int64)
    for i, s in enumerate(seqs): ic[i, :len(s)] = s; am[i, :len(s)] = 1
    with torch.no_grad():
        h = mdl(input_ids=torch.tensor(ic, device=dev), attention_mask=torch.tensor(am, device=dev)).hidden_states[-1].float().cpu().numpy()
    return h

# ---- per-cell baseline ranked token order + baseline gene embeddings ----
cell_order = {}; base_emb = {}
for s in range(0, len(cells), BATCH):
    batch = cells[s:s + BATCH]; seqs = []
    for r in batch:
        x = Xc_g[r].toarray().ravel(); val = x / gf_median; nz = np.where(val > 0)[0]; order = nz[np.argsort(-val[nz])]
        cell_order[r] = order; seqs.append([CLS] + gf_tokens[order].tolist() + [EOS])
    h = embed_seqs(seqs)
    for i, r in enumerate(batch):
        k = len(cell_order[r]); base_emb[r] = h[i, 1:1 + k].copy()      # (k, D) for present genes in rank order
log(f"baseline embeddings for {len(cells)} cells")

# ---- KO each TF; plus a random-gene deletion control ----
infl = np.zeros((Ng, Ng), np.float64); cnt = np.zeros((Ng, Ng), np.float64)
ctrl_infl = np.zeros(Ng, np.float64); ctrl_cnt = np.zeros(Ng, np.float64)   # per-target mean control shift
t0 = time.time()
for ti, g in enumerate(tf_rows):
    cells_with = [r for r in cells if g in set(cell_order[r].tolist())]
    if not cells_with: continue
    for s in range(0, len(cells_with), BATCH):
        batch = cells_with[s:s + BATCH]; seqs = []; metas = []
        for r in batch:
            order = cell_order[r]; pos = int(np.where(order == g)[0][0])
            neworder = np.delete(order, pos)
            seqs.append([CLS] + gf_tokens[neworder].tolist() + [EOS]); metas.append((r, neworder))
        h = embed_seqs(seqs)
        for i, (r, neworder) in enumerate(metas):
            pe = h[i, 1:1 + len(neworder)]                              # perturbed emb, aligned to neworder
            be = base_emb[r]; bpos = {int(gg): p for p, gg in enumerate(cell_order[r])}
            for p2, gg in enumerate(neworder):
                b = be[bpos[int(gg)]]; a = pe[p2]
                sh = 1.0 - float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9))
                infl[g, gg] += sh; cnt[g, gg] += 1
    if (ti + 1) % 50 == 0: log(f"  KO {ti+1}/{len(tf_rows)} TFs ({time.time()-t0:.0f}s)")

# ---- random-gene deletion control (one random present gene per cell) ----
for s in range(0, len(cells), BATCH):
    batch = cells[s:s + BATCH]; seqs = []; metas = []
    for r in batch:
        order = cell_order[r]
        if len(order) < 3: continue
        pos = int(rng.integers(len(order))); neworder = np.delete(order, pos)
        seqs.append([CLS] + gf_tokens[neworder].tolist() + [EOS]); metas.append((r, neworder))
    if not seqs: continue
    h = embed_seqs(seqs)
    for i, (r, neworder) in enumerate(metas):
        pe = h[i, 1:1 + len(neworder)]; be = base_emb[r]; bpos = {int(gg): p for p, gg in enumerate(cell_order[r])}
        for p2, gg in enumerate(neworder):
            b = be[bpos[int(gg)]]; a = pe[p2]
            sh = 1.0 - float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9))
            ctrl_infl[gg] += sh; ctrl_cnt[gg] += 1
ctrl = ctrl_infl / np.maximum(ctrl_cnt, 1)                              # per-target baseline shift

cnt[cnt == 0] = 1; G_ko = (infl / cnt).astype(np.float32)
G_ko_ctrl = np.clip(G_ko - ctrl[None, :], 0, None).astype(np.float32)   # net of position-shift artifact
del mdl; torch.cuda.empty_cache()

ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows)); m = ii != jj; ii, jj = ii[m], jj[m]
a = G_atac[ii, jj]; co = G_co[ii, jj]
def partial(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def rsd(v, c): c1 = np.c_[np.ones_like(c), c]; b = np.linalg.lstsq(c1, v, rcond=None)[0]; return v - c1 @ b
    return float(np.corrcoef(rsd(rx, rz), rsd(ry, rz))[0, 1])
np.savez(f"{OUT}/G_ko_v2.npz", G_ko=G_ko, G_ko_ctrl=G_ko_ctrl, tf_rows=tf_rows)

# ---- Mantel gene-label permutation null for the KO partials (the one positive readout) ----
def mantel_partial(Gpred, nperm=1000):
    x = Gpred[ii, jj]; obs = partial(x, a, co); null = []
    for _ in range(nperm):
        p = rng.permutation(Ng); null.append(partial(x, G_atac[p[ii], p[jj]], co))
    nd = np.array(null)
    return dict(observed=round(obs, 4), null_mean=round(float(nd.mean()), 4), null_sd=round(float(nd.std()), 4),
                z=round(float((obs - nd.mean()) / (nd.std() + 1e-9)), 2),
                p_perm=round(float((np.sum(np.abs(nd) >= abs(obs)) + 1) / (nperm + 1)), 4))
log("Mantel null for KO partials…")
ko_mantel = mantel_partial(G_ko); ko_ctrl_mantel = mantel_partial(G_ko_ctrl)

res = dict(readout="geneformer_insilico_KO", n_cells=int(len(cells)), n_tf=int(len(tf_rows)), n_pairs=int(len(ii)),
           ko_vs_atac=round(float(spearmanr(G_ko[ii, jj], a).statistic), 4),
           ko_partial_given_coexp=round(partial(G_ko[ii, jj], a, co), 4),
           ko_vs_coexp=round(float(spearmanr(G_ko[ii, jj], co).statistic), 4),
           ko_ctrl_vs_atac=round(float(spearmanr(G_ko_ctrl[ii, jj], a).statistic), 4),
           ko_ctrl_partial_given_coexp=round(partial(G_ko_ctrl[ii, jj], a, co), 4),
           coexp_vs_atac=round(float(spearmanr(co, a).statistic), 4),
           mantel_ko_partial=ko_mantel, mantel_ko_ctrl_partial=ko_ctrl_mantel)
json.dump(res, open(f"{OUT}/insilico_ko_v2.json", "w"), indent=2)
log("=== IN-SILICO KO READOUT (vs G_ATAC v2) ===")
for k, v in res.items():
    if isinstance(v, float): log(f"  {k}: {v}")
log(f"SAVED {OUT}/insilico_ko_v2.json")
