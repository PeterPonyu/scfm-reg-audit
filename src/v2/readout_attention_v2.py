#!/usr/bin/env python
"""
scfm-reg-audit v2 — ATTENTION readout vs construct-valid G_ATAC v2.

The embedding-cosine readout is degenerate here (geneformer_vs_coexp=-0.18, a Geneformer anisotropy
artifact), so it can't ground an FM claim. Attention is the interpretability readout both Kendiukhov
scoops used. This recomputes G_FM as Geneformer last-layer gene-gene attention (mean over heads,
pooled), with CORRECTED non-log CP10k ranking + the frozen manifest, and reruns the decisive test
on the regulatory pair set P (TF->target). Reports symmetric and directed (TF->target, target->TF).
"""
import os, json, hashlib, time, numpy as np, torch
from scipy.stats import spearmanr, rankdata
import fm_readout as fr
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
dev = "cuda" if torch.cuda.is_available() else "cpu"; BATCH = int(os.environ.get("BATCH", "8"))
CAP = int(os.environ.get("CELL_CAP", "4000")); NPERM = int(os.environ.get("NPERM", "500"))

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"; OUT = f"{ROOT}/results/v2"
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
GFM = f"{DATA_ROOT}/models/Geneformer/Geneformer-V2-104M"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"

man = json.load(open(MANI)); genes = man["genes"]
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
Ng = len(genes)
Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
types = [str(t) for t in Z["types"]]; tf_rows = np.array(Z["tf_rows"])
G_atac = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
G_co = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")["co"]     # reuse cached co-expression
ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows))
m = ii != jj; ii, jj = ii[m], jj[m]
a = G_atac[ii, jj]; co = G_co[ii, jj]
log(f"manifest OK, P={len(ii)} TF->target pairs, {len(tf_rows)} TFs")

# ---- RNA -> manifest genes (non-log CP10k for GF ranking) ----
A, Xc, Xl, _ = fr.load_norm(RNA)
rsym = {str(s): k for k, s in enumerate(A.var_names)}; ri = np.array([rsym[g] for g in genes])
Xc_g = Xc[:, ri].tocsr()
rng = np.random.default_rng(20260713)
cells = np.arange(Xc_g.shape[0])
if len(cells) > CAP: cells = rng.choice(cells, size=CAP, replace=False)
log(f"attention over {len(cells)} cells")

rd = fr.FMReadout(genes, batch=BATCH)
gf_tokens, gf_median = rd.gf_tokens, rd.gf_median
from transformers.models.bert.modeling_bert import BertModel
mdl = BertModel.from_pretrained(GFM, output_attentions=True, add_pooling_layer=False)
mdl = mdl.to(dev).half().eval() if dev == "cuda" else mdl.eval()
CLS, EOS, PAD = rd.gtok["<cls>"], rd.gtok["<eos>"], rd.gtok["<pad>"]
Asum = np.zeros((Ng, Ng), np.float64); Acnt = np.zeros((Ng, Ng), np.float64)
t0 = time.time()
for s in range(0, len(cells), BATCH):
    rows, orders, lens = [], [], []
    for r in cells[s:s + BATCH]:
        x = Xc_g[r].toarray().ravel(); val = x / gf_median      # non-log CP10k / median (corrected)
        nz = np.where(val > 0)[0]; order = nz[np.argsort(-val[nz])]
        rows.append([CLS] + gf_tokens[order].tolist() + [EOS]); orders.append(order); lens.append(len(order) + 2)
    L = max(lens); ic = np.full((len(rows), L), PAD, np.int64); am = np.zeros((len(rows), L), np.int64)
    for i, r in enumerate(rows): ic[i, :len(r)] = r; am[i, :len(r)] = 1
    with torch.no_grad():
        att = mdl(input_ids=torch.tensor(ic, device=dev), attention_mask=torch.tensor(am, device=dev)).attentions[-1]
        att = att.float().mean(1).cpu().numpy()                  # (b, L, L) mean over heads
    for i, order in enumerate(orders):
        k = len(order); sub = att[i, 1:1 + k, 1:1 + k]
        Asum[np.ix_(order, order)] += sub; Acnt[np.ix_(order, order)] += 1
    if s % (BATCH * 100) == 0: log(f"  {s}/{len(cells)} ({time.time()-t0:.0f}s)")
Acnt[Acnt == 0] = 1; Aatt = (Asum / Acnt).astype(np.float32)     # directional gene-gene attention
G_sym = np.abs((Aatt + Aatt.T) / 2)

def partial(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def rsd(v, c): c1 = np.c_[np.ones_like(c), c]; b = np.linalg.lstsq(c1, v, rcond=None)[0]; return v - c1 @ b
    return float(np.corrcoef(rsd(rx, rz), rsd(ry, rz))[0, 1])

at_sym = G_sym[ii, jj]; at_dir = np.abs(Aatt[ii, jj]); at_rev = np.abs(Aatt[jj, ii])
obs = dict(
    attn_sym_vs_atac=float(spearmanr(at_sym, a).statistic),
    attn_sym_partial_given_coexp=partial(at_sym, a, co),
    attn_sym_vs_coexp=float(spearmanr(at_sym, co).statistic),
    attn_tf2target_vs_atac=float(spearmanr(at_dir, a).statistic),
    attn_tf2target_partial_given_coexp=partial(at_dir, a, co),
    attn_target2tf_vs_atac=float(spearmanr(at_rev, a).statistic),
    coexp_vs_atac=float(spearmanr(co, a).statistic),
)
# Mantel for the primary (symmetric attention partial | coexp)
rng2 = np.random.default_rng(7)
null = []
for _ in range(NPERM):
    perm = rng2.permutation(Ng); null.append(partial(at_sym, G_atac[perm[ii], perm[jj]], co))
nd = np.array(null); o = obs["attn_sym_partial_given_coexp"]
mantel = dict(observed=round(o, 4), null_mean=round(float(nd.mean()), 4), null_sd=round(float(nd.std()), 4),
              z=round(float((o - nd.mean()) / (nd.std() + 1e-9)), 2),
              p_perm=round(float((np.sum(np.abs(nd) >= abs(o)) + 1) / (NPERM + 1)), 4))
res = dict(readout="geneformer_attention_lastlayer", n_cells=int(len(cells)), n_pairs=int(len(ii)),
           observed={k: round(v, 4) for k, v in obs.items()}, mantel_attn_sym_partial=mantel, manifest_sha=man["sha256"])
json.dump(res, open(f"{OUT}/readout_attention_v2.json", "w"), indent=2)
log("=== ATTENTION READOUT v2 (vs construct-valid G_ATAC) ===")
for k, v in obs.items(): log(f"  {k}: {round(v,4)}")
log(f"  MANTEL attn_sym_partial: z={mantel['z']} p={mantel['p_perm']}")
log(f"SAVED {OUT}/readout_attention_v2.json")
