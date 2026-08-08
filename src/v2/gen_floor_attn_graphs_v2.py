#!/usr/bin/env python
"""
scfm-reg-audit v2 -- generate two brain gene-graphs that complete the pooled main table (zero net,
local reuse only):
  (B) random-init Geneformer floor  -- same architecture, SEEDED random weights, embedding-cosine
      readout. The "even an untrained model scores like this" anchor; pairs with the positive control.
  (C) brain last-layer attention    -- was computed inline in readout_attention_v2.py and never
      cached, so it was the one pooled row missing from stats_enhanced_v2 / Fig 2. Cache it here.
Both use the SAME brain corpus, manifest, and cell sample as the trained-Geneformer pooled readout,
so they are directly comparable. Graphs saved to results/v2/; run stats separately (stats_floor_attn_v2.py).
"""
import os, json, hashlib, time, numpy as np, torch
import fm_readout as fr
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUT = f"{ROOT}/results/v2"
MANI = f"{ROOT}/data/manifest/shared_genes.v2.json"
RNA = f"{DATA_ROOT}/datasets/extra_preprocessed/ad_hm_prepped.h5ad"
BATCH = int(os.environ.get("BATCH", "16")); CAP = int(os.environ.get("CELL_CAP", "4000"))
SEED = int(os.environ.get("SEED", "0"))

man = json.load(open(MANI)); genes = man["genes"]
assert hashlib.sha256(("\n".join(genes)).encode()).hexdigest() == man["sha256"]
Ng = len(genes)

A, Xc, Xl, _ = fr.load_norm(RNA)
rsym = {str(s): k for k, s in enumerate(A.var_names)}; ri = np.array([rsym[g] for g in genes])
Xc_g = Xc[:, ri].tocsr()
rng = np.random.default_rng(20260713)
cells = np.arange(Xc_g.shape[0])
if len(cells) > CAP: cells = rng.choice(cells, size=CAP, replace=False)
log(f"brain corpus: {len(cells)} cells, {Ng} genes, batch={BATCH}, dev={fr.DEV}")

rd = fr.FMReadout(genes, batch=BATCH)

# (B) random-init Geneformer floor (embedding-cosine readout, same as trained geneformer_embed) -----
t0 = time.time()
E_rand = rd.geneformer(Xc_g, cells, random_init=True, seed=SEED)
G_rand = fr.FMReadout.cos_graph(E_rand)
np.savez(f"{OUT}/brain_floor_graph_v2.npz", G=G_rand.astype(np.float32), seed=SEED, n_cells=len(cells))
log(f"(B) random-init floor graph done ({time.time()-t0:.0f}s) -> brain_floor_graph_v2.npz")

# (C) brain last-layer attention (symmetrized) -----------------------------------------------------
t0 = time.time()
G_attn = rd.geneformer_attention(Xc_g, cells)
np.savez(f"{OUT}/brain_attention_graph_v2.npz", G_sym=G_attn.astype(np.float32), n_cells=len(cells))
log(f"(C) brain attention graph done ({time.time()-t0:.0f}s) -> brain_attention_graph_v2.npz")
log("DONE gen_floor_attn_graphs_v2")
