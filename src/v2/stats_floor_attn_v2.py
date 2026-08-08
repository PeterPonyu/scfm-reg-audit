#!/usr/bin/env python
"""
Confound-controlled partial rho + TF-block bootstrap CI + Mantel-on-adjusted for the two graphs from
gen_floor_attn_graphs_v2.py, using the IDENTICAL machinery as the pooled main table (imports
stats_enhanced_v2.analyze) so the rows are directly comparable:
  - geneformer_randinit_floor : random-init Geneformer embedding readout (the "untrained" floor)
  - geneformer_attn (brain)   : the previously-missing brain attention pooled row
"""
import os
import json, numpy as np
import stats_enhanced_v2 as se
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))

OUT = se.OUT
ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"

Zb = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
types_b = [str(t) for t in Zb["types"]]; tf_b = np.array(Zb["tf_rows"])
Gb = np.mean([Zb[f"G_{t}"] for t in types_b], axis=0).astype(np.float32)
Fb = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")

G_floor = np.load(f"{OUT}/brain_floor_graph_v2.npz")["G"]
G_attn = np.load(f"{OUT}/brain_attention_graph_v2.npz")["G_sym"]

results = []
se.log("=== random-init Geneformer floor (brain, embedding readout) ===")
results.append(se.analyze("brain", tf_b, Gb, Fb["co"], G_floor, ATAC_B, "geneformer_randinit_floor"))
se.log("=== brain attention (previously-missing pooled row) ===")
results.append(se.analyze("brain", tf_b, Gb, Fb["co"], G_attn, ATAC_B, "geneformer_attn"))

json.dump(results, open(f"{OUT}/stats_floor_attn_v2.json", "w"), indent=2)
se.log("=== SUMMARY ===")
for r in results:
    m = r["mantel_confound_controlled"]
    se.log(f"  {r['label']:28s} obs={r['observed']:+.4f} CI95=[{r['bootstrap_ci95'][0]:+.4f},{r['bootstrap_ci95'][1]:+.4f}] "
           f"z={m['z']:+.2f} p={m['p_perm']:.4f}")
se.log(f"SAVED {OUT}/stats_floor_attn_v2.json")

# compare floor vs trained embedding (from stats_enhanced_v2.json) for the headline message
enh = {r["label"]: r for r in json.load(open(f"{OUT}/stats_enhanced_v2.json")) if r["tag"] == "brain"}
tr = enh.get("geneformer_embed")
fl = results[0]
if tr:
    se.log(f"FLOOR CHECK: trained geneformer_embed obs={tr['observed']:+.4f} CI={tr['bootstrap_ci95']}  vs  "
           f"random-init floor obs={fl['observed']:+.4f} CI={fl['bootstrap_ci95']}")
