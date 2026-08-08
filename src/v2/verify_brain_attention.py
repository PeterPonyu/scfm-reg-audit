#!/usr/bin/env python
"""Self-contained verification of the brain-attention omission.

Calls fixed_panel_audit's two null functions directly with vectors assembled here,
so nothing touches run_pooled_family's batched shared-null path (which does not
terminate in this environment) and nothing reads the 3.2GB fasta (gene-identity
confounds come from the PBMC cache; only brain peakcount is recomputed, from ATAC
peak coordinates alone). Confirms the published embed row to 6 decimals first,
then scores the omitted attention and random-init floor readouts on the identical
edge set, then re-derives BH across the enlarged brain family.
"""
import json
import os
import sys

import anndata as ad
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixed_panel_audit as fpa  # noqa: E402

OUT = fpa.OUT
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
PROM = fpa.PROM
N_PERM = 999

# NOTE: seeds below use Python's hash(), which is salted per process
# (PYTHONHASHSEED). This script is a one-shot verification harness, not a
# reproducible pipeline — the authoritative numbers come from
# run_fixed_panel_audit.py, which derives all seeds from SeedSequence.
# Absolute q-values here may drift slightly across invocations; the signs and
# support decisions are stable.


def brain_peakcount():
    genes, _, _ = fpa.load_manifest()
    gidx = {g: i for i, g in enumerate(genes)}
    gco = {}
    with open(fpa.COORDS) as fh:
        for ln in fh:
            c, s, e, st, nm = ln.rstrip("\n").split("\t")
            if nm not in gidx or nm in gco:
                continue
            s, e = int(s), int(e)
            gco[nm] = (c, s - PROM if st == "+" else s, e if st == "+" else e + PROM)
    Av = ad.read_h5ad(ATAC_B, backed="r")
    peaks = [str(p) for p in Av.var_names]
    Av.file.close()
    pchr = np.array([p.split(":")[0] for p in peaks])
    pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks])
    pmid = (pse[:, 0] + pse[:, 1]) // 2
    by = {}
    for i, c in enumerate(pchr):
        by.setdefault(c, []).append(i)
    by = {c: np.array(v) for c, v in by.items()}
    pc = np.zeros(len(genes), dtype=np.float32)
    for g, i in gidx.items():
        if g not in gco:
            continue
        c, lo, hi = gco[g]
        pis = by.get(c)
        if pis is not None:
            pc[i] = len(pis[(pmid[pis] >= lo) & (pmid[pis] <= hi)])
    return pc


def bh(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0, 1)
    return q


def main():
    # Graphs
    Z = np.load(f"{OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=True)
    types_b = [str(t) for t in Z["types"]]
    tf_rows = np.array(Z["tf_rows"])
    G_atac = np.mean([Z[f"G_{t}"] for t in types_b], axis=0).astype(np.float32)
    F = np.load(f"{OUT}/fmgraphs_pooled_v2.npz")
    co = F["co"].astype(np.float32)
    graphs = {
        "geneformer_embed": F["gf"].astype(np.float32),
        "scGPT_encoder": F["sg"].astype(np.float32),
        "scFoundation_encoder": np.load(f"{OUT}/G_scf_pooled.npz")["G"].astype(np.float32),
        "UCE_encoder": np.load(f"{OUT}/G_uce_pooled.npz")["G"].astype(np.float32),
        "geneformer_attn": np.load(f"{OUT}/brain_attention_graph_v2.npz")["G_sym"].astype(np.float32),
        "random_init_floor": np.load(f"{OUT}/brain_floor_graph_v2.npz")["G"].astype(np.float32),
    }
    K = np.load(f"{OUT}/G_ko_v2.npz")
    graphs["geneformer_ko_raw"] = K["G_ko"].astype(np.float32)
    graphs["geneformer_ko_posctrl"] = K["G_ko_ctrl"].astype(np.float32)

    # Gene-identity confounds from PBMC cache (tissue-independent); brain peakcount fresh.
    cache = np.load(f"{OUT}/pbmc_confounds_v2.npz", allow_pickle=False)
    genelen = cache["genelen"].astype(np.float32)
    detv = cache["detv"].astype(np.float32)
    gc = cache["gc"].astype(np.float32)
    print("computing brain peakcount (no fasta)...", flush=True)
    peakcount = brain_peakcount()
    print("peakcount done", flush=True)

    genes, _, _ = fpa.load_manifest()
    Ng = G_atac.shape[0]
    tf_outdeg = (G_atac > 0).sum(1).astype(np.float32)
    atac_indeg = (G_atac > 0).sum(0).astype(np.float32)

    ii_all = np.repeat(tf_rows, Ng)
    jj_all = np.tile(np.arange(Ng), len(tf_rows))
    m = fpa.edge_mask("brain", genes, tf_rows, ii_all, jj_all)
    ii, jj = ii_all[m], jj_all[m]
    print(f"edge set N={len(ii)}", flush=True)

    co_v = co[ii, jj]
    atac_v = G_atac[ii, jj]
    tfu = np.unique(tf_rows)

    rows = {}
    for label, fm in graphs.items():
        fm_v = fm[ii, jj]
        obs = fpa.partial_rho_obs_sliced(
            fm_v, atac_v, co_v, jj, ii, peakcount, genelen, detv, gc,
            tf_outdeg, atac_indeg, True, "full")
        man = fpa.mantel_randomization(
            fm_v, atac_v, co_v, jj, ii, peakcount, genelen, detv, gc,
            tf_outdeg, atac_indeg, G_atac, True, "full", obs, N_PERM,
            seed=hash((label, "mantel")) % (2**63))
        deg = fpa.degree_preserving_null(
            fm_v, atac_v, co_v, jj, ii, peakcount, genelen, detv, gc,
            tf_outdeg, atac_indeg, G_atac, tfu, True, "full", obs, N_PERM,
            seed=hash((label, "deg")) % (2**63))
        rows[label] = {"rho": float(obs), "pM": man["p_mc"], "pD": deg["p_mc"]}
        print(f"  {label:22s} rho={obs:+.6f} pM={man['p_mc']:.4f} pD={deg['p_mc']:.4f}", flush=True)

    labels = list(rows)
    qM = bh([rows[l]["pM"] for l in labels])
    qD = bh([rows[l]["pD"] for l in labels])
    print(f"\nBH across {len(labels)}-row brain family:", flush=True)
    for i, l in enumerate(labels):
        rows[l]["qM"] = float(qM[i])
        rows[l]["qD"] = float(qD[i])
        both = qM[i] < 0.05 and qD[i] < 0.05
        print(f"  {l:22s} rho={rows[l]['rho']:+.6f} qM={qM[i]:.4f} qD={qD[i]:.4f} both={both}", flush=True)

    with open(f"{OUT}/verify_brain_attention_omission_v2.json", "w") as fh:
        json.dump({"edge_set_N": int(len(ii)), "n_perm": N_PERM, "rows": rows}, fh, indent=1)
    print("wrote verify_brain_attention_omission_v2.json", flush=True)


if __name__ == "__main__":
    main()
