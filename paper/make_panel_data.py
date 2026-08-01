#!/usr/bin/env python3
"""Derived panel data for the extended figure set. Deterministic: every value is
either copied from a pinned results/v2 JSON or recomputed from pinned npz graphs
on the audit's own edge set. Run from paper/:
    python3 make_panel_data.py
Writes paper/panel_data.json, consumed by make_figs.R.
"""
import json
import os
import numpy as np
from scipy.stats import rankdata, spearmanr

base = os.environ.get("SCFM_BASE", "..")
res = os.path.join(base, "results", "v2")
out = {}


def j(name):
    with open(os.path.join(res, name)) as fh:
        return json.load(fh)


def edge_vector(G, tf_rows, n_genes, self_gene_idx, mask_col_idx=None):
    """Flatten TF->target weights over the audit edge set: rows = tf_rows,
    all gene columns, excluding self pairs and (for brain) the masked marker
    column (VCAN -- the only prespecified marker present in the panel)."""
    sub = G[np.asarray(tf_rows)][:, :]
    keep = np.ones(sub.shape, dtype=bool)
    for i, t in enumerate(tf_rows):
        keep[i, self_gene_idx[t]] = False
    if mask_col_idx is not None:
        keep[:, mask_col_idx] = False
    return sub[keep]


def load_panel(npz_name):
    z = np.load(os.path.join(res, npz_name), allow_pickle=False)
    genes = [str(g) for g in z["genes"]]
    return z, genes


# ---------- usability: FM-vs-coexpression Spearman per pooled row ----------
brain_z, brain_genes = load_panel("G_ATAC_v2_GSE174367.npz")
pbmc_z, pbmc_genes = load_panel("G_ATAC_v2_PBMC10k.npz")
brain_gene_idx = {g: i for i, g in enumerate(brain_genes)}
pbmc_gene_idx = {g: i for i, g in enumerate(pbmc_genes)}
# marker mask: VCAN is the only one of the 34 prespecified markers in the panel
brain_mask_col = brain_gene_idx.get("VCAN")

fm = np.load(os.path.join(res, "fmgraphs_pooled_v2.npz"))
ko = np.load(os.path.join(res, "G_ko_v2.npz"))
floor = np.load(os.path.join(res, "brain_floor_graph_v2.npz"))
psg = np.load(os.path.join(res, "pbmc_scgpt_pooled_v2.npz"))
puc = np.load(os.path.join(res, "pbmc_uce_pooled_v2.npz"))

brain_tf_rows = brain_z["tf_rows"]
pbmc_tf_rows = pbmc_z["tf_rows"]
brain_self = {t: t for t in brain_tf_rows}  # self pair = same row/col index
pbmc_self = {t: t for t in pbmc_tf_rows}


def spearman_on_edges(A, B, tf_rows, gene_idx, mask_col=None):
    self_idx = {t: t for t in tf_rows}
    a = edge_vector(A, tf_rows, len(gene_idx), self_idx, mask_col)
    b = edge_vector(B, tf_rows, len(gene_idx), self_idx, mask_col)
    return float(spearmanr(a, b).statistic)


co_b = fm["co"]
computed = {}
computed["brain_random_init_floor"] = spearman_on_edges(
    floor["G"], co_b, brain_tf_rows, brain_gene_idx, brain_mask_col)
computed["brain_geneformer_ko_posctrl"] = spearman_on_edges(
    ko["G_ko_ctrl"], co_b, brain_tf_rows, brain_gene_idx, brain_mask_col)
computed["pbmc_scGPT_encoder"] = spearman_on_edges(
    psg["sg"], psg["co"], pbmc_tf_rows, pbmc_gene_idx, None)
computed["pbmc_UCE_encoder"] = spearman_on_edges(
    puc["uce"], puc["co"], pbmc_tf_rows, pbmc_gene_idx, None)

crossmodal = j("crossmodal_v2.json")["observed"]
crossmodal_scf = j("crossmodal_scf_v2.json")["observed"]
crossmodal_uce = j("crossmodal_uce_v2.json")["observed"]
attn = j("readout_attention_v2.json")["observed"]
kostat = j("insilico_ko_v2.json")
pb_eval = j("pbmc_eval_v2.json")
pb_scf = j("pbmc_eval_scf_v2.json")

usability = {
    "brain": {
        "geneformer_embed": crossmodal["geneformer_vs_coexp"],
        "scFoundation_encoder": crossmodal_scf["scf_vs_coexp"],
        "UCE_encoder": crossmodal_uce["uce_vs_coexp"],
        "scGPT_encoder": crossmodal["scgpt_vs_coexp"],
        "geneformer_attn": attn["attn_sym_vs_coexp"],
        "random_init_floor": computed["brain_random_init_floor"],
        "geneformer_ko_raw": kostat["ko_vs_coexp"],
        "geneformer_ko_posctrl": computed["brain_geneformer_ko_posctrl"],
    },
    "pbmc": {
        "geneformer_embed": pb_eval["embed__fm_vs_coexp"],
        "geneformer_attn": pb_eval["attn__fm_vs_coexp"],
        "scFoundation_encoder": pb_scf["scf_vs_coexp"],
        "scGPT_encoder": computed["pbmc_scGPT_encoder"],
        "UCE_encoder": computed["pbmc_UCE_encoder"],
    },
}
out["usability_fm_vs_coexp"] = usability
out["usability_provenance"] = {
    "json": ["crossmodal_v2", "crossmodal_scf_v2", "crossmodal_uce_v2",
             "readout_attention_v2", "insilico_ko_v2", "pbmc_eval_v2", "pbmc_eval_scf_v2"],
    "recomputed_on_audit_edge_set": list(computed.keys()),
}

# ---------- panel composition: binary target-profile sharing ----------
def composition(npz_name):
    z = np.load(os.path.join(res, npz_name))
    tf_rows = z["tf_rows"]
    type_keys = [k for k in z.keys() if k.startswith("G_")]
    consensus = np.mean([z[k] for k in type_keys], axis=0)
    profiles = (consensus[tf_rows] > 0)
    n_tf = profiles.shape[0]
    # identical profiles
    n_identical = 0
    for i in range(n_tf):
        eq = np.all(profiles == profiles[i], axis=1)
        if eq.sum() > 1:
            n_identical += 1
    # binary cosine partners > 0.8
    P = profiles.astype(float)
    norms = np.linalg.norm(P, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    C = (P @ P.T) / (norms @ norms.T)
    np.fill_diagonal(C, 0.0)
    n_partner = int((C.max(axis=1) > 0.8).sum())
    return {"n_tf": int(n_tf), "n_identical_profile": int(n_identical),
            "n_partner_cosine_gt_0.8": n_partner,
            "definition": "binary (>0) profiles of the per-type-averaged consensus proxy"}


out["panel_composition"] = {
    "brain": composition("G_ATAC_v2_GSE174367.npz"),
    "pbmc": composition("G_ATAC_v2_PBMC10k.npz"),
}

# ---------- motif evidence per tissue ----------
motif = {}
for tag, meta_name in [("Brain", "G_ATAC_v2_meta.json"),
                       ("PBMC", "G_ATAC_v2_PBMC10k_meta.json"),
                       ("Fibroblast mix", "G_ATAC_v2_GSE206767_meta.json")]:
    m = j(meta_name)
    motif[tag] = {"relevant_peaks": m["relevant_peaks"],
                  "peak_motif_hits": m["peak_motif_hits"],
                  "hits_per_peak": m["peak_motif_hits"] / m["relevant_peaks"],
                  "expected_random_per_peak": 5.5,
                  "motif_p": m["motif_p"]}
out["motif_evidence"] = motif

# ---------- per-cell-type cell counts ----------
audit = j("fixed_panel_audit_v2.json")
cells = {}
for tissue in ("brain", "pbmc"):
    seen = {}
    for r in audit["per_cell_type"][tissue]["full_confound"]["rows"]:
        seen[r["cell_type"]] = r["n_cells"]
    cells[tissue] = seen
out["pertype_n_cells"] = cells

# ---------- readout resource / coverage table ----------
puc_cells = int(len(puc["cell_ids"])) if "cell_ids" in puc else None
psg_cells = int(len(psg["cell_ids"])) if "cell_ids" in psg else None
out["readout_qc"] = {
    "cells_per_readout": {
        "brain Geneformer attention (capped)": attn_cells if (attn_cells := j("readout_attention_v2.json")["n_cells"]) else None,
        "brain Geneformer KO": j("insilico_ko_v2.json")["n_cells"],
        "brain UCE": j("crossmodal_uce_v2.json")["n_cells"],
        "brain scFoundation": j("crossmodal_scf_v2.json")["n_cells"],
        "PBMC pool (all readouts)": psg_cells,
    },
    "manifest_gene_coverage": {
        "brain UCE": j("crossmodal_uce_v2.json")["manifest_genes_covered"],
        "brain scFoundation": j("crossmodal_scf_v2.json")["manifest_genes_covered"],
        "PBMC scFoundation": pb_scf["manifest_genes_covered"],
        "PBMC UCE": int(np.sum(puc["covered"])) if "covered" in puc else None,
        "all other readouts": 1200,
    },
    "edge_accounting": {
        "tf_x_gene": 446 * 1200,
        "minus_self_pairs": 446,
        "pbmc_edges": 534754,
        "brain_marker_mask_edges": 446,
        "brain_edges": 534308,
    },
}

with open(os.path.join("panel_data.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out["usability_fm_vs_coexp"], indent=1))
print(json.dumps(out["panel_composition"], indent=1))
print("panel_data.json written")
