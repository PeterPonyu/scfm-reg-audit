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


def require(condition, message):
    """Contract gate for the derived values; unlike assert it survives `python -O`."""
    if not condition:
        raise ValueError(message)


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

# ---------- PBMC scGPT/UCE covariate ladder (rows absent from marginal_vs_adjusted_v2) ----------
# Recomputed with the audit's own method: rank-transform graph weights, residualize
# on [intercept + ranked co-expression + standardized covariates], Pearson of residuals.
# Method is validated in-line against the authoritative audit JSON (binary proxy
# degrees; full and non-degree rungs must match to 1e-4) before any value is emitted.
from scipy.stats import rankdata

_gz = np.load(os.path.join(res, "G_ATAC_v2_PBMC10k.npz"))
_tf_rows = _gz["tf_rows"]
_types = [k for k in _gz.keys() if k.startswith("G_")]
_proxy_full = np.mean([_gz[k] for k in _types], axis=0)[_tf_rows]
_n_tf, _n_g = _proxy_full.shape
_self = np.ones((_n_tf, _n_g), bool)
for _i, _t in enumerate(_tf_rows):
    _self[_i, _t] = False
_proxy = _proxy_full[_self]
_co4 = np.load(os.path.join(res, "pbmc_confounds_v2.npz"))
_gene_covs = {k: _co4[k] for k in ("peakcount", "genelen", "detv", "gc")}
_deg_tf = (_proxy_full > 0).sum(1)
_deg_tg = (_proxy_full > 0).sum(0)


def _vec(x, kind):
    if kind == "tf":
        return np.repeat(x, _n_g)[_self.ravel()]
    return np.tile(x, (_n_tf, 1))[_self]


def _z(x):
    return (x - x.mean()) / x.std()


def _partial(fmv, cols):
    fv = rankdata(fmv)
    pv = rankdata(_proxy)
    X = np.column_stack([np.ones_like(fv)] + cols)
    bf = np.linalg.lstsq(X, fv, rcond=None)[0]
    bp = np.linalg.lstsq(X, pv, rcond=None)[0]
    rf = fv - X @ bf
    rp = pv - X @ bp
    return float(np.corrcoef(rf, rp)[0, 1])


def _ladder(fm_vec, co_vec):
    cols_nondeg = [_z(_vec(_gene_covs[k], "tg")) for k in ("peakcount", "genelen", "detv", "gc")]
    cols_deg = [_z(_vec(_deg_tf, "tf")), _z(_vec(_deg_tg, "tg"))]
    co = rankdata(co_vec)
    return {
        "marginal": float(spearmanr(fm_vec, _proxy).statistic),
        "coexp_only": _partial(fm_vec, [co]),
        "nondegree_only": _partial(fm_vec, cols_nondeg),
        "degree_only": _partial(fm_vec, cols_deg),
        "coexp_plus_nondegree": _partial(fm_vec, [co] + cols_nondeg),
        "coexp_plus_full": _partial(fm_vec, [co] + cols_nondeg + cols_deg),
    }


_psg = np.load(os.path.join(res, "pbmc_scgpt_pooled_v2.npz"))
_puc = np.load(os.path.join(res, "pbmc_uce_pooled_v2.npz"))
lad_sg = _ladder(_psg["sg"][_tf_rows][_self], _psg["co"][_tf_rows][_self])
lad_uc = _ladder(_puc["uce"][_tf_rows][_self], _puc["co"][_tf_rows][_self])

_audit = j("fixed_panel_audit_v2.json")
_want = {}
for _x in _audit["pooled"]["pbmc"]["rows"]:
    if _x["row_type"] == "pooled_fm" and _x["model_label"] in ("scGPT_encoder", "UCE_encoder"):
        _want[(_x["model_label"], _x["confound_spec"])] = _x["observed_partial_rho"]
for _lad, _ml in ((lad_sg, "scGPT_encoder"), (lad_uc, "UCE_encoder")):
    require(abs(_lad["coexp_plus_full"] - _want[(_ml, "full")]) < 1e-4,
            f"{_ml} full-spec ladder disagrees with fixed_panel_audit_v2: {_lad}")
    require(abs(_lad["coexp_plus_nondegree"] - _want[(_ml, "non_degree")]) < 1e-4,
            f"{_ml} non-degree ladder disagrees with fixed_panel_audit_v2: {_lad}")
out["ladder_pbmc_extra"] = {
    "definition": "same method as marginal_vs_adjusted_v2 (binary proxy degrees); "
                  "validated in-line against fixed_panel_audit_v2 full and non-degree rows",
    "scGPT_encoder": lad_sg,
    "UCE_encoder": lad_uc,
}

# ---------- alpha-equivalents for PBMC scGPT/UCE (absent from effect_vs_injection_scale_v2) ----------
# Same linear interpolation on the injection curve (main + subdivided points). Method
# validated in-line against every stored INTERPOLATED alpha_equivalent before use.
_es = j("effect_vs_injection_scale_v2.json")
_sub = j("injection_subdivided_v2.json")


def _curve(tn):
    pts = [(p["alpha"], p["mean_rho"]) for p in _es["alpha_to_rho_curves"][tn]["points"]]
    for r in _sub[tn]["rows"]:
        m = float(np.mean([z["observed_partial_rho_axis_aligned"] for z in r["replicate_runs"]]))
        pts.append((r["alpha"], m))
    pts = sorted(set(pts))
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def _alpha_of(rho, tn):
    a, y = _curve(tn)
    return float(np.interp(rho, y, a))


_errs = []
for _r in _es["observed_effects_as_alpha"]:
    if _r["alpha_equivalent"] is not None:
        _errs.append(abs(_alpha_of(_r["observed_rho"], _r["tissue"]) - _r["alpha_equivalent"]))
require(_errs, "no stored INTERPOLATED alpha_equivalent to validate the curve against")
require(max(_errs) < 2e-3,
        f"alpha-equivalent interpolation disagrees with stored values: max err {max(_errs)}")
out["alpha_equiv_extra"] = {
    "method": "linear interpolation on injection curve incl. subdivided points; "
              "validated against all stored INTERPOLATED values (max err %.4f)" % max(_errs),
    "pbmc_scGPT_encoder": _alpha_of(_want[("scGPT_encoder", "full")], "pbmc"),
    "pbmc_UCE_encoder": _alpha_of(_want[("UCE_encoder", "full")], "pbmc"),
}

# ---------- third-tissue transfer (fig11) ----------
# Consensus proxy per tissue = mean over per-cell-type graphs, the exact
# reduction used by src/v2/cross_tissue_additive_decomp.py (load_consensus).
# Edge set = 446 TF rows x 1200 genes minus self pairs (534,308 edges),
# matching the decomposition's "fixed non-self TF-target edge set".
_tt_meta = {"brain": "G_ATAC_v2_meta.json",
            "pbmc": "G_ATAC_v2_PBMC10k_meta.json",
            "fibro": "G_ATAC_v2_GSE206767_meta.json"}
_tt_npz = {"brain": "G_ATAC_v2_GSE174367.npz",
           "pbmc": "G_ATAC_v2_PBMC10k.npz",
           "fibro": "G_ATAC_v2_GSE206767.npz"}
_tt_names = {"brain": "Brain", "pbmc": "PBMC", "fibro": "Fibroblast mix"}
# total_peaks = n_vars of the pinned source h5ad (SHA-256 in docs/FULL_RERUN.md);
# PBMC total 111,743 is the value independently recomputed by reviewer 2.
_tt_total_peaks = {"brain": 219070, "pbmc": 111743, "fibro": 275448}


def _consensus_support(npz_name):
    z = np.load(os.path.join(res, npz_name), allow_pickle=False)
    types = [str(t) for t in z["types"]]
    G = np.mean([z[f"G_{t}"] for t in types], axis=0).astype(np.float64)
    tr = np.asarray(z["tf_rows"])
    sub = G[tr, :].copy()
    jidx = np.arange(G.shape[1])
    sub[jidx[tr[:, None]] == jidx[None, :]] = 0.0  # self pairs
    return sub != 0


_tt_sup = {k: _consensus_support(v) for k, v in _tt_npz.items()}
_tt_keys = ["brain", "pbmc", "fibro"]
_regions = []
for _r in range(1, 4):
    for _combo in __import__("itertools").combinations(_tt_keys, _r):
        _m = np.logical_and.reduce([_tt_sup[c] for c in _combo])
        for c in _tt_keys:
            if c not in _combo:
                _m &= ~_tt_sup[c]
        _regions.append({"combo": [_tt_names[c] for c in _combo], "n": int(_m.sum())})

_dec = j("cross_tissue_additive_decomp_v2.json")
_pairs = {"GSE174367": "Brain", "PBMC10k": "PBMC", "GSE206767": "Fibroblast mix"}
out["third_tissue"] = {
    "provenance": "consensus = mean over per-cell-type graphs (same as "
                  "cross_tissue_additive_decomp.py); meta from G_ATAC_v2*_meta.json; "
                  "total_peaks = n_vars of pinned h5ad (FULL_RERUN.md hashes)",
    "coverage": [{"tissue": _tt_names[k],
                  "tag": {"brain": "GSE174367", "pbmc": "PBMC10k", "fibro": "GSE206767"}[k],
                  "relevant_peaks": j(_tt_meta[k])["relevant_peaks"],
                  "total_peaks": _tt_total_peaks[k],
                  "motif_hits": j(_tt_meta[k])["peak_motif_hits"]}
                 for k in _tt_keys],
    "edge_overlap": {
        "regions": _regions,
        "supported": {_tt_names[k]: int(_tt_sup[k].sum()) for k in _tt_keys},
        "union": int(sum(r["n"] for r in _regions)),
        "panel_edges_nonself": int(_tt_sup["brain"].size - 446),
    },
    "degree_tf_out": {_tt_names[k]: [int(x) for x in np.sort(_tt_sup[k].sum(1))]
                      for k in _tt_keys},
    "degree_gene_in": {_tt_names[k]: [int(x) for x in np.sort(_tt_sup[k].sum(0))]
                       for k in _tt_keys},
    "rho_phi": [{"pair": "--".join(_pairs[t] for t in row["pair"]),
                 "observed": row["observed_spearman"],
                 "phi": row["binary_support_phi"]}
                for row in _dec["rows"]],
}
require(out["third_tissue"]["edge_overlap"]["panel_edges_nonself"] == 534754,
        f"non-self panel edge count changed: "
        f"{out['third_tissue']['edge_overlap']['panel_edges_nonself']}")
require(abs(sum(r["n"] for r in _regions) - out["third_tissue"]["edge_overlap"]["union"]) < 1e-9,
        "Venn region counts do not sum to the edge-overlap union")

with open(os.path.join("panel_data.json"), "w") as fh:
    json.dump(out, fh, indent=1, allow_nan=False)
print(json.dumps(out["usability_fm_vs_coexp"], indent=1))
print(json.dumps(out["panel_composition"], indent=1))
print("ladder scGPT:", {k: round(v, 5) for k, v in lad_sg.items()})
print("ladder UCE:  ", {k: round(v, 5) for k, v in lad_uc.items()})
print("alpha equiv extra:", out["alpha_equiv_extra"])
print("panel_data.json written")