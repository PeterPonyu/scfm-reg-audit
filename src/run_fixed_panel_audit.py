#!/usr/bin/env python
"""
scfm-reg-audit v2 — driver: produce the fixed-TF-panel audit, replacing the
bootstrap-CI / MDE / superpopulation interpretation retired in LEGACY_INFERENCE_NOTE.md.

Writes (only after full-test pass per blocker gate):
  results/v2/fixed_panel_audit_v2.json         (pooled + per-type + cross-tissue rows)
  results/v2/fixed_panel_signal_injection_v2.json  (axis-aligned sensitivity)
  results/v2/inference_status_v2.json          (machine-readable status + hashes)
"""
import json
import os
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixed_panel_audit as fpa
import pbmc_cache
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))


# ----------------------------- defaults --------------------------------------
SEED_ROOT = int(os.environ.get("SEED_ROOT", "20260724"))
N_PERM_POOLED_MANTEL = int(os.environ.get("N_PERM_POOLED_MANTEL", "999"))
N_PERM_POOLED_DEG = int(os.environ.get("N_PERM_POOLED_DEG", "999"))
N_REPLICATES = int(os.environ.get("N_REPLICATES", "30"))
ALPHAS = [0.0, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.75, 1.0]


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def spawn_int_seeds(ss, n: int) -> list:
    """Spawn n SeedSequence children and return n Python int seeds via generate_state(1)[0]."""
    children = ss.spawn(n)
    return [int(s.generate_state(1, dtype=np.uint64)[0]) for s in children]


# ----------------------------- cached graph loaders --------------------------
def load_pooled_brain():
    Z = np.load(f"{fpa.OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
    types_b = [str(t) for t in Z["types"]]
    tf_b = np.array(Z["tf_rows"])
    G_atac = np.mean([Z[f"G_{t}"] for t in types_b], axis=0).astype(np.float32)
    F = np.load(f"{fpa.OUT}/fmgraphs_pooled_v2.npz")
    Sb = np.load(f"{fpa.OUT}/G_scf_pooled.npz")["G"]
    Ub = np.load(f"{fpa.OUT}/G_uce_pooled.npz")["G"]
    Ab = np.load(f"{fpa.OUT}/brain_attention_graph_v2.npz")["G_sym"]
    FLOORb = np.load(f"{fpa.OUT}/brain_floor_graph_v2.npz")["G"]
    models = {
        "geneformer_embed": F["gf"],
        "scFoundation_encoder": Sb,
        "UCE_encoder": Ub,
        "scGPT_encoder": F["sg"],  # present in pooled brain only
        # Brain Geneformer attention and the random-init floor are cached graphs
        # that were computed but never wired into the audit family. Omitting a
        # computed readout is a selective-inclusion hazard; include both.
        "geneformer_attn": Ab,
        "random_init_floor": FLOORb,
    }
    return G_atac, F["co"], models, tf_b, types_b


def load_pooled_pbmc():
    Z = np.load(f"{fpa.OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=False)
    types_p = [str(t) for t in Z["types"]]
    tf_p = np.array(Z["tf_rows"])
    G_atac = np.mean([Z[f"G_{t}"] for t in types_p], axis=0).astype(np.float32)
    Fp = np.load(f"{fpa.OUT}/pbmc_fmgraphs_pooled.npz")
    Sp = np.load(f"{fpa.OUT}/G_scf_pbmc_pooled.npz")["G"]
    models = {
        "geneformer_embed": Fp["gf"],
        "geneformer_attn": Fp["at"],
        "scFoundation_encoder": Sp,
    }
    return G_atac, Fp["co"], models, tf_p, types_p


def load_optional_pbmc_scgpt():
    path = f"{fpa.OUT}/pbmc_scgpt_pooled_v2.npz"
    if not os.path.exists(path):
        return None
    rna_path = f"{fpa.ROOT}/data/multiome/pbmc10k_rna.h5ad"
    import anndata as ad
    rna = ad.read_h5ad(rna_path, backed="r")
    expected_cell_ids = pbmc_cache.select_pool_cell_ids(rna.n_obs, 4000, 20260713)
    rna.file.close()
    expected_rna_sha = pbmc_cache.sha256_file(rna_path)
    with np.load(path, allow_pickle=False) as cache:
        genes = json.loads(Path(fpa.MANI).read_text())["genes"]
        manifest_sha = json.loads(Path(fpa.MANI).read_text())["sha256"]
        if [str(g) for g in cache["genes"]] != genes:
            raise ValueError("PBMC scGPT graph gene order mismatch")
        if str(cache["manifest_sha"].item()) != manifest_sha:
            raise ValueError("PBMC scGPT graph manifest mismatch")
        if int(cache["selection_seed"].item()) != 20260713:
            raise ValueError("PBMC scGPT graph selection seed mismatch")
        if int(cache["pool_cap"].item()) != 4000:
            raise ValueError("PBMC scGPT graph pool cap mismatch")
        if not np.array_equal(cache["cell_ids"], expected_cell_ids):
            raise ValueError("PBMC scGPT graph selected cells mismatch")
        if str(cache["rna_sha256"].item()) != expected_rna_sha:
            raise ValueError("PBMC scGPT graph RNA input mismatch")
        co, sg = cache["co"].copy(), cache["sg"].copy()
    expected_shape = (len(genes), len(genes))
    if co.shape != expected_shape or sg.shape != expected_shape:
        raise ValueError("PBMC scGPT graph shape mismatch")
    if not np.isfinite(co).all() or not np.isfinite(sg).all():
        raise ValueError("PBMC scGPT graph contains non-finite values")
    return path, co, sg


def load_optional_pbmc_uce():
    """UCE PBMC pooled graph with its matched co-expression control.

    pbmc_uce_pooled_v2.npz carries the same manifest/selection provenance as the
    scGPT cache plus the UCE checkpoint and ESM2 hashes; the graph was computed
    (07-29) but never appended to the audit family, so it is loaded here with the
    same provenance gates as scGPT.
    """
    path = f"{fpa.OUT}/pbmc_uce_pooled_v2.npz"
    if not os.path.exists(path):
        return None
    rna_path = f"{fpa.ROOT}/data/multiome/pbmc10k_rna.h5ad"
    import anndata as ad
    rna = ad.read_h5ad(rna_path, backed="r")
    expected_cell_ids = pbmc_cache.select_pool_cell_ids(rna.n_obs, 4000, 20260713)
    rna.file.close()
    expected_rna_sha = pbmc_cache.sha256_file(rna_path)
    with np.load(path, allow_pickle=False) as cache:
        genes = json.loads(Path(fpa.MANI).read_text())["genes"]
        manifest_sha = json.loads(Path(fpa.MANI).read_text())["sha256"]
        if [str(g) for g in cache["genes"]] != genes:
            raise ValueError("PBMC UCE graph gene order mismatch")
        if str(cache["manifest_sha"].item()) != manifest_sha:
            raise ValueError("PBMC UCE graph manifest mismatch")
        if int(cache["selection_seed"].item()) != 20260713:
            raise ValueError("PBMC UCE graph selection seed mismatch")
        if int(cache["pool_cap"].item()) != 4000:
            raise ValueError("PBMC UCE graph pool cap mismatch")
        if not np.array_equal(cache["cell_ids"], expected_cell_ids):
            raise ValueError("PBMC UCE graph selected cells mismatch")
        if str(cache["rna_sha256"].item()) != expected_rna_sha:
            raise ValueError("PBMC UCE graph RNA input mismatch")
        co, uce = cache["co"].copy(), cache["uce"].copy()
    expected_shape = (len(genes), len(genes))
    if co.shape != expected_shape or uce.shape != expected_shape:
        raise ValueError("PBMC UCE graph shape mismatch")
    if not np.isfinite(co).all() or not np.isfinite(uce).all():
        raise ValueError("PBMC UCE graph contains non-finite values")
    return path, co, uce


def load_geneformer_ko():
    K = np.load(f"{fpa.OUT}/G_ko_v2.npz")
    return {
        "geneformer_ko_raw": K["G_ko"],
        "geneformer_ko_posctrl": K["G_ko_ctrl"],
    }, np.array(K["tf_rows"])


def _report_absent_pertype_caches(tissue: str, missing: list, label: str) -> None:
    """Report cell types whose cached graph is absent.

    Skipping them silently shrinks the per-type readout set without any trace in the
    logs or the artifact, so the omission is announced on both channels.
    """
    if not missing:
        return
    message = (f"{tissue}: {len(missing)} cell types have no {label} cache and are excluded "
               f"from the per-type readouts: {sorted(missing)}")
    log(f"WARNING {message}")
    warnings.warn(message, RuntimeWarning, stacklevel=2)


def load_brain_pertype_models():
    Z = np.load(f"{fpa.OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
    types_b = [str(t) for t in Z["types"]]
    tf_b = np.array(Z["tf_rows"])
    out: dict = {}
    missing_fm, missing_scf = [], []
    for t in types_b:
        fpath = f"{fpa.OUT}/brain_fmgraphs_{t}.npz"
        if not os.path.exists(fpath):
            missing_fm.append(t)
            continue
        F = np.load(fpath)
        spath = f"{fpa.OUT}/brain_scfgraphs_{t}.npz"
        scf = np.load(spath)["scf"] if os.path.exists(spath) else None
        out[t] = {
            "geneformer_embed": F["gf"],
            "geneformer_attn": F["at"],
            "coexp": F["co"],
        }
        if scf is None:
            missing_scf.append(t)
        else:
            out[t]["scFoundation_encoder"] = scf
    _report_absent_pertype_caches("brain", missing_fm, "brain_fmgraphs")
    _report_absent_pertype_caches("brain", missing_scf, "brain_scfgraphs")
    return out, tf_b, types_b


def load_pbmc_pertype_models():
    Z = np.load(f"{fpa.OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=False)
    types_p = [str(t) for t in Z["types"]]
    tf_p = np.array(Z["tf_rows"])
    out: dict = {}
    missing_fm = []
    for t in types_p:
        fpath = f"{fpa.OUT}/pbmc_fmgraphs_{t}.npz"
        if not os.path.exists(fpath):
            missing_fm.append(t)
            continue
        F = np.load(fpath)
        out[t] = {
            "geneformer_embed": F["gf"],
            "geneformer_attn": F["at"],
            "coexp": F["co"],
        }
    _report_absent_pertype_caches("pbmc", missing_fm, "pbmc_fmgraphs")
    return out, tf_p, types_p


def load_cross_tissue():
    Zg = np.load(f"{fpa.OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
    Zp = np.load(f"{fpa.OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=False)
    Zc = np.load(f"{fpa.OUT}/G_ATAC_v2_GSE206767.npz", allow_pickle=False)
    tags = ["GSE174367", "PBMC10k", "GSE206767"]
    consensus = {}
    tfs = {}
    for Z, tag in [(Zg, "GSE174367"), (Zp, "PBMC10k"), (Zc, "GSE206767")]:
        ts = [str(t) for t in Z["types"]]
        consensus[tag] = np.mean([Z[f"G_{t}"] for t in ts], axis=0).astype(np.float32)
        tfs[tag] = np.array(Z["tf_rows"])
    return consensus, tfs, tags


# ----------------------------- KO readout tagging ----------------------------
# Readout labelling is separated from FM family labelling so Geneformer readouts can
# all be BH-adjusted as one family while still carrying an explicit readout label.
_FAMILY_LOOKUP = {
    "geneformer_embed": ("geneformer", "embed"),
    "geneformer_attn": ("geneformer", "attn"),
    "geneformer_ko_raw": ("geneformer", "ko_raw"),
    "geneformer_ko_posctrl": ("geneformer", "ko_posctrl"),
    "scFoundation_encoder": ("scFoundation", "encoder"),
    "UCE_encoder": ("UCE", "encoder"),
    "scGPT_encoder": ("scGPT", "encoder"),
    "co_expression": ("co_expression", "marginal"),
    "random_init_floor": ("random_init", "floor"),
}


# ----------------------------- row builders ---------------------------------
def build_row(
    tissue: str, model_label: str, confound_spec: str,
    fm_v: np.ndarray, atac_v: np.ndarray, co_v: np.ndarray,
    G_atac_full: np.ndarray, tf_rows: np.ndarray,
    ii: np.ndarray, jj: np.ndarray,
    peakcount: np.ndarray, genelen: np.ndarray, gc: np.ndarray, detv: np.ndarray,
    tf_outdeg_full: np.ndarray, atac_indeg_full: np.ndarray,
    use_coexp: bool, n_perm_mantel: int, n_perm_deg: int,
    seed_mantel: int, seed_deg: int,
) -> dict:
    """Compute observed + Mantel + degree-preserving null for one row using the
    BATCHED shared-null helpers (one proxy null per replicate, paired against this row).
    Returns the same row shape as before."""
    family_label, readout = _FAMILY_LOOKUP[model_label]
    observed = fpa.partial_rho_obs_sliced(
        fm_v=fm_v, atac_v=atac_v, co_v=co_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg=tf_outdeg_full, atac_indeg=atac_indeg_full,
        use_coexp=use_coexp, confound_spec=confound_spec,
    )
    # BATCHED Mantel: shared proxy null across rows; this call computes atac_resid per
    # perm ONCE and pairs against fm_v.
    nulls_m, meta_m = fpa.batched_mantel_null(
        fm_vecs=[fm_v], co_v=co_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=tf_outdeg_full, atac_indeg_full=atac_indeg_full,
        G_atac_full=G_atac_full, use_coexp=use_coexp, confound_spec=confound_spec,
        n_perm=n_perm_mantel, seed=seed_mantel,
    )
    mantel = fpa.batched_pvalue_summary(
        nulls_m[0], observed, n_perm=n_perm_mantel,
        seed=seed_mantel, batch_id=meta_m["batch_id"],
        test_type="gene_label_mantel_plus_one_corrected", confound_spec=confound_spec,
    )
    mantel.update({
        "null_columns_perm_recomputed": meta_m["null_columns_perm_recomputed"],
        "null_columns_fixed_under_perm": meta_m["null_columns_fixed_under_perm"],
        "shared_proxy_null_batch_id": meta_m["batch_id"],
        "shared_proxy_null_per_replicate": meta_m["shared_proxy_null_per_replicate"],
        "replicate_seed_stream": meta_m["replicate_seed_stream"],
        "n_rows_in_shared_batch": meta_m["n_rows_in_batch"],
    })
    nulls_d, meta_d = fpa.batched_degree_preserving_null(
        fm_vecs=[fm_v], co_v=co_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=tf_outdeg_full, atac_indeg_full=atac_indeg_full,
        G_atac_full=G_atac_full, tf_rows_unique=tf_rows,
        use_coexp=use_coexp, confound_spec=confound_spec,
        n_perm=n_perm_deg, seed=seed_deg,
    )
    deg = fpa.batched_pvalue_summary(
        nulls_d[0], observed, n_perm=n_perm_deg,
        seed=seed_deg, batch_id=meta_d["batch_id"],
        test_type="degree_preserving_row_shuffle_plus_one_corrected", confound_spec=confound_spec,
    )
    deg.update({
        "null_columns_perm_recomputed": meta_d["null_columns_perm_recomputed"],
        "null_columns_fixed_under_perm": meta_d["null_columns_fixed_under_perm"],
        "shared_proxy_null_batch_id": meta_d["batch_id"],
        "shared_proxy_null_per_replicate": meta_d["shared_proxy_null_per_replicate"],
        "replicate_seed_stream": meta_d["replicate_seed_stream"],
        "n_rows_in_shared_batch": meta_d["n_rows_in_batch"],
    })
    return {
        "tissue": tissue,
        "model_label": model_label,
        "model_family": family_label,
        "readout": readout,
        "confound_spec": confound_spec,
        "n_pairs": int(len(ii)),
        "observed_partial_rho": round(float(observed), 6),
        "use_coexp_in_partial": bool(use_coexp),
        "mantel": mantel,
        "degree_preserving": deg,
    }
def run_pooled_family(ss, atac_file: str, tissue: str, G_atac, co, models: dict,
                      tf_rows: np.ndarray, ko_models: dict | None,
                      n_perm_mantel: int, n_perm_deg: int) -> dict:
    """Returns one tissue's pooled primary + sensitivity rows + provenance.
    Batched: shared proxy null per replicate is computed ONCE per (spec, null type)
    across ALL FM rows in the family (no per-row recomputation of atac_perm or
    atac_resid)."""
    Ng = G_atac.shape[0]
    peakcount, genelen, gc, detv = fpa.build_confounds(atac_file)
    tf_outdeg = (G_atac > 0).sum(1).astype(np.float32)
    atac_indeg = (G_atac > 0).sum(0).astype(np.float32)
    rows: list = []

    ii_all = np.repeat(tf_rows, Ng)
    jj_all = np.tile(np.arange(Ng), len(tf_rows))
    m_all = fpa.edge_mask(tissue, json.loads(Path(fpa.MANI).read_text())["genes"], tf_rows, ii_all, jj_all)
    ii, jj = ii_all[m_all], jj_all[m_all]

    fm_entries: list = list(models.items())
    if ko_models is not None:
        fm_entries.append(("geneformer_ko_raw", ko_models["geneformer_ko_raw"]))
        fm_entries.append(("geneformer_ko_posctrl", ko_models["geneformer_ko_posctrl"]))

    coexp_v = co[ii, jj]
    atac_v = G_atac[ii, jj]
    coexp_baseline_observed = float(spearmanr(rankdata(coexp_v), rankdata(atac_v)).statistic)
    rows.append({
        "row_type": "coexp_baseline_marginal",
        "tissue": tissue,
        "model_label": "co_expression",
        "model_family": "co_expression",
        "readout": "marginal",
        "confound_spec": "n/a",
        "n_pairs": int(len(ii)),
        "observed_marginal_rho": round(coexp_baseline_observed, 6),
        "note": "co-expression is the baseline, NOT a model; reported separately; not in BH.",
    })

    # Build per-row observed + fm_vecs ONCE (shared between primary and sensitivity families)
    fm_vecs_per_row = []
    observed_per_row = []
    for model_label, G_fm in fm_entries:
        fm_v = G_fm[ii, jj]
        fm_vecs_per_row.append(fm_v)
        # Observed for primary (full) and sensitivity (non_degree) are distinct
        obs_full = fpa.partial_rho_obs_sliced(
            fm_v=fm_v, atac_v=atac_v, co_v=coexp_v, jj=jj, ii=ii,
            peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
            tf_outdeg=tf_outdeg, atac_indeg=atac_indeg,
            use_coexp=True, confound_spec="full",
        )
        obs_nd = fpa.partial_rho_obs_sliced(
            fm_v=fm_v, atac_v=atac_v, co_v=coexp_v, jj=jj, ii=ii,
            peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
            tf_outdeg=tf_outdeg, atac_indeg=atac_indeg,
            use_coexp=True, confound_spec="non_degree",
        )
        observed_per_row.append({"full": obs_full, "non_degree": obs_nd})

    # Two seeds per spec × 2 specs × 2 null types = 4 seeds total
    seeds = spawn_int_seeds(ss, 4)
    seed_m_full, seed_m_nd, seed_d_full, seed_d_nd = seeds

    # BATCHED: ONE shared proxy null per spec × null_type across ALL FM rows
    mantel_full_nulls, meta_mf = fpa.batched_mantel_null(
        fm_vecs=fm_vecs_per_row, co_v=coexp_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
        G_atac_full=G_atac, use_coexp=True, confound_spec="full",
        n_perm=n_perm_mantel, seed=seed_m_full,
    )
    mantel_nd_nulls, meta_mn = fpa.batched_mantel_null(
        fm_vecs=fm_vecs_per_row, co_v=coexp_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
        G_atac_full=G_atac, use_coexp=True, confound_spec="non_degree",
        n_perm=n_perm_mantel, seed=seed_m_nd,
    )
    deg_full_nulls, meta_df = fpa.batched_degree_preserving_null(
        fm_vecs=fm_vecs_per_row, co_v=coexp_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
        G_atac_full=G_atac, tf_rows_unique=tf_rows,
        use_coexp=True, confound_spec="full",
        n_perm=n_perm_deg, seed=seed_d_full,
    )
    deg_nd_nulls, meta_dn = fpa.batched_degree_preserving_null(
        fm_vecs=fm_vecs_per_row, co_v=coexp_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
        G_atac_full=G_atac, tf_rows_unique=tf_rows,
        use_coexp=True, confound_spec="non_degree",
        n_perm=n_perm_deg, seed=seed_d_nd,
    )

    # Build rows
    def _summary(null_arr, observed, n_perm, seed, batch_meta, test_type, family_id):
        s = fpa.batched_pvalue_summary(
            null_arr, observed, n_perm=n_perm, seed=seed,
            batch_id=batch_meta["batch_id"],
            test_type=test_type, confound_spec=batch_meta["confound_spec"],
        )
        s.update({
            "family_id": family_id,
            "null_columns_perm_recomputed": batch_meta["null_columns_perm_recomputed"],
            "null_columns_fixed_under_perm": batch_meta["null_columns_fixed_under_perm"],
            "shared_proxy_null_batch_id": batch_meta["batch_id"],
            "shared_proxy_null_per_replicate": batch_meta["shared_proxy_null_per_replicate"],
            "replicate_seed_stream": batch_meta["replicate_seed_stream"],
            "n_rows_in_shared_batch": batch_meta["n_rows_in_batch"],
        })
        return s

    for r_idx, (model_label, G_fm) in enumerate(fm_entries):
        family_label, readout = _FAMILY_LOOKUP[model_label]
        for spec, nulls_arr, null_meta, seed_used, test_type, fam_suffix in [
            ("full", mantel_full_nulls[r_idx], meta_mf, seed_m_full,
             "gene_label_mantel_plus_one_corrected", "mantel"),
            ("non_degree", mantel_nd_nulls[r_idx], meta_mn, seed_m_nd,
             "gene_label_mantel_plus_one_corrected", "mantel"),
        ]:
            family_id_mantel = f"{tissue}_pooled_{spec}_confound_{fam_suffix}"
            mantel = _summary(nulls_arr, observed_per_row[r_idx][spec], n_perm_mantel,
                              seed_used, null_meta, test_type, family_id_mantel)
            deg_arr = deg_full_nulls[r_idx] if spec == "full" else deg_nd_nulls[r_idx]
            deg_meta = meta_df if spec == "full" else meta_dn
            deg_seed = seed_d_full if spec == "full" else seed_d_nd
            family_id_deg = f"{tissue}_pooled_{spec}_confound_degree"
            deg = _summary(deg_arr, observed_per_row[r_idx][spec], n_perm_deg,
                           deg_seed, deg_meta,
                           "degree_preserving_row_shuffle_plus_one_corrected",
                           family_id_deg)
            row = {
                "row_type": "pooled_fm",
                "tissue": tissue,
                "model_label": model_label,
                "model_family": family_label,
                "readout": readout,
                "confound_spec": spec,
                "family_id_mantel": family_id_mantel,
                "family_id_degree": family_id_deg,
                "n_pairs": int(len(ii)),
                "observed_partial_rho": round(float(observed_per_row[r_idx][spec]), 6),
                "use_coexp_in_partial": True,
                "mantel": mantel,
                "degree_preserving": deg,
            }
            rows.append(row)

    # BH within each (tissue, spec) family
    primary = [r for r in rows if r["row_type"] == "pooled_fm" and r["confound_spec"] == "full"]
    sensitivity = [r for r in rows if r["row_type"] == "pooled_fm" and r["confound_spec"] == "non_degree"]
    p_mc_pri = [r["mantel"]["p_mc"] for r in primary]
    p_mc_pri_deg = [r["degree_preserving"]["p_mc"] for r in primary]
    for r, q in zip(primary, fpa.bh_qvalues(p_mc_pri)):
        r["mantel"]["bh_q_family"] = round(q, 6)
    for r, q in zip(primary, fpa.bh_qvalues(p_mc_pri_deg)):
        r["degree_preserving"]["bh_q_family"] = round(q, 6)
    p_mc_sen = [r["mantel"]["p_mc"] for r in sensitivity]
    p_mc_sen_deg = [r["degree_preserving"]["p_mc"] for r in sensitivity]
    for r, q in zip(sensitivity, fpa.bh_qvalues(p_mc_sen)):
        r["mantel"]["bh_q_family"] = round(q, 6)
    for r, q in zip(sensitivity, fpa.bh_qvalues(p_mc_sen_deg)):
        r["degree_preserving"]["bh_q_family"] = round(q, 6)

    return {
        "tissue": tissue,
        "n_pairs": int(len(ii)),
        "panel_n_tf": int(len(tf_rows)),
        "panel_n_genes": int(Ng),
        "marker_mask_applied_to_tissue": tissue in fpa.MARKER_TISSUES,
        "rows": rows,
        "primary_family": {
            "family_id": f"{tissue}_pooled_full_confound_primary",
            "description": (
                "All FM rows at full confound spec for THIS tissue. BH within this family. "
                "Geneformer-embed/attn/ko_raw/ko_posctrl are 4 readouts within the 1 "
                "Geneformer family."
            ),
            "n_rows": len(primary),
            "bh_q_method": "Benjamini-Hochberg, p_mc from plus-one Monte-Carlo randomization",
            "rows": primary,
        },
        "sensitivity_family": {
            "family_id": f"{tissue}_pooled_non_degree_confound_sensitivity",
            "description": (
                "All FM rows at non-degree confound spec for THIS tissue. BH within this "
                "family. Mirrors the preregistered primary/sensitivity separation."
            ),
            "n_rows": len(sensitivity),
            "bh_q_method": "Benjamini-Hochberg, p_mc from plus-one Monte-Carlo randomization",
            "rows": sensitivity,
        },
        "provenance": _pooled_provenance(tissue, atac_file, fm_entries),
    }


def _pooled_provenance(tissue: str, atac_file: str, fm_entries: list) -> dict:
    """Compute full SHA-256 provenance for every cached graph used by a pooled row."""
    tag = "GSE174367" if tissue == "brain" else "PBMC10k"
    G_atac_path = f"{fpa.OUT}/G_ATAC_v2_{tag}.npz"
    co_path = (
        f"{fpa.OUT}/fmgraphs_pooled_v2.npz" if tissue == "brain"
        else f"{fpa.OUT}/pbmc_fmgraphs_pooled.npz"
    )
    row_provenance: dict = {}
    for model_label, _G in fm_entries:
        if model_label == "geneformer_embed":
            if tissue == "brain":
                p, k = co_path, "gf"
            else:
                p, k = co_path, "gf"
        elif model_label == "geneformer_attn":
            if tissue == "brain":
                p, k = f"{fpa.OUT}/brain_attention_graph_v2.npz", "G_sym"
            else:
                p, k = co_path, "at"
        elif model_label == "random_init_floor":
            p, k = f"{fpa.OUT}/brain_floor_graph_v2.npz", "G"
        elif model_label in ("geneformer_ko_raw", "geneformer_ko_posctrl"):
            p, k = f"{fpa.OUT}/G_ko_v2.npz", "G_ko" if model_label == "geneformer_ko_raw" else "G_ko_ctrl"
        elif model_label == "scFoundation_encoder":
            if tissue == "brain":
                p, k = f"{fpa.OUT}/G_scf_pooled.npz", "G"
            else:
                p, k = f"{fpa.OUT}/G_scf_pbmc_pooled.npz", "G"
        elif model_label == "UCE_encoder":
            p, k = f"{fpa.OUT}/G_uce_pooled.npz", "G"
        elif model_label == "scGPT_encoder":
            p, k = co_path, "sg"
        else:
            continue
        row_provenance[model_label] = fpa.matrix_provenance(p, k)
    return {
        "G_atac_pooled_source": G_atac_path,
        "G_atac_pooled_sha256": fpa.sha256_file(G_atac_path),
        "confound_source": atac_file,
        "manifest_sha256": json.loads(Path(fpa.MANI).read_text())["sha256"],
        "model_label_provenance": row_provenance,
    }


def recompute_pooled_bh(result: dict) -> None:
    for family_key in ("primary_family", "sensitivity_family"):
        rows = result[family_key]["rows"]
        for null_key in ("mantel", "degree_preserving"):
            qvalues = fpa.bh_qvalues([row[null_key]["p_mc"] for row in rows])
            for row, qvalue in zip(rows, qvalues):
                row[null_key]["bh_q_family"] = round(qvalue, 6)
        result[family_key]["n_rows"] = len(rows)


def append_independent_control_model(
    result: dict, atac_file: str, tissue: str, G_atac: np.ndarray,
    co: np.ndarray, model: np.ndarray, tf_rows: np.ndarray,
    model_label: str, n_perm_mantel: int, n_perm_deg: int,
    seed_sequence: np.random.SeedSequence, model_path: str,
    model_key: str = "sg", control_key: str = "co",
) -> None:
    Ng = G_atac.shape[0]
    peakcount, genelen, gc, detv = fpa.build_confounds(atac_file)
    tf_outdeg = (G_atac > 0).sum(1).astype(np.float32)
    atac_indeg = (G_atac > 0).sum(0).astype(np.float32)
    genes = json.loads(Path(fpa.MANI).read_text())["genes"]
    ii_all = np.repeat(tf_rows, Ng)
    jj_all = np.tile(np.arange(Ng), len(tf_rows))
    mask = fpa.edge_mask(tissue, genes, tf_rows, ii_all, jj_all)
    ii, jj = ii_all[mask], jj_all[mask]
    seeds = spawn_int_seeds(seed_sequence, 4)
    seed_m_full, seed_m_nd, seed_d_full, seed_d_nd = seeds
    for spec, seed_mantel, seed_deg, family_key in (
        ("full", seed_m_full, seed_d_full, "primary_family"),
        ("non_degree", seed_m_nd, seed_d_nd, "sensitivity_family"),
    ):
        row = build_row(
            tissue=tissue, model_label=model_label, confound_spec=spec,
            fm_v=model[ii, jj], atac_v=G_atac[ii, jj], co_v=co[ii, jj],
            G_atac_full=G_atac, tf_rows=tf_rows, ii=ii, jj=jj,
            peakcount=peakcount, genelen=genelen, gc=gc, detv=detv,
            tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
            use_coexp=True, n_perm_mantel=n_perm_mantel, n_perm_deg=n_perm_deg,
            seed_mantel=seed_mantel, seed_deg=seed_deg,
        )
        family_id_mantel = f"{tissue}_pooled_{spec}_confound_mantel"
        family_id_degree = f"{tissue}_pooled_{spec}_confound_degree"
        row.update({
            "row_type": "pooled_fm",
            "family_id_mantel": family_id_mantel,
            "family_id_degree": family_id_degree,
            "coexp_control_source": f"{model_path}:co",
        })
        row["mantel"]["family_id"] = family_id_mantel
        row["degree_preserving"]["family_id"] = family_id_degree
        result["rows"].append(row)
        result[family_key]["rows"].append(row)
    result["provenance"]["model_label_provenance"][model_label] = fpa.matrix_provenance(
        model_path, model_key)
    result["provenance"].setdefault("model_control_provenance", {})[model_label] = (
        fpa.matrix_provenance(model_path, control_key))
    recompute_pooled_bh(result)


# ----------------------------- pooled provenance -----------------------------
def run_pertype_family(ss, atac_file: str, tissue: str, G_atac_pooled: np.ndarray,
                       co_pooled: np.ndarray, type_models: dict, tf_rows: np.ndarray,
                       confound_spec: str, cell_count_lookup: dict | None) -> dict:
    """Per-cell-type rows for ONE confound spec. DESCRIPTIVE EXPLORATORY ROBUSTNESS
    ONLY: report observed effect sizes, n_pairs, n_cells, marker-mask status, and
    signed-summary stats (range, median, sign counts across readouts). NO
    randomization, NO p_mc, NO q, NO BH, NO 'significance' or 'chance' language.
    This eliminates ~70/81s of N=9 wall time per the runtime gate."""
    Ng = G_atac_pooled.shape[0]
    peakcount, genelen, gc, detv = fpa.build_confounds(atac_file)
    rows: list = []
    genes = json.loads(Path(fpa.MANI).read_text())["genes"]
    atac_path = f"{fpa.OUT}/G_ATAC_v2_{ 'GSE174367' if tissue=='brain' else 'PBMC10k' }.npz"
    Z_atac = np.load(atac_path, allow_pickle=False)

    for ctype, models in type_models.items():
        G_atac_type = Z_atac[f"G_{ctype}"].astype(np.float32)
        tf_outdeg_t = (G_atac_type > 0).sum(1).astype(np.float32)
        atac_indeg_t = (G_atac_type > 0).sum(0).astype(np.float32)
        ii_all = np.repeat(tf_rows, Ng)
        jj_all = np.tile(np.arange(Ng), len(tf_rows))
        m = fpa.edge_mask(tissue, genes, tf_rows, ii_all, jj_all)
        ii, jj = ii_all[m], jj_all[m]
        atac_v = G_atac_type[ii, jj]
        co_v = co_pooled[ii, jj]
        n_cells = int((cell_count_lookup or {}).get(ctype, 0))

        coexp_baseline_observed = float(spearmanr(rankdata(models["coexp"][ii, jj]),
                                                    rankdata(atac_v)).statistic)
        rows.append({
            "row_type": "coexp_baseline_marginal",
            "tissue": tissue, "cell_type": ctype,
            "model_label": "co_expression", "model_family": "co_expression",
            "readout": "marginal", "confound_spec": "n/a",
            "n_pairs": int(len(ii)), "n_cells": n_cells,
            "observed_marginal_rho": round(coexp_baseline_observed, 6),
            "note": "co-expression baseline; not in BH; descriptive only.",
        })

        for model_label, G in models.items():
            if model_label == "coexp":
                continue
            family_label, readout = _FAMILY_LOOKUP[model_label]
            fm_v = G[ii, jj]
            observed = fpa.partial_rho_obs_sliced(
                fm_v=fm_v, atac_v=atac_v, co_v=co_v, jj=jj, ii=ii,
                peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
                tf_outdeg=tf_outdeg_t, atac_indeg=atac_indeg_t,
                use_coexp=True, confound_spec=confound_spec,
            )
            rows.append({
                "row_type": "pertype_fm",
                "tissue": tissue, "cell_type": ctype,
                "model_label": model_label, "model_family": family_label,
                "readout": readout, "confound_spec": confound_spec,
                "n_pairs": int(len(ii)), "n_cells": n_cells,
                "observed_partial_rho": round(float(observed), 6),
                "marker_mask_applied": tissue in fpa.MARKER_TISSUES,
                "use_coexp_in_partial": True,
                "note": "descriptive exploratory robustness only; no randomization, no p_mc.",
            })

    fm_rows = [r for r in rows if r.get("row_type") == "pertype_fm"]
    # Descriptive signed summaries across FM rows in this (tissue, spec) family
    rhos = np.array([r["observed_partial_rho"] for r in fm_rows], dtype=float)
    summary = {
        "n_rows_exploratory": int(len(fm_rows)),
        "rho_min": float(rhos.min()) if len(rhos) else None,
        "rho_max": float(rhos.max()) if len(rhos) else None,
        "rho_median": float(np.median(rhos)) if len(rhos) else None,
        "rho_mean": float(rhos.mean()) if len(rhos) else None,
        "rho_std": float(rhos.std()) if len(rhos) else None,
        "n_positive": int(np.sum(rhos > 0)),
        "n_negative": int(np.sum(rhos < 0)),
        "n_zero": int(np.sum(rhos == 0)),
    }
    return {
        "tissue": tissue,
        "confound_spec": confound_spec,
        "n_rows_total": len(rows),
        "n_rows_exploratory": len(fm_rows),
        "rows": rows,
        "descriptive_summary": summary,
    }


# ----------------------------- cross-tissue ----------------------------------
def run_cross_tissue(consensus: dict, tfs: dict) -> dict:
    """3 cross-tissue fixed-panel OBSERVED-SPEARMAN pairs only. No randomization
    (the legacy cross_tissue_bootstrap_v2.json retains a Mantel p_mc which can be
    cited as a legacy sensitivity analysis but is not regenerated here)."""
    Ng = next(iter(consensus.values())).shape[0]
    rows = []
    pairs = [("GSE174367", "PBMC10k"), ("GSE174367", "GSE206767"), ("PBMC10k", "GSE206767")]
    for (a, b) in pairs:
        tf_common = np.intersect1d(tfs[a], tfs[b])
        ii = np.repeat(tf_common, Ng); jj = np.tile(np.arange(Ng), len(tf_common)); m = ii != jj
        ii, jj = ii[m], jj[m]
        x, y = consensus[a][ii, jj], consensus[b][ii, jj]
        observed = float(spearmanr(x, y).statistic)
        rows.append({
            "pair": [a, b],
            "n_tf_common": int(len(tf_common)),
            "observed_spearman": round(observed, 6),
            "provenance": {
                "consensus_a_source": f"results/v2/G_ATAC_v2_{a}.npz",
                "consensus_a_sha256": fpa.sha256_file(f"{fpa.OUT}/G_ATAC_v2_{a}.npz"),
                "consensus_b_source": f"results/v2/G_ATAC_v2_{b}.npz",
                "consensus_b_sha256": fpa.sha256_file(f"{fpa.OUT}/G_ATAC_v2_{b}.npz"),
            },
            "note": "fixed-panel observed Spearman; no randomization in new audit; legacy cross_tissue_bootstrap_v2.json retains a Mantel p_mc for cross-reference.",
        })
    return {
        "family_id": "cross_tissue_construct_reproducibility",
        "description": "Reproducibility of the regulatory_potential_proxy across the 3 ATAC datasets; fixed-panel observed-Spearman only.",
        "n_rows": len(rows),
        "rows": rows,
    }


# ----------------------------- axis-aligned sensitivity ----------------------
def run_sensitivity(ss, tissue: str, G_atac_pooled: np.ndarray, co_pooled: np.ndarray,
                    tf_rows: np.ndarray, atac_file: str, n_replicates: int) -> dict:
    """Per (tissue, alpha, replicate): inject alpha*z(resid) + (1-alpha)*z(noise). NO CI."""
    Ng = G_atac_pooled.shape[0]
    peakcount, genelen, gc, detv = fpa.build_confounds(atac_file)
    tf_outdeg = (G_atac_pooled > 0).sum(1).astype(np.float32)
    atac_indeg = (G_atac_pooled > 0).sum(0).astype(np.float32)
    ii_all = np.repeat(tf_rows, Ng); jj_all = np.tile(np.arange(Ng), len(tf_rows))
    m = fpa.edge_mask(tissue, json.loads(Path(fpa.MANI).read_text())["genes"], tf_rows, ii_all, jj_all)
    ii, jj = ii_all[m], jj_all[m]
    atac_v = G_atac_pooled[ii, jj]
    co_v = co_pooled[ii, jj]

    pc_z = fpa.zscore(peakcount[jj]); gl_z = fpa.zscore(genelen[jj])
    dv_z = fpa.zscore(detv[jj]); gc_z = fpa.zscore(gc[jj])
    od_z = fpa.zscore(tf_outdeg[ii]); ai_z = fpa.zscore(atac_indeg[jj])
    C_full = np.column_stack([pc_z, gl_z, dv_z, gc_z, od_z, ai_z])

    C_with_co = np.column_stack([rankdata(co_v), C_full])
    atac_resid = fpa.resid(rankdata(atac_v), C_with_co)

    n_seeds = len(ALPHAS) * n_replicates
    seeds = spawn_int_seeds(ss, n_seeds)
    seed_iter = iter(seeds)
    rows = []
    for alpha in ALPHAS:
        per_alpha = []
        for rep in range(n_replicates):
            seed = int(next(seed_iter))
            rng = np.random.default_rng(seed)
            noise = rng.standard_normal(len(atac_resid))
            synth = alpha * fpa.zscore(atac_resid) + (1.0 - alpha) * fpa.zscore(noise)
            observed = fpa.partial_rho_obs_sliced(
                fm_v=synth, atac_v=atac_v, co_v=co_v, jj=jj, ii=ii,
                peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
                tf_outdeg=tf_outdeg, atac_indeg=atac_indeg,
                use_coexp=True, confound_spec="full",
            )
            per_alpha.append({
                "alpha": float(alpha),
                "replicate": int(rep),
                "seed": int(seed),
                "n_pairs": int(len(ii)),
                "observed_partial_rho_axis_aligned": round(float(observed), 6),
                "design": "axis_aligned_pipeline_sensitivity",
            })
        rows.append({"alpha": float(alpha), "replicate_runs": per_alpha})
    return {
        "tissue": tissue,
        "regulatory_potential_proxy_basis": f"results/v2/G_ATAC_v2_{ 'GSE174367' if tissue=='brain' else 'PBMC10k' }.npz",
        "regulatory_potential_proxy_basis_sha256": fpa.sha256_file(
            f"{fpa.OUT}/G_ATAC_v2_{ 'GSE174367' if tissue=='brain' else 'PBMC10k' }.npz"),
        "control_coexp_source": ("fmgraphs_pooled_v2.npz:co" if tissue == "brain" else "pbmc_fmgraphs_pooled.npz:co"),
        "confound_spec": "full",
        "alphas": ALPHAS,
        "replicates_per_alpha": int(n_replicates),
        "design": ("axis_aligned_pipeline_sensitivity; alpha is fraction of regulatory_potential_proxy "
                   "residual injected; fresh N(0,1) noise per (tissue,alpha,replicate). NOT a power "
                   "analysis, no exclusion, no coverage, no MDE, no CI, no implied-alpha."),
        "rows": rows,
    }


# ----------------------------- production preflight -------------------------
LEGACY_HASHES = {
    "stats_enhanced_v2.json": "bd84b5af0d81e74739495231ac8a5774f96197253d0efed69e374e50c948b39a",
    "power_analysis_v2.json": "ab6dd2e384ebc3244531cbedfec4a3b5074934881a624a354b3b42b7a52c0e9f",
    "cross_tissue_bootstrap_v2.json": "7f0ffa5a49196df2e843c3080a326b1edec54fbf63c4b5fea9bbc12e4dfbf750",
    "pertype_stats_enhanced_v2.json": "e50f4a552449c47c3a9e14e5787ac49ed981114c6d7f79f05c58b23db475ac86",
}


def verify_legacy_hashes(out_dir: str = fpa.OUT) -> dict:
    actual = {name: fpa.sha256_file(f"{out_dir}/{name}") for name in LEGACY_HASHES}
    mismatches = {
        name: {"expected": LEGACY_HASHES[name], "actual": actual[name]}
        for name in LEGACY_HASHES if actual[name] != LEGACY_HASHES[name]
    }
    if mismatches:
        raise RuntimeError(f"legacy hash mismatch; aborting before analysis: {mismatches}")
    return actual


def validate_fixed_panel_inputs(tf_arrays: dict, matrices: dict) -> dict:
    genes, _detection, manifest_sha256 = fpa.load_manifest()
    n_genes = len(genes)
    if n_genes != 1200:
        raise ValueError(f"fixed panel requires 1200 manifest genes; got {n_genes}")

    canonical_name, canonical_tf = next(iter(tf_arrays.items()))
    canonical_tf = np.asarray(canonical_tf)

    tf_hashes = {}
    for name, tf_rows in tf_arrays.items():
        tf_rows = np.asarray(tf_rows)
        if tf_rows.ndim != 1 or len(tf_rows) != 446:
            raise ValueError(f"{name} must contain 446 TF indices; got shape {tf_rows.shape}")
        if len(np.unique(tf_rows)) != len(tf_rows):
            raise ValueError(f"{name} TF indices are not unique")
        if not np.issubdtype(tf_rows.dtype, np.integer):
            raise ValueError(f"{name} TF indices must be integers")
        if np.any(tf_rows < 0) or np.any(tf_rows >= n_genes):
            raise ValueError(f"{name} TF indices are out of range for {n_genes} genes")
        if not np.array_equal(tf_rows, canonical_tf):
            raise ValueError(f"TF panel mismatch: {name} differs from {canonical_name}")
        tf_hashes[name] = fpa.sha256_array(tf_rows)

    matrix_hashes = {}
    for name, matrix in matrices.items():
        matrix = np.asarray(matrix)
        if matrix.shape != (n_genes, n_genes):
            raise ValueError(f"{name} must have shape {(n_genes, n_genes)}; got {matrix.shape}")
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} contains non-finite values")
        matrix_hashes[name] = fpa.sha256_array(matrix)

    return {
        "n_tf": int(len(canonical_tf)),
        "n_genes": int(n_genes),
        "manifest_sha256": manifest_sha256,
        "tf_panel_sha256": fpa.sha256_array(canonical_tf),
        "tf_array_sha256": tf_hashes,
        "matrix_sha256": matrix_hashes,
    }


def load_cell_count_lookup(path: str, key: str, expected_cell_types) -> dict:
    document = json.loads(Path(path).read_text())
    rows = document.get(key)
    if not isinstance(rows, list):
        raise ValueError(f"{path}:{key} must be a list")
    counts = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("cell_type"), str):
            raise ValueError(f"{path}:{key} has a row without a string cell_type")
        cell_type = row["cell_type"]
        n_cells = row.get("n")
        if cell_type in counts:
            raise ValueError(f"{path}:{key} has duplicate cell_type {cell_type}")
        if not isinstance(n_cells, int) or isinstance(n_cells, bool) or n_cells < 0:
            raise ValueError(f"{path}:{key}:{cell_type} has invalid n={n_cells!r}")
        counts[cell_type] = n_cells
    missing = sorted(set(expected_cell_types) - set(counts))
    if missing:
        raise ValueError(f"{path}:{key} is missing cell counts for {missing}")
    return counts


def exact_bh_family_definitions() -> dict:
    definitions = {}
    for tissue in ("brain", "pbmc"):
        for confound_spec in ("full", "non_degree"):
            for null_type in ("mantel", "degree"):
                family_id = f"{tissue}_pooled_{confound_spec}_confound_{null_type}"
                definitions[family_id] = {
                    "tissue": tissue,
                    "confound_spec": confound_spec,
                    "null_type": null_type,
                    "members": "all pooled FM readout rows for this tissue/specification/null",
                    "method": "Benjamini-Hochberg over plus-one Monte-Carlo p_mc values",
                }
    return definitions


def _stage_json(path: str, document: dict) -> str:
    destination = Path(path)
    fd, temp_path = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent,
    )
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(document, handle, indent=2, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        with open(temp_path) as handle:
            parsed = json.load(handle)
        if parsed != document:
            raise RuntimeError(f"staged JSON validation mismatch for {destination.name}")
        return temp_path
    except BaseException:
        # Cleanup only; the staging failure itself is re-raised unchanged.
        try:
            os.unlink(temp_path)
        except FileNotFoundError:
            pass
        raise


def publish_authoritative_outputs(audit_path: str, audit_doc: dict,
                                  sens_path: str, sens_doc: dict,
                                  status_path: str, status_builder) -> tuple:
    staged = []
    try:
        audit_temp = _stage_json(audit_path, audit_doc)
        staged.append(audit_temp)
        sens_temp = _stage_json(sens_path, sens_doc)
        staged.append(sens_temp)
        audit_hashes = {
            Path(audit_path).name: fpa.sha256_file(audit_temp),
            Path(sens_path).name: fpa.sha256_file(sens_temp),
        }
        status_doc = status_builder(audit_hashes)
        status_temp = _stage_json(status_path, status_doc)
        staged.append(status_temp)
        staged_hashes = {
            **audit_hashes,
            Path(status_path).name: fpa.sha256_file(status_temp),
        }
        os.replace(audit_temp, audit_path)
        staged.remove(audit_temp)
        os.replace(sens_temp, sens_path)
        staged.remove(sens_temp)
        os.replace(status_temp, status_path)
        staged.remove(status_temp)
        return staged_hashes, status_doc
    finally:
        for temp_path in staged:
            try:
                os.unlink(temp_path)
            except FileNotFoundError:
                pass


# ----------------------------- inference_status_v2.json ----------------------
def build_inference_status(audit_hashes: dict, legacy_actual: dict,
                           panel_validation: dict, has_pbmc_scgpt: bool,
                           has_pbmc_uce: bool = False) -> dict:
    status = {
        "schema_version": 1,
        "decision": "retire_bootstrap_CI_and_MDE_from_authoritative_inference",
        "rationale": ("bootstrap CI for the TF-block confounded-controlled partial-rho was demonstrably "
                       "mis-calibrated; Mantel and degree-preserving nulls gave contradictory results on the "
                       "same panel; 'four models' was claimed without comparable per-tissue coverage. The "
                       "fixed-panel audit treats the preregistered 446x1200 TF->target panel as a fixed finite "
                       "graph and reports finite Monte-Carlo randomization p-values plus the degree-preserving "
                       "null side-by-side. No bootstrap CI, no MDE, no coverage, no exclusion claim, no "
                       "implied-alpha, no generalized population inference."),
        "legacy_files_status": {
            k: {
                "status": "retired_not_authoritative",
                "retention": "kept on disk verbatim for audit; not used for new claims",
                "sha256_pinned": LEGACY_HASHES[k],
                "sha256_actual": legacy_actual[k],
                "hash_check": "pass",
            } for k in LEGACY_HASHES
        },
        "new_authoritative_outputs": audit_hashes,
        "seed_root": SEED_ROOT,
        "n_perm_pooled_mantel": N_PERM_POOLED_MANTEL,
        "n_perm_pooled_degree": N_PERM_POOLED_DEG,
        "n_replicates_axis_aligned": N_REPLICATES,
        "panel": panel_validation,
        "p_value_definition": (
            "plus-one-corrected Monte-Carlo randomization p-value: "
            "p_mc = (count(|null| >= |observed|) + 1) / (N_perm + 1). Resolution = 1/(N_perm + 1)."
        ),
        "BH_family_definitions": exact_bh_family_definitions(),
        "non_inferential_families": {
            "pertype_descriptive_only": "Per-type rows are descriptive exploratory robustness only; NO p_mc, NO q, NO BH.",
            "cross_tissue_construct_reproducibility": "Three fixed-panel observed Spearman comparisons only; no randomization.",
        },
        "model_coverage_table": {
            "brain_pooled_FMs": ["geneformer_embed", "geneformer_attn", "geneformer_ko_raw",
                                 "geneformer_ko_posctrl", "scFoundation_encoder", "UCE_encoder",
                                 "scGPT_encoder", "random_init_floor"],
            "brain_per_type_FMs": ["geneformer_embed", "geneformer_attn", "scFoundation_encoder"],
            "pbmc_pooled_FMs": (["geneformer_embed", "geneformer_attn", "scFoundation_encoder"]
                                 + (["scGPT_encoder"] if has_pbmc_scgpt else [])
                                 + (["UCE_encoder"] if has_pbmc_uce else [])),
            "pbmc_per_type_FMs": ["geneformer_embed", "geneformer_attn"],
            "absent_FM_readouts": [
                "UCE per-type",
                "scGPT per-type (no cached per-cell-type scGPT graph)",
            ] + ([] if has_pbmc_scgpt else [
                "scGPT PBMC pooled (no cached PBMC scGPT graph)",
            ]) + ([] if has_pbmc_uce else [
                "UCE PBMC pooled (no cached PBMC UCE graph)",
            ]),
        },
        "scope_verification": {
            "scgpt_pooled_brain_recompute_command": (
                "python3 -c \"import numpy as np; from scipy.stats import rankdata; "
                "ROOT='" + os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")) + "'; "
                "Z=np.load(ROOT+'/results/v2/G_ATAC_v2_GSE174367.npz', allow_pickle=False); "
                "F=np.load(ROOT+'/results/v2/fmgraphs_pooled_v2.npz'); "
                "types=[str(t) for t in Z['types']]; tf=np.array(Z['tf_rows']); "
                "G=np.mean([Z[f'G_{t}'] for t in types], axis=0).astype(np.float32); "
                "ii=np.repeat(tf, 1200); jj=np.tile(np.arange(1200), len(tf)); m=ii!=jj; ii, jj=ii[m], jj[m]; "
                "a=G[ii, jj]; co=F['co'][ii, jj]; sg=F['sg'][ii, jj]; "
                "def partial(x,y,z): rx,ry,rz=rankdata(x),rankdata(y),rankdata(z); "
                "C=np.c_[np.ones_like(rz),rz]; bx=np.linalg.lstsq(C,rx,rcond=None)[0]; "
                "by=np.linalg.lstsq(C,ry,rcond=None)[0]; return float(np.corrcoef(rx-C@bx, ry-C@by)[0,1]); "
                "print(partial(sg, a, co))\""
            ),
            "crossmodal_v2_json_sha256": fpa.sha256_file(f"{fpa.OUT}/crossmodal_v2.json"),
            "crossmodal_v2_json_path": "results/v2/crossmodal_v2.json",
            "fmgraphs_pooled_v2_npz_sha256": fpa.sha256_file(f"{fpa.OUT}/fmgraphs_pooled_v2.npz"),
            "fmgraphs_pooled_v2_sg_matrix_sha256": fpa.sha256_array(
                np.load(f"{fpa.OUT}/fmgraphs_pooled_v2.npz", allow_pickle=False)["sg"]
            ),
        },
    }
    return status


# ----------------------------- main ------------------------------------------
def main():
    log("=== fixed-panel audit ===")
    log(f"SEED_ROOT={SEED_ROOT}  N_PERM_POOLED_MANTEL={N_PERM_POOLED_MANTEL}  "
        f"N_PERM_POOLED_DEG={N_PERM_POOLED_DEG}  N_REPLICATES={N_REPLICATES}")
    legacy_actual = verify_legacy_hashes()
    G_atac_b, co_b, models_b, tf_b, types_b = load_pooled_brain()
    ko_models, ko_tf = load_geneformer_ko()
    G_atac_p, co_p, models_p, tf_p, types_p = load_pooled_pbmc()
    pbmc_scgpt = load_optional_pbmc_scgpt()
    type_models_b, tf_b2, _ = load_brain_pertype_models()
    type_models_p, tf_p2, _ = load_pbmc_pertype_models()
    consensus, tfs, _ = load_cross_tissue()

    tf_arrays = {
        "brain_pooled": tf_b,
        "pbmc_pooled": tf_p,
        "brain_pertype": tf_b2,
        "pbmc_pertype": tf_p2,
        "geneformer_ko": ko_tf,
        **{f"cross_tissue_{name}": rows for name, rows in tfs.items()},
    }
    matrices = {
        "brain_proxy": G_atac_b,
        "brain_coexp": co_b,
        "pbmc_proxy": G_atac_p,
        "pbmc_coexp": co_p,
        **{f"brain_pooled_{name}": matrix for name, matrix in models_b.items()},
        **{f"pbmc_pooled_{name}": matrix for name, matrix in models_p.items()},
        **{f"ko_{name}": matrix for name, matrix in ko_models.items()},
        **{f"cross_tissue_{name}": matrix for name, matrix in consensus.items()},
    }
    brain_proxy_cache = np.load(f"{fpa.OUT}/G_ATAC_v2_GSE174367.npz", allow_pickle=False)
    pbmc_proxy_cache = np.load(f"{fpa.OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=False)
    for cell_type in types_b:
        matrices[f"brain_proxy_{cell_type}"] = brain_proxy_cache[f"G_{cell_type}"]
    for cell_type in types_p:
        matrices[f"pbmc_proxy_{cell_type}"] = pbmc_proxy_cache[f"G_{cell_type}"]
    for cell_type, model_map in type_models_b.items():
        for model_name, matrix in model_map.items():
            matrices[f"brain_pertype_{cell_type}_{model_name}"] = matrix
    for cell_type, model_map in type_models_p.items():
        for model_name, matrix in model_map.items():
            matrices[f"pbmc_pertype_{cell_type}_{model_name}"] = matrix
    panel_validation = validate_fixed_panel_inputs(tf_arrays, matrices)
    ncells_b = load_cell_count_lookup(
        f"{fpa.OUT}/pertype_fm_v2.json", "per_type", type_models_b,
    )
    ncells_p = load_cell_count_lookup(
        f"{fpa.OUT}/pbmc_eval_v2.json", "per_type_coexp", type_models_p,
    )
    log("preflight passed: legacy hashes, fixed panel, graph shapes, and finite values")

    ATAC_B = f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"
    ATAC_P = f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"
    for input_path in (ATAC_B, ATAC_P, fpa.MANI, fpa.COORDS, fpa.HG38):
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"required fixed-panel input missing: {input_path}")

    ss = np.random.SeedSequence(SEED_ROOT)
    log("=== BRAIN pooled primary + sensitivity ===")
    pooled_brain = run_pooled_family(ss.spawn(1)[0], ATAC_B, "brain", G_atac_b, co_b, models_b, tf_b,
                                     ko_models, N_PERM_POOLED_MANTEL, N_PERM_POOLED_DEG)
    log("=== PBMC pooled primary + sensitivity ===")
    pooled_pbmc = run_pooled_family(ss.spawn(1)[0], ATAC_P, "pbmc", G_atac_p, co_p, models_p, tf_p,
                                    None, N_PERM_POOLED_MANTEL, N_PERM_POOLED_DEG)
    if pbmc_scgpt is not None:
        scgpt_path, scgpt_co, scgpt_graph = pbmc_scgpt
        append_independent_control_model(
            pooled_pbmc, ATAC_P, "pbmc", G_atac_p, scgpt_co, scgpt_graph, tf_p,
            "scGPT_encoder", N_PERM_POOLED_MANTEL, N_PERM_POOLED_DEG,
            np.random.SeedSequence([SEED_ROOT, 20260726, 1]), scgpt_path,
        )
        log("PBMC scGPT appended with its matched co-expression control")
    pbmc_uce = load_optional_pbmc_uce()
    if pbmc_uce is not None:
        uce_path, uce_co, uce_graph = pbmc_uce
        append_independent_control_model(
            pooled_pbmc, ATAC_P, "pbmc", G_atac_p, uce_co, uce_graph, tf_p,
            "UCE_encoder", N_PERM_POOLED_MANTEL, N_PERM_POOLED_DEG,
            np.random.SeedSequence([SEED_ROOT, 20260729, 1]), uce_path,
            model_key="uce",
        )
        log("PBMC UCE appended with its matched co-expression control")

    log("=== BRAIN per-type exploratory (full + non-degree, DESCRIPTIVE ONLY) ===")
    pertype_brain_full = run_pertype_family(
        ss.spawn(1)[0], ATAC_B, "brain", G_atac_b, co_b,
        type_models_b, tf_b, "full", ncells_b,
    )
    pertype_brain_nondeg = run_pertype_family(
        ss.spawn(1)[0], ATAC_B, "brain", G_atac_b, co_b,
        type_models_b, tf_b, "non_degree", ncells_b,
    )

    log("=== PBMC per-type exploratory (full + non-degree, DESCRIPTIVE ONLY) ===")
    pertype_pbmc_full = run_pertype_family(
        ss.spawn(1)[0], ATAC_P, "pbmc", G_atac_p, co_p,
        type_models_p, tf_p, "full", ncells_p,
    )
    pertype_pbmc_nondeg = run_pertype_family(
        ss.spawn(1)[0], ATAC_P, "pbmc", G_atac_p, co_p,
        type_models_p, tf_p, "non_degree", ncells_p,
    )

    # Cross-tissue (observed only, no randomization in new audit)
    log("=== cross-tissue construct reproducibility (observed Spearman only) ===")
    cross = run_cross_tissue(consensus, tfs)

    log("=== axis-aligned pipeline sensitivity ===")
    sens_brain = run_sensitivity(ss.spawn(1)[0], "brain", G_atac_b, co_b, tf_b, ATAC_B, N_REPLICATES)
    sens_pbmc = run_sensitivity(ss.spawn(1)[0], "pbmc", G_atac_p, co_p, tf_p, ATAC_P, N_REPLICATES)

    audit_doc = {
        "schema_version": 1,
        "seed_root": SEED_ROOT,
        "panel": panel_validation,
        "marker_mask_policy": {
            "marker_genes": sorted(fpa.MARKER_GENES),
            "marker_tissues": sorted(fpa.MARKER_TISSUES),
            "description": "marker mask applied symmetrically to every (tissue, cell_type) row for tissues in marker_tissues; tissues not in the set are explicitly unmasked.",
        },
        "pooled": {
            "brain": pooled_brain,
            "pbmc": pooled_pbmc,
        },
        "per_cell_type": {
            "brain": {
                "full_confound": pertype_brain_full,
                "non_degree_confound": pertype_brain_nondeg,
            },
            "pbmc": {
                "full_confound": pertype_pbmc_full,
                "non_degree_confound": pertype_pbmc_nondeg,
            },
            "family_design_note": (
                "Per-type is DESCRIPTIVE EXPLORATORY ROBUSTNESS ONLY. No randomization, "
                "no p_mc, no q, no BH. Each (tissue, spec) row carries observed effect "
                "size + descriptive_summary (range/median/sign counts across FM readouts)."
            ),
        },
        "cross_tissue_construct_reproducibility": cross,
        "model_coverage_table": {
            "brain_pooled_FMs": ["geneformer_embed", "geneformer_attn", "geneformer_ko_raw",
                                 "geneformer_ko_posctrl", "scFoundation_encoder", "UCE_encoder",
                                 "scGPT_encoder", "random_init_floor"],
            "brain_per_type_FMs": ["geneformer_embed", "geneformer_attn", "scFoundation_encoder"],
            "pbmc_pooled_FMs": (["geneformer_embed", "geneformer_attn", "scFoundation_encoder"]
                                 + (["scGPT_encoder"] if pbmc_scgpt is not None else [])
                                 + (["UCE_encoder"] if pbmc_uce is not None else [])),
            "pbmc_per_type_FMs": ["geneformer_embed", "geneformer_attn"],
            "absent_FM_readouts": [
                "UCE per-type",
                "scGPT per-type (no cached per-cell-type scGPT graph)",
            ] + ([] if pbmc_scgpt is not None else [
                "scGPT PBMC pooled (no cached PBMC scGPT graph)",
            ]) + ([] if pbmc_uce is not None else [
                "UCE PBMC pooled (no cached PBMC UCE graph)",
            ]),
        },
        "no_population_inference_contract": (
            "Outputs contain no bootstrap_ci95, no mde_*, no clears_zero, no implied_alpha, "
            "no power, no coverage, no exclusion claims. Every p_mc is a plus-one Monte-Carlo "
            "randomization p-value with explicit seed, N_perm, and resolution 1/(N_perm+1)."
        ),
    }
    audit_path = f"{fpa.OUT}/fixed_panel_audit_v2.json"

    sens_doc = {
        "schema_version": 1,
        "seed_root": SEED_ROOT,
        "design": ("axis_aligned_pipeline_sensitivity; alpha is fraction of the "
                   "regulatory_potential_proxy residual injected; fresh independent N(0,1) "
                   "noise per (tissue, alpha, replicate); NOT a power analysis; no exclusion, "
                   "no coverage, no MDE, no CI, no implied-alpha, no generalised population "
                   "inference."),
        "tissues": {
            "brain": sens_brain,
            "pbmc": sens_pbmc,
        },
    }
    sens_path = f"{fpa.OUT}/fixed_panel_signal_injection_v2.json"
    status_path = f"{fpa.OUT}/inference_status_v2.json"
    legacy_actual = verify_legacy_hashes()
    audit_hashes, status = publish_authoritative_outputs(
        audit_path, audit_doc, sens_path, sens_doc, status_path,
        lambda hashes: build_inference_status(
            hashes, legacy_actual, panel_validation, pbmc_scgpt is not None,
            pbmc_uce is not None),
    )
    log(f"SAVED {audit_path}")
    log(f"SAVED {sens_path}")
    log(f"SAVED {status_path}")
    log("=== summary hashes ===")
    for k, v in audit_hashes.items():
        log(f"  {k}: sha256={v}")
    log("legacy sha256 verified: " + ", ".join(
        f"{k}={status['legacy_files_status'][k]['sha256_actual'][:8]}"
        for k in status["legacy_files_status"]))


if __name__ == "__main__":
    main()