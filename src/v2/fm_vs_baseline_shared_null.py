#!/usr/bin/env python
"""G002 Option B — shared-null randomization of Δρ = ρ_FM − ρ_coexp_baseline.

For each tissue × pooled FM readout under the full confound spec:
  observed Δ = ρ_FM(use_coexp=True, full) − ρ_baseline(use_coexp=False, full)
Under each of N shared gene-label Mantel permutations, recompute BOTH rhos with the
SAME permutation stream and form Δ_null. One-sided plus-one Monte Carlo p-value for
Δ > 0 (FM beats baseline); BH within tissue across FM rows.

Baseline rung matches brain/pbmc_coexp_baseline_null.py (degree_only_no_selfpartial).
FM rows match the audit full-confound partial (use_coexp=True).

Seeds (explicit integer contract):
  BRAIN_SHARED_NULL_SEED = 2026081001
  PBMC_SHARED_NULL_SEED  = 2026081002

Env:
  SCREG_SHARED_NULL_NPERM  default 999; smoke with 49
  SCREG_DATA_ROOT          brain ATAC parent (.../data)
  SCFM_BRAIN_ATAC          override brain ATAC h5ad
  SCREG_PBMC_ATAC          override PBMC ATAC h5ad
  SCREG_GENE_COORDS        override gene_coords_hg38.tsv
  SCREG_SHARED_NULL_OUT_ROOT  output root (default: capsule/workspace root)

Fail-closed: missing required NPZ / ATAC / coords raises FileNotFoundError.
Does not fabricate p-values.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixed_panel_audit as fpa  # noqa: E402
import pbmc_cache  # noqa: E402

# ----------------------------- paths / seeds ---------------------------------
WORKSPACE = Path(__file__).resolve().parents[2]
OUT_ROOT = Path(os.environ.get("SCREG_SHARED_NULL_OUT_ROOT", str(WORKSPACE)))
PUBLIC_JSON = OUT_ROOT / "results" / "fm_vs_baseline_shared_null_v2.public.json"
PRIVATE_JSON = OUT_ROOT / "results" / "v2" / "fm_vs_baseline_shared_null_v2.json"

NPZ_OUT = Path(fpa.OUT)
DATA_ROOT = os.environ.get(
    "SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
# Optional sibling experiment-workspace data root (no machine-specific path literals).
_MONOREPO_CANDIDATES = [
    os.environ.get("SCREG_MONOREPO_DATA", ""),
    os.path.join(str(WORKSPACE.parent), "singlecell-genomics-research",
                 "projects", "scfm-reg-audit", "data"),
]
MONOREPO_DATA = next((p for p in _MONOREPO_CANDIDATES if p and os.path.isdir(p)), "")

ATAC_B = os.environ.get(
    "SCFM_BRAIN_ATAC",
    f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")
_ATAC_P_DEFAULT = f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"
if not os.path.exists(_ATAC_P_DEFAULT) and MONOREPO_DATA:
    _ATAC_P_DEFAULT = f"{MONOREPO_DATA}/multiome/pbmc10k_atac.h5ad"
ATAC_P = os.environ.get("SCREG_PBMC_ATAC", _ATAC_P_DEFAULT)

COORDS = fpa.COORDS
if not os.path.exists(COORDS):
    _coords_fallback = (
        f"{MONOREPO_DATA}/annotation/gene_coords_hg38.tsv" if MONOREPO_DATA else "")
    COORDS = os.environ.get("SCREG_GENE_COORDS", _coords_fallback or fpa.COORDS)

PROM = fpa.PROM
N_PERM = int(os.environ.get("SCREG_SHARED_NULL_NPERM", "999"))
ALPHA = 0.05

BRAIN_SHARED_NULL_SEED = 2026081001
PBMC_SHARED_NULL_SEED = 2026081002
SEED_CONTRACT = "explicit_integer_v1"

# Required NPZ keys relative to results/v2 (fail closed).
REQUIRED_NPZ = {
    "brain": [
        "G_ATAC_v2_GSE174367.npz",
        "fmgraphs_pooled_v2.npz",
        "G_scf_pooled.npz",
        "G_uce_pooled.npz",
        "brain_attention_graph_v2.npz",
        "brain_floor_graph_v2.npz",
        "G_ko_v2.npz",
        "pbmc_confounds_v2.npz",
    ],
    "pbmc": [
        "G_ATAC_v2_PBMC10k.npz",
        "pbmc_fmgraphs_pooled.npz",
        "G_scf_pbmc_pooled.npz",
        "pbmc_scgpt_pooled_v2.npz",
        "pbmc_uce_pooled_v2.npz",
        "pbmc_confounds_v2.npz",
    ],
}


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def require_path(path: str, label: str) -> str:
    if not path or not os.path.exists(path):
        raise FileNotFoundError(
            f"required {label} missing (fail-closed, no fabrication): {path!r}")
    return path


def require_npz(name: str) -> str:
    path = str(NPZ_OUT / name)
    return require_path(path, f"NPZ {name}")


def tissue_peakcount(atac_file: str) -> np.ndarray:
    """Promoter-window peak counts from ATAC peak coordinates (no fasta)."""
    genes, _, _ = fpa.load_manifest()
    gidx = {g: i for i, g in enumerate(genes)}
    gco: Dict[str, Tuple[str, int, int]] = {}
    with open(require_path(COORDS, "gene coords")) as fh:
        for ln in fh:
            c, s, e, st, nm = ln.rstrip("\n").split("\t")
            if nm not in gidx or nm in gco:
                continue
            s, e = int(s), int(e)
            gco[nm] = (c, s - PROM if st == "+" else s, e if st == "+" else e + PROM)
    Av = ad.read_h5ad(require_path(atac_file, "ATAC h5ad"), backed="r")
    peaks = [str(p) for p in Av.var_names]
    Av.file.close()
    pchr = np.array([p.split(":")[0] for p in peaks])
    pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks])
    pmid = (pse[:, 0] + pse[:, 1]) // 2
    by: Dict[str, np.ndarray] = {}
    for i, c in enumerate(pchr):
        by.setdefault(c, []).append(i)
    by = {c: np.asarray(v) for c, v in by.items()}
    pc = np.zeros(len(genes), dtype=np.float32)
    for g, i in gidx.items():
        if g not in gco:
            continue
        c, lo, hi = gco[g]
        pis = by.get(c)
        if pis is not None:
            pc[i] = len(pis[(pmid[pis] >= lo) & (pmid[pis] <= hi)])
    return pc


def load_gene_identity_confounds() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """genelen / detv / gc from PBMC confound cache (tissue-independent)."""
    conf_path = require_npz("pbmc_confounds_v2.npz")
    genes, _, man_sha = fpa.load_manifest()
    _pc, genelen, detv, gc = pbmc_cache.load_confound_cache(conf_path, genes, man_sha)
    return genelen, detv, gc


def make_edge_index(tissue: str, G: np.ndarray, tf: np.ndarray):
    genes, _, _ = fpa.load_manifest()
    Ng = G.shape[0]
    ii0 = np.repeat(tf, Ng)
    jj0 = np.tile(np.arange(Ng), len(tf))
    m = fpa.edge_mask(tissue, genes, tf, ii0, jj0)
    return ii0[m], jj0[m]


def generate_shared_perms(Ng: int, n_perm: int, seed: int) -> np.ndarray:
    """Match batched_mantel_null SeedSequence replicate stream when seed is used."""
    rng_proxy = np.random.SeedSequence(seed)
    perms = np.empty((n_perm, Ng), dtype=np.int64)
    for k in range(n_perm):
        rep_seq = rng_proxy.spawn(1)[0]
        rep_seed = int(rep_seq.generate_state(1, dtype=np.uint64)[0])
        perms[k] = np.random.default_rng(rep_seed).permutation(Ng)
    return perms


def one_sided_plus_one_p(null: np.ndarray, observed: float, n_perm: int) -> dict:
    """Plus-one MC p for H1: Δ > 0  →  p = (#{Δ_null ≥ Δ_obs} + 1) / (N + 1)."""
    count = int(np.sum(null >= observed))
    p_mc = (count + 1) / (n_perm + 1)
    return {
        "p_mc": float(p_mc),
        "N_perm": int(n_perm),
        "resolution": float(1 / (n_perm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "z": float((observed - null.mean()) / (null.std() + 1e-9)),
        "null_obs_count_at_or_above_obs": count,
        "test_type": "shared_mantel_delta_one_sided_plus_one",
        "alternative": "delta_rho > 0 (FM beats coexp baseline)",
    }


def load_brain_bundle() -> dict:
    for name in REQUIRED_NPZ["brain"]:
        require_npz(name)
    Z = np.load(require_npz("G_ATAC_v2_GSE174367.npz"), allow_pickle=False)
    types = [str(t) for t in Z["types"]]
    tf = np.array(Z["tf_rows"])
    G = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
    F = np.load(require_npz("fmgraphs_pooled_v2.npz"))
    co = F["co"].astype(np.float32)
    models = {
        "geneformer_embed": F["gf"].astype(np.float32),
        "scFoundation_encoder": np.load(require_npz("G_scf_pooled.npz"))["G"].astype(np.float32),
        "UCE_encoder": np.load(require_npz("G_uce_pooled.npz"))["G"].astype(np.float32),
        "scGPT_encoder": F["sg"].astype(np.float32),
        "geneformer_attn": np.load(require_npz("brain_attention_graph_v2.npz"))["G_sym"].astype(np.float32),
        "random_init_floor": np.load(require_npz("brain_floor_graph_v2.npz"))["G"].astype(np.float32),
    }
    K = np.load(require_npz("G_ko_v2.npz"))
    models["geneformer_ko_raw"] = K["G_ko"].astype(np.float32)
    models["geneformer_ko_posctrl"] = K["G_ko_ctrl"].astype(np.float32)
    return {"tissue": "brain", "G": G, "co": co, "models": models, "tf": tf,
            "co_by_model": None, "seed": BRAIN_SHARED_NULL_SEED, "atac": ATAC_B}


def load_pbmc_bundle() -> dict:
    for name in REQUIRED_NPZ["pbmc"]:
        require_npz(name)
    Z = np.load(require_npz("G_ATAC_v2_PBMC10k.npz"), allow_pickle=False)
    types = [str(t) for t in Z["types"]]
    tf = np.array(Z["tf_rows"])
    G = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
    Fp = np.load(require_npz("pbmc_fmgraphs_pooled.npz"))
    co = Fp["co"].astype(np.float32)
    models = {
        "geneformer_embed": Fp["gf"].astype(np.float32),
        "geneformer_attn": Fp["at"].astype(np.float32),
        "scFoundation_encoder": np.load(require_npz("G_scf_pbmc_pooled.npz"))["G"].astype(np.float32),
    }
    # Matched co-expression controls for independent-control readouts (audit contract).
    scgpt = np.load(require_npz("pbmc_scgpt_pooled_v2.npz"), allow_pickle=False)
    uce = np.load(require_npz("pbmc_uce_pooled_v2.npz"), allow_pickle=False)
    models["scGPT_encoder"] = scgpt["sg"].astype(np.float32)
    models["UCE_encoder"] = uce["uce"].astype(np.float32)
    co_by_model = {
        "scGPT_encoder": scgpt["co"].astype(np.float32),
        "UCE_encoder": uce["co"].astype(np.float32),
    }
    return {"tissue": "pbmc", "G": G, "co": co, "models": models, "tf": tf,
            "co_by_model": co_by_model, "seed": PBMC_SHARED_NULL_SEED, "atac": ATAC_P}


def run_tissue(bundle: dict, peakcount: np.ndarray, genelen: np.ndarray,
               detv: np.ndarray, gc: np.ndarray, n_perm: int) -> dict:
    tissue = bundle["tissue"]
    G = bundle["G"]
    co = bundle["co"]
    tf = bundle["tf"]
    seed = int(bundle["seed"])
    ii, jj = make_edge_index(tissue, G, tf)
    od = (G > 0).sum(1).astype(np.float32)
    ind = (G > 0).sum(0).astype(np.float32)
    atac_v = G[ii, jj]
    co_v = co[ii, jj]
    Ng = G.shape[0]

    rho_base = fpa.partial_rho_obs_sliced(
        co_v, atac_v, co_v, jj, ii, peakcount, genelen, detv, gc, od, ind,
        False, "full")
    log(f"{tissue}: n_pairs={len(ii)} rho_baseline={rho_base:+.6f} (degree_only_no_selfpartial)")

    log(f"{tissue}: generating {n_perm} shared gene-label perms (seed={seed})...")
    t_perm = time.time()
    perms = generate_shared_perms(Ng, n_perm, seed)
    log(f"{tissue}: perms ready ({time.time() - t_perm:.1f}s)")

    # Baseline null under shared perms (use_coexp=False; predictor = co).
    t0 = time.time()
    base_nulls, base_meta = fpa.batched_mantel_null(
        fm_vecs=[co_v], co_v=co_v, jj=jj, ii=ii,
        peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
        tf_outdeg_full=od, atac_indeg_full=ind, G_atac_full=G,
        use_coexp=False, confound_spec="full",
        n_perm=n_perm, seed=None, perms=perms,
    )
    rho_base_null = base_nulls[0]
    log(f"{tissue}: baseline null done ({time.time() - t0:.1f}s)")

    # Group FM rows by coexp-control source so X_fixed matches the audit partial.
    primary_labels: List[str] = []
    primary_vecs: List[np.ndarray] = []
    matched_groups: Dict[str, List[Tuple[str, np.ndarray, np.ndarray]]] = {}
    co_by_model = bundle.get("co_by_model") or {}

    for label, G_fm in bundle["models"].items():
        fm_v = G_fm[ii, jj].astype(np.float32)
        if label in co_by_model:
            matched_groups.setdefault(label, []).append(
                (label, fm_v, co_by_model[label][ii, jj].astype(np.float32)))
        else:
            primary_labels.append(label)
            primary_vecs.append(fm_v)

    fm_null_by_label: Dict[str, np.ndarray] = {}
    fm_rho_obs: Dict[str, float] = {}
    fm_meta_by_label: Dict[str, dict] = {}

    if primary_vecs:
        t0 = time.time()
        nulls, meta = fpa.batched_mantel_null(
            fm_vecs=primary_vecs, co_v=co_v, jj=jj, ii=ii,
            peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
            tf_outdeg_full=od, atac_indeg_full=ind, G_atac_full=G,
            use_coexp=True, confound_spec="full",
            n_perm=n_perm, seed=None, perms=perms,
        )
        log(f"{tissue}: primary-co FM batch ({len(primary_labels)} rows) "
            f"({time.time() - t0:.1f}s)")
        for label, null, fm_v in zip(primary_labels, nulls, primary_vecs):
            fm_null_by_label[label] = null
            fm_meta_by_label[label] = meta
            fm_rho_obs[label] = fpa.partial_rho_obs_sliced(
                fm_v, atac_v, co_v, jj, ii, peakcount, genelen, detv, gc, od, ind,
                True, "full")

    for label, entries in matched_groups.items():
        # one entry per label by construction
        _lab, fm_v, co_m = entries[0]
        t0 = time.time()
        nulls, meta = fpa.batched_mantel_null(
            fm_vecs=[fm_v], co_v=co_m, jj=jj, ii=ii,
            peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
            tf_outdeg_full=od, atac_indeg_full=ind, G_atac_full=G,
            use_coexp=True, confound_spec="full",
            n_perm=n_perm, seed=None, perms=perms,
        )
        log(f"{tissue}: matched-co FM {label} ({time.time() - t0:.1f}s)")
        fm_null_by_label[label] = nulls[0]
        fm_meta_by_label[label] = meta
        fm_rho_obs[label] = fpa.partial_rho_obs_sliced(
            fm_v, atac_v, co_m, jj, ii, peakcount, genelen, detv, gc, od, ind,
            True, "full")

    rows_pub = []
    rows_priv = []
    # Stable order: dict insertion order from bundle["models"]
    for label in bundle["models"]:
        rho_fm = float(fm_rho_obs[label])
        delta_obs = rho_fm - float(rho_base)
        delta_null = fm_null_by_label[label] - rho_base_null
        summary = one_sided_plus_one_p(delta_null, delta_obs, n_perm)
        summary["seed"] = seed
        summary["shared_proxy_null_batch_id"] = (
            f"shared_delta_seed{seed}_full_n{n_perm}")
        beats = bool(delta_obs > 0)
        rows_pub.append({
            "tissue": tissue,
            "model_label": label,
            "rho_fm_full": round(rho_fm, 6),
            "rho_baseline": float(rho_base),
            "delta_rho": float(delta_obs),
            "beats_baseline_observed": beats,
            "p_mc": summary["p_mc"],
            # bh_q filled after tissue-wide BH
        })
        rows_priv.append({
            **rows_pub[-1],
            "n_pairs": int(len(ii)),
            "null_mean": summary["null_mean"],
            "null_sd": summary["null_sd"],
            "z": summary["z"],
            "null_obs_count_at_or_above_obs": summary["null_obs_count_at_or_above_obs"],
            "test_type": summary["test_type"],
            "alternative": summary["alternative"],
            "resolution": summary["resolution"],
            "coexp_control": ("matched_npz" if label in co_by_model else "tissue_primary"),
            "fm_mantel_batch_id": fm_meta_by_label[label].get("batch_id"),
            "baseline_mantel_batch_id": base_meta.get("batch_id"),
        })

    pvals = [r["p_mc"] for r in rows_pub]
    qvals = fpa.bh_qvalues(pvals)
    for r_pub, r_priv, q in zip(rows_pub, rows_priv, qvals):
        r_pub["bh_q"] = float(q)
        r_pub["shared_null_significant"] = bool(q < ALPHA and r_pub["beats_baseline_observed"])
        r_priv["bh_q"] = float(q)
        r_priv["shared_null_significant"] = r_pub["shared_null_significant"]

    return {
        "tissue": tissue,
        "n_pairs": int(len(ii)),
        "n_perm": int(n_perm),
        "mantel_seed": seed,
        "rho_baseline": float(rho_base),
        "baseline_rung": "degree_only_no_selfpartial",
        "confound_spec": "full",
        "rows_public": rows_pub,
        "rows_private": rows_priv,
        "n_shared_null_significant": sum(1 for r in rows_pub if r["shared_null_significant"]),
        "n_beats_baseline_observed": sum(1 for r in rows_pub if r["beats_baseline_observed"]),
    }


def strip_private_paths(obj):
    """Public JSON must not contain absolute private filesystem paths."""
    if isinstance(obj, dict):
        return {k: strip_private_paths(v) for k, v in obj.items()
                if not (isinstance(v, str) and v.startswith("/home/"))}
    if isinstance(obj, list):
        return [strip_private_paths(x) for x in obj]
    return obj


def main():
    t_wall = time.time()
    log(f"=== FM vs baseline shared-null (N_PERM={N_PERM}) ===")
    log(f"NPZ_OUT={NPZ_OUT}")
    log(f"PUBLIC_JSON={PUBLIC_JSON}")
    log(f"PRIVATE_JSON={PRIVATE_JSON}")

    # Fail closed on required NPZ up front (both tissues).
    missing = []
    for names in REQUIRED_NPZ.values():
        for name in names:
            p = NPZ_OUT / name
            if not p.exists():
                missing.append(str(p))
    if missing:
        raise FileNotFoundError(
            "required NPZ missing (fail-closed, no fabrication):\n  "
            + "\n  ".join(missing))

    genelen, detv, gc = load_gene_identity_confounds()
    log("computing tissue peakcounts...")
    pc_brain = tissue_peakcount(ATAC_B)
    pc_pbmc = tissue_peakcount(ATAC_P)
    log("peakcounts done")

    brain = run_tissue(load_brain_bundle(), pc_brain, genelen, detv, gc, N_PERM)
    pbmc = run_tissue(load_pbmc_bundle(), pc_pbmc, genelen, detv, gc, N_PERM)

    all_pub_rows = brain["rows_public"] + pbmc["rows_public"]
    all_priv_rows = brain["rows_private"] + pbmc["rows_private"]
    wall = time.time() - t_wall

    public_doc = {
        "schema_version": 1,
        "method": "shared_mantel_delta_v1",
        "job": "0x55",
        "goal": "G002_option_B",
        "alternative": "one_sided_delta_gt_0",
        "alternative_note": (
            "p_mc = (#{Δ_null ≥ Δ_obs} + 1) / (N_perm + 1); tests whether FM partial ρ "
            "exceeds the co-expression baseline under shared gene-label Mantel draws."
        ),
        "confound_spec": "full",
        "baseline_rung": "degree_only_no_selfpartial",
        "fm_rung": "full_confound_use_coexp_true",
        "n_perm": int(N_PERM),
        "seeds": {
            "brain": BRAIN_SHARED_NULL_SEED,
            "pbmc": PBMC_SHARED_NULL_SEED,
            "seed_contract": SEED_CONTRACT,
        },
        "bh": {
            "method": "Benjamini-Hochberg",
            "within": "tissue",
            "alpha": ALPHA,
        },
        "wall_time_sec": round(wall, 1),
        "n_fm_rows": len(all_pub_rows),
        "n_beats_baseline_observed": sum(
            1 for r in all_pub_rows if r["beats_baseline_observed"]),
        "n_shared_null_significant": sum(
            1 for r in all_pub_rows if r["shared_null_significant"]),
        "by_tissue": {
            "brain": {
                "rho_baseline": brain["rho_baseline"],
                "n_pairs": brain["n_pairs"],
                "n_shared_null_significant": brain["n_shared_null_significant"],
                "n_beats_baseline_observed": brain["n_beats_baseline_observed"],
            },
            "pbmc": {
                "rho_baseline": pbmc["rho_baseline"],
                "n_pairs": pbmc["n_pairs"],
                "n_shared_null_significant": pbmc["n_shared_null_significant"],
                "n_beats_baseline_observed": pbmc["n_beats_baseline_observed"],
            },
        },
        "summary": (
            f"{sum(1 for r in all_pub_rows if r['shared_null_significant'])}/"
            f"{len(all_pub_rows)} FM rows significant for Δ>0 under shared Mantel "
            f"(BH q<{ALPHA} within tissue); "
            f"{sum(1 for r in all_pub_rows if r['beats_baseline_observed'])}/"
            f"{len(all_pub_rows)} beat baseline on the observed point estimate."
        ),
        "rows": all_pub_rows,
    }
    public_doc = strip_private_paths(public_doc)

    private_doc = {
        **public_doc,
        "schema_version": 1,
        "method": "shared_mantel_delta_v1",
        "npz_root_basename": "results/v2",
        "rows_enriched": all_priv_rows,
    }
    private_doc["commands"] = {
        "smoke_nperm_49": (
            "SCREG_SHARED_NULL_NPERM=49 SCREG_DATA_ROOT=<data> "
            "SCREG_PBMC_ATAC=<pbmc_atac.h5ad> SCREG_GENE_COORDS=<gene_coords_hg38.tsv> "
            "python src/v2/fm_vs_baseline_shared_null.py"
        ),
        "production_nperm_999": (
            "SCREG_SHARED_NULL_NPERM=999 SCREG_DATA_ROOT=<data> "
            "SCREG_PBMC_ATAC=<pbmc_atac.h5ad> SCREG_GENE_COORDS=<gene_coords_hg38.tsv> "
            "python src/v2/fm_vs_baseline_shared_null.py"
        ),
    }

    fpa.write_json_atomic(str(PUBLIC_JSON), public_doc, indent=2)
    fpa.write_json_atomic(str(PRIVATE_JSON), private_doc, indent=2)
    log(f"wrote {PUBLIC_JSON}")
    log(f"wrote {PRIVATE_JSON}")
    log(f"DONE wall={wall:.1f}s  summary={public_doc['summary']}")


if __name__ == "__main__":
    main()
