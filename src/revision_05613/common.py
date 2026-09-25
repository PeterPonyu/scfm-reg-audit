"""Shared inputs for the ARRAY-D-26-05613 revision analyses.

Loads the same cached graphs, edge sets, confounds and seeds as
src/run_fixed_panel_audit.py, without writing anything to results/v2.
Revision outputs go to results/revision_05613/.
"""
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SRC))
import fixed_panel_audit as fpa  # noqa: E402
import run_fixed_panel_audit as rfa  # noqa: E402

OUT_DIR = Path(fpa.ROOT) / "results" / "revision_05613"

def result_json(name, env_var=None):
    """Resolve an explicit input, then local pipeline output, then public JSON.

    An invalid explicit path is an error, never a reason to use another input.
    No sibling workspace is searched. The chosen bytes are hashed by callers.
    """
    configured = os.environ.get(env_var) if env_var else None
    if configured is not None:
        path = Path(configured)
        if not path.is_file():
            raise FileNotFoundError(f"{env_var} input missing: {path}")
        return path
    candidates = [Path(fpa.OUT) / name,
                  Path(fpa.ROOT) / "results" / (Path(name).stem + ".public.json")]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"Result input missing: {name}; checked {candidates}")


AUDIT_JSON = result_json("fixed_panel_audit_v2.json", "SCREG_AUDIT_JSON")
ATAC_FILES = {
    "brain": os.environ.get("SCFM_BRAIN_ATAC",
        f"{rfa.DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"),
    "pbmc": os.environ.get("SCREG_PBMC_ATAC", f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"),
}
SEED_ROOT = 20260724            # the frozen audit's seed root (run_fixed_panel_audit.py:29)
SEED_ROOT_REVISION = 20260924   # new analyses in this revision


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_seeds() -> dict:
    """Re-derive the frozen audit's per-test seeds exactly as main() does.

    main(): ss = SeedSequence(SEED_ROOT); brain uses ss.spawn(1)[0], PBMC the next
    spawn; each family spawns 4 ints (mantel full, mantel non-degree, degree full,
    degree non-degree). PBMC scGPT and UCE are appended with their own sequences.
    """
    ss = np.random.SeedSequence(SEED_ROOT)
    ss_brain = ss.spawn(1)[0]
    ss_pbmc = ss.spawn(1)[0]
    out = {}
    for group, seq in (
        ("brain", ss_brain),
        ("pbmc_shared", ss_pbmc),
        ("pbmc_scgpt", np.random.SeedSequence([SEED_ROOT, 20260726, 1])),
        ("pbmc_uce", np.random.SeedSequence([SEED_ROOT, 20260729, 1])),
    ):
        m_full, m_nd, d_full, d_nd = rfa.spawn_int_seeds(seq, 4)
        out[group] = {("mantel", "full"): m_full, ("mantel", "non_degree"): m_nd,
                      ("degree", "full"): d_full, ("degree", "non_degree"): d_nd}
    return out


def _edges(tissue: str, tf_rows: np.ndarray, Ng: int):
    genes = json.loads(Path(fpa.MANI).read_text())["genes"]
    ii_all = np.repeat(tf_rows, Ng)
    jj_all = np.tile(np.arange(Ng), len(tf_rows))
    m = fpa.edge_mask(tissue, genes, tf_rows, ii_all, jj_all)
    return ii_all[m], jj_all[m]


def load_groups() -> dict:
    """Return the four shared-randomization groups of the frozen pooled audit.

    Each group carries its own co-expression control and edge set; rows inside a
    group share one proxy randomization per replicate (as in the frozen audit).
    """
    groups = {}
    G_b, co_b, models_b, tf_b, _ = rfa.load_pooled_brain()
    ko, _ = rfa.load_geneformer_ko()
    entries_b = list(models_b.items()) + [
        ("geneformer_ko_raw", ko["geneformer_ko_raw"]),
        ("geneformer_ko_posctrl", ko["geneformer_ko_posctrl"]),
    ]
    G_p, co_p, models_p, tf_p, _ = rfa.load_pooled_pbmc()
    scgpt = rfa.load_optional_pbmc_scgpt()
    uce = rfa.load_optional_pbmc_uce()
    if scgpt is None or uce is None:
        raise FileNotFoundError("PBMC scGPT/UCE caches are required for the revision analyses")
    spec = [
        ("brain", "brain", G_b, co_b, tf_b, entries_b),
        ("pbmc_shared", "pbmc", G_p, co_p, tf_p, list(models_p.items())),
        ("pbmc_scgpt", "pbmc", G_p, scgpt[1], tf_p, [("scGPT_encoder", scgpt[2])]),
        ("pbmc_uce", "pbmc", G_p, uce[1], tf_p, [("UCE_encoder", uce[2])]),
    ]
    conf_cache = {}
    for group, tissue, G, co, tf, entries in spec:
        if tissue not in conf_cache:
            conf_cache[tissue] = fpa.build_confounds(ATAC_FILES[tissue])
        peakcount, genelen, gc, detv = conf_cache[tissue]
        Ng = G.shape[0]
        ii, jj = _edges(tissue, tf, Ng)
        groups[group] = {
            "tissue": tissue,
            "G": G, "tf_rows": tf, "ii": ii, "jj": jj,
            "co_v": co[ii, jj], "atac_v": G[ii, jj],
            "labels": [lab for lab, _ in entries],
            "fm_vecs": [np.asarray(M)[ii, jj] for _, M in entries],
            "peakcount": peakcount, "genelen": genelen, "gc": gc, "detv": detv,
            "tf_outdeg": (G > 0).sum(1).astype(np.float32),
            "atac_indeg": (G > 0).sum(0).astype(np.float32),
        }
    return groups


def frozen_rows() -> dict:
    """(tissue, spec, model_label) -> frozen audit row, from fixed_panel_audit_v2.json."""
    doc = json.loads(AUDIT_JSON.read_text())
    out = {}
    for tissue in ("brain", "pbmc"):
        for fam in ("primary_family", "sensitivity_family"):
            for r in doc["pooled"][tissue][fam]["rows"]:
                out[(tissue, r["confound_spec"], r["model_label"])] = r
    return out


def observed_unrounded(groups: dict) -> dict:
    """(tissue, spec, model_label) -> unrounded partial rho, recomputed from the audit inputs.

    The audit JSON stores observed_partial_rho rounded to 6 decimals, but its p-values were
    computed from the unrounded statistic; for rows with |rho| near 0 the rounding moves a few
    null draws across |obs|. Every recomputed value must round to the stored one.
    """
    frozen = frozen_rows()
    out = {}
    for g in groups.values():
        for lab, fm_v in zip(g["labels"], g["fm_vecs"]):
            for spec in ("full", "non_degree"):
                rho = fpa.partial_rho_obs_sliced(
                    fm_v=fm_v, atac_v=g["atac_v"], co_v=g["co_v"], jj=g["jj"], ii=g["ii"],
                    peakcount=g["peakcount"], genelen=g["genelen"], detv=g["detv"], gc=g["gc"],
                    tf_outdeg=g["tf_outdeg"], atac_indeg=g["atac_indeg"],
                    use_coexp=True, confound_spec=spec)
                stored = frozen[(g["tissue"], spec, lab)]["observed_partial_rho"]
                assert round(float(rho), 6) == stored, (g["tissue"], spec, lab, rho, stored)
                out[(g["tissue"], spec, lab)] = float(rho)
    return out


def plus_one_p(null: np.ndarray, observed: float) -> float:
    return (int(np.sum(np.abs(null) >= abs(observed))) + 1) / (len(null) + 1)


def write_json(path: Path, doc: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(doc, indent=1, sort_keys=False) + "\n")
    os.replace(tmp, path)
