#!/usr/bin/env python
"""
scfm-reg-audit v2 — fixed-TF-panel audit (deterministic, no superpopulation inference).

Design contract:
  - The preregistered 446 TF x 1200 target panel is treated as fixed.
  - Confound-controlled partial Spearman rho is the point effect.
  - Uncertainty is a finite Monte-Carlo randomization p-value (p_mc) from BOTH a gene-label
    Mantel permutation AND a degree-preserving row-shuffle null, computed on the same panel.
  - No bootstrap CI; no MDE; no coverage; no exclusion; no implied-alpha; no population
    inference language.
  - All RNG derives from one explicit SEED_ROOT via numpy.random.SeedSequence (driver).
  - Tissues apply the same confounds and the same marker/exclusion mask symmetrically.
  - Signal injection is an axis-aligned pipeline-sensitivity diagnostic, never an
    inferential claim: alpha = fraction of the regulatory_potential_proxy residual
    injected into a synthetic graph built from fresh independent N(0,1) noise per
    (tissue, alpha, replicate).

Null semantics — gene-label Mantel randomization, two confound specs:

  confound_spec = 'full':
    Design matrix columns = [1, rankdata(co), z(peakcount), z(genelen), z(detection),
                              z(GC), z(tf_outdeg_perm), z(atac_indeg_perm)]
    Under each gene-label permutation, ONLY:
      - atac_perm = G_atac[perm[ii], perm[jj]]
      - tf_outdeg_perm = tf_outdeg[perm]
      - atac_indeg_perm = atac_indeg[perm]
    are recomputed. rankdata(fm) is FIXED (the FM labels are not permuted). rankdata(co)
    and the four gene-identity confounds are FIXED (they are properties of gene identity,
    not of G_atac structure). p_mc = (count(|null| >= |observed|) + 1) / (N_perm + 1).

  confound_spec = 'non_degree':
    Design matrix columns = [1, rankdata(co), z(peakcount), z(genelen), z(detection),
                              z(GC)]
    NO degree columns. The null randomizes atac_perm only; no degree columns are
    recomputed under perm (they are absent from the design). The OBSERVED statistic
    for this family is also non-degree; the null statistic matches the observed
    statistic exactly (same column set).

Degree-preserving row-shuffle null:
    For each TF row independently, permute its edge-weight values across the non-self
    targets. Per-row out-degree is preserved by construction; indegree changes. The
    statistic is the same column set as 'full' above. Tests whether the SPECIFIC target
    assignment carries signal beyond aggregate magnitude.

Pure numpy/scipy on cached graphs (results/v2/*.npz). No GPU. No FM re-embedding.
"""
import hashlib
import json
import os
import warnings
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pyfaidx
from scipy.stats import rankdata


# ----------------------------- paths / constants -----------------------------
def _project_root():
    """Resolve the project/capsule root: env override, then the nearest ancestor of
    this file that contains data/manifest/shared_genes.v2.json (works from both the
    development layout src/v2/ and the release capsule layout src/)."""
    env = os.environ.get("SCREG_PROJECT_ROOT")
    if env:
        return os.path.abspath(env)
    here = os.path.abspath(os.path.dirname(__file__))
    for level in range(1, 6):
        candidate = os.path.abspath(os.path.join(here, *[".."] * level))
        if os.path.exists(os.path.join(candidate, "data", "manifest", "shared_genes.v2.json")):
            return candidate
    return os.path.abspath(os.path.join(here, "..", ".."))


ROOT = _project_root()
OUT = f"{ROOT}/results/v2"
DATA = f"{ROOT}/data"
MANI = f"{DATA}/manifest/shared_genes.v2.json"
COORDS = f"{DATA}/annotation/gene_coords_hg38.tsv"
HG38 = f"{DATA}/genome/hg38.fa"
PROM = 2000
W = 500

MARKER_GENES = {
    "MOBP", "MOG", "PLP1", "MBP", "ST18", "CTNNA3", "SLC17A7", "SATB2", "RBFOX3",
    "NRGN", "SLC17A6", "GAD1", "GAD2", "SLC32A1", "DLX1", "AQP4", "GFAP", "SLC1A2",
    "ALDH1L1", "CSF1R", "P2RY12", "CX3CR1", "C1QA", "C1QB", "PDGFRA", "CSPG4",
    "OLIG1", "OLIG2", "VCAN", "CLDN5", "FLT1", "PECAM1", "PDGFRB", "RGS5",
}
MARKER_TISSUES = {"brain"}


# ----------------------------- manifest / confound cache ---------------------
def load_manifest() -> Tuple[List[str], Dict[str, float], str]:
    from pathlib import Path
    man = json.loads(Path(MANI).read_text())
    genes = man["genes"]
    sha = hashlib.sha256(("\n".join(genes)).encode()).hexdigest()
    if sha != man["sha256"]:
        raise ValueError(f"manifest sha mismatch for {MANI}: {sha} vs {man['sha256']}")
    return genes, man["detection"], sha


_CONF_CACHE: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}


def build_confounds(atac_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if atac_file in _CONF_CACHE:
        return _CONF_CACHE[atac_file]
    genes, det, _ = load_manifest()
    gidx = {g: i for i, g in enumerate(genes)}
    Ng = len(genes)
    gco: Dict[str, Tuple[str, int, int, int]] = {}
    from pathlib import Path
    with Path(COORDS).open() as _coords_fh:
        for ln in _coords_fh:
            c, s, e, st, nm = ln.rstrip("\n").split("\t")
            if nm not in gidx or nm in gco:
                continue
            s, e = int(s), int(e)
            lo = s - PROM if st == "+" else s
            hi = e if st == "+" else e + PROM
            gco[nm] = (c, lo, hi, abs(e - s))
    Av = ad.read_h5ad(atac_file, backed="r")
    peaks = [str(p) for p in Av.var_names]
    pchr = np.array([p.split(":")[0] for p in peaks])
    pse = np.array([[int(x) for x in p.split(":")[1].split("-")] for p in peaks])
    pmid = (pse[:, 0] + pse[:, 1]) // 2
    by: Dict[str, np.ndarray] = {}
    for i, c in enumerate(pchr):
        by.setdefault(c, []).append(i)
    for c in list(by.keys()):
        by[c] = np.array(by[c])
    fa = pyfaidx.Fasta(HG38, sequence_always_upper=True)
    peakcount = np.zeros(Ng); genelen = np.zeros(Ng); gc = np.zeros(Ng)
    without_coords: List[str] = []
    without_peak_chromosome: List[str] = []
    without_genome_chromosome: List[str] = []
    for g, i in gidx.items():
        if g not in gco:
            without_coords.append(g)
            continue
        c, lo, hi, ln_ = gco[g]
        genelen[i] = ln_
        pis = by.get(c)
        if pis is None:
            without_peak_chromosome.append(g)
            continue
        sel = pis[(pmid[pis] >= lo) & (pmid[pis] <= hi)]
        peakcount[i] = len(sel)
        if not len(sel):
            continue
        if c not in fa.keys():
            without_genome_chromosome.append(g)
            continue
        ks = []
        for p in sel[:40]:
            s = str(fa[c][max(0, int(pmid[p]) - W // 2): int(pmid[p]) + W // 2])
            ks.append((s.count("G") + s.count("C")) / W)
        gc[i] = float(np.mean(ks))
    _warn_zero_filled_confounds(atac_file, without_coords, without_peak_chromosome,
                               without_genome_chromosome)
    missing_detection = [g for g in genes if g not in det]
    if missing_detection:
        warnings.warn(
            f"{len(missing_detection)} manifest genes have no detection rate in {MANI} and are "
            f"zero-filled (first: {missing_detection[:5]})",
            RuntimeWarning, stacklevel=2,
        )
    detv = np.array([det.get(g, 0.0) for g in genes])
    out = (peakcount.astype(np.float32), genelen.astype(np.float32),
           gc.astype(np.float32), detv.astype(np.float32))
    _CONF_CACHE[atac_file] = out
    return out


def _warn_zero_filled_confounds(atac_file: str, without_coords: List[str],
                                without_peak_chromosome: List[str],
                                without_genome_chromosome: List[str]) -> None:
    """Surface silently zero-filled confound columns.

    A gene with no coordinate record, no peak on its chromosome, or no genome
    contig keeps a zero peakcount/genelen/GC entry. That degradation is a data
    problem, not a statistic, so report it instead of letting the design matrix
    absorb it unnoticed.
    """
    for genes_missing, reason in (
        (without_coords, f"no coordinate record in {COORDS}"),
        (without_peak_chromosome, f"no ATAC peak on their chromosome in {atac_file}"),
        (without_genome_chromosome, f"no matching contig in {HG38}"),
    ):
        if genes_missing:
            warnings.warn(
                f"{len(genes_missing)} manifest genes have {reason}; their confound columns "
                f"are zero-filled (first: {sorted(genes_missing)[:5]})",
                RuntimeWarning, stacklevel=3,
            )


# ----------------------------- core statistics -------------------------------
def zscore(v: np.ndarray) -> np.ndarray:
    v = v.astype(float)
    return (v - v.mean()) / (v.std() + 1e-9)


def resid(y: np.ndarray, C: np.ndarray) -> np.ndarray:
    """OLS residual y - [1|C] @ beta via np.linalg.solve on C1.T @ C1."""
    C1 = np.column_stack([np.ones(len(y)), C])
    CtC = C1.T @ C1
    try:
        beta = np.linalg.solve(CtC, C1.T @ y)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(C1, y, rcond=None)[0]
    return y - C1 @ beta


def pcorr(x: np.ndarray, y: np.ndarray, controls: np.ndarray) -> float:
    """Reference partial correlation.

    ``controls`` must not contain an intercept; this function adds exactly one.
    Call :func:`_pcorr_full_design` only when the design already contains its single
    intercept column.
    """
    if controls.ndim == 1:
        controls = controls.reshape(-1, 1)
    design = np.column_stack([np.ones(len(x)), controls])
    return _pcorr_full_design(x, y, design)


def _stable_solve(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return the least-squares solution without forming ill-conditioned normal equations."""
    try:
        return np.linalg.lstsq(X, y, rcond=None)[0]
    except np.linalg.LinAlgError:
        return np.linalg.pinv(X) @ y


def _column_space_basis(X: np.ndarray) -> np.ndarray:
    """Return an orthonormal basis for the numerical column space of ``X``."""
    U, singular_values, _ = np.linalg.svd(X, full_matrices=False)
    if singular_values.size == 0:
        return np.empty((X.shape[0], 0), dtype=np.float64)
    tolerance = np.finfo(singular_values.dtype).eps * max(X.shape) * singular_values[0]
    return U[:, singular_values > tolerance]


def _rho_from_moments(num: float, den: float, context: str) -> float:
    """Correlation from precomputed cross- and self-moments of centred residuals.

    ``den == 0`` means a genuinely constant residual, for which zero correlation is
    the defined answer. A non-finite moment instead means NaN/inf entered the inputs;
    reporting that as a zero effect would silently bias every downstream p_mc, so it
    is raised.
    """
    if not np.isfinite(num) or not np.isfinite(den):
        raise ValueError(
            f"non-finite partial-correlation moments in {context}: num={num}, den={den}; "
            "check the input graphs and confounds for NaN/inf")
    return num / den if den > 0 else 0.0


def _pcorr_full_design(x: np.ndarray, y: np.ndarray, design: np.ndarray) -> float:
    """Low-level partial correlation for a design containing exactly one intercept."""
    bx = _stable_solve(design, x)
    by = _stable_solve(design, y)
    rx = x - design @ bx
    ry = y - design @ by
    rxn = rx - rx.mean()
    ryn = ry - ry.mean()
    num = float(rxn @ ryn)
    den = float(np.sqrt((rxn @ rxn) * (ryn @ ryn)))
    return _rho_from_moments(num, den, "pcorr")


def pcorr_inplace(x: np.ndarray, y: np.ndarray, C_full: np.ndarray) -> float:
    """Compatibility wrapper for a full design that already contains one intercept."""
    return _pcorr_full_design(x, y, C_full)


def pcorr_fwl(x_r: np.ndarray, y_r: np.ndarray, X_fixed_resid_x: np.ndarray,
              Q: np.ndarray, X_changing: np.ndarray) -> float:
    """Frisch-Waugh-Lovell (FWL) joint OLS for the full design [1|C_fixed|X_changing]
    without forming CtC of the full matrix.

    FWL theorem: the partial residual of y on the full design equals the residual of
    (y - Q Q^T y) on X_changing_orth, where Q is the orthonormal factor of X_fixed and
    X_changing_orth = X_changing - Q Q^T X_changing.

    `X_fixed_resid_x` is the fm_r residual against X_fixed only, precomputed ONCE per row
    (fm is fixed across perms). Falls back to lstsq/pinv on singular 2-col XtX."""
    rx_fixed = X_fixed_resid_x
    ry_fixed = y_r - Q @ (Q.T @ y_r)
    X_ch_orth = X_changing - Q @ (Q.T @ X_changing)
    bx = _stable_solve(X_ch_orth, rx_fixed)
    by = _stable_solve(X_ch_orth, ry_fixed)
    rx = rx_fixed - X_ch_orth @ bx
    ry = ry_fixed - X_ch_orth @ by
    rxn = rx - rx.mean()
    ryn = ry - ry.mean()
    num = float(rxn @ ryn)
    den = float(np.sqrt((rxn @ rxn) * (ryn @ ryn)))
    return _rho_from_moments(num, den, "pcorr_fwl")


def spearman_paired(x: np.ndarray, y: np.ndarray) -> float:
    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        raise ValueError("spearman_paired needs varying ranks on both sides; got a constant input")
    rho = float(np.corrcoef(rx, ry)[0, 1])
    if not np.isfinite(rho):
        raise ValueError(f"spearman_paired produced a non-finite rho ({rho}); inputs carry NaN/inf")
    return rho


# ----------------------------- mask helpers ----------------------------------
def edge_mask(tissue: str, genes: List[str], tf_rows: np.ndarray,
              ii: np.ndarray, jj: np.ndarray) -> np.ndarray:
    """Symmetric marker-gene mask: drops edges involving a marker gene for tissues in
    MARKER_TISSUES. Same policy applied identically to every (tissue, cell_type) row."""
    m = ii != jj
    if tissue in MARKER_TISSUES:
        mgi = np.array([g not in MARKER_GENES for g in genes])
        m = m & mgi[ii] & mgi[jj]
    return m


# ----------------------------- confound-controlled partial rho ---------------
def partial_rho_obs_sliced(
    fm_v: np.ndarray, atac_v: np.ndarray, co_v: np.ndarray,
    jj: np.ndarray, ii: np.ndarray,
    peakcount: np.ndarray, genelen: np.ndarray, detv: np.ndarray, gc: np.ndarray,
    tf_outdeg: np.ndarray, atac_indeg: np.ndarray,
    use_coexp: bool, confound_spec: str,
) -> float:
    """Observed partial Spearman rho under the chosen confound spec. confound_spec in
    {'non_degree', 'full'}. For non_degree, NO degree columns are used; the statistic
    matches the null's statistic exactly."""
    pc_z = zscore(peakcount[jj])
    gl_z = zscore(genelen[jj])
    dv_z = zscore(detv[jj])
    gc_z = zscore(gc[jj])
    conf_cols = [pc_z, gl_z, dv_z, gc_z]
    if confound_spec == "full":
        conf_cols += [zscore(tf_outdeg[ii]), zscore(atac_indeg[jj])]
    C = np.column_stack(conf_cols)
    C_use = np.column_stack([rankdata(co_v), C]) if use_coexp else C
    return pcorr(rankdata(fm_v), rankdata(atac_v), C_use)


# ----------------------------- gene-label Mantel randomization --------------
def mantel_randomization(
    fm_v: np.ndarray, atac_v: np.ndarray, co_v: np.ndarray,
    jj: np.ndarray, ii: np.ndarray,
    peakcount: np.ndarray, genelen: np.ndarray, detv: np.ndarray, gc: np.ndarray,
    tf_outdeg_full: np.ndarray, atac_indeg_full: np.ndarray,
    G_atac_full: np.ndarray,
    use_coexp: bool, confound_spec: str,
    observed: float, n_perm: int, seed: int,
) -> Dict[str, float]:
    """Gene-label Mantel randomization.

    NULL SEMANTICS:
      - fm_r = rankdata(fm_v) is FIXED across perms (FM labels are not randomized).
      - For 'full' spec, the 2 degree columns (zscore(tf_outdeg_perm[ii]),
        zscore(atac_indeg_perm[jj])) are RECOMPUTED under perm because they are
        G_atac-derived.
      - For 'non_degree' spec, NO degree columns are present; the design matrix is
        the same as the observed statistic.
      - co_r = rankdata(co_v) and the 4 non-degree gene-identity confounds
        (peakcount, genelen, detection, GC) are FIXED under perm.
    p_mc = (count(|null| >= |observed|) + 1) / (n_perm + 1).
    """
    Ng = G_atac_full.shape[0]
    N = len(ii)
    fm_r = rankdata(fm_v)
    co_r = rankdata(co_v)
    pc_z = zscore(peakcount[jj])
    gl_z = zscore(genelen[jj])
    dv_z = zscore(detv[jj])
    gc_z = zscore(gc[jj])

    n_fixed = 6 if use_coexp else 5
    n_changing = 2 if confound_spec == "full" else 0
    n_cols = n_fixed + n_changing
    C_full = np.empty((N, n_cols), dtype=np.float64)
    if use_coexp:
        C_full[:, 0] = 1.0
        C_full[:, 1] = co_r
        C_full[:, 2] = pc_z
        C_full[:, 3] = gl_z
        C_full[:, 4] = dv_z
        C_full[:, 5] = gc_z
    else:
        C_full[:, 0] = 1.0
        C_full[:, 1] = pc_z
        C_full[:, 2] = gl_z
        C_full[:, 3] = dv_z
        C_full[:, 4] = gc_z

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        perm = rng.permutation(Ng)
        atac_perm = G_atac_full[perm[ii], perm[jj]]
        ap_r = rankdata(atac_perm)
        if confound_spec == "full":
            outdeg_perm = tf_outdeg_full[perm]
            indeg_perm = atac_indeg_full[perm]
            C_full[:, n_fixed] = zscore(outdeg_perm[ii])
            C_full[:, n_fixed + 1] = zscore(indeg_perm[jj])
        null[k] = pcorr_inplace(fm_r, ap_r, C_full)
    abs_null = np.abs(null)
    p_mc = (int(np.sum(abs_null >= abs(observed))) + 1) / (n_perm + 1)
    return {
        "p_mc": float(p_mc),
        "N_perm": int(n_perm),
        "seed": int(seed),
        "resolution": float(1 / (n_perm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "z": float((observed - null.mean()) / (null.std() + 1e-9)),
        "null_obs_count_at_or_above_obs": int(np.sum(abs_null >= abs(observed))),
        "test_type": "gene_label_mantel_plus_one_corrected",
        "confound_spec": confound_spec,
        "null_columns_perm_recomputed": (["tf_outdeg", "atac_indeg"] if confound_spec == "full" else ["atac_only"]),
        "null_columns_fixed_under_perm": ["fm", "co", "peakcount", "genelen", "detection", "GC"],
    }


# ----------------------------- degree-preserving row-shuffle null ------------
def degree_preserving_null(
    fm_v: np.ndarray, atac_v: np.ndarray, co_v: np.ndarray,
    jj: np.ndarray, ii: np.ndarray,
    peakcount: np.ndarray, genelen: np.ndarray, detv: np.ndarray, gc: np.ndarray,
    tf_outdeg_full: np.ndarray, atac_indeg_full: np.ndarray,
    G_atac_full: np.ndarray, tf_rows_unique: np.ndarray,
    use_coexp: bool, confound_spec: str,
    observed: float, n_perm: int, seed: int,
) -> Dict[str, float]:
    """Row-shuffle null: per TF row, permute edge weights across non-self targets.
    The design uses the same ``confound_spec`` as the observed statistic."""
    Ng = G_atac_full.shape[0]
    N = len(ii)
    fm_r = rankdata(fm_v)
    co_r = rankdata(co_v)
    pc_z = zscore(peakcount[jj])
    gl_z = zscore(genelen[jj])
    dv_z = zscore(detv[jj])
    gc_z = zscore(gc[jj])
    fixed_columns = [co_r, pc_z, gl_z, dv_z, gc_z] if use_coexp else [pc_z, gl_z, dv_z, gc_z]
    if confound_spec == "full":
        fixed_columns.append(zscore(tf_outdeg_full[ii]))
    C_full = np.column_stack([np.ones(N), *fixed_columns])

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    tf_list = np.unique(tf_rows_unique)
    non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, Ng)]) for t in tf_list}
    # Cache row hashes / non-self indices once (used inside hot loop)
    Gp = G_atac_full.copy()
    for k in range(n_perm):
        for t in tf_list:
            nz = non_self[t]
            vals = Gp[t, nz].copy()
            rng.shuffle(vals)
            Gp[t, nz] = vals
        atac_p = Gp[ii, jj]
        ap_r = rankdata(atac_p)
        if confound_spec == "full":
            indeg_p = (Gp > 0).sum(0).astype(np.float32)
            C_perm = np.column_stack([C_full, zscore(indeg_p[jj])])
        else:
            C_perm = C_full
        null[k] = pcorr_inplace(fm_r, ap_r, C_perm)
    abs_null = np.abs(null)
    p_mc = (int(np.sum(abs_null >= abs(observed))) + 1) / (n_perm + 1)
    return {
        "p_mc": float(p_mc),
        "N_perm": int(n_perm),
        "seed": int(seed),
        "resolution": float(1 / (n_perm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "z": float((observed - null.mean()) / (null.std() + 1e-9)),
        "null_obs_count_at_or_above_obs": int(np.sum(abs_null >= abs(observed))),
        "test_type": "degree_preserving_row_shuffle_plus_one_corrected",
        "confound_spec": confound_spec,
        "null_columns_perm_recomputed": (["atac", "atac_indeg"]
                                           if confound_spec == "full" else ["atac"]),
        "null_columns_fixed_under_perm": (["fm", "co", "peakcount", "genelen", "detection", "GC", "tf_outdeg"]
                                            if confound_spec == "full" else
                                            ["fm", "co", "peakcount", "genelen", "detection", "GC"]),
    }


# ----------------------------- batched shared-null helpers -------------------
def batched_mantel_null(
    fm_vecs: List[np.ndarray],
    co_v: np.ndarray,
    jj: np.ndarray,
    ii: np.ndarray,
    peakcount: np.ndarray,
    genelen: np.ndarray,
    detv: np.ndarray,
    gc: np.ndarray,
    tf_outdeg_full: np.ndarray,
    atac_indeg_full: np.ndarray,
    G_atac_full: np.ndarray,
    use_coexp: bool,
    confound_spec: str,
    n_perm: int,
    seed: int,
    perms: Optional[np.ndarray] = None,
) -> Tuple[List[np.ndarray], Dict[str, object]]:
    """Batched gene-label Mantel randomization: one shared proxy null distribution per
    replicate, paired against every FM vector.

    `perms` (optional, mutually exclusive with `seed`): explicit (n_perm, Ng) array of
    gene-label permutations; when provided, used as-is (no RNG). When None, perms are
    generated from a SeedSequence via replicate-spawning uint64 seeds.

    Returns (per-row null arrays) + shared batch metadata."""
    if perms is not None and seed is not None:
        raise ValueError("perms and seed are mutually exclusive")
    if perms is not None and (perms.ndim != 2 or perms.shape[0] != n_perm or perms.shape[1] != G_atac_full.shape[0]):
        raise ValueError(f"perms must have shape (n_perm={n_perm}, Ng={G_atac_full.shape[0]}); got {perms.shape}")

    Ng = G_atac_full.shape[0]
    N = len(ii)
    n_rows = len(fm_vecs)
    fm_r_list = [rankdata(f) for f in fm_vecs]

    pc_z = zscore(peakcount[jj]); gl_z = zscore(genelen[jj])
    dv_z = zscore(detv[jj]); gc_z = zscore(gc[jj])
    co_r = rankdata(co_v)

    if use_coexp:
        X_fixed = np.column_stack([np.ones(N), co_r, pc_z, gl_z, dv_z, gc_z])
    else:
        X_fixed = np.column_stack([np.ones(N), pc_z, gl_z, dv_z, gc_z])
    Q = _column_space_basis(X_fixed)

    X_fixed_resid_per_fm = [f_r - Q @ (Q.T @ f_r) for f_r in fm_r_list]

    nulls: List[np.ndarray] = [np.empty(n_perm, dtype=float) for _ in range(n_rows)]
    if perms is not None:
        # Explicit perms path: derive degenerate replicate_seed metadata
        replicate_seeds: List[int] = []
    else:
        rng_proxy = np.random.SeedSequence(seed)
        replicate_seeds = []

    for k in range(n_perm):
        if perms is not None:
            perm = perms[k]
        else:
            rep_seq = rng_proxy.spawn(1)[0]
            rep_seed = int(rep_seq.generate_state(1, dtype=np.uint64)[0])
            replicate_seeds.append(rep_seed)
            rng = np.random.default_rng(rep_seed)
            perm = rng.permutation(Ng)
        # IMPORTANT: atac is rankdata'd BEFORE residualizing against the rank-scale
        # design. The pcorr statistic is rank-based; mixing raw and rank scales would
        # bias the residualization. (Discovered via the explicit-index equivalence test.)
        atac_perm_k = rankdata(G_atac_full[perm[ii], perm[jj]])
        atac_fixed_resid = atac_perm_k - Q @ (Q.T @ atac_perm_k)
        if confound_spec == "full":
            outdeg_perm = tf_outdeg_full[perm]; indeg_perm = atac_indeg_full[perm]
            X_changing_perm = np.column_stack([zscore(outdeg_perm[ii]),
                                                zscore(indeg_perm[jj])])
        else:
            X_changing_perm = np.empty((N, 0), dtype=np.float64)
        if X_changing_perm.shape[1] > 0:
            X_ch_orth = X_changing_perm - Q @ (Q.T @ X_changing_perm)
            beta_atac = _stable_solve(X_ch_orth, atac_fixed_resid)
            atac_full_resid = atac_fixed_resid - X_ch_orth @ beta_atac
        else:
            atac_full_resid = atac_fixed_resid
        atac_centered = atac_full_resid - atac_full_resid.mean()
        atac_norm_sq = float(atac_centered @ atac_centered)
        for r, X_fixed_resid_x in enumerate(X_fixed_resid_per_fm):
            if X_changing_perm.shape[1] > 0:
                beta_fm = _stable_solve(X_ch_orth, X_fixed_resid_x)
                fm_full_resid = X_fixed_resid_x - X_ch_orth @ beta_fm
            else:
                fm_full_resid = X_fixed_resid_x
            fm_centered = fm_full_resid - fm_full_resid.mean()
            fm_norm_sq = float(fm_centered @ fm_centered)
            num = float(fm_centered @ atac_centered)
            den = float(np.sqrt(fm_norm_sq * atac_norm_sq))
            nulls[r][k] = _rho_from_moments(num, den, f"batched_mantel_null row {r} perm {k}")

    meta = {
        "shared_seed": int(seed) if seed is not None else None,
        "replicate_seed_stream": "explicit perms" if perms is not None
            else "uint64 per perm via numpy.random.SeedSequence.generate_state(1)",
        "replicate_seeds": replicate_seeds,
        "X_fixed_columns": (["1", "rankdata(co)", "z(peakcount)", "z(genelen)", "z(detection)", "z(GC)"]
                            if use_coexp else
                            ["1", "z(peakcount)", "z(genelen)", "z(detection)", "z(GC)"]),
        "X_changing_columns": (["z(tf_outdeg_perm)", "z(atac_indeg_perm)"]
                                if confound_spec == "full" else []),
        "null_columns_fixed_under_perm": ["fm_r (rankdata is fixed)", "co_r (rankdata is fixed)",
                                            "peakcount", "genelen", "detection", "GC"],
        "null_columns_perm_recomputed": (["atac", "tf_outdeg", "atac_indeg"]
                                          if confound_spec == "full" else ["atac"]),
        "confound_spec": confound_spec,
        "n_perm": int(n_perm),
        "batch_id": f"mantel_seed{seed}_spec{confound_spec}_n{n_perm}",
        "shared_proxy_null_per_replicate": True,
        "n_rows_in_batch": n_rows,
        "resolution": float(1 / (n_perm + 1)),
    }
    return nulls, meta


def batched_degree_preserving_null(
    fm_vecs: List[np.ndarray],
    co_v: np.ndarray,
    jj: np.ndarray,
    ii: np.ndarray,
    peakcount: np.ndarray,
    genelen: np.ndarray,
    detv: np.ndarray,
    gc: np.ndarray,
    tf_outdeg_full: np.ndarray,
    atac_indeg_full: np.ndarray,
    G_atac_full: np.ndarray,
    tf_rows_unique: np.ndarray,
    use_coexp: bool,
    confound_spec: str,
    n_perm: int,
    seed: int,
    precomputed_shuffled_graphs: Optional[List[np.ndarray]] = None,
) -> Tuple[List[np.ndarray], Dict[str, object]]:
    """Batched degree-preserving row-shuffle null: one shared row-shuffled graph per
    replicate, paired against every FM vector.

    `precomputed_shuffled_graphs` (optional, mutually exclusive with `seed`): explicit
    list of (Ng, Ng) row-shuffled graphs (one per replicate) the caller wants used
    verbatim. When None, graphs are generated from `seed` via replicate-spawning.

    IMPORTANT (revised semantics, per blocker 1): in the FULL design, the changing
    columns are the perm-recomputed TARGET indegree ONLY. TF outdegree is fixed
    under the row shuffle (per-row outdegree is preserved by construction) and is
    already part of the FIXED design, so it is NOT recomputed per perm. The
    previous version double-counted tf_outdeg in the changing block, which biased
    the statistic. The non-degree specification still uses the row-shuffle null,
    but omits both degree columns."""
    if precomputed_shuffled_graphs is not None and seed is not None:
        raise ValueError("precomputed_shuffled_graphs and seed are mutually exclusive")
    Ng = G_atac_full.shape[0]
    N = len(ii)
    n_rows = len(fm_vecs)
    fm_r_list = [rankdata(f) for f in fm_vecs]
    pc_z = zscore(peakcount[jj]); gl_z = zscore(genelen[jj])
    dv_z = zscore(detv[jj]); gc_z = zscore(gc[jj])
    co_r = rankdata(co_v)
    od_z_ii = zscore(tf_outdeg_full[ii])  # row-shuffle-invariant

    # Non-degree and full specifications share the same row-shuffled proxy null.
    # The full specification additionally controls invariant TF out-degree and
    # perm-recomputed target in-degree.
    if confound_spec == "full":
        fixed_columns = [co_r, pc_z, gl_z, dv_z, gc_z, od_z_ii] if use_coexp else [pc_z, gl_z, dv_z, gc_z, od_z_ii]
    else:
        fixed_columns = [co_r, pc_z, gl_z, dv_z, gc_z] if use_coexp else [pc_z, gl_z, dv_z, gc_z]
    X_fixed = np.column_stack([np.ones(N), *fixed_columns])
    Q = _column_space_basis(X_fixed)
    X_fixed_resid_per_fm = [f_r - Q @ (Q.T @ f_r) for f_r in fm_r_list]

    nulls: List[np.ndarray] = [np.empty(n_perm, dtype=float) for _ in range(n_rows)]
    tf_list = np.unique(tf_rows_unique)
    non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, Ng)]) for t in tf_list}

    if precomputed_shuffled_graphs is not None:
        if len(precomputed_shuffled_graphs) != n_perm:
            raise ValueError(f"precomputed_shuffled_graphs length {len(precomputed_shuffled_graphs)} != n_perm {n_perm}")
        replicate_seeds: List[int] = []
    else:
        rng_proxy = np.random.SeedSequence(seed)
        replicate_seeds = []
    Gp = G_atac_full.copy()
    for k in range(n_perm):
        if precomputed_shuffled_graphs is not None:
            Gp = precomputed_shuffled_graphs[k].copy()
        else:
            rep_seq = rng_proxy.spawn(1)[0]
            rep_seed = int(rep_seq.generate_state(1, dtype=np.uint64)[0])
            replicate_seeds.append(rep_seed)
            rng = np.random.default_rng(rep_seed)
            # Row-shuffle: each TF row's non-self values permuted INDEPENDENTLY (preserves per-row outdegree).
            for t in tf_list:
                nz = non_self[t]
                vals = Gp[t, nz].copy()
                rng.shuffle(vals)
                Gp[t, nz] = vals
        atac_p = Gp[ii, jj]
        indeg_p = (Gp > 0).sum(0).astype(np.float64)
        atac_p_r = rankdata(atac_p)
        atac_fixed_resid = atac_p_r - Q @ (Q.T @ atac_p_r)
        if confound_spec == "full":
            # Only target in-degree changes under a within-row shuffle; TF
            # out-degree is already present in X_fixed and is invariant.
            X_changing_perm = zscore(indeg_p[jj]).reshape(-1, 1)
            X_ch_orth = X_changing_perm - Q @ (Q.T @ X_changing_perm)
            beta_atac = _stable_solve(X_ch_orth, atac_fixed_resid)
            atac_full_resid = atac_fixed_resid - X_ch_orth @ beta_atac
        else:
            X_changing_perm = np.empty((N, 0), dtype=np.float64)
            X_ch_orth = X_changing_perm
            atac_full_resid = atac_fixed_resid
        atac_centered = atac_full_resid - atac_full_resid.mean()
        atac_norm_sq = float(atac_centered @ atac_centered)
        for r, X_fixed_resid_x in enumerate(X_fixed_resid_per_fm):
            if confound_spec == "full":
                beta_fm = _stable_solve(X_ch_orth, X_fixed_resid_x)
                fm_full_resid = X_fixed_resid_x - X_ch_orth @ beta_fm
            else:
                fm_full_resid = X_fixed_resid_x
            fm_centered = fm_full_resid - fm_full_resid.mean()
            fm_norm_sq = float(fm_centered @ fm_centered)
            num = float(fm_centered @ atac_centered)
            den = float(np.sqrt(fm_norm_sq * atac_norm_sq))
            nulls[r][k] = _rho_from_moments(
                num, den, f"batched_degree_preserving_null row {r} perm {k}")

    meta = {
        "shared_seed": int(seed) if seed is not None else None,
        "replicate_seed_stream": "explicit shuffled graphs" if precomputed_shuffled_graphs is not None
            else "uint64 per perm via numpy.random.SeedSequence.generate_state(1)",
        "replicate_seeds": replicate_seeds,
        "X_fixed_columns": ((["1", "rankdata(co)", "z(peakcount)", "z(genelen)", "z(detection)", "z(GC)", "z(tf_outdeg)"]
                             if use_coexp else
                             ["1", "z(peakcount)", "z(genelen)", "z(detection)", "z(GC)", "z(tf_outdeg)"])
                            if confound_spec == "full" else
                            (["1", "rankdata(co)", "z(peakcount)", "z(genelen)", "z(detection)", "z(GC)"]
                             if use_coexp else
                             ["1", "z(peakcount)", "z(genelen)", "z(detection)", "z(GC)"])),
        "X_changing_columns": (["z(atac_indeg)"] if confound_spec == "full" else []),
        "null_columns_fixed_under_perm": (["fm_r", "co_r", "peakcount", "genelen", "detection", "GC", "tf_outdeg"]
                                            if confound_spec == "full" else
                                            ["fm_r", "co_r", "peakcount", "genelen", "detection", "GC"]),
        "null_columns_perm_recomputed": (["atac", "atac_indeg"] if confound_spec == "full" else ["atac"]),
        "confound_spec": confound_spec,
        "n_perm": int(n_perm),
        "batch_id": (f"degree_precomputed_spec{confound_spec}_n{n_perm}" if precomputed_shuffled_graphs is not None
                     else f"degree_seed{seed}_spec{confound_spec}_n{n_perm}"),
        "shared_proxy_null_per_replicate": True,
        "n_rows_in_batch": n_rows,
        "resolution": float(1 / (n_perm + 1)),
    }
    return nulls, meta


def batched_pvalue_summary(null: np.ndarray, observed: float, n_perm: int,
                            seed: int, batch_id: str,
                            test_type: Optional[str] = None,
                            confound_spec: Optional[str] = None) -> Dict[str, float]:
    """Plus-one Monte-Carlo p-value summary for a single null distribution.
    Optionally stamps test_type and confound_spec into the summary so callers don't
    need to overlay those fields separately."""
    abs_null = np.abs(null)
    p_mc = (int(np.sum(abs_null >= abs(observed))) + 1) / (n_perm + 1)
    out = {
        "p_mc": float(p_mc),
        "N_perm": int(n_perm),
        "seed": int(seed),
        "resolution": float(1 / (n_perm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "z": float((observed - null.mean()) / (null.std() + 1e-9)),
        "null_obs_count_at_or_above_obs": int(np.sum(abs_null >= abs(observed))),
        "batch_id": batch_id,
    }
    if test_type is not None:
        out["test_type"] = test_type
    if confound_spec is not None:
        out["confound_spec"] = confound_spec
    return out


# ----------------------------- BH-FDR ---------------------------------------
def bh_qvalues(pvals: List[float]) -> List[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return [round(float(x), 6) for x in out]


# ----------------------------- output writing --------------------------------
def write_json_atomic(path: str, document: dict, indent: int = 1) -> None:
    """Write an artifact JSON via a staged temp file.

    Serialising into a temp file first means a failure part-way through (non-finite
    value, full disk, interrupt) leaves the previous artifact intact instead of a
    truncated file that later steps would read as authoritative.
    """
    import tempfile
    from pathlib import Path
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp",
                                     dir=destination.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(document, handle, indent=indent, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, destination)
    except Exception:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise


# ----------------------------- file hashing ---------------------------------
def sha256_file(p: str) -> str:
    from pathlib import Path
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sha256_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


# ----------------------------- provenance record builder ---------------------
def matrix_provenance(path: str, key: Optional[str] = None) -> Dict[str, object]:
    """Build a provenance record for a (path, optional key) matrix."""
    rec: Dict[str, object] = {"file_sha256": sha256_file(path)}
    if key is not None:
        z = np.load(path, allow_pickle=False)
        arr = z[key]
        rec["key"] = key
        rec["matrix_sha256"] = sha256_array(arr)
        rec["shape"] = list(arr.shape)
        rec["dtype"] = str(arr.dtype)
    return rec