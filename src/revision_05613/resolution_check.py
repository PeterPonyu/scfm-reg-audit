#!/usr/bin/env python
"""R1.1 / R2.4 resolution check: every frozen Monte Carlo test rerun at N = 9,999.

Covers the 13 pooled FM rows (both nulls, both confound specifications), the two
same-edge co-expression baselines, and the TF-disjoint probe (family permutation
and paired sign-flip). Every test reuses its frozen seed. In all three seed
contracts the replicate stream is consumed in order, so the first 999 replicates
of the deeper run are the frozen replicates; the script asserts that the p-values
recomputed from that prefix equal the frozen values before it reports anything.

The frozen N = 999 audit stays the registered primary analysis. This file only
adds resolution. Null arrays are saved so the joint-null (Westfall-Young) and
MDE analyses can reuse them.

Usage:
  python resolution_check.py --n-perm 9999 --workers 8
  python resolution_check.py --n-perm 999 --only pbmc_uce:mantel:full   # prefix pilot
"""
import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import fpa, log  # noqa: E402

BASELINE_SEEDS = {  # src/brain_coexp_baseline_null.py, src/pbmc_coexp_baseline_null.py
    "brain": {"mantel": 2026073101, "degree": 2026073102},
    "pbmc": {"mantel": 2026073103, "degree": 2026073104},
}
PROBE_SEED_ROOT = 20260730  # src/pair_probe_stats.py

GROUPS = None
SEEDS = None
BASE = None


# ----------------------------- inputs ----------------------------------------
def load_baseline_inputs():
    """Rebuild the baseline scripts' inputs exactly (they use the PBMC confound cache
    for gene length, detection and GC in both tissues, and their own peak counts)."""
    import brain_coexp_baseline_null as bb
    import pbmc_coexp_baseline_null as pb
    cache = np.load(f"{fpa.OUT}/pbmc_confounds_v2.npz", allow_pickle=False)
    gl, dv, gc = cache["genelen"], cache["detv"], cache["gc"]
    genes, _, _ = fpa.load_manifest()
    out = {}
    for tissue, npz, fm_npz, pc_fn in (
        ("brain", "G_ATAC_v2_GSE174367.npz", "fmgraphs_pooled_v2.npz", bb.brain_peakcount),
        ("pbmc", "G_ATAC_v2_PBMC10k.npz", "pbmc_fmgraphs_pooled.npz", pb.pbmc_peakcount),
    ):
        Z = np.load(f"{fpa.OUT}/{npz}", allow_pickle=False)
        types = [str(t) for t in Z["types"]]
        tf = np.array(Z["tf_rows"])
        G = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
        co = np.load(f"{fpa.OUT}/{fm_npz}")["co"].astype(np.float32)
        pc = pc_fn()
        Ng = G.shape[0]
        od = (G > 0).sum(1).astype(np.float32)
        ind = (G > 0).sum(0).astype(np.float32)
        ii0 = np.repeat(tf, Ng)
        jj0 = np.tile(np.arange(Ng), len(tf))
        m = fpa.edge_mask(tissue, genes, tf, ii0, jj0)
        ii, jj = ii0[m], jj0[m]
        out[tissue] = dict(G=G, tf=tf, ii=ii, jj=jj, co_v=co[ii, jj], at_v=G[ii, jj],
                           pc=pc, gl=gl, dv=dv, gc=gc, od=od, ind=ind)
    return out


def _init(groups, seeds, base, probe):
    global GROUPS, SEEDS, BASE, PROBE
    GROUPS, SEEDS, BASE, PROBE = groups, seeds, base, probe


# ----------------------------- null arrays -----------------------------------
def baseline_mantel_array(b, n_perm, seed):
    """Same loop as fpa.mantel_randomization (use_coexp=False, full), returning the array."""
    Ng = b["G"].shape[0]
    N = len(b["ii"])
    fm_r = rankdata(b["co_v"])
    C = np.empty((N, 7))
    C[:, 0] = 1.0
    C[:, 1] = fpa.zscore(b["pc"][b["jj"]]); C[:, 2] = fpa.zscore(b["gl"][b["jj"]])
    C[:, 3] = fpa.zscore(b["dv"][b["jj"]]); C[:, 4] = fpa.zscore(b["gc"][b["jj"]])
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(Ng)
        ap_r = rankdata(b["G"][perm[b["ii"]], perm[b["jj"]]])
        C[:, 5] = fpa.zscore(b["od"][perm][b["ii"]])
        C[:, 6] = fpa.zscore(b["ind"][perm][b["jj"]])
        null[k] = fpa.pcorr_inplace(fm_r, ap_r, C)
    return null


def baseline_degree_array(b, n_perm, seed):
    """Same loop as fpa.degree_preserving_null (use_coexp=False, full), returning the array."""
    Ng = b["G"].shape[0]
    N = len(b["ii"])
    fm_r = rankdata(b["co_v"])
    jj = b["jj"]
    C = np.column_stack([np.ones(N), fpa.zscore(b["pc"][jj]), fpa.zscore(b["gl"][jj]),
                         fpa.zscore(b["dv"][jj]), fpa.zscore(b["gc"][jj]),
                         fpa.zscore(b["od"][b["ii"]])])
    rng = np.random.default_rng(seed)
    tf_list = np.unique(b["tf"])
    non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, Ng)]) for t in tf_list}
    Gp = b["G"].copy()
    null = np.empty(n_perm)
    for k in range(n_perm):
        for t in tf_list:
            nz = non_self[t]
            vals = Gp[t, nz].copy()
            rng.shuffle(vals)
            Gp[t, nz] = vals
        ap_r = rankdata(Gp[b["ii"], jj])
        indeg_p = (Gp > 0).sum(0).astype(np.float32)
        null[k] = fpa.pcorr_inplace(fm_r, ap_r, np.column_stack([C, fpa.zscore(indeg_p[jj])]))
    return null


def _fwl_corr(x_fixed_resid, y_fixed_resid, X_changing, Q):
    """Partial correlation on [fixed | changing] via Frisch-Waugh-Lovell, as the batched nulls do."""
    if X_changing.shape[1]:
        Xo = X_changing - Q @ (Q.T @ X_changing)
        U = fpa._column_space_basis(Xo)
        x = x_fixed_resid - U @ (U.T @ x_fixed_resid)
        y = y_fixed_resid - U @ (U.T @ y_fixed_resid)
    else:
        x, y = x_fixed_resid, y_fixed_resid
    x = x - x.mean()
    y = y - y.mean()
    return float(x @ y / np.sqrt((x @ x) * (y @ y)))


def baseline_mantel_array_fast(b, n_perm, seed):
    """Same random stream and design as baseline_mantel_array, computed by FWL (much faster)."""
    Ng = b["G"].shape[0]
    N = len(b["ii"])
    jj, ii = b["jj"], b["ii"]
    Q = fpa._column_space_basis(np.column_stack([
        np.ones(N), fpa.zscore(b["pc"][jj]), fpa.zscore(b["gl"][jj]),
        fpa.zscore(b["dv"][jj]), fpa.zscore(b["gc"][jj])]))
    fm_r = rankdata(b["co_v"])
    fm_res = fm_r - Q @ (Q.T @ fm_r)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(Ng)
        ap_r = rankdata(b["G"][perm[ii], perm[jj]])
        ch = np.column_stack([fpa.zscore(b["od"][perm][ii]), fpa.zscore(b["ind"][perm][jj])])
        null[k] = _fwl_corr(fm_res, ap_r - Q @ (Q.T @ ap_r), ch, Q)
    return null


def baseline_degree_array_fast(b, n_perm, seed):
    """Same random stream and design as baseline_degree_array, computed by FWL."""
    Ng = b["G"].shape[0]
    N = len(b["ii"])
    jj, ii = b["jj"], b["ii"]
    Q = fpa._column_space_basis(np.column_stack([
        np.ones(N), fpa.zscore(b["pc"][jj]), fpa.zscore(b["gl"][jj]),
        fpa.zscore(b["dv"][jj]), fpa.zscore(b["gc"][jj]), fpa.zscore(b["od"][ii])]))
    fm_r = rankdata(b["co_v"])
    fm_res = fm_r - Q @ (Q.T @ fm_r)
    rng = np.random.default_rng(seed)
    tf_list = np.unique(b["tf"])
    non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, Ng)]) for t in tf_list}
    Gp = b["G"].copy()
    null = np.empty(n_perm)
    for k in range(n_perm):
        for t in tf_list:
            nz = non_self[t]
            vals = Gp[t, nz].copy()
            rng.shuffle(vals)
            Gp[t, nz] = vals
        ap_r = rankdata(Gp[ii, jj])
        indeg_p = (Gp > 0).sum(0).astype(np.float32)
        null[k] = _fwl_corr(fm_res, ap_r - Q @ (Q.T @ ap_r),
                            fpa.zscore(indeg_p[jj]).reshape(-1, 1), Q)
    return null


def run_job(job):
    kind, group, null_type, spec, n_perm, out_dir = job
    t0 = time.time()
    path = Path(out_dir) / f"{kind}_{group}_{null_type}_{spec}.npz"
    if kind == "fm":
        g = GROUPS[group]
        seed = SEEDS[group][(null_type, spec)]
        kw = dict(fm_vecs=g["fm_vecs"], co_v=g["co_v"], jj=g["jj"], ii=g["ii"],
                  peakcount=g["peakcount"], genelen=g["genelen"], detv=g["detv"], gc=g["gc"],
                  tf_outdeg_full=g["tf_outdeg"], atac_indeg_full=g["atac_indeg"],
                  G_atac_full=g["G"], use_coexp=True, confound_spec=spec,
                  n_perm=n_perm, seed=seed)
        if null_type == "mantel":
            nulls, meta = fpa.batched_mantel_null(**kw)
        else:
            nulls, meta = fpa.batched_degree_preserving_null(tf_rows_unique=g["tf_rows"], **kw)
        arrays = {lab: arr for lab, arr in zip(g["labels"], nulls)}
        arrays["_replicate_seeds"] = np.array(meta["replicate_seeds"], dtype=np.uint64)
        np.savez(path, **arrays)
    elif kind == "baseline":
        b = BASE[group]
        seed = BASELINE_SEEDS[group][null_type]
        fn = baseline_mantel_array_fast if null_type == "mantel" else baseline_degree_array_fast
        np.savez(path, co_expression=fn(b, n_perm, seed))
    elif kind == "probe":
        fam_index = int(group)
        np.savez(path, null=probe_family_null(PROBE, fam_index, n_perm))
    else:
        raise ValueError(kind)
    return job[:4], round(time.time() - t0, 1), str(path)


# ----------------------------- probe -----------------------------------------
PROBE = None


def probe_context():
    """Inputs of pair_probe_stats.py (same arm, features, confound design, ridge path)."""
    sys.path.insert(0, str(common.SRC / "v2"))
    import pair_probe_stats as pps
    from sklearn.linear_model import RidgeCV
    ROOT = Path(fpa.ROOT)
    pair_dir = ROOT / "results/v2/tf_probe_pair"
    ev = json.loads((ROOT / "results/v2/tf_probe_pair_eval_v2.json").read_text())
    arm = ev["primary_arm"]
    families = sorted(ev["results"][arm])
    n_genes, n_tf = ev["n_genes"], ev["n_test_tfs"]
    prov = json.loads((pair_dir / "provenance.json").read_text())
    cols = [prov["feature_names"].index(f) for f in ev["arms"][arm]]
    tg = np.load(pair_dir / "pair_targets.npz", allow_pickle=False)
    y_train = tg["y_train"].astype(np.float64)
    y_test = tg["y_test"].astype(np.float64).reshape(n_tf, n_genes)
    cg = tg["conf_gene"].astype(np.float64)
    cz = (cg - cg.mean(0)) / np.where(cg.std(0) == 0, 1.0, cg.std(0))
    resid = pps.make_residualiser(np.column_stack([np.ones(n_genes), cz]))
    preds, per_tf, obs = {}, {}, {}
    for fam in families:
        d = np.load(pair_dir / f"{fam}_pairs.npz", allow_pickle=False)
        xtr = d["X_train"].astype(np.float64)[:, cols]
        xte = d["X_test"].astype(np.float64)[:, cols]
        mu, sd = xtr.mean(0), xtr.std(0)
        sd = np.where(sd == 0, 1.0, sd)
        model = RidgeCV(alphas=pps.ALPHAS, scoring="neg_mean_squared_error")
        model.fit((xtr - mu) / sd, y_train)
        preds[fam] = model.predict((xte - mu) / sd).reshape(n_tf, n_genes)
        per_tf[fam] = pps.adjusted_per_tf(y_test, preds[fam], resid)
        obs[fam] = pps.mean_over_valid(per_tf[fam], fam)
    return {"pps": pps, "families": families, "y_test": y_test, "resid": resid,
            "preds": preds, "per_tf": per_tf, "obs": obs, "n_tf": n_tf}


def probe_family_null(ctx, fam_index, n_perm):
    pps, fam = ctx["pps"], ctx["families"][fam_index]
    rng = np.random.default_rng(PROBE_SEED_ROOT * 1000 + fam_index)
    pred = ctx["preds"][fam]
    null = np.empty(n_perm)
    for b in range(n_perm):
        shuffled = np.stack([rng.permutation(pred[i]) for i in range(ctx["n_tf"])])
        null[b] = pps.mean_over_valid(pps.adjusted_per_tf(ctx["y_test"], shuffled, ctx["resid"]), fam)
    return null


def probe_signflip(ctx, n_perm):
    pps = ctx["pps"]
    out = {}
    for fam in ctx["families"]:
        if fam == pps.BASELINE:
            continue
        diff = ctx["per_tf"][fam] - ctx["per_tf"][pps.BASELINE]
        diff = diff[np.isfinite(diff)]
        rng = np.random.default_rng(PROBE_SEED_ROOT + 1)
        signs = rng.choice([-1.0, 1.0], size=(n_perm, diff.size))
        out[fam] = (float(diff.mean()), (signs * diff).mean(axis=1))
    return out


# ----------------------------- assemble --------------------------------------
def summarize(null, observed, prefix=999):
    out = {"p": common.plus_one_p(null, observed), "N": int(len(null)),
           "null_sd": float(null.std()), "null_mean": float(null.mean()),
           "exceed": int(np.sum(np.abs(null) >= abs(observed)))}
    if len(null) >= prefix:
        out["p_prefix999"] = common.plus_one_p(null[:prefix], observed)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=9999)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", default="", help="kind:group:null:spec filter, comma-separated")
    ap.add_argument("--skip-probe", action="store_true")
    args = ap.parse_args()

    tag = f"n{args.n_perm}"
    out_dir = common.OUT_DIR / f"resolution_{tag}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"loading inputs; output -> {out_dir}")
    groups = common.load_groups()
    seeds = common.audit_seeds()
    base = load_baseline_inputs()
    frozen = common.frozen_rows()

    # Seed derivation must reproduce the frozen per-test seeds exactly.
    for group, g in groups.items():
        tissue = g["tissue"]
        for lab in g["labels"]:
            for spec in ("full", "non_degree"):
                r = frozen[(tissue, spec, lab)]
                assert r["mantel"]["seed"] == seeds[group][("mantel", spec)], (group, lab, spec)
                assert r["degree_preserving"]["seed"] == seeds[group][("degree", spec)], (group, lab, spec)
    log("seed derivation reproduces all frozen per-test seeds")

    jobs = [("fm", grp, nt, sp, args.n_perm, str(out_dir))
            for grp in groups for nt in ("mantel", "degree") for sp in ("full", "non_degree")]
    jobs += [("baseline", t, nt, "full", args.n_perm, str(out_dir))
             for t in ("brain", "pbmc") for nt in ("mantel", "degree")]
    probe = None
    if not args.skip_probe and not args.only:
        log("probe context (ridge fits)")
        probe = probe_context()
        jobs += [("probe", str(i), "family", "na", args.n_perm, str(out_dir))
                 for i in range(len(probe["families"]))]
    if args.only:
        keep = {tuple(s.split(":")) for s in args.only.split(",")}
        jobs = [j for j in jobs if (j[0], j[1], j[2], j[3]) in keep or (j[1], j[2], j[3]) in keep]
    # Heaviest (brain, 8 rows) first so the pool drains evenly.
    jobs.sort(key=lambda j: (j[1] != "brain", j[0] != "fm"))
    done = [j for j in jobs if (out_dir / f"{j[0]}_{j[1]}_{j[2]}_{j[3]}.npz").exists()]
    jobs = [j for j in jobs if j not in done]
    log(f"{len(done)} jobs already on disk (resumed); {len(jobs)} null jobs at N={args.n_perm} "
        f"on {args.workers} workers")
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers, initializer=_init, initargs=(groups, seeds, base, probe)) as pool:
        for key, secs, _path in pool.imap_unordered(run_job, jobs):
            log(f"  done {key} in {secs}s")

    if probe is not None:
        arrays = {f"family_{fam}": np.load(out_dir / f"probe_{i}_family_na.npz")["null"]
                  for i, fam in enumerate(probe["families"])}
        flips = probe_signflip(probe, args.n_perm)
        for fam, (_, null) in flips.items():
            arrays[f"signflip_{fam}"] = null
        np.savez(out_dir / "probe_nulls.npz", **arrays)
        probe = {"families": probe["families"], "observed": probe["obs"],
                 "signflip_observed": {fam: v[0] for fam, v in flips.items()}}

    # ---- assemble per-row results and verify the frozen prefix ----
    obs_exact = common.observed_unrounded(groups)
    rows, mismatches = [], []
    for group, g in groups.items():
        tissue = g["tissue"]
        for spec in ("full", "non_degree"):
            for null_type, key in (("mantel", "mantel"), ("degree", "degree_preserving")):
                f = out_dir / f"fm_{group}_{null_type}_{spec}.npz"
                if not f.exists():
                    continue
                Z = np.load(f)
                for lab in g["labels"]:
                    fr = frozen[(tissue, spec, lab)]
                    obs = obs_exact[(tissue, spec, lab)]
                    s = summarize(Z[lab], obs)
                    s.update(tissue=tissue, spec=spec, model_label=lab, null=null_type,
                             group=group, observed_partial_rho=obs,
                             frozen_p=fr[key]["p_mc"], seed=seeds[group][(null_type, spec)])
                    if "p_prefix999" in s and abs(s["p_prefix999"] - fr[key]["p_mc"]) > 1e-12:
                        mismatches.append(s)
                    rows.append(s)
    base_rows = []
    frozen_base = {t: json.loads((Path(fpa.OUT) / f"{t}_coexp_baseline_null_v2.json").read_text())
                   for t in ("brain", "pbmc")}
    for t in ("brain", "pbmc"):
        for null_type, key in (("mantel", "pM"), ("degree", "pD")):
            f = out_dir / f"baseline_{t}_{null_type}_full.npz"
            if not f.exists():
                continue
            obs = frozen_base[t]["observed_rho"]
            s = summarize(np.load(f)["co_expression"], obs)
            s.update(tissue=t, null=null_type, frozen_p=frozen_base[t][key],
                     seed=BASELINE_SEEDS[t][null_type])
            if "p_prefix999" in s and abs(s["p_prefix999"] - frozen_base[t][key]) > 1e-12:
                mismatches.append(s)
            base_rows.append(s)

    probe_doc = None
    if probe is not None:
        frozen_probe = json.loads((Path(fpa.OUT) / "tf_probe_pair_stats_v2.json").read_text())
        P = np.load(out_dir / "probe_nulls.npz")
        probe_doc = {"families": {}, "contrasts": {}}
        for fam in probe["families"]:
            s = summarize(P[f"family_{fam}"], probe["observed"][fam])
            s["frozen_p"] = frozen_probe["families"][fam]["mantel_p"]
            if abs(s["p_prefix999"] - s["frozen_p"]) > 1e-12:
                mismatches.append({"probe_family": fam, **s})
            probe_doc["families"][fam] = s
        for fam, d in probe["signflip_observed"].items():
            s = summarize(P[f"signflip_{fam}"], d)
            s["frozen_p"] = frozen_probe["contrasts_vs_baseline"][fam]["signflip_p"]
            if abs(s["p_prefix999"] - s["frozen_p"]) > 1e-12:
                mismatches.append({"probe_contrast": fam, **s})
            probe_doc["contrasts"][fam] = s

    doc = {
        "schema_version": 1,
        "analysis": "resolution_check",
        "purpose": ("R1.1/R2.4: rerun every frozen Monte Carlo test with its frozen seed at a "
                    "deeper N. The frozen N=999 audit remains the registered primary analysis."),
        "n_perm": args.n_perm,
        "prefix_contract": ("replicates are consumed in seed order, so the first 999 replicates "
                            "equal the frozen replicates; p_prefix999 must equal the frozen p_mc"),
        "prefix_mismatches": mismatches,
        "fm_rows": rows,
        "baseline_rows": base_rows,
        "probe": probe_doc,
        "inputs": {"fixed_panel_audit_v2_sha256": common.sha256_file(common.AUDIT_JSON)},
    }
    common.write_json(out_dir / f"resolution_{tag}.json", doc)
    log(f"wrote {out_dir / f'resolution_{tag}.json'}; prefix mismatches: {len(mismatches)}")
    if mismatches:
        sys.exit(2)


if __name__ == "__main__":
    main()
