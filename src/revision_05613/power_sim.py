#!/usr/bin/env python
"""R1.4: simulation check of the analytic MDE for the full-specification tests.

Uses the frozen audit's own planted-signal construction (run_sensitivity in
run_fixed_panel_audit.py): synth = a * z(proxy residual) + (1 - a) * z(noise),
where the proxy residual is the proxy's rank residual on ranked co-expression
and the six structural covariates. Each synthetic graph is scored with the same
statistic and the same two randomization tests (N = 999, full specification) as
the audit rows. Power at each injected fraction is the share of replicates that
the test declares unusual; the empirical MDE is the observed rho at which power
reaches 0.80. This checks the test, not regulatory biology.
"""
import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path
from statistics import NormalDist

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import fpa, log  # noqa: E402

A_GRID = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.010, 0.015]
FAMILY_M = {"brain": 8, "pbmc": 5}
GROUP = {"brain": "brain", "pbmc": "pbmc_shared"}
DATA = None


def synth_rows(g, n_rep, seed):
    jj, ii = g["jj"], g["ii"]
    C = np.column_stack([rankdata(g["co_v"]), fpa.zscore(g["peakcount"][jj]), fpa.zscore(g["genelen"][jj]),
                         fpa.zscore(g["detv"][jj]), fpa.zscore(g["gc"][jj]),
                         fpa.zscore(g["tf_outdeg"][ii]), fpa.zscore(g["atac_indeg"][jj])])
    za = fpa.zscore(fpa.resid(rankdata(g["atac_v"]), C))
    seeds = [int(s.generate_state(1, dtype=np.uint64)[0])
             for s in np.random.SeedSequence(seed).spawn(len(A_GRID) * n_rep)]
    rows, meta, it = [], [], iter(seeds)
    for a in A_GRID:
        for rep in range(n_rep):
            s = next(it)
            noise = np.random.default_rng(s).standard_normal(len(za))
            rows.append(a * za + (1 - a) * fpa.zscore(noise))
            meta.append((a, rep, s))
    return rows, meta


def run_job(job):
    tissue, null_type, n_perm, n_rep, seed_rows, seed_null, out_dir = job
    t0 = time.time()
    g = DATA[GROUP[tissue]]
    rows, meta = synth_rows(g, n_rep, seed_rows)
    obs = np.array([fpa.partial_rho_obs_sliced(
        fm_v=r, atac_v=g["atac_v"], co_v=g["co_v"], jj=g["jj"], ii=g["ii"],
        peakcount=g["peakcount"], genelen=g["genelen"], detv=g["detv"], gc=g["gc"],
        tf_outdeg=g["tf_outdeg"], atac_indeg=g["atac_indeg"], use_coexp=True,
        confound_spec="full") for r in rows])
    kw = dict(fm_vecs=rows, co_v=g["co_v"], jj=g["jj"], ii=g["ii"], peakcount=g["peakcount"],
              genelen=g["genelen"], detv=g["detv"], gc=g["gc"], tf_outdeg_full=g["tf_outdeg"],
              atac_indeg_full=g["atac_indeg"], G_atac_full=g["G"], use_coexp=True,
              confound_spec="full", n_perm=n_perm, seed=seed_null)
    if null_type == "mantel":
        nulls, _ = fpa.batched_mantel_null(**kw)
    else:
        nulls, _ = fpa.batched_degree_preserving_null(tf_rows_unique=g["tf_rows"], **kw)
    nulls = np.column_stack(nulls)
    np.savez(Path(out_dir) / f"power_{tissue}_{null_type}.npz", observed=obs, nulls=nulls,
             alpha=np.array([m[0] for m in meta]), rep=np.array([m[1] for m in meta]),
             noise_seed=np.array([m[2] for m in meta], dtype=np.uint64))
    return tissue, null_type, round(time.time() - t0, 1)


def _init(data):
    global DATA
    DATA = data


def assemble(out_dir):
    z80 = NormalDist().inv_cdf(0.80)
    out = {}
    for tissue in ("brain", "pbmc"):
        M = np.load(Path(out_dir) / f"power_{tissue}_mantel.npz")
        D = np.load(Path(out_dir) / f"power_{tissue}_degree.npz")
        assert np.allclose(M["observed"], D["observed"])
        obs, a = M["observed"], M["alpha"]
        pM = np.array([common.plus_one_p(M["nulls"][:, r], obs[r]) for r in range(len(obs))])
        pD = np.array([common.plus_one_p(D["nulls"][:, r], obs[r]) for r in range(len(obs))])
        sdM, sdD = M["nulls"].std(axis=0), D["nulls"].std(axis=0)
        m = FAMILY_M[tissue]
        levels = {"a05": 0.05, "a05m": 0.05 / m}
        curve = []
        for av in A_GRID:
            sel = a == av
            rec = {"alpha_injected": av, "n": int(sel.sum()),
                   "mean_rho": float(obs[sel].mean()), "sd_rho": float(obs[sel].std())}
            for tag, lev in levels.items():
                rec[f"power_M_{tag}"] = float(np.mean(pM[sel] <= lev))
                rec[f"power_D_{tag}"] = float(np.mean(pD[sel] <= lev))
                rec[f"power_dual_{tag}"] = float(np.mean((pM[sel] <= lev) & (pD[sel] <= lev)))
            curve.append(rec)
        emp = {}
        for tag, lev in levels.items():
            pw = np.array([c[f"power_dual_{tag}"] for c in curve])
            rho = np.array([c["mean_rho"] for c in curve])
            k = int(np.argmax(pw >= 0.8)) if np.any(pw >= 0.8) else None
            if k is None or k == 0:
                emp[tag] = None
            else:
                t = (0.8 - pw[k - 1]) / (pw[k] - pw[k - 1])
                emp[tag] = float(rho[k - 1] + t * (rho[k] - rho[k - 1]))
        sd_dual = float(np.median(np.maximum(sdM, sdD)))
        analytic = {tag: (NormalDist().inv_cdf(1 - lev / 2) + z80) * sd_dual for tag, lev in levels.items()}
        out[tissue] = {"curve": curve, "empirical_mde_dual": emp, "analytic_mde_dual": analytic,
                       "median_null_sd": {"mantel": float(np.median(sdM)), "degree": float(np.median(sdD))},
                       "false_positive_rate_at_alpha0": {
                           tag: curve[0][f"power_dual_{tag}"] for tag in levels}}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=999)
    ap.add_argument("--n-rep", type=int, default=20)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--assemble-only", action="store_true")
    args = ap.parse_args()
    out_dir = common.OUT_DIR / f"power_sim_n{args.n_perm}_r{args.n_rep}"
    out_dir.mkdir(parents=True, exist_ok=True)
    root = np.random.SeedSequence([common.SEED_ROOT_REVISION, 2])  # R1.4 stream
    s = [int(x.generate_state(1, dtype=np.uint64)[0]) for x in root.spawn(6)]
    seeds = {"brain": (s[0], s[1], s[2]), "pbmc": (s[3], s[4], s[5])}  # rows, mantel, degree
    if not args.assemble_only:
        groups = common.load_groups()
        data = {GROUP[t]: groups[GROUP[t]] for t in ("brain", "pbmc")}
        jobs = []
        for t in ("brain", "pbmc"):
            jobs.append((t, "mantel", args.n_perm, args.n_rep, seeds[t][0], seeds[t][1], str(out_dir)))
            jobs.append((t, "degree", args.n_perm, args.n_rep, seeds[t][0], seeds[t][2], str(out_dir)))
        log(f"{len(jobs)} jobs; {len(A_GRID)} alphas x {args.n_rep} replicates; N={args.n_perm}")
        with mp.get_context("fork").Pool(args.workers, initializer=_init, initargs=(data,)) as pool:
            for t, nt, secs in pool.imap_unordered(run_job, jobs):
                log(f"  done {t} {nt} in {secs}s")
    res = {"schema_version": 1, "analysis": "power_simulation_full_spec",
           "design": ("synth = a*z(proxy residual) + (1-a)*z(noise); same statistic and nulls as "
                      "the audit rows; full specification; checks test sensitivity only"),
           "alpha_grid": A_GRID, "n_rep": args.n_rep, "n_perm": args.n_perm,
           "seed_contract": "SeedSequence([20260924, 2]).spawn(6): brain rows/mantel/degree, pbmc rows/mantel/degree",
           "seeds": seeds, "tissues": assemble(out_dir)}
    common.write_json(out_dir / f"power_sim_n{args.n_perm}_r{args.n_rep}.json", res)
    for t, d in res["tissues"].items():
        log(t, "empirical MDE", d["empirical_mde_dual"], "analytic", d["analytic_mde_dual"],
            "FPR", d["false_positive_rate_at_alpha0"])


if __name__ == "__main__":
    main()
