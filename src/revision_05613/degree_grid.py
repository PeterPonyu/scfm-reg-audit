#!/usr/bin/env python
"""R1.3: sensitivity of the 13 pooled rows to how graph degree is controlled.

Every specification keeps the frozen base design (intercept, ranked co-expression,
peak count, gene length, RNA detection, GC) and changes only the degree block:

  non_degree      none (frozen prespecified sensitivity)
  full            z(TF out-degree), z(target in-degree)   (frozen primary)
  outdeg_only     z(TF out-degree)
  indeg_only      z(target in-degree)
  log_degree      z(log1p out-degree), z(log1p in-degree)
  rank_degree     ranked out-degree, ranked in-degree
  spline_degree   natural cubic splines of log1p out- and in-degree (5 fixed knots each)
  strength        weighted degree: proxy row sums and column sums
  external_fibro  out/in-degree of the fibroblast-mix proxy (GSE206767; not the tested proxy)
  external_other  out/in-degree of the other FM tissue's proxy (brain <-> PBMC)

Degree terms computed from the tested proxy are recomputed inside each
randomization exactly as in the frozen audit (gene-label null: permuted with the
labels; row-shuffle null: out-degree and row sums invariant, in-degree and column
sums recomputed). External degree terms are gene properties and stay fixed, like
peak count or GC. All specifications share one randomization per replicate, so
their p-values are comparable row by row. The two frozen specifications are
recomputed here with new seeds; their observed rho must match the frozen audit.
"""
import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import fpa, log  # noqa: E402

SPECS = ["non_degree", "full", "outdeg_only", "indeg_only", "log_degree", "rank_degree",
         "spline_degree", "strength", "external_fibro", "external_other"]
GROUP_ORDER = ["brain", "pbmc_shared", "pbmc_scgpt", "pbmc_uce"]
KNOT_Q = (0.05, 0.275, 0.5, 0.725, 0.95)
ALPHA = 0.05
GROUPS = None
EXT = None


def ncs_basis(x, knots):
    """Natural cubic spline basis without intercept (Hastie et al. 2009, eq. 5.4-5.5)."""
    K = len(knots)
    xi = np.asarray(knots, float)

    def d(k):
        return (np.clip(x - xi[k], 0, None) ** 3 - np.clip(x - xi[K - 1], 0, None) ** 3) / (xi[K - 1] - xi[k])

    cols = [x] + [d(k) - d(K - 2) for k in range(K - 2)]
    return np.column_stack(cols)


def zcols(M):
    M = np.atleast_2d(np.asarray(M, float))
    if M.shape[0] == 1 and M.shape[1] > 1:
        M = M.T
    sd = M.std(axis=0)
    return (M - M.mean(axis=0)) / np.where(sd > 0, sd, 1.0)


def degree_block(spec, g, od, ind, so, si, knots):
    ii, jj = g["ii"], g["jj"]
    if spec == "non_degree":
        return np.empty((len(ii), 0))
    if spec == "full":
        return zcols(np.column_stack([od[ii], ind[jj]]))
    if spec == "outdeg_only":
        return zcols(od[ii][:, None])
    if spec == "indeg_only":
        return zcols(ind[jj][:, None])
    if spec == "log_degree":
        return zcols(np.column_stack([np.log1p(od[ii]), np.log1p(ind[jj])]))
    if spec == "rank_degree":
        return zcols(np.column_stack([rankdata(od[ii]), rankdata(ind[jj])]))
    if spec == "spline_degree":
        return zcols(np.column_stack([ncs_basis(np.log1p(od[ii]), knots["od"]),
                                      ncs_basis(np.log1p(ind[jj]), knots["ind"])]))
    if spec == "strength":
        return zcols(np.column_stack([so[ii], si[jj]]))
    if spec in ("external_fibro", "external_other"):
        e = EXT[g["tissue"]][spec]
        return zcols(np.column_stack([e["od"][ii], e["ind"][jj]]))
    raise ValueError(spec)


def resid_block(vecs, X):
    """Residualize each vector (already orthogonal to the base design) on block X."""
    if X.shape[1] == 0:
        return vecs
    U = fpa._column_space_basis(X)
    return [v - U @ (U.T @ v) for v in vecs]


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt((a @ a) * (b @ b))
    return float(a @ b / den) if den > 0 else 0.0


def prepare(g):
    N = len(g["ii"])
    jj = g["jj"]
    base = np.column_stack([np.ones(N), rankdata(g["co_v"]), fpa.zscore(g["peakcount"][jj]),
                            fpa.zscore(g["genelen"][jj]), fpa.zscore(g["detv"][jj]),
                            fpa.zscore(g["gc"][jj])])
    Q = fpa._column_space_basis(base)
    fm_res = [r - Q @ (Q.T @ r) for r in (rankdata(v) for v in g["fm_vecs"])]
    G = g["G"]
    od, ind = (G > 0).sum(1).astype(float), (G > 0).sum(0).astype(float)
    knots = {"od": np.quantile(np.log1p(od[g["ii"]]), KNOT_Q),
             "ind": np.quantile(np.log1p(ind[g["jj"]]), KNOT_Q)}
    return Q, fm_res, knots


def stats_for(g, Q, fm_res, knots, atac_ranks, od, ind, so, si):
    a0 = atac_ranks - Q @ (Q.T @ atac_ranks)
    out = np.empty((len(SPECS), len(fm_res)))
    for s, spec in enumerate(SPECS):
        X = degree_block(spec, g, od, ind, so, si, knots)
        if X.shape[1]:
            X = X - Q @ (Q.T @ X)
        res = resid_block([a0] + fm_res, X)
        a = res[0]
        for r, f in enumerate(res[1:]):
            out[s, r] = corr(f, a)
    return out


def run_job(job):
    group, null_type, n_perm, seed, out_dir = job
    t0 = time.time()
    g = GROUPS[group]
    Q, fm_res, knots = prepare(g)
    G = g["G"].astype(np.float64)
    ii, jj = g["ii"], g["jj"]
    Ng = G.shape[0]
    od, ind = (G > 0).sum(1).astype(float), (G > 0).sum(0).astype(float)
    so, si = G.sum(1), G.sum(0)
    obs = stats_for(g, Q, fm_res, knots, rankdata(G[ii, jj]), od, ind, so, si)
    nulls = np.empty((n_perm, len(SPECS), len(fm_res)))
    ss = np.random.SeedSequence(seed)
    tf_list = np.unique(g["tf_rows"])
    non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, Ng)]) for t in tf_list}
    Gp = G.copy()
    for k in range(n_perm):
        rng = np.random.default_rng(int(ss.spawn(1)[0].generate_state(1, dtype=np.uint64)[0]))
        if null_type == "mantel":
            perm = rng.permutation(Ng)
            ranks = rankdata(G[perm[ii], perm[jj]])
            nulls[k] = stats_for(g, Q, fm_res, knots, ranks, od[perm], ind[perm], so[perm], si[perm])
        else:
            for t in tf_list:
                nz = non_self[t]
                vals = Gp[t, nz].copy()
                rng.shuffle(vals)
                Gp[t, nz] = vals
            ranks = rankdata(Gp[ii, jj])
            nulls[k] = stats_for(g, Q, fm_res, knots, ranks, od,
                                 (Gp > 0).sum(0).astype(float), so, Gp.sum(0))
    path = Path(out_dir) / f"grid_{group}_{null_type}.npz"
    np.savez(path, observed=obs, nulls=nulls, specs=np.array(SPECS),
             labels=np.array(g["labels"]), seed=np.uint64(seed))
    return group, null_type, round(time.time() - t0, 1)


def _init(groups, ext):
    global GROUPS, EXT
    GROUPS, EXT = groups, ext


def external_degrees(groups):
    """Degree vectors from proxies that are not the tested one (same frozen panel)."""
    def load(tag):
        Z = np.load(f"{fpa.OUT}/G_ATAC_v2_{tag}.npz", allow_pickle=False)
        G = np.mean([Z[f"G_{t}"] for t in [str(x) for x in Z["types"]]], axis=0)
        return {"od": (G > 0).sum(1).astype(float), "ind": (G > 0).sum(0).astype(float)}
    fib, brain, pbmc = load("GSE206767"), load("GSE174367"), load("PBMC10k")
    return {"brain": {"external_fibro": fib, "external_other": pbmc},
            "pbmc": {"external_fibro": fib, "external_other": brain}}


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = p[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(q, 1)
    return out


def assemble(out_dir, n_perm):
    frozen = common.frozen_rows()
    by_tissue = {"brain": ["brain"], "pbmc": ["pbmc_shared", "pbmc_scgpt", "pbmc_uce"]}
    result = {"specs": SPECS, "families": {}, "checks": []}
    for tissue, grps in by_tissue.items():
        labels, obs, nulls = [], {}, {}
        for grp in grps:
            M = np.load(Path(out_dir) / f"grid_{grp}_mantel.npz")
            D = np.load(Path(out_dir) / f"grid_{grp}_degree.npz")
            labs = [str(x) for x in M["labels"]]
            for r, lab in enumerate(labs):
                labels.append(lab)
                obs[lab] = M["observed"][:, r]
                assert np.allclose(M["observed"][:, r], D["observed"][:, r])
                nulls[lab] = {"mantel": M["nulls"][:, :, r], "degree": D["nulls"][:, :, r]}
        for s, spec in enumerate(SPECS):
            if spec in ("full", "non_degree"):
                for lab in labels:
                    fr = frozen[(tissue, spec, lab)]["observed_partial_rho"]
                    ok = abs(round(float(obs[lab][s]), 6) - fr) <= 1.5e-6
                    result["checks"].append({"tissue": tissue, "spec": spec, "model_label": lab,
                                             "observed": float(obs[lab][s]), "frozen": fr, "match": ok})
            pm = [common.plus_one_p(nulls[lab]["mantel"][:, s], obs[lab][s]) for lab in labels]
            pd = [common.plus_one_p(nulls[lab]["degree"][:, s], obs[lab][s]) for lab in labels]
            qm, qd = bh(pm), bh(pd)
            rows = []
            for r, lab in enumerate(labels):
                rho = float(obs[lab][s])
                full_rho = float(obs[lab][SPECS.index("full")])
                rows.append({"model_label": lab, "rho": rho,
                             "pM": pm[r], "pD": pd[r], "qM": float(qm[r]), "qD": float(qd[r]),
                             "support": bool(qm[r] < ALPHA and qd[r] < ALPHA),
                             "same_sign_as_full": bool(np.sign(rho) == np.sign(full_rho)),
                             "sdM": float(nulls[lab]["mantel"][:, s].std()),
                             "sdD": float(nulls[lab]["degree"][:, s].std())})
            result["families"][f"{tissue}_{spec}"] = rows
    summary = {}
    for spec in SPECS:
        rows = result["families"][f"brain_{spec}"] + result["families"][f"pbmc_{spec}"]
        summary[spec] = {
            "support": sum(r["support"] for r in rows),
            "support_positive": sum(r["support"] and r["rho"] > 0 for r in rows),
            "support_negative": sum(r["support"] and r["rho"] < 0 for r in rows),
            "sign_agrees_with_full": sum(r["same_sign_as_full"] for r in rows),
            "d_only": sum(r["qD"] < ALPHA <= r["qM"] for r in rows),
            "m_only": sum(r["qM"] < ALPHA <= r["qD"] for r in rows),
        }
    result["summary"] = summary
    result["all_frozen_observed_match"] = all(c["match"] for c in result["checks"])
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=999)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--assemble-only", action="store_true")
    args = ap.parse_args()
    out_dir = common.OUT_DIR / f"degree_grid_n{args.n_perm}"
    out_dir.mkdir(parents=True, exist_ok=True)
    if not args.assemble_only:
        groups = common.load_groups()
        ext = external_degrees(groups)
        root = np.random.SeedSequence([common.SEED_ROOT_REVISION, 3])  # R1.3 stream
        seeds = [int(s.generate_state(1, dtype=np.uint64)[0]) for s in root.spawn(8)]
        jobs, i = [], 0
        for grp in GROUP_ORDER:
            for nt in ("mantel", "degree"):
                jobs.append((grp, nt, args.n_perm, seeds[i], str(out_dir)))
                i += 1
        log(f"{len(jobs)} jobs, N={args.n_perm}, specs={len(SPECS)}; seeds {seeds}")
        with mp.get_context("fork").Pool(args.workers, initializer=_init, initargs=(groups, ext)) as pool:
            for grp, nt, secs in pool.imap_unordered(run_job, jobs):
                log(f"  done {grp} {nt} in {secs}s")
    res = assemble(out_dir, args.n_perm)
    res.update({"schema_version": 1, "analysis": "degree_specification_grid",
                "n_perm": args.n_perm, "seed_root": common.SEED_ROOT_REVISION,
                "seed_contract": "SeedSequence([20260924, 3]).spawn(8) -> one seed per (group, null) "
                                 "in order brain, pbmc_shared, pbmc_scgpt, pbmc_uce x (mantel, degree); "
                                 "replicates spawned from each seed in order",
                "inputs": {"fixed_panel_audit_v2_sha256": common.sha256_file(common.AUDIT_JSON)}})
    common.write_json(out_dir / f"degree_grid_n{args.n_perm}.json", res)
    log("frozen observed rho reproduced:", res["all_frozen_observed_match"])
    for spec, s in res["summary"].items():
        log(f"  {spec:15s} {s}")


if __name__ == "__main__":
    main()
