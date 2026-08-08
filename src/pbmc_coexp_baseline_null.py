#!/usr/bin/env python
"""Co-expression baseline through both nulls, PBMC tissue.

Mirrors brain_coexp_baseline_null.py: the baseline is scored at its own primary
rung (no self-partial, use_coexp=False, full confounds) because partialling
co-expression out of co-expression is degenerate by construction. Same edge set,
confounds, and null machinery as the pooled audit rows. Explicit integer seeds
make both Monte Carlo tests reproducible across Python processes.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, ".")
import fixed_panel_audit as fpa  # noqa: E402

OUT = fpa.OUT
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
ATAC_P = os.path.join(fpa.ROOT, "data", "multiome", "pbmc10k_atac.h5ad")
PROM = fpa.PROM
PBMC_MANTEL_SEED = 2026073103
PBMC_DEGREE_SEED = 2026073104


def pbmc_peakcount():
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
    import anndata as ad
    Av = ad.read_h5ad(ATAC_P, backed="r")
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


def main():
    Z = np.load(f"{OUT}/G_ATAC_v2_PBMC10k.npz", allow_pickle=False)
    types = [str(t) for t in Z["types"]]
    tf = np.array(Z["tf_rows"])
    G = np.mean([Z[f"G_{t}"] for t in types], axis=0).astype(np.float32)
    F = np.load(f"{OUT}/pbmc_fmgraphs_pooled.npz")
    co = F["co"].astype(np.float32)
    cache = np.load(f"{OUT}/pbmc_confounds_v2.npz", allow_pickle=False)
    gl, dv, gc = cache["genelen"], cache["detv"], cache["gc"]
    print("pbmc peakcount...", flush=True)
    pc = pbmc_peakcount()
    genes, _, _ = fpa.load_manifest()
    Ng = G.shape[0]
    od = (G > 0).sum(1).astype(np.float32)
    ind = (G > 0).sum(0).astype(np.float32)
    ii0 = np.repeat(tf, Ng)
    jj0 = np.tile(np.arange(Ng), len(tf))
    m = fpa.edge_mask("pbmc", genes, tf, ii0, jj0)
    ii, jj = ii0[m], jj0[m]
    co_v, at_v = co[ii, jj], G[ii, jj]
    print(f"pbmc edges {len(ii)}", flush=True)
    obs = fpa.partial_rho_obs_sliced(co_v, at_v, co_v, jj, ii, pc, gl, dv, gc,
                                     od, ind, False, "full")
    print(f"pbmc coexp baseline obs (degree_only, no self-partial): {obs:.6f}", flush=True)
    t = time.time()
    man = fpa.mantel_randomization(co_v, at_v, co_v, jj, ii, pc, gl, dv, gc,
                                   od, ind, G, False, "full", obs, 999,
                                   seed=PBMC_MANTEL_SEED)
    deg = fpa.degree_preserving_null(co_v, at_v, co_v, jj, ii, pc, gl, dv, gc,
                                     od, ind, G, np.unique(tf), False, "full", obs, 999,
                                     seed=PBMC_DEGREE_SEED)
    print(f"pbmc baseline: rho={obs:+.6f} pM={man['p_mc']} pD={deg['p_mc']} "
          f"zM={man['z']:.2f} zD={deg['z']:.2f}  ({time.time()-t:.0f}s)", flush=True)
    document = {
        "schema_version": 2,
        "tissue": "pbmc", "model_label": "co_expression_baseline",
        "rung": "degree_only_no_selfpartial", "observed_rho": float(obs),
        "pM": man["p_mc"], "pD": deg["p_mc"], "zM": man["z"], "zD": deg["z"],
        "n_perm": 999, "mantel_seed": PBMC_MANTEL_SEED,
        "degree_seed": PBMC_DEGREE_SEED, "seed_contract": "explicit_integer_v1",
    }
    fpa.write_json_atomic(f"{OUT}/pbmc_coexp_baseline_null_v2.json", document)
    print("wrote pbmc_coexp_baseline_null_v2.json", flush=True)


if __name__ == "__main__":
    main()
