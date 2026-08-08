"""
Benchmark pooled-only at N=99 (descriptive per-type, no randomization; cross-tissue
observed-only). Unbuffered progress. Not authoritative; bench-only.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit_utils as au
import fixed_panel_audit as fpa
import run_fixed_panel_audit as drv
from run_fixed_panel_audit import ATAC_BRAIN, ATAC_PBMC

log = au.log

N_PERM = int(os.environ.get("BENCH_N_PERM", "99"))
SEED_ROOT = int(os.environ.get("SEED_ROOT", "20260724"))


def bench_pooled():
    log(f"=== BENCH POOLED N_PERM_MANTEL=N_PERM_DEG={N_PERM} ===")
    ss = np.random.SeedSequence(SEED_ROOT)
    G_atac_b, co_b, models_b, tf_b, types_b = drv.load_pooled_brain()
    ko_models, _ = drv.load_geneformer_ko()
    t0 = time.time()
    pooled_brain = drv.run_pooled_family(ss.spawn(1)[0], ATAC_BRAIN, "brain", G_atac_b, co_b,
                                         models_b, tf_b, ko_models, N_PERM, N_PERM)
    log(f"brain pooled (N={N_PERM}): {time.time() - t0:.2f}s  "
        f"({len(pooled_brain['primary_family']['rows'])} primary, "
        f"{len(pooled_brain['sensitivity_family']['rows'])} sensitivity)")

    G_atac_p, co_p, models_p, tf_p, types_p = drv.load_pooled_pbmc()
    t0 = time.time()
    pooled_pbmc = drv.run_pooled_family(ss.spawn(1)[0], ATAC_PBMC, "pbmc", G_atac_p, co_p,
                                        models_p, tf_p, None, N_PERM, N_PERM)
    log(f"pbmc pooled (N={N_PERM}): {time.time() - t0:.2f}s  "
        f"({len(pooled_pbmc['primary_family']['rows'])} primary, "
        f"{len(pooled_pbmc['sensitivity_family']['rows'])} sensitivity)")
    return pooled_brain, pooled_pbmc


def bench_pertype_descriptive():
    log(f"=== BENCH PER-TYPE DESCRIPTIVE (no randomization, no p_mc) ===")
    ss = np.random.SeedSequence(SEED_ROOT)
    G_atac_b, co_b, _, tf_b, _ = drv.load_pooled_brain()
    type_models_b, _, _ = drv.load_brain_pertype_models()
    ncells_b = json.load(open(f"{fpa.OUT}/pertype_fm_v2.json"))["per_type"]
    ncells_b = {r["cell_type"]: r["n"] for r in ncells_b}
    for cs in ("full", "non_degree"):
        t0 = time.time()
        pt = drv.run_pertype_family(
            ss.spawn(1)[0], ATAC_BRAIN, "brain", G_atac_b, co_b,
            type_models_b, tf_b, cs, ncells_b,
        )
        log(f"brain per-type {cs}: {time.time() - t0:.2f}s  "
            f"({pt['n_rows_exploratory']} FM rows)")

    G_atac_p, co_p, _, tf_p, _ = drv.load_pooled_pbmc()
    type_models_p, _, _ = drv.load_pbmc_pertype_models()
    ncells_p = json.load(open(f"{fpa.OUT}/pbmc_eval_v2.json"))["per_type_coexp"]
    ncells_p = {r["cell_type"]: r["n"] for r in ncells_p}
    for cs in ("full", "non_degree"):
        t0 = time.time()
        pt = drv.run_pertype_family(
            ss.spawn(1)[0], ATAC_PBMC, "pbmc", G_atac_p, co_p,
            type_models_p, tf_p, cs, ncells_p,
        )
        log(f"pbmc per-type {cs}: {time.time() - t0:.2f}s  "
            f"({pt['n_rows_exploratory']} FM rows)")


if __name__ == "__main__":
    t_total = time.time()
    bench_pooled()
    bench_pertype_descriptive()
    log(f"=== TOTAL BENCHMARK: {time.time() - t_total:.2f}s ===")
    log("(production scaling: multiply by N_PERM_production / N_PERM_benchmark; per-type is descriptive so cost is constant)")