#!/usr/bin/env python
import json
import os
import tempfile
import time

import numpy as np

import fixed_panel_audit as fpa
import run_fixed_panel_audit as drv


DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))
N_PERM_MANTEL = int(os.environ.get("N_PERM_POOLED_MANTEL", "999"))
N_PERM_DEG = int(os.environ.get("N_PERM_POOLED_DEG", "999"))
OUT_PATH = os.environ.get("PBMC_SCGPT_STATS_OUT", f"{fpa.OUT}/pbmc_scgpt_stats_v2.json")


def log(*args):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, flush=True)


def write_json_atomic(path, document):
    parent = os.path.dirname(path) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=os.path.basename(path), suffix=".json", dir=parent)
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(document, fh, indent=2)
            fh.write("\n")
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def main():
    G_atac, co, models, tf_rows, types = drv.load_pooled_pbmc()
    scgpt = drv.load_optional_pbmc_scgpt()
    if scgpt is None:
        raise FileNotFoundError(f"required PBMC scGPT graph missing: {fpa.OUT}/pbmc_scgpt_pooled_v2.npz")
    scgpt_path, scgpt_co, scgpt_graph = scgpt
    atac_path = f"{fpa.ROOT}/data/multiome/pbmc10k_atac.h5ad"
    result = drv.run_pooled_family(
        np.random.SeedSequence(drv.SEED_ROOT).spawn(2)[1],
        atac_path, "pbmc", G_atac, co, models, tf_rows, None,
        N_PERM_MANTEL, N_PERM_DEG,
    )
    drv.append_independent_control_model(
        result, atac_path, "pbmc", G_atac, scgpt_co, scgpt_graph, tf_rows,
        "scGPT_encoder", N_PERM_MANTEL, N_PERM_DEG,
        np.random.SeedSequence([drv.SEED_ROOT, 20260726, 1]), scgpt_path,
    )
    document = {
        "schema_version": 1,
        "design": "PBMC pooled fixed-panel update with matched-control scGPT",
        "seed_root": drv.SEED_ROOT,
        "n_perm_mantel": N_PERM_MANTEL,
        "n_perm_degree": N_PERM_DEG,
        "types": types,
        "pooled_pbmc": result,
    }
    write_json_atomic(OUT_PATH, document)
    log(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
