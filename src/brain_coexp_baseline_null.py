#!/usr/bin/env python
"""Co-expression baseline through both nulls, brain tissue."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coexp_baseline_null as baseline  # noqa: E402

DATA_ROOT = baseline.data_root()
ATAC_B = os.environ.get(
    "SCFM_BRAIN_ATAC",
    f"{DATA_ROOT}/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")
BRAIN_MANTEL_SEED = 2026073101
BRAIN_DEGREE_SEED = 2026073102


def main():
    baseline.run_baseline(
        tissue="brain",
        atac_file=ATAC_B,
        proxy_npz="G_ATAC_v2_GSE174367.npz",
        fm_npz="fmgraphs_pooled_v2.npz",
        mantel_seed=BRAIN_MANTEL_SEED,
        degree_seed=BRAIN_DEGREE_SEED,
        out_name="brain_coexp_baseline_null_v2.json",
    )


if __name__ == "__main__":
    main()
