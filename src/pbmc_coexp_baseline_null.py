#!/usr/bin/env python
"""Co-expression baseline through both nulls, PBMC tissue.

Mirrors brain_coexp_baseline_null.py; both delegate to coexp_baseline_null.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coexp_baseline_null as baseline  # noqa: E402
import fixed_panel_audit as fpa  # noqa: E402

ATAC_P = os.path.join(fpa.ROOT, "data", "multiome", "pbmc10k_atac.h5ad")
PBMC_MANTEL_SEED = 2026073103
PBMC_DEGREE_SEED = 2026073104


def main():
    baseline.run_baseline(
        tissue="pbmc",
        atac_file=ATAC_P,
        proxy_npz="G_ATAC_v2_PBMC10k.npz",
        fm_npz="pbmc_fmgraphs_pooled.npz",
        mantel_seed=PBMC_MANTEL_SEED,
        degree_seed=PBMC_DEGREE_SEED,
        out_name="pbmc_coexp_baseline_null_v2.json",
    )


if __name__ == "__main__":
    main()
