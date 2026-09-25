#!/usr/bin/env python
"""Deterministic rerun of the brain attention-omission guard (Figure 6D).

The submitted guard (src/v2/verify_brain_attention.py) drew its seeds from Python's
salted hash() and took GC from the PBMC covariate cache, so it was neither
reproducible nor on the audit's brain design. This rerun uses the audit's brain
design (brain GC), the batched randomization functions, fixed seeds, and N = 999:
  guard A: the seven brain rows without Geneformer attention, BH over seven rows;
  guard B: all eight rows with the new seed stream, BH over eight rows.
The comparison with the frozen audit is by Support decision, as in the paper.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import fpa, log  # noqa: E402
from multiplicity import bh  # noqa: E402


def main():
    g = common.load_groups()["brain"]
    frozen = common.frozen_rows()
    seeds = [int(s.generate_state(1, dtype=np.uint64)[0])
             for s in np.random.SeedSequence([common.SEED_ROOT_REVISION, 4]).spawn(2)]
    kw = dict(co_v=g["co_v"], jj=g["jj"], ii=g["ii"], peakcount=g["peakcount"],
              genelen=g["genelen"], detv=g["detv"], gc=g["gc"], tf_outdeg_full=g["tf_outdeg"],
              atac_indeg_full=g["atac_indeg"], G_atac_full=g["G"], use_coexp=True,
              confound_spec="full", n_perm=999)
    labels = g["labels"]
    obs_exact = common.observed_unrounded({"brain": g})
    obs = [obs_exact[("brain", "full", lab)] for lab in labels]
    nm, _ = fpa.batched_mantel_null(fm_vecs=g["fm_vecs"], seed=seeds[0], **kw)
    nd, _ = fpa.batched_degree_preserving_null(fm_vecs=g["fm_vecs"], tf_rows_unique=g["tf_rows"],
                                               seed=seeds[1], **kw)
    pM = np.array([common.plus_one_p(nm[r], obs[r]) for r in range(len(labels))])
    pD = np.array([common.plus_one_p(nd[r], obs[r]) for r in range(len(labels))])
    audit = {lab: (frozen[("brain", "full", lab)]["mantel"]["bh_q_family"] < 0.05
                   and frozen[("brain", "full", lab)]["degree_preserving"]["bh_q_family"] < 0.05)
             for lab in labels}
    out = {"schema_version": 1, "analysis": "attention_omission_guard_rerun", "n_perm": 999,
           "seeds": {"mantel": seeds[0], "degree": seeds[1]},
           "seed_contract": "SeedSequence([20260924, 4]).spawn(2)", "guards": {}}
    for name, keep in (("without_attention", [l != "geneformer_attn" for l in labels]),
                       ("all_eight", [True] * len(labels))):
        idx = [i for i, k in enumerate(keep) if k]
        qM, qD = bh(pM[idx]), bh(pD[idx])
        rows = []
        for j, i in enumerate(idx):
            sup = bool(qM[j] < 0.05 and qD[j] < 0.05)
            rows.append({"model_label": labels[i], "pM": float(pM[i]), "pD": float(pD[i]),
                         "qM": float(qM[j]), "qD": float(qD[j]), "support": sup,
                         "audit_support": audit[labels[i]], "agrees": sup == audit[labels[i]]})
        out["guards"][name] = {"rows": rows, "all_agree": all(r["agrees"] for r in rows)}
        log(name, "agrees with audit:", out["guards"][name]["all_agree"],
            [(r["model_label"], r["support"]) for r in rows])
    common.write_json(common.OUT_DIR / "attention_guard.json", out)


if __name__ == "__main__":
    main()
