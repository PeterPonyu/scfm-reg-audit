#!/usr/bin/env python
"""Brain co-expression baseline and FM-versus-baseline contrast with brain covariates.

Found while preparing the ARRAY-D-26-05613 revision: brain_coexp_baseline_null.py and
fm_vs_baseline_shared_null.py take gene length, detection and GC from the PBMC
confound cache for both tissues ("tissue-independent"). Gene length and detection
are identical across tissues, but GC is computed from each tissue's own linked
peaks (brain vs PBMC GC correlate at r = 0.50). The brain FM rows of the primary
audit use brain GC, so the brain baseline and the brain rows of the shared-null
contrast used a slightly different design from the rows they are compared with.
PBMC is unaffected (its cache equals its own covariates).

This script recomputes, with the frozen seeds and N:
  1. the brain baseline rho and its two randomization p-values (N = 999 and 9,999);
  2. the brain shared-null contrast (Table 5), after first checking that the old
     inputs reproduce the frozen values;
  3. the "rho > baseline" protocol-pass gate for the brain rows.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from common import fpa, log  # noqa: E402
import resolution_check as rc  # noqa: E402

sys.path.insert(0, str(common.SRC / "v2"))
import fm_vs_baseline_shared_null as svb  # noqa: E402


def main():
    frozen_base = json.loads((Path(fpa.OUT) / "brain_coexp_baseline_null_v2.json").read_text())
    frozen_delta = json.loads((Path(fpa.OUT) / "fm_vs_baseline_shared_null_v2.json").read_text())
    base_inputs = rc.load_baseline_inputs()["brain"]           # PBMC-cache gl/dv/gc (frozen)
    pc, gl, gc, dv = fpa.build_confounds(common.ATAC_FILES["brain"])  # brain covariates
    assert np.allclose(pc, base_inputs["pc"]) and np.allclose(gl, base_inputs["gl"]) \
        and np.allclose(dv, base_inputs["dv"]), "only GC is expected to differ"
    b_old = base_inputs
    b_new = dict(base_inputs, gc=gc)

    out = {"schema_version": 1, "analysis": "brain_baseline_covariate_fix",
           "gc_correlation_brain_vs_pbmc_cache": float(np.corrcoef(gc, b_old["gc"])[0, 1])}
    # ---- 1. baseline rho and randomization p ----
    def obs(b):
        return fpa.partial_rho_obs_sliced(b["co_v"], b["at_v"], b["co_v"], b["jj"], b["ii"],
                                          b["pc"], b["gl"], b["dv"], b["gc"], b["od"], b["ind"],
                                          False, "full")
    rho_old, rho_new = obs(b_old), obs(b_new)
    assert abs(rho_old - frozen_base["observed_rho"]) < 1e-9, (rho_old, frozen_base["observed_rho"])
    base = {"rho_frozen": rho_old, "rho_brain_covariates": rho_new}
    chk_m = common.plus_one_p(rc.baseline_mantel_array_fast(b_old, 999, rc.BASELINE_SEEDS["brain"]["mantel"]), rho_old)
    chk_d = common.plus_one_p(rc.baseline_degree_array_fast(b_old, 999, rc.BASELINE_SEEDS["brain"]["degree"]), rho_old)
    assert abs(chk_m - frozen_base["pM"]) < 1e-12 and abs(chk_d - frozen_base["pD"]) < 1e-12, (chk_m, chk_d)
    log("fast baseline path reproduces the frozen brain pM/pD with the old inputs")
    for n in (999, 9999):
        log(f"baseline nulls N={n}")
        nm = rc.baseline_mantel_array_fast(b_new, n, rc.BASELINE_SEEDS["brain"]["mantel"])
        nd = rc.baseline_degree_array_fast(b_new, n, rc.BASELINE_SEEDS["brain"]["degree"])
        base[f"n{n}"] = {"pM": common.plus_one_p(nm, rho_new), "pD": common.plus_one_p(nd, rho_new),
                         "sdM": float(nm.std()), "sdD": float(nd.std())}
    base["frozen"] = {"pM": frozen_base["pM"], "pD": frozen_base["pD"]}
    out["baseline"] = base
    log("baseline", base)

    # ---- 2. shared-null contrast, old inputs (reproduction) then brain covariates ----
    bundle = svb.load_brain_bundle()
    frozen_rows = {r["model_label"]: r for r in frozen_delta["rows"] if r["tissue"] == "brain"}
    res_old = svb.run_tissue(bundle, b_old["pc"], b_old["gl"], b_old["dv"], b_old["gc"], 999)
    for r in res_old["rows_public"]:
        f = frozen_rows[r["model_label"]]
        assert abs(r["p_mc"] - f["p_mc"]) < 1e-12 and abs(r["rho_fm_full"] - f["rho_fm_full"]) < 1e-6, \
            (r["model_label"], r["p_mc"], f["p_mc"])
    log("shared-null contrast: frozen brain values reproduced with the old inputs")
    res_new = svb.run_tissue(bundle, pc, gl, dv, gc, 999)
    primary = common.frozen_rows()
    rows = []
    for r in res_new["rows_public"]:
        lab = r["model_label"]
        audit_rho = primary[("brain", "full", lab)]["observed_partial_rho"]
        rows.append({"model_label": lab, "rho_fm": r["rho_fm_full"], "audit_rho": audit_rho,
                     "rho_fm_matches_audit": abs(r["rho_fm_full"] - audit_rho) <= 1.5e-6,
                     "delta": r["delta_rho"], "p": r["p_mc"], "q": r["bh_q"],
                     "shared_sig": r["shared_null_significant"],
                     "rho_gt_base": bool(audit_rho > rho_new),
                     "frozen": {k: frozen_rows[lab][k] for k in
                                ("rho_fm_full", "delta_rho", "p_mc", "bh_q", "shared_null_significant")},
                     "frozen_rho_gt_base": bool(audit_rho > rho_old)})
    out["contrast_brain"] = {"rho_baseline": res_new["rho_baseline"], "rows": rows,
                             "n_shared_sig": sum(r["shared_sig"] for r in rows),
                             "n_beats_observed": sum(r["delta"] > 0 for r in rows)}
    out["gate_changes"] = [r["model_label"] for r in rows if r["rho_gt_base"] != r["frozen_rho_gt_base"]]
    out["inputs"] = {"fixed_panel_audit_v2_sha256": common.sha256_file(common.AUDIT_JSON)}
    common.write_json(common.OUT_DIR / "brain_baseline_fix.json", out)
    log("brain baseline rho", rho_old, "->", rho_new, "| gate changes:", out["gate_changes"])
    for r in rows:
        log(f"  {r['model_label']:22s} rho {r['rho_fm']:+.5f} (audit {r['audit_rho']:+.5f}) "
            f"delta {r['delta']:+.5f} q {r['q']:.3f} sig {r['shared_sig']} | frozen q {r['frozen']['bh_q']:.3f}")


if __name__ == "__main__":
    main()
