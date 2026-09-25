#!/usr/bin/env python
"""Write the Appendix B tables of the ARRAY-D-26-05613 revision from the result JSONs.

Every number in paper/revision_array_05613_R1/source/appendix_b_tables.tex comes from
results/v2/ (frozen) or results/revision_05613/ (revision analyses). A missing input
becomes \\PENDING{...}, which stops the LaTeX build instead of printing a gap.
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

R = common.OUT_DIR
_PAPER_OUT = Path(common.fpa.ROOT) / "paper/revision_array_05613_R1/source"
OUT_TEX = Path(os.environ.get("SCREG_APPENDIX_OUT") or (
    _PAPER_OUT / "appendix_b_tables.tex" if _PAPER_OUT.is_dir() else R / "appendix_b_tables.tex"))
LABEL = {
    "geneformer_embed": "Geneformer embed", "geneformer_attn": "Geneformer attention",
    "scFoundation_encoder": "scFoundation encoder", "UCE_encoder": "UCE encoder",
    "scGPT_encoder": "scGPT encoder", "random_init_floor": "Random-init floor",
    "geneformer_ko_raw": "Geneformer KO", "geneformer_ko_posctrl": "Artifact-corrected KO",
}
SPEC_LABEL = {
    "non_degree": "none (non-degree)", "full": "linear out + in (full)",
    "outdeg_only": "linear out only", "indeg_only": "linear in only",
    "log_degree": r"$\log(1+\cdot)$ out + in", "rank_degree": "ranked out + in",
    "spline_degree": "spline out + in", "strength": "weighted out + in",
    "external_fibro": "fibroblast-mix proxy degree", "external_other": "other-tissue proxy degree",
}


def load(name):
    p = R / name
    return json.loads(p.read_text()) if p.exists() else None


def pend(what):
    return "\\PENDING{" + what + "}"


def f(x, d=3):
    if x is None:
        return "--"
    return f"${x:.{d}f}$"


def table(label, caption, colspec, header, rows, resize=True):
    body = "\n".join(" & ".join(r) + r" \\" for r in rows)
    inner = (f"\\begin{{tabular}}{{{colspec}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n"
             "\\bottomrule\n\\end{tabular}")
    wrap = f"\\fitfig{{{inner}}}" if resize else inner
    return (f"\\begin{{table}}[!htbp]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            f"{{\\small {wrap}}}\n\\end{{table}}\n")


def t_panel():
    d = load("panel_cascade.json")
    if d is None or not d.get("reproduces_frozen_manifest"):
        return pend("panel cascade")
    rows = [[s["step"].replace(">=", r"$\ge$").replace("TSS-2 kb", r"TSS$-$2\,kb").replace("%", r"\%"),
             f"{s['genes']:,}".replace(",", "{,}"), str(s["tfs"])] for s in d["steps"]]
    cap = ("Construction of the frozen 1,200-gene panel. Filters are applied in the order listed; "
           "counts are the genes, and the JASPAR TFs among them, remaining after each step. "
           "Re-applying the rules reproduces the frozen manifest and its SHA-256 "
           f"(\\texttt{{{d['manifest_sha256'][:12]}}}\\ldots). The cap of 1,200 is set by the scGPT "
           "input length.")
    return table("tab:panel", cap, "p{0.66\\linewidth}rr", "Step & Genes & TFs", rows)


def t_motif():
    d = load("motif_expectation.json")
    if d is None:
        return pend("motif expectation")
    rows = [
        ["JASPAR 2024 CORE vertebrate matrices, total / used", f"{d['jaspar_motifs_total']} / {d['motifs_used']}"],
        ["Panel TFs with a matrix / with more than one / heterodimer matrices",
         f"{d['panel_tfs_with_motif']} / {d['tfs_with_multiple_motifs']} / {d['heterodimer_motifs_used']}"],
        ["Matrix width, bp (min / median / max)",
         f"{d['motif_width_bp']['min']} / {d['motif_width_bp']['median']:.0f} / {d['motif_width_bp']['max']}"],
        [r"First-order expectation $554\times2\times500\times10^{-5}$", f"{d['expected_first_order']:.2f}"],
        ["Expected hit positions per window (widths accounted for)", f"{d['expected_hit_positions']:.2f}"],
        ["Expected distinct motifs hit per window", f"{d['expected_distinct_motifs_per_peak']:.2f}"],
        ["Expected distinct TFs hit per window (heterodimers count twice)", f"{d['expected_distinct_tfs_per_peak']:.2f}"],
    ]
    for tissue, lab in (("brain", "brain"), ("pbmc", "PBMC"), ("fibroblast_mix", "fibroblast mix")):
        o = d["observed"][tissue]
        rows.append([f"Observed distinct motifs per linked peak, {lab} "
                     f"({o['distinct_motif_hits']:,}/{o['linked_peaks']:,})".replace(",", "{,}"),
                     f"{o['per_peak']:.2f} ({o['ratio_to_expected']:.1f}$\\times$)"])
    cap = ("Motif scan accounting behind the expected 5.5 random hits per peak. Windows are 500\\,bp "
           "around each linked peak midpoint, both strands, $p<10^{-5}$ against a uniform "
           "background. A motif counts once per window however many positions score; a TF counts "
           "once per peak after an OR over its matrices.")
    return table("tab:motif", cap, "p{0.72\\linewidth}r", "Quantity & Value", rows)


def t_seeds():
    doc = json.loads(common.AUDIT_JSON.read_text())
    rows = []
    for tissue in ("brain", "pbmc"):
        seen = {}
        for fam in ("primary_family", "sensitivity_family"):
            for r in doc["pooled"][tissue][fam]["rows"]:
                key = (r["mantel"]["seed"], r["degree_preserving"]["seed"], r["confound_spec"])
                seen.setdefault(key, []).append(LABEL[r["model_label"]])
        for (sm, sd, spec), labs in seen.items():
            who = "all eight rows" if len(labs) == 8 else ", ".join(labs)
            rows.append([f"{'Brain' if tissue == 'brain' else 'PBMC'}, {spec.replace('_', '-')}",
                         who, f"\\texttt{{{sm}}}", f"\\texttt{{{sd}}}"])
    for t, sm, sd in (("Brain", 2026073101, 2026073102), ("PBMC", 2026073103, 2026073104)):
        rows.append([f"{t}, co-expression baseline", "baseline", f"\\texttt{{{sm}}}", f"\\texttt{{{sd}}}"])
    rows.append(["Brain / PBMC, FM versus baseline contrast", "shared gene-label draws",
                 "\\texttt{2026081001} / \\texttt{2026081002}", "--"])
    rows.append(["PBMC, TF-disjoint probe", "families (sorted index $k$)",
                 "\\texttt{20260730000}$+k$", "sign flip \\texttt{20260731}"])
    rows.append(["Revision: degree grid; power simulation", "Appendix B",
                 "\\texttt{SeedSequence([20260924,3])}", "\\texttt{SeedSequence([20260924,2])}"])
    rows.append(["Revision: attention-omission guard", "brain rows",
                 "\\texttt{SeedSequence([20260924,4])}", "(same sequence)"])
    cap = ("Seeds of every randomization test. Audit seeds derive from root \\texttt{20260724}: "
           "\\texttt{SeedSequence(20260724)} spawns one child per tissue, and each child spawns four "
           "integers (gene-label full, gene-label non-degree, row-shuffle full, row-shuffle "
           "non-degree); PBMC scGPT and UCE use \\texttt{SeedSequence([20260724,20260726,1])} and "
           "\\texttt{[20260724,20260729,1]}. Each audit test draws replicate $b$ from the $b$-th child "
           "of \\texttt{SeedSequence(seed)}; baseline and probe tests draw sequentially from "
           "\\texttt{default\\_rng(seed)}. All are recorded in the release JSON files.")
    return table("tab:seeds", cap, "p{0.24\\linewidth}p{0.26\\linewidth}ll",
                 "Tests & Rows & Gene-label ($M$) & Row-shuffle ($D$)", rows)


def t_multiplicity():
    d = load("multiplicity_robustness.json")
    if d is None:
        return pend("multiplicity table")
    fz = d["frozen_n999"]["families"]
    deep = d.get("resolution_n9999", {}).get("families")
    rows = []
    for tissue in ("brain", "pbmc"):
        fam = fz[f"{tissue}_full"]
        for i, lab in enumerate(fam["labels"]):
            r = [("Brain" if tissue == "brain" else "PBMC"), LABEL[lab]]
            for m in ("BH", "BY", "Holm"):
                r.append(f"{fam['q'][m]['mantel'][i]:.3f} / {fam['q'][m]['degree'][i]:.3f}")
            r.append(f"{fam['q']['IUT-BY'][i]:.3f}")
            if deep:
                dfam = deep[f"{tissue}_full"]
                r.append(f"{dfam['q']['WY']['mantel'][i]:.4f} / {dfam['q']['WY']['degree'][i]:.4f}")
            else:
                r.append(pend("Westfall-Young"))
            r.append("yes" if fam["support"]["BH"][i] else "no")
            rows.append(r)
    counts = d["frozen_n999"]["support_counts"]
    extra = ""
    if deep:
        c2 = d["resolution_n9999"]["support_counts"]
        extra = (f" Westfall--Young (joint null, $N=9{{,}}999$): {c2['WY']['full']}/13 full, "
                 f"{c2['WY']['non_degree']}/13 non-degree.")
    cap = ("Dual-null Support under multiplicity procedures that do not assume independence or "
           "positive dependence, full specification. Entries are adjusted $q_M / q_D$ within the "
           "frozen families (BH, Benjamini--Yekutieli, Holm), the intersection--union $q$ from "
           "Benjamini--Yekutieli on $\\max(p_M,p_D)$, and Westfall--Young step-down min$P$ adjusted "
           "$p_M / p_D$ from the joint randomization distribution. Support counts under BH, BY, Holm, "
           f"IUT-BH and IUT-BY: {counts['BH']['full']}, {counts['BY']['full']}, "
           f"{counts['Holm']['full']}, {counts['IUT-BH']['full']} and {counts['IUT-BY']['full']} of 13 "
           f"(full); {counts['BH']['non_degree']}, {counts['BY']['non_degree']}, "
           f"{counts['Holm']['non_degree']}, {counts['IUT-BH']['non_degree']} and "
           f"{counts['IUT-BY']['non_degree']} of 13 (non-degree).{extra}")
    return table("tab:multiplicity", cap, "llccccll",
                 r"Tissue & Readout & BH & BY & Holm & IUT-BY & WY & Support", rows)


def t_resolution():
    d = load("resolution_n9999/resolution_n9999.json")
    if d is None:
        return pend("resolution table")
    if d["prefix_mismatches"]:
        return pend("resolution prefix mismatch")
    rows = []
    by = {}
    for x in d["fm_rows"]:
        by.setdefault((x["tissue"], x["spec"], x["model_label"]), {})[x["null"]] = x
    mult = load("multiplicity_robustness.json")
    deep = (mult or {}).get("resolution_n9999", {}).get("families", {})
    for tissue in ("brain", "pbmc"):
        fam = deep.get(f"{tissue}_full")
        for key in [k for k in by if k[0] == tissue and k[1] == "full"]:
            m, dd = by[key]["mantel"], by[key]["degree"]
            i = fam["labels"].index(key[2]) if fam else None
            sup = ("yes" if fam and fam["support"]["BH"][i] else "no") if fam else pend("support")
            rows.append([("Brain" if tissue == "brain" else "PBMC"), LABEL[key[2]],
                         f"{m['frozen_p']:.3f}", f"{m['p']:.4f}", f"{dd['frozen_p']:.3f}", f"{dd['p']:.4f}", sup])
    base = {}
    for b in d["baseline_rows"]:
        base.setdefault(b["tissue"], {})[b["null"]] = b
    for t, v in base.items():
        rows.append([("Brain" if t == "brain" else "PBMC"), "Co-expression baseline",
                     f"{v['mantel']['frozen_p']:.3f}", f"{v['mantel']['p']:.4f}",
                     f"{v['degree']['frozen_p']:.3f}", f"{v['degree']['p']:.4f}", "--"])
    pr = d.get("probe")
    extra = ""
    if pr:
        cx = pr["families"]["co_expression"]
        extra = (f" TF-disjoint probe: the co-expression family's label-permutation $p$ is "
                 f"{cx['frozen_p']:.3f} at $N=999$ and {cx['p']:.4f} at $N=9{{,}}999$; sign-flip "
                 "contrast $p$ at $N=9{,}999$: " + "; ".join(
                     f"{LABEL.get(k, 'random-init floor' if k == 'random_floor' else k)} {v['p']:.4f}"
                     for k, v in pr["contrasts"].items()) + ".")
    cap = ("Resolution check: every frozen test rerun with its own seed at $N=9{,}999$ "
           "(resolution $10^{-4}$). The first 999 replicates are the frozen replicates, and the "
           "$p$-values recomputed from them equal the frozen values in every row. Support is BH "
           "dual-null Support recomputed at $N=9{,}999$ within the frozen families, full "
           "specification." + extra)
    return table("tab:resolution", cap, "llrrrrl",
                 r"Tissue & Readout & $p_M$ (999) & $p_M$ (9{,}999) & $p_D$ (999) & $p_D$ (9{,}999) & Support", rows)


def t_mde():
    d = load("mde.json")
    ps = None
    for cand in sorted(R.glob("power_sim_n999_r*/power_sim_*.json")):
        ps = json.loads(cand.read_text())
    if d is None:
        return pend("MDE table")
    rows = []
    for r in d["rows"]:
        if r["spec"] != "full":
            continue
        rows.append([("Brain" if r["tissue"] == "brain" else "PBMC"), LABEL[r["model_label"]],
                     f(r["rho"], 5), f"{r['sd_mantel_n999']:.5f}", f"{r['sd_degree_n999']:.5f}",
                     f"{r['mde_dual_a05']:.4f}", f"{r['mde_dual_a05m']:.4f}"])
    sr = d["summary_ranges"]
    nd = (f"Non-degree specification: dual MDE {sr['brain_non_degree']['mde_dual_a05'][0]:.4f}--"
          f"{sr['brain_non_degree']['mde_dual_a05'][1]:.4f} (brain) and "
          f"{sr['pbmc_non_degree']['mde_dual_a05'][0]:.4f}--{sr['pbmc_non_degree']['mde_dual_a05'][1]:.4f} "
          "(PBMC) at $a=0.05$.")
    if ps:
        sim = "; ".join(
            f"{'brain' if t == 'brain' else 'PBMC'}: empirical {v['empirical_mde_dual']['a05']:.4f} versus "
            f"analytic {v['analytic_mde_dual']['a05']:.4f} ($a=0.05$), empirical "
            f"{v['empirical_mde_dual']['a05m']:.4f} versus analytic {v['analytic_mde_dual']['a05m']:.4f} "
            f"($a=0.05/m$), false-positive rate at $\\alpha=0$ {v['false_positive_rate_at_alpha0']['a05']:.2f}"
            for t, v in ps["tissues"].items())
        sim = (f" Planted-signal simulation ({len(ps['alpha_grid'])} injected fractions $\\times$ "
               f"{ps['n_rep']} replicates, $N={ps['n_perm']}$, full specification): " + sim + ".")
    else:
        sim = " " + pend("power simulation")
    cap = ("Sensitivity of the full-specification tests. $\\sigma_M,\\sigma_D$ are the standard "
           "deviations of each row's gene-label and row-shuffle randomization distributions; the "
           "dual MDE is $(z_{1-a/2}+z_{0.80})\\max(\\sigma_M,\\sigma_D)$ at $a=0.05$ and at "
           "$a=0.05/m$. It is the smallest $|\\rho|$ the tests would call unusual with 80\\% "
           "probability on this panel, not a population power or an exclusion bound. " + nd + sim)
    return table("tab:mde", cap, "llrrrrr",
                 r"Tissue & Readout & $\rho$ & $\sigma_M$ & $\sigma_D$ & MDE ($a{=}0.05$) & MDE ($a{=}0.05/m$)", rows)


def t_grid():
    cands = sorted(R.glob("degree_grid_n999/degree_grid_n999.json"))
    if not cands:
        return pend("degree grid")
    d = json.loads(cands[0].read_text())
    if not d.get("all_frozen_observed_match"):
        return pend("degree grid frozen-rho mismatch")
    rows = []
    for spec in d["specs"]:
        s = d["summary"][spec]
        rows.append([SPEC_LABEL[spec], f"{s['support']} ({s['support_positive']}+, {s['support_negative']}$-$)",
                     f"{s['sign_agrees_with_full']}/13", str(s["d_only"]), str(s["m_only"])])
    cap = ("Degree-specification grid. Every specification adds its degree block to the same "
           "non-degree base design; terms from the tested proxy are recomputed inside each "
           "randomization, terms from another proxy stay fixed. All specifications share one "
           "randomization per replicate ($N=999$, new seeds), and the observed $\\rho$ of the two "
           "frozen specifications reproduce the frozen audit. Support is BH dual-null within "
           "tissue $\\times$ specification $\\times$ null; signs are compared with the full "
           "specification row by row.")
    return table("tab:degreegrid", cap, "p{0.34\\linewidth}cccc",
                 r"Degree block & Support (sign) & Same sign as full & D-only & M-only", rows)


def t_basefix():
    d = load("brain_baseline_fix.json")
    if d is None:
        return pend("brain baseline fix")
    b = d["baseline"]
    reso = load("resolution_n9999/resolution_n9999.json")
    old = {x["null"]: x["p"] for x in (reso or {}).get("baseline_rows", []) if x["tissue"] == "brain"}
    old9999 = f"{old['mantel']:.4f} / {old['degree']:.4f}" if len(old) == 2 else pend("frozen brain baseline at N=9,999")
    rows = [["Brain baseline $\\rho$", f"{b['rho_frozen']:.5f}", f"{b['rho_brain_covariates']:.5f}"],
            ["$p_M$ / $p_D$ ($N=999$)", f"{b['frozen']['pM']:.3f} / {b['frozen']['pD']:.3f}",
             f"{b['n999']['pM']:.3f} / {b['n999']['pD']:.3f}"],
            ["$p_M$ / $p_D$ ($N=9{,}999$)", old9999, f"{b['n9999']['pM']:.4f} / {b['n9999']['pD']:.4f}"]]
    for r in d["contrast_brain"]["rows"]:
        rows.append([f"{LABEL[r['model_label']]}: $\\Delta\\rho$, $q_\\Delta$",
                     f"{r['frozen']['delta_rho']:+.5f}, {r['frozen']['bh_q']:.3f}",
                     f"{r['delta']:+.5f}, {r['q']:.3f}"])
    cap = ("Brain co-expression baseline and FM-versus-baseline contrast with GC from the PBMC "
           "covariate cache (earlier capsule versions) and from brain peaks (this version). Gene "
           "length, detection and peak count are the same in both. With brain GC the contrast's "
           "$\\rho_{\\mathrm{FM}}$ equals the primary audit $\\rho$ in every row. Frozen seeds and $N$.")
    return table("tab:basefix", cap, "p{0.46\\linewidth}ll", "Quantity & PBMC-cache GC & Brain GC", rows)


def main():
    parts = [
        "% Generated by src/revision_05613/make_appendix_b.py; do not edit by hand.\n",
        t_panel(), t_motif(), t_seeds(), t_multiplicity(), t_resolution(), t_mde(), t_grid(), t_basefix(),
    ]
    OUT_TEX.write_text("\n".join(parts))
    n_pending = sum(p.count("\\PENDING{") for p in parts)
    print(f"wrote {OUT_TEX} ({n_pending} pending)")


if __name__ == "__main__":
    main()
