#!/usr/bin/env python
"""
PROVENANCE RECORD — exploratory matplotlib-based figures (2024).

The authoritative generator is paper/make_figs.R (tikzDevice, Times-compatible fonts).
This file is preserved for provenance but is NOT used in the current build pipeline.

Historical context: early prototype figures for PAPER_OUTLINE.md, data-driven from
results/v2/*.json. Palette/form choices follow the dataviz skill: categorical hues
in fixed order (never cycled), sequential = one hue light->dark, dumbbell = one hue
two shades for before/after, thin marks, recessive gridlines, direct labels over
legends where the series count is small.

Original run command: python3 make_figures.py   (wrote figs/*.pdf + figs/*.png)
Current pipeline: python3 make_panel_data.py && Rscript make_figs.R
"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "results" / "v2"
FIGS = Path(__file__).resolve().parent / "figs"
FIGS.mkdir(exist_ok=True)

def load(name): return json.load(open(R / name))

# ---- validated palette (dataviz skill, references/palette.md) ----
CAT = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"]
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"]
INK = "#0b0b0b"; INK2 = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"; BASE = "#c3c2b7"
SURFACE = "#fcfcfb"
STATUS_CRIT = "#d03b3b"; STATUS_GOOD = "#0ca30c"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.edgecolor": BASE, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "svg.fonttype": "none",
})

def save(fig, name):
    fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("wrote", name)

# ============================================================ Fig 1 — truth construct
def fig1():
    svd = json.load(open(R / "svd_spectrum.json"))
    tri = {
        ("brain", "fibroblast"): load("cross_tissue_atac_v2.json")["truth_spearman"],
        ("brain", "PBMC"): load("cross_tissue_brain_vs_pbmc.json")["truth_spearman"],
        ("PBMC", "fibroblast"): load("cross_tissue_pbmc_vs_fibroblast.json")["truth_spearman"],
    }
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[1.1, 1]))

    # Panel A: SVD energy spectrum (non-degeneracy) — sequential magnitude, one hue
    x = np.arange(1, 11)
    w = 0.38
    axA.bar(x - w/2, svd["GSE174367"]["top10_energy"], width=w, color=SEQ_BLUE[3], label="Brain")
    axA.bar(x + w/2, svd["PBMC10k"]["top10_energy"], width=w, color=CAT[1], label="PBMC")
    axA.set_xlabel("singular value rank"); axA.set_ylabel("energy fraction")
    axA.set_title("A  Truth graph is high-rank\n(not motif-scan collapse)", loc="left", fontsize=9, color=INK)
    axA.set_xticks(x); axA.legend(frameon=False, loc="upper right", fontsize=7.5)
    axA.spines[["top", "right"]].set_visible(False)
    axA.annotate("rank-1 collapse\nwould put ~100%\nof energy here",
                xy=(1, svd["GSE174367"]["top10_energy"][0]), xytext=(3.2, 0.11),
                fontsize=6.8, color=MUTED, ha="left",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))

    # Panel B: 3-way cross-tissue truth reproducibility — sequential heatmap
    labs = ["Brain", "PBMC", "Fibroblast"]
    M = np.array([[1, tri[("brain","PBMC")], tri[("brain","fibroblast")]],
                  [tri[("brain","PBMC")], 1, tri[("PBMC","fibroblast")]],
                  [tri[("brain","fibroblast")], tri[("PBMC","fibroblast")], 1]])
    im = axB.imshow(M, cmap=matplotlib.colors.LinearSegmentedColormap.from_list("seqblue", SEQ_BLUE),
                    vmin=0, vmax=1)
    for i in range(3):
        for j in range(3):
            txt = "1.00" if i == j else f"{M[i,j]:.2f}"
            col = "white" if M[i, j] > 0.6 else INK
            axB.text(j, i, txt, ha="center", va="center", fontsize=9, color=col)
    axB.set_xticks(range(3)); axB.set_yticks(range(3))
    axB.set_xticklabels(labs, fontsize=8); axB.set_yticklabels(labs, fontsize=8)
    axB.set_title("B  Regulatory truth replicates\nacross 3 independent tissues", loc="left", fontsize=9, color=INK)
    for s in axB.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.06)
    cb.ax.tick_params(labelsize=7); cb.set_label("truth-graph Spearman ρ", fontsize=7.5, color=INK2)
    fig.suptitle("Sequence-grounded regulatory truth: valid and reproducible", fontsize=10.5, y=1.04, x=0.02, ha="left")
    save(fig, "fig1_truth_construct")

# ============================================================ Fig 2 — decisive result (dumbbell)
def fig2():
    conf = load("confound_regression_v2.json")
    ko = load("ko_confound_check.json")
    scf_b = load("scf_confound_check.json")
    scf_p = load("pbmc_eval_scf_v2.json")
    uce = load("uce_confound_check.json")
    pbmc = load("pbmc_eval_v2.json")

    rows = [
        # (label, pre_confound, post_confound, killed_flag)
        ("Co-expression — brain",              conf["coexp_vs_atac_marginal"],      conf["coexp_partial_confounds_only"], False),
        ("Co-expression — PBMC (paired)",       pbmc["embed__coexp_vs_atac"],         pbmc["embed__coexp_partial_confounds"], False),
        ("Geneformer embed — brain",            conf["fm_partial_coexp_only"],        conf["fm_partial_coexp_plus_confounds"], False),
        ("Geneformer embed — PBMC (paired)",    pbmc["embed__fm_partial_coexp"],      pbmc["embed__fm_partial_coexp_confounds"], False),
        ("Geneformer attention — PBMC (paired)", pbmc["attn__fm_partial_coexp"],      pbmc["attn__fm_partial_coexp_confounds"], False),
        ("Geneformer perturbation (KO) — brain", ko["ko_raw"]["partial_coexp"],        ko["ko_raw"]["partial_coexp_confounds"], True),
        ("scFoundation — brain",                scf_b["scf_partial_coexp"],           scf_b["scf_partial_coexp_confounds"], True),
        ("scFoundation — PBMC (paired)",        scf_p["scf_partial_coexp"],           scf_p["scf_partial_coexp_confounds"], False),
        ("UCE — brain",                         uce["uce_partial_coexp"],             uce["uce_partial_coexp_confounds"], True),
    ]
    n = len(rows)
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    y = np.arange(n)[::-1]
    for yi, (label, pre, post, killed) in zip(y, rows):
        ax.plot([pre, post], [yi, yi], color=BASE, lw=1.2, zorder=1)
        ax.scatter([pre], [yi], s=42, color=SEQ_BLUE[1], edgecolor=INK, linewidth=0.4, zorder=2, label="_pre")
        ax.scatter([post], [yi], s=42, color=SEQ_BLUE[4], edgecolor=INK, linewidth=0.4, zorder=3, label="_post")
        if killed:
            ax.annotate("significant\nbefore confound\ncontrol →killed", xy=(pre, yi), xytext=(pre + 0.012, yi + 0.42),
                        fontsize=6.3, color=STATUS_CRIT, ha="left",
                        arrowprops=dict(arrowstyle="-", color=STATUS_CRIT, lw=0.5))
    ax.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlabel("partial Spearman ρ (regulation, controlling co-expression)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    h1 = plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=SEQ_BLUE[1], markeredgecolor=INK, markersize=6, label="raw / partial|coexp")
    h2 = plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=SEQ_BLUE[4], markeredgecolor=INK, markersize=6, label="+ confound control")
    ax.legend(handles=[h1, h2], frameon=False, loc="lower right", fontsize=7.5)
    ax.set_title("No FM adds regulatory signal beyond co-expression, once confounded\n"
                 "4 FMs × 2 tissues × 3 readouts — every apparent positive was a confound artifact",
                 fontsize=9.5, loc="left")
    save(fig, "fig2_decisive_result")

# ============================================================ Fig 3 — confound artifact case study
def fig3():
    ko = load("insilico_ko_v2.json")
    m = ko["mantel_ko_partial"]
    kc = load("ko_confound_check.json")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 2.9))

    x = np.linspace(m["null_mean"] - 4*m["null_sd"], m["null_mean"] + 4*m["null_sd"], 400)
    from scipy.stats import norm
    axA.fill_between(x, norm.pdf(x, m["null_mean"], m["null_sd"]), color=SEQ_BLUE[1], alpha=0.55, lw=0)
    axA.plot(x, norm.pdf(x, m["null_mean"], m["null_sd"]), color=SEQ_BLUE[3], lw=1.2)
    axA.axvline(m["observed"], color=STATUS_CRIT, lw=1.6)
    axA.annotate(f"observed = {m['observed']}\nz = {m['z']}, p = {m['p_perm']}", xy=(m["observed"], 0),
                xytext=(m["observed"] - 0.005, norm.pdf(x, m["null_mean"], m["null_sd"]).max()*0.55),
                fontsize=7.2, color=STATUS_CRIT, ha="right")
    axA.set_yticks([]); axA.set_xlabel("partial ρ (Mantel gene-label permutation null, N=1000)")
    axA.set_title("A  Before confound control:\nlooks like a real signal", loc="left", fontsize=9)
    axA.spines[["top", "right", "left"]].set_visible(False)

    labels = ["raw partial", "+ confound\ncontrol"]
    vals = [kc["ko_raw"]["partial_coexp"], kc["ko_raw"]["partial_coexp_confounds"]]
    cols = [STATUS_CRIT, MUTED]
    axB.bar(labels, vals, color=cols, width=0.55)
    axB.axhline(0, color=BASE, lw=0.8)
    for i, v in enumerate(vals): axB.text(i, v + (0.003 if v >= 0 else -0.006), f"{v:+.4f}", ha="center", fontsize=8, color=INK)
    axB.set_ylabel("partial Spearman ρ")
    axB.set_title("B  After: collapses to noise\n(gene-hubness confound, not regulation)", loc="left", fontsize=9)
    axB.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Confound-artifact case study — in-silico perturbation readout, Geneformer", fontsize=10, y=1.06, x=0.02, ha="left")
    save(fig, "fig3_confound_case_study")

# ============================================================ Fig 4 — paired-cell calibration
def fig4():
    conf = load("confound_regression_v2.json"); pbmc = load("pbmc_eval_v2.json")
    scf_b = load("scf_confound_check.json"); scf_p = load("pbmc_eval_scf_v2.json")
    groups = ["Co-expression", "Geneformer\n(embedding)", "scFoundation"]
    brain = [conf["coexp_partial_confounds_only"], conf["fm_partial_coexp_plus_confounds"], scf_b["scf_partial_coexp_confounds"]]
    pbmcv = [pbmc["embed__coexp_partial_confounds"], pbmc["embed__fm_partial_coexp_confounds"], scf_p["scf_partial_coexp_confounds"]]

    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    x = np.arange(len(groups)); w = 0.34
    b1 = ax.bar(x - w/2, brain, width=w, color=SEQ_BLUE[3], label="Brain (unpaired, cross-study)")
    b2 = ax.bar(x + w/2, pbmcv, width=w, color=CAT[1], label="PBMC Multiome (paired, same cells)")
    ax.axhline(0, color=BASE, lw=0.8)
    for bars in (b1, b2):
        for r in bars: ax.text(r.get_x()+r.get_width()/2, r.get_height()+(0.0004 if r.get_height()>=0 else -0.0012),
                               f"{r.get_height():+.4f}", ha="center", fontsize=7, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8.5)
    ax.set_ylabel("confound-controlled partial ρ")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title("Paired-cell calibration: same near-zero verdict\nwhether RNA/ATAC come from different studies or the same cells", fontsize=9.5, loc="left")
    save(fig, "fig4_paired_calibration")

# ============================================================ Fig 5 — per-type robustness
def fig5():
    brain_co = load("pertype_coexp_v2.json")["per_type"]
    brain_fm = load("pertype_fm_v2.json")["per_type"]
    brain_scf = load("pertype_fm_scf_v2.json")["per_type"]
    pbmc_co = load("pbmc_eval_v2.json")["per_type_coexp"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.4, 3.3), gridspec_kw=dict(wspace=0.45))

    typesA = [r["cell_type"] for r in brain_co] + [r["cell_type"] for r in pbmc_co]
    valsA = [r["coexp_vs_atac"] for r in brain_co] + [r["coexp_vs_atac"] for r in pbmc_co]
    colsA = [SEQ_BLUE[3]] * len(brain_co) + [CAT[1]] * len(pbmc_co)
    yA = np.arange(len(typesA))[::-1]
    axA.barh(yA, valsA, color=colsA, height=0.62)
    axA.axvline(0, color=BASE, lw=0.8)
    axA.set_yticks(yA); axA.set_yticklabels(typesA, fontsize=7.5)
    axA.set_xlabel("co-expression → truth (per cell type)")
    axA.spines[["top", "right"]].set_visible(False)
    h1 = mpatches.Patch(color=SEQ_BLUE[3], label="Brain"); h2 = mpatches.Patch(color=CAT[1], label="PBMC")
    axA.legend(handles=[h1, h2], frameon=False, fontsize=7.5, loc="lower right")
    axA.set_title("A  Trivial co-expression floor,\nconsistent per cell type (both tissues)", fontsize=9, loc="left")

    types = [r["cell_type"] for r in brain_fm if r["n"] >= 150]
    emb = {r["cell_type"]: r["emb_partial"] for r in brain_fm}
    att = {r["cell_type"]: r["attn_partial"] for r in brain_fm}
    scf = {r["cell_type"]: r["scf_partial"] for r in brain_scf}
    x = np.arange(len(types)); w = 0.26
    axB.bar(x - w, [emb[t] for t in types], width=w, color=CAT[0], label="Geneformer embed")
    axB.bar(x, [att[t] for t in types], width=w, color=CAT[2], label="Geneformer attn")
    axB.bar(x + w, [scf.get(t, 0) for t in types], width=w, color=CAT[4], label="scFoundation")
    axB.axhline(0, color=BASE, lw=0.8)
    axB.set_xticks(x); axB.set_xticklabels(types, fontsize=8)
    axB.set_ylabel("FM partial ρ | co-expr.", fontsize=8, labelpad=2)
    axB.spines[["top", "right"]].set_visible(False)
    axB.legend(frameon=False, fontsize=7, loc="lower right")
    axB.set_title("B  FM adds nothing per cell type\n(brain, 3 FM readouts)", fontsize=9, loc="left")
    fig.suptitle("Per-cell-type robustness", fontsize=10.5, y=1.05, x=0.02, ha="left")
    save(fig, "fig5_pertype_robustness")

# ============================================================ Fig 6 (supp) — readout degeneracy
def fig6():
    cross = load("crossmodal_v2.json"); pbmc = load("pbmc_eval_v2.json")
    attn = load("readout_attention_v2.json"); scf_b = load("crossmodal_scf_v2.json")
    scf_p = load("pbmc_eval_scf_v2.json"); uce = load("crossmodal_uce_v2.json")
    rows = [
        ("Geneformer\nembedding", cross["observed"]["geneformer_vs_coexp"], pbmc["embed__fm_vs_coexp"]),
        ("Geneformer\nattention", attn["observed"]["attn_sym_vs_coexp"], pbmc["attn__fm_vs_coexp"]),
        ("scFoundation", scf_b["observed"]["scf_vs_coexp"], scf_p["scf_vs_coexp"]),
        ("UCE", uce["observed"]["uce_vs_coexp"], None),
    ]
    labels = [r[0] for r in rows]; brain = [r[1] for r in rows]; pbmcv = [r[2] for r in rows]
    fig, ax = plt.subplots(figsize=(5.8, 3.0))
    x = np.arange(len(labels)); w = 0.34
    ax.bar(x - w/2, brain, width=w, color=SEQ_BLUE[3], label="Brain")
    ax.bar([xi + w/2 for xi, v in zip(x, pbmcv) if v is not None], [v for v in pbmcv if v is not None],
          width=w, color=CAT[1], label="PBMC")
    ax.axhline(0, color=BASE, lw=0.8)
    ax.axhspan(-0.25, 0, color=STATUS_CRIT, alpha=0.07)
    ax.text(len(labels)-0.5, -0.22, "readout anti-correlates\nwith co-expression\n(degenerate)", fontsize=6.6, color=STATUS_CRIT, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Spearman ρ vs co-expression")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title("Supp. — readout sanity check: does the FM graph even\ntrack the known attention≈co-expression prior?", fontsize=9, loc="left")
    save(fig, "fig6_readout_degeneracy")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print("done —", len(list(FIGS.glob("*.pdf"))), "figures in", FIGS)
