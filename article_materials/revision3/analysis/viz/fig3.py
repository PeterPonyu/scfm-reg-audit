"""Figure 3: degree conditioning, paired effects and proxy margin structure."""
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from . import style as S
from .common import (ADMITTED, FIGDATA, HIST, READOUTS, ROOT, clean, corrected_rows, fmt, group_labels,
                     hbars, legend_above, load, payload, pt, short)
from .layout import Canvas, split

SPECS = ["non_degree", "external_fibro", "external_other", "full", "outdeg_only", "indeg_only", "log_degree",
         "rank_degree", "spline_degree", "strength"]
SPEC_NAMES = ["No degree", "Fibroblast", "Other tissue", "Own linear", "Out only", "In only", "Log", "Rank", "Spline",
              "Weighted"]
CONTEXT_LETTER = {"brain": "M", "pbmc": "P"}
CONTEXT_NAME = {"M": "Mixed", "P": "PBMC", "B": "Brain", "D": "BMMC"}


class Data:
    def __init__(self):
        grid = load(HIST / "analysis/degree_iut_portable.json")
        self.rows, self.summary = grid["rows"], grid["summary"]
        self.full = [r for r in self.rows if r["spec"] == "full"]
        self.corrected = [r for r in corrected_rows() if r["null_mode"] == "restricted"]
        self.bmmc = payload("figure5_independent_validation")["final"]["rows"]
        self.additive = load(ROOT / "inputs/additive-decomposition.json")["rows"]
        self.third = load(FIGDATA / "third_tissue_original_metrics.json")
        self.split = load(FIGDATA / "degree_cell_split.json")
        assert len(self.full) == 13

    def grid_value(self, spec, tissue, model, key="rho"):
        return next(v[key] for v in self.rows if (v["spec"], v["tissue"], v["model_label"]) == (spec, tissue, model))


def panel_a(d):
    def draw(ax):
        mat = np.array([[d.grid_value(s, r["tissue"], r["model_label"]) for s in SPECS] for r in d.full])
        im = ax.imshow(mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.02, vcenter=0, vmax=0.04), aspect="auto",
                       interpolation="nearest")
        for k in range(1, 13):
            ax.axhline(k - 0.5, color="white", lw=0.5)
        for k in range(1, 10):
            ax.axvline(k - 0.5, color="white", lw=0.5)
        ax.axhline(7.5, color=S.GRAY, lw=0.9)
        row_labels = [short(r["model_label"]) for r in d.full]
        ax.set_yticks(range(13))
        ax.set_yticklabels(row_labels)
        ax.set_xticks(range(10)); ax.set_xticklabels(SPEC_NAMES, rotation=90)
        ax.tick_params(length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        counts = [f"{d.summary[s]['IUT_BY_support']}/{d.summary[s]['IUT_Holm_support']}" for s in SPECS]
        for j, t in enumerate(counts):
            ax.text(j, -0.78, t, ha="center", va="bottom", fontsize=S.FS_NOTE)
        ax.text(-0.5, -1.7, "IUT-BY / Holm selections of 13", ha="left", va="bottom", fontsize=S.FS_NOTE,
                fontweight="normal")
        cax = ax.inset_axes([0.22, -0.43, 0.56, 0.028])
        cb = ax.figure.colorbar(im, cax=cax, orientation="horizontal", ticks=[-0.02, 0, 0.02, 0.04])
        cb.ax.set_xticklabels(["−0.02", "0", "0.02", "0.04"], fontsize=S.FS_TICK)
        cb.ax.tick_params(length=2, width=0.5, pad=1.5)
        cb.outline.set_linewidth(0.4)
        cb.set_label("Partial ρ", fontsize=S.FS_LABEL, labelpad=1.5)
        group_labels(ax, [("Mixed", 0, 7), ("PBMC", 8, 12)], row_labels)
        return {"rows": [f"{CONTEXT_LETTER[r['tissue']]} {short(r['model_label'])}" for r in d.full],
                "specifications": SPEC_NAMES, "rho": mat.tolist(), "iut_by_holm_selected": counts}
    return draw


def paired_rows(d):
    out = []
    for r in d.full:
        none = d.grid_value("non_degree", r["tissue"], r["model_label"])
        out.append({"context": CONTEXT_LETTER[r["tissue"]], "readout": r["model_label"], "none": none, "own": r["rho"]})
    for m in ADMITTED:
        out.append({"context": "B", "readout": m, **{deg: next(x["observed"] for x in d.corrected
                                                                if x["readout"] == m and x["degree_mode"] == deg)
                                                      for deg in ("none", "own")}})
    for m in READOUTS:
        if any(x["readout"] == m for x in d.bmmc):
            out.append({"context": "D", "readout": m, **{deg: next(x["observed"] for x in d.bmmc if x["readout"] == m
                        and x["degree_mode"] == deg and x["null_mode"] == "restricted") for deg in ("none", "own")}})
    return out


def panel_b(d):
    def draw(ax):
        rows = paired_rows(d)
        for i, r in enumerate(rows):
            ax.plot([r["none"], r["own"]], [i, i], color="#C3C8CC", lw=0.8)
            ax.scatter(r["none"], i - 0.1, c=S.ORANGE, marker="^", s=17, zorder=3, lw=0)
            ax.scatter(r["own"], i + 0.1, c=S.TEAL, marker="s", s=15, zorder=3, lw=0)
        counts = {k: sum(r["context"] == k for r in rows) for k in "MPBD"}
        assert counts == {"M": 8, "P": 5, "B": 4, "D": 5}
        for y in (7.5, 12.5, 16.5):
            ax.axhline(y, color=S.RULE, lw=0.7)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(ylim=(21.7, -0.7), xlim=(-0.03, 0.03), xticks=[-0.02, 0, 0.02], xlabel="Partial Spearman ρ")
        ax.set_xticklabels(["−0.02", "0", "0.02"])
        labels = [short(r["readout"]) for r in rows]
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(labels)
        clean(ax)
        group_labels(ax, [("Mixed", 0, 7), ("PBMC", 8, 12), ("Brain", 13, 16), ("BMMC", 17, 21)], labels)
        legend_above(ax, [pt("^", S.ORANGE, "No degree"), pt("s", S.TEAL, "Own degree")])
        return {"rows": [{"context": CONTEXT_NAME[r["context"]], "readout": short(r["readout"]), "none": r["none"],
                          "own": r["own"]} for r in rows]}
    return draw


def panel_c(d):
    def draw(ax):
        ratios = []
        for r in d.full:
            none_sd = d.grid_value("non_degree", r["tissue"], r["model_label"], key="sdM")
            ratio = none_sd / r["sdM"]
            group = 2 if r["model_label"] == "random_init_floor" else (0 if r["tissue"] == "brain" else 1)
            ax.scatter(ratio, group, s=20, color=[S.BLUE, S.TEAL, S.GRAY][group], zorder=3, lw=0)
            ratios.append({"tissue": r["tissue"], "readout": r["model_label"], "ratio": ratio})
        ax.axvline(1, color=S.GRAY, ls=":", lw=0.8)
        ax.set(ylim=(2.6, -0.6), xlim=(0, 17.5), xticks=[1, 6, 11, 16], xlabel="Null SD ratio (none / own)")
        ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["M trained", "P trained", "Random GF"])
        clean(ax)
        return {"ratios": ratios}
    return draw


def panel_d(d):
    def draw(ax):
        names = {("GSE174367", "PBMC10k"): "Brain–PBMC", ("GSE174367", "GSE206767"): "Brain–fibro.",
                 ("PBMC10k", "GSE206767"): "PBMC–fibro."}
        rec = []
        for i, r in enumerate(d.additive):
            lab = names[tuple(r["pair"])]
            obs, res = r["observed_spearman"], r["residual_spearman_after_own_additive_fits"]
            ax.plot([obs, res], [i, i], color=S.GRAY, lw=0.8)
            ax.scatter(obs, i, c=S.BLUE, s=22, zorder=3, lw=0)
            ax.scatter(res, i, c=S.TEAL, marker="s", s=20, zorder=3, lw=0)
            ax.text(0.585, i, fmt(r["fraction_explained_by_additive_marginals"], 3), va="center", fontsize=S.FS_NOTE)
            rec.append({"pair": lab, "observed": obs, "residual": res,
                        "additive_ratio": r["fraction_explained_by_additive_marginals"]})
        ax.set(ylim=(2.6, -0.6), xlim=(0.40, 0.67), xticks=[0.4, 0.5, 0.6], xlabel="Observed / residual ρ  |  ratio")
        ax.set_yticks(range(3)); ax.set_yticklabels([names[tuple(r["pair"])] for r in d.additive])
        clean(ax)
        legend_above(ax, [pt("o", S.BLUE, "Observed"), pt("s", S.TEAL, "Residual")])
        return {"rows": rec}
    return draw


def panel_e(d):
    def draw(ax):
        for (name, vals), color in zip(d.third["degrees"].items(), (S.BLUE, S.TEAL, S.ORANGE)):
            v = np.sort(vals)
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=color, lw=1.0,
                    label=name.replace(" mix", ""))
        ax.set(xlim=(0, 1200), ylim=(0, 1.02), xticks=[0, 600, 1200], yticks=[0, 0.5, 1],
               xlabel="Supported targets per TF", ylabel="Cumulative fraction")
        ax.set_xticklabels(["0", "600", "1,200"])
        ax.legend(loc="lower right", handlelength=1.2, labelspacing=0.25, borderpad=0.2)
        clean(ax, grid="both")
        degrees = d.third["degrees"]
        return {"contexts": list(degrees), "n_tfs": [len(v) for v in degrees.values()],
                "median_targets_per_tf": [float(np.median(v)) for v in degrees.values()],
                "sorted_targets_per_tf": {k: sorted(int(x) for x in v) for k, v in degrees.items()}}
    return draw


def panel_f(d):
    def draw(ax):
        items = list(d.third["linked_peaks"].items())
        fr = [n / total for _, (n, total) in items]
        hbars(ax, [k for k, _ in items], fr, [S.BLUE, S.TEAL, S.ORANGE], xlim=(0, 0.1), xticks=[0, 0.04, 0.08],
              xticklabels=["0", "0.04", "0.08"], xlabel="Linked / accessible peaks", height=0.58,
              annotate=[f"{n:,}" for _, (n, _) in items], annotate_dx=0.002)
        return {"contexts": [k for k, _ in items], "linked": [v[0] for _, v in items], "accessible": [v[1] for _, v in items]}
    return draw


def panel_g(d):
    def draw(ax):
        items = list(d.third["overlap"].items())
        pretty = {"B only": "Brain only", "P only": "PBMC only", "F only": "Fibroblast only", "B + P": "Brain + PBMC",
                  "B + F": "Brain + fibro.", "P + F": "PBMC + fibro.", "All three": "All three"}
        hbars(ax, [pretty[k] for k, _ in items], [v for _, v in items], S.TEAL, xlim=(0, 56000), xticks=[0, 20000, 40000],
              xticklabels=["0", "20k", "40k"], xlabel="Supported TF–target edges", height=0.62,
              annotate=[f"{v:,}" for _, v in items], annotate_dx=800)
        return {"sets": [k for k, _ in items], "edges": [v for _, v in items]}
    return draw


def panel_h(d):
    def draw(ax):
        rows = d.split["leave_edge_out"]
        for i, r in enumerate(rows):
            ax.plot([0, r["difference"] * 1e4], [i, i], color="#C3C8CC", lw=0.8)
            ax.scatter(r["difference"] * 1e4, i, c=S.TEAL, s=17, zorder=3, lw=0)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        labels = [short(r["readout"]) for r in rows]
        ax.set(ylim=(len(rows) - 0.4, -0.6), xlim=(-8.5, 10.5), xticks=[-6, 0, 6], xlabel=r"Δρ ($\times10^{-4}$)")
        ax.set_xticklabels(["−6", "0", "6"])
        ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels)
        clean(ax)
        n_mixed = sum(r["tissue"] == "brain" for r in rows)
        group_labels(ax, [("Mixed", 0, n_mixed - 1), ("PBMC", n_mixed, len(rows) - 1)], labels)
        return {"rows": [{"tissue": r["tissue"], "readout": short(r["readout"]), "own_rho": r["own_rho"],
                          "leave_edge_out_rho": r["leave_edge_out_rho"], "difference": r["difference"]} for r in rows]}
    return draw


def build(out_height=156.0):
    d = Data()
    c = Canvas("Figure3", out_height)
    r1, r2, r3 = split(c.rect, "v", [78.0, 27.0, 46.0])
    a_cell, b_cell = split(r1, "h", [1.62, 1.0])
    c_cell, d_cell, f_cell = split(r2, "h", [1.0, 1.1, 1.0])
    e_cell, g_cell, h_cell = split(r3, "h", [1.0, 1.1, 1.0])
    c.panel(a_cell, "A", "Ten degree designs", panel_a(d))
    c.panel(b_cell, "B", "Paired no / own degree", panel_b(d))
    c.panel(c_cell, "C", "Gene-label null width", panel_c(d), align_y="r2")
    c.panel(d_cell, "D", "Additive margins", panel_d(d), align_y="r2")
    c.panel(f_cell, "E", "Linked-peak admission", panel_f(d), align_y="r2")
    c.panel(e_cell, "F", "TF out-degree", panel_e(d), align_y="r3")
    c.panel(g_cell, "G", "Edge support by context set", panel_g(d), align_y="r3")
    c.panel(h_cell, "H", "Leave-edge-out", panel_h(d), align_y="r3")
    return c, {}
