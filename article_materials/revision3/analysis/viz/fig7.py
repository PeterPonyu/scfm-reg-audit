"""Figure 7: descriptive cell-type effects, TF-disjoint probe and construct-only transfer."""
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from statsmodels.stats.multitest import multipletests

from . import style as S
from .common import FIGDATA, HIST, clean, fmt, fmt_rho, hbars, heat, legend_above, load, pt, short
from .layout import Canvas, split

PROBE_NAMES = {"co_expression": "Co-expr.", "random_floor": "Random GF"}


def public(name):
    return load(HIST / "public_results" / f"{name}_v2.public.json")


def lab(m):
    return PROBE_NAMES.get(m, short(m))


class Data:
    def __init__(self):
        self.fixed = load(HIST / "analysis/inputs/v2/fixed_panel_audit_v2.json")
        self.split = load(FIGDATA / "degree_cell_split.json")
        self.stats = public("tf_probe_pair_stats")
        self.eval = public("tf_probe_pair_eval")
        self.sens = public("tf_probe_pair_sensitivity")
        self.transfer = load(FIGDATA / "construct_transfer_v2.json")
        order = ["geneformer_embed", "geneformer_attn", "random_floor", "scGPT_encoder", "UCE_encoder", "scFoundation_encoder",
                 "co_expression"]
        present = list(self.stats["families"])
        self.families = [m for m in order if m in present]
        assert sorted(self.families) == sorted(present) and len(self.families) == 6


def pertype_panel(d, tissue):
    def draw(ax):
        per = d.fixed["per_cell_type"][tissue]
        rows = [r for r in per["full_confound"]["rows"] if r["row_type"] == "pertype_fm"]
        types = list(dict.fromkeys(r["cell_type"] for r in rows))
        models = list(dict.fromkeys(r["model_label"] for r in rows))
        other = next(k for k in per if k != "full_confound")
        columns = [(deg, m) for deg in ("full_confound", other) for m in models]
        cells = {r["cell_type"]: r["n_cells"] for r in per["full_confound"]["rows"]}
        mat = np.array([[next(r["observed_partial_rho"] for r in per[deg]["rows"] if r["row_type"] == "pertype_fm"
                              and r["cell_type"] == ty and r["model_label"] == m) for deg, m in columns] for ty in types])
        labels = [f"{ty} ({cells[ty]:,})" for ty in types]
        cols = [short(m) + (" O" if deg == "full_confound" else " N") for deg, m in columns]
        heat(ax, mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.06, vcenter=0, vmax=0.06), fmt=fmt_rho,
             dark=lambda v: abs(v) > 0.03, row_labels=labels, col_labels=cols, col_rotation=90)
        for k in range(len(models), len(columns), len(models)):
            ax.axvline(k - 0.5, color=S.GRAY, lw=0.9)
        return {"cell_types": types, "cells": [cells[t] for t in types], "columns": cols, "partial_rho": mat.tolist()}
    return draw


def panel_c(d):
    def draw(ax):
        diffs = np.array([r["difference"] for r in d.split["rows"]])
        unchanged = int((np.abs(diffs) < 1e-12).sum())
        assert (unchanged, len(diffs), len(d.split["splits"])) == (100, 100, 10)
        hbars(ax, ["Unchanged", "Changed"], [unchanged, len(diffs) - unchanged], [S.TEAL, S.GRAY], xlim=(0, 150),
              xticks=[0, 50, 100], xlabel="Readout / half comparisons", height=0.55,
              annotate=[f"{unchanged}/{len(diffs)}", "0"], annotate_dx=3)
        return {"unchanged": unchanged, "changed": 0, "comparisons": len(diffs), "cell_splits": 10}
    return draw


def panel_e(d):
    def draw(ax):
        for i, m in enumerate(d.families):
            ax.scatter(d.stats["families"][m]["adjusted_rho_mean"], i, c=S.GRAY if m == "co_expression" else S.BLUE, s=20, zorder=3, lw=0)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(ylim=(5.6, -0.6), xlim=(-0.026, 0.011), xticks=[-0.02, -0.01, 0], xlabel="Mean adjusted ρ, 135 TFs")
        ax.set_xticklabels(["−0.02", "−0.01", "0"])
        ax.set_yticks(range(6)); ax.set_yticklabels([lab(m) for m in d.families])
        clean(ax)
        return {"families": [lab(m) for m in d.families],
                "adjusted_rho_mean": [d.stats["families"][m]["adjusted_rho_mean"] for m in d.families],
                "train_tfs": d.stats["n_train_tfs"], "test_tfs": d.stats["n_test_tfs"]}
    return draw


def panel_f(d):
    def draw(ax):
        contr = d.stats["contrasts_vs_baseline"]
        real = [m for m in d.families if m not in ("co_expression", "random_floor")]
        four = dict(zip(real, multipletests([contr[m]["signflip_p"] for m in real], method="fdr_bh")[1]))
        assert abs(four["UCE_encoder"] - 0.052) < 1e-3
        for i, m in enumerate(d.families):
            if m == "co_expression":
                continue
            ax.scatter(contr[m]["signflip_q"], i - 0.17, c=S.BLUE, s=19, zorder=3, lw=0)
            if m in four:
                ax.scatter(four[m], i + 0.17, marker="s", c=S.ORANGE, s=17, zorder=3, lw=0)
        ax.axvline(0.05, color=S.RED, ls="--", lw=0.7)
        ax.set(ylim=(5.6, -0.6), xscale="log", xlim=(0.02, 1.3), xticks=[0.05, 1], xlabel="BH q vs co-expr.")
        ax.set_xticklabels(["0.05", "1"])
        ax.set_yticks(range(6)); ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        ax.minorticks_off(); clean(ax)
        legend_above(ax, [pt("o", S.BLUE, "All 5"), pt("s", S.ORANGE, "4 trained")])
        return {"contrasts": [lab(m) for m in d.families if m != "co_expression"],
                "bh_q_five_contrasts": [contr[m]["signflip_q"] for m in d.families if m != "co_expression"],
                "bh_q_four_trained": four}
    return draw


def panel_g(d):
    def draw(ax):
        sub = ["none", "detv_only", "atac_construction", "full"]
        mat = np.array([[d.sens["grid"][m][k]["adjusted_rho_mean"] for k in sub] for m in d.families])
        heat(ax, mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.10, vcenter=0, vmax=0.10), fmt=fmt_rho,
             dark=lambda v: abs(v) > 0.05, row_labels=None,
             col_labels=["None", "Detect.", "Constr.", "Full"], col_rotation=90)
        return {"families": [lab(m) for m in d.families], "designs": ["None", "Detection", "Construction", "Full"],
                "adjusted_rho_mean": mat.tolist()}
    return draw


def panel_h(d):
    def draw(ax):
        for i, m in enumerate(d.families):
            r = d.eval["results"]["edge_only"][m]
            ax.errorbar(r["adjusted_rho_mean"], i, xerr=r["adjusted_rho_std"], fmt="o", color=S.BLUE, ms=3.2, capsize=1.8, lw=0.8)
        ax.axvline(0, color=S.GRAY, lw=0.7)
        ax.set(ylim=(5.6, -0.6), xlim=(-0.095, 0.085), xticks=[-0.08, 0, 0.08], xlabel="Mean ± 1 SD")
        ax.set_xticklabels(["−0.08", "0", "0.08"])
        ax.set_yticks(range(6)); ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        clean(ax)
        return {"families": [lab(m) for m in d.families],
                "mean": [d.eval["results"]["edge_only"][m]["adjusted_rho_mean"] for m in d.families],
                "sd": [d.eval["results"]["edge_only"][m]["adjusted_rho_std"] for m in d.families]}
    return draw


def panel_i(d):
    def draw(ax):
        t = d.transfer
        mat = np.array([t["observed_rho"][o] for o in t["overlays"]])
        heat(ax, mat, cmap="Blues", vmin=0, vmax=1, fmt="{:.3f}", dark=lambda v: v > 0.6, row_labels=t["overlays"],
             col_labels=t["references"], col_rotation=0)
        return {"overlays": t["overlays"], "references": t["references"], "observed_spearman": mat.tolist()}
    return draw


def panel_j(d):
    def draw(ax):
        t = d.transfer
        cr, bc = t["collectri"], t["bmmc_coverage"]
        items = [("Prior TFs", cr["tf_hit"]), ("Prior genes", cr["gene_hit"]), ("Prior edges", cr["edge_keep"]),
                 ("BMMC genes", bc["genes"]), ("BMMC TFs", bc["tf"])]
        vals = [v["numerator"] / v["denominator"] for _, v in items]
        ax.barh(range(5), vals, height=0.58, color=[S.BLUE] * 3 + [S.TEAL] * 2)
        ax.axvline(bc["predeclared_gate"], color=S.RED, ls="--", lw=0.8)
        for i, v in enumerate(vals):
            ax.text(v + 0.02, i, fmt(v, 3), va="center", fontsize=S.FS_NOTE)
        ax.set(ylim=(4.55, -0.55), xlim=(0, 1.22), xticks=[0, 0.5, 0.9], xlabel="Frozen-panel coverage")
        ax.set_yticks(range(5)); ax.set_yticklabels([k for k, _ in items])
        clean(ax)
        return {"items": [k for k, _ in items], "fraction": vals, "gate": bc["predeclared_gate"]}
    return draw


def build(out_height=122.0):
    d = Data()
    c = Canvas("Figure7", out_height)
    r1, r2, r3 = split(c.rect, "v", [1.25, 1.2, 0.92])
    a_cell, b_cell = split(r1, "h", [1.0, 0.9])
    e_cell, f_cell, g_cell, h_cell = split(r2, "h", [1.0, 0.6, 1.0, 0.7])
    c_cell, i_cell, j_cell = split(r3, "h", [1.0, 1.0, 1.0])
    c.panel(a_cell, "A", "Mixed per-type effects", pertype_panel(d, "brain"), align_y="r1")
    c.panel(b_cell, "B", "PBMC per-type effects", pertype_panel(d, "pbmc"), align_y="r1")
    c.panel(e_cell, "C", "Adjusted probe", panel_e(d), align_y="r2")
    c.panel(f_cell, "D", "Probe contrasts", panel_f(d), align_y="r2")
    c.panel(g_cell, "E", "Covariate sensitivity", panel_g(d), align_y="r2")
    c.panel(h_cell, "F", "Per-TF spread", panel_h(d), align_y="r2")
    c.panel(c_cell, "G", "Cell-half degree", panel_c(d), align_y="r3")
    c.panel(i_cell, "H", "Construct-only transfer", panel_i(d), align_y="r3")
    c.panel(j_cell, "I", "Prior / overlay coverage", panel_j(d), align_y="r3")
    return c, {}
