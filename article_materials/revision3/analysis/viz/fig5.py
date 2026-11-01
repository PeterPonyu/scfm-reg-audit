"""Figure 5: component-test validity, restricted transformations and synthetic sensitivity."""
import numpy as np
from scipy.stats import beta

from . import style as S
from .common import (FIGDATA, clean, fmt, legend_above, load, payload, pt, short, group_labels)
from .layout import Canvas, split

SCENARIOS = [("M", S.BLUE, "o", "Gene label M"), ("D", S.ORANGE, "s", "Target shuffle D"), ("both", S.TEAL, "D", "Joint")]


class Data:
    def __init__(self):
        self.rand = payload("figure6_randomization")
        self.restricted = load(FIGDATA / "restricted_randomization.json")
        self.ladder = load(FIGDATA / "injection_ladder_v2.json")
        self.power = load(FIGDATA / "power_extension_n999_r100.json")
        self.stress = self.rand["stress"]["tissues"]
        assert self.rand["stress"]["n_replicates"] == 200


def stress_panel(d, tissue, labels):
    def draw(ax):
        rows = d.stress[tissue]["rows"]
        assert [r["scenario"] for r in rows[:2]] == ["iid", "heteroskedastic_row_column"]
        rec = []
        for j, (key, color, marker, name) in enumerate(SCENARIOS):
            for i, row in enumerate(rows[:2]):
                r = row[key]
                lo, hi = r["wilson95"]
                y = i * 1.4 + (j - 1) * 0.4
                ax.errorbar(r["rate"], y, xerr=[[r["rate"] - lo], [hi - r["rate"]]], fmt=marker, color=color, ms=3.6,
                            capsize=1.8, lw=0.8)
                ax.text(0.2, y, f"{r['count']}/200", va="center", fontsize=S.FS_NOTE)
                rec.append({"scenario": row["scenario"], "component": key, "count": r["count"], "rate": r["rate"],
                            "wilson95": r["wilson95"]})
        ax.axvline(0.05, color=S.RED, ls="--", lw=0.8)
        ax.set(ylim=(2.2, -0.7), xlim=(0, 0.275), xticks=[0, 0.05, 0.1, 0.15, 0.2], xlabel="Rejection rate")
        ax.set_xticklabels(["0", "0.05", "0.10", "0.15", "0.20"])
        ax.set_yticks([0, 1.4]); ax.set_yticklabels(["Independent\nnoise", "Shared\nshocks"] if labels else [])
        if not labels:
            ax.tick_params(axis="y", length=0)
        clean(ax)
        if tissue == "brain":
            legend_above(ax, [pt(m, col, n.replace("Gene label ", "").replace("Target shuffle ", "")) for _, col, m, n in SCENARIOS])
        return {"rows": rec}
    return draw


def panel_c(d):
    def draw(ax):
        rows = [(t, d.stress[t]["rows"][2]) for t in ("brain", "pbmc")]
        assert all(r["scenario"] == "opportunity_mean_stress" and r["M"]["count"] == r["D"]["count"] == r["both"]["count"] == 200
                   for _, r in rows)
        rec, labels, y = [], [], 0
        for t, row in rows:
            for key, color, marker, name in SCENARIOS:
                r = row[key]
                lo, hi = r["wilson95"]
                ax.errorbar(r["rate"], y, xerr=[[r["rate"] - lo], [hi - r["rate"]]], fmt=marker, color=color, ms=3.8,
                            capsize=1.8, lw=0.8)
                ax.text(lo - 0.0012, y, f"{r['count']}/200", ha="right", va="center", fontsize=S.FS_NOTE)
                labels.append({"M": "M", "D": "D", "both": "Joint"}[key])
                rec.append({"context": "Mixed" if t == "brain" else "PBMC", "component": key, "count": r["count"],
                            "wilson95": r["wilson95"]})
                y += 1
        ax.axhline(2.5, color=S.RULE, lw=0.6)
        ax.set(xlim=(0.93, 1.004), ylim=(5.6, -0.6), xticks=[0.95, 1.0], xlabel="Rejection rate")
        ax.set_xticklabels(["0.95", "1.00"])
        ax.set_yticks(range(6)); ax.set_yticklabels(labels)
        clean(ax)
        group_labels(ax, [("Mixed", 0, 2), ("PBMC", 3, 5)], labels)
        return {"rows": rec}
    return draw


def panel_d(d):
    def draw(ax):
        pairs, y = [], 0
        for t in ("brain", "pbmc"):
            fam = d.rand["multiplicity"]["frozen_n999"]["families"][t + "_full"]
            names, old = fam["labels"], fam["q"]["IUT-BY"]
            for r in d.restricted["tissues"][t]["rows"]:
                q_old, q_new = old[names.index(r["readout"])], r["q_IUT_BY"]
                color = S.RED if r["readout"] == "geneformer_embed" else S.GRAY
                ax.plot([q_old, q_new], [y, y], color=color, lw=0.9)
                ax.scatter(q_old, y - 0.12, facecolor="white", edgecolor=color, s=17, zorder=3, lw=0.9)
                ax.scatter(q_new, y + 0.12, color=color, s=17, zorder=3, lw=0)
                pairs.append({"tissue": t, "readout": r["readout"], "unrestricted": q_old, "restricted": q_new})
                if r["readout"] == "geneformer_embed":
                    ax.text(q_new * 1.25, y - 0.35, fmt(q_new, 3), ha="left", va="center", color=S.RED, fontsize=S.FS_NOTE)
                y += 1
        ax.axhline(7.5, color=S.RULE, lw=0.6)
        ax.axvline(0.05, color=S.RED, ls="--", lw=0.8)
        labels = [short(p["readout"]) for p in pairs]
        ax.set(xscale="log", xlim=(0.003, 1.5), ylim=(12.6, -0.6), xticks=[0.01, 0.1, 1], xlabel="IUT-BY q")
        ax.set_xticklabels(["0.01", "0.1", "1"])
        ax.set_yticks(range(13)); ax.set_yticklabels(labels)
        ax.minorticks_off(); clean(ax)
        group_labels(ax, [("Mixed", 0, 7), ("PBMC", 8, 12)], labels)
        legend_above(ax, [pt("o", S.GRAY, "Unrestr.", hollow=True), pt("o", S.GRAY, "Restr.")])
        return {"trajectories": pairs}
    return draw


def _lines(ax, d):
    for t, color, name in (("brain", S.BLUE, "Mixed"), ("pbmc", S.TEAL, "PBMC")):
        pts = d.ladder["tissues"][t]["points"]
        x = np.array([r["alpha"] for r in pts]); y = np.array([r["mean_rho"] for r in pts]); sd = np.array([r["sd_rho"] for r in pts])
        ax.plot(x, y, color=color, lw=1.0, label=name)
        ax.fill_between(x, y - sd, y + sd, color=color, alpha=0.14, lw=0)


def panel_e(d):
    def draw(ax):
        _lines(ax, d)
        ax.set(xlim=(0, 1), ylim=(-0.02, 1.02), xticks=[0, 0.5, 1], yticks=[0, 0.5, 1], xlabel="Injected residual fraction",
               ylabel="Recovered ρ")
        ax.legend(loc="upper left", handlelength=1.2, labelspacing=0.25, borderpad=0.2)
        clean(ax, grid="both")
        inset = ax.inset_axes([0.60, 0.14, 0.36, 0.36])
        _lines(inset, d)
        inset.set(xlim=(0, 0.022), ylim=(-0.004, 0.026), xticks=[0, 0.02], yticks=[0, 0.02])
        inset.set_title("Low dose", fontsize=S.FS_NOTE, pad=1.5)
        inset.tick_params(labelsize=S.FS_TICK, length=1.8, pad=1.2)
        clean(inset, grid="both")
        pts = {t: d.ladder["tissues"][t]["points"] for t in ("brain", "pbmc")}
        return {"tissues": ["Mixed", "PBMC"], "alpha": [r["alpha"] for r in pts["brain"]],
                "mean_rho": {t: [r["mean_rho"] for r in v] for t, v in pts.items()},
                "sd_rho": {t: [r["sd_rho"] for r in v] for t, v in pts.items()},
                "graphs_per_dose": sorted({r["n_replicates"] for v in pts.values() for r in v})}
    return draw


def power_panel(d, tissue):
    def draw(ax):
        curves = d.power["tissues"][tissue]["curves"]
        rec = {}
        for key, color, name in (("a05", S.BLUE, "α = 0.05"), ("a05m", S.ORANGE, "α = 0.05 / m")):
            x = [r["alpha_injected"] for r in curves]
            y = [r[key]["proportion"] for r in curves]
            lo = [r[key]["simultaneous95_per_curve"][0] for r in curves]
            hi = [r[key]["simultaneous95_per_curve"][1] for r in curves]
            ax.plot(x, y, "o-", color=color, ms=2.6, lw=0.9, label=name)
            ax.fill_between(x, lo, hi, color=color, alpha=0.14, lw=0)
            rec[key] = {"alpha": x, "proportion": y, "simultaneous95": [lo, hi], "noise_draws_per_dose": curves[0]["n"]}
        ax.axhline(0.8, color=S.GRAY, ls="--", lw=0.7)
        ax.set(xlim=(0, 0.0155), ylim=(0, 1.04), xticks=[0, 0.005, 0.01, 0.015], yticks=[0, 0.5, 1],
               xlabel="Injected residual fraction", ylabel="Joint rejection")
        ax.set_xticklabels(["0", "0.005", "0.010", "0.015"])
        if tissue == "brain":
            ax.legend(loc="lower right", handlelength=1.2, labelspacing=0.25, borderpad=0.2)
        clean(ax, grid="both")
        return rec
    return draw


def panel_i(d):
    def draw(ax):
        n = np.array([20, 100])
        lower = beta.ppf(0.0025, n, 1)
        ax.errorbar([0, 1], [1, 1], yerr=[1 - lower, [0, 0]], fmt="o", color=S.TEAL, capsize=3, ms=4)
        for x, v in zip([0, 1], lower):
            ax.text(x + 0.12, v + 0.012, fmt(v, 3), fontsize=S.FS_NOTE, va="bottom")
        ax.set(xlim=(-0.6, 1.9), ylim=(0.65, 1.03), xticks=[0, 1], xticklabels=["20", "100"], yticks=[0.7, 0.8, 0.9, 1],
               xlabel="Successes / tests", ylabel="Lower bound")
        clean(ax, grid="y")
        return {"successes_of": n.tolist(), "lower_bound_99_5_clopper_pearson": lower.tolist()}
    return draw


def build(out_height=104.0):
    d = Data()
    c = Canvas("Figure5", out_height)
    r1, r2 = split(c.rect, "v", [1.45, 1.0])
    ab_cell, c_cell, d_cell = split(r1, "h", [1.3, 0.55, 1.0])
    a_cell, b_cell = split(ab_cell, "v", [1.0, 1.0])
    e_cell, g_cell, h_cell, i_cell = split(r2, "h", [1.15, 1.0, 1.0, 0.72])
    c.panel(a_cell, "A", "Mixed-source stress", stress_panel(d, "brain", True), align_left="ab")
    c.panel(b_cell, "B", "PBMC stress", stress_panel(d, "pbmc", True), align_left="ab")
    c.panel(c_cell, "C", "Opportunity mean", panel_c(d))
    c.panel(d_cell, "D", "Restricted IUT-BY", panel_d(d))
    c.panel(e_cell, "E", "Synthetic ladder", panel_e(d), align_y="r2")
    c.panel(g_cell, "F", "Mixed joint rejection", power_panel(d, "brain"), align_y="r2")
    c.panel(h_cell, "G", "PBMC joint rejection", power_panel(d, "pbmc"), align_y="r2")
    c.panel(i_cell, "H", "Precision", panel_i(d), align_y="r2")
    return c, {}
