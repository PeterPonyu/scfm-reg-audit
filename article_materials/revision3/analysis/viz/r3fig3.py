"""Figure 3 (R3): source-matched alignment on identical edge sets (common support, 9,999 draws)."""
import json

import numpy as np

from . import style as S
from .common import ROOT, clean, corrected_rows, legend_above, pt, short
from .layout import Canvas, split

RES = ROOT / "evidence" / "R3_RESULTS.json"
JOBS = ["brain_selected_panel-brain", "brain_selected_panel-pbmc", "historical_panel_corrected_source-brain"]
CTX = {"brain_selected_panel-brain": "New / brain", "brain_selected_panel-pbmc": "New / PBMC",
       "historical_panel_corrected_source-brain": "Old / brain"}
ORDER = ["geneformer_embed", "geneformer_attention", "random_init_geneformer", "scgpt_encoder", "scfoundation_encoder", "uce_encoder"]
DEG = {"none": (S.ORANGE, "^"), "own": (S.TEAL, "s")}


class Data:
    def __init__(self):
        r = json.loads(RES.read_text())
        cs = r["common_support"]
        self.six = cs["six"]["rows"]
        self.gl = cs["gene_level"]["rows"]
        self.fam = cs["six"]["family"]
        assert len(self.six) == 72 and len(self.gl) == 60
        self.registered = corrected_rows()

    def pick(self, rows, job, m, deg, null):
        return next(x for x in rows if (x["job"], x["readout"], x["degree_mode"], x["null_mode"]) == (job, m, deg, null))


def effect_panel(d):
    def draw(ax):
        rec, ticks, names, y = [], [], [], 0.0
        for job in JOBS:
            ax.axhspan(y - 0.6, y + 5.6, color=S.LIGHT if JOBS.index(job) % 2 else "white", lw=0, zorder=0)
            for k, m in enumerate(ORDER):
                for deg, off in (("none", -0.16), ("own", 0.16)):
                    r = d.pick(d.six, job, m, deg, "restricted")
                    col, mk = DEG[deg]
                    ax.scatter(r["observed"], y + k + off, marker=mk, color=col, s=13, zorder=3, lw=0)
                    rec.append({"context": CTX[job], "readout": m, "degree": deg, "rho": r["observed"]})
                ticks.append(y + k); names.append(short(m))
            ax.text(0.047, y + 2.5, CTX[job], rotation=270, ha="center", va="center", fontsize=S.FS_NOTE)
            y += 6.6
        ax.axvline(0, color=S.GRAY, lw=0.6)
        ax.set(xlim=(-0.03, 0.052), xticks=[-0.02, 0, 0.02, 0.04], ylim=(y - 1.0, -0.8), xlabel="Incremental partial ρ")
        ax.set_xticklabels(["−0.02", "0", "0.02", "0.04"])
        ax.set_yticks(ticks); ax.set_yticklabels(names)
        clean(ax)
        legend_above(ax, [pt("^", S.ORANGE, "No degree"), pt("s", S.TEAL, "Own degree")])
        return {"rows": rec}
    return draw


def pvalue_panel(d):
    def draw(ax):
        rec = []
        y = {m: i for i, m in enumerate(ORDER)}
        for r in d.six:
            col, mk = DEG[r["degree_mode"]]
            jit = {"brain_selected_panel-brain": -0.24, "brain_selected_panel-pbmc": 0.0,
                   "historical_panel_corrected_source-brain": 0.24}[r["job"]]
            filled = r["null_mode"] == "restricted"
            ax.scatter(-np.log10(r["pIUT"]), y[r["readout"]] + jit, marker=mk, s=12, lw=0.7, zorder=3,
                       facecolor=col if filled else "white", edgecolor=col)
            rec.append({k: r[k] for k in ("job", "readout", "degree_mode", "null_mode", "pIUT", "q_BY", "p_Holm")})
        m = len(d.six)
        holm1 = -np.log10(0.05 / m)
        by1 = -np.log10(0.05 / (m * np.sum(1 / np.arange(1, m + 1))))
        ax.axvline(-np.log10(0.05), color=S.GRAY, ls=":", lw=0.7)
        ax.axvline(holm1, color=S.RED, ls="--", lw=0.7)
        ax.axvline(by1, color=S.PURPLE, ls="-.", lw=0.7)
        ax.text(holm1 - 0.06, -0.75, "Holm", fontsize=S.FS_NOTE, color=S.RED, va="center", ha="right")
        ax.text(by1 + 0.06, 5.35, "BY", fontsize=S.FS_NOTE, color=S.PURPLE, va="center", ha="left")
        ax.set(xlim=(0, 4.3), xticks=[0, 1, 2, 3, 4], ylim=(5.6, -1.1), xlabel="−log10 joint p (9,999 draws)")
        ax.set_yticks(range(6)); ax.set_yticklabels([short(m) for m in ORDER])
        clean(ax)
        return {"rows": rec, "holm_first_step": 0.05 / m, "family_slots": m,
                "min_BY": d.fam["min_BY"], "min_Holm": d.fam["min_Holm"]}
    return draw


def component_panel(d):
    def draw(ax):
        rec = []
        for r in d.six:
            col = S.TEAL if r["null_mode"] == "restricted" else S.ORANGE
            ax.scatter(-np.log10(r["pM"]), -np.log10(r["pD"]), s=10, color=col, lw=0, alpha=0.9, zorder=3)
            rec.append({"pM": r["pM"], "pD": r["pD"], "null": r["null_mode"]})
        ax.plot([0, 4.2], [0, 4.2], color=S.GRAY, lw=0.6, ls=":")
        ax.set(xlim=(0, 4.2), ylim=(0, 4.2), xticks=[0, 2, 4], yticks=[0, 2, 4],
               xlabel="−log10 p, gene label M", ylabel="−log10 p, target shuffle D")
        clean(ax, grid="both")
        ax.legend(handles=[pt("o", S.TEAL, "Restricted"), pt("o", S.ORANGE, "Unrestricted")], loc="lower right",
                  frameon=False, fontsize=S.FS_NOTE, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.3)
        return {"rows": rec}
    return draw


def registered_panel(d):
    def draw(ax):
        rec = []
        for r in d.registered:
            c = d.pick(d.six, "brain_selected_panel-brain", r["readout"], r["degree_mode"], r["null_mode"])
            col, mk = DEG[r["degree_mode"]]
            ax.scatter(r["observed"], c["observed"], marker=mk, color=col, s=14, lw=0, zorder=3)
            rec.append({"readout": r["readout"], "degree": r["degree_mode"], "null": r["null_mode"],
                        "registered": r["observed"], "common_support": c["observed"]})
        ax.plot([-0.012, 0.006], [-0.012, 0.006], color=S.GRAY, lw=0.6, ls=":")
        ax.set(xlim=(-0.012, 0.006), ylim=(-0.012, 0.006), xticks=[-0.01, -0.005, 0, 0.005], yticks=[-0.01, -0.005, 0, 0.005],
               xlabel="Registered ρ (4 readouts)", ylabel="Common-support ρ")
        ax.set_xticklabels(["−0.01", "", "0", ""]); ax.set_yticklabels(["−0.01", "", "0", ""])
        clean(ax, grid="both")
        return {"rows": rec}
    return draw


def sensitivity_panel(d):
    def draw(ax):
        rec = []
        for r in d.gl:
            c = d.pick(d.six, r["job"], r["readout"], r["degree_mode"], r["null_mode"])
            if r["null_mode"] != "restricted":
                continue
            col, mk = DEG[r["degree_mode"]]
            ax.scatter(c["observed"], r["observed"], marker=mk, color=col, s=12, lw=0, zorder=3)
            rec.append({"job": r["job"], "readout": r["readout"], "degree": r["degree_mode"], "six": c["observed"],
                        "gene_level": r["observed"]})
        ax.plot([-0.03, 0.04], [-0.03, 0.04], color=S.GRAY, lw=0.6, ls=":")
        ax.set(xlim=(-0.025, 0.03), ylim=(-0.025, 0.03), xticks=[-0.02, 0, 0.02], yticks=[-0.02, 0, 0.02],
               xlabel="ρ, six-readout edges", ylabel="ρ, all gene-level edges")
        ax.set_xticklabels(["−0.02", "0", "0.02"]); ax.set_yticklabels(["−0.02", "0", "0.02"])
        clean(ax, grid="both")
        return {"rows": rec}
    return draw


def build(out_height=118.0):
    d = Data()
    c = Canvas("Figure3", out_height)
    left, right = split(c.rect, "h", [1.0, 1.55])
    rt, rb = split(right, "v", [1.0, 1.0])
    b_cell, c_cell = split(rt, "h", [1.35, 1.0])
    d_cell, e_cell = split(rb, "h", [1.0, 1.0])
    c.panel(left, "A", "Effects on identical edges", effect_panel(d))
    c.panel(b_cell, "B", "Joint p, 72 hypotheses", pvalue_panel(d), align_y="t")
    c.panel(c_cell, "C", "Component probabilities", component_panel(d), align_y="t")
    c.panel(d_cell, "D", "Registered vs common support", registered_panel(d), align_y="b")
    c.panel(e_cell, "E", "Attention-free sensitivity", sensitivity_panel(d), align_y="b")
    return c, {"family": d.fam}
