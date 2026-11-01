"""Figure 5 (R3): symmetric comparison with matched co-expression and the CTCF occupancy reference."""
import json

from . import style as S
from .common import ROOT, clean, short
from .layout import Canvas, split
from . import fig6 as old

RES = ROOT / "evidence" / "R3_RESULTS.json"
JOBS = ["brain_selected_panel-brain", "brain_selected_panel-pbmc", "historical_panel_corrected_source-brain"]
CTX = {"brain_selected_panel-brain": "New / brain", "brain_selected_panel-pbmc": "New / PBMC",
       "historical_panel_corrected_source-brain": "Old / brain"}
ORDER = ["geneformer_embed", "geneformer_attention", "random_init_geneformer", "scgpt_encoder", "scfoundation_encoder", "uce_encoder"]
DEG = {"none": (S.ORANGE, "^", "No degree"), "own": (S.TEAL, "s", "Own degree")}


class Data:
    def __init__(self):
        r = json.loads(RES.read_text())
        self.boot = {(b["job"], b["degree_mode"]): b for b in r["symmetric_bootstrap"]}
        self.old = old.Data()


def alignment_panel(d, degree, labels):
    def draw(ax):
        rec, ticks, names, y = [], [], [], 0.0
        for job in JOBS:
            b = d.boot[(job, degree)]
            co = b["coexpression_rho"]
            lo, hi = b["coexpression_ci95"]
            ax.axhspan(y - 0.5, y + 6.5, color=S.LIGHT if JOBS.index(job) % 2 else "white", lw=0, zorder=0)
            ax.errorbar(co, y, xerr=[[co - lo], [hi - co]], fmt="D", color=S.GRAY, ms=3.2, capsize=1.4, lw=0.75)
            ticks.append(y); names.append("Co-expr.")
            rows = {r["readout"]: r for r in b["readouts"]}
            for k, m in enumerate(ORDER, start=1):
                r = rows[m]
                v = r["fm_rho"]
                a, z = r["fm_ci95"]
                col, mk, _ = DEG[degree]
                ax.errorbar(v, y + k, xerr=[[v - a], [z - v]], fmt=mk, color=col, ms=3.0, capsize=1.4, lw=0.75)
                ticks.append(y + k); names.append(short(m))
                rec.append({"context": CTX[job], "readout": m, "rho": v, "ci95": r["fm_ci95"]})
            rec.append({"context": CTX[job], "readout": "co-expression", "rho": co, "ci95": [lo, hi]})
            ax.text(0.0915, y + 3.0, CTX[job], rotation=270, ha="center", va="center", fontsize=S.FS_NOTE)
            y += 7.6
        ax.axvline(0, color=S.GRAY, lw=0.6)
        if degree == "none":
            ax.set(xlim=(-0.032, 0.098), xticks=[0, 0.04, 0.08], xlabel="Partial ρ")
            ax.set_xticklabels(["0", "0.04", "0.08"])
        else:
            ax.set(xlim=(-0.014, 0.0185), xticks=[-0.01, 0, 0.01], xlabel="Partial ρ")
            ax.set_xticklabels(["−0.01", "0", "0.01"])
            for t in ax.texts:
                t.set_x(0.0168)
        ax.set_ylim(y - 1.1, -0.8)
        ax.set_yticks(ticks)
        ax.set_yticklabels(names if labels else [])
        if not labels:
            ax.tick_params(axis="y", length=0)
        clean(ax)
        return {"degree": degree, "rows": rec}
    return draw


def difference_panel(d):
    def draw(ax):
        rec, y, ticks, names = [], 0.0, [], []
        for job in JOBS:
            ax.axhspan(y - 0.6, y + 5.6, color=S.LIGHT if JOBS.index(job) % 2 else "white", lw=0, zorder=0)
            for k, m in enumerate(ORDER):
                for degree, off in (("none", -0.17), ("own", 0.17)):
                    r = next(v for v in d.boot[(job, degree)]["readouts"] if v["readout"] == m)
                    v, (a, z) = r["difference"], r["difference_ci95"]
                    col, mk, _ = DEG[degree]
                    ax.errorbar(v, y + k + off, xerr=[[v - a], [z - v]], fmt=mk, color=col, ms=2.8, capsize=1.2, lw=0.7)
                    rec.append({"context": CTX[job], "readout": m, "degree": degree, "difference": v, "ci95": r["difference_ci95"]})
                ticks.append(y + k); names.append(short(m))
            ax.text(0.0505, y + 2.5, CTX[job], rotation=270, ha="center", va="center", fontsize=S.FS_NOTE)
            y += 6.6
        ax.axvline(0, color=S.GRAY, lw=0.6)
        ax.set(xlim=(-0.095, 0.056), ylim=(y - 1.0, -0.8), xticks=[-0.08, -0.04, 0, 0.04],
               xlabel="Model minus co-expression ρ")
        ax.set_xticklabels(["−0.08", "−0.04", "0", "0.04"])
        ax.set_yticks(ticks); ax.set_yticklabels(names)
        clean(ax)
        return {"rows": rec}
    return draw


def build(out_height=152.0):
    d = Data()
    c = Canvas("Figure5", out_height)
    top, ctcf_strip, bottom = split(c.rect, "v", [1.0, 0.075, 0.42])
    a_cell, b_cell, c_cell = split(top, "h", [1.12, 0.82, 1.15])
    f_cell, g_cell = split(bottom, "h", [1.0, 1.0])
    c.panel(a_cell, "A", "Structural covariates only", alignment_panel(d, "none", True), align_y="top")
    c.panel(b_cell, "B", "Plus own degree", alignment_panel(d, "own", False), align_y="top")
    c.panel(c_cell, "C", "Symmetric difference", difference_panel(d), align_y="top")
    c.panel(ctcf_strip, "D", "CTCF targets (new panel, genuine brain)", old.panel_f(d.old))
    c.panel(f_cell, "E", "CTCF ROC-AUC", old.metric_panel(d.old, "auc", 0, (0.38, 0.72), [0.4, 0.5, 0.6], True), align_y="bot")
    c.panel(g_cell, "F", "CTCF average precision", old.metric_panel(d.old, "ap", 0, (0.79, 0.955), [0.8, 0.85, 0.9], True), align_y="bot")
    return c, {"source": "evidence/R3_RESULTS.json symmetric_bootstrap; inputs/ctcf-compact.json"}
