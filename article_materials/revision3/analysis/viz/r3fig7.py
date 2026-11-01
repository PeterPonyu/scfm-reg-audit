"""Figure 7 (R3): five bone-marrow evaluation donors on the common finite gene set."""
import json

import numpy as np

from . import style as S
from .common import ROOT, clean, legend_above, pt, short
from .layout import Canvas, split
from . import fig6 as old

RES = ROOT / "evidence" / "R3_RESULTS.json"
FIVE = ["geneformer_embed", "random_init_geneformer", "scgpt_encoder", "scfoundation_encoder", "uce_encoder"]
DEGC = {"reference": (S.BLUE, "o"), "own": (S.TEAL, "s"), "none": (S.ORANGE, "^")}
DONOR_COL = [S.BLUE, S.TEAL, S.ORANGE, S.PURPLE, S.RED]


class Data:
    def __init__(self):
        b = json.loads(RES.read_text())["bmmc"]
        assert b["streams_found"] == 30
        self.b = b
        self.rows = b["rows"]
        self.summary = {s["readout"]: s for s in b["donor_summary"]}
        self.donors = sorted({r["donor"] for r in self.rows})
        self.old = old.Data()


def coverage_panel(d):
    def draw(ax):
        ex = d.b["common"]["excluded_by_donor"]
        donors = d.old.donors
        vals = [len(ex[x]) for x in donors]
        ax.bar(range(len(donors)), vals, color=S.ORANGE, width=0.6)
        for k, v in enumerate(vals):
            ax.text(k, v + 0.4, str(v), ha="center", va="bottom", fontsize=S.FS_NOTE)
        ax.set(ylim=(0, 18), yticks=[0, 5, 10, 15], ylabel="Genes excluded")
        ax.set_xticks(range(len(donors))); ax.set_xticklabels(donors, rotation=0)
        ax.set_xlabel(f"Common set {d.b['common']['genes_common']:,} genes, {d.b['common']['tfs_common']} TFs")
        clean(ax, grid="y")
        return {"excluded_by_donor": {x: len(ex[x]) for x in donors}, "common_genes": d.b["common"]["genes_common"]}
    return draw


def donor_panel(d, key, xlabel):
    def draw(ax):
        rec = []
        for k, m in enumerate(FIVE):
            s = d.summary[m]
            vals = s["observed"] if key == "observed" else s["differences"]
            inter = s["donor_interval"] if key == "observed" else s["paired_baseline_interval"]
            for j, (dn, v) in enumerate(zip(s["donors"], vals)):
                ax.scatter(v, k - 0.25 + 0.1 * j, s=9, color=DONOR_COL[j], lw=0, zorder=3)
            lo, hi = inter["ci95"]
            pnt = inter["point"]
            ax.errorbar(pnt, k + 0.32, xerr=[[pnt - lo], [hi - pnt]], fmt="D", color="black", ms=3.0, capsize=1.6, lw=0.8)
            rec.append({"readout": m, "values": vals, "mean": pnt, "ci95": inter["ci95"],
                        "signflip_p": s.get("paired_baseline_signflip_p") if key != "observed" else None})
            if key != "observed":
                ax.text(1.02, k, f"p = {s['paired_baseline_signflip_p']:.3g}", transform=ax.get_yaxis_transform(),
                        ha="left", va="center", fontsize=S.FS_NOTE)
        ax.axvline(0, color=S.GRAY, lw=0.6)
        if key == "observed":
            ax.set(xlim=(-0.016, 0.038), xticks=[-0.01, 0, 0.01, 0.02, 0.03])
            ax.set_xticklabels(["−0.01", "0", "0.01", "0.02", "0.03"])
            legend_above(ax, [pt("o", c, dn) for c, dn in zip(DONOR_COL, d.summary[FIVE[0]]["donors"])], ncol=3)
        else:
            ax.set(xlim=(-0.035, 0.026), xticks=[-0.03, -0.02, -0.01, 0, 0.01, 0.02])
            ax.set_xticklabels(["−0.03", "−0.02", "−0.01", "0", "0.01", "0.02"])
        ax.set(ylim=(4.7, -0.6), xlabel=xlabel)
        ax.set_yticks(range(5)); ax.set_yticklabels([short(m) for m in FIVE])
        clean(ax)
        return {"rows": rec}
    return draw


def pvalue_panel(d):
    def draw(ax):
        rec = []
        ypos = {m: i for i, m in enumerate(FIVE)}
        offs = {"reference": -0.24, "own": 0.0, "none": 0.24}
        for r in d.rows:
            col, mk = DEGC[r["degree_mode"]]
            ax.scatter(-np.log10(r["pIUT"]), ypos[r["readout"]] + offs[r["degree_mode"]], marker=mk, s=10, lw=0.7, zorder=3,
                       facecolor=col if r["null_mode"] == "restricted" else "white", edgecolor=col)
            rec.append({k: r[k] for k in ("donor", "readout", "degree_mode", "null_mode", "observed", "pIUT", "q_BY", "p_Holm")})
        ax.axvline(-np.log10(0.05 / 150), color=S.RED, ls="--", lw=0.7)
        ax.axvline(-np.log10(0.05), color=S.GRAY, ls=":", lw=0.7)
        ax.set(xlim=(0, 4.3), xticks=[0, 1, 2, 3, 4], ylim=(4.6, -0.6), xlabel="−log10 joint p (150 hypotheses)")
        ax.set_yticks(range(5)); ax.set_yticklabels([short(m) for m in FIVE])
        clean(ax)
        legend_above(ax, [pt("o", S.BLUE, "Reference"), pt("s", S.TEAL, "Own"), pt("^", S.ORANGE, "None")])
        return {"rows": rec}
    return draw


def build(out_height=104.0):
    d = Data()
    c = Canvas("Figure7", out_height)
    roles, top, bot = split(c.rect, "v", [0.32, 1.0, 1.0])
    b_cell, e_cell = split(top, "h", [0.8, 1.2])
    c_cell, d_cell = split(bot, "h", [1.0, 1.0])
    c.panel(roles, "A", "Frozen donor roles", old.panel_a(d.old))
    c.panel(b_cell, "B", "Common finite gene set", coverage_panel(d), align_y="t")
    c.panel(e_cell, "C", "Joint p by donor and design", pvalue_panel(d), align_y="t")
    c.panel(c_cell, "D", "Donor effects (reference degree)", donor_panel(d, "observed", "Partial ρ; diamond, donor mean and 95% CI"), align_y="b")
    c.panel(d_cell, "E", "Model minus co-expression", donor_panel(d, "difference", "Δρ; diamond, donor mean and 95% CI"), align_y="b")
    return c, {"family_150": d.b.get("family_150"), "donor_summary": d.b.get("donor_summary")}
