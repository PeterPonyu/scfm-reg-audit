"""Figure 1: source correction, panel construction, coverage and admission."""
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from . import style as S
from .common import (HIST, JOBS, READOUTS, ROOT, bar_labels, clean, boolean_heat, hbars, load, payload,
                     short)
from .layout import Canvas, fixed, split
from .schematics import workflow, proximal

CONTEXT = {"brain_selected_panel-brain": "New / brain", "brain_selected_panel-pbmc": "New / PBMC",
           "historical_panel_corrected_source-brain": "Old / brain"}


class Data:
    def __init__(self):
        self.source = payload("figure1_source")
        self.design = payload("figure2_design")
        self.coverage = payload("figure3_coverage")
        self.cascade = load(HIST / "analysis/inputs/revision_05613/panel_cascade.json")
        self.contexts = load(ROOT / "inputs/corrective-contexts.json")["jobs"]
        self.fixed_audit = load(HIST / "analysis/inputs/v2/fixed_panel_audit_v2.json")
        self.gate_rows = load(HIST / "figures/data/gate_external_mapping.json")["rows"]
        undefined = {j: len(v["undefined_coexpression_genes"]) for j, v in self.design["prep"]["jobs"].items()}
        assert undefined == {j: self.contexts[j]["undefined_coexpression_genes"] for j in JOBS}
        self.undefined = undefined


def panel_a(d):
    def draw(ax):
        layers = d.source["audit_layers"]
        names = ["Local", "AD1", "Normal1"]
        frac = [r["matched_cells"] / r["audited_cells"] for r in layers]
        ax.barh(range(3), frac, height=0.62, color=S.BLUE)
        for i, r in enumerate(layers):
            ax.text(0.985, i, f"{r['matched_cells']:,}/{r['audited_cells']:,}", ha="right", va="center",
                    color="white", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, 1.0), ylim=(2.55, -0.55), xticks=[0, 0.5, 1], xlabel="Matched / checked cells")
        ax.set_yticks(range(3)); ax.set_yticklabels(names)
        clean(ax)
        return {"labels": names, "matched": [r["matched_cells"] for r in layers],
                "checked": [r["audited_cells"] for r in layers], "features": layers[0]["features"]}
    return draw


def panel_b(d):
    def draw(ax):
        sel = d.source["prep"]["selection"]["brain"]
        types = sel["types"]
        quota = [sel["per_type_quota"][t] for t in types]
        ax.barh(range(len(types)), quota, height=0.66, color=S.CELLTYPES[:len(types)])
        bar_labels(ax, range(len(types)), quota, [str(q) for q in quota], 0.8)
        total = sum(quota)
        ax.set(xlim=(0, 25), ylim=(len(types) - 0.45, -0.55), xticks=[0, 9, 18], xlabel="Nuclei per sample")
        ax.text(0.98, 0.02, f"Σ = {total}", transform=ax.transAxes, ha="right", va="bottom", fontsize=S.FS_NOTE)
        ax.set_yticks(range(len(types))); ax.set_yticklabels(types)
        clean(ax)
        return {"types": types, "quota": quota, "per_sample": total, "samples": 18, "total": total * 18}
    return draw


def panel_c(d):
    def draw(ax):
        tf = d.design["tf_counts"]
        shared, unique, total = d.design["shared_genes"], d.design["panel_specific_genes"], d.design["frozen_genes"]
        assert shared + unique == total
        for i, n in enumerate(tf):
            ax.barh(i, shared, color=S.BLUE, height=0.6)
            ax.barh(i, unique, left=shared, color=S.LIGHT, edgecolor=S.GRAY, lw=0.5, height=0.6)
            ax.text(shared / 2, i, str(shared), ha="center", va="center", color="white", fontsize=S.FS_NOTE)
            ax.text(shared + unique / 2, i, str(unique), ha="center", va="center", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, total), ylim=(1.55, -0.55), xticks=[0, 600, 1200], xlabel="Genes: shared | unique")
        ax.set_yticks([0, 1]); ax.set_yticklabels([f"New · {tf[0]} TF", f"Old · {tf[1]} TF"])
        clean(ax)
        return {"shared": shared, "unique": unique, "tfs": {"new": tf[0], "old": tf[1]}}
    return draw


def panel_d(d):
    def draw(ax):
        steps = d.cascade["steps"]
        labels = ["GENCODE", "Geneformer", "scGPT", "ATAC locus", "Detection", "Cap"]
        genes = [s["genes"] for s in steps]
        ax.barh(range(len(steps)), genes, height=0.66, color=S.BLUE)
        for i, s in enumerate(steps):
            ax.text(genes[i] * 1.08, i, f"{s['genes']:,} · {s['tfs']} TF", va="center", fontsize=S.FS_NOTE)
        ax.set(xscale="log", xlim=(900, 4.2e5), ylim=(len(steps) - 0.45, -0.55), xticks=[1e3, 1e4, 1e5],
               xlabel="Genes retained")
        ax.set_xticklabels(["1,000", "10,000", "100,000"])
        ax.set_yticks(range(len(steps))); ax.set_yticklabels(labels)
        ax.minorticks_off(); clean(ax)
        return {"steps": steps, "labels": labels}
    return draw


def panel_e(d, labels=True):
    def draw(ax):
        peaks = [d.contexts[j]["linked_peaks"] for j in JOBS]
        hbars(ax, [CONTEXT[j] for j in JOBS], peaks, [S.TEAL, S.BLUE, S.GRAY], xlim=(0, 50000),
              xticks=[0, 15000, 30000], xticklabels=["0", "15k", "30k"], xlabel="Linked peaks", height=0.6,
              annotate=[f"{p:,}" for p in peaks], annotate_dx=700, show_labels=labels)
        return {"contexts": [CONTEXT[j] for j in JOBS], "linked_peaks": peaks,
                "tfs": [d.contexts[j]["tfs"] for j in JOBS]}
    return draw


def panel_f(d, labels=True):
    def draw(ax):
        values = [d.undefined[j] for j in JOBS]
        hbars(ax, [CONTEXT[j] for j in JOBS], values, [S.TEAL, S.ORANGE, S.ORANGE], xlim=(0, 58),
              xticks=[0, 20, 40], xlabel="Zero-variance genes", height=0.6, annotate=[str(v) for v in values],
              annotate_dx=1.0, show_labels=labels)
        return {"contexts": [CONTEXT[j] for j in JOBS], "undefined_comparator_genes": values}
    return draw


def panel_g(d, legend=False):
    rows = d.coverage["inferred"]["rows"]

    def draw(ax):
        states = np.zeros((6, 3), int)
        gene = np.zeros((6, 3))
        pair = np.zeros((6, 3))
        for i, m in enumerate(READOUTS):
            for j, job in enumerate(JOBS):
                r = next(v for v in rows if v["job"] == job and v["readout"] == m)
                gene[i, j], pair[i, j] = r["gene_coverage_fraction"], r["evaluated_pair_coverage_fraction"]
                states[i, j] = 2 if r["primary_statistics_eligible"] else (1 if j == 0 else 0)
        cmap = ListedColormap([S.LIGHT, S.AMBER, S.TEAL])
        ax.imshow(states, cmap=cmap, vmin=0, vmax=2, aspect="auto")
        for i in range(6):
            for j in range(3):
                ax.text(j, i, f"{100 * gene[i, j]:.1f}/{100 * pair[i, j]:.1f}", ha="center", va="center",
                        fontsize=S.FS_NOTE, color="white" if states[i, j] == 2 else "black")
        for k in range(1, 6):
            ax.axhline(k - 0.5, color="white", lw=0.7)
        for k in range(1, 3):
            ax.axvline(k - 0.5, color="white", lw=0.7)
        ax.set_yticks(range(6)); ax.set_yticklabels([short(m) for m in READOUTS])
        ax.set_xticks(range(3)); ax.set_xticklabels([CONTEXT[j] for j in JOBS])
        ax.tick_params(length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        if legend:
            ax.legend(handles=[Patch(facecolor=S.TEAL, label="Admitted"), Patch(facecolor=S.AMBER, label="Model gap"),
                               Patch(facecolor=S.LIGHT, edgecolor=S.GRAY, lw=0.5, label="No comparator")],
                      loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, handlelength=1.0,
                      handletextpad=0.4, columnspacing=1.0, borderaxespad=0.0, borderpad=0.0, fontsize=S.FS_NOTE)
        return {"rows": [short(m) for m in READOUTS], "columns": [CONTEXT[j] for j in JOBS],
                "gene_coverage_pct": (100 * gene).tolist(), "pair_coverage_pct": (100 * pair).tolist(),
                "state": states.tolist()}
    return draw


def panel_h(d, labels=True):
    cov = d.coverage["results"]["rows"]

    def draw(ax):
        counts = [sum(r["measured"] and r["job"] == job for r in cov) for job in JOBS]
        totals = [sum(r["job"] == job for r in cov) for job in JOBS]
        assert totals == [24, 24, 24] and sum(counts) == 16
        ax.barh(range(3), counts, height=0.6, color=S.TEAL)
        ax.barh(range(3), np.array(totals) - counts, left=counts, color=S.LIGHT, edgecolor=S.GRAY, lw=0.5, height=0.6)
        for i, v in enumerate(counts):
            ax.text(25, i, f"{v}/24", va="center", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, 31), ylim=(2.55, -0.55), xticks=[0, 12, 24], xlabel="Measured slots")
        ax.set_yticks(range(3)); ax.set_yticklabels([CONTEXT[j] for j in JOBS] if labels else [])
        if not labels:
            ax.tick_params(axis="y", length=0)
        clean(ax)
        return {"contexts": [CONTEXT[j] for j in JOBS], "measured": counts, "registered": totals}
    return draw


def panel_i(d):
    def draw(ax):
        models = list(dict.fromkeys(r["readout"] for r in d.gate_rows))
        scopes = [("brain", "pooled"), ("pbmc", "pooled"), ("brain", "per_type"), ("pbmc", "per_type")]
        table = d.fixed_audit["model_coverage_table"]
        mat = np.array([[m in table[f"{t}_{s}_FMs"] for t, s in scopes] for m in models], int)
        boolean_heat(ax, mat, colors=(S.LIGHT, S.BLUE), row_labels=[short(m) for m in models],
                     col_labels=["M\npool", "P\npool", "M\ntype", "P\ntype"], col_rotation=0)
        return {"rows": [short(m) for m in models], "columns": ["M pooled", "P pooled", "M per type", "P per type"],
                "available": mat.tolist()}
    return draw


def panel_j(d):
    def draw(ax):
        panel_tfs = 446
        total = panel_tfs * 1200
        pbmc = d.fixed_audit["pooled"]["pbmc"]["n_pairs"]
        brain = d.fixed_audit["pooled"]["brain"]["n_pairs"]
        assert pbmc == panel_tfs * 1199
        values = [total, pbmc, brain]
        ax.barh(range(3), values, height=0.6, color=[S.GRAY, S.TEAL, S.BLUE])
        for i, v in enumerate(values):
            ax.text(v * 0.5, i, f"{v:,}", ha="center", va="center", color="white", fontsize=S.FS_NOTE)
        ax.set(xlim=(0, 600000), ylim=(2.55, -0.55), xticks=[0, 300000, 600000], xlabel="TF–gene pairs")
        ax.set_xticklabels(["0", "300k", "600k"])
        ax.set_yticks(range(3)); ax.set_yticklabels(["All", "No self", "M mask"])
        clean(ax)
        return {"labels": ["All pairs", "Non-self pairs", "Mixed-source mask"], "pairs": values}
    return draw


def build(out_height=158.0):
    d = Data()
    c = Canvas("Figure1", out_height)
    schematic_row, results = split(c.rect, "v", [fixed(48.0), 1.0])
    workflow_cell, proximal_cell = split(schematic_row, "h", [1.0, 1.0])
    c.panel(workflow_cell, "A", "Source-to-audit workflow", workflow)
    c.panel(proximal_cell, "B", "Proximal proxy construction", proximal)
    top, mid, bot = split(results, "v", [4.1, 2.55, 3.9])
    a_c, b_cell, d_cell = split(top, "h", [1.0, 0.9, 1.45])
    a_cell, c_cell = split(a_c, "v", [1.05, 1.0])
    e_cell, f_cell, g_cell = split(mid, "h", [1.35, 0.85, 1.0])
    i_cell, j_cell, h_cell = split(bot, "h", [0.92, 1.6, 0.78])
    c.panel(a_cell, "C", "Source fingerprints", panel_a(d))
    c.panel(c_cell, "D", "Panel overlap", panel_c(d))
    c.panel(b_cell, "E", "Selection quota", panel_b(d))
    c.panel(d_cell, "F", "Historical panel cascade", panel_d(d))
    c.panel(e_cell, "G", "Linked peaks", panel_e(d), align_y="r2")
    c.panel(f_cell, "H", "Undefined genes", panel_f(d, labels=False), align_y="r2")
    c.panel(g_cell, "I", "Measured slots", panel_h(d, labels=False), align_y="r2")
    c.panel(i_cell, "J", "Graph-cache availability", panel_i(d), align_y="r3")
    c.panel(j_cell, "K", "Gene / pair coverage (%)", panel_g(d, legend=True), align_y="r3")
    c.panel(h_cell, "L", "Edge accounting", panel_j(d), align_y="r3")
    return c, {}
