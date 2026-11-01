"""Build R3 figures: ``python -m viz.build_r3 r3fig2 ... [--preview]``."""
import importlib
import json
import sys
from pathlib import Path

from . import qa
from .common import ROOT
from .export import export

FIGURES = {"r3fig1": "Figure1", "r3fig2": "Figure2", "r3fig3": "Figure3", "r3fig4": "Figure4", "r3fig5": "Figure5",
           "r3fig6": "Figure6", "r3fig7": "Figure7", "r3figA4": "FigureA4"}
OUT = ROOT / "figures"
QA = ROOT / "qa"
PREVIEW = ROOT / "qa" / "preview"


def run(key, preview):
    module = importlib.import_module(f"viz.{key}")
    canvas, data = module.build()
    canvas.fit()
    report = qa.full_report(canvas)
    stem = FIGURES[key]
    if preview:
        PREVIEW.mkdir(parents=True, exist_ok=True)
        canvas.fig.savefig(PREVIEW / f"{stem}.png", dpi=200)
    else:
        export(canvas, OUT, stem, data, report)
        QA.mkdir(exist_ok=True)
        (QA / f"layout-{stem}.json").write_text(json.dumps(report, indent=1) + "\n")
    brief = {k: report[k] for k in ("figure", "size_mm", "passed", "minimum_font_pt", "white_fraction")}
    brief.update({"overlaps": report["overlaps"][:4], "outside_canvas": report["outside_canvas"][:4],
                  "outside_cell": report["outside_cell"][:4], "axes_overlaps": report["inner_axes_overlaps"][:4]})
    print(json.dumps(brief))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    for k in (args or list(FIGURES)):
        run(k, "--preview" in sys.argv)
