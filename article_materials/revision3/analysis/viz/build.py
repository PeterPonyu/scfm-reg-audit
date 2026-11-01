"""Build figures: ``python -m viz.build fig1 fig2 ...`` (add --preview to skip the large exports)."""
import importlib
import json
import sys
import tempfile
from pathlib import Path

from . import qa
from .common import ROOT
from .export import export

FIGURES = {"fig1": "Figure1", "fig2": "Figure2", "fig3": "Figure3", "fig4": "Figure4", "fig5": "Figure5",
           "fig6": "Figure6", "fig7": "Figure7"}
OUT = ROOT / "figures"
QA = ROOT / "qa"
PREVIEW = Path(tempfile.gettempdir()) / "screg-preview"


def run(key: str, preview: bool):
    module = importlib.import_module(f"viz.{key}")
    canvas, data = module.build()
    canvas.fit()
    report = qa.full_report(canvas)
    stem = FIGURES[key]
    if preview:
        PREVIEW.mkdir(exist_ok=True)
        canvas.fig.savefig(PREVIEW / f"{stem}.png", dpi=220)
    else:
        export(canvas, OUT, stem, data, report)
        QA.mkdir(exist_ok=True)
        (QA / f"layout-{stem}.json").write_text(json.dumps(report, indent=1) + "\n")
    brief = {k: report[k] for k in ("figure", "size_mm", "passed", "minimum_font_pt", "white_fraction",
                                    "thickest_blank_row_band_mm", "thickest_blank_column_corridor_mm",
                                    "largest_empty_rectangle_fraction", "largest_empty_rectangle_mm", "largest_empty_rectangle_panel")}
    brief.update({"overlaps": report["overlaps"][:6], "outside_canvas": report["outside_canvas"][:6],
                  "outside_cell": report["outside_cell"][:6], "axes_overlaps": report["inner_axes_overlaps"][:6],
                  "min_panel_ink": report["minimum_panel_ink"]})
    print(json.dumps(brief))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    for k in (args or list(FIGURES)):
        run(k, "--preview" in sys.argv)
