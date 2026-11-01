"""Appendix Figure A4 (R3): PBMC descriptive effects, TF-disjoint probe and construct-only transfer.

Identical panels to the Revision 2 Figure 7 except that the mixed-source per-cell-type panel
(aortic RNA cells stratified by brain cell-type labels) is removed.
"""
import inspect

from . import fig7 as old
from .layout import Canvas, split

# Same panel code as Revision 2, with horizontal column labels for the full-width layout.
_src = inspect.getsource(old.pertype_panel).replace("col_rotation=90", "col_rotation=0")
_ns = dict(vars(old))
exec(compile(_src, "fig7_pertype_horizontal", "exec"), _ns)
pertype_panel = _ns["pertype_panel"]


def build(out_height=104.0):
    d = old.Data()
    c = Canvas("FigureA4", out_height)
    r1, r2, r3 = split(c.rect, "v", [0.78, 1.2, 0.92])
    a_cell, _pad = split(r1, "h", [1.0, 0.0001])
    e_cell, f_cell, g_cell, h_cell = split(r2, "h", [1.0, 0.6, 1.0, 0.7])
    c_cell, i_cell, j_cell = split(r3, "h", [1.0, 1.0, 1.0])
    c.panel(a_cell, "A", "PBMC per-type effects", pertype_panel(d, "pbmc"), align_y="r1")
    c.panel(e_cell, "B", "Adjusted probe", old.panel_e(d), align_y="r2")
    c.panel(f_cell, "C", "Probe contrasts", old.panel_f(d), align_y="r2")
    c.panel(g_cell, "D", "Covariate sensitivity", old.panel_g(d), align_y="r2")
    c.panel(h_cell, "E", "Per-TF spread", old.panel_h(d), align_y="r2")
    c.panel(c_cell, "F", "Cell-half degree", old.panel_c(d), align_y="r3")
    c.panel(i_cell, "G", "Construct-only transfer", old.panel_i(d), align_y="r3")
    c.panel(j_cell, "H", "Prior / overlay coverage", old.panel_j(d), align_y="r3")
    return c, {"derived_from": "Revision 2 Figure 7 without its mixed-source per-cell-type panel"}
