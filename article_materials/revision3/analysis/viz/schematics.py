"""Restore the R1 workflow and proximal-linkage illustrations with current labels.

The brain contour and matrix/network glyphs retain the original drawing grammar;
the source hierarchy and audit labels follow the corrected manuscript. All shapes
are editable vector artists. Coordinates are schematic, never measured data.
"""
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Rectangle

from . import style as S


def setup(ax):
    ax.set(xlim=(0, 100), ylim=(0, 100))
    ax.axis("off")


def text(ax, x, y, value, **kwargs):
    return ax.text(x, y, value, ha="center", va="center", fontsize=S.FS_NOTE,
                   fontweight="normal", linespacing=1.12, **kwargs)


def route(ax, points, color=S.GRAY, arrow=True):
    """Every inter-element connector consists of horizontal/vertical segments."""
    for a, b in zip(points[:-2], points[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=0.7)
    a, b = points[-2:]
    if arrow:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=6,
                                    color=color, lw=0.7, shrinkA=0, shrinkB=0))
    else:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=0.7)


def brain(ax, x, y, size, color=S.BLUE):
    # R1 build_overview.py brain contour, translated from its top-down coordinates.
    commands = [("M",48,17),("C",39,4,20,10,20,24),("C",5,24,2,43,12,52),
                ("C",0,66,11,82,25,82),("C",24,99,43,101,49,88),("L",49,23),
                ("C",53,4,77,8,80,25),("C",95,23,100,44,89,53),
                ("C",102,67,87,86,75,81),("C",74,100,54,99,49,87)]
    verts, codes = [], []
    for cmd, *coords in commands:
        for i in range(0, len(coords), 2):
            verts.append((x + coords[i] * size / 100, y + (100 - coords[i+1]) * size / 100))
            codes.append({"M":Path.MOVETO,"L":Path.LINETO,"C":Path.CURVE4}[cmd])
    ax.add_patch(PathPatch(Path(verts,codes), fc="#E8EFF6", ec=color, lw=0.75))


def cells(ax, x, y, size, color=S.TEAL):
    for a,b,r in ((.25,.65,.20),(.72,.72,.19),(.53,.24,.22)):
        ax.add_patch(Circle((x+a*size,y+b*size),r*size,fc="#E9F4F1",ec=color,lw=.65))
        ax.add_patch(Circle((x+a*size,y+b*size),r*size*.48,fc=color,ec="none",alpha=.5))


def matrix(ax, x, y, w, h, color=S.BLUE):
    ax.add_patch(Rectangle((x,y),w,h,fc="white",ec=color,lw=.65))
    for j in range(3):
        for k in range(5):
            ax.add_patch(Rectangle((x+.5+k*w/5,y+.5+j*h/3),w/5-.9,h/3-.9,
                                   fc=color if (j+k)%3==0 else "#DCE7EC",ec="none"))


def graph(ax, x, y, w, h, color):
    nodes=[(.12,.18),(.12,.82),(.48,.50),(.86,.16),(.86,.84)]
    for points in (((.12,.18),(.48,.18),(.48,.50)),((.12,.82),(.48,.82),(.48,.50)),
                   ((.48,.50),(.86,.50),(.86,.16)),((.48,.50),(.86,.50),(.86,.84))):
        route(ax,[(x+a*w,y+b*h) for a,b in points],color,arrow=False)
    for a,b in nodes:
        ax.add_patch(Circle((x+a*w,y+b*h),.8,fc="white",ec=color,lw=.75))


def workflow(ax):
    setup(ax)
    brain(ax,11,83,11,S.GRAY)
    brain(ax,44,83,11,S.BLUE)
    cells(ax,78,83,11)
    text(ax,16.5,72,"Aorta RNA +\nbrain ATAC",color=S.GRAY)
    text(ax,50,72,"Brain correction\nsame RNA pool",color=S.BLUE)
    text(ax,83.5,72,"BMMC donors\nseparate roles",color=S.TEAL)
    for x in (16.5,50,83.5):
        route(ax,[(x,63),(x,59),(50,59)],arrow=False)
    route(ax,[(50,59),(50,55)])
    matrix(ax,8,44,11,9)
    text(ax,60,48,"Fixed panel / non-self pairs")
    route(ax,[(50,42),(50,38),(25,38),(25,34)])
    route(ax,[(50,38),(77,38),(77,34)])
    graph(ax,8,23,14,10,S.BLUE)
    graph(ax,61,23,14,10,S.TEAL)
    text(ax,38,28,"RNA graph",color=S.BLUE)
    text(ax,87,28,"ATAC\nproxy",color=S.TEAL)
    route(ax,[(25,21),(25,17),(50,17)],arrow=False)
    route(ax,[(77,21),(77,17),(50,17),(50,13)])
    text(ax,50,8,"Degree contrasts → M/D → IUT-BY / Holm")
    return {"schematic":True,"restored_from":"R1 workflow illustration",
            "sources":["historical aorta RNA / brain ATAC","same-pool genuine brain","donor-separated BMMC"],
            "interpretation":"Named IUT-BY/Holm families require valid components; diagrams are illustrative."}


def proximal(ax):
    setup(ax)
    ax.add_patch(Rectangle((15,66),70,24,fc="#ECF2F7",ec="none"))
    ax.plot([2,98],[72,72],color=S.GRAY,lw=.65)
    for x,h in ((6,6),(22,9),(38,15),(55,8),(72,12),(82,7),(95,7)):
        xx=np.linspace(x-3,x+3,35)
        ax.plot(xx,72+h*np.exp(-((xx-x)/1.35)**2),color=S.BLUE if 15<=x<=85 else S.GRAY,lw=.7)
    for x in (15,35,85):ax.plot([x,x],[65,89],color=S.BLUE,lw=.55,ls=(0,(3,2)))
    route(ax,[(15,66),(87,66)],S.BLUE)
    for x,w in ((35,8),(54,10),(72,11)):
        ax.add_patch(Rectangle((x,64),w,4,fc=S.BLUE,ec="none"))
    text(ax,14,95,"TSS − 2 kb"); text(ax,40,95,"TSS"); text(ax,87,95,"TES")
    route(ax,[(27,62),(27,56)],S.BLUE)
    route(ax,[(72,62),(72,56)],S.TEAL)
    for x,color in ((3,S.BLUE),(53,S.TEAL)):
        ax.add_patch(Rectangle((x,38),44,17,fc="#F4F7F8",ec=color,lw=.55))
    text(ax,25,50,"Accessibility",color=S.BLUE)
    text(ax,25,42,"Equal-type mean")
    text(ax,75,50,"Motif OR",color=S.TEAL)
    text(ax,75,42,"500-bp windows")
    route(ax,[(25,37),(25,32),(50,32)],arrow=False)
    route(ax,[(75,37),(75,32),(50,32),(50,27)])
    graph(ax,8,13,16,10,S.TEAL)
    text(ax,63,21,"Directed TF → target weights")
    text(ax,63,12,"Self-edges = 0")
    text(ax,50,3,"Brain panel: 442 TFs × 1,200 genes")
    return {"schematic":True,"restored_from":"R1 proximal-linkage illustration",
            "window":"strand-aware TSS minus 2 kb through gene body",
            "brain_tfs":442,"brain_genes":1200,"causal_ground_truth":False}
