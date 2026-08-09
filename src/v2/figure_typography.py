#!/usr/bin/env python3
"""Figure typography and layout contract for scfm-reg-audit paper.

Defines explicit, measurable layout constraints that unify the scReg-Eval
manuscript figures and standalone figure builds. All dimensions are in points
unless noted otherwise.

Typography contract:
- Base font size: 11pt (manuscript body is 11pt)
- Font packages: newtxtext, newtxmath (Times-compatible, matches NimbusRomNo9L in PeerJ PDF)
- LaTeX engine: pdftex (for tikzDevice compatibility)
- Smallest rendered text: ~10pt after \fitfig downscale to 6.5in PeerJ text block

Layout contract:
- Figure width: 6.8in design width (downscales to ~6.5in in manuscript)
- Panel margins: computed from content, not hardcoded
- Legend ownership: one legend per composite (enforced structurally)
"""

from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class TypographyContract:
    """Immutable typography specification."""

    base_font_size: int = 11
    manuscript_font_size: int = 11
    font_packages: Tuple[str, ...] = ("\\usepackage{amsmath}",
                                       "\\usepackage{newtxtext}",
                                       "\\usepackage{newtxmath}")
    latex_engine: str = "pdftex"
    target_width_inches: float = 6.8
    peerj_text_block_inches: float = 6.5
    min_rendered_size_pt: float = 10.0

    def latex_document_declaration(self) -> str:
        """Return the LaTeX document class declaration."""
        return f"\\documentclass[{self.base_font_size}pt]{{article}}"

    def font_packages_list(self) -> Tuple[str, ...]:
        """Return the font packages for tikzDevice metrics and rendering."""
        return self.font_packages

    def expected_downscale_factor(self) -> float:
        """Return the expected downscale when fitting to PeerJ text block."""
        return self.peerj_text_block_inches / self.target_width_inches


@dataclass
class LayoutMeasurements:
    """Measured layout dimensions for a figure panel.

    Margins are computed to accommodate the actual content (axis labels, tick
    labels, titles) rather than using fixed values. All measurements are in
    points unless specified.
    """

    # Content measurements (computed from actual text/labels)
    max_ytick_label_width_pt: Optional[float] = None
    max_xtick_label_height_pt: Optional[float] = None
    ylabel_width_pt: Optional[float] = None
    xlabel_height_pt: Optional[float] = None
    title_height_pt: Optional[float] = None

    # Legend measurements
    legend_width_pt: Optional[float] = None
    legend_height_pt: Optional[float] = None
    legend_owner: Optional[str] = None  # Panel ID that owns this legend

    # Computed margins (right margin for forest plots with q-value columns)
    right_margin_pt: Optional[float] = None

    def compute_right_margin_for_text_column(
        self,
        column_text_width_pt: float,
        padding_pt: float = 12.0
    ) -> float:
        """Compute right margin needed to fit a text column (e.g., q-values).

        Args:
            column_text_width_pt: Width of the widest text in the column
            padding_pt: Padding between plot region and text column

        Returns:
            Right margin in points
        """
        return column_text_width_pt + padding_pt


@dataclass
class FigureComposite:
    """A composite figure (e.g., Fig3 with panels A, B, C, D).

    Enforces structural constraints:
    - Exactly one legend owner per composite
    - All panels share the same typography contract
    """

    figure_id: str
    panels: Tuple[str, ...]
    legend_owner: Optional[str] = None
    typography: TypographyContract = TypographyContract()

    def __post_init__(self):
        """Validate legend ownership."""
        if self.legend_owner is not None and self.legend_owner not in self.panels:
            raise ValueError(
                f"Legend owner '{self.legend_owner}' not in panels {self.panels}"
            )

    def validate_legend_ownership(self, measurements: dict) -> None:
        """Validate that only the designated panel owns a legend.

        Args:
            measurements: Dict mapping panel IDs to LayoutMeasurements

        Raises:
            ValueError: If multiple panels claim legend ownership
        """
        owners = [
            panel_id
            for panel_id in self.panels
            if measurements.get(panel_id) and measurements[panel_id].legend_owner
        ]
        if len(owners) > 1:
            raise ValueError(
                f"Multiple panels claim legend ownership in {self.figure_id}: {owners}"
            )
        if owners and self.legend_owner and owners[0] != self.legend_owner:
            raise ValueError(
                f"Legend owner mismatch in {self.figure_id}: "
                f"declared={self.legend_owner}, actual={owners[0]}"
            )


# Figure map: canonical R figure names -> submission package names
FIGURE_MAP = (
    ("fig10_coverage_qc", "Figure1"),
    ("fig1_truth_construct", "Figure2"),
    ("fig2_cross_tissue_decomp", "Figure3"),
    ("fig3_primary_audit", "Figure4"),
    ("fig4_usability_check", "Figure5"),
    ("fig5_null_diagnostics", "Figure6"),
    ("fig6_spec_sensitivity", "Figure7"),
    ("fig7_pertype_descriptive", "Figure8"),
    ("fig8_injection_ladder", "Figure9"),
    ("fig9_tf_probe", "Figure10"),
    ("fig11_third_tissue_transfer", "Figure11"),
    ("fig12_protocol_pass_matrix", "Figure12"),
    ("fig13_scope_card", "Figure13"),
)

# Composite structure for each figure (panel tags)
FIGURE_COMPOSITES = {
    "fig1_truth_construct": FigureComposite("fig1", ("A", "B", "C", "D")),
    "fig2_cross_tissue_decomp": FigureComposite("fig2", ("A", "B", "C", "D")),
    "fig3_primary_audit": FigureComposite("fig3", ("A", "B", "C", "D")),
    "fig4_usability_check": FigureComposite("fig4", ("A", "B", "C", "D")),
    "fig5_null_diagnostics": FigureComposite("fig5", ("A", "B", "C", "D")),
    "fig6_spec_sensitivity": FigureComposite("fig6", ("A", "B", "C", "D")),
    "fig7_pertype_descriptive": FigureComposite("fig7", ("A", "B", "C", "D")),
    "fig8_injection_ladder": FigureComposite("fig8", ("A", "B", "C", "D")),
    "fig9_tf_probe": FigureComposite("fig9", ("A", "B", "C", "D")),
    "fig10_coverage_qc": FigureComposite("fig10", ("A", "B", "C", "D")),
    "fig11_third_tissue_transfer": FigureComposite("fig11", ("A", "B", "C", "D")),
    # Single-panel gate heatmap (0x10); not a 2x2 multipanel.
    "fig12_protocol_pass_matrix": FigureComposite("fig12", ("A",)),
    "fig13_scope_card": FigureComposite("fig13", ("A",)),
}


def validate_figure_map(tex_path: str) -> None:
    """Validate that manuscript references match FIGURE_MAP.

    Args:
        tex_path: Path to manuscript.tex

    Raises:
        RuntimeError: If referenced figures don't match the map
    """
    import re
    from pathlib import Path

    text = Path(tex_path).read_text()
    referenced = set(re.findall(r"figs/(fig[0-9]+_[A-Za-z0-9_]+)\.tex", text))
    mapped = {name for name, _ in FIGURE_MAP}

    if referenced != mapped:
        raise RuntimeError(
            f"Figure map mismatch: "
            f"manuscript-only={sorted(referenced - mapped)}; "
            f"map-only={sorted(mapped - referenced)}"
        )


def compute_forest_plot_margin(
    max_qvalue_text_width_pt: float = 72.0,
    padding_pt: float = 12.0
) -> float:
    """Compute right margin for forest plots with q-value columns.

    Default width (72pt) accommodates "$q_M=0.000;\\ q_D=0.000$" at 11pt base.

    Args:
        max_qvalue_text_width_pt: Width of widest q-value text
        padding_pt: Padding between plot and text

    Returns:
        Right margin in points
    """
    # Legacy values for reference:
    # - fig3_primary_audit panel A: margin(4, 118, 3, 4)
    # - fig9_tf_probe panel A: margin(4, 132, 3, 4)
    # These were measured empirically for specific text widths
    return max_qvalue_text_width_pt + padding_pt
