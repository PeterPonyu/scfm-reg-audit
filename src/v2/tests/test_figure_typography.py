#!/usr/bin/env python3
"""Tests for figure typography and layout contract."""
import pytest
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from figure_typography import (
    TypographyContract,
    LayoutMeasurements,
    FigureComposite,
    FIGURE_MAP,
    FIGURE_COMPOSITES,
    validate_figure_map,
    compute_forest_plot_margin,
)


class TestTypographyContract:
    """Test typography contract immutability and calculations."""

    def test_default_values(self):
        """Test default contract values match specification."""
        typo = TypographyContract()
        assert typo.base_font_size == 11
        assert typo.manuscript_font_size == 11
        assert typo.latex_engine == "pdftex"
        assert typo.target_width_inches == 6.8
        assert typo.peerj_text_block_inches == 6.5

    def test_font_packages(self):
        """Test font packages include newtx families."""
        typo = TypographyContract()
        packages = typo.font_packages_list()
        assert "\\usepackage{newtxtext}" in packages
        assert "\\usepackage{newtxmath}" in packages
        assert "\\usepackage{amsmath}" in packages

    def test_latex_document_declaration(self):
        """Test LaTeX document class declaration."""
        typo = TypographyContract()
        decl = typo.latex_document_declaration()
        assert "11pt" in decl
        assert "article" in decl

    def test_downscale_factor(self):
        """Test expected downscale factor calculation."""
        typo = TypographyContract()
        factor = typo.expected_downscale_factor()
        expected = 6.5 / 6.8
        assert abs(factor - expected) < 0.001

    def test_immutability(self):
        """Test that TypographyContract is immutable."""
        typo = TypographyContract()
        with pytest.raises(AttributeError):
            typo.base_font_size = 12  # type: ignore


class TestLayoutMeasurements:
    """Test layout measurement calculations."""

    def test_compute_right_margin_default(self):
        """Test right margin computation with default padding."""
        layout = LayoutMeasurements()
        margin = layout.compute_right_margin_for_text_column(72.0)
        assert margin == 72.0 + 12.0  # default padding

    def test_compute_right_margin_custom_padding(self):
        """Test right margin computation with custom padding."""
        layout = LayoutMeasurements()
        margin = layout.compute_right_margin_for_text_column(72.0, padding_pt=20.0)
        assert margin == 72.0 + 20.0

    def test_forest_plot_margin_legacy_values(self):
        """Test forest plot margin matches legacy empirical values."""
        # Legacy fig3 panel A: margin(4, 118, 3, 4)
        # This was ~72pt text + ~46pt padding (empirical)
        margin = compute_forest_plot_margin(max_qvalue_text_width_pt=72.0, padding_pt=46.0)
        assert margin == 118.0

        # Legacy fig9 panel A: margin(4, 132, 3, 4)
        # This was ~86pt text + ~46pt padding (empirical)
        margin = compute_forest_plot_margin(max_qvalue_text_width_pt=86.0, padding_pt=46.0)
        assert margin == 132.0


class TestFigureComposite:
    """Test figure composite structure and validation."""

    def test_valid_composite(self):
        """Test creating a valid composite."""
        comp = FigureComposite("fig1", ("A", "B", "C", "D"), legend_owner="A")
        assert comp.figure_id == "fig1"
        assert comp.panels == ("A", "B", "C", "D")
        assert comp.legend_owner == "A"

    def test_invalid_legend_owner(self):
        """Test that invalid legend owner raises error."""
        with pytest.raises(ValueError, match="Legend owner 'Z' not in panels"):
            FigureComposite("fig1", ("A", "B", "C", "D"), legend_owner="Z")

    def test_validate_single_legend(self):
        """Test validation passes with single legend owner."""
        comp = FigureComposite("fig1", ("A", "B", "C", "D"), legend_owner="A")
        measurements = {
            "A": LayoutMeasurements(legend_owner="A"),
            "B": LayoutMeasurements(),
            "C": LayoutMeasurements(),
            "D": LayoutMeasurements(),
        }
        # Should not raise
        comp.validate_legend_ownership(measurements)

    def test_validate_multiple_legends_fails(self):
        """Test validation fails with multiple legend owners."""
        comp = FigureComposite("fig1", ("A", "B", "C", "D"), legend_owner="A")
        measurements = {
            "A": LayoutMeasurements(legend_owner="A"),
            "B": LayoutMeasurements(legend_owner="B"),  # conflict
            "C": LayoutMeasurements(),
            "D": LayoutMeasurements(),
        }
        with pytest.raises(ValueError, match="Multiple panels claim legend ownership"):
            comp.validate_legend_ownership(measurements)

    def test_validate_legend_mismatch(self):
        """Test validation fails when declared and actual owners differ."""
        comp = FigureComposite("fig1", ("A", "B", "C", "D"), legend_owner="A")
        measurements = {
            "A": LayoutMeasurements(),
            "B": LayoutMeasurements(legend_owner="B"),  # actual owner differs
            "C": LayoutMeasurements(),
            "D": LayoutMeasurements(),
        }
        with pytest.raises(ValueError, match="Legend owner mismatch"):
            comp.validate_legend_ownership(measurements)


class TestFigureMap:
    """Test figure map consistency."""

    def test_figure_map_structure(self):
        """Test figure map has expected structure."""
        assert len(FIGURE_MAP) == 12
        for r_name, sub_name in FIGURE_MAP:
            assert r_name.startswith("fig")
            assert sub_name.startswith("Figure")
            assert sub_name[6:].isdigit()

    def test_figure_map_unique_names(self):
        """Test figure map has unique R and submission names."""
        r_names = [name for name, _ in FIGURE_MAP]
        sub_names = [name for _, name in FIGURE_MAP]
        assert len(r_names) == len(set(r_names))
        assert len(sub_names) == len(set(sub_names))

    def test_figure_composites_match_map(self):
        """Test that all composites have entries in the map."""
        map_names = {name for name, _ in FIGURE_MAP}
        composite_names = set(FIGURE_COMPOSITES.keys())
        assert map_names == composite_names

    def test_all_composites_have_declared_panels(self):
        """Most figures are 2x2 (A–D); fig12 is protocol-pass (A) + scope (B)."""
        for name, comp in FIGURE_COMPOSITES.items():
            assert len(comp.panels) >= 1
            if name == "fig12_protocol_pass_matrix":
                assert comp.panels == ("A", "B")
            else:
                assert comp.panels == ("A", "B", "C", "D")


class TestFigureMapValidation:
    """Test figure map validation against manuscript."""

    def test_validate_requires_all_figures(self, tmp_path):
        """Test validation fails if manuscript is missing figures."""
        manuscript = tmp_path / "manuscript.tex"
        # Only reference first figure
        manuscript.write_text(r"\input{figs/fig1_truth_construct.tex}")

        with pytest.raises(RuntimeError, match="Figure map mismatch"):
            validate_figure_map(str(manuscript))

    def test_validate_rejects_extra_figures(self, tmp_path):
        """Test validation fails if manuscript has extra figures."""
        manuscript = tmp_path / "manuscript.tex"
        # Reference all figures plus a fake one
        content = "\n".join(
            [rf"\input{{figs/{name}.tex}}" for name, _ in FIGURE_MAP]
        )
        content += r"\input{figs/fig99_fake.tex}"
        manuscript.write_text(content)

        with pytest.raises(RuntimeError, match="Figure map mismatch"):
            validate_figure_map(str(manuscript))

    def test_validate_accepts_exact_match(self, tmp_path):
        """Test validation passes when manuscript matches map exactly."""
        manuscript = tmp_path / "manuscript.tex"
        content = "\n".join(
            [rf"\input{{figs/{name}.tex}}" for name, _ in FIGURE_MAP]
        )
        manuscript.write_text(content)

        # Should not raise
        validate_figure_map(str(manuscript))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
