#!/usr/bin/env python3
"""Integration tests for PeerJ package build with typography contract."""
import pytest
import subprocess
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from build_peerj_package import (
    PAPER, PKG, FIG_NAMES, TABLE_FRAGMENTS,
    validate_figure_map, TYPOGRAPHY
)


class TestBuildPipeline:
    """Test build pipeline integration."""

    def test_paper_directory_exists(self):
        """Test that paper directory exists."""
        assert PAPER.exists()
        assert PAPER.is_dir()

    def test_manuscript_exists(self):
        """Test that manuscript.tex exists."""
        manuscript = PAPER / "manuscript.tex"
        assert manuscript.exists()

    def test_references_exist(self):
        """Test that references.bib exists."""
        refs = PAPER / "references.bib"
        assert refs.exists()

    def test_wlpeerj_class_exists(self):
        """Test that wlpeerj.cls exists."""
        cls_file = PAPER / "wlpeerj.cls"
        assert cls_file.exists()

    def test_figure_map_validates(self):
        """Test that figure map validates against manuscript."""
        # Should not raise
        validate_figure_map(str(PAPER / "manuscript.tex"))

    def test_all_figure_fragments_exist(self):
        """Test that all figure .tex fragments exist."""
        for name, _ in FIG_NAMES:
            fig_path = PAPER / "figs" / f"{name}.tex"
            assert fig_path.exists(), f"Missing figure fragment: {name}.tex"

    def test_all_table_fragments_exist(self):
        """Test that all table .tex fragments exist."""
        for name in TABLE_FRAGMENTS:
            table_path = PAPER / "figs" / f"{name}.tex"
            assert table_path.exists(), f"Missing table fragment: {name}.tex"

    def test_typography_contract_accessible(self):
        """Test that typography contract is accessible from build script."""
        assert TYPOGRAPHY.base_font_size == 11
        assert TYPOGRAPHY.target_width_inches == 6.8


class TestManuscriptConsistency:
    """Test manuscript consistency with contract."""

    def test_manuscript_references_all_figures(self):
        """Test that manuscript references all figures in map."""
        manuscript_text = (PAPER / "manuscript.tex").read_text()
        for name, _ in FIG_NAMES:
            assert f"figs/{name}.tex" in manuscript_text, \
                f"Manuscript missing reference to {name}.tex"

    def test_manuscript_references_all_tables(self):
        """Test that manuscript references all tables."""
        manuscript_text = (PAPER / "manuscript.tex").read_text()
        for name in TABLE_FRAGMENTS:
            assert f"figs/{name}.tex" in manuscript_text, \
                f"Manuscript missing reference to {name}.tex"

    def test_no_orphan_figure_references(self):
        """Test that manuscript doesn't reference unlisted figures."""
        import re
        manuscript_text = (PAPER / "manuscript.tex").read_text()
        referenced = set(re.findall(r"figs/(fig[0-9]+_[A-Za-z0-9_]+)\.tex", manuscript_text))
        mapped = {name for name, _ in FIG_NAMES}
        orphans = referenced - mapped
        assert not orphans, f"Manuscript references unlisted figures: {orphans}"


class TestFigureFragments:
    """Test figure fragment structure."""

    def test_figure_fragments_are_tex(self):
        """Test that all figure fragments are .tex files."""
        for name, _ in FIG_NAMES:
            fig_path = PAPER / "figs" / f"{name}.tex"
            content = fig_path.read_text()
            # tikzDevice fragments should contain tikzpicture
            assert "tikzpicture" in content or "pgf" in content, \
                f"{name}.tex doesn't appear to be a tikzDevice output"

    def test_table_fragments_are_tex(self):
        """Test that all table fragments are valid LaTeX."""
        for name in TABLE_FRAGMENTS:
            table_path = PAPER / "figs" / f"{name}.tex"
            content = table_path.read_text()
            # Tables should contain tabular or similar environment
            assert "tabular" in content or "longtable" in content, \
                f"{name}.tex doesn't appear to be a LaTeX table"


class TestBuildScriptImports:
    """Test that build script imports work correctly."""

    def test_can_import_build_script(self):
        """Test that build script can be imported."""
        try:
            import build_peerj_package
            assert True
        except ImportError as e:
            pytest.fail(f"Cannot import build_peerj_package: {e}")

    def test_can_import_typography_contract(self):
        """Test that typography contract can be imported."""
        try:
            from figure_typography import TypographyContract
            assert True
        except ImportError as e:
            pytest.fail(f"Cannot import TypographyContract: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
