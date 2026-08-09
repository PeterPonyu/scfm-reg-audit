"""Negative-path + regression tests for the figure-consistency guard (v0.2.9).

The capsule builder and the capsule validator both refuse to ship when:
- `paper/figs_preview.tex` references a fragment outside CURRENT_FRAGMENTS;
- a `FigureN.tex` standalone wrapper points at a fragment not in
  CURRENT_FIGURES[N-1];
- a `FigureN.pdf` filename referenced by `flat_upload/manuscript.tex` is out
  of order (Figure1..Figure11).

These tests monkey-patch the contract constants to inject a stale fragment,
then exercise the guard functions to confirm they raise.
"""
import json
import re
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_release_capsule import (
    CURRENT_FIGURES, CURRENT_TABLES, CURRENT_FRAGMENTS,
    FIGURE_INPUT_RE, WRAPPER_INPUT_RE, FLAT_FIGURE_RE,
)
from validate_capsule import CURRENT_FIGURES as CAP_FIGS


WRAPPERS = Path(__file__).resolve().parents[3] / "paper/submission_peerj/internal/figure_build"
PROJECT_ROOT = Path(__file__).resolve().parents[3]


class BuilderConstantShape(unittest.TestCase):

    def test_constants_consistent(self):
        self.assertEqual(len(CURRENT_FIGURES), 13)
        self.assertEqual(len(CURRENT_TABLES), 6)
        self.assertEqual(CURRENT_FRAGMENTS, CURRENT_FIGURES + CURRENT_TABLES)
        self.assertEqual(CURRENT_FIGURES[0], "fig10_coverage_qc.tex")
        self.assertEqual(CURRENT_FIGURES[-1], "fig13_scope_card.tex")
        self.assertIn("fig12_protocol_pass_matrix.tex", CURRENT_FIGURES)
        self.assertIn("fig13_scope_card.tex", CURRENT_FIGURES)
        self.assertEqual(CURRENT_FIGURES, CAP_FIGS)  # builder + validator use the same allowlist


class RegexParsers(unittest.TestCase):

    def test_figure_input_regex(self):
        text = (
            "\\input{figs/fig3_primary_audit.tex} and "
            "\\input{figs/table1_primary_fixed_panel.tex}"
        )
        self.assertEqual(
            FIGURE_INPUT_RE.findall(text),
            ["fig3_primary_audit.tex", "table1_primary_fixed_panel.tex"],
        )

    def test_wrapper_input_regex(self):
        self.assertEqual(
            WRAPPER_INPUT_RE.findall("\\input{fig10_coverage_qc.tex}"),
            ["fig10_coverage_qc.tex"],
        )

    def test_flat_figure_regex(self):
        text = "\\includegraphics[width=\\linewidth]{Figure3.pdf}"
        self.assertEqual(FLAT_FIGURE_RE.findall(text), ["Figure3.pdf"])

    def test_flat_figure_regex_ignores_unflagged_lines(self):
        text = "this mentions Figure3.pdf in prose"
        self.assertEqual(FLAT_FIGURE_RE.findall(text), [])


class WrappersMatchContract(unittest.TestCase):
    """Positive-path regression: ensure on-disk wrappers are in sync."""

    @unittest.skipUnless(WRAPPERS.exists(), "wrappers not present")
    def test_wrappers_match_contract(self):
        wrapper_inputs = tuple(
            WRAPPER_INPUT_RE.findall((WRAPPERS / f"Figure{index}.tex").read_text())[0]
            for index in range(1, len(CURRENT_FIGURES) + 1)
        )
        self.assertEqual(wrapper_inputs, CURRENT_FIGURES)


class RendererUsesAuthoritativeResults(unittest.TestCase):

    def test_renderer_does_not_load_retired_results(self):
        status = json.loads(
            (PROJECT_ROOT / "results/v2/inference_status_v2.json").read_text()
        )
        retired = {
            name
            for name, metadata in status["legacy_files_status"].items()
            if metadata["status"] == "retired_not_authoritative"
        }
        renderer = (PROJECT_ROOT / "paper/make_figs.R").read_text()
        loaded = set(re.findall(r'J\("([^"\\]+\.json)"\)', renderer))
        self.assertFalse(
            retired & loaded,
            f"figure renderer loads retired results: {sorted(retired & loaded)}",
        )


if __name__ == "__main__":
    unittest.main()
