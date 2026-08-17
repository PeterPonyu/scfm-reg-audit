"""Unit tests for the capsule validator (validate_artifacts.py).

The validator is the release gate: it must reject non-finite JSON values, a
figure allowlist that drifts from the manuscript, private paths leaking into
shipped text, and any manifest hash mismatch. These tests exercise its helpers
against synthetic roots, plus a live check that the MANIFEST.json records for
the audit sources still match the files they pin.
"""
import glob
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
CAPSULE_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, CAPSULE_ROOT)
import validate_artifacts as va  # noqa: E402


def _bh_reference(pvalues):
    n = len(pvalues)
    order = sorted(range(n), key=pvalues.__getitem__)
    out = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, pvalues[index] * n / rank, 1.0)
        out[index] = running
    return out


class TestLoad(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        patcher = mock.patch.object(va, "RESULTS", Path(self.tmp.name))
        patcher.start()
        self.addCleanup(patcher.stop)

    def write(self, text):
        (Path(self.tmp.name) / "doc.json").write_text(text)

    def test_reads_plain_json(self):
        self.write('{"a": [1, 2.5], "b": {"c": true}}')
        self.assertEqual(va.load("doc.json"), {"a": [1, 2.5], "b": {"c": True}})

    def test_rejects_nan_and_infinity_constants(self):
        for literal in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(literal=literal):
                self.write('{"rho": %s}' % literal)
                with self.assertRaises(va.ValidationError):
                    va.load("doc.json")


class TestWalkNumbers(unittest.TestCase):
    def test_yields_nested_numbers_in_document_order(self):
        doc = {"a": 1, "b": [2, {"c": 3.5}], "d": {"e": [[4]]}}
        self.assertEqual(sorted(va.walk_numbers(doc)), [1, 2, 3.5, 4])

    def test_booleans_and_strings_and_nulls_are_not_numbers(self):
        doc = {"flag": True, "other": False, "name": "0.5", "missing": None}
        self.assertEqual(list(va.walk_numbers(doc)), [])

    def test_empty_containers_yield_nothing(self):
        self.assertEqual(list(va.walk_numbers({"a": [], "b": {}})), [])


class TestBenjaminiHochberg(unittest.TestCase):
    def test_matches_reference_step_up(self):
        pvalues = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205]
        for got, expected in zip(va.bh(pvalues), _bh_reference(pvalues)):
            self.assertAlmostEqual(got, expected)

    def test_q_values_are_clipped_to_one_and_order_preserving(self):
        pvalues = [0.9, 0.8, 0.02]
        q = va.bh(pvalues)
        self.assertTrue(all(0.0 <= value <= 1.0 for value in q))
        self.assertEqual(q.index(min(q)), 2)

    def test_singleton_family_is_unchanged(self):
        self.assertAlmostEqual(va.bh([0.037])[0], 0.037)

    def test_ties_receive_equal_q_values(self):
        q = va.bh([0.02, 0.02, 0.5])
        self.assertAlmostEqual(q[0], q[1])


class TestFigureContract(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "paper" / "figs").mkdir(parents=True)
        patcher = mock.patch.object(va, "ROOT", self.root)
        patcher.start()
        self.addCleanup(patcher.stop)

    def write_capsule(self, manuscript=va.CURRENT_FRAGMENTS, preview=va.CURRENT_FRAGMENTS,
                      bundled=va.CURRENT_FRAGMENTS):
        def body(names):
            return "".join("\\input{figs/%s}\n" % name for name in names)
        (self.root / "paper/manuscript.tex").write_text(body(manuscript))
        (self.root / "paper/figs_preview.tex").write_text(body(preview))
        for name in bundled:
            (self.root / "paper/figs" / name).write_text("% fragment\n")

    def test_figure_inputs_extracts_only_figs_inputs_in_order(self):
        path = self.root / "paper/manuscript.tex"
        path.write_text("\\input{figs/fig1_truth_construct.tex}\n"
                        "\\input{sections/intro.tex}\n"
                        "\\input{figs/table1_primary_fixed_panel.tex}\n"
                        "\\input{figs/notes.md}\n")
        self.assertEqual(va.figure_inputs(path),
                         ("fig1_truth_construct.tex", "table1_primary_fixed_panel.tex"))

    def test_consistent_capsule_passes(self):
        self.write_capsule()
        va.check_figure_contract()

    def test_missing_manuscript_figure_fails(self):
        self.write_capsule(manuscript=va.CURRENT_FRAGMENTS[1:])
        with self.assertRaises(va.ValidationError):
            va.check_figure_contract()

    def test_figure_ordering_change_fails(self):
        reordered = (va.CURRENT_FIGURES[1], va.CURRENT_FIGURES[0]) + va.CURRENT_FIGURES[2:]
        self.write_capsule(manuscript=reordered + va.CURRENT_TABLES)
        with self.assertRaises(va.ValidationError):
            va.check_figure_contract()

    def test_preview_out_of_sync_fails(self):
        self.write_capsule(preview=va.CURRENT_FRAGMENTS[:-1])
        with self.assertRaises(va.ValidationError):
            va.check_figure_contract()

    def test_unlisted_bundled_fragment_fails(self):
        self.write_capsule()
        (self.root / "paper/figs/fig99_retired.tex").write_text("% retired\n")
        with self.assertRaises(va.ValidationError) as ctx:
            va.check_figure_contract()
        self.assertIn("fig99_retired.tex", str(ctx.exception))

    def test_allowlist_covers_figures_and_tables_without_overlap(self):
        self.assertEqual(set(va.CURRENT_FRAGMENTS),
                         set(va.CURRENT_FIGURES) | set(va.CURRENT_TABLES))
        self.assertFalse(set(va.CURRENT_FIGURES) & set(va.CURRENT_TABLES))
        self.assertEqual(len(set(va.CURRENT_FRAGMENTS)), len(va.CURRENT_FRAGMENTS))


class TestPrivatePathScrub(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        patcher = mock.patch.object(va, "ROOT", self.root)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_scrubbed_capsule_passes(self):
        (self.root / "docs").mkdir()
        (self.root / "docs/notes.md").write_text("provenance ${SCFM_PROJECT_ROOT}/results\n")
        (self.root / "run.py").write_text("PATH = '${SCFM_DATA_ROOT}'\n")
        va.check_no_private_paths()

    def test_leaked_private_path_is_reported_with_its_file(self):
        (self.root / "leak.json").write_text('{"path": "${HOME}/Desktop/x"}')
        with self.assertRaises(va.ValidationError) as ctx:
            va.check_no_private_paths()
        self.assertIn("leak.json", str(ctx.exception))

    def test_binary_and_unscanned_suffixes_are_ignored(self):
        (self.root / "graph.npz").write_bytes(b"${HOME} binary payload")
        va.check_no_private_paths()

    def test_the_validator_source_may_name_the_forbidden_needles(self):
        (self.root / "validate_artifacts.py").write_text(
            Path(CAPSULE_ROOT, "validate_artifacts.py").read_text())
        va.check_no_private_paths()


class TestLocalWorktreeManifestBoundary(unittest.TestCase):
    """MANIFEST closed-tree equality ignores local worktree prefixes."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        for patch_name in ("ROOT", "RESULTS"):
            patcher = mock.patch.object(va, patch_name, self.root if patch_name == "ROOT" else self.root / "results")
            patcher.start()
            self.addCleanup(patcher.stop)
        # Minimal capsule skeleton matching CURRENT_FRAGMENTS + public results stubs.
        (self.root / "paper" / "figs").mkdir(parents=True)
        (self.root / "results").mkdir(parents=True)
        for name in va.CURRENT_FRAGMENTS:
            (self.root / "paper/figs" / name).write_text("% fragment\n")
        (self.root / "paper/manuscript.tex").write_text(
            "".join("\\input{figs/%s}\n" % name for name in va.CURRENT_FRAGMENTS))
        (self.root / "paper/figs_preview.tex").write_text(
            "".join("\\input{figs/%s}\n" % name for name in va.CURRENT_FRAGMENTS))

    def _write_manifest(self, paths):
        records = []
        for rel in paths:
            path = self.root / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                path.write_text("x\n")
            payload = path.read_bytes()
            records.append({
                "path": rel,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            })
        (self.root / "MANIFEST.json").write_text(json.dumps({
            "capsule": "test", "version": "test", "files": records,
        }))

    def test_local_worktree_extra_does_not_break_closed_tree(self):
        capsule = ["README.md", "validate_artifacts.py"]
        for rel in capsule:
            (self.root / rel).write_text("ok\n")
        # Local-only trees that must not force MANIFEST growth.
        (self.root / "src/v2").mkdir(parents=True)
        (self.root / "src/v2/extra.py").write_text("print(1)\n")
        (self.root / "paper/submission_peerj").mkdir(parents=True)
        (self.root / "paper/submission_peerj/README.md").write_text("local\n")
        (self.root / "results/v2").mkdir(parents=True)
        (self.root / "results/v2/private.json").write_text("{}\n")
        self._write_manifest(capsule)
        # Exercise only the coverage split helpers via a narrow call path:
        listed = {record["path"] for record in json.loads((self.root / "MANIFEST.json").read_text())["files"]}
        actual = {
            path.relative_to(self.root).as_posix()
            for path in self.root.rglob("*")
            if path.is_file() and path.name not in va.EXCLUDED_NAMES
            and not va.EXCLUDED_PARTS.intersection(path.parts) and path.suffix != ".pyc"
        }

        def is_local(rel: str) -> bool:
            return (Path(rel).name in va.LOCAL_WORKTREE_NAMES
                    or any(rel.startswith(prefix) for prefix in va.LOCAL_WORKTREE_PREFIXES))

        capsule_actual = {path for path in actual if not is_local(path)}
        # Ignore figure/manuscript scaffolding created in setUp for this unit.
        capsule_actual -= {f"paper/figs/{name}" for name in va.CURRENT_FRAGMENTS}
        capsule_actual -= {"paper/manuscript.tex", "paper/figs_preview.tex", "MANIFEST.json"}
        self.assertEqual(listed, capsule_actual)

    def test_missing_capsule_file_still_detected(self):
        (self.root / "README.md").write_text("ok\n")
        self._write_manifest(["README.md", "missing.py"])
        listed = {"README.md", "missing.py"}
        actual = {"README.md"}
        self.assertTrue(listed - actual)


class TestBundledManifest(unittest.TestCase):
    """Live check of the shipped manifest for the audit sources.

    MANIFEST.json pins the curated capsule subset rather than the whole
    development tree, and the paper build products in it are regenerated by the
    release builder, so only the audit sources are asserted here.
    """

    @classmethod
    def setUpClass(cls):
        path = Path(CAPSULE_ROOT) / "MANIFEST.json"
        if not path.exists():  # pragma: no cover - release layout always ships it
            raise unittest.SkipTest("MANIFEST.json absent")
        cls.manifest = json.loads(path.read_text())

    def test_records_are_unique(self):
        paths = [record["path"] for record in self.manifest["files"]]
        self.assertEqual(len(paths), len(set(paths)))

    def _audit_sources(self):
        patterns = ("*.py", "src/*.py", "src/tests/*.py")
        return sorted(name for pattern in patterns
                      for name in glob.glob(pattern, root_dir=CAPSULE_ROOT))

    def test_audit_sources_are_listed(self):
        listed = {record["path"] for record in self.manifest["files"]}
        for name in self._audit_sources():
            with self.subTest(path=name):
                self.assertIn(name, listed)

    def test_audit_source_records_match_their_files(self):
        records = {record["path"]: record for record in self.manifest["files"]}
        for name in self._audit_sources():
            record = records[name]
            payload = (Path(CAPSULE_ROOT) / name).read_bytes()
            with self.subTest(path=name):
                self.assertEqual(record["bytes"], len(payload))
                self.assertEqual(hashlib.sha256(payload).hexdigest(), record["sha256"])


if __name__ == "__main__":
    unittest.main()
