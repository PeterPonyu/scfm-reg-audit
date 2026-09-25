"""Small data fixtures for the independent public release entry points."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src/v2"))
import screg_paths
from split_pbmc_multiome import split_multiome
sys.path.insert(0, str(ROOT / "src/revision_05613"))
import common
import multiplicity


class TestIndependentPaths(unittest.TestCase):
    def test_metadata_is_required_unless_pooled_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {}, clear=True):
            root = Path(tmp)
            with self.assertRaisesRegex(FileNotFoundError, "META_FILE"):
                screg_paths.proxy_metadata_path(root)
            local = root / "data/annotation/atac_cell_meta.csv.gz"
            local.parent.mkdir(parents=True)
            local.write_bytes(b"fixture")
            self.assertEqual(Path(screg_paths.proxy_metadata_path(root)), local)
            for invalid in (str(root / "missing.csv.gz"), ""):
                with mock.patch.dict(os.environ, {"META_FILE": invalid}):
                    with self.assertRaises(FileNotFoundError):
                        screg_paths.proxy_metadata_path(root)
            with mock.patch.dict(os.environ, {"META_FILE": "none"}):
                self.assertIsNone(screg_paths.proxy_metadata_path(root))

    def test_shared_null_does_not_search_sibling_and_honors_override(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {}, clear=True):
            root = Path(tmp) / "project"
            sibling = Path(tmp) / "other_study"
            for rel in ("multiome/pbmc10k_atac.h5ad", "annotation/gene_coords_hg38.tsv"):
                p = sibling / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("fixture")
            coords = root / "data/annotation/gene_coords_hg38.tsv"
            coords.parent.mkdir(parents=True)
            coords.write_text("local")
            explicit = Path(tmp) / "configured.tsv"
            explicit.write_text("explicit")
            with mock.patch.dict(os.environ, {"SCREG_MONOREPO_DATA": str(sibling),
                                               "SCREG_GENE_COORDS": str(explicit)}), \
                 mock.patch.object(common.fpa, "ROOT", str(root)), \
                 mock.patch.object(common.fpa, "COORDS", str(coords)):
                spec = importlib.util.spec_from_file_location("shared_null_fixture", ROOT / "src/v2/fm_vs_baseline_shared_null.py")
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                self.assertEqual(module.ATAC_P, str(root / "data/multiome/pbmc10k_atac.h5ad"))
                self.assertEqual(module.COORDS, str(explicit))
                with self.assertRaisesRegex(FileNotFoundError, "fail-closed"):
                    module.require_path(module.ATAC_P, "PBMC")
                with self.assertRaises(FileNotFoundError):
                    module.require_path(str(sibling), "directory is not an input")

    def test_result_fallback_and_invalid_explicit_input(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {}, clear=True):
            root = Path(tmp)
            public = root / "results/example.public.json"
            public.parent.mkdir()
            public.write_text("{}")
            local = root / "results/v2/example.json"
            with mock.patch.object(common.fpa, "ROOT", str(root)), mock.patch.object(common.fpa, "OUT", str(local.parent)):
                self.assertEqual(common.result_json("example.json", "SCREG_AUDIT_JSON"), public)
                local.parent.mkdir()
                local.write_text("{}")
                self.assertEqual(common.result_json("example.json"), local)
                with mock.patch.dict(os.environ, {"SCREG_AUDIT_JSON": str(root / "missing.json")}):
                    with self.assertRaisesRegex(FileNotFoundError, "SCREG_AUDIT_JSON"):
                        common.result_json("example.json", "SCREG_AUDIT_JSON")

    def test_copied_driver_default_stays_in_project(self):
        self.assertEqual(Path(common.rfa.DATA_ROOT).resolve(), ROOT / "data")
        self.assertEqual(Path(common.ATAC_FILES["brain"]).resolve(), ROOT / "data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")


class TestPBMCConversion(unittest.TestCase):
    def make_input(self, path, include_peaks=True):
        matrix = sp.csc_matrix([[1, 0], [2, 3], [0, 4]], dtype=np.int32)
        with h5py.File(path, "w") as handle:
            group = handle.create_group("matrix")
            for name, value in (("data", matrix.data), ("indices", matrix.indices),
                                ("indptr", matrix.indptr), ("shape", matrix.shape)):
                group.create_dataset(name, data=value)
            group.create_dataset("barcodes", data=[b"cell-a", b"cell-b"])
            features = group.create_group("features")
            features.create_dataset("feature_type", data=[b"Gene Expression", b"Gene Expression", b"Peaks" if include_peaks else b"Gene Expression"])
            features.create_dataset("name", data=[b"GENE", b"GENE", b"peak-name"])
            features.create_dataset("id", data=[b"ENSG1", b"ENSG2", b"chr1:10-20"])

    def test_matched_barcodes_counts_and_feature_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.h5"
            self.make_input(source)
            rna_path, atac_path = split_multiome(source, Path(tmp) / "out")
            rna, atac = ad.read_h5ad(rna_path), ad.read_h5ad(atac_path)
            self.assertEqual(rna.obs_names.tolist(), ["cell-a", "cell-b"])
            self.assertEqual(atac.obs_names.tolist(), rna.obs_names.tolist())
            self.assertEqual(rna.var_names.tolist(), ["GENE", "GENE-1"])
            self.assertEqual(atac.var_names.tolist(), ["chr1:10-20"])
            np.testing.assert_array_equal(rna.X.toarray(), [[1, 2], [0, 3]])
            np.testing.assert_array_equal(atac.X.toarray(), [[0], [4]])

    def test_invalid_inputs_produce_no_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, out = Path(tmp) / "input.h5", Path(tmp) / "out"
            with self.assertRaises(FileNotFoundError):
                split_multiome(source, out)
            self.make_input(source, include_peaks=False)
            with self.assertRaisesRegex(ValueError, "both Gene Expression and Peaks"):
                split_multiome(source, out)
            self.assertFalse(out.exists())


class TestPublicSummaryEntryPoints(unittest.TestCase):
    def assert_numeric_record(self, actual, expected):
        if isinstance(expected, dict):
            self.assertEqual(set(actual), set(expected))
            for key in expected:
                self.assert_numeric_record(actual[key], expected[key])
        elif isinstance(expected, list):
            self.assertEqual(len(actual), len(expected))
            for a, b in zip(actual, expected):
                self.assert_numeric_record(a, b)
        elif isinstance(expected, float):
            self.assertAlmostEqual(actual, expected, delta=1e-12)
        else:
            self.assertEqual(actual, expected)

    def run_summary(self, script):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "check.json"
            env = {k: v for k, v in os.environ.items() if not k.startswith(("SCREG_", "SCFM_"))}
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            proc = subprocess.run([sys.executable, "-B", str(ROOT / f"src/revision_05613/{script}.py"),
                                   "--summary-only", "--output", str(output)],
                                  cwd=tmp, env=env, text=True, capture_output=True, timeout=45)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            return json.loads(output.read_text())

    def test_corrections_known_unsorted_values_and_invalid_values(self):
        p = [0.04, 0.01, 0.03, 0.002]
        for method, expected in (("BH", [0.04, 0.02, 0.04, 0.008]),
                                 ("BY", [1/12, 1/24, 1/12, 1/60]),
                                 ("Holm", [0.06, 0.03, 0.06, 0.008])):
            np.testing.assert_allclose(multiplicity.adjust(p, method), expected, rtol=0, atol=1e-12)
        for invalid in ([float("nan")], [-0.1], [1.1], [[0.2]]):
            with self.assertRaises(ValueError):
                multiplicity.adjust(invalid, "BH")

    def test_public_multiplicity_matches_archived_frozen_results(self):
        actual = self.run_summary("multiplicity")
        expected = json.loads((ROOT / "results/revision_05613/multiplicity_robustness.json").read_text())
        self.assertNotIn("resolution_n9999", actual)
        self.assert_numeric_record(actual["frozen_n999"], expected["frozen_n999"])
        self.assertEqual(actual["inputs"]["fixed_panel_audit_v2_sha256"], common.sha256_file(ROOT / "results/fixed_panel_audit_v2.public.json"))

    def test_public_analytic_mde_matches_archived_analytic_fields(self):
        actual = self.run_summary("mde")
        expected = json.loads((ROOT / "results/revision_05613/mde.json").read_text())
        self.assertFalse(actual["inputs"]["n9999_nulls_used"])
        self.assertEqual(len(actual["rows"]), len(expected["rows"]))
        for row, archived in zip(actual["rows"], expected["rows"]):
            self.assertFalse(any(k.endswith(("_emp", "_n9999")) for k in row))
            self.assert_numeric_record(row, {k: archived[k] for k in row})
        for key in ("summary_ranges", "probe"):
            self.assert_numeric_record(actual[key], expected[key])


if __name__ == "__main__":
    unittest.main()
