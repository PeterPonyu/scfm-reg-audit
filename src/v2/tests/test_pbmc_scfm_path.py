import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, V2)

import pbmc_cache


class TestPbmcCache(unittest.TestCase):
    def setUp(self):
        self.genes = ["A", "B", "C"]
        self.manifest_sha = "abc123"
        self.cell_ids = np.array([4, 7, 9])
        self.co = np.eye(3, dtype=np.float32)
        self.sg = np.full((3, 3), 0.25, dtype=np.float32)
        self.rna_sha256 = "rna123"

    def test_pool_cell_selection_is_deterministic(self):
        first = pbmc_cache.select_pool_cell_ids(10, 4, 20260713)
        second = pbmc_cache.select_pool_cell_ids(10, 4, 20260713)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(first), 4)
        np.testing.assert_array_equal(
            pbmc_cache.select_pool_cell_ids(3, 4, 20260713), np.arange(3))

    def test_scgpt_cache_round_trip_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cache.npz")
            pbmc_cache.write_scgpt_cache(
                path, self.co, self.sg, self.cell_ids, self.genes,
                self.manifest_sha, 20260713, 4000, self.rna_sha256)
            co, sg = pbmc_cache.load_scgpt_cache(
                path, self.cell_ids, self.genes, self.manifest_sha, 20260713, 4000,
                self.rna_sha256)
            np.testing.assert_array_equal(co, self.co)
            np.testing.assert_array_equal(sg, self.sg)
            with self.assertRaisesRegex(ValueError, "provenance mismatch"):
                pbmc_cache.load_scgpt_cache(
                    path, self.cell_ids[::-1], self.genes, self.manifest_sha, 20260713,
                    4000, self.rna_sha256)
            with self.assertRaisesRegex(ValueError, "provenance mismatch"):
                pbmc_cache.load_scgpt_cache(
                    path, self.cell_ids, self.genes, self.manifest_sha, 20260713,
                    4000, "different-rna")

    def test_atomic_write_failure_preserves_existing_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cache.npz")
            pbmc_cache.write_scgpt_cache(
                path, self.co, self.sg, self.cell_ids, self.genes,
                self.manifest_sha, 20260713, 4000, self.rna_sha256)
            with open(path, "rb") as fh:
                before = fh.read()
            with mock.patch.object(pbmc_cache.os, "replace", side_effect=OSError("injected")):
                with self.assertRaisesRegex(OSError, "injected"):
                    pbmc_cache.write_scgpt_cache(
                        path, self.co * 2, self.sg * 2, self.cell_ids, self.genes,
                        self.manifest_sha, 20260713, 4000, self.rna_sha256)
            with open(path, "rb") as fh:
                self.assertEqual(fh.read(), before)
            self.assertEqual(os.listdir(tmp), ["cache.npz"])

    def test_preflight_report_records_versions_hash_and_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference_path = os.path.join(tmp, "reference.json")
            report_path = os.path.join(tmp, "report.json")
            with open(reference_path, "w") as fh:
                fh.write('{"n_pairs": 3}\n')
            pbmc_cache.write_preflight_report(
                report_path, {"n_pairs": 3}, reference_path, self.manifest_sha)
            import hashlib
            import json
            with open(reference_path, "rb") as fh:
                reference_sha = hashlib.sha256(fh.read()).hexdigest()
            with open(report_path) as fh:
                report = json.load(fh)
            self.assertEqual(report["status"], "passed")
            self.assertEqual(report["manifest_sha"], self.manifest_sha)
            self.assertEqual(report["reference_metrics_sha256"], reference_sha)
            self.assertEqual(report["metrics"], {"n_pairs": 3})
            for key in ("python", "numpy", "scipy", "pandas", "anndata"):
                self.assertTrue(report[key])

    def test_reference_metric_gate_rejects_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "reference.json")
            current = {
                "n_pairs": 3,
                "n_tf": 2,
                "tissue": "PBMC",
                "types": ["T"],
                "embed__fm_vs_atac": 0.1,
                "attn__fm_vs_atac": 0.2,
            }
            import json
            with open(path, "w") as fh:
                json.dump(current, fh)
            pbmc_cache.verify_reference_metrics(current, path)
            drifted = dict(current, embed__fm_vs_atac=0.1001)
            with self.assertRaisesRegex(RuntimeError, "reference metric drift"):
                pbmc_cache.verify_reference_metrics(drifted, path)

    def test_confound_cache_validates_manifest_shape_and_finiteness(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "confounds.npz")
            vectors = {key: np.arange(3, dtype=np.float32)
                       for key in ("peakcount", "genelen", "detv", "gc")}
            np.savez(path, genes=np.asarray(self.genes), manifest_sha=np.asarray(self.manifest_sha), **vectors)
            loaded = pbmc_cache.load_confound_cache(path, self.genes, self.manifest_sha)
            for vector in loaded:
                np.testing.assert_array_equal(vector, np.arange(3, dtype=np.float32))
            with self.assertRaisesRegex(ValueError, "manifest mismatch"):
                pbmc_cache.load_confound_cache(path, self.genes, "wrong")
            vectors["gc"] = np.array([0.0, np.nan, 1.0], dtype=np.float32)
            np.savez(path, genes=np.asarray(self.genes), manifest_sha=np.asarray(self.manifest_sha), **vectors)
            with self.assertRaisesRegex(ValueError, "non-finite"):
                pbmc_cache.load_confound_cache(path, self.genes, self.manifest_sha)


class TestPbmcRunnerContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(os.path.join(V2, "pbmc_eval_v2.py")) as fh:
            cls.source = fh.read()

    def test_scgpt_uses_matched_coexpression_and_separate_cache(self):
        self.assertIn('build_pooled_scgpt_graphs(allcells, f"{OUT}/pbmc_scgpt_pooled_v2.npz")', self.source)
        self.assertIn('run_test(a_p, Gco_sg[iiv, jjv], Gsg[iiv, jjv], conf, "scgpt")', self.source)
        self.assertNotIn('sg=Gsg', self.source)

    def test_preflight_exits_before_scgpt_checkpoint_load(self):
        self.assertIn('os.environ.get("PBMC_PREFLIGHT_ONLY", "0") == "1"', self.source)
        gate = self.source.index("pbmc_cache.verify_reference_metrics")
        preflight = self.source.index("if PREFLIGHT_ONLY:")
        scgpt = self.source.index("Gco_sg, Gsg = build_pooled_scgpt_graphs")
        self.assertLess(gate, preflight)
        self.assertLess(preflight, scgpt)
        self.assertIn('raise RuntimeError("PBMC_PREFLIGHT_ONLY requires PBMC_REFERENCE_METRICS")', self.source)
        self.assertIn("raise SystemExit(0)", self.source)

    def test_remote_run_can_skip_per_type_graphs(self):
        self.assertIn('os.environ.get("PBMC_POOLED_ONLY", "0") == "1"', self.source)
        self.assertIn("for t in ([] if POOLED_ONLY else types):", self.source)

    def test_scgpt_cache_is_bound_to_rna_content(self):
        build = self.source[self.source.index("def build_pooled_scgpt_graphs"):]
        self.assertIn("rna_sha256 = pbmc_cache.sha256_file(RNA)", build)
        self.assertIn("POOL_CAP, rna_sha256", build)


if __name__ == "__main__":
    unittest.main()
