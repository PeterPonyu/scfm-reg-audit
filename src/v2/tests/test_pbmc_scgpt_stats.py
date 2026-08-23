import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, V2)

import pbmc_scgpt_stats_v2 as stats


class TestPbmcScgptStatsEntryPoint(unittest.TestCase):
    def test_requires_scgpt_graph(self):
        graph = np.zeros((3, 3), dtype=np.float32)
        with mock.patch.object(stats.drv, "load_pooled_pbmc", return_value=(
                graph, graph, {"geneformer_embed": graph}, np.array([0]), ["T"])), \
             mock.patch.object(stats.drv, "load_optional_pbmc_scgpt", return_value=None):
            with self.assertRaisesRegex(FileNotFoundError, "required PBMC scGPT graph missing"):
                stats.main()

    def test_atomic_output_failure_preserves_existing_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = os.path.join(tmp, "stats.json")
            with open(output, "w") as fh:
                fh.write("old-result\n")
            with mock.patch.object(stats.os, "replace", side_effect=OSError("injected")):
                with self.assertRaisesRegex(OSError, "injected"):
                    stats.write_json_atomic(output, {"new": True})
            with open(output) as fh:
                self.assertEqual(fh.read(), "old-result\n")
            self.assertEqual(os.listdir(tmp), ["stats.json"])

    def test_writes_independent_output(self):
        graph = np.zeros((3, 3), dtype=np.float32)
        base_result = {
            "rows": [],
            "primary_family": {"rows": [], "n_rows": 0},
            "sensitivity_family": {"rows": [], "n_rows": 0},
            "provenance": {"model_label_provenance": {}},
        }
        with tempfile.TemporaryDirectory() as tmp:
            output = os.path.join(tmp, "stats.json")
            with mock.patch.object(stats, "OUT_PATH", output), \
                 mock.patch.object(stats.drv, "load_pooled_pbmc", return_value=(
                     graph, graph, {"geneformer_embed": graph}, np.array([0]), ["T"])), \
                 mock.patch.object(stats.drv, "load_optional_pbmc_scgpt", return_value=(
                     "scgpt.npz", graph, graph)), \
                 mock.patch.object(stats.drv, "run_pooled_family", return_value=base_result), \
                 mock.patch.object(stats.drv, "append_independent_control_model") as append:
                stats.main()
            append.assert_called_once()
            self.assertTrue(os.path.exists(output))
            import json
            with open(output) as fh:
                document = json.load(fh)
            self.assertEqual(document["design"], "PBMC pooled fixed-panel update with matched-control scGPT")
            self.assertEqual(document["types"], ["T"])
            self.assertEqual(document["pooled_pbmc"], base_result)


if __name__ == "__main__":
    unittest.main()
