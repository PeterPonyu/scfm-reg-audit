import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, V2)

import pbmc_uce_eval_v2 as eval_uce
import pbmc_uce_stats_v2 as stats_uce


class TestPbmcUceCache(unittest.TestCase):
    def test_round_trip_and_provenance_rejection(self):
        genes = ["A", "B", "C"]
        cell_ids = np.array([2, 5], dtype=np.int64)
        co = np.eye(3, dtype=np.float32)
        uce = np.ones((3, 3), dtype=np.float32)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cache.npz"
            eval_uce.write_uce_cache(
                path, co, uce, 3, cell_ids, genes, "manifest", 2,
                "rna", "checkpoint", "esm2",
            )
            loaded = eval_uce.load_uce_cache(
                path, cell_ids, genes, "manifest", 2,
                "rna", "checkpoint", "esm2",
            )
            np.testing.assert_array_equal(loaded[0], co)
            np.testing.assert_array_equal(loaded[1], uce)
            self.assertEqual(loaded[2], 3)
            with self.assertRaisesRegex(ValueError, "provenance mismatch"):
                eval_uce.load_uce_cache(
                    path, cell_ids, genes, "manifest", 2,
                    "changed", "checkpoint", "esm2",
                )

    def test_normalization_is_library_size_invariant(self):
        matrix = np.array([[1, 3], [10, 30]], dtype=np.float64)
        transformed = eval_uce.normalized_log_counts(matrix).toarray()
        np.testing.assert_allclose(transformed[0], transformed[1], rtol=1e-12, atol=1e-12)

    @unittest.skipUnless(
        Path(stats_uce.UCE_PATH).exists() and Path(eval_uce.RNA).exists(),
        "installed PBMC UCE cache or RNA h5ad not present (fresh checkout)",
    )
    def test_formal_graph_loader_accepts_installed_cache(self):
        co, uce = stats_uce.load_uce_graph(stats_uce.UCE_PATH)
        self.assertEqual(co.shape, (1200, 1200))
        self.assertEqual(uce.shape, (1200, 1200))
        self.assertTrue(np.isfinite(co).all())
        self.assertTrue(np.isfinite(uce).all())


if __name__ == "__main__":
    unittest.main(verbosity=2)
