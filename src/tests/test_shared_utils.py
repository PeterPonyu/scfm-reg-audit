"""Unit tests for the shared utilities extracted from the audit scripts.

Layer 1 only: pure numpy/JSON behaviour, no cached graphs and no subprocesses.
The contracts pinned here are the ones the audit's published numbers depend on:
the plus-one Monte Carlo p-value, BH q-value order, provenance gating of cached
graphs, and crash-safe writes.
"""
import json
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, SRC)
import audit_utils as au  # noqa: E402
import pair_probe_common as ppc  # noqa: E402
import pbmc_cache  # noqa: E402


class TestHashing(unittest.TestCase):
    def test_chunked_file_hash_matches_whole_file_hash(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "payload.bin")
            payload = os.urandom(3 * 1024 * 1024 + 17)
            with open(path, "wb") as fh:
                fh.write(payload)
            self.assertEqual(au.sha256_file(path, chunk_size=1024),
                             hashlib.sha256(payload).hexdigest())

    def test_array_hash_is_layout_independent(self):
        base = np.arange(12, dtype=np.float64).reshape(3, 4)
        self.assertEqual(au.sha256_array(base), au.sha256_array(base[:, :]))
        self.assertNotEqual(au.sha256_array(base), au.sha256_array(base.T))


class TestAtomicWrites(unittest.TestCase):
    def test_json_write_replaces_and_leaves_no_temporary(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "doc.json")
            au.write_json_atomic(path, {"b": 1, "a": 2})
            self.assertEqual(os.listdir(tmp), ["doc.json"])
            with open(path) as fh:
                text = fh.read()
            self.assertEqual(json.loads(text), {"a": 2, "b": 1})
            self.assertLess(text.index('"a"'), text.index('"b"'))
            self.assertTrue(text.endswith("\n"))

    def test_json_write_rejects_non_finite_values_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "doc.json")
            with self.assertRaises(ValueError):
                au.write_json_atomic(path, {"x": float("nan")})
            self.assertFalse(os.path.exists(path))

    def test_npz_write_keeps_requested_filename(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cache.npz")
            au.write_npz_atomic(path, a=np.arange(3), label=np.asarray("v1"))
            self.assertEqual(os.listdir(tmp), ["cache.npz"])
            with np.load(path, allow_pickle=False) as cache:
                np.testing.assert_array_equal(cache["a"], np.arange(3))
                self.assertEqual(str(cache["label"].item()), "v1")


class TestMonteCarloAndBH(unittest.TestCase):
    def test_plus_one_two_sided_p_value(self):
        null = np.array([-3.0, -1.0, 0.5, 2.0])
        p, count = au.mc_pvalue(null, observed=1.0, n_perm=4)
        self.assertEqual(count, 3)
        self.assertAlmostEqual(p, 4 / 5)

    def test_null_summary_fields(self):
        null = np.array([0.1, -0.2, 0.3])
        summary = au.mc_null_summary(null, observed=0.25, n_perm=3, seed=7)
        self.assertEqual(summary["N_perm"], 3)
        self.assertEqual(summary["seed"], 7)
        self.assertAlmostEqual(summary["resolution"], 0.25)
        self.assertAlmostEqual(
            summary["p_mc"],
            (summary["null_obs_count_at_or_above_obs"] + 1) / 4)
        self.assertAlmostEqual(summary["null_mean"], float(null.mean()))
        self.assertAlmostEqual(summary["z"],
                               (0.25 - null.mean()) / (null.std() + 1e-9))

    def test_bh_is_order_preserving_and_monotone(self):
        pvalues = [0.04, 0.01, 0.9, 0.03]
        q = au.bh_qvalues(pvalues)
        self.assertEqual(len(q), len(pvalues))
        self.assertAlmostEqual(q[1], 0.04)
        self.assertAlmostEqual(q[2], 0.9)
        by_rank = [q[i] for i in np.argsort(pvalues)]
        self.assertEqual(by_rank, sorted(by_rank))
        self.assertTrue(all(0.0 <= value <= 1.0 for value in q))

    def test_bh_rounding_is_opt_in(self):
        self.assertEqual(au.bh_qvalues([1 / 3, 1 / 3], round_to=6), [0.333333, 0.333333])


class TestSeeds(unittest.TestCase):
    def test_spawned_seeds_are_distinct_ints_and_reproducible(self):
        first = au.spawn_int_seeds(np.random.SeedSequence(20260724), 5)
        second = au.spawn_int_seeds(np.random.SeedSequence(20260724), 5)
        self.assertEqual(first, second)
        self.assertEqual(len(set(first)), 5)
        self.assertTrue(all(isinstance(seed, int) for seed in first))


def _write_cache(path, **arrays):
    pbmc_cache.write_graph_cache(path, **arrays)


class TestGraphCacheProvenance(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "graphs.npz")
        self.genes = ["A", "B"]
        self.expectations = {
            "genes": self.genes,
            "manifest_sha": "abc",
            "selection_seed": 11,
            "cell_ids": np.array([3, 5]),
        }
        _write_cache(self.path, co=np.eye(2), sg=np.eye(2), covered=2,
                     genes=self.genes, manifest_sha="abc", selection_seed=11,
                     cell_ids=np.array([3, 5]))

    def tearDown(self):
        self.tmp.cleanup()

    def load(self):
        return pbmc_cache.load_graph_cache(
            self.path, ("co", "sg"), self.expectations, "TEST", len(self.genes),
            scalar_keys=("covered",))

    def test_missing_cache_returns_none(self):
        self.assertIsNone(pbmc_cache.load_graph_cache(
            os.path.join(self.tmp.name, "absent.npz"), ("co",),
            self.expectations, "TEST", 2))

    def test_matching_cache_round_trips(self):
        cache = self.load()
        np.testing.assert_array_equal(cache["co"], np.eye(2))
        self.assertEqual(cache["covered"], 2)

    def test_each_provenance_field_is_gated(self):
        for key, wrong in (("genes", ["A", "C"]), ("manifest_sha", "def"),
                           ("selection_seed", 12), ("cell_ids", np.array([3, 6]))):
            self.expectations[key] = wrong
            with self.assertRaises(ValueError) as raised:
                self.load()
            self.assertIn(key, str(raised.exception))
            self.expectations[key] = {
                "genes": self.genes, "manifest_sha": "abc",
                "selection_seed": 11, "cell_ids": np.array([3, 5])}[key]

    def test_shape_and_finiteness_are_gated(self):
        _write_cache(self.path, co=np.eye(3), sg=np.eye(3), covered=2,
                     genes=self.genes, manifest_sha="abc", selection_seed=11,
                     cell_ids=np.array([3, 5]))
        with self.assertRaises(ValueError):
            self.load()
        broken = np.eye(2)
        broken[0, 0] = np.inf
        _write_cache(self.path, co=broken, sg=np.eye(2), covered=2,
                     genes=self.genes, manifest_sha="abc", selection_seed=11,
                     cell_ids=np.array([3, 5]))
        with self.assertRaises(ValueError):
            self.load()


class TestPairProbeHelpers(unittest.TestCase):
    def test_residualiser_removes_design_span(self):
        rng = np.random.default_rng(0)
        design = np.column_stack([np.ones(50), rng.normal(size=50)])
        residuals = ppc.rank_residualiser(design)(rng.normal(size=50))
        np.testing.assert_allclose(design.T @ residuals, np.zeros(2), atol=1e-8)

    def test_confound_design_is_standardised_with_intercept(self):
        conf = np.column_stack([np.arange(10.0), np.ones(10)])
        design = ppc.confound_design(conf, 10)
        self.assertEqual(design.shape, (10, 3))
        np.testing.assert_allclose(design[:, 0], np.ones(10))
        self.assertAlmostEqual(float(design[:, 1].mean()), 0.0)
        self.assertAlmostEqual(float(design[:, 1].std()), 1.0)
        np.testing.assert_allclose(design[:, 2], np.zeros(10))

    def test_blocks_splits_pairs_per_tf(self):
        block = ppc.blocks(np.arange(6), 2, 3)
        self.assertEqual(block.shape, (2, 3))
        np.testing.assert_array_equal(block[1], [3, 4, 5])

    def test_per_tf_rho_is_nan_for_degenerate_tfs(self):
        y_true = np.array([[1.0, 2.0, 3.0], [1.0, 1.0, 1.0]])
        rhos = ppc.per_tf_rho(y_true, y_true)
        self.assertAlmostEqual(rhos[0], 1.0)
        self.assertTrue(np.isnan(rhos[1]))

    def test_permutation_p_matches_plus_one_formula(self):
        null = np.array([0.1, -0.4, 0.6])
        self.assertAlmostEqual(ppc.permutation_p(null, 0.5, 3), 2 / 4)


if __name__ == "__main__":
    unittest.main()
