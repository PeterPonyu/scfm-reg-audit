"""Unit tests for the PBMC cache provenance helpers (src/pbmc_cache.py).

Pure unit layer: temporary files only, no real caches, no model inference. The
contracts under test are the ones the audit relies on when it reuses a cached
graph instead of recomputing it — any provenance field that changes must make
the loader refuse the cache rather than silently return stale numbers.
"""
import hashlib
import json
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, SRC)
import pbmc_cache  # noqa: E402

GENES = ["A", "B", "C"]
MANIFEST_SHA = "0" * 64
SEED = 20260713
POOL_CAP = 2
RNA_SHA = "1" * 64


def _confound_cache(path, genes=GENES, manifest_sha=MANIFEST_SHA, vectors=None, n=None):
    n = len(genes) if n is None else n
    peakcount, genelen, detv, gc = vectors or (
        np.arange(n, dtype=np.float32),
        np.full(n, 1000.0, dtype=np.float32),
        np.linspace(0.1, 0.9, n).astype(np.float32),
        np.full(n, 0.45, dtype=np.float32),
    )
    np.savez(path, genes=np.asarray(genes), manifest_sha=np.asarray(manifest_sha),
             peakcount=peakcount, genelen=genelen, detv=detv, gc=gc)
    return path


class TestSelectPoolCellIds(unittest.TestCase):
    def test_no_subsampling_below_cap(self):
        for n_cells in (0, 1, 5):
            ids = pbmc_cache.select_pool_cell_ids(n_cells, 5, SEED)
            np.testing.assert_array_equal(ids, np.arange(n_cells))

    def test_subsamples_without_replacement_above_cap(self):
        ids = pbmc_cache.select_pool_cell_ids(100, 7, SEED)
        self.assertEqual(ids.size, 7)
        self.assertEqual(np.unique(ids).size, 7)
        self.assertTrue(((ids >= 0) & (ids < 100)).all())

    def test_selection_is_seed_deterministic(self):
        first = pbmc_cache.select_pool_cell_ids(100, 7, SEED)
        np.testing.assert_array_equal(first, pbmc_cache.select_pool_cell_ids(100, 7, SEED))
        self.assertFalse(np.array_equal(
            first, pbmc_cache.select_pool_cell_ids(100, 7, SEED + 1)))


class TestSha256File(unittest.TestCase):
    def test_matches_hashlib_and_is_chunk_size_invariant(self):
        payload = os.urandom(4096) + b"tail"
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "blob.bin")
            with open(path, "wb") as fh:
                fh.write(payload)
            expected = hashlib.sha256(payload).hexdigest()
            self.assertEqual(pbmc_cache.sha256_file(path), expected)
            self.assertEqual(pbmc_cache.sha256_file(path, chunk_size=7), expected)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "empty.bin")
            open(path, "wb").close()
            self.assertEqual(pbmc_cache.sha256_file(path), hashlib.sha256(b"").hexdigest())


class TestVerifyReferenceMetrics(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = os.path.join(self.tmp.name, "reference.json")
        self.reference = {"n_pairs": 10, "n_tf": 2, "tissue": "pbmc", "types": ["T", "B"],
                          "embed__fm_vs_coexp": 0.5, "attn__fm_vs_coexp": 0.25}
        with open(self.path, "w") as fh:
            json.dump(self.reference, fh)

    def test_identical_metrics_pass(self):
        pbmc_cache.verify_reference_metrics(dict(self.reference), self.path)

    def test_untracked_extra_keys_are_ignored(self):
        current = dict(self.reference, unrelated_note="anything")
        pbmc_cache.verify_reference_metrics(current, self.path)

    def test_drifted_readout_metric_raises(self):
        current = dict(self.reference, embed__fm_vs_coexp=0.5001)
        with self.assertRaises(RuntimeError) as ctx:
            pbmc_cache.verify_reference_metrics(current, self.path)
        self.assertIn("embed__fm_vs_coexp", str(ctx.exception))

    def test_missing_tracked_key_raises(self):
        current = {k: v for k, v in self.reference.items() if k != "n_pairs"}
        with self.assertRaises(RuntimeError) as ctx:
            pbmc_cache.verify_reference_metrics(current, self.path)
        self.assertIn("n_pairs", str(ctx.exception))

    def test_new_readout_key_absent_from_reference_raises(self):
        current = dict(self.reference, embed__new_readout=0.1)
        with self.assertRaises(RuntimeError):
            pbmc_cache.verify_reference_metrics(current, self.path)


class TestLoadConfoundCache(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = os.path.join(self.tmp.name, "confounds.npz")

    def test_returns_four_vectors_in_contract_order(self):
        _confound_cache(self.path)
        peakcount, genelen, detv, gc = pbmc_cache.load_confound_cache(
            self.path, GENES, MANIFEST_SHA)
        np.testing.assert_allclose(peakcount, [0.0, 1.0, 2.0])
        np.testing.assert_allclose(genelen, [1000.0] * 3)
        self.assertEqual(detv.shape, (3,))
        np.testing.assert_allclose(gc, [0.45] * 3)

    def test_gene_order_mismatch_raises(self):
        _confound_cache(self.path, genes=["A", "C", "B"])
        with self.assertRaises(ValueError):
            pbmc_cache.load_confound_cache(self.path, GENES, MANIFEST_SHA)

    def test_manifest_sha_mismatch_raises(self):
        _confound_cache(self.path, manifest_sha="9" * 64)
        with self.assertRaises(ValueError):
            pbmc_cache.load_confound_cache(self.path, GENES, MANIFEST_SHA)

    def test_vector_shape_mismatch_raises(self):
        _confound_cache(self.path, n=4)
        with self.assertRaises(ValueError):
            pbmc_cache.load_confound_cache(self.path, GENES, MANIFEST_SHA)

    def test_non_finite_values_raise(self):
        _confound_cache(self.path, vectors=(
            np.array([0.0, np.nan, 2.0], dtype=np.float32),
            np.ones(3, dtype=np.float32),
            np.ones(3, dtype=np.float32),
            np.ones(3, dtype=np.float32),
        ))
        with self.assertRaises(ValueError):
            pbmc_cache.load_confound_cache(self.path, GENES, MANIFEST_SHA)


class TestScgptCacheRoundTrip(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = os.path.join(self.tmp.name, "scgpt.npz")
        self.cell_ids = np.arange(POOL_CAP)
        self.co = np.arange(9, dtype=np.float32).reshape(3, 3)
        self.sg = self.co.T.copy()

    def write(self, co=None, sg=None, **overrides):
        kwargs = {"cell_ids": self.cell_ids, "genes": GENES, "manifest_sha": MANIFEST_SHA,
                  "selection_seed": SEED, "pool_cap": POOL_CAP, "rna_sha256": RNA_SHA}
        kwargs.update(overrides)
        pbmc_cache.write_scgpt_cache(
            self.path, self.co if co is None else co, self.sg if sg is None else sg, **kwargs)

    def load(self, **overrides):
        kwargs = {"cell_ids": self.cell_ids, "genes": GENES, "manifest_sha": MANIFEST_SHA,
                  "selection_seed": SEED, "pool_cap": POOL_CAP, "rna_sha256": RNA_SHA}
        kwargs.update(overrides)
        return pbmc_cache.load_scgpt_cache(self.path, **kwargs)

    def test_missing_cache_returns_none(self):
        self.assertIsNone(self.load())

    def test_write_then_load_returns_copies_of_graphs(self):
        self.write()
        co, sg = self.load()
        np.testing.assert_allclose(co, self.co)
        np.testing.assert_allclose(sg, self.sg)
        co[0, 0] = -1.0
        np.testing.assert_allclose(self.load()[0], self.co)

    def test_write_leaves_no_temporary_files(self):
        self.write()
        self.assertEqual(os.listdir(self.tmp.name), [os.path.basename(self.path)])

    def test_every_provenance_field_is_enforced(self):
        self.write()
        for override in ({"cell_ids": np.arange(POOL_CAP + 1)},
                         {"genes": ["A", "B", "D"]},
                         {"manifest_sha": "9" * 64},
                         {"selection_seed": SEED + 1},
                         {"pool_cap": POOL_CAP + 1},
                         {"rna_sha256": "2" * 64}):
            with self.subTest(field=next(iter(override))):
                with self.assertRaises(ValueError):
                    self.load(**override)

    def test_graph_shape_mismatch_raises(self):
        self.write(co=np.zeros((3, 4), dtype=np.float32), sg=np.zeros((3, 4), dtype=np.float32))
        with self.assertRaises(ValueError):
            self.load()

    def test_non_finite_graph_raises(self):
        bad = self.co.copy()
        bad[1, 1] = np.inf
        self.write(co=bad)
        with self.assertRaises(ValueError):
            self.load()


class TestWritePreflightReport(unittest.TestCase):
    def test_report_records_environment_and_reference_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference_path = os.path.join(tmp, "reference.json")
            with open(reference_path, "w") as fh:
                fh.write("{}")
            report_path = os.path.join(tmp, "preflight.json")
            metrics = {"n_pairs": 10}
            pbmc_cache.write_preflight_report(report_path, metrics, reference_path, MANIFEST_SHA)

            with open(report_path) as fh:
                report = json.load(fh)
            self.assertEqual(report["status"], "passed")
            self.assertEqual(report["metrics"], metrics)
            self.assertEqual(report["manifest_sha"], MANIFEST_SHA)
            self.assertEqual(report["reference_metrics_sha256"],
                             hashlib.sha256(b"{}").hexdigest())
            self.assertEqual(report["numpy"], np.__version__)
            for key in ("python", "scipy", "pandas", "anndata"):
                self.assertTrue(report[key])
            self.assertEqual(sorted(os.listdir(tmp)),
                             ["preflight.json", "reference.json"])

    def test_report_is_replaced_atomically_on_rewrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference_path = os.path.join(tmp, "reference.json")
            with open(reference_path, "w") as fh:
                fh.write("{}")
            report_path = os.path.join(tmp, "preflight.json")
            pbmc_cache.write_preflight_report(report_path, {"n_pairs": 1},
                                              reference_path, MANIFEST_SHA)
            pbmc_cache.write_preflight_report(report_path, {"n_pairs": 2},
                                              reference_path, MANIFEST_SHA)
            with open(report_path) as fh:
                self.assertEqual(json.load(fh)["metrics"], {"n_pairs": 2})
            self.assertEqual(sorted(os.listdir(tmp)),
                             ["preflight.json", "reference.json"])


if __name__ == "__main__":
    unittest.main()
