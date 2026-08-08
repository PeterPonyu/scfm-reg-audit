"""Unit tests for the TF-disjoint pair-level probe statistics.

Covers src/run_pair_probe.py (rank residualisation, per-TF Spearman readouts,
partial rho given the co-expression probe, summary reduction) and
src/pair_probe_stats.py (BH q-values, the fixed-design residualiser, and the
confound-adjusted per-TF readout the permutation null is built on).

Both modules import scikit-learn at module scope; the capsule's declared test
dependencies do not include it, so these tests skip when it is absent.
"""
import os
import sys
import unittest

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, SRC)

try:
    import pair_probe_stats as pps
    import run_pair_probe as rpp
    SKLEARN_MISSING = None
except ImportError as exc:  # pragma: no cover - environment dependent
    pps = rpp = None
    SKLEARN_MISSING = str(exc)


def _design(n_genes, n_conf=2, seed=0):
    rng = np.random.default_rng(seed)
    return np.column_stack([np.ones(n_genes), rng.normal(size=(n_genes, n_conf))])


def _bh_reference(pvalues):
    """Textbook step-up BH, computed independently of both implementations."""
    n = len(pvalues)
    order = np.argsort(pvalues, kind="stable")
    q = np.empty(n, dtype=float)
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, pvalues[index] * n / rank, 1.0)
        q[index] = running
    return q


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestResidualise(unittest.TestCase):
    def test_residuals_are_orthogonal_to_design(self):
        design = _design(40)
        vec = np.random.default_rng(1).normal(size=40)
        resid = rpp.residualise(vec, design)
        np.testing.assert_allclose(design.T @ resid, np.zeros(design.shape[1]), atol=1e-8)

    def test_vector_enters_only_through_its_ranks(self):
        design = _design(30)
        vec = np.random.default_rng(2).normal(size=30)
        monotone = np.exp(vec) + 3.0
        np.testing.assert_allclose(rpp.residualise(vec, design),
                                   rpp.residualise(monotone, design), atol=1e-9)

    def test_intercept_only_design_centres_ranks(self):
        design = np.ones((5, 1))
        resid = rpp.residualise(np.array([10.0, 20.0, 30.0, 40.0, 50.0]), design)
        np.testing.assert_allclose(resid, np.array([-2.0, -1.0, 0.0, 1.0, 2.0]), atol=1e-12)

    def test_design_column_is_fully_removed(self):
        column = np.arange(20, dtype=float)
        design = np.column_stack([np.ones(20), rankdata(column)])
        np.testing.assert_allclose(rpp.residualise(column, design),
                                   np.zeros(20), atol=1e-8)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestPerTfRho(unittest.TestCase):
    def test_monotone_prediction_scores_one_per_tf(self):
        n_tf, n_genes = 3, 8
        y_true = np.arange(n_tf * n_genes, dtype=float)
        rhos = rpp.per_tf_rho(y_true, np.exp(y_true / 10.0), n_tf, n_genes)
        np.testing.assert_allclose(rhos, np.ones(n_tf))

    def test_reversed_prediction_scores_minus_one(self):
        n_tf, n_genes = 2, 6
        y_true = np.arange(n_tf * n_genes, dtype=float)
        y_pred = np.concatenate([y_true[:n_genes][::-1], y_true[n_genes:][::-1]])
        np.testing.assert_allclose(rpp.per_tf_rho(y_true, y_pred, n_tf, n_genes),
                                   -np.ones(n_tf))

    def test_degenerate_tf_rows_are_nan_not_zero(self):
        n_tf, n_genes = 2, 5
        y_true = np.concatenate([np.arange(n_genes, dtype=float), np.ones(n_genes)])
        y_pred = np.concatenate([np.arange(n_genes, dtype=float), np.arange(n_genes, dtype=float)])
        rhos = rpp.per_tf_rho(y_true, y_pred, n_tf, n_genes)
        self.assertAlmostEqual(rhos[0], 1.0)
        self.assertTrue(np.isnan(rhos[1]))
        rhos = rpp.per_tf_rho(y_pred, y_true, n_tf, n_genes)  # degenerate prediction
        self.assertTrue(np.isnan(rhos[1]))

    def test_design_argument_matches_manual_residual_spearman(self):
        n_tf, n_genes = 4, 25
        rng = np.random.default_rng(3)
        y_true = rng.normal(size=n_tf * n_genes)
        y_pred = rng.normal(size=n_tf * n_genes)
        design = _design(n_genes, seed=4)
        got = rpp.per_tf_rho(y_true, y_pred, n_tf, n_genes, design=design)
        for i in range(n_tf):
            a = rpp.residualise(y_true.reshape(n_tf, n_genes)[i], design)
            b = rpp.residualise(y_pred.reshape(n_tf, n_genes)[i], design)
            self.assertAlmostEqual(got[i], spearmanr(a, b).statistic)

    def test_adjustment_removes_a_confound_driven_association(self):
        """Two vectors that only share the confound must lose their marginal rho."""
        n_genes = 200
        rng = np.random.default_rng(20)
        confound = np.arange(n_genes, dtype=float)
        design = np.column_stack([np.ones(n_genes), rankdata(confound)])
        y_true = confound + rng.normal(scale=1.0, size=n_genes)
        y_pred = confound + rng.normal(scale=1.0, size=n_genes)
        marginal = rpp.per_tf_rho(y_true, y_pred, 1, n_genes)[0]
        adjusted = rpp.per_tf_rho(y_true, y_pred, 1, n_genes, design=design)[0]
        self.assertGreater(marginal, 0.9)
        self.assertLess(abs(adjusted), 0.2)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestPartialRho(unittest.TestCase):
    def test_signal_shared_with_the_baseline_is_removed(self):
        """A probe that only reproduces the co-expression baseline gains nothing."""
        n_genes = 200
        rng = np.random.default_rng(5)
        y_ctrl = np.arange(n_genes, dtype=float)
        y_true = y_ctrl + rng.normal(size=n_genes)
        y_pred = y_ctrl + rng.normal(size=n_genes)
        marginal = rpp.per_tf_rho(y_true, y_pred, 1, n_genes)[0]
        partial = rpp.partial_rho(y_true, y_pred, y_ctrl, 1, n_genes)[0]
        self.assertGreater(marginal, 0.9)
        self.assertLess(abs(partial), 0.2)

    def test_control_orthogonal_prediction_keeps_signal(self):
        n_genes = 40
        y_true = np.arange(n_genes, dtype=float)
        y_ctrl = np.zeros(n_genes)
        y_ctrl[::2] = 1.0
        out = rpp.partial_rho(y_true, y_true, y_ctrl, 1, n_genes)
        self.assertGreater(out[0], 0.9)

    def test_matches_manual_partial_spearman(self):
        n_tf, n_genes = 3, 30
        rng = np.random.default_rng(7)
        y_true = rng.normal(size=n_tf * n_genes)
        y_pred = rng.normal(size=n_tf * n_genes)
        y_ctrl = rng.normal(size=n_tf * n_genes)
        got = rpp.partial_rho(y_true, y_pred, y_ctrl, n_tf, n_genes)
        for i in range(n_tf):
            block = slice(i * n_genes, (i + 1) * n_genes)
            design = np.column_stack([np.ones(n_genes), rankdata(y_ctrl[block])])
            a = rpp.residualise(y_true[block], design)
            b = rpp.residualise(y_pred[block], design)
            self.assertAlmostEqual(got[i], spearmanr(a, b).statistic)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestSummarise(unittest.TestCase):
    def test_ignores_degenerate_tfs_and_reports_their_count(self):
        rhos = np.array([0.2, np.nan, 0.4, 0.6, np.nan])
        got = rpp.summarise("adjusted_rho", rhos)
        self.assertEqual(got["adjusted_rho_n_ok"], 3)
        self.assertAlmostEqual(got["adjusted_rho_mean"], 0.4)
        self.assertAlmostEqual(got["adjusted_rho_median"], 0.4)
        self.assertAlmostEqual(got["adjusted_rho_std"], np.std([0.2, 0.4, 0.6]))

    def test_keys_are_namespaced_and_values_are_plain_python(self):
        got = rpp.summarise("marginal_rho", np.array([0.1, 0.3]))
        self.assertEqual(sorted(got), ["marginal_rho_mean", "marginal_rho_median",
                                       "marginal_rho_n_ok", "marginal_rho_std"])
        self.assertIsInstance(got["marginal_rho_mean"], float)
        self.assertIsInstance(got["marginal_rho_n_ok"], int)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestArmContract(unittest.TestCase):
    def test_edge_only_arm_withholds_the_degree_columns(self):
        self.assertEqual(rpp.PRIMARY_ARM, "edge_only")
        self.assertEqual(set(rpp.ARMS["all"]) - set(rpp.ARMS["edge_only"]), {2, 3})
        self.assertTrue(set(rpp.ARMS["edge_only"]).issubset(rpp.ARMS["all"]))
        self.assertEqual(rpp.BASELINE, pps.BASELINE)

    def test_seed_contract_is_explicit_and_reproducible(self):
        self.assertEqual(pps.N_PERM, 999)
        self.assertIsInstance(pps.SEED_ROOT, int)
        family_seeds = [pps.SEED_ROOT * 1000 + i for i in range(6)]
        self.assertEqual(len(set(family_seeds)), 6)
        self.assertNotIn(pps.SEED_ROOT + 1, family_seeds)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestBenjaminiHochberg(unittest.TestCase):
    def test_matches_reference_step_up_on_random_pvalues(self):
        p = np.random.default_rng(8).uniform(size=50)
        np.testing.assert_allclose(pps.bh(p), _bh_reference(p), atol=1e-12)

    def test_order_is_preserved_and_qs_are_bounded(self):
        p = np.array([0.9, 0.001, 0.5, 0.02])
        q = pps.bh(p)
        self.assertEqual(np.argmin(q), np.argmin(p))
        self.assertTrue(((q >= 0) & (q <= 1)).all())
        self.assertTrue(np.all(q >= p - 1e-12))

    def test_monotone_in_the_sorted_p_order(self):
        p = np.random.default_rng(9).uniform(size=30)
        q = pps.bh(p)[np.argsort(p)]
        self.assertTrue(np.all(np.diff(q) >= -1e-12))

    def test_single_pvalue_is_unchanged(self):
        np.testing.assert_allclose(pps.bh([0.037]), [0.037])

    def test_all_ones_stay_at_one(self):
        np.testing.assert_allclose(pps.bh([1.0, 1.0, 1.0]), np.ones(3))

    def test_agrees_with_the_validator_implementation(self):
        sys.path.insert(0, os.path.abspath(os.path.join(SRC, "..")))
        import validate_artifacts as va
        p = list(np.random.default_rng(10).uniform(size=17))
        np.testing.assert_allclose(pps.bh(p), va.bh(p), atol=1e-12)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestMakeResidualiser(unittest.TestCase):
    def test_matches_the_lstsq_residualiser_used_by_the_probe(self):
        design = _design(35, seed=11)
        resid = pps.make_residualiser(design)
        vec = np.random.default_rng(12).normal(size=35)
        np.testing.assert_allclose(resid(vec), rpp.residualise(vec, design), atol=1e-8)

    def test_residuals_are_orthogonal_to_the_fixed_design(self):
        design = _design(24, n_conf=3, seed=13)
        resid = pps.make_residualiser(design)
        for seed in (14, 15):
            out = resid(np.random.default_rng(seed).normal(size=24))
            np.testing.assert_allclose(design.T @ out, np.zeros(design.shape[1]), atol=1e-8)

    def test_rank_deficient_design_still_projects(self):
        base = np.arange(12, dtype=float)
        design = np.column_stack([np.ones(12), base, 2.0 * base])  # collinear
        resid = pps.make_residualiser(design)
        out = resid(np.random.default_rng(16).normal(size=12))
        np.testing.assert_allclose(design.T @ out, np.zeros(3), atol=1e-8)


@unittest.skipIf(SKLEARN_MISSING, f"scikit-learn unavailable: {SKLEARN_MISSING}")
class TestAdjustedPerTf(unittest.TestCase):
    def setUp(self):
        self.n_tf, self.n_genes = 4, 30
        self.design = _design(self.n_genes, seed=17)
        self.resid = pps.make_residualiser(self.design)

    def test_shape_and_agreement_with_the_probe_readout(self):
        rng = np.random.default_rng(18)
        y_true = rng.normal(size=(self.n_tf, self.n_genes))
        y_pred = rng.normal(size=(self.n_tf, self.n_genes))
        got = pps.adjusted_per_tf(y_true, y_pred, self.resid)
        self.assertEqual(got.shape, (self.n_tf,))
        expected = rpp.per_tf_rho(y_true.ravel(), y_pred.ravel(), self.n_tf,
                                  self.n_genes, design=self.design)
        np.testing.assert_allclose(got, expected, atol=1e-8)

    def test_constant_rows_are_nan_under_the_intercept_residualiser(self):
        resid = pps.make_residualiser(np.ones((self.n_genes, 1)))
        y_true = np.tile(np.arange(self.n_genes, dtype=float), (self.n_tf, 1))
        y_pred = y_true.copy()
        y_pred[2] = 7.0
        got = pps.adjusted_per_tf(y_true, y_pred, resid)
        self.assertTrue(np.isnan(got[2]))
        np.testing.assert_allclose(np.delete(got, 2), np.ones(self.n_tf - 1))

    def test_within_tf_permutation_destroys_the_pairing(self):
        """The permutation null's mean must sit near zero while the observed does not."""
        rng = np.random.default_rng(19)
        y_true = rng.normal(size=(self.n_tf, self.n_genes))
        y_pred = y_true + 0.1 * rng.normal(size=(self.n_tf, self.n_genes))
        observed = np.nanmean(pps.adjusted_per_tf(y_true, y_pred, self.resid))
        null = [np.nanmean(pps.adjusted_per_tf(
            y_true, np.stack([rng.permutation(row) for row in y_pred]), self.resid))
            for _ in range(30)]
        self.assertGreater(observed, 0.8)
        self.assertLess(abs(float(np.mean(null))), 0.2)


if __name__ == "__main__":
    unittest.main()
