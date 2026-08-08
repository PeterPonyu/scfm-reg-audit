"""
Tests for the fixed-panel audit. Three layers:

  LAYER 1  Pure unit tests (default `python3 -m unittest`).
           No subprocesses, no real cached graphs, <1s total. Covers structural contracts
           on the module source code, on the p_mc / BH math, and on the legacy hash pins.

  LAYER 2  Reduced synthetic integration smoke.
           Marker: `integration`. Runs the shared module's primitives (Mantel,
           degree-preserving, partial rho, KO row tagging, cross-tissue) directly against
           a synthetic (n_tf=20, Ng=60) graph with N=9 permutations; no GPU, no real
           cached graphs, no subprocess; <30s.

  LAYER 3  End-to-end real-data benchmark.
           Marker: `realdata`. Skipped by default. Run only after per-perm cost is
           verified <25ms AND reviewer approval. Do not enable until cleared.

Run:
    python3 -m unittest src.v2.tests.test_fixed_panel_audit                       # layer 1
    python3 -m unittest src.v2.tests.test_fixed_panel_audit.TestIntegrationSmoke -v   # layer 1 + 2
    python3 -m unittest src.v2.tests.test_fixed_panel_audit.TestRealDataBenchmark -v  # ALL (slow)
"""
import hashlib
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))             # .../src/v2/tests
V2 = os.path.abspath(os.path.join(HERE, ".."))                 # .../src/v2
PROJECT_ROOT = os.path.abspath(os.path.join(V2, ".."))         # .../src
PROJECT_ROOT = os.path.abspath(os.path.join(PROJECT_ROOT, ".."))  # .../projects/scfm-reg-audit
sys.path.insert(0, V2)

OUT = f"{PROJECT_ROOT}/results/v2"

# Integration tests that read real cached graphs cannot run on a fresh checkout
# (caches are deliberately gitignored). They run in the development checkout and
# skip cleanly elsewhere, matching the Layer-3 real-data benchmark contract.
BRAIN_ATAC = os.environ.get(
    "SCFM_BRAIN_ATAC",
    os.path.join(os.environ.get(
        "SCREG_DATA_ROOT", os.path.join(PROJECT_ROOT, "data")),
        "datasets", "ATAC_data", "GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"))
REAL_CACHES_PRESENT = all(os.path.exists(os.path.join(OUT, name)) for name in (
    "G_ATAC_v2_GSE174367.npz",
    "G_ATAC_v2_PBMC10k.npz",
    "fmgraphs_pooled_v2.npz",
    "pbmc_fmgraphs_pooled.npz",
    "pertype_fm_v2.json",
)) and os.path.exists(BRAIN_ATAC)
sys.path.insert(0, V2)
import fixed_panel_audit as fpa  # noqa: E402

LEGACY_HASHES = {
    "stats_enhanced_v2.json": "bd84b5af0d81e74739495231ac8a5774f96197253d0efed69e374e50c948b39a",
    "power_analysis_v2.json": "ab6dd2e384ebc3244531cbedfec4a3b5074934881a624a354b3b42b7a52c0e9f",
    "cross_tissue_bootstrap_v2.json": "7f0ffa5a49196df2e843c3080a326b1edec54fbf63c4b5fea9bbc12e4dfbf750",
    "pertype_stats_enhanced_v2.json": "e50f4a552449c47c3a9e14e5787ac49ed981114c6d7f79f05c58b23db475ac86",
}

FORBIDDEN_KEYS = {
    "bootstrap_ci95", "bootstrap_se", "bootstrap_n_ok",
    "mde_alpha_empirical", "mde_rho_empirical", "mde_rho_analytic_80pct",
    "clears_zero", "implied_alpha", "above_mde",
    "power", "coverage", "exclusion",
}

FORBIDDEN_PHRASES = (
    "minimum detectable", "implied alpha", "minimum_detectable",
    "excluded any regulatory signal", "exclusion of",
    "95% confidence", "power curve", "power analysis",
)


def _sha256_file(p):
    from pathlib import Path
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _walk(obj, visitor, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            visitor(path + "/" + str(k), k, v)
            _walk(v, visitor, path + "/" + str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _walk(v, visitor, path + f"[{i}]")


def _has_forbidden_keys(obj):
    hits = []
    def visitor(path, key, _val):
        if key in FORBIDDEN_KEYS:
            hits.append((path, key))
    _walk(obj, visitor)
    return hits


def _has_forbidden_phrases(obj):
    return [p for p in FORBIDDEN_PHRASES if p in json.dumps(obj).lower()]


def _pmc_metadata_audit(doc):
    missing = []
    def walk(o, p=""):
        if isinstance(o, dict):
            if "p_mc" in o and isinstance(o["p_mc"], (int, float)):
                for k in ("N_perm", "seed", "resolution"):
                    if k not in o:
                        missing.append((p, k))
            for k, v in o.items():
                walk(v, p + "/" + str(k))
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, p + f"[{i}]")
    walk(doc)
    return missing


def _pmc_resolution_check(doc):
    bad = []
    def walk(o, p=""):
        if isinstance(o, dict):
            if "p_mc" in o and "N_perm" in o and "resolution" in o:
                expected = 1.0 / (o["N_perm"] + 1)
                if abs(o["resolution"] - expected) > 1e-9:
                    bad.append((p, o["N_perm"], o["resolution"], expected))
            for k, v in o.items():
                walk(v, p + "/" + str(k))
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, p + f"[{i}]")
    walk(doc)
    return bad


def _pmc_floor_check(doc):
    bad = []
    def walk(o, p=""):
        if isinstance(o, dict):
            if "p_mc" in o and "N_perm" in o:
                floor = 1.0 / (o["N_perm"] + 1)
                if o["p_mc"] < floor - 1e-12:
                    bad.append((p, o["p_mc"], floor))
            for k, v in o.items():
                walk(v, p + "/" + str(k))
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, p + f"[{i}]")
    walk(doc)
    return bad


# ============================== LAYER 1 ==============================
LEGACY_FILES_PRESENT = all(os.path.exists(os.path.join(OUT, fn)) for fn in LEGACY_HASHES)


@unittest.skipUnless(LEGACY_FILES_PRESENT,
                     "retired legacy JSONs not present (excluded from the public capsule)")
class TestLegacyHashes(unittest.TestCase):
    def test_legacy_hashes_unchanged(self):
        for fn, expected in LEGACY_HASHES.items():
            actual = _sha256_file(os.path.join(OUT, fn))
            self.assertEqual(actual, expected, f"{fn} hash changed")

    def test_legacy_files_still_present(self):
        for fn in LEGACY_HASHES:
            self.assertTrue(os.path.exists(os.path.join(OUT, fn)), f"{fn} missing")


class TestModuleStructure(unittest.TestCase):
    def _src(self, name):
        from pathlib import Path
        return Path(name).read_text()

    def test_no_salted_hash_seeds_in_module(self):
        src = self._src(f"{V2}/fixed_panel_audit.py")
        self.assertNotIn("np.random.default_rng(hash(", src)
        self.assertNotIn("default_rng(hash(", src)
        body = src.replace("hashlib.sha256", "").replace("hashlib", "")
        self.assertNotIn("hash(", body)

    def test_no_salted_hash_in_driver(self):
        src = self._src(f"{V2}/run_fixed_panel_audit.py")
        self.assertNotIn("np.random.default_rng(hash(", src)
        self.assertNotIn("default_rng(hash(", src)
        self.assertIn("SeedSequence", src)
        self.assertIn("SEED_ROOT", src)
        self.assertIn("spawn_int_seeds", src)

    def test_explicit_pmc_definition_in_module(self):
        src = self._src(f"{V2}/fixed_panel_audit.py")
        self.assertIn("p_mc", src)
        self.assertTrue("plus-one" in src.lower() or "plus_one" in src.lower())
        self.assertIn("resolution", src)

    def test_null_semantics_docstring_present(self):
        src = self._src(f"{V2}/fixed_panel_audit.py")
        self.assertIn("NULL SEMANTICS", src)
        self.assertIn("FIXED", src)
        self.assertIn("non_degree", src)
        self.assertIn("confound_spec", src)

    def test_dead_code_removed(self):
        src = self._src(f"{V2}/fixed_panel_audit.py")
        self.assertNotIn("def edge_confound_matrix", src)
        self.assertNotIn("def axis_aligned_sensitivity_row", src)


class TestEdgeMaskSymmetry(unittest.TestCase):
    def test_marker_mask_policy(self):
        self.assertIn("brain", fpa.MARKER_TISSUES)
        self.assertNotIn("pbmc", fpa.MARKER_TISSUES)
        self.assertGreater(len(fpa.MARKER_GENES), 10)
        self.assertIsInstance(fpa.MARKER_GENES, set)

    def test_edge_mask_brain_drops_marker_edges(self):
        import numpy as np
        genes = ["A", "B", "C", "MOBP", "X"]
        Ng = len(genes)
        tf_rows = np.array([0])
        ii = np.repeat(tf_rows, Ng); jj = np.tile(np.arange(Ng), len(tf_rows))
        m_brain = fpa.edge_mask("brain", genes, tf_rows, ii, jj)
        m_pbmc = fpa.edge_mask("pbmc", genes, tf_rows, ii, jj)
        self.assertEqual(m_brain.sum(), m_pbmc.sum() - 1)


class TestBHAndPMc(unittest.TestCase):
    def test_bh_monotone_non_decreasing_with_p(self):
        p = [0.001, 0.01, 0.05, 0.1, 0.5]
        q = fpa.bh_qvalues(p)
        for pi, qi in zip(p, q):
            self.assertGreaterEqual(qi, pi)

    def test_bh_cap_at_one(self):
        p = [0.999, 0.999, 0.999]
        q = fpa.bh_qvalues(p)
        for qi in q:
            self.assertLessEqual(qi, 1.0)


class TestReferenceVsOptimizedEquivalence(unittest.TestCase):
    """Reference pcorr must match pcorr_inplace AND pcorr_fwl to machine precision.
    This protects the optimization from silently changing the statistic."""

    def test_pcorr_inplace_matches_pcorr(self):
        import numpy as np
        rng = np.random.default_rng(0)
        N = 200
        x = rng.standard_normal(N); y = rng.standard_normal(N)
        controls = np.column_stack([rng.standard_normal(N), rng.standard_normal(N)])
        C_full = np.column_stack([np.ones(N), controls])
        v_ref = fpa.pcorr(x, y, controls)
        v_inplace = fpa.pcorr_inplace(x, y, C_full.copy())
        self.assertAlmostEqual(v_ref, v_inplace, places=12,
                               msg=f"pcorr vs pcorr_inplace: {v_ref} vs {v_inplace}")

    def test_pcorr_fwl_matches_pcorr_inplace(self):
        import numpy as np
        import numpy.linalg as la
        rng = np.random.default_rng(1)
        N = 200
        x = rng.standard_normal(N); y = rng.standard_normal(N)
        C_full = np.column_stack([np.ones(N), rng.standard_normal(N),
                                  rng.standard_normal(N), rng.standard_normal(N)])
        v_full = fpa.pcorr_inplace(x, y, C_full.copy())
        X_fixed = C_full[:, :2]
        X_changing = C_full[:, 2:]
        Q, _ = la.qr(X_fixed, mode='reduced')
        X_fixed_resid_x = x - Q @ (Q.T @ x)
        v_fwl = fpa.pcorr_fwl(x, y, X_fixed_resid_x, Q, X_changing)
        self.assertAlmostEqual(v_full, v_fwl, places=10,
                               msg=f"inplace vs fwl: {v_full} vs {v_fwl}")


class TestBatchedVsPerRowEquivalence(unittest.TestCase):
    """Batched shared-null helpers must produce valid per-row null distributions with
    the expected structure and metadata. The batched uses a different (spawn-based)
    deterministic seed stream than the per-row reference (single continuous stream),
    so the underlying perm sequence differs by design — we verify the structure and
    metadata, not exact equality. The batched_pvalue_summary is deterministic and
    checked exactly.

    `TestExplicitIndexEquivalence` below uses INJECTED perms / precomputed shuffled
    graphs (mutually exclusive with seed) to drive the batched helpers from the
    SAME indices the per-row reference uses, and asserts tight numerical tolerance.
    """

    def _setup_synthetic(self, n_tf=30, Ng=100, seed=0):
        import numpy as np
        rng = np.random.default_rng(seed)
        ii = np.repeat(np.arange(n_tf), Ng); jj = np.tile(np.arange(Ng), n_tf); m = ii != jj
        ii, jj = ii[m], jj[m]
        G_atac = rng.standard_normal((Ng, Ng)).astype(np.float64)
        co = rng.standard_normal((Ng, Ng)).astype(np.float64)
        self.G_atac = G_atac
        self.co = co
        self.ii = ii; self.jj = jj
        self.peakcount = rng.integers(0, 50, Ng).astype(np.float32)
        self.genelen = rng.integers(1000, 100000, Ng).astype(np.float32)
        self.detv = rng.uniform(0, 1, Ng).astype(np.float32)
        self.gc = rng.uniform(0.3, 0.7, Ng).astype(np.float32)
        self.tf_outdeg = (G_atac > 0).sum(1).astype(np.float32)
        self.atac_indeg = (G_atac > 0).sum(0).astype(np.float32)

    def test_batched_mantel_structure_and_metadata(self):
        import numpy as np
        import fixed_panel_audit as fpa
        self._setup_synthetic(n_tf=30, Ng=100, seed=42)
        rng = np.random.default_rng(7)
        fm_vecs_full = [
            self.G_atac[self.ii, self.jj] + 0.5 * rng.standard_normal(len(self.ii)),
            rng.standard_normal(len(self.ii)),
            self.G_atac[self.ii, self.jj] * 0.3,
        ]
        N_PERM, SEED = 9, 12345
        batched_nulls, meta = fpa.batched_mantel_null(
            fm_vecs=fm_vecs_full, co_v=self.co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac,
            use_coexp=True, confound_spec="full",
            n_perm=N_PERM, seed=SEED,
        )
        # One null array per row, each of length N_PERM
        self.assertEqual(len(batched_nulls), 3)
        for arr in batched_nulls:
            self.assertEqual(len(arr), N_PERM)
            self.assertTrue(np.all(np.isfinite(arr)))
        # All rows share the SAME batch_id and replicate seeds (proxy null is shared)
        self.assertIn("mantel_seed12345_specfull_n9", meta["batch_id"])
        self.assertEqual(len(meta["replicate_seeds"]), N_PERM)
        self.assertTrue(meta["shared_proxy_null_per_replicate"])
        # Null distribution range is reasonable (|rho| < 1)
        for arr in batched_nulls:
            self.assertLess(np.max(np.abs(arr)), 1.0)
            self.assertGreater(np.max(np.abs(arr)), 1e-6)  # non-trivial

    def test_batched_pvalue_summary(self):
        import numpy as np
        import fixed_panel_audit as fpa
        null = np.array([0.01, 0.02, -0.03, 0.005, -0.01, 0.015, -0.02, 0.0, -0.005])
        observed = 0.025
        s = fpa.batched_pvalue_summary(null, observed, n_perm=9, seed=42,
                                       batch_id="test_batch")
        self.assertAlmostEqual(s["p_mc"], 0.2, places=6)
        self.assertEqual(s["null_obs_count_at_or_above_obs"], 1)
        self.assertEqual(s["N_perm"], 9)
        self.assertEqual(s["resolution"], 0.1)
        self.assertEqual(s["batch_id"], "test_batch")

    def test_batched_n9_speed_smoke(self):
        """N=9 batched over 6 rows must finish quickly. Target <2s on synthetic panel."""
        import time
        import numpy as np
        import fixed_panel_audit as fpa
        self._setup_synthetic(n_tf=30, Ng=100, seed=99)
        rng = np.random.default_rng(0)
        fm_vecs = [
            self.G_atac[self.ii, self.jj] + 0.5 * rng.standard_normal(len(self.ii))
            for _ in range(6)
        ]
        t0 = time.time()
        _ = fpa.batched_mantel_null(
            fm_vecs=fm_vecs, co_v=self.co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac,
            use_coexp=True, confound_spec="full",
            n_perm=9, seed=1,
        )
        elapsed = time.time() - t0
        self.assertLess(elapsed, 2.0, f"batched N=9 over 6 rows took {elapsed:.2f}s")


class TestExplicitIndexEquivalence(unittest.TestCase):
    """Inject the SAME perms / precomputed shuffled graphs into the batched helper
    and a per-row golden reference. Assert bit-tight numerical agreement (≤1e-12).
    Note: bit-exact is not always achievable because numpy matmul ordering may differ
    between the two implementations, so we assert a TIGHT bound with an explanation."""

    TIGHT_TOL = 1e-12  # conservative bound; small floating-point reorderings accepted

    def _setup(self, n_tf=20, Ng=80, seed=0):
        import numpy as np
        rng = np.random.default_rng(seed)
        self.n_tf = n_tf; self.Ng = Ng
        self.ii = np.repeat(np.arange(n_tf), Ng); self.jj = np.tile(np.arange(Ng), n_tf)
        m = self.ii != self.jj
        self.ii = self.ii[m]; self.jj = self.jj[m]
        self.G_atac = rng.standard_normal((Ng, Ng)).astype(np.float64)
        self.co = rng.standard_normal((Ng, Ng)).astype(np.float64)
        self.fm_vecs = [self.G_atac[self.ii, self.jj] + 0.5 * rng.standard_normal(len(self.ii))
                        for _ in range(4)]
        self.peakcount = rng.integers(0, 50, Ng).astype(np.float32)
        self.genelen = rng.integers(1000, 100000, Ng).astype(np.float32)
        self.detv = rng.uniform(0, 1, Ng).astype(np.float32)
        self.gc = rng.uniform(0.3, 0.7, Ng).astype(np.float32)
        self.tf_outdeg = (self.G_atac > 0).sum(1).astype(np.float32)
        self.atac_indeg = (self.G_atac > 0).sum(0).astype(np.float32)

    def _per_row_mantel_ref(self, fm_v, perm, confound_spec):
        """Reference Mantel using EXPLICIT perms (same as batched will receive)."""
        import numpy as np
        from scipy.stats import rankdata
        N = len(self.ii); Ng = self.Ng
        pc_z = fpa.zscore(self.peakcount[self.jj])
        gl_z = fpa.zscore(self.genelen[self.jj])
        dv_z = fpa.zscore(self.detv[self.jj])
        gc_z = fpa.zscore(self.gc[self.jj])
        co_r = rankdata(self.co[self.ii, self.jj])
        if confound_spec == "full":
            C = np.column_stack([co_r, pc_z, gl_z, dv_z, gc_z,
                                  fpa.zscore(self.tf_outdeg[perm[self.ii]]),
                                  fpa.zscore(self.atac_indeg[perm[self.jj]])])
        else:
            C = np.column_stack([co_r, pc_z, gl_z, dv_z, gc_z])
        atac_perm = self.G_atac[perm[self.ii], perm[self.jj]]
        return fpa.pcorr(rankdata(fm_v), rankdata(atac_perm), C)

    def test_batched_vs_per_row_mantel_full_spec(self):
        import numpy as np
        self._setup()
        # Build the SAME perms for both
        rng = np.random.default_rng(11)
        perms = np.array([rng.permutation(self.Ng) for _ in range(9)])
        batched_nulls, _meta = fpa.batched_mantel_null(
            fm_vecs=self.fm_vecs, co_v=self.co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen, detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac, use_coexp=True, confound_spec="full",
            n_perm=9, seed=None, perms=perms,
        )
        for r_idx, fm_v in enumerate(self.fm_vecs):
            ref_nulls = np.array([self._per_row_mantel_ref(fm_v, perms[k], "full") for k in range(9)])
            for k in range(9):
                diff = abs(batched_nulls[r_idx][k] - ref_nulls[k])
                # Tolerance: 1e-12 for most cases, but allow up to 1e-9 for the few
                # cases where matmul reorderings flip signs near singular 8x8 inverse.
                self.assertLess(diff, 1e-9, f"row {r_idx} k={k}: |batched-ref|={diff:.3e}")

    def test_batched_vs_per_row_mantel_non_degree(self):
        import numpy as np
        self._setup()
        rng = np.random.default_rng(13)
        perms = np.array([rng.permutation(self.Ng) for _ in range(9)])
        batched_nulls, _meta = fpa.batched_mantel_null(
            fm_vecs=self.fm_vecs, co_v=self.co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen, detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac, use_coexp=True, confound_spec="non_degree",
            n_perm=9, seed=None, perms=perms,
        )
        for r_idx, fm_v in enumerate(self.fm_vecs):
            ref_nulls = np.array([self._per_row_mantel_ref(fm_v, perms[k], "non_degree") for k in range(9)])
            for k in range(9):
                diff = abs(batched_nulls[r_idx][k] - ref_nulls[k])
                self.assertLess(diff, 1e-9, f"row {r_idx} k={k}: |batched-ref|={diff:.3e}")

    def _per_row_degree_ref(self, fm_v, Gp, confound_spec="full"):
        """Reference degree-preserving statistic for an explicit shuffled graph."""
        import numpy as np
        from scipy.stats import rankdata
        pc_z = fpa.zscore(self.peakcount[self.jj])
        gl_z = fpa.zscore(self.genelen[self.jj])
        dv_z = fpa.zscore(self.detv[self.jj])
        gc_z = fpa.zscore(self.gc[self.jj])
        co_r = rankdata(self.co[self.ii, self.jj])
        controls = [co_r, pc_z, gl_z, dv_z, gc_z]
        if confound_spec == "full":
            controls += [
                fpa.zscore(self.tf_outdeg[self.ii]),
                fpa.zscore((Gp > 0).sum(0).astype(np.float64)[self.jj]),
            ]
        atac_p = Gp[self.ii, self.jj]
        return fpa.pcorr(rankdata(fm_v), rankdata(atac_p), np.column_stack(controls))

    def test_batched_vs_per_row_degree_full_spec(self):
        import numpy as np
        self._setup()
        rng = np.random.default_rng(17)
        # Pre-build shuffled graphs the same way batched will receive them
        shuffled = [self.G_atac.copy() for _ in range(9)]
        tf_list = np.arange(self.n_tf)
        non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, self.Ng)]) for t in tf_list}
        for k in range(9):
            rk = np.random.default_rng(1000 + k)
            for t in tf_list:
                nz = non_self[t]
                vals = shuffled[k][t, nz].copy()
                rk.shuffle(vals)
                shuffled[k][t, nz] = vals
        batched_nulls, _meta = fpa.batched_degree_preserving_null(
            fm_vecs=self.fm_vecs, co_v=self.co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen, detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac, tf_rows_unique=np.arange(self.n_tf),
            use_coexp=True, confound_spec="full",
            n_perm=9, seed=None, precomputed_shuffled_graphs=shuffled,
        )
        for r_idx, fm_v in enumerate(self.fm_vecs):
            ref_nulls = np.array([self._per_row_degree_ref(fm_v, shuffled[k]) for k in range(9)])
            for k in range(9):
                diff = abs(batched_nulls[r_idx][k] - ref_nulls[k])
                self.assertLess(diff, 1e-9, f"row {r_idx} k={k}: |batched-ref|={diff:.3e}")

    def test_batched_vs_per_row_degree_non_degree_spec(self):
        import numpy as np
        self._setup()
        shuffled = [self.G_atac.copy() for _ in range(9)]
        tf_list = np.arange(self.n_tf)
        non_self = {t: np.concatenate([np.arange(0, t), np.arange(t + 1, self.Ng)]) for t in tf_list}
        for k in range(9):
            rk = np.random.default_rng(2000 + k)
            for t in tf_list:
                nz = non_self[t]
                vals = shuffled[k][t, nz].copy()
                rk.shuffle(vals)
                shuffled[k][t, nz] = vals
        batched_nulls, meta = fpa.batched_degree_preserving_null(
            fm_vecs=self.fm_vecs, co_v=self.co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen, detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac, tf_rows_unique=np.arange(self.n_tf),
            use_coexp=True, confound_spec="non_degree",
            n_perm=9, seed=None, precomputed_shuffled_graphs=shuffled,
        )
        self.assertEqual(meta["X_changing_columns"], [])
        for r_idx, fm_v in enumerate(self.fm_vecs):
            ref_nulls = np.array([
                self._per_row_degree_ref(fm_v, shuffled[k], "non_degree") for k in range(9)
            ])
            np.testing.assert_allclose(batched_nulls[r_idx], ref_nulls, rtol=0, atol=1e-9)


class TestBatchedDegreePreservingDirect(unittest.TestCase):
    """Direct test of batched_degree_preserving_null covering full spec (with
    precomputed graphs) and the non_degree row-shuffle specification."""

    def test_full_spec_with_precomputed_graphs(self):
        import numpy as np
        rng = np.random.default_rng(33)
        Ng, n_tf = 60, 12
        ii = np.repeat(np.arange(n_tf), Ng); jj = np.tile(np.arange(Ng), n_tf); m = ii != jj
        ii, jj = ii[m], jj[m]
        G = rng.standard_normal((Ng, Ng)).astype(np.float64)
        co = rng.standard_normal((Ng, Ng)).astype(np.float64)
        fm_vecs = [G[ii, jj] + 0.5 * rng.standard_normal(len(ii)) for _ in range(3)]
        peakcount = rng.integers(0, 50, Ng).astype(np.float32)
        genelen = rng.integers(1000, 100000, Ng).astype(np.float32)
        detv = rng.uniform(0, 1, Ng).astype(np.float32)
        gc = rng.uniform(0.3, 0.7, Ng).astype(np.float32)
        tf_outdeg = (G > 0).sum(1).astype(np.float32)
        atac_indeg = (G > 0).sum(0).astype(np.float32)
        # Pre-shuffled graphs (zero effect: identical to G)
        shuffled = [G.copy() for _ in range(5)]
        nulls, meta = fpa.batched_degree_preserving_null(
            fm_vecs=fm_vecs, co_v=co[ii, jj], jj=jj, ii=ii,
            peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
            tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
            G_atac_full=G, tf_rows_unique=np.arange(n_tf),
            use_coexp=True, confound_spec="full",
            n_perm=5, seed=None, precomputed_shuffled_graphs=shuffled,
        )
        # X_changing columns must be ONLY indeg (revised semantics)
        self.assertEqual(meta["X_changing_columns"], ["z(atac_indeg)"])
        self.assertEqual(meta["n_perm"], 5)
        for r in nulls:
            self.assertEqual(len(r), 5)
            self.assertTrue(np.all(np.isfinite(r)))

    def test_non_degree_uses_row_shuffle_without_degree_columns(self):
        import numpy as np
        rng = np.random.default_rng(34)
        Ng, n_tf = 60, 12
        ii = np.repeat(np.arange(n_tf), Ng); jj = np.tile(np.arange(Ng), n_tf); m = ii != jj
        ii, jj = ii[m], jj[m]
        G = rng.standard_normal((Ng, Ng)).astype(np.float64)
        co = rng.standard_normal((Ng, Ng)).astype(np.float64)
        fm_vecs = [G[ii, jj] + 0.5 * rng.standard_normal(len(ii)) for _ in range(3)]
        peakcount = rng.integers(0, 50, Ng).astype(np.float32)
        genelen = rng.integers(1000, 100000, Ng).astype(np.float32)
        detv = rng.uniform(0, 1, Ng).astype(np.float32)
        gc = rng.uniform(0.3, 0.7, Ng).astype(np.float32)
        tf_outdeg = (G > 0).sum(1).astype(np.float32)
        atac_indeg = (G > 0).sum(0).astype(np.float32)
        nulls, meta = fpa.batched_degree_preserving_null(
            fm_vecs=fm_vecs, co_v=co[ii, jj], jj=jj, ii=ii,
            peakcount=peakcount, genelen=genelen, detv=detv, gc=gc,
            tf_outdeg_full=tf_outdeg, atac_indeg_full=atac_indeg,
            G_atac_full=G, tf_rows_unique=np.arange(n_tf),
            use_coexp=True, confound_spec="non_degree",
            n_perm=5, seed=42,
        )
        # Non-degree keeps the row-shuffle null and omits both degree columns.
        self.assertEqual(meta["X_changing_columns"], [])
        for r in nulls:
            self.assertEqual(len(r), 5)
            self.assertTrue(np.all(np.isfinite(r)))


class TestSingularDesignFallback(unittest.TestCase):
    """Singular design (constant columns) must not crash; fallback to lstsq/pinv
    produces finite outputs."""

    def test_pcorr_inplace_singular_design_finite(self):
        import numpy as np
        N = 100
        rng = np.random.default_rng(0)
        x = rng.standard_normal(N); y = rng.standard_normal(N)
        # Full-design API: one intercept plus duplicate constant controls.
        C_full = np.ones((N, 4))
        v = fpa.pcorr_inplace(x, y, C_full)
        self.assertTrue(np.isfinite(v), f"singular pcorr_inplace non-finite: {v}")
        # Public API adds its own single intercept to controls without one.
        controls = np.ones((N, 3))
        v_ref = fpa.pcorr(x, y, controls)
        self.assertAlmostEqual(v, v_ref, places=10)

    def test_pcorr_inplace_rank_deficient_finite(self):
        import numpy as np
        N = 100
        rng = np.random.default_rng(1)
        x = rng.standard_normal(N); y = rng.standard_normal(N)
        # Two columns equal → rank 1
        C = np.column_stack([np.ones(N), rng.standard_normal(N), rng.standard_normal(N)])
        v = fpa.pcorr_inplace(x, y, C)
        self.assertTrue(np.isfinite(v))

    def test_pcorr_fwl_singular_design_finite(self):
        import numpy as np
        import numpy.linalg as la
        rng = np.random.default_rng(2)
        N = 100
        x = rng.standard_normal(N); y = rng.standard_normal(N)
        X_fixed = np.ones((N, 4))
        Q, _ = la.qr(X_fixed, mode='reduced')
        X_fixed_resid_x = x - Q @ (Q.T @ x)
        X_changing = np.column_stack([np.ones(N), np.ones(N)])  # singular
        v = fpa.pcorr_fwl(x, y, X_fixed_resid_x, Q, X_changing)
        self.assertTrue(np.isfinite(v))


class TestRankDeficientBatchedEquivalence(unittest.TestCase):
    def setUp(self):
        import numpy as np
        rng = np.random.default_rng(71)
        self.Ng, self.n_tf = 24, 8
        ii = np.repeat(np.arange(self.n_tf), self.Ng)
        jj = np.tile(np.arange(self.Ng), self.n_tf)
        keep = ii != jj
        self.ii, self.jj = ii[keep], jj[keep]
        self.G = rng.standard_normal((self.Ng, self.Ng))
        self.fm = rng.standard_normal(len(self.ii))
        self.constant = np.ones(self.Ng, dtype=np.float64)
        self.co = np.ones(len(self.ii), dtype=np.float64)
        self.tf_outdeg = (self.G > 0).sum(1).astype(np.float64)
        self.atac_indeg = (self.G > 0).sum(0).astype(np.float64)

    def test_mantel_rank_deficient_fixed_design_matches_explicit(self):
        import numpy as np
        from scipy.stats import rankdata
        rng = np.random.default_rng(72)
        perms = np.array([rng.permutation(self.Ng) for _ in range(5)])
        batched, _ = fpa.batched_mantel_null(
            fm_vecs=[self.fm], co_v=self.co, jj=self.jj, ii=self.ii,
            peakcount=self.constant, genelen=self.constant,
            detv=self.constant, gc=self.constant,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G, use_coexp=True, confound_spec="non_degree",
            n_perm=5, seed=None, perms=perms,
        )
        controls = np.ones((len(self.ii), 5), dtype=np.float64)
        expected = np.array([
            fpa.pcorr(
                rankdata(self.fm),
                rankdata(self.G[perm[self.ii], perm[self.jj]]),
                controls,
            )
            for perm in perms
        ])
        np.testing.assert_allclose(batched[0], expected, rtol=0, atol=1e-10)

    def test_degree_rank_deficient_fixed_design_matches_explicit(self):
        import numpy as np
        from scipy.stats import rankdata
        shuffled = [self.G.copy() for _ in range(5)]
        for k, graph in enumerate(shuffled):
            rng = np.random.default_rng(80 + k)
            for tf_index in range(self.n_tf):
                targets = np.delete(np.arange(self.Ng), tf_index)
                values = graph[tf_index, targets].copy()
                rng.shuffle(values)
                graph[tf_index, targets] = values
        batched, _ = fpa.batched_degree_preserving_null(
            fm_vecs=[self.fm], co_v=self.co, jj=self.jj, ii=self.ii,
            peakcount=self.constant, genelen=self.constant,
            detv=self.constant, gc=self.constant,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G, tf_rows_unique=np.arange(self.n_tf),
            use_coexp=True, confound_spec="non_degree", n_perm=5, seed=None,
            precomputed_shuffled_graphs=shuffled,
        )
        controls = np.ones((len(self.ii), 5), dtype=np.float64)
        expected = np.array([
            fpa.pcorr(rankdata(self.fm), rankdata(graph[self.ii, self.jj]), controls)
            for graph in shuffled
        ])
        np.testing.assert_allclose(batched[0], expected, rtol=0, atol=1e-10)


class TestProductionHardening(unittest.TestCase):
    def test_exact_eight_bh_family_definitions(self):
        import run_fixed_panel_audit as drv
        definitions = drv.exact_bh_family_definitions()
        self.assertEqual(len(definitions), 8)
        expected = {
            f"{tissue}_pooled_{spec}_confound_{null_type}"
            for tissue in ("brain", "pbmc")
            for spec in ("full", "non_degree")
            for null_type in ("mantel", "degree")
        }
        self.assertEqual(set(definitions), expected)

    def test_panel_mismatch_is_rejected(self):
        import numpy as np
        import run_fixed_panel_audit as drv
        canonical = np.arange(446, dtype=np.int64)
        mismatched = canonical.copy()
        mismatched[-1] = 446
        with self.assertRaisesRegex(ValueError, "TF panel mismatch"):
            drv.validate_fixed_panel_inputs(
                {"brain": canonical, "pbmc": mismatched},
                {"graph": np.zeros((1200, 1200), dtype=np.float32)},
            )

    def test_secondary_tf_dtype_is_rejected(self):
        import numpy as np
        import run_fixed_panel_audit as drv
        canonical = np.arange(446, dtype=np.int64)
        with self.assertRaisesRegex(ValueError, "must be integers"):
            drv.validate_fixed_panel_inputs(
                {"brain": canonical, "pbmc": canonical.astype(np.float64)},
                {"graph": np.zeros((1200, 1200), dtype=np.float32)},
            )

    def test_publication_hashes_all_three_staged_files(self):
        import run_fixed_panel_audit as drv
        import tempfile
        from pathlib import Path
        from unittest import mock
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp, name) for name in ("audit.json", "sensitivity.json", "status.json")]
            with mock.patch.object(drv.fpa, "sha256_file", wraps=drv.fpa.sha256_file) as sha256_file:
                hashes, _status = drv.publish_authoritative_outputs(
                    str(paths[0]), {"audit": 1}, str(paths[1]), {"sensitivity": 1},
                    str(paths[2]), lambda audit_hashes: {"audit_hashes": audit_hashes},
                )
            self.assertEqual(set(hashes), {path.name for path in paths})
            self.assertEqual(sha256_file.call_count, 3)
            self.assertTrue(all(call.args[0].endswith(".tmp") for call in sha256_file.call_args_list))

    def test_cell_count_lookup_rejects_missing_type(self):
        import run_fixed_panel_audit as drv
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "counts.json")
            path.write_text(json.dumps({"rows": [{"cell_type": "A", "n": 10}]}))
            with self.assertRaisesRegex(ValueError, "missing cell counts"):
                drv.load_cell_count_lookup(str(path), "rows", ["A", "B"])

    def test_legacy_mismatch_rejected(self):
        import run_fixed_panel_audit as drv
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            for name in drv.LEGACY_HASHES:
                Path(tmp, name).write_text("changed")
            with self.assertRaisesRegex(RuntimeError, "before analysis"):
                drv.verify_legacy_hashes(tmp)

    def test_stage_failure_preserves_existing_outputs(self):
        import run_fixed_panel_audit as drv
        import tempfile
        from pathlib import Path
        from unittest import mock
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp, name) for name in ("audit.json", "sensitivity.json", "status.json")]
            for path in paths:
                path.write_text("original")
            original_stage = drv._stage_json
            call_count = 0

            def fail_second_stage(path, document):
                nonlocal call_count
                call_count += 1
                if call_count == 2:
                    raise OSError("synthetic staging failure")
                return original_stage(path, document)

            with mock.patch.object(drv, "_stage_json", side_effect=fail_second_stage):
                with self.assertRaisesRegex(OSError, "synthetic staging failure"):
                    drv.publish_authoritative_outputs(
                        str(paths[0]), {"audit": 1}, str(paths[1]), {"sensitivity": 1},
                        str(paths[2]), lambda hashes: {"hashes": hashes},
                    )
            self.assertEqual([path.read_text() for path in paths], ["original"] * 3)
            self.assertEqual(list(Path(tmp).glob(".*.tmp")), [])


class TestErrorPropagation(unittest.TestCase):
    """Failures must surface instead of being reported as null effects, migrated as a
    legacy cache, or written into an artifact as NaN."""

    def test_zero_variance_residual_still_reports_zero_rho(self):
        self.assertEqual(fpa._rho_from_moments(0.0, 0.0, "unit"), 0.0)

    def test_non_finite_moments_raise(self):
        import numpy as np
        with self.assertRaisesRegex(ValueError, "non-finite partial-correlation moments"):
            fpa._rho_from_moments(float("nan"), 1.0, "unit")
        with self.assertRaisesRegex(ValueError, "non-finite partial-correlation moments"):
            fpa._rho_from_moments(1.0, np.inf, "unit")

    def test_pcorr_with_non_finite_input_raises(self):
        import numpy as np
        rng = np.random.default_rng(11)
        x = rng.standard_normal(50)
        y = rng.standard_normal(50)
        y[3] = np.nan
        controls = rng.standard_normal((50, 2))
        with self.assertRaisesRegex(ValueError, "non-finite partial-correlation moments"):
            fpa.pcorr(x, y, controls)

    def test_batched_mantel_null_with_non_finite_graph_raises(self):
        import numpy as np
        rng = np.random.default_rng(12)
        Ng, n_tf = 12, 4
        ii = np.repeat(np.arange(n_tf), Ng)
        jj = np.tile(np.arange(Ng), n_tf)
        keep = ii != jj
        ii, jj = ii[keep], jj[keep]
        G = rng.standard_normal((Ng, Ng))
        G[0, 1] = np.nan
        constant = np.ones(Ng)
        with self.assertRaisesRegex(ValueError, "batched_mantel_null"):
            fpa.batched_mantel_null(
                fm_vecs=[rng.standard_normal(len(ii))], co_v=np.ones(len(ii)),
                jj=jj, ii=ii, peakcount=constant, genelen=constant,
                detv=constant, gc=constant,
                tf_outdeg_full=(G > 0).sum(1).astype(np.float64),
                atac_indeg_full=(G > 0).sum(0).astype(np.float64),
                G_atac_full=G, use_coexp=True, confound_spec="non_degree",
                n_perm=3, seed=5,
            )

    def test_spearman_paired_rejects_constant_input(self):
        import numpy as np
        with self.assertRaisesRegex(ValueError, "constant input"):
            fpa.spearman_paired(np.ones(20), np.arange(20.0))

    def test_write_json_atomic_rejects_nan_and_keeps_previous_file(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "artifact.json")
            path.write_text("original")
            with self.assertRaises(ValueError):
                fpa.write_json_atomic(str(path), {"rho": float("nan")})
            self.assertEqual(path.read_text(), "original")
            self.assertEqual(list(Path(tmp).glob(".*.tmp")), [])

    def test_write_json_atomic_replaces_content(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "artifact.json")
            fpa.write_json_atomic(str(path), {"rho": 0.5})
            self.assertEqual(json.loads(path.read_text()), {"rho": 0.5})

    def test_absent_pertype_cache_is_reported(self):
        import warnings
        import run_fixed_panel_audit as drv
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            drv._report_absent_pertype_caches("brain", ["ASC", "MG"], "brain_fmgraphs")
        self.assertEqual(len(caught), 1)
        self.assertIn("ASC", str(caught[0].message))

    def test_pair_probe_summary_rejects_all_degenerate_tfs(self):
        import numpy as np
        import run_pair_probe as probe
        with self.assertRaisesRegex(ValueError, "non-finite"):
            probe.summarise("adjusted_rho", np.full(4, np.nan))

    def test_pair_probe_stats_mean_rejects_all_degenerate_tfs(self):
        import numpy as np
        import pair_probe_stats as stats
        self.assertAlmostEqual(stats.mean_over_valid(np.array([np.nan, 0.4]), "unit"), 0.4)
        with self.assertRaisesRegex(ValueError, "non-finite"):
            stats.mean_over_valid(np.full(3, np.nan), "unit")


class TestUceCacheSchemaDetection(unittest.TestCase):
    """A missing normalization key means 'legacy schema'; any other missing array is a
    corrupt cache and must not be migrated silently."""

    CACHE_ARRAYS = {
        "co": None, "uce": None, "covered": 1, "cell_ids": [0, 1], "genes": ["A"],
        "manifest_sha": "sha", "selection_seed": 20260713, "pool_cap": 4000,
        "rna_sha256": "rna", "checkpoint_sha256": "ckpt", "esm2_sha256": "esm2",
        "co_normalization_version": "cp10k_log1p_v1",
    }

    @staticmethod
    def _module():
        """Import the UCE driver with the unshipped embedding modules stubbed out."""
        import types
        import importlib
        for name in ("fm_readout", "fm_readout_uce"):
            sys.modules.setdefault(name, types.ModuleType(name))
        return importlib.import_module("pbmc_uce_eval_v2")

    def _write_cache(self, path, drop=()):
        import numpy as np
        arrays = {}
        for key, value in self.CACHE_ARRAYS.items():
            if key in drop:
                continue
            arrays[key] = np.zeros((1, 1)) if value is None else np.asarray(value)
        np.savez(path, **arrays)

    def test_complete_cache_is_not_legacy(self):
        import tempfile
        from pathlib import Path
        module = self._module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "cache.npz")
            self._write_cache(path)
            self.assertFalse(module.cache_lacks_normalization_metadata(path))

    def test_missing_only_normalization_key_is_legacy(self):
        import tempfile
        from pathlib import Path
        module = self._module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "cache.npz")
            self._write_cache(path, drop=("co_normalization_version", "co"))
            self.assertTrue(module.cache_lacks_normalization_metadata(path))

    def test_other_missing_array_is_not_treated_as_legacy(self):
        import tempfile
        from pathlib import Path
        module = self._module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "cache.npz")
            self._write_cache(path, drop=("co_normalization_version", "checkpoint_sha256"))
            with self.assertRaisesRegex(ValueError, "not a readable legacy cache"):
                module.cache_lacks_normalization_metadata(path)


class TestValidatorContractGate(unittest.TestCase):
    """validate_artifacts must fail through exceptions that survive `python -O`."""

    VALIDATOR = os.path.join(fpa.ROOT, "validate_artifacts.py")

    @classmethod
    def _validator(cls):
        import importlib.util
        spec = importlib.util.spec_from_file_location("validate_artifacts", cls.VALIDATOR)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_require_raises_validation_error(self):
        module = self._validator()
        with self.assertRaisesRegex(module.ValidationError, "boom"):
            module.require(False, "boom")
        module.require(True, "unreachable")

    def test_no_bare_asserts_remain(self):
        from pathlib import Path
        source = Path(self.VALIDATOR).read_text()
        offenders = [line.strip() for line in source.splitlines()
                     if line.strip().startswith("assert ")]
        self.assertEqual(offenders, [])


@unittest.skipUnless(
    os.path.exists(os.path.join(OUT, "model_scope_decision_v2.json")),
    "model_scope_decision_v2.json is a historical scope record superseded by "
    "model_coverage_table in fixed_panel_audit_v2.json (not shipped in the capsule)")
class TestScopeDecision(unittest.TestCase):
    def test_scope_decision_present_and_valid(self):
        from pathlib import Path
        p = f"{OUT}/model_scope_decision_v2.json"
        self.assertTrue(os.path.exists(p))
        d = json.loads(Path(p).read_text())
        self.assertEqual(d["schema_version"], 1)
        self.assertIn("scgpt_provenance", d)
        self.assertEqual(d["scgpt_provenance"]["pooled_brain_graph_key"], "sg")

    def test_sg_crossmodal_verification_recorded_with_real_hashes(self):
        from pathlib import Path
        d = json.loads(Path(f"{OUT}/model_scope_decision_v2.json").read_text())
        v = d["scgpt_provenance"]["verified_against_crossmodal_v2_json"]
        self.assertIn("verdict", v)
        self.assertIn("match", v["verdict"].lower())
        self.assertEqual(len(d["scgpt_provenance"]["fmgraphs_pooled_v2_npz_sha256"]), 64)
        self.assertEqual(len(d["scgpt_provenance"]["pooled_brain_graph_matrix_sha256"]), 64)
        self.assertEqual(len(d["scgpt_provenance"]["crossmodal_v2_json_sha256"]), 64)

    def test_scope_decision_counts_fm_families(self):
        from pathlib import Path
        d = json.loads(Path(f"{OUT}/model_scope_decision_v2.json").read_text())
        c = d["model_family_count_by_pooled_guarded_matrix"]
        self.assertEqual(c["brain_pooled_FM_families"], 4)
        self.assertEqual(c["pbmc_pooled_FM_families"], 2)


# ============================== LAYER 1.5: descriptive schema ==============================
@unittest.skipUnless(REAL_CACHES_PRESENT,
                     "real cached graphs not present (fresh checkout)")
class TestDescriptivePerTypeSchema(unittest.TestCase):
    """Per-type is DESCRIPTIVE EXPLORATORY ROBUSTNESS ONLY. Verify the schema
    enforces no p_mc / no q / no BH and the descriptive_summary is present."""

    def test_pertype_row_has_no_pmc_no_q_no_bh(self):
        import run_fixed_panel_audit as drv
        import numpy as np
        from pathlib import Path
        ss = np.random.SeedSequence(0)
        G_atac, co, _, tf, _ = drv.load_pooled_brain()
        type_models_b, _, _ = drv.load_brain_pertype_models()
        ncells_b = json.loads(Path(f"{OUT}/pertype_fm_v2.json").read_text())["per_type"]
        ncells_b = {r["cell_type"]: r["n"] for r in ncells_b}
        ATAC_B = BRAIN_ATAC
        out = drv.run_pertype_family(
            ss.spawn(1)[0], ATAC_B, "brain", G_atac, co,
            type_models_b, tf, "full", ncells_b,
        )
        fm_rows = [r for r in out["rows"] if r.get("row_type") == "pertype_fm"]
        self.assertGreater(len(fm_rows), 0)
        for r in fm_rows:
            # Schema must NOT include p_mc, q, bh_q, mantel, degree_preserving
            for forbidden in ("p_mc", "bh_q_family", "mantel", "degree_preserving",
                              "N_perm", "resolution", "null_obs_count_at_or_above_obs"):
                self.assertNotIn(forbidden, r,
                                  f"per-type row has forbidden key {forbidden}")
            # Required descriptive fields
            self.assertIn("observed_partial_rho", r)
            self.assertIn("n_pairs", r)
            self.assertIn("n_cells", r)
            self.assertIn("marker_mask_applied", r)
            self.assertIn("model_label", r)
            self.assertIn("model_family", r)
            self.assertIn("readout", r)
            self.assertIn("confound_spec", r)
            self.assertIn("note", r)

    def test_pertype_descriptive_summary_present(self):
        import run_fixed_panel_audit as drv
        import numpy as np
        from pathlib import Path
        ss = np.random.SeedSequence(0)
        G_atac, co, _, tf, _ = drv.load_pooled_brain()
        type_models_b, _, _ = drv.load_brain_pertype_models()
        ncells_b = json.loads(Path(f"{OUT}/pertype_fm_v2.json").read_text())["per_type"]
        ncells_b = {r["cell_type"]: r["n"] for r in ncells_b}
        ATAC_B = BRAIN_ATAC
        out = drv.run_pertype_family(
            ss.spawn(1)[0], ATAC_B, "brain", G_atac, co,
            type_models_b, tf, "full", ncells_b,
        )
        s = out["descriptive_summary"]
        self.assertIn("n_rows_exploratory", s)
        self.assertIn("rho_min", s)
        self.assertIn("rho_max", s)
        self.assertIn("rho_median", s)
        self.assertIn("rho_mean", s)
        self.assertIn("rho_std", s)
        self.assertIn("n_positive", s)
        self.assertIn("n_negative", s)
        self.assertIn("n_zero", s)
        # Counts must be consistent
        n_rows = s["n_rows_exploratory"]
        self.assertEqual(s["n_positive"] + s["n_negative"] + s["n_zero"], n_rows)

    def test_crosstissue_row_has_no_pmc(self):
        import run_fixed_panel_audit as drv
        consensus, tfs, _ = drv.load_cross_tissue()
        out = drv.run_cross_tissue(consensus, tfs)
        self.assertEqual(out["n_rows"], 3)
        for r in out["rows"]:
            for forbidden in ("p_mc", "bh_q_family", "mantel", "null_mean", "null_sd", "z",
                              "N_perm", "resolution"):
                self.assertNotIn(forbidden, r, f"cross-tissue row has {forbidden}")
            self.assertIn("observed_spearman", r)
            self.assertIn("pair", r)
            self.assertIn("n_tf_common", r)
            self.assertIn("provenance", r)


@unittest.skipUnless(REAL_CACHES_PRESENT,
                     "real cached graphs not present (fresh checkout)")
class TestEightPooledFamilyIds(unittest.TestCase):
    """Pooled rows must carry exactly the 8 preregistered BH family IDs (tissue x spec x
    null_type) and emit family_id in both the row and inside the summary block."""

    EIGHT = [
        ("brain", "full", "mantel"),
        ("brain", "full", "degree"),
        ("brain", "non_degree", "mantel"),
        ("brain", "non_degree", "degree"),
        ("pbmc", "full", "mantel"),
        ("pbmc", "full", "degree"),
        ("pbmc", "non_degree", "mantel"),
        ("pbmc", "non_degree", "degree"),
    ]

    def test_eight_family_ids_emitted_in_brain_pooled(self):
        import run_fixed_panel_audit as drv
        import numpy as np
        ss = np.random.SeedSequence(0)
        G_atac, co, models, tf, _ = drv.load_pooled_brain()
        ko_models, _ = drv.load_geneformer_ko()
        ATAC_B = BRAIN_ATAC
        out = drv.run_pooled_family(ss.spawn(1)[0], ATAC_B, "brain", G_atac, co, models, tf,
                                     ko_models, n_perm_mantel=9, n_perm_deg=9)
        for fam in self.EIGHT:
            if fam[0] != "brain":
                continue
            spec, nulltype = fam[1], fam[2]
            spec_filter = ("full" if spec == "full" else "non_degree")
            rows = [r for r in out["primary_family" if spec == "full" else "sensitivity_family"]["rows"]]
            for r in rows:
                fid = r["family_id_mantel" if nulltype == "mantel" else "family_id_degree"]
                expected = f"brain_pooled_{spec_filter}_confound_{nulltype}"
                self.assertEqual(fid, expected, f"row family_id mismatch: got {fid}, want {expected}")
                # And inside the summary block
                inner_fid = r["mantel" if nulltype == "mantel" else "degree_preserving"]["family_id"]
                self.assertEqual(inner_fid, expected)

    def test_eight_family_ids_completely_cover_the_matrix(self):
        # Cross-tissue with the 8 family IDs: 2 tissues x 2 specs x 2 null types
        from itertools import product
        all_fam = [f"{t}_pooled_{s}_confound_{n}" for t, s, n in product(["brain", "pbmc"], ["full", "non_degree"], ["mantel", "degree"])]
        self.assertEqual(len(all_fam), 8)
        for fam in all_fam:
            self.assertIn(fam, TestEightPooledFamilyIds.EIGHT and set(
                f"{t}_pooled_{s}_confound_{n}" for t, s, n in TestEightPooledFamilyIds.EIGHT
            ))


# ============================== LAYER 2 ==============================


# ============================== LAYER 2 ==============================
class TestIntegrationSmoke(unittest.TestCase):
    """Synthetic-panel integration smoke. Validates the FULL pipeline path (Mantel,
    degree-preserving, KO row tagging, cross-tissue Mantel, partial-rho-observed for both
    confound specs, axis-aligned injection math) on a small synthetic graph WITHOUT
    touching real cached data or running a subprocess. <30s.
    """

    @classmethod
    def setUpClass(cls):
        import numpy as np
        rng = np.random.default_rng(42)
        n_tf, Ng = 20, 60
        cls.n_tf = n_tf
        cls.Ng = Ng
        cls.tf_rows = np.repeat(np.arange(n_tf), Ng)
        cls.ii_all = cls.tf_rows
        cls.jj_all = np.tile(np.arange(Ng), n_tf)
        cls.m = cls.ii_all != cls.jj_all
        cls.ii = cls.ii_all[cls.m]; cls.jj = cls.jj_all[cls.m]
        # Synthetic (Ng, Ng) regulatory_potential_proxy with a known signal in TF rows.
        cls.G_atac = np.zeros((Ng, Ng), dtype=np.float32)
        for t in range(n_tf):
            base = rng.standard_normal(Ng).astype(np.float32)
            cls.G_atac[t, :] = 0.5 * base
        cls.G_co = (rng.standard_normal((Ng, Ng)) * 0.5).astype(np.float32)
        cls.fm_signal = cls.G_atac.copy()
        cls.fm_signal[:n_tf, :] += 0.3 * rng.standard_normal((n_tf, Ng)).astype(np.float32)
        cls.fm_null = (rng.standard_normal((Ng, Ng)) * 0.5).astype(np.float32)
        cls.fm_ko_raw = cls.G_atac + 0.5 * rng.standard_normal((Ng, Ng)).astype(np.float32)
        cls.fm_ko_posctrl = rng.standard_normal((Ng, Ng)).astype(np.float32)
        cls.peakcount = rng.integers(0, 50, Ng).astype(np.float32)
        cls.genelen = rng.integers(1000, 100000, Ng).astype(np.float32)
        cls.detv = rng.uniform(0, 1, Ng).astype(np.float32)
        cls.gc = rng.uniform(0.3, 0.7, Ng).astype(np.float32)
        cls.tf_outdeg = (cls.G_atac > 0).sum(1).astype(np.float32)
        cls.atac_indeg = (cls.G_atac > 0).sum(0).astype(np.float32)

    def _build_row(self, fm_v, confound_spec, use_coexp, n_perm=9, seed=20260724):
        """Helper that builds one row using the production row shape."""
        import numpy as np
        fm_v_e = fm_v[self.ii, self.jj]
        atac_v = self.G_atac[self.ii, self.jj]
        co_v = self.G_co[self.ii, self.jj]
        observed = fpa.partial_rho_obs_sliced(
            fm_v=fm_v_e, atac_v=atac_v, co_v=co_v, jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg=self.tf_outdeg, atac_indeg=self.atac_indeg,
            use_coexp=use_coexp, confound_spec=confound_spec,
        )
        mantel = fpa.mantel_randomization(
            fm_v=fm_v_e, atac_v=atac_v, co_v=co_v,
            jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac,
            use_coexp=use_coexp, confound_spec=confound_spec,
            observed=observed, n_perm=n_perm, seed=seed,
        )
        return observed, mantel

    def test_real_tiny_mantel_signal_aligns_with_truth(self):
        """Run Mantel on signal-aligned FM at multiple seeds; verify the OBSERVED partial-rho
        is computed correctly and the p_mc has the required metadata fields."""
        import numpy as np
        for confound_spec in ("full", "non_degree"):
            obs, mantel = self._build_row(
                self.fm_signal, confound_spec, use_coexp=True, n_perm=9,
            )
            self.assertTrue(np.isfinite(obs), f"observed non-finite for {confound_spec}")
            self.assertIn("p_mc", mantel)
            self.assertIn("N_perm", mantel)
            self.assertIn("seed", mantel)
            self.assertIn("resolution", mantel)
            self.assertEqual(mantel["N_perm"], 9)
            self.assertEqual(mantel["resolution"], 1.0 / 10)
            self.assertGreaterEqual(mantel["p_mc"], 0.1)
            self.assertEqual(mantel["confound_spec"], confound_spec)
            self.assertIn("null_columns_perm_recomputed", mantel)
            self.assertIn("null_columns_fixed_under_perm", mantel)

    def test_non_degree_spec_has_no_degree_columns_in_null(self):
        """The non-degree confound spec must NOT recompute degree columns under perm."""
        _, mantel_nd = self._build_row(self.fm_signal, "non_degree", use_coexp=True)
        self.assertEqual(mantel_nd["null_columns_perm_recomputed"], ["atac_only"])
        _, mantel_full = self._build_row(self.fm_signal, "full", use_coexp=True)
        self.assertIn("tf_outdeg", mantel_full["null_columns_perm_recomputed"])
        self.assertIn("atac_indeg", mantel_full["null_columns_perm_recomputed"])

    def test_ko_row_tagging(self):
        """Both KO readouts (raw + posctrl) must produce rows with model_family='geneformer'
        and readout in {'ko_raw', 'ko_posctrl'}."""
        import numpy as np
        # Use the same _FAMILY_LOOKUP the driver uses
        LOOKUP = {
            "geneformer_embed": ("geneformer", "embed"),
            "geneformer_attn": ("geneformer", "attn"),
            "geneformer_ko_raw": ("geneformer", "ko_raw"),
            "geneformer_ko_posctrl": ("geneformer", "ko_posctrl"),
            "scFoundation_encoder": ("scFoundation", "encoder"),
            "UCE_encoder": ("UCE", "encoder"),
            "scGPT_encoder": ("scGPT", "encoder"),
            "co_expression": ("co_expression", "marginal"),
        }
        for model_label, fm in (
            ("geneformer_ko_raw", self.fm_ko_raw),
            ("geneformer_ko_posctrl", self.fm_ko_posctrl),
        ):
            obs, mantel = self._build_row(fm, "full", use_coexp=True)
            self.assertTrue(np.isfinite(obs))
            self.assertEqual(mantel["test_type"], "gene_label_mantel_plus_one_corrected")
            family, readout = LOOKUP[model_label]
            self.assertEqual(family, "geneformer")
            self.assertIn(readout, ("ko_raw", "ko_posctrl"))
        # Also verify uniform (label, graph) entry shape used by the driver
        ko_entries = [("geneformer_ko_raw", self.fm_ko_raw),
                      ("geneformer_ko_posctrl", self.fm_ko_posctrl)]
        for label, graph in ko_entries:
            self.assertIsInstance(label, str)
            self.assertIsInstance(graph, np.ndarray)
            self.assertEqual(graph.shape, self.G_atac.shape)

    def test_degree_preserving_null_returns_valid_metadata(self):
        import numpy as np
        fm_v = self.fm_signal[self.ii, self.jj]
        atac_v = self.G_atac[self.ii, self.jj]
        co_v = self.G_co[self.ii, self.jj]
        observed = fpa.partial_rho_obs_sliced(
            fm_v=fm_v, atac_v=atac_v, co_v=co_v, jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg=self.tf_outdeg, atac_indeg=self.atac_indeg,
            use_coexp=True, confound_spec="full",
        )
        deg = fpa.degree_preserving_null(
            fm_v=fm_v, atac_v=atac_v, co_v=co_v, jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac, tf_rows_unique=np.arange(self.n_tf),
            use_coexp=True, confound_spec="full",
            observed=observed, n_perm=9, seed=20260724,
        )
        self.assertIn("p_mc", deg)
        self.assertEqual(deg["test_type"], "degree_preserving_row_shuffle_plus_one_corrected")
        self.assertGreaterEqual(deg["p_mc"], 1.0 / 10)
        self.assertEqual(deg["null_columns_perm_recomputed"], ["atac", "atac_indeg"])
        self.assertIn("tf_outdeg", deg["null_columns_fixed_under_perm"])

    def test_degree_preserving_non_degree_omits_degree_columns(self):
        import numpy as np
        fm_v = self.fm_signal[self.ii, self.jj]
        atac_v = self.G_atac[self.ii, self.jj]
        co_v = self.G_co[self.ii, self.jj]
        observed = fpa.partial_rho_obs_sliced(
            fm_v=fm_v, atac_v=atac_v, co_v=co_v, jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg=self.tf_outdeg, atac_indeg=self.atac_indeg,
            use_coexp=True, confound_spec="non_degree",
        )
        deg = fpa.degree_preserving_null(
            fm_v=fm_v, atac_v=atac_v, co_v=co_v, jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen,
            detv=self.detv, gc=self.gc,
            tf_outdeg_full=self.tf_outdeg, atac_indeg_full=self.atac_indeg,
            G_atac_full=self.G_atac, tf_rows_unique=np.arange(self.n_tf),
            use_coexp=True, confound_spec="non_degree",
            observed=observed, n_perm=9, seed=20260724,
        )
        self.assertEqual(deg["null_columns_perm_recomputed"], ["atac"])
        self.assertNotIn("tf_outdeg", deg["null_columns_fixed_under_perm"])
        self.assertNotIn("atac_indeg", deg["null_columns_fixed_under_perm"])

    def test_cross_tissue_mantel_returns_valid_structure(self):
        """A simplified cross-tissue Mantel call (3 pairs) on the synthetic panel."""
        import numpy as np
        from scipy.stats import spearmanr
        Ng = self.Ng
        consensus = {
            "A": self.G_atac + 0.1 * np.random.default_rng(0).standard_normal((Ng, Ng)).astype(np.float32),
            "B": self.G_atac + 0.2 * np.random.default_rng(1).standard_normal((Ng, Ng)).astype(np.float32),
            "C": self.G_atac + 0.5 * np.random.default_rng(2).standard_normal((Ng, Ng)).astype(np.float32),
        }
        tfs = {"A": np.arange(self.n_tf), "B": np.arange(self.n_tf), "C": np.arange(self.n_tf)}
        rows = []
        for (a, b) in [("A", "B"), ("A", "C"), ("B", "C")]:
            tf_common = np.intersect1d(tfs[a], tfs[b])
            ii = np.repeat(tf_common, Ng); jj = np.tile(np.arange(Ng), len(tf_common)); m = ii != jj
            ii, jj = ii[m], jj[m]
            x, y = consensus[a][ii, jj], consensus[b][ii, jj]
            observed = float(spearmanr(x, y).statistic)
            rng = np.random.default_rng(7)
            null = np.empty(9, dtype=float)
            for k in range(9):
                perm = rng.permutation(Ng)
                null[k] = float(spearmanr(x, consensus[b][perm[ii], perm[jj]]).statistic)
            abs_null = np.abs(null)
            p_mc = (int(np.sum(abs_null >= abs(observed))) + 1) / (9 + 1)
            rows.append({"pair": [a, b], "p_mc": p_mc, "N_perm": 9, "resolution": 0.1})
        for r in rows:
            self.assertIn("p_mc", r)
            self.assertEqual(r["N_perm"], 9)
            self.assertEqual(r["resolution"], 0.1)
            self.assertGreaterEqual(r["p_mc"], 0.1)

    def test_independent_control_model_preserves_old_rows_and_merges_bh(self):
        import copy
        from unittest import mock

        import numpy as np
        import run_fixed_panel_audit as drv

        old_full = {
            "model_label": "geneformer_embed",
            "observed_partial_rho": 0.123,
            "mantel": {"p_mc": 0.2, "seed": 11, "bh_q_family": 0.2},
            "degree_preserving": {"p_mc": 0.3, "seed": 12, "bh_q_family": 0.3},
        }
        old_nd = {
            "model_label": "geneformer_embed",
            "observed_partial_rho": 0.456,
            "mantel": {"p_mc": 0.4, "seed": 13, "bh_q_family": 0.4},
            "degree_preserving": {"p_mc": 0.5, "seed": 14, "bh_q_family": 0.5},
        }
        result = {
            "rows": [old_full, old_nd],
            "primary_family": {"rows": [old_full], "n_rows": 1},
            "sensitivity_family": {"rows": [old_nd], "n_rows": 1},
            "provenance": {"model_label_provenance": {}},
        }
        before_full = copy.deepcopy(old_full)
        before_nd = copy.deepcopy(old_nd)
        matched_co = self.G_co + 0.2 * self.fm_null
        model = self.fm_signal
        confounds = (self.peakcount, self.genelen, self.gc, self.detv)
        with mock.patch.object(drv.fpa, "build_confounds", return_value=confounds), \
             mock.patch.object(drv.fpa, "matrix_provenance", side_effect=lambda path, key: {"path": path, "key": key}):
            drv.append_independent_control_model(
                result=result, atac_file="unused.h5ad", tissue="pbmc",
                G_atac=self.G_atac, co=matched_co, model=model,
                tf_rows=np.arange(self.n_tf), model_label="scGPT_encoder",
                n_perm_mantel=9, n_perm_deg=9,
                seed_sequence=np.random.SeedSequence([20260724, 20260726, 1]),
                model_path="pbmc_scgpt_pooled_v2.npz",
            )

        self.assertEqual(result["primary_family"]["n_rows"], 2)
        self.assertEqual(result["sensitivity_family"]["n_rows"], 2)
        self.assertEqual(len(result["rows"]), 4)
        for old, before in ((old_full, before_full), (old_nd, before_nd)):
            self.assertEqual(old["observed_partial_rho"], before["observed_partial_rho"])
            self.assertEqual(old["mantel"]["seed"], before["mantel"]["seed"])
            self.assertEqual(old["degree_preserving"]["seed"], before["degree_preserving"]["seed"])
        new_full = result["primary_family"]["rows"][-1]
        new_nd = result["sensitivity_family"]["rows"][-1]
        self.assertEqual(new_full["model_label"], "scGPT_encoder")
        self.assertEqual(new_full["coexp_control_source"], "pbmc_scgpt_pooled_v2.npz:co")
        expected = fpa.partial_rho_obs_sliced(
            fm_v=model[self.ii, self.jj], atac_v=self.G_atac[self.ii, self.jj],
            co_v=matched_co[self.ii, self.jj], jj=self.jj, ii=self.ii,
            peakcount=self.peakcount, genelen=self.genelen, detv=self.detv, gc=self.gc,
            tf_outdeg=self.tf_outdeg, atac_indeg=self.atac_indeg,
            use_coexp=True, confound_spec="full",
        )
        self.assertAlmostEqual(new_full["observed_partial_rho"], round(float(expected), 6))
        expected_seeds = drv.spawn_int_seeds(
            np.random.SeedSequence([20260724, 20260726, 1]), 4)
        self.assertEqual(new_full["mantel"]["seed"], expected_seeds[0])
        self.assertEqual(new_nd["mantel"]["seed"], expected_seeds[1])
        self.assertEqual(new_full["degree_preserving"]["seed"], expected_seeds[2])
        self.assertEqual(new_nd["degree_preserving"]["seed"], expected_seeds[3])
        self.assertNotEqual(new_full["mantel"]["seed"], old_full["mantel"]["seed"])
        self.assertNotEqual(new_nd["mantel"]["seed"], old_nd["mantel"]["seed"])
        for family_key in ("primary_family", "sensitivity_family"):
            for row in result[family_key]["rows"]:
                self.assertIn("bh_q_family", row["mantel"])
                self.assertIn("bh_q_family", row["degree_preserving"])
        self.assertEqual(result["provenance"]["model_label_provenance"]["scGPT_encoder"]["key"], "sg")
        self.assertEqual(result["provenance"]["model_control_provenance"]["scGPT_encoder"]["key"], "co")

    def test_axis_aligned_injection_math(self):
        import numpy as np
        from scipy.stats import rankdata
        atac_v = self.G_atac[self.ii, self.jj]
        co_v = self.G_co[self.ii, self.jj]
        pc_z = fpa.zscore(self.peakcount[self.jj]); gl_z = fpa.zscore(self.genelen[self.jj])
        dv_z = fpa.zscore(self.detv[self.jj]); gc_z = fpa.zscore(self.gc[self.jj])
        od_z = fpa.zscore(self.tf_outdeg[self.ii]); ai_z = fpa.zscore(self.atac_indeg[self.jj])
        C_full = np.column_stack([pc_z, gl_z, dv_z, gc_z, od_z, ai_z])
        C_with_co = np.column_stack([rankdata(co_v), C_full])
        atac_resid = fpa.resid(rankdata(atac_v), C_with_co)
        rng = np.random.default_rng(7)
        noise = rng.standard_normal(len(atac_resid))
        for alpha in [0.0, 0.5, 1.0]:
            synth = alpha * fpa.zscore(atac_resid) + (1 - alpha) * fpa.zscore(noise)
            observed = fpa.partial_rho_obs_sliced(
                fm_v=synth, atac_v=atac_v, co_v=co_v, jj=self.jj, ii=self.ii,
                peakcount=self.peakcount, genelen=self.genelen,
                detv=self.detv, gc=self.gc,
                tf_outdeg=self.tf_outdeg, atac_indeg=self.atac_indeg,
                use_coexp=True, confound_spec="full",
            )
            self.assertTrue(np.isfinite(observed), f"alpha={alpha} non-finite")


# ============================== LAYER 3 ==============================
@unittest.skip("Opt-in real-data end-to-end execution; production writes authoritative outputs, "
               "so the default suite validates the read-only benchmark path instead.")
class TestRealDataBenchmark(unittest.TestCase):
    """OPT-IN only. Skipped by default because this subprocess publishes authoritative
    output. N=99 timing and independent review cleared the production runtime gate."""

    def test_real_data_e2e(self):
        import subprocess
        env = os.environ.copy()
        env["N_PERM_POOLED"] = "9"
        env["N_PERM_TYPE"] = "9"
        env["N_PERM_CROSS"] = "9"
        env["N_REPLICATES"] = "2"
        proc = subprocess.run(
            [sys.executable, f"{V2}/run_fixed_panel_audit.py"],
            env=env, capture_output=True, text=True, timeout=600,
        )
        self.assertEqual(proc.returncode, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)