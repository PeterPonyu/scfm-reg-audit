"""Unit tests for the co-expression baseline drivers.

src/brain_coexp_baseline_null.py and src/pbmc_coexp_baseline_null.py score the
baseline at its own primary rung. Their only tissue-specific logic outside
main() is the promoter peak-count confound, so these tests drive
brain_peakcount / pbmc_peakcount against a synthetic gene-coordinate table and
a synthetic ATAC h5ad, and pin the seed contract both drivers publish.
"""
import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, SRC)
import brain_coexp_baseline_null as brain  # noqa: E402
import fixed_panel_audit as fpa  # noqa: E402
import pbmc_coexp_baseline_null as pbmc  # noqa: E402

try:
    import anndata as ad
    import scipy.sparse as sp
    ANNDATA_MISSING = None
except ImportError as exc:  # pragma: no cover - environment dependent
    ANNDATA_MISSING = str(exc)

GENES = ["PLUS", "MINUS", "OTHERCHR", "NOCOORD"]
# name -> (chrom, start, end, strand); PROM padding is applied by the drivers.
COORDS = {
    "PLUS": ("chr1", 10_000, 12_000, "+"),
    "MINUS": ("chr1", 50_000, 52_000, "-"),
    "OTHERCHR": ("chr2", 10_000, 12_000, "+"),
}
PEAKS = [
    "chr1:10500-10600",    # inside PLUS body
    "chr1:11900-12100",    # midpoint 12000, PLUS end boundary
    "chr1:8000-8100",      # PLUS promoter window only (upstream of start)
    "chr1:100-200",        # outside every window
    "chr1:52500-52600",    # MINUS downstream promoter window
    "chr1:60000-60100",    # beyond the MINUS window
    "chr2:10500-10600",    # inside OTHERCHR body
    "chr3:10500-10600",    # chromosome with no panel gene
]


def _write_coords(path, extra_lines=()):
    with open(path, "w") as fh:
        for name, (chrom, start, end, strand) in COORDS.items():
            fh.write(f"{chrom}\t{start}\t{end}\t{strand}\t{name}\n")
        for line in extra_lines:
            fh.write(line)


def _write_atac(path, peaks=PEAKS):
    matrix = sp.csr_matrix(np.ones((2, len(peaks)), dtype=np.float32))
    adata = ad.AnnData(X=matrix)
    adata.var_names = peaks
    adata.write_h5ad(path)


@unittest.skipIf(ANNDATA_MISSING, f"anndata unavailable: {ANNDATA_MISSING}")
class PeakCountMixin:
    """Shared contract: both drivers count promoter-window peaks per panel gene."""

    module = None
    function_name = None
    atac_attribute = None

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.coords = os.path.join(self.tmp.name, "gene_coords.tsv")
        self.atac = os.path.join(self.tmp.name, "atac.h5ad")
        _write_coords(self.coords)
        _write_atac(self.atac)

    def peakcount(self, genes=GENES, coords=None, atac=None):
        with mock.patch.object(fpa, "COORDS", coords or self.coords), \
                mock.patch.object(fpa, "load_manifest",
                                  return_value=(list(genes), {}, "sha")), \
                mock.patch.object(self.module, self.atac_attribute, atac or self.atac):
            return getattr(self.module, self.function_name)()

    def test_promoter_window_counts_per_gene(self):
        counts = self.peakcount()
        self.assertEqual(counts.dtype, np.float32)
        self.assertEqual(counts.shape, (len(GENES),))
        # PLUS: body peak, boundary peak at its end, and the upstream promoter peak.
        # MINUS: the downstream promoter peak only. OTHERCHR: its own chr2 peak.
        # NOCOORD has no annotation row, so it stays at zero.
        np.testing.assert_array_equal(counts, np.array([3, 1, 1, 0], dtype=np.float32))

    def test_promoter_padding_is_strand_aware(self):
        """The window is padded upstream on + strand and downstream on - strand."""
        self.assertEqual(self.module.PROM, fpa.PROM)

        def window(lo, hi):
            atac = os.path.join(self.tmp.name, f"peaks_{lo}_{hi}.h5ad")
            _write_atac(atac, peaks=[f"chr1:{lo}-{hi}"])
            return self.peakcount(atac=atac)

        plus_start, minus_end = COORDS["PLUS"][1], COORDS["MINUS"][2]
        inside_plus_promoter = window(plus_start - fpa.PROM + 100, plus_start - fpa.PROM + 300)
        outside_plus_promoter = window(plus_start - fpa.PROM - 300, plus_start - fpa.PROM - 100)
        inside_minus_promoter = window(minus_end + fpa.PROM - 300, minus_end + fpa.PROM - 100)
        outside_minus_promoter = window(minus_end + fpa.PROM + 100, minus_end + fpa.PROM + 300)

        self.assertEqual(inside_plus_promoter[0], 1.0)
        self.assertEqual(outside_plus_promoter[0], 0.0)
        self.assertEqual(inside_minus_promoter[1], 1.0)
        self.assertEqual(outside_minus_promoter[1], 0.0)
        # Padding never leaks onto the opposite side of the gene.
        self.assertEqual(inside_plus_promoter[1], 0.0)
        self.assertEqual(inside_minus_promoter[0], 0.0)

    def test_genes_absent_from_the_annotation_are_zero(self):
        counts = self.peakcount(genes=["NOCOORD", "PLUS"])
        self.assertEqual(counts[0], 0.0)
        self.assertGreater(counts[1], 0.0)

    def test_first_annotation_row_per_gene_wins(self):
        coords = os.path.join(self.tmp.name, "duplicated.tsv")
        _write_coords(coords, extra_lines=["chr2\t10000\t12000\t+\tPLUS\n"])
        np.testing.assert_array_equal(self.peakcount(coords=coords), self.peakcount())

    def test_annotation_rows_outside_the_panel_are_ignored(self):
        coords = os.path.join(self.tmp.name, "extra_gene.tsv")
        _write_coords(coords, extra_lines=["chr1\t10000\t12000\t+\tNOTINPANEL\n"])
        np.testing.assert_array_equal(self.peakcount(coords=coords), self.peakcount())

    def test_counts_are_independent_of_gene_order(self):
        reordered = list(reversed(GENES))
        counts = self.peakcount()
        reordered_counts = self.peakcount(genes=reordered)
        np.testing.assert_array_equal(reordered_counts, counts[::-1])


class TestBrainPeakCount(PeakCountMixin, unittest.TestCase):
    module = brain
    function_name = "brain_peakcount"
    atac_attribute = "ATAC_B"


class TestPbmcPeakCount(PeakCountMixin, unittest.TestCase):
    module = pbmc
    function_name = "pbmc_peakcount"
    atac_attribute = "ATAC_P"


class TestSeedContract(unittest.TestCase):
    def test_all_four_baseline_seeds_are_distinct_explicit_integers(self):
        seeds = [brain.BRAIN_MANTEL_SEED, brain.BRAIN_DEGREE_SEED,
                 pbmc.PBMC_MANTEL_SEED, pbmc.PBMC_DEGREE_SEED]
        for seed in seeds:
            self.assertIsInstance(seed, int)
        self.assertEqual(len(set(seeds)), 4)

    def test_published_seeds_match_the_shipped_baseline_json(self):
        import json
        from pathlib import Path
        capsule = Path(SRC).parent
        for name, mantel, degree in (
            ("brain_coexp_baseline_null_v2.public.json",
             brain.BRAIN_MANTEL_SEED, brain.BRAIN_DEGREE_SEED),
            ("pbmc_coexp_baseline_null_v2.public.json",
             pbmc.PBMC_MANTEL_SEED, pbmc.PBMC_DEGREE_SEED),
        ):
            path = capsule / "results" / name
            if not path.exists():  # pragma: no cover - development layout
                self.skipTest(f"{name} absent")
            doc = json.loads(path.read_text())
            with self.subTest(artifact=name):
                self.assertEqual(doc["seed_contract"], "explicit_integer_v1")
                self.assertEqual(doc["mantel_seed"], mantel)
                self.assertEqual(doc["degree_seed"], degree)
                self.assertEqual(doc["n_perm"], 999)

    def test_both_drivers_write_into_the_shared_results_directory(self):
        self.assertEqual(brain.OUT, fpa.OUT)
        self.assertEqual(pbmc.OUT, fpa.OUT)


if __name__ == "__main__":
    unittest.main()
