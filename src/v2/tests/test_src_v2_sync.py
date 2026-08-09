"""Guard: capsule COPY_FILES sources in src/v2 must match shipped src/ twins."""
import os
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
ROOT = V2.parent.parent

COPY_PAIRS = [
    ("src/v2/fixed_panel_audit.py", "src/fixed_panel_audit.py"),
    ("src/v2/run_fixed_panel_audit.py", "src/run_fixed_panel_audit.py"),
    ("src/v2/pbmc_uce_eval_v2.py", "src/pbmc_uce_eval_v2.py"),
    ("src/v2/brain_coexp_baseline_null.py", "src/brain_coexp_baseline_null.py"),
    ("src/v2/pbmc_coexp_baseline_null.py", "src/pbmc_coexp_baseline_null.py"),
    ("src/v2/pair_probe_stats.py", "src/pair_probe_stats.py"),
    ("src/v2/run_pair_probe.py", "src/run_pair_probe.py"),
    ("src/v2/pbmc_cache.py", "src/pbmc_cache.py"),
    ("src/v2/benchmark_n99.py", "src/benchmark_n99.py"),
    ("src/v2/tests/test_fixed_panel_audit.py", "src/tests/test_fixed_panel_audit.py"),
]


class TestSrcV2Sync(unittest.TestCase):
    def test_copy_files_pairs_are_byte_identical(self):
        for src_rel, dst_rel in COPY_PAIRS:
            with self.subTest(pair=f"{src_rel}->{dst_rel}"):
                left = (ROOT / src_rel).read_bytes()
                right = (ROOT / dst_rel).read_bytes()
                self.assertEqual(left, right)


if __name__ == "__main__":
    unittest.main()
