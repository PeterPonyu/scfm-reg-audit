"""Smoke tests for extension-lane local compute overlay (Option B-prime)."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
ROOT = V2.parents[1]
EXT = V2 / "extension"
sys.path.insert(0, str(EXT))

import emit_claim_pack  # noqa: E402
import registry as regmod  # noqa: E402


class TestExtensionRegistry(unittest.TestCase):
    def test_rna_lakes_denied(self):
        reg = regmod.load_extension_registry()
        for tid in ("cancer_rna_lakes", "development_rna_lakes"):
            with self.assertRaises(PermissionError):
                reg.assert_may_emit_g_atac(tid)

    def test_bmmc_construct_candidate(self):
        reg = regmod.load_extension_registry()
        self.assertEqual(reg.get_tissue("bmmc")["role"], "construct_candidate")
        self.assertNotIn("bmmc", reg.summary()["extension_audit_candidates"])
        # construct lane still may emit G_ATAC (construct NPZ), but not FM Support
        reg.assert_may_emit_g_atac("bmmc")

    def test_methods_status_synced(self):
        reg = regmod.load_extension_registry()
        self.assertNotEqual(reg.get_method("degree_matched_random")["status"], "stub")
        self.assertNotEqual(reg.get_method("motif_only_rp")["status"], "stub")
        self.assertEqual(
            reg.get_method("collectri_prior")["status"], "skipped_no_local_cache"
        )


class TestClaimPackSmoke(unittest.TestCase):
    def test_emit_claim_pack_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            paths = emit_claim_pack.emit(out)
            dual = json.loads(paths["dual_null"].read_text())
            index = json.loads(paths["index"].read_text())
            self.assertEqual(dual["n_full_rows"], 13)
            self.assertEqual(dual["n_dual_null_full"], 7)
            self.assertEqual(index["counts"]["dual_null_full"], 7)
            self.assertEqual(index["counts"]["protocol_pass_frozen"], 0)
            self.assertFalse(index["peerj_support_rows_touched"])


class TestFetchDemoted(unittest.TestCase):
    def test_fetch_script_no_curl(self):
        script = EXT / "scripts" / "fetch_optional_pilots.sh"
        text = script.read_text()
        # Forbid executable downloader invocations (comments mentioning the word are OK).
        for line in text.splitlines():
            stripped = line.split("#", 1)[0].strip()
            if not stripped:
                continue
            self.assertNotRegex(stripped, r"\bcurl\b")
            self.assertNotRegex(stripped, r"\bwget\b")
        proc = subprocess.run(
            ["bash", str(script)],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(proc.returncode, 0)
        self.assertIn("DISABLED", proc.stdout)


class TestConstructFibroSmoke(unittest.TestCase):
    def test_fibro_execute_reads_locked_npz(self):
        locked = ROOT / "results" / "v2" / "G_ATAC_v2_GSE206767.npz"
        if not locked.exists():
            self.skipTest("locked fibro G_ATAC not present in this checkout")
        import construct_hooks as ch

        plan = ch.run_construct("fibroblast", execute=True)
        self.assertEqual(plan["status"], "executed")
        out = ROOT / "results" / "v2" / "extension" / "construct" / "GSE206767"
        self.assertTrue((out / "mantel_vs_locked.json").exists())
        self.assertTrue((out / "additive_decomp_row.json").exists())
        self.assertEqual(plan.get("peerj_support_rows_unchanged"), 13)


class TestBaselineEmitters(unittest.TestCase):
    def test_motif_only_and_degree_emit(self):
        locked = ROOT / "results" / "v2" / "G_ATAC_v2_GSE174367.npz"
        if not locked.exists():
            self.skipTest("locked brain G_ATAC not present in this checkout")
        import baseline_stubs as bs

        motif = bs.run_baseline("motif_only_rp", execute=True, proxy_tag="GSE174367")
        deg = bs.run_baseline("degree_matched_random", execute=True, proxy_tag="GSE174367")
        self.assertEqual(motif["run_status"], "emitted")
        self.assertEqual(deg["run_status"], "emitted")
        encode = bs.run_baseline("encode_chip_binding", execute=True)
        self.assertEqual(encode["run_status"], "emitted")


class TestCliEntrypoint(unittest.TestCase):
    def test_cli_registry_and_claim_pack(self):
        cli = EXT / "cli.py"
        reg = subprocess.run(
            [sys.executable, str(cli), "registry", "--json"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(reg.returncode, 0, reg.stderr)
        payload = json.loads(reg.stdout)
        self.assertIn("cancer_rna_lakes", payload["summary"]["forbidden_g_atac"])
        with tempfile.TemporaryDirectory() as tmp:
            pack = subprocess.run(
                [sys.executable, str(cli), "claim-pack", "--out-dir", tmp],
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(pack.returncode, 0, pack.stderr)
            self.assertTrue((Path(tmp) / "index.json").exists())


if __name__ == "__main__":
    unittest.main()
