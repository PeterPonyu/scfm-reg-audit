"""Smoke tests for extension-lane local compute overlay (Option B-prime)."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path
from uuid import uuid4

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
ROOT = V2.parents[1]
EXT = V2 / "extension"
sys.path.insert(0, str(EXT))

import download_gate  # noqa: E402
import emit_claim_pack  # noqa: E402
import paths as pathmod  # noqa: E402
import registry as regmod  # noqa: E402


def _wipe_dir(path: Path) -> None:
    """Recursively remove ``path`` if it exists (files then dirs)."""
    if not path.exists():
        return
    for p in sorted(path.rglob("*"), reverse=True):
        if p.is_file():
            p.unlink()
        elif p.is_dir():
            p.rmdir()
    if path.is_dir():
        path.rmdir()


def _backup_download_env() -> dict[str, str | None]:
    return {
        k: os.environ.get(k)
        for k in (download_gate.APPROVAL_ENV, download_gate.PLAN_MATCH_ENV)
    }


def _restore_download_env(backup: dict[str, str | None]) -> None:
    for k, v in backup.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


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
        out = ROOT / "results" / "v2" / "extension" / f"_pytest_claim_pack_{uuid4().hex}"
        _wipe_dir(out)
        try:
            paths = emit_claim_pack.emit(out)
            dual = json.loads(paths["dual_null"].read_text())
            index = json.loads(paths["index"].read_text())
            self.assertEqual(dual["n_full_rows"], 13)
            self.assertEqual(dual["n_dual_null_full"], 7)
            self.assertEqual(index["counts"]["dual_null_full"], 7)
            self.assertEqual(index["counts"]["protocol_pass_frozen"], 0)
            self.assertFalse(index["peerj_support_rows_touched"])
        finally:
            _wipe_dir(out)


class TestPathConfinement(unittest.TestCase):
    def test_traversal_rejected(self):
        with self.assertRaises(ValueError):
            pathmod.assert_confined_write_path(
                "docs/reports/extension-claim-pack/../../secrets",
                label="--out-dir",
            )
        with self.assertRaises(ValueError):
            pathmod.assert_confined_write_path("/tmp/escape", label="--out-dir")
        with self.assertRaises(ValueError):
            pathmod.assert_safe_tag("../evil", label="tag")
        with self.assertRaises(ValueError):
            pathmod.assert_safe_tag("a/b", label="tag")

    def test_allowlisted_roots_ok(self):
        claim = pathmod.assert_confined_write_path(
            ROOT / "docs" / "reports" / "extension-claim-pack",
            label="--out-dir",
        )
        self.assertTrue(claim.is_absolute())
        ext = pathmod.assert_confined_write_path(
            "results/v2/extension/download-plans",
            label="out",
        )
        self.assertTrue(str(ext).endswith("download-plans"))

    def test_emit_claim_pack_rejects_escape(self):
        with self.assertRaises(ValueError):
            emit_claim_pack.emit(Path("/tmp/not-allowed-claim-pack"))


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
        self.assertEqual(proc.returncode, 2)
        self.assertIn("DISABLED", proc.stdout)


class TestDownloadGate(unittest.TestCase):
    def test_refuses_without_approval(self):
        env_backup = _backup_download_env()
        for k in (download_gate.APPROVAL_ENV, download_gate.PLAN_MATCH_ENV):
            os.environ.pop(k, None)
        try:
            code, payload = download_gate.run_download_gate("D1", write_plan=False)
            self.assertEqual(code, 2)
            self.assertEqual(payload["status"], "approval_required")
            self.assertFalse(payload["network_fetch_performed"])
        finally:
            _restore_download_env(env_backup)

    def test_approved_still_no_network(self):
        env_backup = _backup_download_env()
        os.environ[download_gate.APPROVAL_ENV] = "1"
        os.environ[download_gate.PLAN_MATCH_ENV] = "D1"
        try:
            with unittest.mock.patch("urllib.request.urlopen") as urlopen:
                with unittest.mock.patch("subprocess.run") as sprout:
                    code, payload = download_gate.run_download_gate(
                        "D1", write_plan=False
                    )
                    self.assertEqual(code, 0)
                    self.assertEqual(payload["status"], "approved_dry_run_no_fetch")
                    self.assertFalse(payload["network_fetch_performed"])
                    urlopen.assert_not_called()
                    sprout.assert_not_called()
        finally:
            _restore_download_env(env_backup)

    def test_approved_writes_dry_run_plan(self):
        env_backup = _backup_download_env()
        os.environ[download_gate.APPROVAL_ENV] = "1"
        os.environ[download_gate.PLAN_MATCH_ENV] = "D1"
        plan_path = (
            ROOT / "results" / "v2" / "extension" / "download-plans" / "D1.dry_run.json"
        )
        try:
            code, payload = download_gate.run_download_gate("D1", write_plan=True)
            self.assertEqual(code, 0)
            self.assertTrue(plan_path.exists())
            doc = json.loads(plan_path.read_text())
            self.assertFalse(doc["network_fetch_performed"])
            self.assertEqual(payload.get("dry_run_plan_path"), str(plan_path.relative_to(ROOT)))
        finally:
            plan_path.unlink(missing_ok=True)
            plans_dir = plan_path.parent
            if plans_dir.is_dir() and not any(plans_dir.iterdir()):
                plans_dir.rmdir()
            _restore_download_env(env_backup)

    def test_rejected_plan_fails_closed(self):
        env_backup = _backup_download_env()
        os.environ[download_gate.APPROVAL_ENV] = "1"
        os.environ[download_gate.PLAN_MATCH_ENV] = "D4"
        try:
            code, payload = download_gate.run_download_gate("D4", write_plan=False)
            self.assertEqual(code, 2)
            self.assertEqual(payload["status"], "rejected_plan")
            self.assertFalse(payload["network_fetch_performed"])
        finally:
            _restore_download_env(env_backup)

    def test_plan_mismatch_fails_closed(self):
        env_backup = _backup_download_env()
        os.environ[download_gate.APPROVAL_ENV] = "1"
        os.environ[download_gate.PLAN_MATCH_ENV] = "D2"
        try:
            with unittest.mock.patch("urllib.request.urlopen") as urlopen:
                code, payload = download_gate.run_download_gate("D1", write_plan=False)
                self.assertEqual(code, 2)
                self.assertEqual(payload["status"], "approval_required")
                self.assertFalse(payload["network_fetch_performed"])
                urlopen.assert_not_called()
        finally:
            _restore_download_env(env_backup)

    def test_invalid_plan_id_no_traceback(self):
        env_backup = _backup_download_env()
        try:
            code, payload = download_gate.run_download_gate("../evil", write_plan=False)
            self.assertEqual(code, 2)
            self.assertEqual(payload["status"], "invalid_plan_id")
            self.assertFalse(payload["network_fetch_performed"])
            # CLI path must also stay structured (no SystemExit traceback).
            rc = download_gate.main(["--plan-id", "a/b", "--no-write"])
            self.assertEqual(rc, 2)
        finally:
            _restore_download_env(env_backup)


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
        # Plan env stays PeerJ-redacted (no /home leak via ATAC_FILE).
        self.assertNotIn("/home/", plan["env"].get("ATAC_FILE", ""))

    def test_atac_file_machine_usable_when_local_exists(self):
        import construct_hooks as ch

        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as fh:
            tmp = Path(fh.name)
            fh.write(b"placeholder")
        try:
            env = ch.fibro_style_env("fibroblast", atac_file=str(tmp))
            self.assertEqual(Path(env["ATAC_FILE"]).resolve(), tmp.resolve())
            self.assertTrue(Path(env["ATAC_FILE"]).exists())
            self.assertEqual(env["ATAC_FILE_REDACTED"], tmp.name)
            safe = ch.log_safe_env(env)
            self.assertEqual(safe["ATAC_FILE"], env["ATAC_FILE_REDACTED"])
            self.assertNotEqual(safe["ATAC_FILE"], env["ATAC_FILE"])
        finally:
            tmp.unlink(missing_ok=True)


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
        out = ROOT / "results" / "v2" / "extension" / f"_pytest_cli_claim_pack_{uuid4().hex}"
        _wipe_dir(out)
        try:
            pack = subprocess.run(
                [sys.executable, str(cli), "claim-pack", "--out-dir", str(out)],
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(pack.returncode, 0, pack.stderr)
            self.assertTrue((out / "index.json").exists())
        finally:
            _wipe_dir(out)

    def test_cli_download_refuses_without_approval(self):
        cli = EXT / "cli.py"
        env = {k: v for k, v in os.environ.items() if k not in (
            download_gate.APPROVAL_ENV,
            download_gate.PLAN_MATCH_ENV,
        )}
        proc = subprocess.run(
            [sys.executable, str(cli), "download", "--plan-id", "D1", "--no-write"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        self.assertEqual(proc.returncode, 2)
        self.assertIn("approval", proc.stdout.lower() + proc.stderr.lower())


if __name__ == "__main__":
    unittest.main()
