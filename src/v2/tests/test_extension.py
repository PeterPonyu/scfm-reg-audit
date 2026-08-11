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


class TestBuilderOutResolution(unittest.TestCase):
    def test_default_out_when_unset(self):
        out = pathmod.resolve_builder_out_dir(
            extension_out=None,
            peerj_lock=False,
            env={},
            create=False,
        )
        self.assertEqual(out, pathmod.CANONICAL_V2_RESULTS.resolve())

    def test_extension_out_honored(self):
        # Use repo-relative overlay path (confinement allowlist).
        rel = f"results/v2/extension/_pytest_out_{uuid4().hex}"
        target = ROOT / rel
        try:
            out = pathmod.resolve_builder_out_dir(
                extension_out=rel,
                peerj_lock=True,
                create=True,
            )
            self.assertEqual(out, target.resolve())
            self.assertTrue(out.is_dir())
            self.assertTrue(pathmod.under_extension_overlay(out))
        finally:
            _wipe_dir(target)

    def test_peerj_lock_refuses_canonical_without_extension_out(self):
        with self.assertRaises(ValueError) as ctx:
            pathmod.resolve_builder_out_dir(
                extension_out=None,
                peerj_lock=True,
                env={"SCREG_PEERJ_SUPPORT_LOCK": "1"},
                create=False,
            )
        self.assertIn("SCREG_EXTENSION_OUT", str(ctx.exception))

    def test_peerj_lock_refuses_canonical_explicit_out(self):
        with self.assertRaises(ValueError):
            pathmod.resolve_builder_out_dir(
                extension_out="results/v2",
                peerj_lock=True,
                create=False,
            )

    def test_peak_name_normalize(self):
        self.assertEqual(
            pathmod.normalize_peak_name("chr1-9776-10668"), "chr1:9776-10668"
        )
        self.assertTrue(pathmod.peak_name_is_builder_ready("chr1:9776-10668"))
        self.assertFalse(pathmod.peak_name_is_builder_ready("chr1-9776-10668"))


class TestDescartesBridge(unittest.TestCase):
    def test_absent_fail_closed(self):
        import descartes_bridge as db

        with tempfile.TemporaryDirectory() as td:
            status = db.bridge_status(pilot_dir=Path(td))
            self.assertEqual(status["status"], "absent_local")
            self.assertFalse(status["network_fetch_performed"])
            self.assertIn("place", status["message"].lower())

    def test_ready_h5ad(self):
        import anndata as ad
        import descartes_bridge as db
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp

        with tempfile.TemporaryDirectory() as td:
            pilot = Path(td)
            X = sp.csr_matrix(np.ones((4, 3), dtype=np.float32))
            A = ad.AnnData(
                X=X,
                obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]),
                var=pd.DataFrame(
                    index=["chr1:100-200", "chr1:300-400", "chr2:50-80"]
                ),
            )
            h5 = pilot / db.EXPECTED_H5AD_NAME
            A.write_h5ad(h5)
            status = db.bridge_status(pilot_dir=pilot)
            self.assertEqual(status["status"], "ready_atac_file")
            self.assertIn("build_command", status)
            self.assertIn("SCREG_EXTENSION_OUT", status["build_command"])
            self.assertIn("SCREG_PEERJ_SUPPORT_LOCK=1", status["build_command"])


class TestBmmcPrepareTiny(unittest.TestCase):
    def test_extract_renames_peaks(self):
        import anndata as ad
        import bmmc_prepare as bp
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            X = sp.csr_matrix(np.ones((5, 4), dtype=np.float32))
            var = pd.DataFrame(
                {
                    "feature_types": ["GEX", "ATAC", "ATAC", "GEX"],
                },
                index=["GENEA", "chr1-10-20", "chr2-30-40", "GENEB"],
            )
            obs = pd.DataFrame(
                {"cell_type": ["T", "T", "B", "B", "Mono"]},
                index=[f"bc{i}" for i in range(5)],
            )
            src = td_path / "multiome_tiny.h5ad"
            ad.AnnData(X=X, obs=obs, var=var).write_h5ad(src)
            out_atac = td_path / "peaks.h5ad"
            out_meta = td_path / "meta.csv.gz"
            written = bp.extract_atac_peak_matrix(
                src, out_atac=out_atac, out_meta=out_meta
            )
            self.assertEqual(written["n_peaks"], 2)
            B = ad.read_h5ad(out_atac)
            self.assertEqual(list(B.var_names), ["chr1:10-20", "chr2:30-40"])
            self.assertTrue(out_meta.exists())


class TestConstructBmmcDryRun(unittest.TestCase):
    def test_bmmc_next_steps_mention_prepare_or_build(self):
        import construct_hooks as ch

        plan = ch.run_construct("bmmc", execute=False)
        self.assertEqual(plan["status"], "awaiting_g_atac")
        joined = " ".join(plan["next_steps"])
        self.assertIn("P3", joined)
        # Either prepare path (raw multiome) or build path (prepared peaks).
        self.assertTrue(
            "prepare-bmmc" in joined or "build overlay G_ATAC" in joined,
            joined,
        )
        self.assertEqual(plan["env"]["SCREG_PEERJ_SUPPORT_LOCK"], "1")
        self.assertTrue(
            plan["env"]["SCREG_EXTENSION_OUT"].startswith("results/v2/extension/")
        )


if __name__ == "__main__":
    unittest.main()
