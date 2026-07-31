import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, V2)

from audit_pbmc_remote_results import audit


class TestAuditPbmcRemoteResults(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.reference = {
            "n_pairs": 3,
            "n_tf": 2,
            "tissue": "PBMC",
            "types": ["T"],
            **{f"embed__m{i}": float(i) for i in range(6)},
            **{f"attn__m{i}": float(i) for i in range(6)},
        }
        self.reference_path = self.root / "reference.json"
        self.reference_path.write_text(json.dumps(self.reference))
        self.preflight_path = self.root / "pbmc_preflight_remote.json"
        self.preflight_path.write_text(json.dumps({
            "status": "passed",
            "manifest_sha": "manifest",
            "metrics": self.reference,
            "reference_metrics_sha256": hashlib.sha256(
                self.reference_path.read_bytes()).hexdigest(),
        }))
        self.eval_path = self.root / "pbmc_eval_v2.json"
        self.new_eval = dict(
            self.reference, **{f"scgpt__m{i}": float(i) for i in range(6)})
        self.eval_path.write_text(json.dumps(self.new_eval))
        self.graph_path = self.root / "pbmc_scgpt_pooled_v2.npz"
        np.savez(
            self.graph_path,
            co=np.eye(1200, dtype=np.float32),
            sg=np.eye(1200, dtype=np.float32),
            cell_ids=np.arange(4000),
            genes=np.asarray([f"g{i}" for i in range(1200)]),
            manifest_sha=np.asarray("manifest"),
            selection_seed=np.asarray(20260713),
            pool_cap=np.asarray(4000),
            rna_sha256=np.asarray("rna"),
        )

    def tearDown(self):
        self.tmp.cleanup()

    def create_evidence(self, execution_id="20260728-120000-123"):
        manifest = {"version": "1", "files": {}}
        manifest_sha = hashlib.sha256(
            json.dumps(manifest, sort_keys=True).encode()).hexdigest()
        manifest["manifest_sha256"] = manifest_sha
        run_id = manifest_sha[:16]
        manifest_path = self.root / f"manifest.{run_id}.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

        verification_path = self.root / (
            f"manifest_verification.{run_id}.{execution_id}.json")
        verification_path.write_text(json.dumps({
            "status": "verified",
            "manifest_sha256": manifest_sha,
            "verified_count": 1,
        }))
        preflight_log = self.root / f"preflight_run.{run_id}.{execution_id}.log"
        scgpt_log = self.root / f"scgpt_run.{run_id}.{execution_id}.log"
        preflight_log.write_text("preflight complete\n")
        scgpt_log.write_text("scGPT complete\n")

        receipt_name = f"run_receipt.{run_id}.{execution_id}.json"
        marker_fields = {
            "run_id": run_id,
            "execution_id": execution_id,
            "manifest_sha256": manifest_sha,
            "status": "success",
        }
        preflight_marker = self.root / f"PREFLIGHT_DONE.{run_id}.{execution_id}"
        preflight_marker.write_text("\n".join(
            f"{key}={value}" for key, value in {
                **marker_fields, "phase": "preflight"}.items()) + "\n")
        scgpt_marker = self.root / f"SCGPT_DONE.{run_id}.{execution_id}"
        scgpt_marker.write_text("\n".join(
            f"{key}={value}" for key, value in {
                **marker_fields, "receipt": receipt_name,
                "phase": "scgpt"}.items()) + "\n")

        targets = {
            "manifest_json": manifest_path,
            "manifest_verification_json": verification_path,
            "preflight_json": self.preflight_path,
            "eval_json": self.eval_path,
            "graph_npz": self.graph_path,
            "reference_json": self.reference_path,
            "preflight_log": preflight_log,
            "scgpt_log": scgpt_log,
        }
        receipt_path = self.root / receipt_name
        receipt_path.write_text(json.dumps({
            "run_id": run_id,
            "execution_id": execution_id,
            "manifest_sha256": manifest_sha,
            "phases": {"preflight": "success", "scgpt": "success"},
            "hashes": {
                key: hashlib.sha256(path.read_bytes()).hexdigest()
                for key, path in targets.items()
            },
        }))
        return run_id, execution_id, receipt_path, preflight_marker

    def run_evidence_audit(self, run_id, execution_id):
        return audit(
            self.preflight_path, self.graph_path, self.eval_path,
            self.reference_path, expected_manifest_sha="manifest",
            expected_rna_sha="rna", evidence_dir=self.root,
            run_id=run_id, execution_id=execution_id)

    def run_audit(self, expected_rna_sha="rna"):
        return audit(
            self.preflight_path, self.graph_path, self.eval_path,
            self.reference_path, expected_manifest_sha="manifest",
            expected_rna_sha=expected_rna_sha)

    def test_accepts_complete_consistent_results(self):
        result = self.run_audit()
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["historical_fields_unchanged"], 16)
        self.assertEqual(len(result["scgpt_metric_keys"]), 6)
        self.assertNotIn("evidence_manifest_sha256", result)

    def test_accepts_valid_evidence(self):
        run_id, execution_id, _receipt, _marker = self.create_evidence()
        result = self.run_evidence_audit(run_id, execution_id)
        self.assertEqual(result["run_id"], run_id)
        self.assertEqual(result["execution_id"], execution_id)
        self.assertTrue(result["evidence_manifest_sha256"].startswith(run_id))

    def test_rejects_tampered_output_hash(self):
        run_id, execution_id, _receipt, _marker = self.create_evidence()
        self.new_eval["scgpt__m0"] = 99.0
        self.eval_path.write_text(json.dumps(self.new_eval))
        with self.assertRaisesRegex(ValueError, "receipt hash mismatch: eval_json"):
            self.run_evidence_audit(run_id, execution_id)

    def test_rejects_tampered_receipt_hash(self):
        run_id, execution_id, receipt_path, _marker = self.create_evidence()
        receipt = json.loads(receipt_path.read_text())
        receipt["hashes"]["graph_npz"] = "0" * 64
        receipt_path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "receipt hash mismatch: graph_npz"):
            self.run_evidence_audit(run_id, execution_id)

    def test_rejects_missing_marker(self):
        run_id, execution_id, _receipt, marker = self.create_evidence()
        marker.unlink()
        with self.assertRaises(FileNotFoundError):
            self.run_evidence_audit(run_id, execution_id)

    def test_rejects_wrong_execution_id(self):
        run_id, _execution_id, _receipt, _marker = self.create_evidence()
        with self.assertRaises(FileNotFoundError):
            self.run_evidence_audit(run_id, "wrong-execution")

    def test_rejects_artifact_outside_evidence_dir(self):
        run_id, execution_id, _receipt, _marker = self.create_evidence()
        outside = self.root.parent / "outside-pbmc-eval.json"
        outside.write_bytes(self.eval_path.read_bytes())
        try:
            with self.assertRaisesRegex(ValueError, "outside evidence_dir"):
                audit(
                    self.preflight_path, self.graph_path, outside,
                    self.reference_path, expected_manifest_sha="manifest",
                    expected_rna_sha="rna", evidence_dir=self.root,
                    run_id=run_id, execution_id=execution_id)
        finally:
            outside.unlink(missing_ok=True)

    def test_rejects_changed_historical_metric(self):
        self.new_eval["embed__m0"] = 99.0
        self.eval_path.write_text(json.dumps(self.new_eval))
        with self.assertRaisesRegex(ValueError, "changed historical metric"):
            self.run_audit()

    def test_rejects_rna_hash_mismatch(self):
        with self.assertRaisesRegex(ValueError, "RNA hash mismatch"):
            self.run_audit(expected_rna_sha="different")

    def test_rejects_nonfinite_graph(self):
        graph = np.eye(1200, dtype=np.float32)
        graph[0, 0] = np.nan
        np.savez(
            self.graph_path,
            co=graph,
            sg=np.eye(1200, dtype=np.float32),
            cell_ids=np.arange(4000),
            genes=np.asarray([f"g{i}" for i in range(1200)]),
            manifest_sha=np.asarray("manifest"),
            selection_seed=np.asarray(20260713),
            pool_cap=np.asarray(4000),
            rna_sha256=np.asarray("rna"),
        )
        with self.assertRaisesRegex(ValueError, "non-finite"):
            self.run_audit()


if __name__ == "__main__":
    unittest.main()
