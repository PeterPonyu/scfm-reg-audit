#!/usr/bin/env python
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256_file(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path):
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        raise ValueError(f"invalid JSON evidence file: {path.name}") from exc


def _load_marker(path):
    marker = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            raise ValueError(f"invalid evidence marker line: {path.name}")
        key, value = line.split("=", 1)
        if not key or key in marker:
            raise ValueError(f"invalid evidence marker key: {path.name}")
        marker[key] = value
    return marker


def _require_marker(path, run_id, execution_id, manifest_sha, phase,
                    receipt_name=None):
    marker = _load_marker(path)
    expected = {
        "run_id": run_id,
        "execution_id": execution_id,
        "manifest_sha256": manifest_sha,
        "phase": phase,
        "status": "success",
    }
    if receipt_name is not None:
        expected["receipt"] = receipt_name
    for key, value in expected.items():
        if marker.get(key) != value:
            raise ValueError(f"evidence marker mismatch: {path.name}: {key}")


def _audit_evidence(evidence_dir, run_id, execution_id, preflight_file,
                    graph_file, eval_file, reference_file):
    if not run_id or not execution_id:
        raise ValueError("run_id and execution_id are required with evidence_dir")
    evidence = Path(evidence_dir)
    if not evidence.is_dir():
        raise FileNotFoundError(evidence)

    expected_artifacts = {
        preflight_file.resolve(): (evidence / "pbmc_preflight_remote.json").resolve(),
        graph_file.resolve(): (evidence / "pbmc_scgpt_pooled_v2.npz").resolve(),
        eval_file.resolve(): (evidence / "pbmc_eval_v2.json").resolve(),
    }
    for supplied, expected in expected_artifacts.items():
        if supplied != expected:
            raise ValueError(
                f"strict evidence artifact is outside evidence_dir: {supplied.name}")

    names = {
        "manifest_json": f"manifest.{run_id}.json",
        "manifest_verification_json":
            f"manifest_verification.{run_id}.{execution_id}.json",
        "preflight_marker": f"PREFLIGHT_DONE.{run_id}.{execution_id}",
        "scgpt_marker": f"SCGPT_DONE.{run_id}.{execution_id}",
        "receipt": f"run_receipt.{run_id}.{execution_id}.json",
        "preflight_log": f"preflight_run.{run_id}.{execution_id}.log",
        "scgpt_log": f"scgpt_run.{run_id}.{execution_id}.log",
    }
    evidence_paths = {key: evidence / name for key, name in names.items()}
    for path in evidence_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest = _load_json(evidence_paths["manifest_json"])
    manifest_sha = manifest.get("manifest_sha256")
    unhashed_manifest = dict(manifest)
    unhashed_manifest.pop("manifest_sha256", None)
    computed_manifest_sha = hashlib.sha256(
        json.dumps(unhashed_manifest, sort_keys=True).encode()).hexdigest()
    if manifest_sha != computed_manifest_sha:
        raise ValueError("evidence manifest self-hash mismatch")
    if len(run_id) != 16 or manifest_sha[:16] != run_id:
        raise ValueError("RUN_ID does not match evidence manifest hash prefix")

    verification = _load_json(evidence_paths["manifest_verification_json"])
    if verification.get("status") != "verified":
        raise ValueError("evidence manifest verification did not pass")
    if verification.get("manifest_sha256") != manifest_sha:
        raise ValueError("evidence manifest verification hash mismatch")

    receipt_path = evidence_paths["receipt"]
    _require_marker(
        evidence_paths["preflight_marker"], run_id, execution_id,
        manifest_sha, "preflight")
    _require_marker(
        evidence_paths["scgpt_marker"], run_id, execution_id,
        manifest_sha, "scgpt", receipt_path.name)

    receipt = _load_json(receipt_path)
    if receipt.get("run_id") != run_id:
        raise ValueError("evidence receipt run_id mismatch")
    if receipt.get("execution_id") != execution_id:
        raise ValueError("evidence receipt execution_id mismatch")
    if receipt.get("manifest_sha256") != manifest_sha:
        raise ValueError("evidence receipt manifest hash mismatch")
    if receipt.get("phases") != {"preflight": "success", "scgpt": "success"}:
        raise ValueError("evidence receipt phases did not pass")

    hash_targets = {
        "manifest_json": evidence_paths["manifest_json"],
        "manifest_verification_json": evidence_paths["manifest_verification_json"],
        "preflight_json": preflight_file,
        "eval_json": eval_file,
        "graph_npz": graph_file,
        "reference_json": reference_file,
        "preflight_log": evidence_paths["preflight_log"],
        "scgpt_log": evidence_paths["scgpt_log"],
    }
    receipt_hashes = receipt.get("hashes")
    if not isinstance(receipt_hashes, dict):
        raise ValueError("evidence receipt hashes missing")
    missing_hashes = set(hash_targets).difference(receipt_hashes)
    unknown_hashes = set(receipt_hashes).difference(hash_targets)
    if missing_hashes or unknown_hashes:
        raise ValueError(
            "evidence receipt hash set mismatch: "
            f"missing={sorted(missing_hashes)}, unknown={sorted(unknown_hashes)}")
    for key, path in hash_targets.items():
        if receipt_hashes[key] != sha256_file(path):
            raise ValueError(f"evidence receipt hash mismatch: {key}")
    return manifest_sha


def audit(preflight_path, graph_path, eval_path, reference_path,
          expected_manifest_sha=None, expected_rna_sha=None,
          evidence_dir=None, run_id=None, execution_id=None):
    paths = [Path(path) for path in (preflight_path, graph_path, eval_path, reference_path)]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    preflight_file, graph_file, eval_file, reference_file = paths
    preflight = json.loads(preflight_file.read_text())
    reference = json.loads(reference_file.read_text())
    new_eval = json.loads(eval_file.read_text())
    if preflight.get("status") != "passed":
        raise ValueError("remote preflight did not pass")
    if (expected_manifest_sha is not None
            and preflight.get("manifest_sha") != expected_manifest_sha):
        raise ValueError("remote preflight manifest hash mismatch")
    reference_sha = sha256_file(reference_file)
    if preflight.get("reference_metrics_sha256") != reference_sha:
        raise ValueError("remote preflight reference hash mismatch")
    keys = ["n_pairs", "n_tf", "tissue", "types"]
    keys.extend(sorted(key for key in reference if key.startswith(("embed__", "attn__"))))
    for key in keys:
        if preflight["metrics"].get(key) != reference.get(key):
            raise ValueError(f"remote preflight metric mismatch: {key}")
        if new_eval.get(key) != reference.get(key):
            raise ValueError(f"new evaluation changed historical metric: {key}")

    with np.load(graph_file, allow_pickle=False) as cache:
        required = {
            "co", "sg", "cell_ids", "genes", "manifest_sha", "selection_seed",
            "pool_cap", "rna_sha256",
        }
        missing = required.difference(cache.files)
        if missing:
            raise ValueError(f"PBMC scGPT cache missing keys: {sorted(missing)}")
        co, sg = cache["co"], cache["sg"]
        cell_ids = cache["cell_ids"]
        genes = cache["genes"]
        if co.shape != (1200, 1200) or sg.shape != (1200, 1200):
            raise ValueError("PBMC scGPT graph shape mismatch")
        if len(cell_ids) != 4000 or len(np.unique(cell_ids)) != 4000:
            raise ValueError("PBMC scGPT selected-cell IDs invalid")
        if len(genes) != 1200 or len(np.unique(genes)) != 1200:
            raise ValueError("PBMC scGPT gene IDs invalid")
        if int(cache["selection_seed"].item()) != 20260713:
            raise ValueError("PBMC scGPT selection seed mismatch")
        if int(cache["pool_cap"].item()) != 4000:
            raise ValueError("PBMC scGPT pool cap mismatch")
        if (expected_manifest_sha is not None
                and str(cache["manifest_sha"].item()) != expected_manifest_sha):
            raise ValueError("PBMC scGPT manifest hash mismatch")
        if (expected_rna_sha is not None
                and str(cache["rna_sha256"].item()) != expected_rna_sha):
            raise ValueError("PBMC scGPT RNA hash mismatch")
        if not np.isfinite(co).all() or not np.isfinite(sg).all():
            raise ValueError("PBMC scGPT graph contains non-finite values")

    scgpt_keys = sorted(key for key in new_eval if key.startswith("scgpt__"))
    if len(scgpt_keys) != 6:
        raise ValueError(f"expected 6 scGPT metrics, found {scgpt_keys}")
    if any(not np.isfinite(float(new_eval[key])) for key in scgpt_keys):
        raise ValueError("new scGPT metrics contain non-finite values")
    evidence_manifest_sha = None
    if evidence_dir is not None:
        evidence_manifest_sha = _audit_evidence(
            evidence_dir, run_id, execution_id, preflight_file, graph_file,
            eval_file, reference_file)
    elif run_id is not None or execution_id is not None:
        raise ValueError("evidence_dir is required with run_id or execution_id")

    result = {
        "status": "passed",
        "historical_fields_unchanged": len(keys),
        "scgpt_metric_keys": scgpt_keys,
        "graph_sha256": sha256_file(graph_file),
    }
    if evidence_manifest_sha is not None:
        result["evidence_manifest_sha256"] = evidence_manifest_sha
        result["run_id"] = run_id
        result["execution_id"] = execution_id
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("preflight")
    parser.add_argument("graph")
    parser.add_argument("evaluation")
    parser.add_argument("reference")
    parser.add_argument("--evidence-dir")
    parser.add_argument("--run-id")
    parser.add_argument("--execution-id")
    args = parser.parse_args()
    if bool(args.evidence_dir) != bool(args.run_id and args.execution_id):
        parser.error("--evidence-dir, --run-id, and --execution-id must be supplied together")
    root = Path(__file__).resolve().parents[2]
    manifest_path = root / "data/manifest/shared_genes.v2.json"
    rna_path = root / "data/multiome/pbmc10k_rna.h5ad"
    manifest_sha = json.loads(manifest_path.read_text())["sha256"]
    rna_sha = sha256_file(rna_path)
    print(json.dumps(audit(
        args.preflight, args.graph, args.evaluation, args.reference,
        expected_manifest_sha=manifest_sha, expected_rna_sha=rna_sha,
        evidence_dir=args.evidence_dir, run_id=args.run_id,
        execution_id=args.execution_id), indent=2))


if __name__ == "__main__":
    main()
