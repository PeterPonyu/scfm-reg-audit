#!/usr/bin/env python3
"""Single entry point: offline inventory by default; later authorized work explicitly."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

from screg_validation.common import ContractError, authorize, digest, load_json, sha256, write_json
from screg_validation.hardware import gpu_readiness
from screg_validation.models import READOUTS, READOUT_RESOURCES, check_resources
from screg_validation.replication import context, preparation_resources

HERE = Path(__file__).resolve().parent
TASKS = ("source-match", "training-overlap", "external-reference", "gate-validation",
         "randomization-stress", "replication-prepare", "replication-infer", "replication-analyze")
OPERATIONAL_TASKS = ("replication-pilot",)
EXECUTION_TASKS = TASKS + OPERATIONAL_TASKS


def input_path(resources, name):
    if name not in resources or not Path(resources[name].get("path", "")).is_file():
        raise ContractError(f"required explicit resource absent: {name}")
    resource = resources[name]
    if not resource.get("sha256") or sha256(resource["path"]) != resource["sha256"]:
        raise ContractError(f"resource hash missing or changed: {name}")
    return Path(resource["path"])


def validate_plan(plan):
    if plan.get("schema") != "screg-validation.v1" or set(plan.get("tasks", [])) != set(TASKS):
        raise ContractError("task registry schema incomplete")
    if plan.get("registered_before_new_scientific_outcomes") is not True:
        raise ContractError("new scientific outcomes require a separately frozen plan")
    if set(plan["replication"]["models"]["readouts"]) - set(READOUTS):
        raise ContractError("unsupported model readout")
    if plan["replication"]["dataset"]["scale"] != "counts":
        raise ContractError("independent replication including UCE requires explicit raw counts")
    for r in plan["downloads"]:
        if not 0 < r["max_bytes"] <= 5 * 1024**3 or not r["url"].startswith("https://"):
            raise ContractError("download resource lacks a bounded HTTPS contract")
    if len({r["id"] for r in plan["downloads"]}) != len(plan["downloads"]):
        raise ContractError("duplicate download IDs")
    stats = plan["replication"]["statistics"]
    if not (0 < stats["permutations"] <= 9999 and 100 <= stats["resamples"] <= 10000):
        raise ContractError("unbounded statistical budget")
    for section in (stats, plan["replication"]["models"], plan["randomization_stress"]):
        if not 0 < section.get("max_seconds", 0) <= 28800:
            raise ContractError("stage time budget missing or exceeds eight hours")


def preflight(plan, resources, *, verify_hashes=False, gpu=False):
    validate_plan(plan)
    checks = check_resources(resources, list(resources), verify_hashes=verify_hashes)
    tasks = []
    requirements = {"source-match": ["source_left", "source_right"] if {"source_left", "source_right"} <= resources.keys() else ["source_candidate"],
                    "training-overlap": ["training_inventory", "evaluation_identity"],
                    "external-reference": ["external_labels"], "gate-validation": ["biological_readout_classes"],
                    "randomization-stress": ["randomization_input"],
                    "replication-prepare": preparation_resources(plan),
                    "replication-infer": sorted({n for r in plan["replication"]["models"]["readouts"] for n in READOUT_RESOURCES[r]}),
                    "replication-analyze": []}
    for task in TASKS:
        missing = [name for name in requirements[task] if name not in checks or checks[name]["blockers"]]
        blockers = ["required_input:" + name for name in missing]
        if task == "replication-prepare":
            evidence = plan["replication"]["dataset"].get("genome_evidence")
            if not evidence:
                blockers.append("explicit_genome_build_evidence_not_yet_recorded")
            elif "multiome" in checks and not checks["multiome"]["blockers"]:
                from screg_validation.h5ad import genome_evidence
                import h5py
                try:
                    with h5py.File(resources["multiome"]["path"], "r") as f:
                        genome_evidence(f, plan["replication"]["dataset"])
                except (ContractError, OSError) as exc:
                    blockers.append("genome_evidence:" + str(exc))
        if task == "external-reference" and plan["external_reference"].get("score_replication_graphs"):
            blockers.append("verified_prepared_and_model_graph_receipts_required")
        if task == "replication-infer":
            from screg_validation.models import validate_adapter_layout
            if not blockers:
                try:
                    for readout in plan["replication"]["models"]["readouts"]:
                        validate_adapter_layout(resources, readout)
                except ContractError as exc:
                    blockers.append(str(exc))
            blockers.append("prepared_stage_receipt_required_before_inference")
            blockers.append("bounded_pilot_receipt_required_before_full_inference")
        if task == "replication-analyze":
            blockers.append("complete_prepared_and_inference_receipts_required")
        tasks.append({"task": task, "implemented": True,
                      "technical_inputs": "blocked" if blockers else "present_hash_verified" if verify_hashes else "present_hash_unchecked",
                      "blockers": blockers, "execution_authorized": False})
    tasks.append({"task": "replication-pilot", "implemented": True, "technical_inputs": "blocked",
                  "blockers": ["prepared_stage_receipt_required"], "execution_authorized": False})
    return {"status": "offline_preflight_completed", "plan_sha256": digest(plan), "resources": checks,
            "tasks": tasks, "free_bytes": shutil.disk_usage(HERE).free,
            "dependencies": {name: importlib.util.find_spec(name) is not None for name in
                             ("numpy", "scipy", "h5py", "pandas", "torch", "transformers", "sklearn", "MOODS", "pyfaidx", "einops")},
            "gpu": gpu_readiness(plan["gpu"]) if gpu else {"status": "not_queried"},
            "execution": {"data_downloads": 0, "model_inference": 0, "real_data_analyses": 0,
                          "reason": "this command inspects resources only; downloads, inference and analyses need a run record"}}


def run_task(plan, resources, task, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if task.startswith("replication-"):
        from screg_validation.replication import prepare_replication, pilot_replication, infer_replication, analyze_replication
        function = {"replication-prepare": prepare_replication, "replication-pilot": pilot_replication, "replication-infer": infer_replication,
                    "replication-analyze": analyze_replication}[task]
        return function(plan, resources, out / "replication")
    if task == "source-match":
        from screg_validation.provenance import compare_candidate, donor_join, validate_chain
        if "source_left" in resources and "source_right" in resources:
            from screg_validation.source_adapter import compare_h5ad_sources
            from screg_validation.h5ad import index
            import h5py
            result = compare_h5ad_sources(resources["source_left"], resources["source_right"],
                                         **plan["source_h5ad_match"])
            with h5py.File(input_path(resources, "source_left"), "r") as f:
                source_cells = index(f, "obs")
        else:
            arrays = np.load(input_path(resources, "source_candidate"), allow_pickle=False)
            source_cells = arrays["left_cells"].tolist()
            result = compare_candidate(source_cells, arrays["left_genes"].tolist(), arrays["left_matrix"],
                                       arrays["right_cells"].tolist(), arrays["right_genes"].tolist(), arrays["right_matrix"],
                                       **plan["source_match"])
        if "source_chain" in resources:
            result["chain_structure"] = validate_chain(load_json(input_path(resources, "source_chain")))
        if "source_donors" in resources:
            with input_path(resources, "source_donors").open() as f:
                joined = donor_join(source_cells, list(csv.DictReader(f)))
            result["donor_join"] = {"donors": len(set(joined)), "mapping_sha256": digest(joined)}
        result["provenance_closed"] = False
    elif task == "training-overlap":
        from screg_validation.provenance import training_overlap
        training = load_json(input_path(resources, "training_inventory"))
        evaluation = load_json(input_path(resources, "evaluation_identity"))
        result = {"rows": [training_overlap(row, evaluation[row["model"]]) for row in training]}
    elif task == "external-reference":
        if plan["external_reference"].get("score_replication_graphs"):
            from screg_validation.external_adapter import score_labels
            result = score_labels(input_path(resources, "external_labels"), out / "replication",
                                plan["replication"]["models"]["readouts"],
                                expected_tissue=plan["external_reference"]["expected_tissue"],
                                min_each=plan["external_reference"]["min_each"],
                                resamples=plan["replication"]["statistics"]["resamples"], seed=plan["seed"],
                                expected_context=context(plan, resources))
            result["plan_sha256"] = digest(plan)
            write_json(out / "external-reference.json", result)
            return result
        from screg_validation.statistics import external_metrics
        with input_path(resources, "external_labels").open() as f:
            rows = list(csv.DictReader(f))
        metrics = external_metrics(rows, **{k: plan["external_reference"][k] for k in ("expected_tissue", "min_each")})
        result = {"metrics": metrics, "scope": "paired external labels; verified negative and unlabeled tasks remain separate"}
    elif task == "gate-validation":
        from screg_validation.statistics import gate_validation
        result = gate_validation(load_json(input_path(resources, "biological_readout_classes")))
    elif task == "randomization-stress":
        from screg_validation.calibration import stress
        arrays = np.load(input_path(resources, "randomization_input"), allow_pickle=False)
        rule = {**plan["randomization_stress"]}
        budget = rule.pop("max_seconds")
        started = time.monotonic()
        def check_time():
            if time.monotonic() - started > budget:
                raise ContractError("registered stress time budget reached")
        result = stress(arrays["graph"], arrays["tf_indices"], arrays["controls"], arrays["strata"],
                        **rule, gpu_policy=plan["gpu"], check=check_time)
    else:
        raise ContractError("unregistered task")
    result["plan_sha256"] = digest(plan)
    write_json(out / (task + ".json"), result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=HERE / "plan.json")
    parser.add_argument("--resources", type=Path, default=HERE / "local_resources.json")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("--verify-hashes", action="store_true")
    p.add_argument("--gpu", action="store_true")
    p.add_argument("--out", type=Path)
    p = sub.add_parser("run")
    p.add_argument("--task", choices=EXECUTION_TASKS, required=True)
    p.add_argument("--execute", action="store_true")
    p.add_argument("--authorization", type=Path)
    p.add_argument("--out", type=Path, default=HERE / "runs")
    p = sub.add_parser("download")
    p.add_argument("--resource", required=True)
    p.add_argument("--execute", action="store_true")
    p.add_argument("--authorization", type=Path)
    p.add_argument("--out", type=Path, default=HERE / "inputs")
    args = parser.parse_args()
    plan, resources = load_json(args.plan), load_json(args.resources)
    validate_plan(plan)
    if args.command == "preflight":
        result = preflight(plan, resources, verify_hashes=args.verify_hashes, gpu=args.gpu)
        if args.out:
            write_json(args.out, result)
        print(json.dumps(result, ensure_ascii=False, indent=2))
    elif args.command == "run":
        if not args.execute:
            if args.authorization is not None:
                authorize(plan, args.authorization, "experiment", args.task)
            print(json.dumps({"status": "planned_only", "task": args.task, "plan_sha256": digest(plan),
                              "execution": False, "authorization_needed": True}, indent=2))
            return
        auth = authorize(plan, args.authorization, "experiment", args.task)
        if args.task in OPERATIONAL_TASKS:
            from screg_validation.replication import operational_contract
            contract = operational_contract()
            if args.task not in contract["extra_tasks"] or auth.get("operational_contract_sha256") != digest(contract):
                raise ContractError("operational pilot authorization absent or drifted")
        core_context = context(plan, resources)
        run_root = args.out.resolve() / (digest(core_context)[:20])
        # Never write experiment outputs into the manuscript or its deliverables.
        if not run_root.is_relative_to(HERE.resolve()):
            raise ContractError("experiment outputs must remain within this validation package")
        run_root.mkdir(parents=True, exist_ok=True)
        # Separate action receipts retain the original plan approval and all
        # extra inputs without changing the shared replication namespace.
        action_context = {"plan": digest(plan), "resources": digest(resources), "core": core_context}
        write_json(run_root / "authorizations" / (args.task + "-" + digest(auth)[:20] + ".json"),
                   {"authorization_sha256": digest(auth), **action_context})
        result = run_task(plan, resources, args.task, run_root)
        print(json.dumps({"task": args.task, "status": result.get("status", "complete"), "out": str(run_root)}, indent=2))
    else:
        resource = next((r for r in plan["downloads"] if r["id"] == args.resource), None)
        if resource is None:
            raise ContractError("download resource not registered")
        if not args.execute:
            if args.authorization is not None:
                authorize(plan, args.authorization, "download", args.resource)
            print(json.dumps({"status": "planned_only", "resource": resource, "execution": False}, indent=2))
            return
        if not args.out.resolve().is_relative_to(HERE.resolve()):
            raise ContractError("downloads must remain within this validation package")
        from screg_validation.fetch import fetch_registered, load_fetch_gate
        gate = load_fetch_gate(input_path(resources, "download_gate") if "download_gate" in resources else None)
        print(json.dumps(fetch_registered(plan, resource, args.out, args.authorization, gate=gate), indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ContractError, OSError, KeyError, ValueError) as exc:
        print(json.dumps({"status": "blocked", "reason": str(exc)}, ensure_ascii=False), file=sys.stderr)
        sys.exit(2)
