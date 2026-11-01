"""Staged donor/panel validation. Every expensive stage is separately authorized."""
from __future__ import annotations

import hashlib
from pathlib import Path
import time

import numpy as np
from scipy import sparse
from scipy.stats import rankdata, spearmanr

from .common import ContractError, cache_receipt, digest, load_json, sha256, write_json
from .design import donor_split, select_panel, validate_split, verify_panel
from .graphs import (build_motif_hits, donor_graph, gene_coexpression, independent_degree,
                     locus_links, read_coordinates, read_meme)
from .h5ad import counts_normalization, inspect, read_rows
from .models import PANEL_RESOURCES, READOUT_RESOURCES, require_resources, run_readout, tokenizable_genes
from .statistics import (adjust, component_randomization, donor_interval,
                         exact_paired_signflip, partial_spearman)
from .batched_null import batched_components
from .hardware import require_gpu


def source_signature():
    root = Path(__file__).resolve().parent
    files = sorted(root.rglob("*.py")) + [root.parent / "run.py", root.parent / "operational_contract.json"]
    return digest({p.relative_to(root.parent).as_posix(): sha256(p) for p in files})


def preparation_resources(plan):
    names = ["multiome", "gene_coordinates", "genome_fasta", "genome_fai", "meme", "original_panel",
             *PANEL_RESOURCES]
    if any(name in plan["replication"]["models"]["readouts"] for name in ("scfoundation_encoder", "uce_encoder")):
        names += ["scf_gene_index", "uce_chrom"]
    return names


def context(plan, resources):
    # Optional assay labels, source candidates and download settings must not
    # invalidate the frozen panel or force expensive model inference to repeat.
    names = set(preparation_resources(plan))
    for readout in plan["replication"]["models"]["readouts"]:
        names.update(READOUT_RESOURCES[readout])
    core = {"replication": plan["replication"], "seed": plan["seed"], "gpu": plan["gpu"]}
    return {"plan": digest(core), "resources": digest({n: resources.get(n) for n in sorted(names)}),
            "code": source_signature()}


def finish(directory, name, ctx, outputs):
    write_json(directory / name, {"context_sha256": digest(ctx),
                                  "outputs": {p: sha256(directory / p) for p in outputs}})


def require_prepared(plan, resources, directory):
    """Current input hashes and complete preparation manifest, never a subset."""
    directory = Path(directory)
    require_resources(resources, preparation_resources(plan))
    ctx = context(plan, resources)
    receipt = cache_receipt(directory / "prepared.receipt.json", ctx)
    if not receipt:
        raise ContractError("verified prepared stage receipt required")
    prepared, panel = load_json(directory / "prepared.json"), load_json(directory / "panel.json")
    verify_panel(panel)
    if (prepared.get("status") != "prepared" or prepared.get("panel_sha256") != panel["panel_sha256"]
            or prepared.get("split") != panel["split"] or prepared.get("source_context") != ctx):
        raise ContractError("prepared metadata incomplete or drifted")
    expected = {"panel.json", "links.npz", "motif_hits.npz", "prepared.json", "confounds.npz"}
    for donor in panel["split"]["degree_reference"] + panel["split"]["evaluation"]:
        stem = prepared["donor_stems"].get(donor)
        if stem != hashlib.sha256(donor.encode()).hexdigest()[:16]:
            raise ContractError("prepared donor identity changed")
        expected.update(stem + suffix for suffix in ("-raw.npz", "-cp10k.npz", "-graphs.npz", "-graph.json"))
    if set(receipt.get("outputs", {})) != expected:
        raise ContractError("prepared receipt output manifest incomplete")
    return prepared, panel


def select_cells(info, split, rule, seed):
    donors = np.asarray(info["donors"])
    labels = np.asarray(info["cell_types"])
    common = set(labels)
    for donor in sum(split.values(), []):
        supported = {label for label in set(labels[donors == donor])
                     if int(((labels == label) & (donors == donor)).sum()) >= rule["min_cells_per_type"]}
        common &= supported
    if rule.get("cell_types"):
        declared = set(rule["cell_types"])
        if not declared <= common:
            raise ContractError("declared cell types lack support in a frozen donor")
        common = declared
    if not common:
        raise ContractError("no cell type shared by every frozen donor; no outcome-driven donor exclusion")
    quota = max(1, rule["max_cells_per_donor"] // len(common))
    if quota < rule["min_cells_per_type"]:
        raise ContractError("cell cap too small for common-type balancing")
    # Match type composition exactly across discovery/reference/evaluation. The
    # quota uses metadata counts only, never RNA values or a successful effect.
    type_quota = {label: min(quota, min(int(((labels == label) & (donors == donor)).sum())
                                      for donor in sum(split.values(), []))) for label in common}
    result = {}
    for donor in sum(split.values(), []):
        rows = []
        for label in sorted(common):
            eligible = np.flatnonzero((donors == donor) & (labels == label))
            selected = sorted(eligible, key=lambda i: hashlib.sha256(f"{seed}:{info['cells'][i]}".encode()).hexdigest())[:type_quota[label]]
            rows.extend(selected)
        if len(rows) < rule["min_cells_per_donor"]:
            raise ContractError("frozen donor below minimum cell support")
        result[donor] = sorted(rows)
    return result, sorted(common)


def prepare_replication(plan, resources, directory, *, motif_builder=build_motif_hits):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    ctx = context(plan, resources)
    require_resources(resources, preparation_resources(plan))
    if cache_receipt(directory / "prepared.receipt.json", ctx):
        return require_prepared(plan, resources, directory)[0]
    spec = plan["replication"]["dataset"]
    path = resources["multiome"]["path"]
    info = inspect(path, spec)
    split_rule = plan["replication"]["split"]
    split = donor_split(info["donors"], **split_rule)
    validate_split(split, info["donors"])
    selection, common_types = select_cells(info, split, plan["replication"]["cells"], plan["seed"])
    genes = [g for g, t in zip(info["features"], info["feature_types"]) if t == spec["rna_value"]]
    rna_cols = np.flatnonzero(np.array(info["feature_types"]) == spec["rna_value"])
    atac_cols = np.flatnonzero(np.array(info["feature_types"]) == spec["atac_value"])
    peaks = [info["features"][i] for i in atac_cols]
    coordinates = read_coordinates(resources["gene_coordinates"]["path"], plan["replication"]["promoter_bp"])
    supported_genes = sorted(set(genes) & set(coordinates))
    candidate_links = locus_links(peaks, supported_genes, coordinates)
    has_peak = {g for g, count in zip(supported_genes, candidate_links.getnnz(0)) if count > 0}
    all_models = any(name in plan["replication"]["models"]["readouts"]
                     for name in ("scfoundation_encoder", "uce_encoder"))
    eligible = has_peak & tokenizable_genes(resources, require_all_models=all_models)
    tf_symbols = {g for m in read_meme(resources["meme"]["path"]).values() for g in m["symbols"]}
    discovery_rows = sum([selection[d] for d in split["discovery"]], [])
    discovery = read_rows(path, spec["matrix"], discovery_rows, rna_cols)
    counts_normalization(discovery, spec["scale"])
    donor_labels = [info["donors"][i] for i in discovery_rows]
    original = load_json(resources["original_panel"]["path"])["genes"]
    panel = select_panel(discovery, genes, donor_labels, split, eligible_genes=eligible,
                         tf_symbols=tf_symbols, original_genes=original, **plan["replication"]["panel"])
    panel["selection_receipt"] = {"source_sha256": resources["multiome"]["sha256"],
                                  "discovery_cell_ids_sha256": digest([info["cells"][i] for i in discovery_rows])}
    panel["panel_sha256"] = digest({k: v for k, v in panel.items() if k != "panel_sha256"})
    write_json(directory / "panel.json", panel)
    gene_lookup = {g: i for i, g in enumerate(genes)}
    panel_cols = [gene_lookup[g] for g in panel["genes"]]
    links = locus_links(peaks, panel["genes"], coordinates)
    keep = np.flatnonzero(links.getnnz(1))
    if not len(keep):
        raise ContractError("no peaks in frozen panel loci")
    kept_peaks = [peaks[i] for i in keep]
    links = links[keep].tocsr()
    hits, peak_gc, peak_length = motif_builder(kept_peaks, panel["tf_genes"],
                                             resources["genome_fasta"]["path"], resources["meme"]["path"],
                                             pvalue=plan["replication"]["motif_pvalue"])
    sparse.save_npz(directory / "links.npz", links)
    sparse.save_npz(directory / "motif_hits.npz", hits)
    tf_indices = [panel["genes"].index(t) for t in panel["tf_genes"]]
    outputs = ["panel.json", "links.npz", "motif_hits.npz"]
    # Reference donors are processed before evaluation donors; no held-out links feed back.
    donors_to_graph = split["degree_reference"] + split["evaluation"]
    for donor in donors_to_graph:
        stem = hashlib.sha256(donor.encode()).hexdigest()[:16]
        full_rna = read_rows(path, spec["matrix"], selection[donor], rna_cols)
        cp, log = counts_normalization(full_rna, spec["scale"])
        raw, cp, log = full_rna[:, panel_cols], cp[:, panel_cols], log[:, panel_cols]
        atac = read_rows(path, spec["matrix"], selection[donor], atac_cols[keep])
        graph, graph_info = donor_graph(atac, log, links, hits, tf_indices,
                                       accessibility_floor=plan["replication"]["accessibility_floor"],
                                       linkage=plan["replication"]["linkage"],
                                       correlation_floor=plan["replication"]["link_correlation_floor"])
        co = gene_coexpression(log)
        sparse.save_npz(directory / (stem + "-raw.npz"), raw)
        sparse.save_npz(directory / (stem + "-cp10k.npz"), cp)
        np.savez_compressed(directory / (stem + "-graphs.npz"), proxy=graph, coexpression=co)
        write_json(directory / (stem + "-graph.json"), {**graph_info, "donor": donor,
                    "role": "degree_reference" if donor in split["degree_reference"] else "evaluation",
                    "cell_ids_sha256": digest([info["cells"][i] for i in selection[donor]])})
        outputs += [stem + suffix for suffix in ("-raw.npz", "-cp10k.npz", "-graphs.npz", "-graph.json")]
    n = len(panel["genes"])
    pc = np.asarray(links.sum(0)).ravel()
    gc = np.asarray(links.T @ (peak_gc * peak_length)).ravel() / np.maximum(1, np.asarray(links.T @ peak_length).ravel())
    lengths = [coordinates[g][3] for g in panel["genes"]]
    np.savez_compressed(directory / "confounds.npz", peak_count=pc, gene_length=lengths,
                        detection=[panel["detection"][g] for g in panel["genes"]], gc=gc,
                        opportunity_degree=np.asarray(((hits.T @ links) > 0).sum(0)).ravel())
    prepared = {"status": "prepared", "split": split, "genes": n, "tfs": len(tf_indices),
                "genome_evidence": info["genome_evidence"],
                "peaks_scanned": len(keep), "common_cell_types": common_types,
                "cell_composition_rule": "identical metadata-only per-type quotas in every discovery/reference/evaluation donor",
                "cells_per_donor": {d: len(selection[d]) for d in selection},
                "donor_stems": {d: hashlib.sha256(d.encode()).hexdigest()[:16] for d in donors_to_graph},
                "panel_sha256": panel["panel_sha256"], "source_context": ctx,
                "scope": "new discovery-selected panel, separate degree-reference donors and held-out evaluation donors; RNA/ATAC proxy remains operational"}
    write_json(directory / "prepared.json", prepared)
    outputs += ["prepared.json", "confounds.npz"]
    finish(directory, "prepared.receipt.json", ctx, outputs)
    return prepared


def operational_contract():
    path = Path(__file__).resolve().parents[1] / "operational_contract.json"
    contract = load_json(path)
    rule = contract.get("pilot", {})
    if (contract.get("schema") != "screg-operational.v1"
            or contract.get("scientific_plan_changed") is not False
            or rule.get("max_batches_per_readout") != 1
            or not 0 < rule.get("max_seconds", 0) <= 1800
            or rule.get("compare_real_cpu_forward") is not True
            or not 0 < rule.get("cpu_cuda_atol", 0) <= .02
            or not 0 < rule.get("cpu_cuda_rtol", 0) <= .02
            or rule.get("min_comparable_pairs") != 1):
        raise ContractError("invalid bounded pilot operational contract")
    return contract


def pilot_context(plan, resources, panel, directory=None):
    ctx = {**context(plan, resources), "stage": "pilot",
           "panel_sha256": panel["panel_sha256"],
           "operational_contract_sha256": digest(operational_contract())}
    if directory is not None:
        ctx["prepared_receipt_sha256"] = sha256(Path(directory) / "prepared.receipt.json")
    return ctx


def require_pilot(plan, resources, directory):
    directory = Path(directory)
    prepared, panel = require_prepared(plan, resources, directory)
    pctx = pilot_context(plan, resources, panel, directory)
    pilot_dir = directory / "inference_pilot"
    pilot_receipt = cache_receipt(pilot_dir / "pilot.receipt.json", pctx)
    if not pilot_receipt:
        raise ContractError("successful bounded pilot receipt required before full inference")
    report = load_json(pilot_dir / "pilot.json")
    rule = plan["replication"]["models"]
    if (report.get("status") != "pilot_complete" or report.get("readouts") != rule["readouts"]
            or report.get("donor") != panel["split"]["evaluation"][0]
            or report.get("device") != rule["device"]
            or report.get("max_batches_per_readout") != 1
            or report.get("panel_sha256") != panel["panel_sha256"]
            or report.get("source_context") != pctx
            or report.get("donor_cells") != prepared["cells_per_donor"][report.get("donor")]):
        raise ContractError("pilot accounting does not match frozen readouts/device/donor")
    expected = {"pilot.json"}
    for readout in rule["readouts"]:
        require_resources(resources, READOUT_RESOURCES[readout])
        expected.update(readout + "/" + name for name in ("result.receipt.json", "pilot.json", "model_graph.npz", "cpu_graph.npz"))
        child = pilot_dir / readout
        child_receipt = cache_receipt(child / "result.receipt.json", {**pctx, "readout": readout})
        if not child_receipt:
            raise ContractError("pilot child receipt missing: " + readout)
        if set(child_receipt["outputs"]) != {"pilot.json", "model_graph.npz", "cpu_graph.npz"}:
            raise ContractError("pilot child receipt manifest incomplete")
        row = load_json(child / "pilot.json")
        if (row.get("status") != "pilot_success" or row.get("readout") != readout
                or row.get("device") != rule["device"]
                or row.get("cells") != min(rule["batch"], report["donor_cells"])
                or row.get("batches") != 1 or not row.get("finite_observed_output")
                or not row.get("cpu_cuda_comparison", {}).get("passed")
                or row.get("panel_sha256") != panel["panel_sha256"]
                or row.get("donor") != report["donor"]
                or not 0 < row.get("elapsed_seconds", 0) <= operational_contract()["pilot"]["max_seconds"]):
            raise ContractError("pilot child incomplete or unsuccessful: " + readout)
        required_hash = digest({n: resources[n]["sha256"] for n in READOUT_RESOURCES[readout]})
        if row.get("resources_sha256") != required_hash:
            raise ContractError("pilot checkpoint identity drifted")
        comparison = row["cpu_cuda_comparison"]
        thresholds = operational_contract()["pilot"]
        if (comparison.get("cpu_device") != "cpu"
                or comparison.get("comparable_pairs", 0) < thresholds["min_comparable_pairs"]
                or comparison.get("atol") != thresholds["cpu_cuda_atol"]
                or comparison.get("rtol") != thresholds["cpu_cuda_rtol"]
                or not np.isfinite(comparison.get("max_absolute_difference", np.nan))
                or not 0 <= comparison["max_absolute_difference"] <= thresholds["cpu_cuda_atol"] + thresholds["cpu_cuda_rtol"]):
            raise ContractError("pilot numeric comparison incomplete or failed")
        graph = np.load(child / "model_graph.npz", allow_pickle=False)["graph"]
        if graph.shape != (len(panel["genes"]), len(panel["genes"])) or np.isinf(graph).any() or not np.isfinite(graph[~np.eye(len(panel["genes"]), dtype=bool)]).any():
            raise ContractError("pilot output shape or finite health invalid")
        cpu_graph = np.load(child / "cpu_graph.npz", allow_pickle=False)["graph"]
        common = np.isfinite(graph) & ~np.eye(len(panel["genes"]), dtype=bool)
        if (cpu_graph.shape != graph.shape or not np.array_equal(np.isfinite(cpu_graph), np.isfinite(graph))
                or not np.allclose(graph[common], cpu_graph[common], atol=thresholds["cpu_cuda_atol"], rtol=thresholds["cpu_cuda_rtol"])
                or int(common.sum()) != comparison["comparable_pairs"]
                or not np.isclose(float(np.max(np.abs(graph[common] - cpu_graph[common]))), comparison["max_absolute_difference"], atol=1e-12, rtol=1e-12)):
            raise ContractError("pilot saved CPU/CUDA graph comparison failed")
        if rule["device"] == "cuda" and (not row.get("gpu_before", {}).get("ready")
                or not row.get("gpu_after", {}).get("ready")
                or type(row.get("peak_memory_bytes")) is not int or row["peak_memory_bytes"] <= 0):
            raise ContractError("pilot CUDA telemetry incomplete")
    if set(pilot_receipt["outputs"]) != expected:
        raise ContractError("pilot receipt manifest incomplete")
    return report


def pilot_replication(plan, resources, directory, *, readout_runner=run_readout):
    """One real batch per readout, separate from donor inference and statistics."""
    directory = Path(directory)
    ctx = context(plan, resources)
    prepared, panel = require_prepared(plan, resources, directory)
    pctx = pilot_context(plan, resources, panel, directory)
    pilot_dir = directory / "inference_pilot"
    pilot_dir.mkdir(exist_ok=True)
    if cache_receipt(pilot_dir / "pilot.receipt.json", pctx):
        return require_pilot(plan, resources, directory)
    rule = plan["replication"]["models"]
    contract = operational_contract()["pilot"]
    donor = panel["split"]["evaluation"][0]
    stem = prepared["donor_stems"][donor]
    all_cp = sparse.load_npz(directory / (stem + "-cp10k.npz"))
    all_raw = sparse.load_npz(directory / (stem + "-raw.npz"))
    cp, raw = all_cp[:rule["batch"]], all_raw[:rule["batch"]]
    if not len(cp.shape) == 2 or cp.shape[0] == 0:
        raise ContractError("empty pilot batch")
    started = time.monotonic()
    outputs = ["pilot.json"]
    def check():
        if time.monotonic() - started >= contract["max_seconds"]:
            raise ContractError("bounded pilot time budget reached")
        if rule["device"] == "cuda":
            require_gpu(plan["gpu"])
    for readout in rule["readouts"]:
        child = pilot_dir / readout
        child.mkdir(exist_ok=True)
        child_context = {**pctx, "readout": readout}
        if not cache_receipt(child / "result.receipt.json", child_context):
            check()
            tick = time.monotonic()
            try:
                graph, telemetry = readout_runner(readout, panel["genes"], cp, raw, resources,
                    device=rule["device"], batch=rule["batch"], seed=plan["seed"], gpu_policy=plan["gpu"],
                    inference_context=rule["inference_context"], check=check)
                check()
                # Same real weights, cells and deterministic seed on CPU; this is
                # numeric health evidence, not an extra biological observation.
                cpu_graph, cpu_report = readout_runner(readout, panel["genes"], cp, raw, resources,
                    device="cpu", batch=rule["batch"], seed=plan["seed"], gpu_policy=plan["gpu"],
                    inference_context=rule["inference_context"], check=check)
                check()
                finite = np.isfinite(graph)
                common = finite & np.isfinite(cpu_graph)
                off_diagonal = ~np.eye(len(panel["genes"]), dtype=bool)
                common &= off_diagonal
                pair_count = int(common.sum())
                if (np.isinf(graph).any() or pair_count < contract["min_comparable_pairs"]
                        or not np.array_equal(finite, np.isfinite(cpu_graph))):
                    raise ContractError("pilot has undefined observed output or CPU/CUDA coverage differs")
                error = np.abs(graph[common] - cpu_graph[common])
                passed = np.allclose(graph[common], cpu_graph[common],
                                     atol=contract["cpu_cuda_atol"], rtol=contract["cpu_cuda_rtol"])
                row = {**telemetry, "status": "pilot_success" if passed else "pilot_failed",
                       "readout": readout, "donor": donor, "cells": cp.shape[0], "batches": 1,
                       "device": rule["device"], "elapsed_seconds": time.monotonic() - tick,
                       "panel_sha256": panel["panel_sha256"], "finite_observed_output": True,
                       "missing_gene_count": len(telemetry["unobserved_genes"]),
                       "cpu_cuda_comparison": {"passed": bool(passed), "comparable_pairs": pair_count,
                           "max_absolute_difference": float(error.max()), "atol": contract["cpu_cuda_atol"],
                           "rtol": contract["cpu_cuda_rtol"], "cpu_device": cpu_report.get("device"),
                           "scope": "same checkpoint, panel, cells and seed; actual CPU model forward"}}
                np.savez_compressed(child / "model_graph.npz", graph=graph)
                np.savez_compressed(child / "cpu_graph.npz", graph=cpu_graph)
                write_json(child / "pilot.json", row)
                if not passed:
                    raise ContractError("pilot real CPU/CUDA model graph comparison failed")
                finish(child, "result.receipt.json", child_context, ["pilot.json", "model_graph.npz", "cpu_graph.npz"])
            except Exception as exc:
                write_json(child / "failure.json", {"status": "pilot_failed", "readout": readout,
                    "reason": str(exc), "context": child_context})
                raise
        outputs += [readout + "/result.receipt.json", readout + "/pilot.json", readout + "/model_graph.npz", readout + "/cpu_graph.npz"]
    report = {"status": "pilot_complete", "readouts": rule["readouts"], "donor": donor,
              "donor_cells": all_cp.shape[0], "device": rule["device"], "max_batches_per_readout": 1,
              "panel_sha256": panel["panel_sha256"], "source_context": pctx,
              "scope": "bounded engineering pilot only; never donor inference or scientific evidence"}
    write_json(pilot_dir / "pilot.json", report)
    finish(pilot_dir, "pilot.receipt.json", pctx, outputs)
    return require_pilot(plan, resources, directory)


def coverage_status(graph, report, panel, rule):
    genes = panel["genes"]
    n = len(genes)
    graph = np.asarray(graph)
    counts = np.asarray(report.get("coverage_per_gene", []))
    if (graph.shape != (n, n) or np.isinf(graph).any() or counts.shape != (n,)
            or not np.isfinite(counts).all() or np.any(counts < 0)):
        raise ContractError("invalid graph or gene coverage")
    missing = [g for g, c in zip(genes, counts) if c == 0]
    if report.get("unobserved_genes") != missing:
        raise ContractError("unobserved gene list contradicts coverage")
    unobserved = np.flatnonzero(counts == 0)
    if len(unobserved) and (np.isfinite(graph[unobserved, :]).any() or np.isfinite(graph[:, unobserved]).any()):
        raise ContractError("unobserved graph scores were finite; forbidden zero imputation")
    ii, jj = np.repeat([genes.index(t) for t in panel["tf_genes"]], n), np.tile(np.arange(n), len(panel["tf_genes"]))
    ii, jj = ii[ii != jj], jj[ii != jj]
    gene_coverage = 1 - len(missing) / n
    pair_coverage = float(np.isfinite(graph[ii, jj]).mean())
    return {"gene_coverage_fraction": gene_coverage, "evaluated_pair_coverage_fraction": pair_coverage,
            "coverage_threshold_pass": gene_coverage >= rule["min_gene_coverage"],
            "primary_statistics_eligible": gene_coverage >= rule["min_gene_coverage"] and pair_coverage == 1.}


def require_inference(plan, resources, directory):
    directory = Path(directory)
    prepared, panel = require_prepared(plan, resources, directory)
    require_pilot(plan, resources, directory)
    ctx = context(plan, resources)
    receipt = cache_receipt(directory / "inference.receipt.json", ctx)
    if not receipt:
        raise ContractError("complete inference stage receipt required before analysis")
    report = load_json(directory / "inference.json")
    tasks = [prepared["donor_stems"][d] + "-" + r
             for d in panel["split"]["evaluation"] for r in plan["replication"]["models"]["readouts"]]
    expected = {"inference.json", *(t + "/result.receipt.json" for t in tasks)}
    if (set(receipt["outputs"]) != expected or report.get("tasks") != tasks
            or report.get("status") != "inference_complete"):
        raise ContractError("inference receipt manifest or task accounting incomplete")
    for donor in panel["split"]["evaluation"]:
        for readout in plan["replication"]["models"]["readouts"]:
            child = directory / (prepared["donor_stems"][donor] + "-" + readout)
            tctx = {**ctx, "donor": donor, "readout": readout, "panel_sha256": panel["panel_sha256"]}
            record = cache_receipt(child / "result.receipt.json", tctx)
            if not record or set(record["outputs"]) != {"model_graph.npz", "coverage.json"}:
                raise ContractError("inference child receipt incomplete")
            graph = np.load(child / "model_graph.npz", allow_pickle=False)["graph"]
            coverage = load_json(child / "coverage.json")
            status = coverage_status(graph, coverage, panel, plan["replication"]["models"])
            if any(coverage.get(k) != v for k, v in status.items()):
                raise ContractError("inference coverage report inconsistent with graph")
    return report


def infer_replication(plan, resources, directory, *, readout_runner=run_readout):
    directory = Path(directory)
    ctx = context(plan, resources)
    prepared, panel = require_prepared(plan, resources, directory)
    require_pilot(plan, resources, directory)
    tasks, rejected = [], []
    rule = plan["replication"]["models"]
    started = time.monotonic()
    def inference_check():
        if time.monotonic() - started >= rule["max_seconds"]:
            raise ContractError("registered inference time budget reached; completed receipts preserved")
        if rule["device"] == "cuda":
            require_gpu(plan["gpu"])
    n = len(panel["genes"])
    ii, jj = np.repeat([panel["genes"].index(t) for t in panel["tf_genes"]], n), np.tile(np.arange(n), len(panel["tf_genes"]))
    ii, jj = ii[ii != jj], jj[ii != jj]
    for donor in panel["split"]["evaluation"]:
        stem = prepared["donor_stems"][donor]
        cp = sparse.load_npz(directory / (stem + "-cp10k.npz"))
        raw = sparse.load_npz(directory / (stem + "-raw.npz"))
        for readout in rule["readouts"]:
            task = directory / (stem + "-" + readout)
            task.mkdir(exist_ok=True)
            task_context = {**ctx, "donor": donor, "readout": readout, "panel_sha256": panel["panel_sha256"]}
            if cache_receipt(task / "result.receipt.json", task_context):
                report = load_json(task / "coverage.json")
            else:
                if time.monotonic() - started >= rule["max_seconds"]:
                    raise ContractError("registered inference time budget reached; completed model receipts preserved")
                # Readiness is checked before every model task. Errors preserve previous receipts.
                graph, report = readout_runner(readout, panel["genes"], cp, raw, resources,
                        device=rule["device"], batch=rule["batch"], seed=plan["seed"], gpu_policy=plan["gpu"],
                        inference_context=rule["inference_context"], check=inference_check)
                inference_check()
                report.update(coverage_status(graph, report, panel, rule))
                # NaN entries represent unavailable genes and remain explicit in NPZ.
                np.savez_compressed(task / "model_graph.npz", graph=graph)
                write_json(task / "coverage.json", report)
                finish(task, "result.receipt.json", task_context, ["model_graph.npz", "coverage.json"])
            if not report["primary_statistics_eligible"]:
                rejected.append({"donor": donor, "readout": readout,
                    "reason": "coverage_below_frozen_threshold" if not report["coverage_threshold_pass"] else "incomplete_frozen_pair_coverage",
                    "gene_coverage_fraction": report["gene_coverage_fraction"],
                    "evaluated_pair_coverage_fraction": report["evaluated_pair_coverage_fraction"]})
            tasks.append(str(task.relative_to(directory)))
    report = {"status": "inference_complete", "tasks": tasks, "coverage_rejections": rejected,
              "scope": "90% token coverage is a minimum check; primary statistics additionally require every frozen TF-target pair finite; exclusions persist on resume"}
    write_json(directory / "inference.json", report)
    finish(directory, "inference.receipt.json", ctx, ["inference.json"] + [t + "/result.receipt.json" for t in tasks])
    return report


def restricted_strata(panel, confounds):
    n = len(panel["genes"])
    def quantile_bin(values, points):
        values = np.asarray(values)
        return np.digitize(values, np.unique(np.quantile(values, points)))
    opportunity = quantile_bin(confounds["opportunity_degree"], [.2, .4, .6, .8])
    gc = quantile_bin(confounds["gc"], [1/3, 2/3])
    return opportunity * 3 + gc + 15 * np.isin(np.arange(n), [panel["genes"].index(g) for g in panel["tf_genes"]])


def analyze_replication(plan, resources, directory):
    directory = Path(directory)
    ctx = context(plan, resources)
    prepared, panel = require_prepared(plan, resources, directory)
    require_inference(plan, resources, directory)
    panel, prepared = load_json(directory / "panel.json"), load_json(directory / "prepared.json")
    verify_panel(panel)
    tf_indices = [panel["genes"].index(t) for t in panel["tf_genes"]]
    n = len(panel["genes"])
    ii, jj = np.repeat(tf_indices, n), np.tile(np.arange(n), len(tf_indices))
    ii, jj = ii[ii != jj], jj[ii != jj]
    confounds = dict(np.load(directory / "confounds.npz", allow_pickle=False))
    common = np.column_stack([confounds[k][jj] for k in ("peak_count", "gene_length", "detection", "gc")])
    reference = [np.load(directory / (prepared["donor_stems"][d] + "-graphs.npz"))["proxy"]
                 for d in panel["split"]["degree_reference"]]
    rows, exclusions = [], []
    started = time.monotonic()
    statistics = plan["replication"]["statistics"]
    for donor in panel["split"]["evaluation"]:
        stem = prepared["donor_stems"][donor]
        arrays = dict(np.load(directory / (stem + "-graphs.npz"), allow_pickle=False))
        g, co = arrays["proxy"], arrays["coexpression"]
        od, ind = independent_degree(reference, panel["split"]["degree_reference"], donor, panel)
        baseline_controls = np.column_stack([common, od[ii], ind[jj]])
        baseline_rho = partial_spearman(co[ii, jj], g[ii, jj], baseline_controls)
        eligible_names, eligible_scores = [], []
        for readout in plan["replication"]["models"]["readouts"]:
            task = directory / (stem + "-" + readout)
            task_context = {**ctx, "donor": donor, "readout": readout, "panel_sha256": panel["panel_sha256"]}
            if not cache_receipt(task / "result.receipt.json", task_context):
                raise ContractError(f"inference not complete: {donor}/{readout}")
            scores = np.load(task / "model_graph.npz", allow_pickle=False)["graph"]
            coverage = load_json(task / "coverage.json")
            coverage_check = coverage_status(scores, coverage, panel, plan["replication"]["models"])
            if not coverage_check["primary_statistics_eligible"]:
                exclusions.append({"donor": donor, "readout": readout, "reason": "coverage_below_frozen_threshold" if not coverage_check["coverage_threshold_pass"] else "incomplete_frozen_pair_coverage; not zero imputed"})
                continue
            eligible_names.append(readout)
            eligible_scores.append(scores)
        if not eligible_names:
            continue
        controls = np.column_stack([common, rankdata(co[ii, jj])])
        def time_check():
            if time.monotonic() - started > statistics["max_seconds"]:
                raise ContractError("registered analysis time budget reached; preserve completed inference")
        for degree_mode in statistics["degree_modes"]:
            for null_mode in statistics["null_modes"]:
                strata = (np.zeros(n, int) if null_mode == "unrestricted"
                          else restricted_strata(panel, confounds))
                # Each donor/specification gets a receipt, so interrupted analysis resumes.
                statistics_dir = directory / (stem + "-statistics-" + degree_mode + "-" + null_mode)
                statistics_dir.mkdir(exist_ok=True)
                statistics_context = {**ctx, "donor": donor, "degree": degree_mode, "null": null_mode,
                                      "eligible_readouts": eligible_names}
                if cache_receipt(statistics_dir / "result.receipt.json", statistics_context):
                    completed_rows = load_json(statistics_dir / "rows.json")
                else:
                    tested = batched_components(np.stack(eligible_scores), g, tf_indices, controls,
                             strata=strata, permutations=statistics["permutations"], seed=plan["seed"],
                             degree_mode=degree_mode, reference_degree=(od, ind),
                             backend=statistics.get("backend", "cpu"), gpu_policy=plan["gpu"], check=time_check,
                             checkpoint=statistics_dir / "draw-progress.json",
                             checkpoint_context=statistics_context)
                    completed_rows = []
                    for readout, scores, row in zip(eligible_names, eligible_scores, tested):
                        row.pop("nulls")
                        # Identical outcome/pairs/nuisance, not mismatched adjusted statistics.
                        fm_rho = partial_spearman(scores[ii, jj], g[ii, jj], baseline_controls)
                        row.update(donor=donor, readout=readout, null_mode=null_mode,
                                   comparable_marginal_fm=fm_rho, comparable_marginal_baseline=baseline_rho,
                                   comparable_rho_difference=fm_rho - baseline_rho,
                                   concordance=float(spearmanr(scores[ii, jj], co[ii, jj]).statistic))
                        completed_rows.append(row)
                    write_json(statistics_dir / "rows.json", completed_rows)
                    finish(statistics_dir, "result.receipt.json", statistics_context, ["rows.json"])
                rows.extend(completed_rows)
    if not rows:
        raise ContractError("no eligible model/donor rows")
    registered_count = (len(panel["split"]["evaluation"]) * len(plan["replication"]["models"]["readouts"])
                        * len(statistics["degree_modes"]) * len(statistics["null_modes"]))
    family_p = [r["pIUT"] for r in rows] + [1.] * (registered_count - len(rows))
    # An untestable registered readout cannot shrink the correction family.
    qby, pholm = adjust(family_p)[:len(rows)], adjust(family_p, "Holm")[:len(rows)]
    for r, q, p in zip(rows, qby, pholm):
        r.update(q_IUT_BY=float(q), p_IUT_Holm=float(p))
    summary = []
    expected = panel["split"]["evaluation"]
    for readout in plan["replication"]["models"]["readouts"]:
        base = [r for r in rows if r["readout"] == readout and r["degree_mode"] == "reference"
                and r["null_mode"] == "restricted"]
        if set(r["donor"] for r in base) != set(expected):
            summary.append({"readout": readout, "donor_interval": None,
                            "reason": "missing frozen donor/readout coverage; no selective donor removal"})
            continue
        summary.append({"readout": readout,
                        "donor_interval": donor_interval([r["observed"] for r in base], [r["donor"] for r in base],
                             expected_donors=expected, resamples=statistics["resamples"], seed=plan["seed"]),
                        "paired_baseline_interval": donor_interval([r["comparable_rho_difference"] for r in base],
                             [r["donor"] for r in base], expected_donors=expected,
                             resamples=statistics["resamples"], seed=plan["seed"]),
                        "paired_baseline_signflip_p": exact_paired_signflip([r["comparable_rho_difference"] for r in base])})
    valid = [r for r in summary if r.get("paired_baseline_signflip_p") is not None]
    if valid:
        p_family = [r["paired_baseline_signflip_p"] for r in valid] + [1.] * (len(summary) - len(valid))
        for r, p in zip(valid, adjust(p_family, "Holm")[:len(valid)]):
            r["paired_baseline_signflip_Holm"] = float(p)
    report = {"status": "complete", "rows": rows, "donor_summary": summary, "exclusions": exclusions,
              "registered_hypotheses": registered_count, "untestable_hypotheses_retained_at_p1": registered_count - len(rows),
              "family": "all registered donor x readout x degree-mode x null-mode rows; untestable rows retained at p=1 in pooled IUT-BY/Holm",
              "scope": "new held-out-donor panel extension; conditional randomization and donor bootstrap; no claim of causal recovery or unseen pretraining",
              "panel_sha256": panel["panel_sha256"], "source_context": ctx}
    write_json(directory / "replication.json", report)
    finish(directory, "analysis.receipt.json", ctx, ["replication.json"])
    return report
