"""Explicit local resources, lazy adapters and recorded per-gene token coverage."""
from __future__ import annotations

import importlib
import json
from pathlib import Path
import pickle
import time

import numpy as np
from scipy import sparse

from .common import ContractError, digest, sha256
from .hardware import require_gpu


READOUTS = ("geneformer_embed", "geneformer_attention", "random_init_geneformer",
            "scgpt_encoder", "scfoundation_encoder", "uce_encoder")
PANEL_RESOURCES = ("gf_tokens", "gf_medians", "gf_name_ids", "scgpt_vocab")
READOUT_RESOURCES = {
    "geneformer_embed": ("gf_checkpoint", "gf_config", *PANEL_RESOURCES),
    "geneformer_attention": ("gf_checkpoint", "gf_config", *PANEL_RESOURCES),
    "random_init_geneformer": ("gf_config", *PANEL_RESOURCES),
    "scgpt_encoder": ("scgpt_checkpoint", "scgpt_args", "gf_config", *PANEL_RESOURCES),
    "scfoundation_encoder": ("scf_checkpoint", "scf_gene_index"),
    "uce_encoder": ("uce_checkpoint", "uce_esm2", "uce_chrom", "uce_offsets"),
}


def check_resources(resources, required, *, verify_hashes):
    result = {}
    for name in required:
        resource = resources.get(name, {})
        path = Path(resource.get("path", "__absent__"))
        issues = []
        if not path.is_file():
            issues.append("missing_file")
        if not resource.get("sha256"):
            issues.append("hash_not_frozen")
        actual = sha256(path) if verify_hashes and path.is_file() else None
        if actual is not None and actual != resource.get("sha256"):
            issues.append("resource_hash_drift")
        result[name] = {"present": path.is_file(), "bytes": path.stat().st_size if path.is_file() else None,
                        "sha256_verified": actual is not None and actual == resource.get("sha256"),
                        "blockers": issues}
    return result


def require_resources(resources, required):
    checks = check_resources(resources, required, verify_hashes=True)
    if any(row["blockers"] for row in checks.values()):
        raise ContractError("missing or changed model/input resources: " + str(checks))


def resource_path(resources, name):
    return resources[name]["path"]


def validate_adapter_layout(resources, readout):
    """Legacy loaders use fixed basenames; bind every file they will actually open."""
    groups = []
    if readout in READOUTS[:4]:
        groups = [{"gf_tokens": "token_dictionary_gc104M.pkl", "gf_medians": "gene_median_dictionary_gc104M.pkl",
                   "gf_name_ids": "gene_name_id_dict_gc104M.pkl"}]
        if readout == "scgpt_encoder":
            groups.append({"scgpt_vocab": "vocab.json", "scgpt_args": "args.json", "scgpt_checkpoint": "best_model.pt"})
        else:
            groups.append({"scgpt_vocab": "vocab.json"})
        if readout in ("geneformer_embed", "geneformer_attention"):
            groups.append({"gf_config": "config.json", "gf_checkpoint": "model.safetensors"})
        else:
            groups.append({"gf_config": "config.json"})
    for group in groups:
        paths = [Path(resources[name]["path"]).resolve() for name in group]
        if len({p.parent for p in paths}) != 1 or any(p.name != expected for p, expected in zip(paths, group.values())):
            raise ContractError("adapter resource layout differs from explicit loader files")


def tokenizable_genes(resources, *, require_all_models=False):
    required = list(PANEL_RESOURCES)
    if require_all_models:
        required += ["scf_gene_index", "uce_chrom"]
    require_resources(resources, required)
    def pkl(name):
        with Path(resource_path(resources, name)).open("rb") as f:
            return pickle.load(f)
    tokens, medians, names = pkl("gf_tokens"), pkl("gf_medians"), pkl("gf_name_ids")
    vocab = json.loads(Path(resource_path(resources, "scgpt_vocab")).read_text())
    genes = {s for s, e in names.items() if s in vocab and e in tokens and e in medians and float(medians[e]) > 0}
    if require_all_models:
        import pandas as pd
        scf = set(pd.read_csv(resource_path(resources, "scf_gene_index"), sep="\t")["gene_name"])
        chromosomes = pd.read_csv(resource_path(resources, "uce_chrom"))
        uce = set(chromosomes.loc[chromosomes["species"] == "human", "gene_symbol"].astype(str).str.upper())
        genes = {g for g in genes if g in scf and g.upper() in uce}
    return genes


def embedding_graph(embedding, coverage):
    e = np.asarray(embedding, float)
    count = np.asarray(coverage)
    if e.ndim != 2 or count.shape != (len(e),) or not np.isfinite(e).all():
        raise ContractError("invalid gene embedding/coverage")
    observed = count > 0
    graph = np.full((len(e), len(e)), np.nan)
    if observed.any():
        centered = e[observed] - e[observed].mean(0, keepdims=True)
        norms = np.linalg.norm(centered, axis=1)
        valid = norms > 1e-12
        positions = np.flatnonzero(observed)[valid]
        unit = centered[valid] / norms[valid, None]
        graph[np.ix_(positions, positions)] = np.abs(unit @ unit.T)
    return graph


def run_readout(readout, genes, cp10k, raw, resources, *, device, batch, seed, gpu_policy,
                inference_context="frozen_panel", check=None):
    if readout not in READOUTS or inference_context != "frozen_panel":
        raise ContractError("readout/context not registered; full-context is a distinct future estimand")
    if device not in {"cpu", "cuda"}:
        raise ContractError("explicit CPU/CUDA device required")
    started = time.monotonic()
    before = require_gpu(gpu_policy) if device == "cuda" else None
    def check_batch():
        if check is not None:
            check()
        if device == "cuda":
            require_gpu(gpu_policy)
    required = READOUT_RESOURCES[readout]
    require_resources(resources, required)
    validate_adapter_layout(resources, readout)
    cp10k, raw = sparse.csr_matrix(cp10k), sparse.csr_matrix(raw)
    if cp10k.shape != raw.shape or cp10k.shape[1] != len(genes):
        raise ContractError("readout panel/count shape mismatch")
    import torch
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(4)
    if device == "cuda":
        torch.cuda.set_device(int(gpu_policy.get("index", 0)))
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.cuda.reset_peak_memory_stats()
    cells = np.arange(cp10k.shape[0])
    if readout in READOUTS[:4]:
        module = importlib.import_module(".vendor.fm_readout", __package__)
        module.DEV = device
        module.READINESS_CHECK = check_batch
        module.GFM = str(Path(resource_path(resources, "gf_config")).parent)
        module.DCT = str(Path(resource_path(resources, "gf_tokens")).parent)
        module.CK = str(Path(resource_path(resources, "scgpt_vocab")).parent)
        adapter = module.FMReadout(genes, batch=batch)
        if readout == "geneformer_attention":
            graph = adapter.geneformer_attention(cp10k, cells)
            if not np.isfinite(graph).all():
                raise ContractError("nonfinite raw attention output")
            coverage = np.asarray((cp10k > 0).sum(0)).ravel()
            graph[coverage == 0, :] = np.nan
            graph[:, coverage == 0] = np.nan
            graph[adapter.last_pair_counts == 0] = np.nan
        else:
            if readout == "scgpt_encoder":
                log = cp10k.copy()
                log.data = np.log1p(log.data)
                embedding = adapter.scgpt(log, cells)
            else:
                embedding = adapter.geneformer(cp10k, cells, random_init=readout == "random_init_geneformer", seed=seed)
            coverage = adapter.last_gene_counts
            graph = embedding_graph(embedding, coverage)
    elif readout == "scfoundation_encoder":
        module = importlib.import_module(".vendor.fm_readout_scf", __package__)
        module.DEV, module.CKPT, module.GIDX = device, resource_path(resources, "scf_checkpoint"), resource_path(resources, "scf_gene_index")
        module.READINESS_CHECK = check_batch
        adapter = module.SCFReadout(genes, batch=batch)
        graph = embedding_graph(adapter.gene_embed(cp10k, cells), adapter.last_gene_counts)
        coverage = adapter.last_gene_counts
    else:
        module = importlib.import_module(".vendor.fm_readout_uce", __package__)
        module.DEV = device
        module.READINESS_CHECK = check_batch
        module.CKPT, module.ESM2_HUMAN = resource_path(resources, "uce_checkpoint"), resource_path(resources, "uce_esm2")
        module.CHROM_CSV, module.OFFSETS_PKL = resource_path(resources, "uce_chrom"), resource_path(resources, "uce_offsets")
        adapter = module.UCEReadout(genes, batch=batch, seed=seed)
        # Preserve original gene-sentence universe; only frozen-panel input genes carry counts.
        lookup = {g: i for i, g in enumerate(adapter.uni_syms)}
        mapping = [lookup.get(g.upper(), -1) for g in genes]
        values = raw.tocoo()
        columns = np.asarray(mapping)[values.col]
        ok = columns >= 0
        universe = sparse.csr_matrix((values.data[ok], (values.row[ok], columns[ok])), shape=(len(cells), len(lookup)))
        embedding = adapter.gene_embed(universe, cells, None)
        coverage = adapter.last_gene_counts
        graph = embedding_graph(embedding, coverage)
    if np.isinf(graph).any() or not np.isfinite(graph).any():
        raise ContractError("undefined model graph")
    del adapter
    after = None
    peak = None
    if device == "cuda":
        torch.cuda.synchronize()
        peak = int(torch.cuda.max_memory_allocated())
        torch.cuda.empty_cache()
        after = require_gpu(gpu_policy)
    return graph, {"cells": int(cp10k.shape[0]), "batches": int(np.ceil(cp10k.shape[0] / batch)),
                   "peak_memory_bytes": peak, "elapsed_seconds": time.monotonic() - started,
                   "gpu_before": before, "gpu_after": after,
                   "numeric_dtype": "float16 model, float32 readout" if device == "cuda" and readout in READOUTS[:3] else "float32 model/readout",
                   ** {"readout": readout, "coverage_per_gene": np.asarray(coverage).astype(int).tolist(),
                   "unobserved_genes": [g for g, c in zip(genes, coverage) if c == 0],
                   "context": inference_context, "device": device,
                   "resources_sha256": digest({n: resources[n]["sha256"] for n in required}),
                   "scope": "new-panel, same encoder/readout definitions; not official end-to-end regulatory prediction"}}
