#!/usr/bin/env python3
"""Independent, local, read-only source-correction evidence verification.

No producer imports, network, model forwards, motif scans, random draws or GPU
queries. A receipt certifies file identity, not the truth of its metadata: source
selection, counts, normalization, coverage and fixed families are replayed here.
Nothing is written unless --out is requested. Exit 0 means the requested stage
is verified (possibly statistics untestable), 1 partial, and 2 invalid evidence.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import pickle
import re
import sys
from urllib.parse import urlsplit
from zipfile import BadZipFile

import h5py
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

from verify_execution import (HASH, InvalidEvidence, PANEL_RESOURCES,
                              READOUT_RESOURCES, READOUT_SCOPE, ROW_SCOPE,
                              _column, _text, close, correction, digest,
                              load_json, number, require, sha256, unique,
                              verify_receipt)


HERE = Path(__file__).resolve().parent
ARMS = {"brain_selected_panel": ["brain", "pbmc"],
        "historical_panel_corrected_source": ["brain"]}
READOUTS = list(READOUT_RESOURCES)
PREPARED_FILES = ("raw.npz", "cp10k.npz", "graphs.npz", "graph.json")
PREPARED_SCOPE = "metadata-restricted source-corrected motif-accessibility proxy; no causal reference or same-cell brain pairing"
ANALYSIS_SCOPE = "Post-review corrective fixed-panel conditional analysis; no independent donor CI, causal ground truth, gate validity or pretraining exclusion"
PAIRING = "same biological sample identifiers; separate nuclei, no cell pairing or verified donor crosswalk"
FEATURE_NAMING = "exact official name order with standard uniqueness replay, no suffix stripping"
EXCLUSION_REASON = "incomplete_model_or_baseline_frozen_pair_coverage; no zero imputation"
OFFICIAL_FILES = {
    "brain_rna": "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
    "brain_rna_metadata": "GSE174367_snRNA-seq_cell_meta.csv.gz",
    "brain_atac_metadata": "GSE174367_snATAC-seq_cell_meta.csv.gz",
}
FM_ADAPTER_SHA256 = "01b3f951b0e22a66234ce4b3bf37696a69c0d3f116eb61d371fd147e989814e3"
UCE_ADAPTER_SHA256 = "9df70d6f1e161be91cb2579d6a1bdc99e36256317376fefe6c5695fee780c439"
UCE_TOKENS_PER_CELL = 1024
PEAK = re.compile(r"^([^:\s-]+)[:-](\d+)-(\d+)$")


def integer(value, label, low=0, high=math.inf):
    require(type(value) is int and low <= value <= high, f"invalid integer {label}")
    return value


def array_close(actual, expected, label, *, atol=1e-10, rtol=1e-9):
    a, b = np.asarray(actual), np.asarray(expected)
    require(a.shape == b.shape and np.issubdtype(a.dtype, np.number)
            and not np.isinf(a).any() and np.array_equal(np.isnan(a), np.isnan(b))
            and np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=True),
            f"incorrect {label}")


def sparse_close(actual, expected, label, *, exact=False):
    require(sparse.issparse(actual) and actual.shape == expected.shape
            and np.isfinite(actual.data).all(), f"invalid {label} shape/values")
    delta = (actual.astype(float) - expected).tocsr()
    delta.eliminate_zeros()
    require(not delta.nnz or not exact and float(np.max(np.abs(delta.data))) <= 1e-8,
            f"incorrect {label}; source row/feature alignment or normalization differs")


def core_context(plan, resources, code_root=HERE / "screg_validation"):
    """Reconstruct core + runner + ops + corrective/provenance source digest."""
    root = Path(code_root)
    require(root.is_dir(), "core source directory missing")
    files = sorted(root.rglob("*.py"))
    require(bool(files), "core Python source absent")
    files += [root.parent / "run.py", root.parent / "operational_contract.json"]
    core = digest({p.relative_to(root.parent).as_posix(): sha256(p) for p in files})
    code = {"core": core, "source_correction.py": sha256(root.parent / "source_correction.py"),
            "provenance_fingerprint.py": sha256(root.parent / "provenance_fingerprint.py")}
    require(plan.get("resources_sha256") == digest(resources), "plan/resources registry digest drift")
    historical = root.parent / "plan.json"
    require(sha256(historical) == plan.get("historical_plan_sha256"), "historical scientific plan bytes drift")
    return {"plan": digest(plan), "resources": digest(resources), "code": digest(code)}, code


def contract(plan):
    require(plan.get("schema") == "screg-source-correction.v1", "unsupported corrective schema")
    integer(plan.get("seed"), "seed")
    require(plan.get("registered_before_corrective_model_outcomes") is True,
            "corrective preregistration declaration absent")
    require(isinstance(plan.get("arms"), dict) and set(plan["arms"]) == set(ARMS), "registered two arms differ")
    for arm, tissues in ARMS.items():
        require(plan["arms"][arm].get("tissues") == tissues, "registered arm/tissue design differs")
    models, stats, pilot = plan["models"], plan["statistics"], plan["pilot"]
    require(models.get("readouts") == READOUTS and models.get("registered_tasks") == 18,
            "registered six readouts / eighteen tasks differ")
    require(models.get("excluded_historical_readouts") == ["geneformer_ko_raw", "geneformer_ko_posctrl"],
            "withdrawn historical knockout scope differs")
    require(models.get("device") == "cuda" and models.get("inference_context") == "frozen_panel"
            and models.get("batch") == 2, "registered inference device/context/batch differs")
    require(models.get("min_gene_coverage") == .9 and models.get("pair_coverage_required") == 1.,
            "registered coverage gates differ")
    require(stats.get("degree_modes") == ["own", "none"]
            and stats.get("null_modes") == ["restricted", "unrestricted"]
            and stats.get("registered_hypotheses") == 72 and stats.get("donor_intervals") is False,
            "registered seventy-two-slot family/donor scope differs")
    integer(stats.get("permutations"), "registered draws", 1)
    require(stats.get("backend") in {"cpu", "cuda"}, "explicit statistics backend required")
    require(pilot.get("first_tissue") == "brain" and pilot.get("max_batches_per_readout") == 1
            and pilot.get("finite_mask_equal") is True, "two-panel pilot contract differs")
    for key in ("cpu_cuda_atol", "cpu_cuda_rtol"):
        number(pilot.get(key), key, 1e-8, .02)
    integer(pilot.get("min_comparable_pairs"), "pilot comparable pairs", 1)
    require(plan["panel"].get("all_model_intersection") is False,
            "panel cannot be selected on all-model outcome coverage")
    for key in ("max_genes", "min_genes", "min_tfs"):
        integer(plan["panel"].get(key), key, 1)
    require(plan["panel"]["min_genes"] <= plan["panel"]["max_genes"], "invalid panel limits")
    number(plan["panel"].get("detection_floor"), "panel detection threshold", 0, 1)
    integer(plan["panel"].get("promoter_bp"), "promoter interval")
    require(plan["brain"].get("exclude_unannotated") is True
            and plan["brain"].get("sample_key_rna") == "SampleID"
            and plan["brain"].get("sample_key_atac") == "Sample.ID"
            and plan["brain"].get("cell_type_key") == "Cell.Type"
            and plan["brain"].get("reference") == "GRCh38.p12.premrna2"
            and plan["brain"].get("rna_scale") == "integer_counts", "brain metadata/source contract differs")
    integer(plan["brain"].get("max_cells_per_sample"), "sample cell cap", 1)
    integer(plan["pbmc"].get("max_cells"), "PBMC cell cap", 1)
    require(plan["pbmc"].get("rna_scale") == "counts", "PBMC raw-count scale required")
    require(plan["proxy"].get("background") == [.25] * 4
            and plan["proxy"].get("self_edges") == 0, "proxy background/self-edge contract differs")
    integer(plan["proxy"].get("peak_width_bp"), "peak sequence width", 2)
    number(plan["proxy"].get("motif_pvalue"), "motif threshold", 1e-12, 1-1e-12)
    require(isinstance(plan.get("limits"), str) and bool(plan["limits"]), "scientific scope limits absent")


def receipt(path, ctx, outputs):
    """Reject symlinked execution artifacts before the reusable strict checker."""
    path = Path(path)
    for relative in [path.name, *outputs]:
        p = path.parent / relative
        require(not p.is_symlink(), f"symlinked execution evidence: {relative}")
        for parent in p.parents:
            if parent == path.parent:
                break
            require(not parent.is_symlink(), f"symlinked execution directory: {relative}")
    return verify_receipt(path, ctx, outputs)


def csv_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames))
                and {"Barcode", "Cell.Type"} <= set(reader.fieldnames), "metadata CSV columns invalid")
        records = list(reader)
    require(bool(records) and all(None not in row and all(isinstance(v, str) for v in row.values())
                                  for row in records), "malformed metadata CSV row")
    unique([row["Barcode"] for row in records], "metadata barcode")
    return records


def unique_symbols(names):
    """Exact make_index_unique(join='-'), preserving existing suffixed symbols."""
    require(bool(names) and all(isinstance(n, str) and n.strip() == n and n for n in names),
            "missing official feature symbol")
    occupied, seen, counters, result = set(names), set(), {}, []
    for name in names:
        if name not in seen:
            seen.add(name)
            result.append(name)
            continue
        suffix = counters.get(name, 0) + 1
        while f"{name}-{suffix}" in occupied:
            suffix += 1
        counters[name] = suffix
        value = f"{name}-{suffix}"
        occupied.add(value)
        result.append(value)
    return result


def h5ad_info(path):
    with h5py.File(path, "r") as handle:
        cells = unique(_column(handle, "obs", _text(handle["obs"].attrs.get("_index", "_index"))), "source cell")
        genes = unique(_column(handle, "var", _text(handle["var"].attrs.get("_index", "_index"))), "source feature")
        node = handle.get("X")
        require(node is not None, "source X matrix absent")
        shape = node.shape if isinstance(node, h5py.Dataset) else node.attrs.get("shape")
        require(shape is not None and tuple(shape) == (len(cells), len(genes)), "H5AD dimensions disagree with indices")
        genome = _column(handle, "var", "genome") if "genome" in handle["var"] else None
    return cells, genes, genome


def csr_rows(node, rows, shape):
    ptr = np.asarray(node["indptr"][:])
    require(np.issubdtype(ptr.dtype, np.integer) and ptr.shape == (shape[0] + 1,)
            and ptr[0] == 0 and np.all(np.diff(ptr) >= 0)
            and ptr[-1] == len(node["data"]) == len(node["indices"]), "invalid source sparse pointers")
    data, indices, selected_ptr = [], [], [0]
    for row in rows:
        lo, hi = int(ptr[row]), int(ptr[row+1])
        ix, values = np.asarray(node["indices"][lo:hi]), np.asarray(node["data"][lo:hi])
        require(np.issubdtype(ix.dtype, np.integer) and np.all((ix >= 0) & (ix < shape[1]))
                and len(np.unique(ix)) == len(ix), "invalid/duplicate source sparse feature index")
        data.append(values)
        indices.append(ix)
        selected_ptr.append(selected_ptr[-1] + len(ix))
    return sparse.csr_matrix((np.concatenate(data), np.concatenate(indices), selected_ptr),
                             shape=(len(rows), shape[1]), dtype=float)


def read_h5ad_rows(path, rows):
    with h5py.File(path, "r") as handle:
        node = handle["X"]
        shape = tuple(node.shape if isinstance(node, h5py.Dataset) else node.attrs["shape"])
        require(len(rows) == len(set(rows)) and bool(rows)
                and min(rows) >= 0 and max(rows) < shape[0], "invalid source cell selection")
        if isinstance(node, h5py.Dataset):
            require(len(rows) * shape[1] * node.dtype.itemsize <= 512 * 1024**2,
                    "selected dense source exceeds bounded memory")
            order = np.argsort(rows)
            return sparse.csr_matrix(node[np.asarray(rows)[order], :][np.argsort(order)], dtype=float)
        encoding = _text(node.attrs.get("encoding-type", ""))
        require(encoding in {"csr_matrix", "csc_matrix"}, "unknown sparse source encoding")
        if encoding == "csr_matrix":
            return csr_rows(node, rows, shape)
        # Rare CSC inputs are read by feature, keeping only selected cells.
        ptr = np.asarray(node["indptr"][:])
        require(np.issubdtype(ptr.dtype, np.integer) and ptr.shape == (shape[1]+1,)
                and ptr[0] == 0 and np.all(np.diff(ptr) >= 0)
                and ptr[-1] == len(node["data"]) == len(node["indices"]), "invalid CSC pointers")
        lookup = {r: i for i, r in enumerate(rows)}
        rr, cc, dd = [], [], []
        for col in range(shape[1]):
            lo, hi = int(ptr[col]), int(ptr[col+1])
            ix, values = np.asarray(node["indices"][lo:hi]), node["data"][lo:hi]
            require(np.issubdtype(ix.dtype, np.integer) and np.all((ix >= 0) & (ix < shape[0]))
                    and len(np.unique(ix)) == len(ix), "invalid/duplicate CSC cell index")
            for row, value in zip(ix, values):
                if int(row) in lookup:
                    rr.append(lookup[int(row)]); cc.append(col); dd.append(value)
        return sparse.csr_matrix((dd, (rr, cc)), shape=(len(rows), shape[1]), dtype=float)


def normalize(raw, label):
    require(np.isfinite(raw.data).all() and np.all(raw.data >= 0)
            and np.allclose(raw.data, np.round(raw.data), atol=1e-6, rtol=0), f"noninteger/negative {label} counts")
    sums = np.asarray(raw.sum(axis=1)).ravel()
    require(np.isfinite(sums).all() and np.all(sums > 0), f"zero/invalid {label} full-feature library")
    cp = raw.multiply((10000. / sums)[:, None]).tocsr()
    log = cp.copy()
    log.data = np.log1p(log.data)
    return cp, log, sums


def subset(values, universe, genes):
    lookup = {g: i for i, g in enumerate(universe)}
    mapping = np.array([lookup.get(g, -1) for g in genes])
    present = np.flatnonzero(mapping >= 0)
    source = values[:, mapping[present]].tocoo()
    return sparse.csr_matrix((source.data, (source.row, present[source.col])), shape=(len(values.indptr)-1, len(genes)))


class SafeDictionary(pickle.Unpickler):
    def find_class(self, module, name):
        # The pinned Geneformer medians contain NumPy numeric scalars. Permit
        # only their data reconstruction; object/structured dtypes and all
        # other globals remain forbidden.
        if module == "numpy" and name == "dtype":
            def numeric_dtype(*args):
                value = np.dtype(*args)
                require(value.kind in "biuf" and value.fields is None and value.itemsize <= 16,
                        "nonprimitive pickle dtype")
                return value
            return numeric_dtype
        if module in {"numpy.core.multiarray", "numpy._core.multiarray"} and name == "scalar":
            def numeric_scalar(dtype, payload):
                require(isinstance(dtype, np.dtype) and dtype.kind in "biuf" and dtype.fields is None
                        and isinstance(payload, bytes) and len(payload) == dtype.itemsize,
                        "nonprimitive pickle scalar")
                return np.frombuffer(payload, dtype=dtype, count=1)[0].item()
            return numeric_scalar
        raise InvalidEvidence(f"executable/global pickle forbidden in panel dictionary: {module}.{name}")


def panel_dictionaries(resources):
    result = []
    for name in ("gf_tokens", "gf_medians", "gf_name_ids"):
        with Path(resources[name]["path"]).open("rb") as handle:
            value = SafeDictionary(handle).load()
        require(isinstance(value, dict), "panel pickle is not a primitive dictionary")
        result.append(value)
    tokens, medians, names = result
    vocab = load_json(resources["scgpt_vocab"]["path"])
    require(isinstance(vocab, dict), "scGPT vocabulary is not a dictionary")
    allowed = set()
    for name, identifier in names.items():
        require(isinstance(name, str) and isinstance(identifier, str), "invalid gene name/ID dictionary")
        if name in vocab and identifier in tokens and identifier in medians:
            require(type(tokens[identifier]) is int and tokens[identifier] >= 0, "invalid model token")
            if number(medians[identifier], "Geneformer median", 0) > 0:
                allowed.add(name)
    return allowed


def coordinates(path, promoter):
    result = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        require(len(fields) == 5, "malformed coordinate TSV")
        chrom, lo, hi, strand, gene = fields
        lo, hi = int(lo), int(hi)
        require(chrom and gene and 0 <= lo < hi and strand in {"+", "-"}, "invalid locus coordinate")
        if gene not in result:
            result[gene] = (chrom, max(0, lo-promoter) if strand == "+" else lo,
                            hi if strand == "+" else hi+promoter, hi-lo)
    return result


def motif_symbols(path):
    result, motif_ids = set(), set()
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        if not line.startswith("MOTIF "):
            continue
        parts = line.split()
        require(len(parts) >= 3 and parts[1] not in motif_ids, "invalid/duplicate MEME motif identity")
        motif_ids.add(parts[1])
        j = i + 1
        while j < len(lines) and not lines[j].startswith("letter-probability matrix"):
            require(not lines[j].startswith("MOTIF "), "missing MEME PWM")
            j += 1
        require(j < len(lines), "missing MEME PWM")
        width = re.search(r"\bw\s*=\s*(\d+)", lines[j])
        require(width is not None, "MEME width absent")
        pwm = np.array([[float(v) for v in row.split()] for row in lines[j+1:j+1+int(width[1])]])
        require(pwm.shape == (int(width[1]), 4) and np.isfinite(pwm).all()
                and np.all(pwm >= 0) and np.allclose(pwm.sum(1), 1, atol=2e-5), "invalid MEME probabilities")
        result.update(re.sub(r"\(var\.\d+\)", "", parts[2]).upper().split("::"))
    require(bool(motif_ids), "no pinned MEME motifs")
    return result


def peak(value):
    match = PEAK.fullmatch(value)
    require(match is not None and int(match[2]) < int(match[3]), f"invalid peak {value}")
    return match[1], int(match[2]), int(match[3])


def links_for(peaks, genes, loci):
    by_chrom = {}
    for i, value in enumerate(peaks):
        chrom, start, end = peak(value)
        by_chrom.setdefault(chrom, []).append(((start+end)//2, i))
    ordered = {}
    for chrom, values in by_chrom.items():
        entries = sorted(values)
        ordered[chrom] = (np.asarray([v[0] for v in entries]), np.asarray([v[1] for v in entries], dtype=int))
    rr, cc = [], []
    for j, gene in enumerate(genes):
        require(gene in loci, f"panel gene lacks coordinate: {gene}")
        chrom, lo, hi, _ = loci[gene]
        positions, indices = ordered.get(chrom, (np.array([]), np.array([], dtype=int)))
        left, right = np.searchsorted(positions, [lo, hi])
        selected = indices[left:right].tolist()
        rr.extend(selected); cc.extend([j] * len(selected))
    return sparse.csr_matrix((np.ones(len(rr)), (rr, cc)), shape=(len(peaks), len(genes)))


def sequence_gc(fasta, fai, windows):
    """Read only registered windows through the existing FAI; never build an index."""
    entries = {}
    for line in Path(fai).read_text().splitlines():
        name, length, offset, bases, width, *_ = line.split("\t")
        require(name not in entries, "duplicate FASTA index contig")
        entries[name] = tuple(map(int, (length, offset, bases, width)))
    values = []
    with Path(fasta).open("rb") as handle:
        for window in windows:
            chrom, lo, hi = peak(window)
            require(chrom in entries, "peak chromosome absent from pinned FASTA")
            length, offset, bases, width = entries[chrom]
            require(0 <= lo < hi <= length and 0 < bases <= width and offset >= 0, "invalid FAI/window bounds")
            first = offset + (lo//bases)*width + lo % bases
            last = offset + ((hi-1)//bases)*width + (hi-1) % bases + 1
            handle.seek(first)
            seq = handle.read(last-first).replace(b"\n", b"").replace(b"\r", b"").upper()
            require(len(seq) == hi-lo and set(seq) <= set(b"ACGTRYSWKMBDHVN"), "FAI sequence replay mismatch")
            denominator = sum(seq.count(c) for c in (b"A", b"C", b"G", b"T"))
            values.append((seq.count(b"G")+seq.count(b"C"))/denominator if denominator else 0.)
    return np.asarray(values)


class CorrectionVerifier:
    def __init__(self, run_dir, plan, resources, *, code_root=HERE / "screg_validation"):
        contract(plan)
        self.root, self.plan, self.resources = Path(run_dir).resolve(), plan, resources
        self.ctx, self.code = core_context(plan, resources, code_root)
        self.code_root = Path(code_root)
        self.jobs = [arm+"-"+t for arm, rule in plan["arms"].items() for t in rule["tissues"]]
        self.stages, self.errors, self.panels, self.data = {}, [], {}, {}
        self.pilot_rows, self.inference_rows, self.stat_rows = {}, {}, {}
        self.arrays = {}
        self.input_hashes, self.artifact_hashes = {}, {}
        self.receipt_contexts, self.input_stats = {}, {}

    def resources_check(self):
        required = {"brain_rna", "brain_atac", "brain_rna_metadata", "brain_atac_metadata",
                    "pbmc_rna", "pbmc_atac", "pbmc_metadata", "gene_coordinates", "genome_fasta",
                    "genome_fai", "meme", "original_panel", *PANEL_RESOURCES}
        required.update(itertools.chain.from_iterable(READOUT_RESOURCES.values()))
        require(required <= set(self.resources), "incomplete corrective resource registry")
        for name, record in self.resources.items():
            require(isinstance(record, dict) and HASH.fullmatch(record.get("sha256", "")), f"unfrozen resource {name}")
            p = Path(record["path"])
            integer(record.get("bytes"), name+" bytes", 1)
            require(p.is_file() and p.stat().st_size == record["bytes"] and sha256(p) == record["sha256"],
                    f"resource missing/length/hash drift: {name}")
            self.input_hashes[p] = record["sha256"]
            stat = p.stat()
            self.input_stats[p] = (stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)
        downloads = []
        for name, filename in OFFICIAL_FILES.items():
            record = self.resources[name]
            require(record.get("source") == "official_GEO_download", "official RNA/metadata provenance absent")
            side = load_json(str(record["path"])+".download.json")
            head = side.get("head", {})
            require(side.get("status") == "downloaded_verified" and side.get("sha256") == record["sha256"]
                    and side.get("bytes") == record["bytes"] and head.get("length") == record["bytes"],
                    f"official download sidecar mismatch: {name}")
            url = urlsplit(head.get("url", ""))
            require(url.scheme == "https" and url.hostname == "ftp.ncbi.nlm.nih.gov"
                    and url.path == "/geo/series/GSE174nnn/GSE174367/suppl/"+filename,
                    f"official source URL mismatch: {name}")
            require(HASH.fullmatch(side.get("authorization_sha256", "")) is not None
                    and HASH.fullmatch(side.get("resource_sha256", "")) is not None
                    and side["resource_sha256"] == head.get("resource_sha256")
                    and bool(head.get("etag") or head.get("last_modified")), "unpinned download provenance sidecar")
            if filename.endswith(".gz"):
                with gzip.open(record["path"], "rb") as handle:
                    while handle.read(2**20):
                        pass
            downloads.append(name)
        return {"state": "verified", "verified_files": len(self.resources),
                "official_local_downloads_verified": downloads, "network_used": False,
                "download_authorization_scope_replayed": False}

    def source_check(self):
        paths = {n: r["path"] for n, r in self.resources.items()}
        rmeta, ameta = csv_records(paths["brain_rna_metadata"]), csv_records(paths["brain_atac_metadata"])
        with h5py.File(paths["brain_rna"], "r") as handle:
            m = handle["matrix"]
            bc = unique(m["barcodes"].asstr()[:].tolist(), "official brain barcode")
            names = m["features/name"].asstr()[:].tolist()
            ids = unique(m["features/id"].asstr()[:].tolist(), "official feature ID")
            shape = np.asarray(m["shape"][:])
            require(shape.shape == (2,) and np.issubdtype(shape.dtype, np.integer)
                    and tuple(shape) == (len(names), len(bc)) and len(ids) == len(names), "official RNA dimensions/features differ")
            require(set(m["features/feature_type"].asstr()[:]) == {"Gene Expression"}
                    and len(m["features/feature_type"]) == len(names), "non-Gene-Expression official RNA feature")
            genome = sorted(set(m["features/genome"].asstr()[:]))
            require(genome == [self.plan["brain"]["reference"]]
                    and len(m["features/genome"]) == len(names), "official brain genome identity differs")
        allgenes = unique_symbols(names)
        ac, peaks, genome_atac = h5ad_info(paths["brain_atac"])
        require(genome_atac is not None and len(genome_atac) == len(peaks)
                and all("GRCh38" in g or "hg38" in g for g in genome_atac), "brain ATAC reference absent/different")
        require({r["Barcode"] for r in rmeta} <= set(bc)
                and {r["Barcode"] for r in ameta} <= set(ac), "complete assay-specific metadata barcode mismatch")
        rule, seed = self.plan["brain"], self.plan["seed"]
        exclude = {"", "na", "nan", "none", "unknown", "null"}
        pools = {}
        for assay, records, sk in (("rna", rmeta, rule["sample_key_rna"]), ("atac", ameta, rule["sample_key_atac"])):
            require(all(sk in r for r in records), "sample metadata column absent")
            pools[assay] = [r for r in records if r["Cell.Type"].strip().casefold() not in exclude]
            require(all(r[sk].strip() for r in pools[assay]), "blank annotated sample ID")
        samples = sorted({r[rule["sample_key_rna"]] for r in pools["rna"]}
                         & {r[rule["sample_key_atac"]] for r in pools["atac"]})
        require(len(samples) >= 3, "fewer than three shared sample blocks")
        grouped = {}
        for assay, sk in (("rna", rule["sample_key_rna"]), ("atac", rule["sample_key_atac"])):
            grouped[assay] = {}
            for row in pools[assay]:
                grouped[assay].setdefault((row[sk], row["Cell.Type"]), []).append(row)
        types = sorted(set.intersection(*[{label for sample, label in grouped[a] if sample == s}
                                         for a in ("rna", "atac") for s in samples]))
        require(bool(types), "no common annotated sample/assay types")
        quota = {t: min(rule["max_cells_per_sample"]//len(types),
                        min(len(grouped[a][(s, t)]) for a in grouped for s in samples)) for t in types}
        require(min(quota.values()) >= 1, "empty per-type quota")
        selection = {}
        for assay in ("rna", "atac"):
            selection[assay] = []
            for sample in samples:
                for label in types:
                    chosen = sorted(grouped[assay][sample, label], key=lambda row:
                                    hashlib.sha256(f'{seed}:{assay}:{sample}:{label}:{row["Barcode"]}'.encode()).hexdigest())[:quota[label]]
                    selection[assay].extend(chosen)
        composition = {"samples": samples, "types": types, "per_type_quota": quota,
                       "cells_per_sample": sum(quota.values()), "sample_count": len(samples), "pairing": PAIRING}
        rlookup, alookup = {c: i for i, c in enumerate(bc)}, {c: i for i, c in enumerate(ac)}
        rrows = [rlookup[r["Barcode"]] for r in selection["rna"]]
        arows = [alookup[r["Barcode"]] for r in selection["atac"]]
        with h5py.File(paths["brain_rna"], "r") as handle:
            raw = csr_rows(handle["matrix"], rrows, (len(bc), len(names)))
        cp, log, library = normalize(raw, "official brain full-feature")
        self.data["brain"] = {"raw": raw, "cp": cp, "log": log, "universe": allgenes,
                              "library": library, "peaks": peaks, "atac_rows": arows,
                              "labels": [r["Cell.Type"] for r in selection["atac"]]}
        pbc, pgenes, _ = h5ad_info(paths["pbmc_rna"])
        pac, ppeaks, _ = h5ad_info(paths["pbmc_atac"])
        require(pbc == pac, "PBMC RNA/ATAC complete cell order differs")
        pmeta = csv_records(paths["pbmc_metadata"])
        pbc_set = set(pbc)
        valid = [r for r in pmeta if r["Barcode"] in pbc_set and r["Cell.Type"] not in {"", "NA", "nan", "NaN"}]
        chosen = sorted(valid, key=lambda r: hashlib.sha256(f'{seed}:pbmc:{r["Barcode"]}'.encode()).hexdigest())[:self.plan["pbmc"]["max_cells"]]
        require(bool(chosen), "no annotated paired PBMC cells")
        lookup = {c: i for i, c in enumerate(pbc)}
        prows = [lookup[r["Barcode"]] for r in chosen]
        praw = read_h5ad_rows(paths["pbmc_rna"], prows)
        pcp, plog, plibrary = normalize(praw, "PBMC full-feature")
        self.data["pbmc"] = {"raw": praw, "cp": pcp, "log": plog, "universe": pgenes,
                             "library": plibrary, "peaks": ppeaks, "atac_rows": prows,
                             "labels": [r["Cell.Type"] for r in chosen]}
        self.loci = coordinates(paths["gene_coordinates"], self.plan["panel"]["promoter_bp"])
        tokenizable = panel_dictionaries(self.resources)
        candidates = sorted(set(self.loci) & tokenizable)
        opportunity = links_for(peaks, candidates, self.loci).getnnz(0)
        detection = dict(zip(allgenes, np.asarray((raw > 0).mean(0)).ravel()))
        passing = [g for g, count in zip(candidates, opportunity)
                   if count > 0 and detection.get(g, 0) >= self.plan["panel"]["detection_floor"]]
        tfs = sorted(set(passing) & motif_symbols(paths["meme"]))
        require(len(tfs) <= self.plan["panel"]["max_genes"], "TF count exceeds panel cap")
        non_tfs = sorted(set(passing)-set(tfs), key=lambda g: (-detection[g], g))
        genes = sorted(tfs+non_tfs[:self.plan["panel"]["max_genes"]-len(tfs)])
        require(len(genes) >= self.plan["panel"]["min_genes"] and len(tfs) >= self.plan["panel"]["min_tfs"],
                "genuine brain-selected panel below registered support")
        original = load_json(paths["original_panel"])
        historical = unique(original["genes"], "historical panel gene")
        require(set(historical) <= set(original["tf_flags"])
                and all(type(original["tf_flags"][g]) is bool for g in historical), "historical TF flags invalid")
        htfs = unique([g for g in historical if original["tf_flags"][g]], "historical TF")
        for arm, gs, ts in (("brain_selected_panel", genes, tfs), ("historical_panel_corrected_source", historical, htfs)):
            self.panels[arm] = {"genes": gs, "tf_genes": ts, "panel_sha256": digest({"genes": gs, "tf_genes": ts})}
        self.selection = {"brain": composition, "rna_cells": selection["rna"], "atac_cells": selection["atac"],
                          "brain_rna_metadata_missing_cells": len(bc)-len(rmeta), "brain_atac_metadata_missing_cells": len(ac)-len(ameta),
                          "new_panel_selection_scale": "detection over the registered metadata-balanced genuine brain RNA cell pool; no model outcomes",
                          "historical_overlap": len(set(genes) & set(historical)),
                          "withdrawn_readouts": self.plan["models"]["excluded_historical_readouts"],
                          "source": {"features": len(names), "cells": len(bc), "feature_ids_sha256": digest(ids),
                                     "genome": genome, "feature_names_sha256": digest(names), "feature_naming": FEATURE_NAMING},
                          "brain_rna_sha256": self.resources["brain_rna"]["sha256"], "brain_atac_sha256": self.resources["brain_atac"]["sha256"],
                          "pbmc_cells": chosen}
        return {"state": "verified", "official_brain": self.selection["source"], "brain": composition,
                "brain_rna_selected_cells": len(rrows), "brain_atac_selected_cells": len(arows),
                "unannotated_rna_cells_excluded": len(bc)-len(rmeta), "unannotated_atac_cells_excluded": len(ac)-len(ameta),
                "pbmc_selected_cells": len(prows), "pbmc_original_types_absent_from_selection": sorted({r["Cell.Type"] for r in valid}-{r["Cell.Type"] for r in chosen}),
                "panels": {arm: {"genes": len(p["genes"]), "tfs": len(p["tf_genes"]), "sha256": p["panel_sha256"]} for arm, p in self.panels.items()},
                "historical_panel_overlap": self.selection["historical_overlap"], "metadata_selection_replayed": True,
                "normalization_universe": {t: len(d["universe"]) for t, d in self.data.items()},
                "upstream_local_ATAC_and_PBMC_processing_chain_verified": False,
                "panel_selection_independent_or_held_out": False, "brain_same_cell_pairing": False}

    def prepared_check(self):
        expected = ["panels.json", "selection.json", "prepared.json"] + [j+"/"+f for j in self.jobs for f in PREPARED_FILES]
        if not (self.root / "prepared.receipt.json").is_file():
            details = {"state": "partial" if self.root.exists() and any(self.root.iterdir()) else "not_started",
                       "missing_receipt": "prepared.receipt.json", "present_outputs": [p for p in expected if (self.root/p).is_file()]}
            # Stable selection/panel JSON can be reviewed during preparation;
            # it never upgrades the incomplete stage to prepared.
            for filename, wanted in (("panels.json", self.panels), ("selection.json", self.selection)):
                p = self.root/filename
                if p.is_file():
                    require(load_json(p) == wanted, f"partial {filename} differs from independent source replay")
            return details
        self.check_receipt(self.root / "prepared.receipt.json", self.ctx, expected)
        require(load_json(self.root/"panels.json") == self.panels, "panel selection/history differs from independent replay")
        require(load_json(self.root/"selection.json") == self.selection, "metadata/feature selection differs from independent replay")
        report = load_json(self.root/"prepared.json")
        require(report.get("status") == "prepared" and report.get("jobs") == self.jobs
                and report.get("source_context") == self.ctx
                and report.get("panels_sha256") == sha256(self.root/"panels.json")
                and report.get("selection_sha256") == sha256(self.root/"selection.json")
                and report.get("claims") == "Separate corrective analysis; historical biological scope withdrawn, historical outputs untouched",
                "prepared report identity/scope differs")
        summaries = []
        for job in self.jobs:
            arm, tissue = job.rsplit("-", 1)
            panel, data = self.panels[arm], self.data[tissue]
            genes, tfs = panel["genes"], panel["tf_genes"]
            n = len(genes)
            raw, cp, log = [subset(data[k], data["universe"], genes) for k in ("raw", "cp", "log")]
            sparse_close(sparse.load_npz(self.root/job/"raw.npz"), raw, job+" raw counts", exact=True)
            sparse_close(sparse.load_npz(self.root/job/"cp10k.npz"), cp, job+" full-library CP10k")
            with np.load(self.root/job/"graphs.npz", allow_pickle=False) as archive:
                require(set(archive.files) == {"proxy", "coexpression", "rna_variance", "peak_count", "gene_length", "gc", "detection"},
                        "prepared graph keys differ")
                a = {k: archive[k] for k in archive.files}
            self.arrays[job] = a
            proxy = a["proxy"]
            require(proxy.shape == (n,n) and np.isfinite(proxy).all() and np.all(proxy >= 0)
                    and np.all(np.diag(proxy) == 0), "prepared proxy shape/finiteness/self edges invalid")
            tfpos = [genes.index(g) for g in tfs]
            require(np.all(proxy[~np.isin(np.arange(n), tfpos)] == 0), "non-TF proxy rows contain hits")
            means = np.asarray(log.mean(0)).ravel()
            variance = np.asarray(log.multiply(log).mean(0)).ravel()-means**2
            covariance = (log.T@log).toarray()/log.shape[0]-np.outer(means, means)
            std = np.sqrt(np.maximum(0, np.diag(covariance)))
            den = np.outer(std, std)
            co = np.divide(covariance, den, out=np.zeros_like(covariance), where=den > 0)
            co = np.clip(np.abs(co), 0, 1)
            co[variance <= 0, :] = np.nan; co[:, variance <= 0] = np.nan
            array_close(a["rna_variance"], variance, "RNA variance")
            array_close(a["coexpression"], co, "full-library log1p coexpression")
            link = links_for(data["peaks"], genes, self.loci)
            keep = np.flatnonzero(link.getnnz(1))
            linked = link[keep].tocsr()
            kept = [data["peaks"][i] for i in keep]
            windows = []
            for value in kept:
                chrom, lo, hi = peak(value)
                mid, half = (lo+hi)//2, self.plan["proxy"]["peak_width_bp"]//2
                windows.append(f"{chrom}:{max(0, mid-half)}-{mid+half}")
            gc = sequence_gc(self.resources["genome_fasta"]["path"], self.resources["genome_fai"]["path"], windows)
            peak_counts = linked.getnnz(0)
            expected_gc = np.asarray(linked.T@gc).ravel()/np.maximum(1, peak_counts)
            expected_detection = np.asarray((raw > 0).mean(0)).ravel()
            array_close(a["peak_count"], peak_counts, "strand-aware half-open peak counts", rtol=0, atol=0)
            array_close(a["gene_length"], [self.loci[g][3] for g in genes], "gene lengths", rtol=0, atol=0)
            array_close(a["gc"], expected_gc, "pinned FASTA linked-peak sequence GC")
            array_close(a["detection"], expected_detection, "tissue RNA detection")
            atac = read_h5ad_rows(self.resources[tissue+"_atac"]["path"], data["atac_rows"])[:, keep].tocsr()
            require(np.isfinite(atac.data).all() and np.all(atac.data >= 0)
                    and np.allclose(atac.data, np.round(atac.data), atol=1e-6, rtol=0), "ATAC values are not raw integer counts")
            labels = np.asarray(data["labels"])
            types = sorted(set(data["labels"]))
            access = np.mean([np.log1p(np.asarray(atac[labels == t].mean(0)).ravel()) for t in types], axis=0)
            bound = np.asarray(linked.T@access).ravel()
            require(np.all(proxy <= bound[None, :]+1e-9), "motif proxy exceeds all-hit equal-type accessibility bound")
            positions = {g: i for i, g in enumerate(data["universe"])}
            expected_report = {"arm": arm, "tissue": tissue, "cells": len(data["atac_rows"]), "types": types,
                               "rna_cells": raw.shape[0], "panel_sha256": panel["panel_sha256"],
                               "peak_ids_sha256": digest(kept), "window_ids_sha256": digest(windows), "linked_peaks": len(kept),
                               "motif_parameters": self.plan["proxy"], "ordered_tf_sha256": digest(tfs),
                               "source_feature_position": [positions.get(g) for g in genes],
                               "missing_rna_features": [g for g in genes if g not in positions],
                               "undefined_coexpression_genes": [g for g,v in zip(genes, variance) if v <= 0], "scope": PREPARED_SCOPE}
            require(load_json(self.root/job/"graph.json") == expected_report, "prepared graph/source metadata differs")
            ii, jj = self.pairs(panel)
            summaries.append({"job": job, "cells": raw.shape[0], "genes": n, "tfs": len(tfs), "linked_peaks": len(kept),
                              "missing_rna_features": expected_report["missing_rna_features"],
                              "undefined_coexpression_genes": expected_report["undefined_coexpression_genes"],
                              "baseline_all_frozen_pairs_finite": bool(np.isfinite(co[ii,jj]).all()),
                              "proxy_nonconstant_frozen_pairs": bool(np.ptp(proxy[ii,jj]) > 0),
                              "full_library_sizes_min_max": [float(data["library"].min()), float(data["library"].max())]})
        return {"state": "verified", "outputs_verified": len(expected), "jobs": summaries,
                "counts_and_full_feature_normalization_replayed": True,
                "motif_hit_scan_reexecuted": False,
                "motif_proxy_verified_scope": "hash-bound output, TF/self-edge rules, linked peaks, source sequence GC and accessibility upper bound; individual motif hits were not independently rescanned"}

    def check_receipt(self, path, ctx, outputs):
        record = receipt(path, ctx, outputs)
        self.artifact_hashes[Path(path)] = sha256(path)
        self.receipt_contexts[Path(path)] = record["context_sha256"]
        for name, value in record["outputs"].items():
            self.artifact_hashes[Path(path).parent/name] = value
        return record

    @staticmethod
    def pairs(panel):
        n = len(panel["genes"])
        ii = np.repeat([panel["genes"].index(g) for g in panel["tf_genes"]], n)
        jj = np.tile(np.arange(n), len(panel["tf_genes"]))
        return ii[ii != jj], jj[ii != jj]

    def graph_file(self, path, panel):
        with np.load(path, allow_pickle=False) as archive:
            require(set(archive.files) == {"graph"}, "model graph archive keys differ")
            g = archive["graph"]
        n = len(panel["genes"])
        require(g.shape == (n,n) and np.issubdtype(g.dtype, np.number)
                and not np.isinf(g).any() and np.isfinite(g).any(), "model graph shape/finiteness differs")
        require(np.all(g[np.isfinite(g)] >= -1e-9) and np.all(g[np.isfinite(g)] <= 1+1e-6),
                "model affinity graph outside [0,1]")
        return g

    def telemetry(self, row, readout, panel, cells, device, graph, *, source_raw=None):
        require(row.get("readout") == readout and row.get("cells") == cells
                and type(row.get("cells")) is int and row.get("batches") == math.ceil(cells/self.plan["models"]["batch"])
                and type(row.get("batches")) is int and row.get("device") == device
                and row.get("context") == "frozen_panel", "readout telemetry execution identity differs")
        require(row.get("resources_sha256") == digest({n: self.resources[n]["sha256"] for n in READOUT_RESOURCES[readout]})
                and row.get("scope") == READOUT_SCOPE, "readout resource identity/scope differs")
        counts, genes = row.get("coverage_per_gene"), panel["genes"]
        require(isinstance(counts, list) and len(counts) == len(genes)
                and all(type(c) is int and c >= 0 for c in counts), "noninteger/negative gene coverage")
        # UCE samples gene tokens with replacement; its observation count is
        # not a number of distinct cells. Other adapters count observed cells.
        if readout == "uce_encoder":
            adapter = self.code_root/"vendor/fm_readout_uce.py"
            source_receipt = load_json(self.code_root/"vendor/source_receipt.json")
            require(source_receipt.get("files",{}).get("fm_readout_uce.py",{}).get("adapted_sha256") == UCE_ADAPTER_SHA256
                    and sha256(adapter) == UCE_ADAPTER_SHA256, "UCE token multiplicity requires exact independently pinned adapter")
            require(max(counts) <= cells*UCE_TOKENS_PER_CELL, "UCE token count exceeds hash-bound 1024 draws per selected cell")
        else:
            require(max(counts) <= cells, "cell-observation coverage exceeds selected cells")
        if source_raw is not None:
            expressed = np.asarray((source_raw > 0).sum(0)).ravel().astype(int)
            require(expressed.shape == (len(genes),), "coverage source panel dimensions differ")
            require(all(c == 0 for c,e in zip(counts,expressed) if e == 0), "tokens claim source-zero gene observations")
            if readout in READOUTS[:4]:
                adapter = self.code_root/"vendor/fm_readout.py"
                source_receipt = load_json(self.code_root/"vendor/source_receipt.json")
                require(source_receipt.get("files",{}).get("fm_readout.py",{}).get("adapted_sha256") == FM_ADAPTER_SHA256
                        and sha256(adapter) == FM_ADAPTER_SHA256, "positive-cell coverage requires pinned FM adapter")
                require(counts == expressed.tolist(), "deterministic positive-cell token counts differ from source panel")
            elif readout == "uce_encoder":
                require(np.all(np.asarray(counts) <= expressed*UCE_TOKENS_PER_CELL),
                        "UCE counts exceed available source-positive cells times 1024")
        missing = [g for g,c in zip(genes,counts) if c == 0]
        require(row.get("unobserved_genes") == missing, "unobserved genes contradict token counts")
        for i, count in enumerate(counts):
            if count == 0:
                require(np.isnan(graph[i,:]).all() and np.isnan(graph[:,i]).all(), "missing gene zero/finite imputation")
        number(row.get("elapsed_seconds"), "readout elapsed seconds", 1e-12,
               self.plan["models"]["max_seconds"] if cells > 2 else self.plan["pilot"]["max_seconds"])
        require(isinstance(row.get("numeric_dtype"), str) and bool(row["numeric_dtype"]), "numeric dtype telemetry absent")
        if device == "cpu":
            require(row.get("gpu_before") is None and row.get("gpu_after") is None
                    and row.get("peak_memory_bytes") is None, "CPU readout claims CUDA telemetry")
        else:
            integer(row.get("peak_memory_bytes"), "CUDA peak memory", 1)
            for key in ("gpu_before", "gpu_after"):
                gpu = row.get(key)
                require(isinstance(gpu, dict) and gpu.get("ready") is True
                        and gpu.get("blockers") == [], "CUDA readiness telemetry absent/blocked")
        ii,jj = self.pairs(panel)
        return {"gene_coverage_fraction": 1-len(missing)/len(genes),
                "evaluated_pair_coverage_fraction": float(np.isfinite(graph[ii,jj]).mean())}

    def pilot_check(self):
        directory = self.root/"pilot"
        expected, missing = ["pilot.json"], []
        for arm, readout in itertools.product(self.plan["arms"], READOUTS):
            names = ["cuda.npz", "cpu.npz", "pilot.json"]
            expected.extend(f"{arm}/{readout}/{f}" for f in [*names,"result.receipt.json"])
            td = directory/arm/readout
            if not (td/"result.receipt.json").is_file():
                missing.append({"job": arm+"-brain", "readout": readout,
                                "state": "partial" if td.exists() else "not_started"})
                continue
            panel = self.panels[arm]
            context = {**self.ctx,"job": arm+"-brain", "readout": readout, "panel_sha256": panel["panel_sha256"]}
            self.check_receipt(td/"result.receipt.json", context, names)
            row = load_json(td/"pilot.json")
            require(row.get("status") == "pilot_success" and row.get("job") == arm+"-brain"
                    and row.get("readout") == readout and row.get("panel_sha256") == panel["panel_sha256"]
                    and row.get("scope") == "engineering pilot, not biological result", "pilot identity/status/scope differs")
            cuda, cpu = self.graph_file(td/"cuda.npz", panel), self.graph_file(td/"cpu.npz", panel)
            cells = min(self.plan["models"]["batch"], self.data["brain"]["raw"].shape[0])
            raw = subset(self.data["brain"]["raw"],self.data["brain"]["universe"],panel["genes"])[:cells]
            for device, graph in (("cpu", cpu), ("cuda", cuda)):
                self.telemetry(row[device], readout, panel, cells, device, graph,source_raw=raw)
            require(row["cpu"]["coverage_per_gene"] == row["cuda"]["coverage_per_gene"],
                    "same seeded CPU/CUDA pilot token coverage differs")
            require(np.array_equal(np.isfinite(cuda),np.isfinite(cpu)), "pilot CPU/CUDA finite masks differ")
            mask = np.isfinite(cuda) & ~np.eye(len(cuda),dtype=bool)
            count = int(mask.sum())
            require(count >= self.plan["pilot"]["min_comparable_pairs"], "pilot lacks comparable nonself pairs")
            require(np.allclose(cuda[mask],cpu[mask],atol=self.plan["pilot"]["cpu_cuda_atol"],
                                rtol=self.plan["pilot"]["cpu_cuda_rtol"]), "actual saved CPU/CUDA pilot numeric comparison fails")
            integer(row.get("pairs"), "pilot pair count", 1)
            require(row["pairs"] == count, "pilot measured pair count differs")
            close(row.get("max_difference"), float(np.max(np.abs(cuda[mask]-cpu[mask]))), "pilot actual maximum difference")
            self.pilot_rows[arm,readout] = row
        exists = (directory/"pilot.receipt.json").is_file()
        if exists:
            require(not missing, "global pilot receipt claims incomplete children")
            self.check_receipt(directory/"pilot.receipt.json", self.ctx, expected)
            report = load_json(directory/"pilot.json")
            require(report.get("status") == "pilot_complete" and report.get("source_context") == self.ctx
                    and report.get("readouts") == list(self.pilot_rows.values()), "global pilot exact twelve rows differ")
        return {"state": "verified" if exists else "partial" if directory.exists() else "not_started",
                "expected_children": 12, "verified_children": len(self.pilot_rows), "missing": missing,
                "global_receipt_verified": exists, "global_outputs": len(expected) if exists else 0,
                "actual_saved_cpu_cuda_comparison_recomputed": True, "scientific_results": False,
                "per_gene_coverage_bound":{"uce_encoder":"source-positive cells * 1024 sampled tokens, with replacement",
                                           "others":"selected cells; positive-cell counts replayed for Geneformer/scGPT"},
                "uce_sampling_adapter_sha256":UCE_ADAPTER_SHA256}

    def inference_check(self):
        require(self.stages["pilot"]["state"] == "verified", "inference cannot be certified without global two-panel pilot")
        expected, missing, eligibility = ["inference.json"], [], []
        for job, readout in itertools.product(self.jobs,READOUTS):
            directory, panel = self.root/job/readout, self.panels[job.rsplit("-",1)[0]]
            names = ["graph.npz","coverage.json"]
            expected.extend(f"{job}/{readout}/{f}" for f in [*names,"result.receipt.json"])
            if not (directory/"result.receipt.json").is_file():
                missing.append({"job": job,"readout": readout,"state": "partial" if directory.exists() else "not_started"})
                continue
            self.check_receipt(directory/"result.receipt.json", {**self.ctx,"job":job,"readout":readout,"panel":panel["panel_sha256"]}, names)
            row = load_json(directory/"coverage.json")
            require(row.get("job") == job and row.get("panel_sha256") == panel["panel_sha256"], "inference job/panel differs")
            graph = self.graph_file(directory/"graph.npz",panel)
            cells = self.data[job.rsplit("-",1)[1]]["raw"].shape[0]
            data = self.data[job.rsplit("-",1)[1]]
            raw = subset(data["raw"],data["universe"],panel["genes"])
            fractions = self.telemetry(row,readout,panel,cells,"cuda",graph,source_raw=raw)
            for key,value in fractions.items():
                close(row.get(key),value,key)
            passed = fractions["gene_coverage_fraction"] >= self.plan["models"]["min_gene_coverage"]
            eligible = passed and fractions["evaluated_pair_coverage_fraction"] == 1.
            require(row.get("coverage_threshold_pass") is passed and row.get("primary_statistics_eligible") is eligible,
                    "inference eligibility contradicts actual frozen pairs")
            self.inference_rows[job,readout] = row
            eligibility.append({"job":job,"readout":readout,**fractions,"eligible":eligible})
        for p in self.root.glob("*/*/coverage.json"):
            require((p.parent.parent.name,p.parent.name) in set(itertools.product(self.jobs,READOUTS)), "unregistered inference output")
        exists = (self.root/"inference.receipt.json").is_file()
        if exists:
            require(not missing, "global inference claims incomplete eighteen children")
            self.check_receipt(self.root/"inference.receipt.json",self.ctx,expected)
            report = load_json(self.root/"inference.json")
            require(report.get("status") == "inference_complete" and report.get("source_context") == self.ctx
                    and type(report.get("tasks")) is int and report["tasks"] == 18
                    and report.get("rows") == list(self.inference_rows.values())
                    and report.get("eligible") == sum(r["primary_statistics_eligible"] for r in self.inference_rows.values()),
                    "global inference task/row/eligible accounting differs")
        return {"state":"verified" if exists else "partial" if self.inference_rows or missing and any(v["state"] == "partial" for v in missing) else "not_started",
                "expected_tasks":18,"verified_tasks":len(self.inference_rows),"missing":missing,
                "global_receipt_verified":exists,"global_outputs":len(expected) if exists else 0,
                "eligible_tasks":sum(r["eligible"] for r in eligibility),"coverage":eligibility,
                "inference_json_sha256":self.artifact_hashes.get(self.root/"inference.json"),
                "inference_receipt_sha256":self.artifact_hashes.get(self.root/"inference.receipt.json"),
                "context_sha256":digest(self.ctx)}

    def strata_sizes(self,job,null):
        panel,a = self.panels[job.rsplit("-",1)[0]],self.arrays[job]
        if null == "unrestricted":
            return [len(panel["genes"])]
        op = np.digitize(a["peak_count"],np.unique(np.quantile(a["peak_count"],[.2,.4,.6,.8])))
        gc = np.digitize(a["gc"],np.unique(np.quantile(a["gc"],[1/3,2/3])))
        tf = np.isin(np.arange(len(panel["genes"])),[panel["genes"].index(g) for g in panel["tf_genes"]])
        labels = 15*tf+3*op+gc
        return [int(np.sum(labels == l)) for l in np.unique(labels)]

    @staticmethod
    def partial(x,y,controls):
        require(np.isfinite(x).all() and np.isfinite(y).all() and np.ptp(x) > 0 and np.ptp(y) > 0,
                "constant/undefined measured statistic")
        x,y = rankdata(x),rankdata(y)
        x -= x.mean(); y -= y.mean()
        design = np.column_stack([np.ones(len(x)),controls])
        u,s,_ = np.linalg.svd(design,full_matrices=False)
        q = u[:,s > np.finfo(float).eps*max(design.shape)*s[0]]
        xr,yr = x-q@(q.T@x),y-q@(q.T@y)
        nx,ny = np.linalg.norm(xr),np.linalg.norm(yr)
        tolerance = np.finfo(float).eps*max(design.shape)*8
        require(nx > tolerance*np.linalg.norm(x) and ny > tolerance*np.linalg.norm(y), "fully conditioned measured statistic")
        return float(np.clip(np.dot(xr,yr)/(nx*ny),-1,1))

    def checkpoint(self,directory,ctx,names):
        p = directory/"draw-progress.json"
        if not p.is_file():
            return None
        before = sha256(p)
        state = load_json(p)
        body = {k:v for k,v in state.items() if k != "state_sha256"}
        require(state.get("state_sha256") == digest(body) and state.get("context_sha256") == digest(ctx), "draw checkpoint hash/context differs")
        done = integer(state.get("next_draw"),"checkpoint draw position",0,self.plan["statistics"]["permutations"])
        require(isinstance(state.get("nulls"),dict) and set(state["nulls"]) == {"M","D"}, "checkpoint component accounting differs")
        for kind in ("M","D"):
            values = np.asarray(state["nulls"][kind],float)
            require(values.size == done*len(names) and (done == 0 or values.shape == (done,len(names)))
                    and np.isfinite(values).all() and np.all(np.abs(values) <= 1+1e-9), "invalid saved checkpoint null values")
        rng = state.get("rng_state")
        require(isinstance(rng,dict) and rng.get("bit_generator") == "PCG64", "checkpoint RNG identity missing")
        bitgen = np.random.PCG64(0)
        bitgen.state = rng  # Validate serialized state; no random values drawn.
        if done == self.plan["statistics"]["permutations"]:
            require(sha256(p) == before, "completed draw checkpoint changed during read")
            self.artifact_hashes[p] = before
        return state

    def measured_check(self,row,job,readout,degree,null,state,index):
        rule,panel,a = self.plan["statistics"],self.panels[job.rsplit("-",1)[0]],self.arrays[job]
        require(tuple(row.get(k) for k in ("job","readout","degree_mode","null_mode")) == (job,readout,degree,null), "statistics slot identity differs")
        require(row.get("measured") is True and row.get("permutations") == rule["permutations"]
                and type(row.get("permutations")) is int and row.get("backend") == rule["backend"]
                and row.get("scope") == ROW_SCOPE, "statistics measurement/draw/backend/scope differs")
        require(row.get("stratum_sizes") == self.strata_sizes(job,null) and max(row["stratum_sizes"]) > 1,
                "registered randomization strata differ")
        for key in ("observed","comparable_fm","comparable_baseline","concordance"):
            number(row.get(key),key,-1,1)
        number(row.get("comparable_difference"),"comparable difference",-2,2)
        require(state is not None and state["next_draw"] == rule["permutations"], "measured statistics lacks complete actual saved draw checkpoint")
        denominator = rule["permutations"]+1
        for key in ("pM","pD"):
            p = number(row.get(key),key,1/denominator,1)
            close(p*denominator,round(p*denominator),"plus-one Monte Carlo p grid")
            null_values = np.asarray(state["nulls"][key[-1]],float)[:,index]
            actual = (1+int(np.sum(np.abs(null_values) >= abs(row["observed"]))))/denominator
            close(p,actual,"actual saved draws "+key)
        close(row.get("pIUT"),max(row["pM"],row["pD"]),"IUT component maximum")
        close(row.get("correction_input_p"),row["pIUT"],"measured correction input")
        ii,jj = self.pairs(panel)
        proxy,co = a["proxy"][ii,jj],a["coexpression"][ii,jj]
        graph = self.graph_file(self.root/job/readout/"graph.npz",panel)[ii,jj]
        controls = np.column_stack([a[k][jj] for k in ("peak_count","gene_length","detection","gc")])
        if degree == "own":
            controls = np.column_stack([controls,(a["proxy"]>0).sum(1)[ii],(a["proxy"]>0).sum(0)[jj]])
        for key,value in (("observed",self.partial(graph,proxy,np.column_stack([controls,rankdata(co)]))),
                          ("comparable_fm",self.partial(graph,proxy,controls)),
                          ("comparable_baseline",self.partial(co,proxy,controls))):
            require(math.isclose(row[key],value,rel_tol=1e-6,abs_tol=1e-8), f"independent measured {key} differs")
        close(row["comparable_difference"],row["comparable_fm"]-row["comparable_baseline"],"paired comparable difference")
        close(row["concordance"],float(np.corrcoef(rankdata(graph),rankdata(co))[0,1]),"FM/baseline concordance")
        require(row.get("ci95") is None and row.get("donor_intervals") is None, "unregistered donor interval")

    def statistics_check(self):
        require(self.stages["infer"]["state"] == "verified", "statistics cannot be certified without eighteen-task global inference")
        groups,missing,untestable = [],[],[]
        for job in self.jobs:
            panel = self.panels[job.rsplit("-",1)[0]]
            ii,jj = self.pairs(panel)
            baseline_ok = bool(np.isfinite(self.arrays[job]["coexpression"][ii,jj]).all())
            names = [r for r in READOUTS if self.inference_rows[job,r]["primary_statistics_eligible"] and baseline_ok]
            for degree,null in itertools.product(self.plan["statistics"]["degree_modes"],self.plan["statistics"]["null_modes"]):
                directory = self.root/job/("statistics-"+degree+"-"+null)
                ctx = {**self.ctx,"job":job,"degree":degree,"null":null,"readouts":names}
                if not names:
                    require(not (directory/"result.receipt.json").exists(), "untestable specification claims measured receipt")
                    state = None
                    group_state = "untestable"
                else:
                    state = self.checkpoint(directory,ctx,names)
                    if (directory/"result.receipt.json").is_file():
                        self.check_receipt(directory/"result.receipt.json",ctx,["rows.json"])
                        rows = load_json(directory/"rows.json")
                        require(isinstance(rows,list) and [r.get("readout") for r in rows] == names, "statistics names duplicated/missing/reordered")
                        for i,(readout,row) in enumerate(zip(names,rows)):
                            self.measured_check(row,job,readout,degree,null,state,i)
                            self.stat_rows[job,readout,degree,null] = row
                        group_state = "verified"
                    else:
                        missing.append({"job":job,"degree":degree,"null":null,"readouts":names})
                        group_state = "partial" if directory.exists() else "not_started"
                groups.append({"job":job,"degree":degree,"null":null,"state":group_state,
                               "directory":directory.relative_to(self.root).as_posix(),
                               "context":ctx,"context_sha256":digest(ctx),
                               "measured_readouts":names,"draws_saved":state["next_draw"] if state else 0,
                               "rows_sha256":self.artifact_hashes.get(directory/"rows.json"),
                               "result_receipt_sha256":self.artifact_hashes.get(directory/"result.receipt.json"),
                               "completed_checkpoint_sha256":self.artifact_hashes.get(directory/"draw-progress.json")})
                for readout in READOUTS:
                    if readout not in names:
                        c = self.inference_rows[job,readout]
                        untestable.append((job,readout,degree,null))
                        self.stat_rows[job,readout,degree,null] = {"job":job,"readout":readout,"degree_mode":degree,"null_mode":null,
                            "measured":False,"observed":None,"pM":None,"pD":None,"pIUT":None,"permutations":0,"correction_input_p":1.,
                            "reason":EXCLUSION_REASON,"coverage":{k:c[k] for k in ("gene_coverage_fraction","evaluated_pair_coverage_fraction")}}
        return {"state":"partial" if missing else "verified","registered_slots":72,
                "measured_slots_verified":sum(r["measured"] for r in self.stat_rows.values()),
                "untestable_slots":len(untestable),"pending_slots":72-len(self.stat_rows),
                "groups":groups,"missing":missing,"null_randomizations_reexecuted":False,
                "p_values_recomputed_from_saved_draws":True,
                "checkpoint_rng_draw_sequence_replayed":False,
                "untestable_correction_input_p_role":"administrative p=1 family placeholder; actual MC p, effects and CI absent; zero draws",
                "pending_is_not_untestable":True}

    def analyze_check(self):
        if not (self.root/"analysis.receipt.json").is_file():
            return {"state":"partial" if (self.root/"results.json").exists() or self.stages["statistics"]["state"] == "verified" else "not_started",
                    "missing_receipt":"analysis.receipt.json"}
        require(self.stages["statistics"]["state"] == "verified", "complete analysis claims pending statistics")
        self.check_receipt(self.root/"analysis.receipt.json",self.ctx,["results.json"])
        report = load_json(self.root/"results.json")
        rows = report.get("rows")
        require(isinstance(rows,list) and len(rows) == 72, "analysis registered seventy-two rows incomplete")
        identities = [tuple(r.get(k) for k in ("job","readout","degree_mode","null_mode")) for r in rows]
        require(len(set(identities)) == 72 and set(identities) == set(self.stat_rows), "analysis exact Cartesian family differs")
        by,holm = correction([r["correction_input_p"] for r in rows],"BY"),correction([r["correction_input_p"] for r in rows],"Holm")
        for row,key,q,p in zip(rows,identities,by,holm):
            body = {k:v for k,v in row.items() if k not in {"q_IUT_BY","p_IUT_Holm"}}
            require(body == self.stat_rows[key], "analysis slot contradicts measured/excluded evidence; administrative p=1 is not measured MC p")
            close(row.get("q_IUT_BY"),q,"pooled 72-slot BY")
            close(row.get("p_IUT_Holm"),p,"pooled 72-slot Holm")
        measured = sum(r["measured"] for r in rows)
        require(report.get("status") == "complete" and report.get("source_context") == self.ctx
                and report.get("registered_hypotheses") == 72 and report.get("measured_hypotheses") == measured
                and report.get("limits") == self.plan["limits"] and report.get("historical_results_replaced") is False
                and report.get("donor_intervals") is None and report.get("scope") == ANALYSIS_SCOPE,
                "analysis terminal count/context/scientific scope differs")
        groups = {(g["job"],g["degree"],g["null"]):g for g in self.stages["statistics"]["groups"]}
        row_bindings = []
        for index,(row,key) in enumerate(zip(rows,identities)):
            job,readout,degree,null = key
            group = groups[job,degree,null]
            row_bindings.append({"row_index":index,
                "identity":dict(zip(("job","readout","degree_mode","null_mode"),key)),
                "row_sha256":digest(row),"evidence_body_sha256":digest(self.stat_rows[key]),
                "statistics_group":group["directory"],"statistics_context_sha256":group["context_sha256"],
                "statistics_child_row_index":group["measured_readouts"].index(readout) if row["measured"] else None,
                "state":"measured" if row["measured"] else "untestable","draws":row["permutations"],
                "correction_input_p_role":"actual_MC_IUT" if row["measured"] else "administrative_family_placeholder"})
        return {"state":"verified","registered_slots":72,"measured_slots":measured,"untestable_slots":72-measured,
                "statistics_state":"measured" if measured else "untestable","corrections_recomputed":True,
                "donor_intervals":None,"historical_results_replaced":False,
                "evidence_binding":{"results_json_sha256":self.artifact_hashes[self.root/"results.json"],
                    "analysis_receipt_sha256":self.artifact_hashes[self.root/"analysis.receipt.json"],
                    "context_sha256":digest(self.ctx),"ordered_row_identity_sha256":digest(identities),
                    "ordered_rows_sha256":digest(rows),"rows":row_bindings}}

    def budget_check(self):
        records = {}
        for stage in ("prepare","pilot","infer","analyze"):
            path = self.root.parent/(self.root.name+"-"+stage+"-budget.json")
            if not path.is_file():
                records[stage] = {"state":"not_recorded"}
                continue
            state = load_json(path)
            maximum = self.plan["pilot"]["max_seconds"] if stage == "pilot" else self.plan["models"]["max_seconds"] if stage == "infer" else self.plan["statistics"]["max_seconds"]
            elapsed = number(state.get("elapsed_seconds"),"cumulative stage elapsed",0)
            integer(state.get("invocations"),"stage invocations",1)
            require(state.get("maximum") == maximum and type(state.get("running")) is bool, "cumulative stage budget identity differs")
            records[stage] = {"state":"running" if state["running"] else "finished", "elapsed_seconds":elapsed,
                              "maximum":maximum,"invocations":state["invocations"],"resume_budget_available":not state["running"] and elapsed < maximum}
        return {"state":"verified","records":records,"authorization_inferred":False}

    def verify(self,stage="analyze"):
        order = ["resources","source","prepare","pilot","infer","statistics","analyze"]
        stop = order.index(stage)
        for name in order[:stop+1]:
            if self.errors or name in {"pilot","infer","statistics","analyze"} and any(
                    self.stages.get(n,{}).get("state") != "verified" for n in order[2:order.index(name)]):
                self.stages[name] = {"state":"blocked","reason":"upstream evidence incomplete/invalid"}
                continue
            try:
                fn = self.prepared_check if name == "prepare" else self.inference_check if name == "infer" else getattr(self,name+"_check")
                self.stages[name] = fn()
            except (ValueError,OSError,KeyError,TypeError,EOFError,IndexError,BadZipFile,pickle.UnpicklingError) as exc:
                self.stages[name] = {"state":"invalid","reason":str(exc)}
                self.errors.append({"stage":name,"reason":str(exc)})
        snapshot_stable = True
        try:
            self.stages["budgets"] = self.budget_check()
            # Receipt hashes are rechecked after parsing to detect replacement
            # during a live run. Inputs also need stable file metadata; their
            # hashes were computed before any source replay.
            for path,wanted in self.artifact_hashes.items():
                require(path.is_file() and sha256(path) == wanted, "execution artifacts changed during verification")
            for path,wanted in self.input_stats.items():
                stat = path.stat()
                require((stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns) == wanted,
                        "audited source/resource changed during verification")
            ctx,_ = core_context(self.plan,self.resources,self.code_root)
            require(ctx == self.ctx, "source signature changed during verification")
        except (ValueError,OSError,KeyError,TypeError) as exc:
            snapshot_stable = False
            self.errors.append({"stage":"snapshot","reason":str(exc)})
        state = self.stages[stage]["state"]
        all_complete = stage == "analyze" and state == "verified" and not self.errors
        return {"schema":"screg-independent-correction-verification.v1","run_dir":str(self.root),
                "verifier_sha256":sha256(__file__),"shared_read_only_helpers_sha256":sha256(HERE/"verify_execution.py"),
                "source_context":self.ctx,"code_identity":self.code,"requested_stage":stage,
                "requested_stage_state":"invalid" if self.errors else state,
                "status":"invalid" if self.errors else "complete" if all_complete else "partial",
                "execution_complete":all_complete,"statistics_state":self.stages.get("analyze",{}).get("statistics_state","pending"),
                "stages":self.stages,"errors":self.errors,
                "evidence_hashes":{"snapshot_stable":snapshot_stable,"source_context_sha256":digest(self.ctx),
                    "artifact_sha256":{p.relative_to(self.root).as_posix():h for p,h in sorted(self.artifact_hashes.items())},
                    "receipt_context_sha256":{p.relative_to(self.root).as_posix():h for p,h in sorted(self.receipt_contexts.items())}},
                "read_only":True,"network_used":False,"model_forward_executed":False,
                "scope":"Separate post-review conditional panel/source case study; no causal truth, independent gate validity, donor CI or pretraining exclusion"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--resources", type=Path, required=True)
    parser.add_argument("--stage", choices=("prepare", "pilot", "infer", "analyze"), default="analyze")
    parser.add_argument("--out", type=Path, help="Explicit optional report; never overwrite audited inputs or an existing report")
    args = parser.parse_args(argv)
    try:
        plan, resources = load_json(args.plan), load_json(args.resources)
        verifier = CorrectionVerifier(args.run_dir, plan, resources)
        result = verifier.verify(args.stage)
    except (ValueError, OSError, KeyError, TypeError, EOFError, IndexError, BadZipFile) as exc:
        result = {"status": "invalid", "requested_stage": args.stage, "errors": [{"reason": str(exc)}]}
    exit_code = 2 if result["status"] == "invalid" else 0 if result.get("requested_stage_state") == "verified" else 1
    text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n"
    if args.out:
        target = args.out.resolve()
        protected = {args.plan.resolve(), args.resources.resolve()}
        if isinstance(locals().get("resources"), dict):
            protected.update(Path(r["path"]).resolve() for r in resources.values() if isinstance(r, dict) and "path" in r)
        require(not target.is_relative_to(args.run_dir.resolve()) and target not in protected,
                "report path would modify audited inputs/run")
        require(not args.out.is_symlink() and not args.out.exists(), "refuse report overwrite/symlink")
        with args.out.open("x") as handle:
            handle.write(text)
    sys.stdout.write(text)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
