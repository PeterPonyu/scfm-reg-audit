"""Explicit H5AD metadata and bounded CSR reads, old and new encodings."""
from __future__ import annotations

from pathlib import Path
import re

import h5py
import numpy as np
from scipy import sparse

from .common import ContractError, identities, unique


def _str(item):
    return item.decode("utf-8") if isinstance(item, bytes) else str(item)


def column(f: h5py.File, group: str, key: str) -> list[str]:
    node = f.get(f"{group}/{key}")
    if node is None:
        raise ContractError(f"missing H5AD column: {group}/{key}")
    if isinstance(node, h5py.Group):
        if "codes" not in node or "categories" not in node:
            raise ContractError(f"unsupported H5AD column encoding: {group}/{key}")
        codes, categories = np.asarray(node["codes"][:]), node["categories"][:]
    else:
        values = node[:]
        old = f.get(f"{group}/__categories/{key}")
        if old is None:
            return [_str(v) for v in values]
        codes, categories = np.asarray(values), old[:]
    if not np.issubdtype(codes.dtype, np.integer) or np.any(codes < 0) or np.any(codes >= len(categories)):
        raise ContractError(f"missing/invalid categorical codes: {group}/{key}")
    return [_str(v) for v in categories[codes]]


def index(f, group):
    return unique(column(f, group, _str(f[group].attrs.get("_index", "_index"))), group + " index")


def matrix_shape(f, key):
    node = f.get(key)
    if node is None:
        raise ContractError(f"missing matrix: {key}; no fallback to X")
    shape = node.shape if isinstance(node, h5py.Dataset) else node.attrs.get("shape")
    if shape is None or len(shape) != 2:
        raise ContractError("matrix has no 2D shape")
    return tuple(map(int, shape))


def read_rows(path, matrix, rows, columns=None, *, max_dense_bytes=512 * 1024 ** 2):
    """Read requested rows in caller order. Never materialize the full atlas."""
    rows = np.asarray(rows, dtype=np.int64)
    if rows.ndim != 1 or not len(rows) or len(rows) != len(set(rows)):
        raise ContractError("nonempty unique row indices required")
    with h5py.File(path, "r") as f:
        nr, nc = matrix_shape(f, matrix)
        if rows.min() < 0 or rows.max() >= nr:
            raise ContractError("row index out of bounds")
        if columns is not None:
            columns = np.asarray(columns, dtype=np.int64)
            if columns.ndim != 1 or not len(columns) or len(set(columns)) != len(columns) or columns.min() < 0 or columns.max() >= nc:
                raise ContractError("invalid column indices")
        node = f[matrix]
        if isinstance(node, h5py.Dataset):
            if len(rows) * nc * node.dtype.itemsize > max_dense_bytes:
                raise ContractError("dense read exceeds memory budget")
            order = np.argsort(rows)
            values = node[rows[order], :][np.argsort(order)]
            out = sparse.csr_matrix(values)
        else:
            encoding = _str(node.attrs.get("encoding-type", ""))
            if encoding == "csc_matrix":
                # Source matching may encounter old CSC RNA. Read selected columns
                # one at a time; keep only requested cells, never transpose the atlas.
                wanted_cols = np.arange(nc) if columns is None else columns
                row_lookup = {int(r): i for i, r in enumerate(rows)}
                rr, cc, dd = [], [], []
                for j, col in enumerate(wanted_cols):
                    start, end = map(int, node["indptr"][col:col + 2])
                    values = node["data"][start:end]
                    indices = node["indices"][start:end]
                    for k, row in enumerate(indices):
                        if int(row) in row_lookup:
                            rr.append(row_lookup[int(row)])
                            cc.append(j)
                            dd.append(values[k])
                return sparse.csr_matrix((dd, (rr, cc)), shape=(len(rows), len(wanted_cols)))
            if encoding != "csr_matrix":
                raise ContractError("unsupported sparse matrix encoding")
            sizes, data, indices = [0], [], []
            for row in rows:
                start, end = map(int, node["indptr"][row:row + 2])
                data.append(node["data"][start:end])
                indices.append(node["indices"][start:end])
                sizes.append(sizes[-1] + end - start)
            out = sparse.csr_matrix((np.concatenate(data), np.concatenate(indices), sizes), shape=(len(rows), nc))
    return out if columns is None else out[:, columns].tocsr()


def counts_normalization(values, scale):
    x = sparse.csr_matrix(values, dtype=np.float64)
    x.sum_duplicates()
    x.eliminate_zeros()
    if not np.isfinite(x.data).all() or np.any(x.data < 0):
        raise ContractError("nonfinite or negative RNA values")
    if scale == "counts":
        if not np.allclose(x.data, np.round(x.data), atol=1e-6, rtol=0):
            raise ContractError("declared counts are not integer-valued")
        rs = np.asarray(x.sum(1)).ravel()
        if np.any(rs == 0):
            raise ContractError("zero-library cell; explicitly exclude before sampling")
        cp = x.multiply((10000 / rs)[:, None]).tocsr()
    elif scale == "cp10k":
        rs = np.asarray(x.sum(1)).ravel()
        if not np.allclose(rs, 10000, rtol=.02):
            raise ContractError("declared CP10k must use the full RNA feature universe")
        cp = x
    elif scale == "log1p_cp10k":
        cp = x.copy()
        cp.data = np.expm1(cp.data)
        if not np.isfinite(cp.data).all() or not np.allclose(np.asarray(cp.sum(1)).ravel(), 10000, rtol=.02):
            raise ContractError("declared log1p CP10k scale inconsistent")
    else:
        raise ContractError("explicit scale must be counts, cp10k or log1p_cp10k")
    log = cp.copy()
    log.data = np.log1p(log.data)
    return cp, log


PEAK = re.compile(r"^([^:\s-]+)[:-](\d+)-(\d+)$")


def parse_peak(value):
    m = PEAK.fullmatch(str(value))
    if not m or int(m[2]) >= int(m[3]):
        raise ContractError(f"invalid peak coordinate: {value}")
    return m[1], int(m[2]), int(m[3])


def genome_evidence(f, spec):
    """Verify the declared processed-file reference without inferring it from peaks."""
    evidence = spec.get("genome_evidence")
    if spec.get("genome_build") != "GRCh38" or not evidence:
        raise ContractError("GRCh38 requires explicit provenance; coordinate syntax alone cannot identify build")
    if isinstance(evidence, dict):
        if (evidence.get("kind") != "hash_bound_h5ad_metadata" or evidence.get("resource") != "multiome"
                or evidence.get("expected_value") != spec["genome_build"]):
            raise ContractError("unsupported genome evidence contract")
        node = f.get(evidence.get("metadata_path", ""))
        if not isinstance(node, h5py.Dataset) or node.shape != ():
            raise ContractError("declared genome metadata absent or nonscalar")
        value = _str(node[()])
        if value != evidence["expected_value"]:
            raise ContractError("input genome metadata disagrees with registered reference")
        return {**evidence, "observed_value": value, "metadata_verified": True}
    # Explicit processing records and synthetic fixtures retain their exact text;
    # they are not elevated to independently verified upstream evidence.
    return {"record": str(evidence), "metadata_verified": False}


def inspect(path, spec):
    with h5py.File(path, "r") as f:
        cells, features = index(f, "obs"), index(f, "var")
        shape = matrix_shape(f, spec["matrix"])
        if shape != (len(cells), len(features)):
            raise ContractError("indices/matrix shape mismatch")
        types = column(f, "var", spec["feature_type_key"])
        donors = identities(column(f, "obs", spec["donor_key"]), "donor")
        labels = identities(column(f, "obs", spec["cell_type_key"]), "cell type")
        if not (len(types) == len(features) and len(donors) == len(labels) == len(cells)):
            raise ContractError("metadata lengths do not match matrix")
        if spec.get("modality_key"):
            modality = set(column(f, "obs", spec["modality_key"]))
            if modality != {spec["paired_modality"]}:
                raise ContractError("same-cell modality contract not met")
        peak_names = [n for n, t in zip(features, types) if t == spec["atac_value"]]
        for name in peak_names:
            parse_peak(name)
        if not peak_names or spec["rna_value"] not in types:
            raise ContractError("paired RNA and ATAC features absent")
        reference = genome_evidence(f, spec)
        matrix_node = f[spec["matrix"]]
        return {"cells": cells, "features": features, "feature_types": types,
                "donors": donors, "cell_types": labels, "shape": list(shape),
                "matrix": spec["matrix"], "scale": spec["scale"],
                "matrix_encoding": _str(matrix_node.attrs.get("encoding-type", "dense")),
                "n_donors": len(set(donors)), "n_rna_features": types.count(spec["rna_value"]),
                "n_peaks": len(peak_names), "genome_evidence": reference,
                "paired_scope": "same-row processed multiome; donor identities from explicit metadata"}
