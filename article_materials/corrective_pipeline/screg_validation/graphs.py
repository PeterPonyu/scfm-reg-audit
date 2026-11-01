"""Independently constructed donor graphs with explicit locus and linkage rules."""
from __future__ import annotations

from pathlib import Path
import re

import numpy as np
from scipy import sparse

from .common import ContractError, identities, unique
from .h5ad import parse_peak


def read_coordinates(path, promoter_bp):
    if promoter_bp < 0:
        raise ContractError("negative promoter interval")
    out = {}
    for line in Path(path).read_text().splitlines():
        chrom, start, end, strand, gene = line.split("\t")
        start, end = int(start), int(end)
        if start < 0 or end <= start or strand not in {"+", "-"}:
            raise ContractError("invalid gene coordinate")
        if gene not in out:
            out[gene] = (chrom, max(0, start - promoter_bp) if strand == "+" else start,
                         end if strand == "+" else end + promoter_bp, end - start)
    return out


def locus_links(peaks, genes, coordinates):
    """Peak midpoint within a fixed half-open promoter+body interval, never nearest gene."""
    unique(peaks, "peak")
    unique(genes, "gene")
    by_chrom = {}
    for i, peak in enumerate(peaks):
        chrom, lo, hi = parse_peak(peak)
        by_chrom.setdefault(chrom, []).append((i, (lo + hi) // 2))
    ordered = {}
    for chrom, entries in by_chrom.items():
        pairs = sorted(entries, key=lambda x: x[1])
        ordered[chrom] = (np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]))
    rows, cols = [], []
    for j, gene in enumerate(genes):
        if gene not in coordinates:
            raise ContractError(f"gene has no pinned coordinate: {gene}")
        chrom, lo, hi, _ = coordinates[gene]
        if chrom in ordered:
            indices, midpoints = ordered[chrom]
            selected = indices[np.searchsorted(midpoints, lo):np.searchsorted(midpoints, hi)]
            rows.extend(selected.tolist())
            cols.extend([j] * len(selected))
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(peaks), len(genes)))


def read_meme(path):
    motifs = {}
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        if not line.startswith("MOTIF "):
            continue
        parts = line.split()
        j = i + 1
        while j < len(lines) and not lines[j].startswith("letter-probability matrix"):
            j += 1
        if j == len(lines):
            raise ContractError("missing MEME PWM")
        width = re.search(r"\bw\s*=\s*(\d+)", lines[j])
        if width is None:
            raise ContractError("missing MEME width")
        pwm = np.array([[float(v) for v in row.split()] for row in lines[j + 1:j + 1 + int(width[1])]])
        if pwm.shape != (int(width[1]), 4) or np.any(pwm < 0) or not np.allclose(pwm.sum(1), 1, atol=2e-5):
            raise ContractError("invalid MEME probabilities")
        symbols = [s.upper() for s in re.sub(r"\(var\.\d+\)", "", parts[2]).split("::")]
        motifs[parts[1]] = {"symbols": symbols, "pwm": pwm.T}
    if not motifs:
        raise ContractError("no MEME motifs")
    return motifs


def build_motif_hits(peaks, tf_genes, fasta_path, meme_path, *, pvalue, check=None):
    """Scan new peak sequences with installed MOODS; no FIMO/bedtools dependency.

    Only peaks in frozen panel loci are passed by the runner. Ambiguous windows
    do not acquire hits by replacing N with A. Both strands are scanned.
    """
    import MOODS.scan
    import MOODS.tools
    from pyfaidx import Fasta

    if not 0 < pvalue < 1:
        raise ContractError("invalid motif threshold")
    motifs = read_meme(meme_path)
    tf_index = {g: i for i, g in enumerate(unique(tf_genes, "TF"))}
    matrices, thresholds, owners = [], [], []
    bg = [.25] * 4
    for motif in motifs.values():
        owner = [tf_index[g] for g in motif["symbols"] if g in tf_index]
        if not owner:
            continue
        lo = MOODS.tools.log_odds(motif["pwm"].tolist(), bg, .01)
        threshold = MOODS.tools.threshold_from_p(lo, bg, pvalue)
        matrices += [lo, MOODS.tools.reverse_complement(lo)]
        thresholds += [threshold] * 2
        owners += [owner] * 2
    if not matrices:
        raise ContractError("no motifs mapped to TF panel")
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)
    # Existing FAI only; no index writes into the resource directory.
    fasta = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True,
                  rebuild=False, build_index=False)
    row, col, gc, lengths = [], [], [], []
    try:
        for i, peak in enumerate(peaks):
            if check and i % 256 == 0:
                check()
            chrom, start, end = parse_peak(peak)
            if chrom not in fasta or end > len(fasta[chrom]):
                raise ContractError("peak outside explicit genome reference")
            seq = str(fasta[chrom][start:end]).upper()
            acgt = sum(seq.count(c) for c in "ACGT")
            gc.append((seq.count("G") + seq.count("C")) / acgt if acgt else 0.)
            lengths.append(acgt)
            hit = set()
            for segment in re.findall("[ACGT]+", seq):
                for k, matches in enumerate(scanner.scan(segment)):
                    if matches:
                        hit.update(owners[k])
            row += [i] * len(hit)
            col += sorted(hit)
    finally:
        fasta.close()
    return (sparse.csr_matrix((np.ones(len(row)), (row, col)), shape=(len(peaks), len(tf_genes))),
            np.array(gc), np.array(lengths))


def donor_graph(atac, rna_log, links, motif_hits, tf_indices, *, accessibility_floor,
                linkage="locus", correlation_floor=.1):
    """No test-donor information is consumed when constructing reference-donor degree.

    `locus` uses fixed windows and donor-specific accessible support. The separately
    registered `paired_correlation` mode estimates links anew within each donor;
    its RNA dependence makes it an exploratory reference, not a causal truth label.
    """
    atac = sparse.csr_matrix(atac, dtype=float)
    rna_log = sparse.csr_matrix(rna_log, dtype=float)
    links = sparse.csr_matrix(links, dtype=float)
    motif_hits = sparse.csr_matrix(motif_hits, dtype=float)
    if atac.shape[0] != rna_log.shape[0] or atac.shape[1] != links.shape[0] or links.shape[1] != rna_log.shape[1] or motif_hits.shape != (links.shape[0], len(tf_indices)):
        raise ContractError("donor graph feature/cell alignment mismatch")
    if not 0 < accessibility_floor <= 1 or linkage not in {"locus", "paired_correlation"}:
        raise ContractError("unknown donor graph rule")
    for x in (atac, rna_log, links, motif_hits):
        if not np.isfinite(x.data).all() or np.any(x.data < 0):
            raise ContractError("graph inputs must be finite and nonnegative")
    n = atac.shape[0]
    frequency = np.asarray((atac > 0).mean(0)).ravel()
    support = frequency >= accessibility_floor
    weighted = links.multiply((support * frequency)[:, None]).tocsr()
    if linkage == "paired_correlation":
        if not 0 <= correlation_floor < 1:
            raise ContractError("invalid preregistered link correlation threshold")
        coo = weighted.tocoo()
        weights = np.zeros(coo.nnz)
        x = (atac > 0).astype(float).tocsc()
        y = rna_log.tocsc()
        for k, (p, g) in enumerate(zip(coo.row, coo.col)):
            xp, yg = x[:, p], y[:, g]
            xm, ym = float(xp.sum()) / n, float(yg.sum()) / n
            sx = float(xp.multiply(xp).sum()) / n - xm ** 2
            sy = float(yg.multiply(yg).sum()) / n - ym ** 2
            rho = ((float(xp.multiply(yg).sum()) / n - xm * ym) / np.sqrt(sx * sy)) if sx > 0 and sy > 0 else 0.
            weights[k] = coo.data[k] * rho if rho >= correlation_floor else 0.
        weighted = sparse.csr_matrix((weights, (coo.row, coo.col)), shape=links.shape)
        weighted.eliminate_zeros()
    values = (motif_hits.T @ weighted).toarray()
    graph = np.zeros((rna_log.shape[1], rna_log.shape[1]))
    graph[np.asarray(tf_indices)] = values
    np.fill_diagonal(graph, 0)
    return graph, {"accessible_peak_count": int(support.sum()), "nonzero_links": weighted.nnz,
                   "linkage": linkage, "cells": n,
                   "scope": "motif family and donor-accessibility proxy; no verified causal edge labels"}


def gene_coexpression(log_values):
    x = sparse.csr_matrix(log_values, dtype=float)
    n = x.shape[0]
    mean = np.asarray(x.mean(0)).ravel()
    covariance = (x.T @ x).toarray() / n - mean[:, None] * mean[None, :]
    sd = np.sqrt(np.maximum(0, np.diag(covariance)))
    denominator = sd[:, None] * sd[None, :]
    out = np.divide(covariance, denominator, out=np.zeros_like(covariance), where=denominator > 0)
    return np.clip(np.abs(out), 0, 1)


def independent_degree(reference_graphs, reference_donors, evaluation_donor, panel):
    reference_donors = unique(reference_donors, "reference donor")
    if evaluation_donor in reference_donors or set(reference_donors) != set(panel["split"]["degree_reference"]):
        raise ContractError("degree reference includes evaluation or wrong donor set")
    graphs = [np.asarray(g, float) for g in reference_graphs]
    n = len(panel["genes"])
    if len(graphs) != len(reference_donors) or any(g.shape != (n, n) or not np.isfinite(g).all() for g in graphs):
        raise ContractError("reference graphs do not match frozen panel/donor identity")
    # Mean donor degree, rather than union support that saturates with donor count.
    return (np.mean([(g > 0).sum(1) for g in graphs], axis=0),
            np.mean([(g > 0).sum(0) for g in graphs], axis=0))
