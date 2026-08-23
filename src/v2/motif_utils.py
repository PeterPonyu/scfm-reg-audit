#!/usr/bin/env python
"""Shared motif helpers: parse JASPAR MEME PWMs and scan peak sequences with MOODS."""
import numpy as np

ACGT = "ACGT"


def parse_meme(path):
    """Parse a JASPAR MEME file -> {motif_id: {'name': TF_symbol, 'pwm': (4, w) probs}}.
    Dimer names (A::B) are kept as-is; caller splits into constituent TFs."""
    motifs = {}; mid = None; name = None; rows = []; reading = False; w = 0
    def flush():
        if mid and rows:
            pwm = np.array(rows, float).T                     # (w,4) -> (4,w) rows=A,C,G,T
            motifs[mid] = {"name": name, "pwm": pwm}
    for ln in open(path):
        s = ln.strip()
        if s.startswith("MOTIF"):
            flush(); rows = []; reading = False
            parts = s.split()
            mid = parts[1]; name = parts[2] if len(parts) > 2 else parts[1]
        elif s.startswith("letter-probability matrix"):
            reading = True
            try: w = int(s.split("w=")[1].split()[0])
            except Exception: w = 0
        elif reading:
            vals = s.split()
            if len(vals) == 4:
                try: rows.append([float(x) for x in vals])
                except ValueError: reading = False
            else:
                reading = False
    flush()
    return motifs


def tf_symbols(name):
    """JASPAR motif name -> list of uppercased gene symbols (splits dimers on '::')."""
    return [t.strip().upper() for t in name.replace("(var.2)", "").replace("(var.3)", "").split("::") if t.strip()]


def _logodds(pwm, bg=(0.25, 0.25, 0.25, 0.25), pseudo=0.01):
    import MOODS.tools
    mat = [pwm[i].tolist() for i in range(4)]                  # 4 rows A,C,G,T (probabilities)
    return MOODS.tools.log_odds(mat, list(bg), pseudo)         # MOODS applies the pseudocount


def build_scanner(motifs, p=1e-4, bg=(0.25, 0.25, 0.25, 0.25)):
    """Return (matrices, thresholds, order) incl. reverse-complement matrices.
    order[k] = motif_id for matrix k (fwd and rev share the same id)."""
    import MOODS.tools
    mats, thr, order = [], [], []
    for mid, d in motifs.items():
        lo = _logodds(d["pwm"], bg)
        rc = MOODS.tools.reverse_complement(lo)
        t = MOODS.tools.threshold_from_p(lo, list(bg), p)
        mats.append(lo); thr.append(t); order.append(mid)
        mats.append(rc); thr.append(t); order.append(mid)
    return mats, thr, order


def scan_seq_hits(seq, mats, thr, order):
    """Return set of motif_ids with >=1 hit in seq (either strand)."""
    import MOODS.scan
    seq = seq.upper()
    if not seq or set(seq) - set("ACGTN"):
        seq = "".join(c if c in "ACGT" else "A" for c in seq)
    res = MOODS.scan.scan_dna(seq, mats, list((0.25, 0.25, 0.25, 0.25)), thr, 7)
    hit = set()
    for k, r in enumerate(res):
        if r: hit.add(order[k])
    return hit
