#!/usr/bin/env python
"""
scfm-reg-audit v2 — UCE (Universal Cell Embedding) gene-level readout (4th FM).

Reimplements the official gene-sentence construction (chromosome-grouped, expression-weighted
sampling with replacement — snap-stanford/UCE eval_data.py::sample_cell_sentences) and captures
the per-token encoder output (model.py::TransformerModel.forward already returns it as the first
value, `gene_output` — the official driver discards it via `_`, keeping only the CLS-pooled cell
embedding). Vendored model.py (BSD-licensed, snap-stanford/UCE) in uce_vendor/.

Gene identity -> embedding table row is reconstructed EXACTLY per the official pipeline
(data_proc/data_utils.py::adata_path_to_prot_chrom_starts + generate_idxs):
  pe_row_idx[gene] = spec_pe_genes.index(gene.upper()) + species_offset
where spec_pe_genes = list(ESM2_human_dict.keys()) in FILE INSERTION ORDER (uppercased) — verified
against evaluate.py::run_eval, which shows the checkpoint's own baked-in pe_embedding.weight is what
is actually used at inference (the all_tokens.torch overwrite path is a no-op once vstacked with the
seeded chrom tensors, since the shape then always matches 145469 and the overwrite is skipped).
Chromosome tokens: pandas Categorical codes over the full (multi-species) species_chrom.csv table,
offset by CHROM_TOKEN_OFFSET=143574 (matching get_spec_chrom_csv exactly) — these are themselves
random (torch.manual_seed(23)) vectors per the official code, not gene-specific, so an approximate
reconstruction only needs internal self-consistency, which categorical-code reuse guarantees.

Deliberate deviation from the vendored code: mask is explicitly cast to bool before being passed to
nn.TransformerEncoder (the original repo relies on an implicit float-mask coercion whose semantics
are torch-version-dependent — this repo is pinned to torch==2.1.1, we run torch 2.12).
"""
import os, sys, json, pickle, time, numpy as np, pandas as pd, torch, torch.nn as nn
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "uce_vendor"))
from model import TransformerModel
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
CKPT = f"{DATA_ROOT}/models/UCE/4layer_model.torch"
ESM2_HUMAN = f"{ROOT}/data/uce/human_esm2.pt"
CHROM_CSV = f"{ROOT}/data/uce/species_chrom.csv"
OFFSETS_PKL = f"{ROOT}/data/uce/species_offsets.pkl"
DEV = "cuda" if torch.cuda.is_available() else "cpu"

PAD_TOK, CHROM_L_TOK, CHROM_R_TOK, CLS_TOK = 0, 1, 2, 3
CHROM_TOKEN_OFFSET = 143574
N_ROWS = 145469


class UCEReadout:
    def __init__(self, gsyms, sample_size=1024, pad_length=1536, batch=4, seed=20260713):
        self.gsyms = list(gsyms); self.Ng = len(gsyms); self.batch = batch
        self.sample_size = sample_size; self.pad_length = pad_length
        self.rng = np.random.default_rng(seed)

        # ---- reconstruct spec_pe_genes (insertion order, uppercased) + offset ----
        esm2 = torch.load(ESM2_HUMAN, map_location="cpu")
        spec_pe_genes = [str(k).upper() for k in esm2.keys()]
        self.sym2row_local = {s: i for i, s in enumerate(spec_pe_genes)}   # local index within human block
        offsets = pickle.load(open(OFFSETS_PKL, "rb"))
        self.offset = int(offsets["human"])
        del esm2

        # ---- chromosome categorical codes (full multi-species table, matches get_spec_chrom_csv) ----
        chrom_df = pd.read_csv(CHROM_CSV)
        chrom_df["spec_chrom"] = pd.Categorical(chrom_df["species"] + "_" + chrom_df["chromosome"])
        human = chrom_df[chrom_df["species"] == "human"].copy()
        human["code"] = human["spec_chrom"].cat.codes
        human = human.drop_duplicates(subset="gene_symbol", keep="first").set_index("gene_symbol")
        self.chrom_code = human["code"].to_dict()          # gene symbol (as in csv) -> chrom category code
        self.gene_start = human["start"].to_dict()

        # ---- sampling universe: genes with BOTH an ESM2 row AND chrom/start info ----
        uni = [g for g in self.sym2row_local if g in self.chrom_code]
        self.uni_syms = uni
        self.uni_pe_row = np.array([self.sym2row_local[g] + self.offset for g in uni], dtype=np.int64)
        self.uni_chrom = np.array([int(self.chrom_code[g]) for g in uni], dtype=np.int64)
        self.uni_start = np.array([self.gene_start[g] for g in uni], dtype=np.int64)
        self.uni_idx = {g: i for i, g in enumerate(uni)}    # symbol -> index into uni_* arrays
        print(f"[UCEReadout] sampling universe: {len(uni)} genes with ESM2+chrom/start info", flush=True)

        # ---- manifest gene -> universe index (None if not in universe) ----
        self.manifest_uni = np.array([self.uni_idx.get(g.upper(), -1) for g in self.gsyms])
        cov = int((self.manifest_uni >= 0).sum())
        print(f"[UCEReadout] manifest genes present in sampling universe: {cov}/{self.Ng}", flush=True)

        # ---- model + checkpoint's own pe_embedding.weight (authoritative, see module docstring) ----
        self.model = TransformerModel(token_dim=5120, d_model=1280, nhead=20, d_hid=5120,
                                      nlayers=4, output_dim=1280, dropout=0.0)
        self.model.pe_embedding = nn.Embedding.from_pretrained(torch.zeros(N_ROWS, 5120))
        sd = torch.load(CKPT, map_location="cpu")
        self.model.load_state_dict(sd, strict=True)
        self.model = self.model.to(DEV).eval()
        self.pe_weight = self.model.pe_embedding.weight.data  # [145469, 5120] on DEV

    def _sample_sentence(self, gene_idx_present, weights):
        """gene_idx_present: indices into uni_* for this cell's nonzero, universe-covered genes.
        weights: matching normalized sampling probability (log1p-count based). Returns ordered_idx
        (pad_length,) of universe-array indices or -1 for special/pad tokens, plus a parallel
        token_id array (row into pe_weight) and the valid sequence length."""
        choice = self.rng.choice(gene_idx_present, size=self.sample_size, p=weights, replace=True)
        chrom = self.uni_chrom[choice]
        order = np.argsort(chrom, kind="stable")
        choice = choice[order]; chrom = chrom[order]; starts = self.uni_start[choice]

        ordered_uni = np.full(self.pad_length, -1, dtype=np.int64)
        ordered_tok = np.full(self.pad_length, PAD_TOK, dtype=np.int64)
        ordered_uni[0] = -1; ordered_tok[0] = CLS_TOK
        i = 1
        uq = np.unique(chrom); self.rng.shuffle(uq)
        for c in uq:
            if i >= self.pad_length - 1: break
            ordered_tok[i] = int(c) + CHROM_TOKEN_OFFSET; i += 1
            loc = np.where(chrom == c)[0]
            loc = loc[np.argsort(starts[loc], kind="stable")]
            take = choice[loc]
            n = min(len(take), self.pad_length - i - 1)
            ordered_uni[i:i + n] = take[:n]
            ordered_tok[i:i + n] = self.uni_pe_row[take[:n]]
            i += n
            if i < self.pad_length: ordered_tok[i] = CHROM_R_TOK; i += 1
        return ordered_uni, ordered_tok, i

    def gene_embed(self, Xraw_genes, cell_ids, uni_col_idx):
        """Xraw_genes: csr (cells x len(uni_syms)) RAW counts, aligned to self.uni_syms (via uni_col_idx
        upstream). Returns [Ng, 1280] averaged embedding for manifest genes covered."""
        gsum = np.zeros((self.Ng, 1280), np.float64); gcnt = np.zeros(self.Ng, np.float64)
        uni2manifest = {}
        for mi, ui in enumerate(self.manifest_uni):
            if ui >= 0: uni2manifest.setdefault(ui, []).append(mi)

        for s in range(0, len(cell_ids), self.batch):
            batch_ids = cell_ids[s:s + self.batch]
            uni_arrs, tok_arrs, lens = [], [], []
            for r in batch_ids:
                x = Xraw_genes[r].toarray().ravel()
                nz = np.where(x > 0)[0]
                if len(nz) == 0: uni_arrs.append(None); continue
                w = np.log1p(x[nz]); w = w / w.sum()
                ou, ot, ln = self._sample_sentence(nz, w)
                uni_arrs.append(ou); tok_arrs.append(ot); lens.append(ln)
            valid = [i for i, u in enumerate(uni_arrs) if u is not None]
            if not valid: continue
            L = max(lens)
            tok_mat = np.stack([tok_arrs[i][:L] for i in valid])            # [b, L]
            uni_mat = np.stack([uni_arrs[i][:L] for i in valid])            # [b, L]
            mask = np.zeros((len(valid), L), dtype=bool)                    # True = PAD (ignore)
            for bi, i in enumerate(valid): mask[bi, lens[i]:] = True

            tok_t = torch.tensor(tok_mat, device=DEV).long()
            src = self.pe_weight[tok_t]                                    # [b, L, 5120]
            src = nn.functional.normalize(src, dim=2)
            src = src.permute(1, 0, 2)                                     # [L, b, 5120] (seq_len, batch, dim)
            pad_mask = torch.tensor(mask, device=DEV).bool()               # [b, L] True=ignore
            with torch.no_grad():
                # bypass TransformerModel.forward()'s internal `(1-mask)` (relies on version-fragile
                # float-mask coercion) — call the submodules directly with an unambiguous bool mask.
                m = self.model
                h = m.encoder(src) * (1280 ** 0.5)
                h = m.pos_encoder(h)
                enc_out = m.transformer_encoder(h, src_key_padding_mask=pad_mask)
                gene_output = m.decoder(enc_out)                            # [L, b, 1280]
            go = gene_output.float().cpu().numpy()
            for bi, i in enumerate(valid):
                for p in range(1, lens[i]):                                 # skip CLS at 0
                    ui = uni_mat[bi, p]
                    if ui < 0: continue                                     # chrom bracket token
                    for mi in uni2manifest.get(int(ui), ()):
                        gsum[mi] += go[p, bi]; gcnt[mi] += 1
            del src, gene_output, tok_t, pad_mask
            if DEV == "cuda": torch.cuda.empty_cache()
        gcnt[gcnt == 0] = 1
        return (gsum / gcnt[:, None]).astype(np.float32)
