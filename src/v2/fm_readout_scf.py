#!/usr/bin/env python
"""
scfm-reg-audit v2 — scFoundation gene-embedding readout (3rd FM, per DESIGN's pre-registered order).

Architecture/loading code vendored from research/sc-fm-benchmark/scripts/scf_repo/ (itself a vendored
copy of biomap-research/scFoundation's model/ folder) and scfoundation_embed.py (the validated,
dependency-free reimplementation used in that project — the official decoder path needs an
uninstallable `pretrainmodels` package + missing reversible.py, so only the ENCODER-ONLY forward is
runnable here). This mirrors exactly what scfoundation_embed.py already validated for CELL embeddings;
the new part is capturing PER-GENE encoder hidden states (before pooling) and scattering them onto our
own frozen 1200-gene manifest, matching the FMReadout.geneformer()/scgpt() pattern in fm_readout.py so
this plugs directly into the existing crossmodal/confound-regression test suite.

Caveat (report alongside any result): this is scFoundation's ENCODER-ONLY embedding (no performer
decoder refinement) — analogous to using Geneformer's raw encoder hidden states, not an official
"gene embedding mode". Per-cell gene set is capped at top-CAP nonzero genes by expression (default
2046, matching the Geneformer token budget) — genes outside a cell's top-CAP are simply absent from
that cell's forward pass, structurally analogous to Geneformer's own rank-truncation behavior.
"""
import os, sys, json, numpy as np, pandas as pd, torch
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scf_vendor"))
import mae_autobin, transformer
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
CKPT = f"{DATA_ROOT}/models/scFoundation-cell/model.pt"
GIDX = f"{ROOT}/data/scfoundation/scf_gene_index.tsv"
DEV = "cuda" if torch.cuda.is_available() else "cpu"


def gatherData(data, labels, pad_token_id):
    value_nums = labels.sum(1); max_num = int(max(value_nums))
    fake_data = torch.full((data.shape[0], max_num), pad_token_id, device=data.device)
    data = torch.hstack([data, fake_data])
    fake_label = torch.full((labels.shape[0], max_num), 1, device=labels.device)
    none_labels = ~labels; labels = labels.float()
    labels[none_labels] = torch.tensor(-float("inf"), device=labels.device)
    tmp = torch.tensor([(i + 1) * 20000 for i in range(labels.shape[1], 0, -1)], device=labels.device)
    labels = labels + tmp; labels = torch.hstack([labels, fake_label])
    idx = labels.topk(max_num).indices
    new_data = torch.gather(data, 1, idx); padding = (new_data == pad_token_id)
    return new_data, padding


class SCFReadout:
    """Frozen-manifest gene embeddings via scFoundation's encoder (cell-only checkpoint)."""

    def __init__(self, gsyms, batch=2, cap=2046):
        self.gsyms = list(gsyms); self.Ng = len(gsyms); self.batch = batch; self.CAP = cap
        glist = list(pd.read_csv(GIDX, header=0, delimiter="\t")["gene_name"])
        self.sym2i = {g: i for i, g in enumerate(glist)}; self.NG = len(glist); assert self.NG == 19264
        self.PAD_ID = 103
        gi = [i for i, s in enumerate(self.gsyms) if s in self.sym2i]
        self.manifest_ok = np.array(gi)                                  # index into gsyms
        self.manifest_scf_idx = np.array([self.sym2i[self.gsyms[i]] for i in gi])  # index into 19264 vocab
        assert len(gi) > len(self.gsyms) * 0.9, f"low manifest overlap: {len(gi)}/{len(self.gsyms)}"

        sd = torch.load(CKPT, map_location="cpu", weights_only=True)
        if isinstance(sd, dict) and "cell" in sd: sd = sd["cell"]
        SEQ = sd["pos_emb.weight"].shape[0] - 1
        m = mae_autobin.MaeAutobin(num_tokens=self.NG, max_seq_len=SEQ, embed_dim=768, decoder_embed_dim=512,
                                   bin_alpha=1.0, bin_num=100, pad_token_id=self.PAD_ID, mask_token_id=102)
        m.encoder = transformer.pytorchTransformerModule(max_seq_len=SEQ, dim=768, depth=12, heads=12)
        m.decoder = torch.nn.Identity()
        miss, unexp = m.load_state_dict(sd, strict=False)
        miss = [k for k in miss if not k.startswith("decoder.")]
        assert not miss, f"scFoundation missing (non-decoder) keys: {miss}"
        self.HALF = DEV == "cuda"
        self.model = m.to(DEV).eval()
        if self.HALF: self.model = self.model.half()
        self.SEQ = SEQ
        self.gene_ids_row = torch.arange(SEQ, device=DEV)                # [SEQ]

    def _vlab(self, mat):
        gene = mat[:, :self.NG]; vlab = torch.zeros_like(mat, dtype=torch.bool)
        if self.NG > self.CAP:
            topv, topi = gene.topk(self.CAP, dim=1); vlab[:, :self.NG].scatter_(1, topi, topv > 0)
        else:
            vlab[:, :self.NG] = gene > 0
        vlab[:, self.NG:] = mat[:, self.NG:] > 0
        return vlab

    def embed(self, Xc19264_cells, gsum, gcnt):
        """Xc19264_cells: (n_cells, NG) dense np.float64 CP10k (NOT log). Accumulates into gsum/gcnt
        [Ng_manifest, 768] in-place for manifest genes present in each cell's gathered token set."""
        n = Xc19264_cells.shape[0]
        for s in range(0, n, self.batch):
            rows = []
            for r in range(s, min(s + self.batch, n)):
                v = Xc19264_cells[r]; tot = v.sum()
                if tot <= 0: continue
                logn = np.log1p(v / tot * 1e4); lt = np.log10(tot)
                rows.append(np.concatenate([logn, [lt, lt]]))
            if not rows: continue
            mat = torch.tensor(np.stack(rows), dtype=torch.float32, device=DEV)
            vlab = self._vlab(mat)
            x, xpad = gatherData(mat, vlab, self.PAD_ID)
            gene_ids = self.gene_ids_row.unsqueeze(0).repeat(mat.shape[0], 1)
            pos, _ = gatherData(gene_ids, vlab, self.PAD_ID)
            xt = x.unsqueeze(2).half() if self.HALF else x.unsqueeze(2).float()
            with torch.no_grad():
                h = self.model.token_emb(xt, output_weight=0) + self.model.pos_emb(pos)
                g = self.model.encoder(h, padding_mask=xpad)             # [b, L, 768]
            g = g.float().cpu().numpy(); posn = pos.cpu().numpy()
            for bi in range(g.shape[0]):
                # map this cell's gathered token positions -> manifest-gene rows (skip resolution/pad tokens)
                tok2man = {int(self.manifest_scf_idx[k]): k for k in range(len(self.manifest_scf_idx))}
                for li in range(g.shape[1]):
                    gid = int(posn[bi, li])
                    if gid in tok2man:
                        k = tok2man[gid]; gsum[self.manifest_ok[k]] += g[bi, li]; gcnt[self.manifest_ok[k]] += 1
            del h, g, xt, x, pos, xpad, vlab, mat
            if DEV == "cuda": torch.cuda.empty_cache()

    def gene_embed(self, X_manifest_cp10k, cell_ids):
        """X_manifest_cp10k: csr (cells x Ng_manifest) on CP10k (NOT log), aligned to self.gsyms.
        Returns [Ng, 768] averaged embedding (zero rows for genes never gathered)."""
        Xd = np.zeros((len(cell_ids), self.NG), np.float64)
        cols = self.manifest_scf_idx; rows_in_manifest = self.manifest_ok
        sub = X_manifest_cp10k[cell_ids][:, rows_in_manifest].toarray()
        Xd[:, cols] = sub
        gsum = np.zeros((self.Ng, 768), np.float64); gcnt = np.zeros(self.Ng, np.float64)
        self.embed(Xd, gsum, gcnt)
        gcnt[gcnt == 0] = 1
        return (gsum / gcnt[:, None]).astype(np.float32)
