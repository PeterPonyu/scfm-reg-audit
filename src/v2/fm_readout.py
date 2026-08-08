#!/usr/bin/env python
"""
scfm-reg-audit v2 — corrected FM gene-graph readout (shared module).

Fixes vs pilot (peer review 2026-07-13):
  P1  Geneformer ranks NON-LOG CP10k / gene_median (pilot ranked log1p values -> token-order bug).
  P2  Never/low-expressed genes are dropped upstream (frozen manifest); std==0 guarded here.
  P4  random-init Geneformer is SEEDED (reproducible floor).
  P5  scGPT bins log1p CP10k (its intended input); order-invariant model.
  P6  strict=False loads assert a matched-parameter count so silent weight drops surface.

Normalization contract (set by the caller, see load_norm()):
  X_cp10k = counts / rowsum * 1e4            (NO log)  -> Geneformer ranking
  X_log   = log1p(X_cp10k)                              -> scGPT value bins + co-expression
Geneformer gene_median_dictionary is on the CP10k scale, so CP10k/median is the correct rank value.
"""
import os, json, pickle, numpy as np, scipy.sparse as sp, torch, torch.nn as nn
from scipy.stats import spearmanr
DATA_ROOT = os.environ.get("SCREG_DATA_ROOT", os.path.join(os.path.dirname(__file__), "..", "..", "data"))

GF = f"{DATA_ROOT}/models/Geneformer"; GFM = f"{GF}/Geneformer-V2-104M"; DCT = f"{GF}/geneformer"
CK = f"{DATA_ROOT}/models/scGPT-human"
DEV = "cuda" if torch.cuda.is_available() else "cpu"


def load_norm(path, ctcol="cell_type"):
    """Return (A, X_cp10k csr, X_log csr, labels). Detects raw-count vs already-normalized input."""
    import anndata as ad
    A = ad.read_h5ad(path)
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    mx = float(X[:min(1000, X.shape[0])].max())
    if mx > 30 and mx == int(mx):                       # raw counts
        rs = np.asarray(X.sum(1)).ravel(); rs[rs == 0] = 1
        Xc = X.multiply(1e4 / rs[:, None]).tocsr()
    else:                                               # already normalized-ish: reconstruct CP10k scale
        # if it looks like log1p already, invert; else assume CP10k
        Xc = X.copy().tocsr()
        if mx < 15:                                     # plausibly log1p
            Xc.data = np.expm1(Xc.data)
            rs = np.asarray(Xc.sum(1)).ravel(); rs[rs == 0] = 1
            Xc = Xc.multiply(1e4 / rs[:, None]).tocsr()
    Xl = Xc.copy(); Xl.data = np.log1p(Xl.data)
    labels = A.obs[ctcol].astype(str).values if ctcol in A.obs else None
    return A, Xc, Xl, labels


def gene_coexp(Xd):
    """|Pearson| gene-gene graph from a (cells x genes) dense log-normalized block."""
    M = Xd - Xd.mean(0, keepdims=True); sd = M.std(0, keepdims=True); sd[sd == 0] = 1
    C = ((M / sd).T @ (M / sd)) / M.shape[0]
    return np.abs(C).astype(np.float32)


class FMReadout:
    """Frozen-manifest gene embeddings. gsyms is the pre-registered ordered gene list."""

    def __init__(self, gsyms, batch=6):
        self.gsyms = list(gsyms); self.Ng = len(gsyms); self.batch = batch
        self.gtok = pickle.load(open(f"{DCT}/token_dictionary_gc104M.pkl", "rb"))
        self.gmed = pickle.load(open(f"{DCT}/gene_median_dictionary_gc104M.pkl", "rb"))
        self.gn2i = pickle.load(open(f"{DCT}/gene_name_id_dict_gc104M.pkl", "rb"))
        self.svoc = json.load(open(f"{CK}/vocab.json"))
        for s in self.gsyms:                            # manifest guarantees membership; assert anyway
            assert s in self.gn2i and self.gn2i[s] in self.gtok and s in self.svoc, f"gene {s} not tokenizable"
        self.gf_tokens = np.array([self.gtok[self.gn2i[s]] for s in self.gsyms])
        self.gf_median = np.array([self.gmed[self.gn2i[s]] for s in self.gsyms], np.float32)
        self.sg_ids = np.array([self.svoc[s] for s in self.gsyms])

    # ---------- Geneformer (P1: rank NON-LOG CP10k / median) ----------
    def _gf_model(self, random_init=False, seed=0):
        from transformers.models.bert.modeling_bert import BertModel, BertConfig
        if random_init:
            torch.manual_seed(seed); np.random.seed(seed)                       # P4
            cfg = BertConfig.from_pretrained(GFM); cfg.output_hidden_states = True
            m = BertModel(cfg, add_pooling_layer=False)
        else:
            m = BertModel.from_pretrained(GFM, output_hidden_states=True, add_pooling_layer=False)
        return m.to(DEV).half().eval() if DEV == "cuda" else m.eval()

    def geneformer(self, Xcp10k_genes, cell_ids, random_init=False, seed=0):
        """Xcp10k_genes: csr (cells x Ng) on CP10k (NO log). Returns E [Ng, HID]."""
        m = self._gf_model(random_init, seed); HID = int(m.config.hidden_size)
        CLS, EOS, PAD = self.gtok["<cls>"], self.gtok["<eos>"], self.gtok["<pad>"]
        gsum = np.zeros((self.Ng, HID), np.float32); gcnt = np.zeros(self.Ng, np.float32)
        for s in range(0, len(cell_ids), self.batch):
            rows, orders, lens = [], [], []
            for r in cell_ids[s:s + self.batch]:
                x = Xcp10k_genes[r].toarray().ravel()
                val = x / self.gf_median                 # <-- CP10k / median, NOT log1p
                nz = np.where(val > 0)[0]
                order = nz[np.argsort(-val[nz])]         # Geneformer rank encoding
                rows.append([CLS] + self.gf_tokens[order].tolist() + [EOS]); orders.append(order); lens.append(len(order) + 2)
            L = max(lens); ic = np.full((len(rows), L), PAD, np.int64); am = np.zeros((len(rows), L), np.int64)
            for i, r in enumerate(rows): ic[i, :len(r)] = r; am[i, :len(r)] = 1
            with torch.no_grad():
                h = m(input_ids=torch.tensor(ic, device=DEV), attention_mask=torch.tensor(am, device=DEV)).hidden_states[-1].float().cpu().numpy()
            for i, o in enumerate(orders): gsum[o] += h[i, 1:1 + len(o)]; gcnt[o] += 1
        del m; torch.cuda.empty_cache() if DEV == "cuda" else None
        gcnt[gcnt == 0] = 1
        return gsum / gcnt[:, None]

    # ---------- scGPT (P5: bin log1p CP10k) ----------
    def scgpt(self, Xlog_genes, cell_ids):
        args = json.load(open(f"{CK}/args.json")); PAD = self.svoc[args["pad_token"]]
        CLS = self.svoc.get("<cls>", self.svoc.get("<CLS>")); D = args["embsize"]; NBIN = args.get("n_bins", 51)

        class CVE(nn.Module):
            def __init__(s, d): super().__init__(); s.linear1 = nn.Linear(1, d); s.activation = nn.ReLU(); s.linear2 = nn.Linear(d, d); s.norm = nn.LayerNorm(d)
            def forward(s, x): x = torch.clamp(x.unsqueeze(-1), max=512); return s.norm(s.linear2(s.activation(s.linear1(x))))

        class Model(nn.Module):
            def __init__(s):
                super().__init__(); s.gene_emb = nn.Embedding(len(self.svoc), D, padding_idx=PAD); s.enc_norm = nn.LayerNorm(D); s.value_encoder = CVE(D)
                s.transformer_encoder = nn.TransformerEncoder(nn.TransformerEncoderLayer(D, args["nheads"], args["d_hid"], dropout=0.0, batch_first=True), args["nlayers"])
            def encode(s, src, val, mask): return s.transformer_encoder(s.enc_norm(s.gene_emb(src)) + s.value_encoder(val), src_key_padding_mask=mask)

        m = Model().eval(); sd = torch.load(f"{CK}/best_model.pt", map_location="cpu", weights_only=False)
        if isinstance(sd, dict) and "model_state_dict" in sd: sd = sd["model_state_dict"]
        rm = {}
        for k, v in sd.items():
            nk = k.replace("encoder.embedding.", "gene_emb.").replace("encoder.enc_norm.", "enc_norm.").replace(".self_attn.Wqkv.weight", ".self_attn.in_proj_weight").replace(".self_attn.Wqkv.bias", ".self_attn.in_proj_bias")
            rm[nk] = v
        info = m.load_state_dict(rm, strict=False)       # P6: assert nothing important dropped
        loaded = len(m.state_dict()) - len(info.missing_keys)
        assert len(info.missing_keys) <= 2, f"scGPT missing {len(info.missing_keys)} params: {info.missing_keys[:8]}"
        m.to(DEV)

        def binv(row):
            if len(row) == 0: return row.astype(int)
            b = np.quantile(row, np.linspace(0, 1, NBIN - 1)); return np.clip(np.digitize(row, b), 1, NBIN - 1)
        gsum = np.zeros((self.Ng, D), np.float32); gcnt = np.zeros(self.Ng, np.float32)
        for s in range(0, len(cell_ids), self.batch):
            rows, vals, orders = [], [], []
            for r in cell_ids[s:s + self.batch]:
                x = Xlog_genes[r].toarray().ravel(); nz = np.where(x > 0)[0]
                rows.append([CLS] + self.sg_ids[nz].tolist()); vals.append([0.0] + binv(x[nz]).astype(float).tolist()); orders.append(nz)
            L = max(len(r) for r in rows); ic = np.full((len(rows), L), PAD, np.int64); vv = np.zeros((len(rows), L), np.float32); am = np.ones((len(rows), L), bool)
            for i, (rr, vl) in enumerate(zip(rows, vals)): ic[i, :len(rr)] = rr; vv[i, :len(vl)] = vl; am[i, :len(rr)] = False
            with torch.no_grad():
                h = m.encode(torch.tensor(ic, device=DEV), torch.tensor(vv, device=DEV), torch.tensor(am, device=DEV)).float().cpu().numpy()
            for i, o in enumerate(orders): gsum[o] += h[i, 1:1 + len(o)]; gcnt[o] += 1
        del m; torch.cuda.empty_cache() if DEV == "cuda" else None
        gcnt[gcnt == 0] = 1
        return gsum / gcnt[:, None]

    # ---------- Geneformer attention (trustworthy readout; last layer, mean heads, symmetrized) ----------
    def geneformer_attention(self, Xcp10k_genes, cell_ids):
        from transformers.models.bert.modeling_bert import BertModel
        m = BertModel.from_pretrained(GFM, output_attentions=True, add_pooling_layer=False)
        m = m.to(DEV).half().eval() if DEV == "cuda" else m.eval()
        CLS, EOS, PAD = self.gtok["<cls>"], self.gtok["<eos>"], self.gtok["<pad>"]
        Asum = np.zeros((self.Ng, self.Ng), np.float64); Acnt = np.zeros((self.Ng, self.Ng), np.float64)
        for s in range(0, len(cell_ids), self.batch):
            rows, orders, lens = [], [], []
            for r in cell_ids[s:s + self.batch]:
                x = Xcp10k_genes[r].toarray().ravel(); val = x / self.gf_median
                nz = np.where(val > 0)[0]; order = nz[np.argsort(-val[nz])]
                rows.append([CLS] + self.gf_tokens[order].tolist() + [EOS]); orders.append(order); lens.append(len(order) + 2)
            L = max(lens); ic = np.full((len(rows), L), PAD, np.int64); am = np.zeros((len(rows), L), np.int64)
            for i, r in enumerate(rows): ic[i, :len(r)] = r; am[i, :len(r)] = 1
            with torch.no_grad():
                att = m(input_ids=torch.tensor(ic, device=DEV), attention_mask=torch.tensor(am, device=DEV)).attentions[-1]
                att = att.float().mean(1).cpu().numpy()
            for i, order in enumerate(orders):
                k = len(order); sub = att[i, 1:1 + k, 1:1 + k]
                Asum[np.ix_(order, order)] += sub; Acnt[np.ix_(order, order)] += 1
        del m; torch.cuda.empty_cache() if DEV == "cuda" else None
        Acnt[Acnt == 0] = 1; Aatt = (Asum / Acnt).astype(np.float32)
        return np.abs((Aatt + Aatt.T) / 2)

    @staticmethod
    def cos_graph(E, center=True):
        if center: E = E - E.mean(0, keepdims=True)
        En = E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-8)
        return np.abs(En @ En.T).astype(np.float32)
