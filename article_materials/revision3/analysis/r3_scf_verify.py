#!/usr/bin/env python3
"""R3-7: verify the evaluated scFoundation weights against the developers' published checkpoint.

Step 1 (no download): the local file's SHA-256 equals the Hugging Face LFS SHA-256 of
`perturblab/scfoundation-cell` model.pt at commit 7b06458e, the mirror it was obtained from.
Step 2 (needs the official `models.ckpt`, 1,432,587,886 bytes, from `genbio-ai/scFoundation`,
LFS SHA-256 9f40bf32...0ea9): hash the official file, then compare every tensor of each
sub-model (`cell`, `gene`, `rde`) with the local state dict by exact equality.

Usage: python r3_scf_verify.py [--official PATH_TO_models.ckpt]
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import torch

W = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from r3paths import need  # noqa: E402

LOCAL = need('scfoundation_weights')
MIRROR = {'repo': 'perturblab/scfoundation-cell', 'commit': '7b06458e108e0690fdd7f2ece61a4658069b801a',
          'file': 'model.pt', 'bytes': 477485967,
          'lfs_sha256': '2446c2fb99a0e183cdc9872b7e9e55ac74f7ec82db7ba0799c83be75cde66dc6'}
OFFICIAL = {'repo': 'genbio-ai/scFoundation', 'commit': 'cb434153a1acfacd215eefc956ea445f7cc39cc3',
            'file': 'models.ckpt', 'bytes': 1432587886,
            'lfs_sha256': '9f40bf324d3d0084c4b288d06f5af4fddd12206e2a3f022551d12e89e33a0ea9',
            'distributor': 'developers of scFoundation (BioMap); code at github.com/biomap-research/scFoundation'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def tensors(obj, prefix=''):
    """Flatten nested dicts to {dotted_key: tensor} for comparison."""
    out = {}
    if isinstance(obj, torch.Tensor):
        out[prefix.rstrip('.')] = obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            out.update(tensors(v, f'{prefix}{k}.'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--official', type=Path)
    a = ap.parse_args()
    report = {'local_file': str(LOCAL), 'local_bytes': LOCAL.stat().st_size, 'local_sha256': sha256(LOCAL), 'mirror': MIRROR}
    report['local_equals_mirror_file'] = report['local_sha256'] == MIRROR['lfs_sha256'] and report['local_bytes'] == MIRROR['bytes']
    local = torch.load(LOCAL, map_location='cpu', weights_only=True)
    report['local_tensors'] = len(local)
    if a.official:
        report['official'] = OFFICIAL
        report['official_sha256'] = sha256(a.official)
        report['official_bytes'] = a.official.stat().st_size
        report['official_file_matches_published_hash'] = report['official_sha256'] == OFFICIAL['lfs_sha256']
        ckpt = torch.load(a.official, map_location='cpu', weights_only=True)
        report['official_top_keys'] = sorted(ckpt) if isinstance(ckpt, dict) else str(type(ckpt))
        comparisons = {}
        for key in (ckpt if isinstance(ckpt, dict) else {}):
            sub = ckpt[key]
            sd = sub.get('model_state_dict', sub) if isinstance(sub, dict) else None
            if not isinstance(sd, dict):
                continue
            flat = tensors(sd)
            common = sorted(set(flat) & set(local))
            equal = [k for k in common if flat[k].shape == local[k].shape and flat[k].dtype == local[k].dtype and torch.equal(flat[k], local[k])]
            maxdiff = max((float((flat[k].double() - local[k].double()).abs().max()) for k in common
                           if flat[k].shape == local[k].shape), default=None)
            comparisons[key] = {'official_tensors': len(flat), 'shared_names': len(common), 'exactly_equal': len(equal),
                                'only_official': len(set(flat) - set(local)), 'only_local': len(set(local) - set(flat)),
                                'max_abs_difference': maxdiff,
                                'identical_state_dict': len(common) == len(flat) == len(local) == len(equal)}
        report['comparisons'] = comparisons
        report['local_identical_to_official_submodel'] = [k for k, v in comparisons.items() if v['identical_state_dict']]
    out = W / 'evidence/scfoundation_verification.json'
    out.write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('mirror', 'official')}, indent=1, default=str))


if __name__ == '__main__':
    main()
