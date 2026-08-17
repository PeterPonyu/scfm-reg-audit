# scReg-Eval repository topology

Last updated: 2026-08-09.

## Why there are two local trees

This research lives in **two intentional places** with different jobs. They are
not accidental duplicates.

| Location | Role | Size (approx.) | Git |
| --- | --- | --- | --- |
| `~/Desktop/scfm-reg-audit` | **Publish capsule** (PeerJ CS source of truth for manuscript + validated code/artifacts) | ~34 MB | Independent repo `main` |
| `~/Desktop/singlecell-genomics-research/projects/scfm-reg-audit` | **Full research workspace** (compute, heavy results, caches, drafts) | ~5 GB | Folder inside monorepo `singlecell-genomics-research` (no nested `.git`) |

### Independent GitHub repo (already exists)

- URL: https://github.com/PeterPonyu/scfm-reg-audit
- Default branch: `main`
- Description: scReg-Eval fixed-panel audit capsule (code MIT; manuscript/results CC BY 4.0)
- Local main worktree: `~/Desktop/scfm-reg-audit`

The monorepo also has a **named remote** `scfm-reg-audit` pointing at the same
GitHub URL. That remote is only a convenience for cherry-picks/sync; the monorepo
folder itself is **not** a checkout of that repo.

### Why Desktop “suddenly” looked independent

The Desktop tree is the **sanitized audit capsule** extracted for submission and
public validation (`README.md` calls it “audit capsule v0.4.0”). The monorepo
path remains the heavy research home (GPU plans, raw-ish derivatives, enhancement
runs). Paper edits in this session were made on the **capsule** path, so
`paper/manuscript.*` can diverge from the monorepo copy until explicitly synced.

## Canonical rules

1. **Manuscript / figure source of truth for PeerJ:** `~/Desktop/scfm-reg-audit/paper/`
2. **Heavy compute / large results:** monorepo `projects/scfm-reg-audit/`
3. **Do not** put multi‑GB NPZ/H5AD into the capsule repo.
4. After paper polish on the capsule, either:
   - copy `paper/` into the monorepo path for research continuity, or
   - open a PR/commit on `PeterPonyu/scfm-reg-audit` and pull there when needed.

## Git worktrees (optional extra checkouts)

From the capsule main worktree you can add sibling worktrees without cloning
again:

```bash
# example: paper-only branch checkout elsewhere
git -C ~/Desktop/scfm-reg-audit worktree add \
  -b paper/peerj-visual-qa \
  ~/Desktop/labs/active/scReg-Eval-paper \
  main

git -C ~/Desktop/scfm-reg-audit worktree list
# remove when finished:
# git -C ~/Desktop/scfm-reg-audit worktree remove ~/Desktop/labs/active/scReg-Eval-paper
```

## Citation style (PeerJ CS)

PeerJ Computer Science uses **Name–Year** with an **alphabetized** reference
list (`apalike` + `\setcitestyle{authoryear,...}`). Not numbered
citation-sequence.

## Related monorepo docs

See monorepo `projects/scfm-reg-audit/README.md` for research-scope boundaries
relative to other single-cell projects under `singlecell-genomics-research`.
