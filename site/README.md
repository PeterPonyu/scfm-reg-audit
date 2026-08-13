# GitHub Pages source (`site/`)

H1+G1 static materials site. Public routes: `/` `/peerj/` `/frontiers/` `/figures/` `/reproducibility/`. Visual tokens follow `DESIGN.md`; routes follow `.omx/specs/github-pages-ia.md`.

Intended production URL (after a human enables Pages): `https://PeterPonyu.github.io/scfm-reg-audit/`. That URL is **not** live until Settings → Pages → Source = GitHub Actions succeeds. Do not set `homepageUrl` before the first 200.

Header may also include outbound GitHub blob links labeled Protocol and Cite (approved hub raster); those are not Pages routes.

## Local assemble (PR-review path)

```bash
bash site/assemble.sh
python3 -m http.server 4173 --directory site/_site
```

Output directory is `site/_site/` (gitignored). It must contain `.nojekyll`, the five HTML routes, `styles.css`, and **exactly two** PDFs:

- `products/peerj-manuscript.pdf` ← `paper/manuscript.pdf`
- `products/frontiers-manuscript.pdf` ← `paper/submission_frontiers_genetics/manuscript.pdf`

Missing SoT PDF → non-zero exit. Do not commit PDFs under `site/`.

Hub screenshot (1280×800):

```bash
npx --yes playwright screenshot --viewport-size=1280,800 \
  http://127.0.0.1:4173/ \
  .omx/artifacts/visual-ralph/github-pages-hub/screenshot-1280x800.png
```

## Human enablement (not done by git)

1. Merge `site/` + `.github/workflows/pages.yml` to `main` (user commit/push).
2. Settings → Pages → Source = **GitHub Actions**.
3. Approve the `github-pages` environment if GitHub requires reviewers. Restrict it to default branch `main`.
4. Confirm `GET https://PeterPonyu.github.io/scfm-reg-audit/` → 200.
5. Then optionally set repo `homepageUrl` to that URL.

## H2 fallback (not v1)

If Actions cannot be enabled: run `bash site/assemble.sh` and publish the contents of `site/_site/` as the root of a `gh-pages` branch. Do not mix `deploy-pages` with a force-push to `gh-pages`.
