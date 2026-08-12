#!/usr/bin/env Rscript
# Extension SI figures (post-PeerJ construct wave).
# Statistical visualizations ONLY via R → tikzDevice → paper/figs_extension/.
# Does NOT write paper/figs/ (PeerJ fig1–12 freeze) or run the PeerJ package builder.
#
# Usage (from repo root or paper/):
#   Rscript paper/make_figs_extension.R
#   SCFM_BASE=.. Rscript make_figs_extension.R   # if cwd is paper/
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(jsonlite)
  library(tikzDevice)
})

# Resolve repo root: allow running from paper/ or repo root
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
if (basename(script_dir) == "paper") {
  default_base <- dirname(script_dir)
} else if (file.exists(file.path(script_dir, "paper", "make_figs.R"))) {
  default_base <- script_dir
} else {
  default_base <- ".."
}
base <- Sys.getenv("SCFM_BASE", default_base)
base <- normalizePath(base, mustWork = TRUE)

source(file.path(base, "src", "v2", "figure_helpers.R"))
setup_tikz_options()

figs_ext <- file.path(base, "paper", "figs_extension")
dir.create(figs_ext, showWarnings = FALSE, recursive = TRUE)

BLUE <- "#2a78d6"; AQUA <- "#1b8f75"; YELLOW <- "#d99a00"
VIOLET <- "#5946b2"; RED <- "#d84a4a"; MUTED <- "grey45"; LIGHT <- "#d9e8f8"

theme_set(theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  plot.title = element_text(size = 11, face = "plain", hjust = 0, margin = margin(b = 3)),
  strip.text = element_text(size = 10), legend.title = element_text(size = 9.5),
  legend.text = element_text(size = 9.5), axis.title = element_text(size = 10),
  axis.text = element_text(size = 9.5), plot.margin = margin(4, 6, 3, 4)))

J <- function(rel) fromJSON(file.path(base, rel), simplifyVector = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- emit into figs_extension only (never paper/figs) ----
emit_ext <- function(name, plot, width, height, tags = LETTERS[1:4]) {
  path <- file.path(figs_ext, paste0(name, ".tex"))
  tikz(path, width = width, height = height, standAlone = FALSE, sanitize = FALSE)
  print(plot)
  dev.off()
  validate_panel_tags(path, tags)
  cat(name, "ok ->", path, "\n")
  invisible(path)
}

# ============== Load + freeze guards ==============
si <- J("docs/reports/extension-claim-pack/construct_si_mantel.json")
if (isTRUE(si$peerj_support_rows_touched))
  stop("construct_si_mantel.json reports peerj_support_rows_touched=TRUE — refuse SI figures")
if (length(si$rows) != 9L)
  stop("expected 9 construct SI rows, got ", length(si$rows))
if (as.integer(si$n_pair_rows_ok) < 9L)
  stop("n_pair_rows_ok < 9")

hon <- J("docs/reports/extension-claim-pack/honesty_policy.json")
if (isTRUE(hon$peerj_support_rows_touched))
  stop("honesty_policy.json peerj_support_rows_touched=TRUE — refuse")
stopifnot(isTRUE(abs(hon$bmmc$gene_coverage - 0.8967) < 1e-4))
stopifnot(isTRUE(abs(hon$bmmc$tf_coverage - 0.7646) < 1e-4))

htan <- J("results/v2/extension/construct/HTAN_GBM_C3N01334/htan_prepare_status.json")
if (!identical(htan$status, "blocked"))
  stop("HTAN status expected blocked, got: ", htan$status)

panel <- J("paper/panel_data.json")
tt_cov <- panel$third_tissue$coverage

# tidy Mantel rows
proxy_lab <- c(GSE174367 = "Brain", PBMC10k = "PBMC", GSE206767 = "Fibroblast mix")
tissue_lab <- c(
  descartes_spleen = "DESCARTES spleen",
  bmmc = "BMMC",
  orphan_treg_gse211155 = "Treg pilot"
)
rows <- bind_rows(lapply(si$rows, function(r) {
  data.frame(
    tissue_id = r$tissue_id,
    tissue = unname(tissue_lab[[r$tissue_id]] %||% r$tissue_id),
    g_atac_tag = r$g_atac_tag,
    locked_proxy = r$locked_proxy,
    proxy = unname(proxy_lab[[r$locked_proxy]] %||% r$locked_proxy),
    rho = as.numeric(r$observed_spearman),
    frac_add = as.numeric(r$fraction_explained_by_additive_marginals),
    add_pred = as.numeric(r$additive_pred_spearman),
    residual = as.numeric(r$residual_spearman_after_own_additive_fits),
    n_tf = as.integer(r$n_tf_common),
    stringsAsFactors = FALSE
  )
}))
rows$tissue <- factor(rows$tissue, levels = c("DESCARTES spleen", "BMMC", "Treg pilot"))
rows$proxy <- factor(rows$proxy, levels = c("Brain", "PBMC", "Fibroblast mix"))

# relevant peaks from extension meta
meta_tags <- c(
  "DESCARTES spleen" = "DESCARTES_spleen",
  "BMMC" = "GSE194122",
  "Treg pilot" = "GSE211155_treg"
)
peak_df <- bind_rows(lapply(names(meta_tags), function(lab) {
  tag <- meta_tags[[lab]]
  m <- J(file.path("results/v2/extension/construct", tag,
                   paste0("G_ATAC_v2_", tag, "_meta.json")))
  data.frame(tissue = lab, source = "extension",
             relevant_peaks = as.numeric(m$relevant_peaks),
             stringsAsFactors = FALSE)
}))
peerj_peaks <- bind_rows(lapply(tt_cov, function(c) {
  data.frame(tissue = c$tissue, source = "PeerJ locked (ref)",
             relevant_peaks = as.numeric(c$relevant_peaks),
             stringsAsFactors = FALSE)
}))
# normalize PeerJ labels already Brain/PBMC/Fibroblast mix

# ============== E1: construct Mantel ==============
# Short labels for PeerJ-scale panels
rows$proxy_s <- factor(
  dplyr::recode(as.character(rows$proxy),
                Brain = "Brain", PBMC = "PBMC", `Fibroblast mix` = "Fibro"),
  levels = c("Brain", "PBMC", "Fibro"))
rows$tissue_s <- factor(
  dplyr::recode(as.character(rows$tissue),
                `DESCARTES spleen` = "Spleen", BMMC = "BMMC", `Treg pilot` = "Treg"),
  levels = c("Spleen", "BMMC", "Treg"))

# Continuous geom_tile fill → tikzDevice raster (leaks .png paths into view).
# Precompute hex colors; scale_fill_identity keeps vector tiles.
fill_hex <- function(x, lo, hi, c_lo, c_hi) {
  t <- pmin(1, pmax(0, (as.numeric(x) - lo) / (hi - lo)))
  a <- grDevices::col2rgb(c_lo)[, 1] / 255
  b <- grDevices::col2rgb(c_hi)[, 1] / 255
  grDevices::rgb((1 - t) * a[1] + t * b[1],
                 (1 - t) * a[2] + t * b[2],
                 (1 - t) * a[3] + t * b[3])
}
rows$fill_rho <- fill_hex(rows$rho, 0.35, 0.92, LIGHT, BLUE)
rows$fill_add <- fill_hex(rows$frac_add, 0.40, 0.90, "#f0f0f0", AQUA)

f1a <- ggplot(rows, aes(proxy_s, tissue_s, fill = fill_rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.0) +
  scale_fill_identity() +
  labs(x = "locked proxy", y = NULL, title = "Mantel rho vs locked G-ATAC") +
  theme(axis.text.x = element_text(size = 9))

f1b <- ggplot(rows, aes(proxy_s, tissue_s, fill = fill_add)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", frac_add)), size = 3.0) +
  scale_fill_identity() +
  labs(x = "locked proxy", y = NULL, title = "Additive fraction of rho") +
  theme(axis.text.x = element_text(size = 9))

rows$pair_s <- paste(as.character(rows$tissue_s), as.character(rows$proxy_s), sep = "-")
rows$pair_s <- factor(rows$pair_s, levels = rows$pair_s[order(rows$rho)])
f1c <- ggplot(rows, aes(y = pair_s)) +
  geom_segment(aes(x = 0.35, xend = rho, yend = pair_s), color = "grey70", linewidth = 0.5) +
  geom_point(aes(x = rho, color = tissue_s), size = 2.5) +
  scale_color_manual(values = c(Spleen = VIOLET, BMMC = YELLOW, Treg = BLUE)) +
  labs(x = "observed rho", y = NULL, color = NULL,
       title = "Transfer rho (BMMC-PBMC highest)") +
  theme(legend.position = "top", legend.text = element_text(size = 8.5))

# D: collapse air between title and panel top spine; D tag right+slightly down
f1d <- ggplot(peak_df, aes(reorder(tissue, relevant_peaks), relevant_peaks)) +
  geom_col(fill = BLUE, width = 0.65) +
  geom_text(aes(label = format(relevant_peaks, big.mark = ",")),
            hjust = -0.05, size = 2.9) +
  coord_flip(ylim = c(0, max(peak_df$relevant_peaks) * 1.22)) +
  labs(x = NULL, y = "linked peaks",
       title = "Motif-linked peaks by tissue") +
  theme(
    # title sits on the panel box (top spine), not floating in outer margin
    plot.title.position = "panel",
    plot.title = element_text(
      size = 9.5,
      margin = margin(t = 1, r = 0, b = 2, l = 0),
      hjust = 0, vjust = 1
    ),
    plot.margin = margin(t = 2, r = 8, b = 2, l = 2),
    # D: slightly right and down vs A/B/C top-left tags
    plot.tag.position = c(0.15, 0.94)
  )

f1a <- f1a + theme(plot.tag.position = c(0.02, 1.0))
f1b <- f1b + theme(plot.tag.position = c(0.02, 1.0))
f1c <- f1c + theme(plot.tag.position = c(0.02, 1.0))

emit_ext("fig_ext1_construct_mantel",
         (f1a | f1b) / (f1c | f1d) +
           plot_annotation(tag_levels = "A") +
           plot_layout(axis_titles = "collect"),
         6.8, 5.8, tags = LETTERS[1:4])

# ============== E2: literature prior + construct transfer stats ONLY ==============
# NO emitter status, freeze badges, file names, or inventory counts.
coll <- J("results/v2/extension/baselines/collectri_prior/summary.json")
cr <- coll$result
stopifnot(!is.null(cr$tf_coverage))
motif <- J("results/v2/extension/baselines/motif_only_rp/summary.json")$result

# A: CollecTRI fraction of panel TFs/genes hit + edge keep rate (scientific coverage)
cov_df <- data.frame(
  metric = factor(c("TF hit", "Gene hit", "Edge keep"),
                  levels = c("TF hit", "Gene hit", "Edge keep")),
  value = c(as.numeric(cr$tf_coverage),
            as.numeric(cr$gene_coverage),
            as.numeric(cr$n_edges_on_panel) / as.numeric(cr$n_edges_raw))
)
f2a <- ggplot(cov_df, aes(metric, value, fill = metric)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.35, size = 3.1) +
  scale_fill_manual(values = c(BLUE, AQUA, YELLOW)) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(x = NULL, y = "fraction", title = "CollecTRI panel hit rates")

# B: edge support scale — motif binary vs literature prior on panel (scientific)
edge_cmp <- data.frame(
  source = factor(c("Motif binary", "CollecTRI prior"),
                  levels = c("Motif binary", "CollecTRI prior")),
  n = c(as.numeric(motif$n_edges), as.numeric(cr$n_edges_on_panel))
)
f2b <- ggplot(edge_cmp, aes(source, n, fill = source)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = format(as.integer(n), big.mark = ",")), vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c(`Motif binary` = MUTED, `CollecTRI prior` = BLUE)) +
  labs(x = NULL, y = "TF-gene edges on panel",
       title = "Edge support: motif vs literature prior") +
  coord_cartesian(ylim = c(0, max(edge_cmp$n) * 1.14))

# C: mean observed Mantel rho per extension tissue (scientific transfer)
mean_rho <- rows %>%
  group_by(tissue_s) %>%
  summarise(rho = mean(rho), .groups = "drop")
f2c <- ggplot(mean_rho, aes(tissue_s, rho, fill = tissue_s)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", rho)), vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c(Spleen = VIOLET, BMMC = YELLOW, Treg = BLUE)) +
  coord_cartesian(ylim = c(0, max(mean_rho$rho) * 1.18)) +
  labs(x = NULL, y = "mean Mantel rho",
       title = "Mean transfer rho (3 locked proxies)")

# D: mean residual rho after additive (scientific decomp remainder)
mean_res <- rows %>%
  group_by(tissue_s) %>%
  summarise(resid = mean(residual), .groups = "drop")
f2d <- ggplot(mean_res, aes(tissue_s, resid, fill = tissue_s)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", resid)), vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c(Spleen = VIOLET, BMMC = YELLOW, Treg = BLUE)) +
  coord_cartesian(ylim = c(0, max(mean_res$resid) * 1.18)) +
  labs(x = NULL, y = "mean residual rho",
       title = "Residual after additive fit")

emit_ext("fig_ext2_baselines_collectri",
         (f2a | f2b) / (f2c | f2d) + plot_annotation(tag_levels = "A") +
           plot_layout(axis_titles = "collect"),
         6.8, 5.6, tags = LETTERS[1:4])

# ============== E3: tissue-level scientific stats ONLY ==============
# NO tar/inventory/filename/lake/freeze-badge packaging content.
treg_rows <- rows[rows$tissue_id == "orphan_treg_gse211155", , drop = FALSE]
if (nrow(treg_rows) == 0)
  treg_rows <- rows[as.character(rows$tissue_s) == "Treg", , drop = FALSE]
treg_rows$proxy_s <- factor(as.character(treg_rows$proxy_s), levels = c("Brain", "PBMC", "Fibro"))

# A: residual Spearman for all tissues x proxies (scientific)
f3a <- ggplot(rows, aes(proxy_s, tissue_s, fill = residual)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", residual)), size = 3.0) +
  scale_fill_gradient(low = LIGHT, high = VIOLET, name = NULL) +
  # continuous residual heatmap risks raster; use identity hex
  labs(x = "locked proxy", y = NULL, title = "Residual rho after additive") +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

# rebuild residual with vector fills (no continuous gradient raster)
rows$fill_res <- fill_hex(rows$residual, 0.20, 0.96, LIGHT, VIOLET)
f3a <- ggplot(rows, aes(proxy_s, tissue_s, fill = fill_res)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", residual)), size = 3.0) +
  scale_fill_identity() +
  labs(x = "locked proxy", y = NULL, title = "Residual rho after additive") +
  theme(axis.text.x = element_text(size = 9))

# B: Treg Mantel rho vs locked proxies
f3b <- ggplot(treg_rows, aes(proxy_s, rho, fill = proxy_s)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", rho)), vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c(Brain = BLUE, PBMC = YELLOW, Fibro = VIOLET)) +
  coord_cartesian(ylim = c(0, max(treg_rows$rho, na.rm = TRUE) * 1.18)) +
  labs(x = "locked proxy", y = "Mantel rho",
       title = "Treg transfer rho vs locked proxies") +
  theme(axis.text.x = element_text(size = 9))

# C: BMMC frozen-panel gene/TF overlap fractions (estimand coverage, not ops)
bmc <- data.frame(
  axis = factor(c("Genes", "TFs"), levels = c("Genes", "TFs")),
  coverage = c(as.numeric(hon$bmmc$gene_coverage), as.numeric(hon$bmmc$tf_coverage)),
  lab = c(
    sprintf("%d/%d", hon$bmmc$gene_overlap, hon$bmmc$gene_panel),
    sprintf("%d/%d", hon$bmmc$tf_overlap, hon$bmmc$tf_panel)
  )
)
f3c <- ggplot(bmc, aes(axis, coverage)) +
  geom_col(fill = BLUE, width = 0.55) +
  geom_text(aes(label = sprintf("%s (%.3f)", lab, coverage)),
            vjust = -0.35, size = 2.8) +
  coord_cartesian(ylim = c(0, 1.12)) +
  labs(x = NULL, y = "panel overlap",
       title = "BMMC overlap with frozen 446 x 1200 panel")

# D: Treg additive fraction
f3d <- ggplot(treg_rows, aes(proxy_s, frac_add, fill = proxy_s)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", frac_add)), vjust = -0.35, size = 3.0) +
  scale_fill_manual(values = c(Brain = AQUA, PBMC = YELLOW, Fibro = VIOLET)) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = "locked proxy", y = "additive fraction",
       title = "Treg additive fraction of rho") +
  theme(axis.text.x = element_text(size = 9))

emit_ext("fig_ext3_honesty_policy",
         (f3a | f3b) / (f3c | f3d) + plot_annotation(tag_levels = "A") +
           plot_layout(axis_titles = "collect"),
         6.8, 5.6, tags = LETTERS[1:4])

# isolation guard: never leave fig_ext* under paper/figs
bad <- list.files(file.path(base, "paper", "figs"), pattern = "^fig_ext")
if (length(bad))
  stop("isolation broken: fig_ext* found under paper/figs: ", paste(bad, collapse = ","))

cat("make_figs_extension.R complete. Outputs in paper/figs_extension/\n")
cat("peerj_support_rows_touched remains false; no package rebuild.\n")
