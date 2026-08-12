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
f1a <- ggplot(rows, aes(proxy, tissue, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 3.2) +
  scale_fill_gradient(low = LIGHT, high = BLUE, limits = c(0.35, 0.92),
                      name = "Spearman $\\rho$") +
  labs(x = "locked proxy", y = NULL,
       title = "Construct SI Mantel vs locked $G_{\\mathrm{ATAC}}$") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 18, hjust = 1))

f1b <- ggplot(rows, aes(proxy, tissue, fill = frac_add)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", frac_add)), size = 3.1) +
  scale_fill_gradient(low = "#f0f0f0", high = AQUA, limits = c(0.4, 0.9),
                      name = "additive frac.") +
  labs(x = "locked proxy", y = NULL,
       title = "Fraction explained by additive marginals") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 18, hjust = 1))

rows$pair_lab <- paste(rows$tissue, rows$proxy, sep = " vs ")
# order by rho
rows$pair_lab <- factor(rows$pair_lab, levels = rows$pair_lab[order(rows$rho)])
f1c <- ggplot(rows, aes(y = pair_lab)) +
  geom_segment(aes(x = 0.35, xend = rho, yend = pair_lab), color = "grey70",
               linewidth = 0.5) +
  geom_point(aes(x = rho, color = tissue), size = 2.4) +
  scale_color_manual(values = c(`DESCARTES spleen` = VIOLET, BMMC = YELLOW,
                                `Treg pilot` = BLUE)) +
  labs(x = "observed Spearman $\\rho$", y = NULL, color = NULL,
       title = "Transfer strength (BMMC--PBMC stands out)") +
  theme(legend.position = "top")

peak_plot <- bind_rows(
  peak_df %>% mutate(tissue = factor(tissue, levels = levels(rows$tissue))),
  peerj_peaks %>% mutate(tissue = factor(tissue, levels = c("Brain", "PBMC", "Fibroblast mix")))
)
# two-panel: extension only for D to avoid scale confusion; PeerJ as ref annotation
f1d <- ggplot(peak_df, aes(reorder(tissue, relevant_peaks), relevant_peaks)) +
  geom_col(fill = BLUE, width = 0.65) +
  geom_text(aes(label = format(relevant_peaks, big.mark = ",")),
            hjust = -0.08, size = 3.0) +
  coord_flip(ylim = c(0, max(peak_df$relevant_peaks) * 1.18)) +
  labs(x = NULL, y = "motif-linked peaks (extension meta)",
       title = "Linkage peaks (extension); PeerJ ref: 6.6k--11.6k") +
  theme(plot.title = element_text(size = 10))

emit_ext("fig_ext1_construct_mantel",
         (f1a | f1b) / (f1c | f1d) + plot_annotation(tag_levels = "A") +
           plot_layout(axis_titles = "collect"),
         6.8, 5.8, tags = LETTERS[1:4])

# ============== E2: CollecTRI + baselines ==============
coll <- J("results/v2/extension/baselines/collectri_prior/summary.json")
cr <- coll$result
stopifnot(identical(cr$status, "projected_partial") || identical(cr$status, "emitter_ready") ||
            !is.null(cr$tf_coverage))
cov_df <- data.frame(
  metric = c("TF hit rate", "Gene hit rate", "Edge keep rate"),
  value = c(as.numeric(cr$tf_coverage),
            as.numeric(cr$gene_coverage),
            as.numeric(cr$n_edges_on_panel) / as.numeric(cr$n_edges_raw)),
  stringsAsFactors = FALSE
)
cov_df$metric <- factor(cov_df$metric, levels = cov_df$metric)
f2a <- ggplot(cov_df, aes(metric, value, fill = metric)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(BLUE, AQUA, YELLOW)) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = NULL, y = "fraction",
       title = "CollecTRI panel projection (edge keep / TF / gene hit rates)") +
  theme(axis.text.x = element_text(angle = 12, hjust = 1))

# baseline status as bar of "ready score" 1=ok 0=not
base_methods <- c("motif_only_rp", "degree_matched_random", "encode_chip_binding", "collectri_prior")
base_rows <- bind_rows(lapply(base_methods, function(mid) {
  s <- J(file.path("results/v2/extension/baselines", mid, "summary.json"))
  st <- as.character(s$result$status %||% s$status %||% "unknown")
  data.frame(method = mid, status = st, stringsAsFactors = FALSE)
}))
base_rows$label <- c(
  motif_only_rp = "motif-only RP",
  degree_matched_random = "degree-matched",
  encode_chip_binding = "ENCODE ChIP",
  collectri_prior = "CollecTRI prior"
)[base_rows$method]
base_rows$ready <- as.numeric(grepl("emit|project|ready|summary", base_rows$status, ignore.case = TRUE))
base_rows$label <- factor(base_rows$label, levels = rev(base_rows$label))
# tikzDevice cannot measure unescaped underscores in text labels
base_rows$status_tex <- gsub("_", " ", base_rows$status, fixed = TRUE)
f2b <- ggplot(base_rows, aes(label, ready, fill = factor(ready))) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = status_tex), vjust = -0.35, size = 2.7) +
  scale_fill_manual(values = c(`0` = RED, `1` = AQUA)) +
  scale_y_continuous(limits = c(0, 1.35), breaks = c(0, 1),
                     labels = c("fail", "ready")) +
  labs(x = NULL, y = NULL,
       title = "Tier-A/C baseline emitter status (extension)") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

edge_df <- data.frame(
  kind = factor(c("raw edges", "on panel"), levels = c("raw edges", "on panel")),
  n = c(as.numeric(cr$n_edges_raw), as.numeric(cr$n_edges_on_panel))
)
f2c <- ggplot(edge_df, aes(kind, n, fill = kind)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.4, size = 3.1) +
  scale_fill_manual(values = c(`raw edges` = MUTED, `on panel` = BLUE)) +
  labs(x = NULL, y = "edges",
       title = "Literature prior: raw vs frozen-panel support") +
  coord_cartesian(ylim = c(0, max(edge_df$n) * 1.12))

freeze_df <- data.frame(
  item = factor(c("Support rows", "BH families", "Panel TFs", "Panel genes"),
                levels = c("Support rows", "BH families", "Panel TFs", "Panel genes")),
  value = c(as.numeric(hon$freeze$peerj_support_rows), 8, 446, 1200)
)
f2d <- ggplot(freeze_df, aes(item, value, fill = item)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = value), vjust = -0.35, size = 3.1) +
  scale_fill_manual(values = c(BLUE, AQUA, YELLOW, VIOLET)) +
  labs(x = NULL, y = "count",
       title = "Freeze badge: Support=13; rows\\_touched=false") +
  coord_cartesian(ylim = c(0, max(freeze_df$value) * 1.15)) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

emit_ext("fig_ext2_baselines_collectri",
         (f2a | f2b) / (f2c | f2d) + plot_annotation(tag_levels = "A") +
           plot_layout(axis_titles = "collect"),
         6.8, 5.6, tags = LETTERS[1:4])

# ============== E3: honesty / policy ==============
htan_df <- data.frame(
  step = factor(c("Local tar", "Inventory", "Peak matrix", "SI outcome"),
                levels = c("Local tar", "Inventory", "Peak matrix", "SI outcome")),
  score = c(1, 1, 0, 0),
  state = c("present", "fragments-only", "absent", "blocked")
)
f3a <- ggplot(htan_df, aes(step, score, fill = factor(score))) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = state), vjust = -0.35, size = 2.9) +
  scale_fill_manual(values = c(`0` = RED, `1` = AQUA)) +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 1),
                     labels = c("no", "yes")) +
  labs(x = NULL, y = NULL,
       title = "HTAN D3 dual-path: blocked (fragments only)") +
  theme(axis.text.x = element_text(angle = 12, hjust = 1))

orphan_df <- data.frame(
  metric = factor(c("h5ad matrices", "with cell-type obs", "pilot types"),
                  levels = c("h5ad matrices", "with cell-type obs", "pilot types")),
  value = c(as.numeric(hon$orphan_lake$n_h5ad),
            as.numeric(hon$orphan_lake$n_with_celltype_obs),
            1)
)
f3b <- ggplot(orphan_df, aes(metric, value, fill = metric)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = value), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(BLUE, RED, YELLOW)) +
  labs(x = NULL, y = "count",
       title = "Orphan lake + GSM6449881 filename meta pilot") +
  coord_cartesian(ylim = c(0, max(orphan_df$value) * 1.15)) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

bmc <- data.frame(
  axis = factor(c("Genes", "TFs"), levels = c("Genes", "TFs")),
  coverage = c(as.numeric(hon$bmmc$gene_coverage), as.numeric(hon$bmmc$tf_coverage)),
  label = c(
    sprintf("%d/%d", hon$bmmc$gene_overlap, hon$bmmc$gene_panel),
    sprintf("%d/%d", hon$bmmc$tf_overlap, hon$bmmc$tf_panel)
  )
)
f3c <- ggplot(bmc, aes(axis, coverage)) +
  geom_col(fill = BLUE, width = 0.55) +
  geom_hline(yintercept = as.numeric(hon$bmmc$gate_threshold), color = RED,
             linetype = 2, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%s (%.3f)", label, coverage)),
            vjust = -0.35, size = 2.9) +
  annotate("text", x = 1.5, y = as.numeric(hon$bmmc$gate_threshold) + 0.04,
           label = "0.8967 / 0.7646 vs 0.90 gate", color = RED, size = 2.8) +
  coord_cartesian(ylim = c(0, 1.12)) +
  labs(x = NULL, y = "frozen-panel coverage",
       title = "BMMC P3: disclosed coverage (construct SI only)")

fr <- data.frame(
  item = factor(c("Support rows", "rows touched", "primary GATAC mutable", "ext tissues"),
                levels = c("Support rows", "rows touched", "primary GATAC mutable", "ext tissues")),
  value = c(13, 0, 0, 3)
)
f3d <- ggplot(fr, aes(item, value, fill = item)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = value), vjust = -0.35, size = 3.1) +
  scale_fill_manual(values = c(BLUE, AQUA, MUTED, YELLOW)) +
  labs(x = NULL, y = "count / flag",
       title = "Freeze: Support=13; peerj rows touched=false") +
  coord_cartesian(ylim = c(0, 15)) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

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
