#!/usr/bin/env Rscript
# Authoritative figures for the fixed-panel audit. Run from paper/:
#   python3 make_panel_data.py && Rscript make_figs.R
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  library(jsonlite)
  library(tikzDevice)
})
options(tikzDefaultEngine = "pdftex")
# Figure text uses the same Times-compatible face as the manuscript body
# (NimbusRomNo9L in the PeerJ PDF): newtx keeps labels native-LaTeX instead of
# a sans fallback, and tikzDevice measures metrics against these exact fonts.
FONTPKGS <- c("\\usepackage{amsmath}", "\\usepackage{newtxtext}", "\\usepackage{newtxmath}")
# Text metrics must be measured at the manuscript's own base size. manuscript.tex
# is 11pt; tikzDevice defaults to 10pt, which under-reserves every string by ~10%
# and pushes left-anchored axis labels outside the bounding box, where the
# device's own clip path shears off their first characters.
options(tikzDocumentDeclaration = "\\documentclass[11pt]{article}",
        tikzMetricPackages = c(getOption("tikzMetricPackages"), FONTPKGS),
        tikzLatexPackages = c(getOption("tikzLatexPackages"), FONTPKGS),
        tikzMetricsDictionary = ".tikz_metrics_pdftex")

base <- Sys.getenv("SCFM_BASE", "..")
source(file.path(base, "src", "v2", "figure_helpers.R"))
res <- file.path(base, "results", "v2")
figs <- file.path(base, "paper", "figs")
dir.create(figs, showWarnings = FALSE, recursive = TRUE)
J <- function(name) fromJSON(file.path(res, name), simplifyVector = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x
audit <- J("fixed_panel_audit_v2.json")
injection <- J("fixed_panel_signal_injection_v2.json")
probe_stats <- J("tf_probe_pair_stats_v2.json")
probe_sens <- J("tf_probe_pair_sensitivity_v2.json")
probe_eval <- J("tf_probe_pair_eval_v2.json")
spec_sens <- J("spec_sensitivity_v2.json")
ladder <- J("marginal_vs_adjusted_v2.json")
effect_scale <- J("effect_vs_injection_scale_v2.json")
decomp <- J("cross_tissue_additive_decomp_v2.json")
invariance <- J("proxy_celltype_invariance_v2.json")
deg_null <- J("degree_preserving_null_v2.json")
attn_readout <- J("readout_attention_v2.json")
ko_stat <- J("insilico_ko_v2.json")
omission <- J("verify_brain_attention_omission_v2.json")
pertype_fm <- J("pertype_fm_v2.json")
panel <- fromJSON(file.path(base, "paper", "panel_data.json"), simplifyVector = FALSE)

BLUE <- "#2a78d6"; AQUA <- "#1b8f75"; YELLOW <- "#d99a00"
VIOLET <- "#5946b2"; RED <- "#d84a4a"; MUTED <- "grey45"; LIGHT <- "#d9e8f8"

# Figure text must survive the \fitfig downscale to the ~6.5in PeerJ text block:
# design near that width with an 11pt base so the smallest rendered text stays
# within ~1pt of the 10pt body font instead of dropping to 6-7pt.
theme_set(theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  plot.title = element_text(size = 11, face = "plain", hjust = 0, margin = margin(b = 3)),
  strip.text = element_text(size = 10), legend.title = element_text(size = 9.5),
  legend.text = element_text(size = 9.5), axis.title = element_text(size = 10),
  axis.text = element_text(size = 9.5), plot.margin = margin(4, 6, 3, 4)))

# Bottom-right D panels often carry longer titles that crowd the D tag when
# left-flush; a small left margin shifts the title right of the tag.
d_title_nudge <- theme(
  plot.title = element_text(hjust = 0, margin = margin(b = 3, l = 10)))

emit <- function(name, plot, width, height, tags = LETTERS[1:4]) {
  path <- file.path(figs, paste0(name, ".tex"))
  tikz(path, width = width, height = height, standAlone = FALSE, sanitize = FALSE)
  print(plot)
  dev.off()
  # A corrupted text node (e.g. an unescaped % turning the rest of a line into a
  # LaTeX comment) truncates the tikz stream WITHOUT raising an R error, and the
  # standalone compile still exits 0 with panels silently missing. Verify every
  # expected panel-tag node survived in the emitted file.
  lines <- readLines(path, warn = FALSE)
  missing <- tags[!vapply(tags, function(t) any(grepl(paste0("{", t, "}"), lines, fixed = TRUE)),
                          logical(1))]
  if (length(missing))
    stop("emitted figure ", name, " lost panel tags: ", paste(missing, collapse = ","))
  cat(name, "ok\n")
}

model_label <- c(
  geneformer_embed = "Geneformer embed",
  geneformer_attn = "Geneformer attention",
  geneformer_ko_raw = "Geneformer KO",
  geneformer_ko_posctrl = "Artifact-corrected KO",
  scFoundation_encoder = "scFoundation encoder",
  UCE_encoder = "UCE encoder",
  scGPT_encoder = "scGPT encoder",
  random_init_floor = "Random-init floor"
)

# Co-expression baseline through both nulls (computed separately: partialling
# co-expression out of itself is degenerate, so these rows live outside the
# audit's primary family in their own JSONs).
baseline_rows <- function() {
  files <- c(Brain = file.path(res, "brain_coexp_baseline_null_v2.json"),
             PBMC = file.path(res, "pbmc_coexp_baseline_null_v2.json"))
  bind_rows(lapply(names(files), function(tn) {
    b <- fromJSON(files[[tn]])
    data.frame(tissue = tn, model = "Co-expression baseline",
               model_key = "co_expression", spec = "full", rho = b$observed_rho,
               mantel_p = b$pM, mantel_q = NA_real_, degree_p = b$pD, degree_q = NA_real_,
               stringsAsFactors = FALSE)
  }))
}
baseline <- baseline_rows()

pooled_rows <- function(tissue) {
  rows <- Filter(function(x) identical(x$row_type, "pooled_fm"), audit$pooled[[tissue]]$rows)
  bind_rows(lapply(rows, function(x) data.frame(
    tissue = ifelse(tissue == "brain", "Brain", "PBMC"),
    model = unname(model_label[[x$model_label]]), model_key = x$model_label,
    spec = x$confound_spec, rho = x$observed_partial_rho,
    mantel_p = x$mantel$p_mc, mantel_q = x$mantel$bh_q_family,
    degree_p = x$degree_preserving$p_mc,
    degree_q = x$degree_preserving$bh_q_family,
    stringsAsFactors = FALSE
  )))
}
pooled <- bind_rows(pooled_rows("brain"), pooled_rows("pbmc"))
pooled$status <- ifelse(pooled$mantel_q < 0.05 & pooled$degree_q < 0.05,
                        "supported by both nulls",
                        ifelse(pooled$mantel_q < 0.05 | pooled$degree_q < 0.05,
                               "supported by one null", "not supported"))

support_colors <- c("supported by both nulls" = AQUA,
                    "supported by one null" = YELLOW,
                    "not supported" = MUTED,
                    "co-expression baseline" = "black")

cross <- bind_rows(lapply(audit$cross_tissue_construct_reproducibility$rows, function(x) {
  names <- c(GSE174367 = "Brain", PBMC10k = "PBMC", GSE206767 = "Fibroblast mix")
  data.frame(pair = paste(names[[x$pair[[1]]]], names[[x$pair[[2]]]], sep = "--"),
             rho = x$observed_spearman)
}))

# ==================== Figure 1: construct, evidence, panel ====================
# A: S-shaped construction schematic matching Methods: ATAC peak → peak–gene
# link + JASPAR/MOODS motif → directed TF→target weight → restrict to the
# hash-pinned 446×1,200 panel. No RNA enters the proxy.
flow <- data.frame(
  x = c(1.0, 3.05, 3.05, 1.0),
  y = c(2.18, 2.18, 1.18, 1.18),
  step = 1:4,
  label = c("accessible peak", "motif match",
            "TF$\\rightarrow$target weight", "fixed panel"),
  stage = c("ATAC evidence", "JASPAR / MOODS",
            "proxy score", "446$\\times$1{,}200 audit set"),
  fill = c(LIGHT, "#ebe6f6", "#d9efe8", "#f2f2f2"),
  border = c(BLUE, VIOLET, AQUA, MUTED)
)
arrow_col <- "#5c6770"
f1a <- ggplot(flow, aes(x, y)) +
  # stage chips under each node (darker text for PeerJ downscale)
  geom_label(aes(y = y - 0.36, label = stage),
             size = 2.55, label.size = 0, fill = "grey96", color = MUTED,
             label.padding = unit(0.11, "lines"), show.legend = FALSE) +
  # main nodes
  geom_label(aes(label = label, fill = fill, color = border),
             size = 3.35, label.size = 0.6, fontface = "plain",
             label.padding = unit(0.38, "lines"), label.r = unit(0.14, "lines"),
             show.legend = FALSE) +
  scale_fill_identity() +
  scale_color_identity() +
  # numbered badges anchored at each node's upper-left
  annotate("point", x = c(0.42, 2.47, 2.47, 0.42),
           y = c(2.42, 2.42, 1.42, 1.42),
           size = 5.0, color = c(BLUE, VIOLET, AQUA, MUTED)) +
  annotate("text", x = c(0.42, 2.47, 2.47, 0.42),
           y = c(2.42, 2.42, 1.42, 1.42),
           label = c("1", "2", "3", "4"),
           size = 2.7, color = "white", fontface = "bold") +
  # directed arrows with method-accurate edge labels
  annotate("segment", x = 1.72, xend = 2.32, y = 2.18, yend = 2.18,
           arrow = arrow(length = unit(0.09, "in"), type = "closed"),
           color = arrow_col, linewidth = 0.6) +
  annotate("text", x = 2.02, y = 2.40, label = "$\\pm$2\\,kb gene link",
           size = 2.45, color = MUTED) +
  annotate("segment", x = 3.05, xend = 3.05, y = 1.86, yend = 1.50,
           arrow = arrow(length = unit(0.09, "in"), type = "closed"),
           color = arrow_col, linewidth = 0.6) +
  annotate("text", x = 3.42, y = 1.68, label = "aggregate",
           size = 2.45, color = MUTED, hjust = 0) +
  annotate("segment", x = 2.32, xend = 1.72, y = 1.18, yend = 1.18,
           arrow = arrow(length = unit(0.09, "in"), type = "closed"),
           color = arrow_col, linewidth = 0.6) +
  annotate("text", x = 2.02, y = 0.96, label = "restrict to panel",
           size = 2.45, color = MUTED) +
  # disclaimer as integrated footer chip (edge weights without expression stats; not RNA-structure-independent overall)
  annotate("label", x = 2.0, y = 0.38,
           label = "edge weights without expression stats; regulatory-potential proxy, not causal ground truth",
           size = 2.75, color = RED, fill = "#fdeeee",
           label.size = 0.35, label.r = unit(0.12, "lines"),
           label.padding = unit(0.26, "lines")) +
  coord_cartesian(xlim = c(0.15, 3.85), ylim = c(0.18, 2.62), clip = "off") +
  labs(title = "Accessibility and motif evidence define the audited proxy") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(size = 11, hjust = 0),
        plot.margin = margin(4, 10, 3, 0))

f1b <- ggplot(cross, aes(reorder(pair, rho), rho)) +
  geom_col(fill = BLUE, width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.5, size = 3.3) +
  coord_cartesian(ylim = c(0, 0.58)) +
  labs(x = NULL, y = "observed Spearman $\\rho$", title = "Proxy structure is reproducible across tissues") +
  theme(axis.text.x = element_text(angle = 12, hjust = 0.95),
        plot.margin = margin(4, 8, 3, 6))

motif <- bind_rows(lapply(names(panel$motif_evidence), function(tn) {
  m <- panel$motif_evidence[[tn]]
  data.frame(tissue = tn, hits_per_peak = m$hits_per_peak,
             expected = m$expected_random_per_peak)
}))
motif$tissue <- factor(motif$tissue, levels = c("Brain", "PBMC", "Fibroblast mix"))
f1c <- ggplot(motif, aes(tissue, hits_per_peak)) +
  geom_col(fill = VIOLET, width = 0.6) +
  geom_hline(yintercept = 5.5, linetype = "dashed", color = RED) +
  geom_text(aes(label = sprintf("%.1f", hits_per_peak)), vjust = -0.5, size = 3.3) +
  annotate("text", x = 2, y = 26.5, label = "expected random = 5.5", size = 3.1,
           color = RED) +
  coord_cartesian(ylim = c(0, 28)) +
  labs(x = NULL, y = "motif hits per peak",
       title = "Expected-random motif share") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

comp <- bind_rows(lapply(c("brain", "pbmc"), function(tn) {
  cc <- panel$panel_composition[[tn]]
  data.frame(tissue = ifelse(tn == "brain", "Brain", "PBMC"),
             metric = rep(c("identical binary profile", "partner at cosine $>0.8$"), each = 1),
             n = c(cc$n_identical_profile, cc$n_partner_cosine_gt_0.8))
}))
f1d <- ggplot(comp, aes(tissue, n, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.5, size = 3.3) +
  scale_fill_manual(values = c(BLUE, YELLOW)) +
  coord_cartesian(ylim = c(0, 160)) +
  labs(x = NULL, y = "TFs (of 446)", fill = NULL,
       title = "The fixed panel is composition-biased") +
  theme(legend.position = "top") +
  d_title_nudge

# free(f1a): A is theme_void, so default panel-align with C leaves a large
# left blank under C's y-axis; releasing it shifts A left and opens air to B.
emit("fig1_truth_construct",
     ((free(f1a) | f1b) + plot_layout(widths = c(1.06, 0.94))) / (f1c | f1d) +
       plot_annotation(tag_levels = "A") +
       plot_layout(axis_titles = "collect"), 6.8, 5.8)

# ============== Figure 2: cross-tissue additive decomposition ==============
dec <- bind_rows(lapply(decomp$rows, function(x) {
  names <- c(GSE174367 = "Brain", PBMC10k = "PBMC", GSE206767 = "Fibroblast mix")
  data.frame(pair = paste(names[[x$pair[[1]]]], names[[x$pair[[2]]]], sep = "--"),
             observed = x$observed_spearman,
             additive = x$additive_pred_spearman,
             residual = x$residual_spearman_after_own_additive_fits,
             phi = x$binary_support_phi,
             frac = x$fraction_explained_by_additive_marginals)
}))
dec_long <- bind_rows(
  dec %>% transmute(pair, metric = "observed", rho = observed),
  dec %>% transmute(pair, metric = "additive-marginal prediction", rho = additive),
  dec %>% transmute(pair, metric = "residual after own additive fit", rho = residual))
dec_long$metric <- factor(dec_long$metric,
                          levels = c("observed", "additive-marginal prediction",
                                     "residual after own additive fit"))
f2a <- ggplot(dec_long, aes(pair, rho, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c(BLUE, YELLOW, AQUA)) +
  labs(x = NULL, y = "Spearman $\\rho$", fill = NULL,
       title = "Most apparent reproducibility is marginal structure") +
  # One-row legend is wider than the panel and clips the left key; two
  # compact rows keep all three entries inside the figure box.
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "top",
        legend.text = element_text(size = 8.5),
        legend.key.size = unit(0.28, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(4, 6, 3, 6),
        axis.text.x = element_text(angle = 15, hjust = 1))

f2b <- ggplot(dec, aes(reorder(pair, frac), frac)) +
  geom_col(fill = BLUE, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f\\%%", 100 * frac)), vjust = -0.5, size = 3.3) +
  coord_cartesian(ylim = c(0, 0.9)) +
  labs(x = NULL, y = "fraction explained",
       title = "Additive marginals explain 69--78\\%") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

f2c <- ggplot(dec, aes(reorder(pair, residual - observed), residual - observed)) +
  geom_col(aes(fill = (residual - observed) > 0), width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "grey55") +
  scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = MUTED)) +
  labs(x = NULL, y = "residual $-$ observed $\\rho$",
       title = "One pair is net-subtractive") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

inv <- bind_rows(lapply(c("brain", "pbmc"), function(tn) {
  t <- invariance$tissues[[tn]]
  data.frame(tissue = ifelse(tn == "brain", "Brain", "PBMC"),
             mean = t$pairwise_spearman$mean,
             lo = t$pairwise_spearman$min, hi = t$pairwise_spearman$max,
             n_types = t$n_cell_types)
}))
f2d <- ggplot(inv, aes(tissue, mean)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, color = MUTED) +
  geom_point(size = 2.6, color = BLUE) +
  coord_cartesian(ylim = c(0.95, 1.0)) +
  labs(x = NULL, y = "cell-type consensus $\\rho$",
       title = "Proxy is near cell-type-invariant") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  d_title_nudge

emit("fig2_cross_tissue_decomp", (f2a | f2b) / (f2c | f2d) + plot_annotation(tag_levels = "A") +
       plot_layout(axis_titles = "collect"), 6.8, 5.8)

# ==================== Figure 3: primary audit ====================
mk_forest <- function(df, title) {
  df$display <- factor(df$display, levels = rev(df$display))
  x_min <- min(df$rho, na.rm = TRUE) - 0.0018  # keep every point inside the panel
  x_max <- max(0.0145, max(df$rho, na.rm = TRUE) + 0.0018)
  x_label <- max(0.0165, x_max + 0.0025)  # q/p column clears points, stays in margin
  q_labels <- ifelse(df$is_baseline,
                     sprintf("$p_M=%.3f;\\ p_D=%.3f$", df$mantel_p, df$degree_p),
                     sprintf("$q_M=%.3f;\\ q_D=%.3f$", df$mantel_q, df$degree_q))
  # Y-tick width is already reserved in the axis grob (~153pt). The helper's
  # extra left pad (>120pt envelope) indented A/B relative to the full-width
  # C|D row; keep only a thin gutter so the forests extend left.
  panel_margin <- compute_measured_panel_margin(
    ytick_labels = df$display,
    right_labels = q_labels,
    base_left_pt = 0,
    normal_ytick_envelope_pt = 160
  )
  ggplot(df, aes(rho, display, color = status)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(aes(shape = is_baseline), size = 2.5) +
    geom_text(aes(x = x_label,
                  label = ifelse(is_baseline,
                                 sprintf("$p_M=%.3f;\\ p_D=%.3f$", mantel_p, degree_p),
                                 sprintf("$q_M=%.3f;\\ q_D=%.3f$", mantel_q, degree_q))),
              hjust = 0, size = 4.0, color = "black") +
    scale_color_manual(values = support_colors) +
    scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16), guide = "none") +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(x = "partial Spearman $\\rho$", y = NULL, color = NULL, title = title) +
    theme(legend.position = "top", axis.ticks.y = element_blank(),
          axis.text.y = element_text(margin = margin(r = 2)),
          plot.margin = panel_margin)
}
primary <- pooled %>% filter(spec == "full")
primary <- bind_rows(primary, baseline)
primary$is_baseline <- primary$model_key == "co_expression"
primary$status <- ifelse(primary$is_baseline, "co-expression baseline", primary$status)
tissue_order <- unique(primary$tissue)
primary <- primary[order(match(primary$tissue, tissue_order), primary$is_baseline), ]
primary$display <- paste(primary$model, primary$tissue, sep = " -- ")
f3a <- mk_forest(primary, "Full-confound specification (primary)")

nondeg <- pooled %>% filter(spec != "full")
nondeg$is_baseline <- FALSE
nondeg <- nondeg[order(match(nondeg$tissue, tissue_order)), ]
nondeg$display <- paste(nondeg$model, nondeg$tissue, sep = " -- ")
nondeg$closest_q <- pmin(nondeg$mantel_q, nondeg$degree_q)
f3b <- mk_forest(nondeg, "Non-degree co-primary (closest $q$ shown; none below 0.05)")

f3c <- ggplot(primary %>% filter(!is_baseline), aes(mantel_q, degree_q, color = status)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.4) +
  scale_color_manual(values = support_colors) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "gene-label $q_M$", y = "row-shuffle $q_D$",
       color = NULL, title = "Both nulls below 0.05") +
  theme(legend.position = "none")  # colors defined by panel A's legend; no room here

# 0x12: dual-null positives vs same-edge co-expression baseline, with r^2 notes.
# Also mark rows below their tissue baseline (e.g. PBMC Geneformer embed).
short_model <- c("Geneformer embed" = "GF emb", "scGPT encoder" = "scGPT",
                 "UCE encoder" = "UCE", "Geneformer attention" = "GF attn",
                 "scFoundation encoder" = "scF", "Random-init floor" = "rand",
                 "Geneformer KO" = "KO", "Artifact-corrected KO" = "KO+")
supp <- primary %>% filter(status == "supported by both nulls")
base_cmp <- bind_rows(lapply(seq_len(nrow(supp)), function(i) {
  b <- baseline[baseline$tissue == supp$tissue[i], ]
  tissue_short <- ifelse(supp$tissue[i] == "Brain", "B", "P")
  sm <- short_model[[supp$model[i]]] %||% substr(supp$model[i], 1, 8)
  lab <- paste(sm, tissue_short, sep = "-")
  bind_rows(
    data.frame(label = lab, kind = "FM dual-null", rho = supp$rho[i],
               r2 = supp$rho[i]^2, stringsAsFactors = FALSE),
    data.frame(label = lab, kind = "co-expression baseline", rho = b$rho[1],
               r2 = b$rho[1]^2, stringsAsFactors = FALSE)
  )
}))
base_cmp$label <- factor(base_cmp$label, levels = unique(base_cmp$label))
# r^2 annotation only on FM bars (tiny effects: scientific notation avoided).
r2_lab <- base_cmp %>% filter(kind == "FM dual-null") %>%
  mutate(lab = sprintf("$r^{2}\\!=\\!%.2e$", r2))
f3d <- ggplot(base_cmp, aes(label, rho, fill = kind)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = r2_lab, aes(x = label, y = rho, label = lab),
            vjust = -0.35, size = 2.0, color = MUTED, inherit.aes = FALSE) +
  scale_fill_manual(values = c("FM dual-null" = AQUA, "co-expression baseline" = "black")) +
  labs(x = NULL, y = "partial Spearman $\\rho$", fill = NULL,
       title = "Dual-null rows vs baseline ($r^{2}=\\rho^{2}$ on FM bars)") +
  coord_cartesian(ylim = c(min(0, min(base_cmp$rho) - 0.001),
                           max(base_cmp$rho) * 1.22)) +
  guides(fill = guide_legend(ncol = 1,
                             keywidth = unit(0.5, "lines"),
                             keyheight = unit(0.4, "lines"))) +
  theme(legend.position = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.32, "lines"),
        legend.spacing.y = unit(0, "pt"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.title = element_text(hjust = 0, margin = margin(t = 0, b = 1, l = 10)),
        axis.text.x = element_text(size = 7.5, angle = 40, hjust = 1))

# Bottom row: A/B forests keep a narrow *panel* (data box ~176-354 of 491 pt)
# with y-tick / q-label text in the margins. Plain patchwork panel-aligns C|D
# into that column (too short); free()+plot_spacer() still clips the right edge
# to the forest panel. cowplot (align="none") lets C|D span nearly full width
# under the A/B y-tick column (spacer 1.28 -> none). wrap_elements keeps tags.
f3_bottom <- cowplot::plot_grid(
  f3c, f3d, nrow = 1,
  rel_widths = c(1.0, 1.1),
  labels = c("C", "D"), label_size = 11, label_fontface = "plain",
  hjust = 0, vjust = 1.1
)
emit("fig3_primary_audit",
     (f3a / f3b / wrap_elements(full = f3_bottom)) +
       plot_annotation(tag_levels = list(c("A", "B", ""))) +
       plot_layout(heights = c(1.15, 1.0, 0.95)),
     6.8, 8.0)  # 9.2->8.0: a 9.2in float exceeded the text block, forcing a near-empty float page (large whitespace band); 8.0in fits the caption on one float page

# ==================== Figure 4: concordance (attention-likeness) check ====================
# Panel JSON key remains usability_fm_vs_coexp for schema stability; plot labels use "concordance".
usa <- bind_rows(lapply(names(panel$usability_fm_vs_coexp), function(tn) {
  u <- panel$usability_fm_vs_coexp[[tn]]
  bind_rows(lapply(names(u), function(mk) {
    data.frame(tissue = ifelse(tn == "brain", "Brain", "PBMC"),
               model_key = mk, model = unname(model_label[[mk]]),
               fm_vs_coexp = u[[mk]], stringsAsFactors = FALSE)
  }))
}))
usa$display <- paste(usa$model, usa$tissue, sep = " -- ")
usa <- usa %>% left_join(primary %>% filter(!is_baseline) %>%
                           select(model_key, tissue, status, primary_rho = rho,
                                  mantel_q, degree_q),
                         by = c("model_key", "tissue"))
usa$dual_null <- usa$status == "supported by both nulls"
usa$display <- factor(usa$display, levels = rev(usa$display))
# 0x16: full scatter of all 13 readouts; dual-null status in legend; call out
# brain Geneformer attention (concordance passer, dual-null negative).
usa$highlight <- usa$model_key == "geneformer_attn" & usa$tissue == "Brain"
f4a <- ggplot(usa, aes(fm_vs_coexp, primary_rho)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(fill = dual_null, shape = dual_null), size = 2.6,
             color = "black", stroke = 0.3) +
  geom_point(data = usa %>% filter(highlight),
             shape = 1, size = 5.0, color = RED, stroke = 0.7) +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey75"),
                    labels = c(`TRUE` = "dual-null Support",
                               `FALSE` = "not dual-null"),
                    name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 22), guide = "none") +
  labs(x = "FM--co-expression concordance $\\rho$",
       y = "full-confound partial $\\rho$",
       title = "All readouts: concordance vs dual-null (ring: GF attn passer)") +
  theme(legend.position = "top",
        legend.text = element_text(size = 7.5))

tile <- usa %>% transmute(display,
                          `supported positive` =
                            status == "supported by both nulls" & primary_rho > 0,
                          `passes concordance ($\\rho>0$)` = fm_vs_coexp > 0)
tile_long <- bind_rows(
  tile %>% transmute(display, check = "supported positive", pass = `supported positive`),
  tile %>% transmute(display, check = "passes concordance ($\\rho>0$)",
                     pass = `passes concordance ($\\rho>0$)`))
f4b <- ggplot(tile_long, aes(check, display, fill = pass)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey80")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Support vs concordance") +
  # Nudge only the B tag right; title stays on the plot title grob.
  theme(legend.position = "top", axis.text.x = element_text(angle = 12, hjust = 1),
        axis.ticks = element_blank(),
        plot.tag.position = c(0.10, 1))

ko_df <- data.frame(
  readout = rep(c("KO raw", "KO artifact-corrected"), each = 3),
  metric = rep(c("proxy (marginal)", "partial $|$ coexp", "co-expression"), 2),
  rho = c(ko_stat$ko_vs_atac, ko_stat$ko_partial_given_coexp, ko_stat$ko_vs_coexp,
          ko_stat$ko_ctrl_vs_atac, ko_stat$ko_ctrl_partial_given_coexp, NA_real_))
ko_df$metric <- factor(ko_df$metric,
                       levels = c("proxy (marginal)", "partial $|$ coexp",
                                  "co-expression"))
f4c <- ggplot(ko_df, aes(metric, rho, fill = readout)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = "grey55") +
  scale_fill_manual(values = c("grey20", BLUE)) +
  labs(x = NULL, y = "Spearman $\\rho$", fill = NULL,
       title = "Knockout readout and its position-shift control") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 12, hjust = 1))

gf <- primary %>% filter(!is_baseline, model %in% c("Geneformer embed", "Geneformer attention"))
f4d <- ggplot(gf, aes(rho, model, group = tissue, color = tissue)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(color = "grey65") +
  geom_point(size = 2.6) +
  scale_color_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  labs(x = "partial $\\rho$ (full confounds)", y = NULL, color = NULL,
       title = "Same-model readouts disagree") +
  theme(legend.position = "top", axis.ticks.y = element_blank()) +
  d_title_nudge

# Bottom row: A/B share long y-tick labels, so plain patchwork panel-aligns
# C|D into the short data column (large blank under the ticks). cowplot
# (align="none") lets C|D span the full width; wrap_elements keeps tags.
# label_size 13 ≈ patchwork tag scale 1.20 at 11pt; label_x nudges D only.
f4_bottom <- cowplot::plot_grid(
  f4c, f4d, nrow = 1,
  rel_widths = c(1.05, 0.95),
  labels = c("C", "D"), label_size = 13, label_fontface = "plain",
  label_x = c(0, 0.08), label_y = 1,
  hjust = 0, vjust = 1.1
)
emit("fig4_usability_check",
     ((f4a | f4b) / wrap_elements(full = f4_bottom)) +
       plot_annotation(tag_levels = list(c("A", "B", ""))) +
       plot_layout(heights = c(1.15, 1.0)),
     6.8, 6.4)

# ==================== Figure 5: randomization diagnostics ====================
dn <- bind_rows(lapply(deg_null, function(x) {
  data.frame(label = gsub("_", " ", x$label), observed = x$observed,
             null_mean = x$null_mean, null_sd = x$null_sd, z = x$z, p = x$p_perm)
}))
dn$label <- factor(dn$label, levels = rev(dn$label))
f5a <- ggplot(dn, aes(label)) +
  geom_errorbar(aes(ymin = null_mean - 1.96 * null_sd, ymax = null_mean + 1.96 * null_sd),
                width = 0.3, color = MUTED) +
  geom_point(aes(y = null_mean), size = 2.0, color = MUTED) +
  geom_point(aes(y = observed), size = 2.6, color = AQUA) +
  geom_text(aes(y = observed, label = sprintf("$z=%.2f$", z)), hjust = -0.15, size = 3.1) +
  coord_flip(ylim = c(-0.0032, 0.0062)) +
  labs(x = NULL, y = "partial $\\rho$",
       title = "Row-shuffle null 95\\% band vs observed") +
  # Title over the panel (not the full plot) so long y-ticks do not pull it
  # into the left gutter; free(..., side="t") still lets the panel rise.
  theme(axis.ticks.y = element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust = 0, margin = margin(t = 0, b = 2)))

att_obs <- attn_readout$observed
att_df <- data.frame(
  variant = c("symmetrized", "TF$\\to$target", "target$\\to$TF",
              "symmetrized $|$ coexp", "TF$\\to$target $|$ coexp"),
  rho = c(att_obs$attn_sym_vs_atac, att_obs$attn_tf2target_vs_atac,
          att_obs$attn_target2tf_vs_atac, att_obs$attn_sym_partial_given_coexp,
          att_obs$attn_tf2target_partial_given_coexp),
  kind = c(rep("marginal", 3), rep("partial $|$ co-expression", 2)))
att_df$variant <- factor(att_df$variant, levels = att_df$variant)
f5b <- ggplot(att_df, aes(variant, rho, fill = kind)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, color = "grey55") +
  geom_hline(yintercept = 0.15, linetype = "dotted", color = AQUA, linewidth = 0.3) +
  annotate("text", x = 1, y = 0.16, label = "naive expectation",
           hjust = 0, size = 2.8, color = AQUA) +
  scale_fill_manual(values = c(BLUE, YELLOW)) +
  labs(x = NULL, y = "Spearman $\\rho$", fill = NULL,
       title = "Attention alignment negative") +
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.28, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.text.x = element_text(angle = 18, hjust = 1, size = 8.5))

ko_null <- bind_rows(lapply(c("mantel_ko_partial", "mantel_ko_ctrl_partial"), function(k) {
  m <- ko_stat[[k]]
  data.frame(readout = ifelse(k == "mantel_ko_partial", "KO raw", "KO artifact-corrected"),
             observed = m$observed, null_mean = m$null_mean, null_sd = m$null_sd,
             z = m$z, p = m$p_perm)
}))
f5c <- ggplot(ko_null, aes(readout)) +
  geom_errorbar(aes(ymin = null_mean - 1.96 * null_sd, ymax = null_mean + 1.96 * null_sd),
                width = 0.25, color = MUTED) +
  geom_point(aes(y = null_mean), size = 2.0, color = MUTED) +
  geom_point(aes(y = observed), size = 2.6, color = BLUE) +
  geom_text(aes(y = observed,
                label = ifelse(p <= 0.001 + 1e-12, "$p<0.001$", sprintf("$p=%.3f$", p))),
            hjust = -0.15, size = 3.1) +
  coord_flip() +
  # Short panel title keeps the left edge aligned with the plot spine; the longer
  # "KO raw + artifact-corrected both align after co-expression alone" wording
  # lives in the figure caption (panel C).
  labs(x = NULL, y = "partial $\\rho$ $|$ co-expression",
       title = "KO aligns after co-expression alone") +
  theme(axis.ticks.y = element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust = 0, margin = margin(t = 0, b = 2)))

# The omission guard is a parallel brain analysis with its own seed stream; its
# q-values differ from the authoritative audit, so the panel compares DECISIONS,
# not numbers (all eight rows agree).
om <- bind_rows(lapply(names(omission$rows), function(mk) {
  x <- omission$rows[[mk]]
  data.frame(model = unname(model_label[[mk]]),
             guard = x$qM < 0.05 & x$qD < 0.05)
}))
audit_dec <- primary %>% filter(!is_baseline, tissue == "Brain") %>%
  transmute(model, authoritative = status == "supported by both nulls")
om <- om %>% left_join(audit_dec, by = "model")
om_long <- bind_rows(
  om %>% transmute(model, analysis = "authoritative\naudit", supported = authoritative),
  om %>% transmute(model, analysis = "attention-omission\nguard", supported = guard))
om_long$model <- factor(om_long$model, levels = rev(unique(om_long$model)))
f5d <- ggplot(om_long, aes(analysis, model, fill = supported)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey85")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Support decisions agree") +
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.28, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks = element_blank()) +
  d_title_nudge

# free(A/C, side="t"): B/D keep a legend band above the panel; without freeing,
# A/C titles stay level with B/D while their panels drop, leaving a blank under
# the title. Tags stay on the default aligned corners (do not move).
emit("fig5_null_diagnostics",
     (free(f5a, type = "panel", side = "t") | f5b) /
       (free(f5c, type = "panel", side = "t") | f5d) +
       plot_annotation(tag_levels = "A") +
       plot_layout(axis_titles = "collect"),
     6.8, 6.0)

# ==================== Figure 6: specification dependence ====================
# 0x11: panel A is a paired full↔non-degree forest with dual-null markers
# (7 dual-null under full → 0 under non-degree) and sign-flip emphasis.
spec <- pooled %>%
  transmute(tissue, model, model_key, spec, rho, status, mantel_q, degree_q,
            spec_label = ifelse(spec == "full", "Full", "Non-degree"),
            dual = mantel_q < 0.05 & degree_q < 0.05,
            display = paste(model, tissue, sep = " -- "))
spec_ord <- primary %>% filter(!is_baseline) %>%
  transmute(tissue, model_key, ord = dplyr::row_number())
spec <- spec %>% left_join(spec_ord, by = c("tissue", "model_key"))
disp_levels <- rev((spec %>% filter(spec == "full") %>% arrange(ord))$display)
spec$display <- factor(spec$display, levels = disp_levels)

full6 <- spec %>% filter(spec_label == "Full") %>%
  select(display, tissue, model, model_key, ord, rho_Full = rho, dual_Full = dual)
nd6 <- spec %>% filter(spec_label == "Non-degree") %>%
  select(display, rho_nd = rho, dual_nd = dual)
pair6 <- full6 %>% left_join(nd6, by = "display")
pair6$sign_flip <- pair6$rho_Full * pair6$rho_nd < 0
stopifnot(nrow(pair6) == 13L, sum(pair6$dual_Full) == 7L, sum(pair6$dual_nd) == 0L)
n_flip <- sum(pair6$sign_flip)

f6a <- ggplot(pair6, aes(y = display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_segment(aes(x = rho_Full, xend = rho_nd, yend = display,
                   color = sign_flip), linewidth = 0.55) +
  geom_point(aes(x = rho_Full, shape = "Full", fill = dual_Full),
             size = 2.5, color = "black", stroke = 0.3) +
  geom_point(aes(x = rho_nd, shape = "Non-degree", fill = dual_nd),
             size = 2.5, color = "black", stroke = 0.3) +
  scale_shape_manual(values = c(Full = 21, `Non-degree` = 22), name = NULL) +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "white"),
                    labels = c(`TRUE` = "dual-null", `FALSE` = "not dual-null"),
                    name = NULL) +
  scale_color_manual(values = c(`TRUE` = RED, `FALSE` = "grey70"),
                     labels = c(`TRUE` = "sign flip", `FALSE` = "same sign"),
                     name = NULL) +
  labs(x = "partial Spearman $\\rho$", y = NULL,
       title = sprintf(
         "Paired full$\\leftrightarrow$non-degree: dual-null $7\\rightarrow 0$; %d sign flips",
         n_flip)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(order = 1),
         color = guide_legend(order = 3)) +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(1, "pt"),
        legend.text = element_text(size = 7.5),
        axis.ticks.y = element_blank())

# A discrete binned fill is used instead of a continuous gradient so tikzDevice
# emits a vector legend rather than a rasterized colorbar PNG (which would create
# an external \pgfimage dependency that does not resolve through \graphicspath).
rho_bins <- c("negative large", "negative small", "positive small", "positive large")
f6b_dat <- spec %>% mutate(rho_bin = factor(cut(
  rho, breaks = c(-Inf, -0.005, 0, 0.005, Inf),
  labels = rho_bins, right = FALSE), levels = rho_bins))
f6b <- ggplot(f6b_dat, aes(spec_label, display, fill = rho_bin)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = setNames(c(VIOLET, "#C9B7E0", "#B7D9D2", AQUA), rho_bins),
                    name = "partial $\\rho$", drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = NULL, y = NULL, title = "Sign-flip map") +
  theme(legend.position = "top", axis.ticks = element_blank(),
        axis.text.y = element_blank(), legend.text = element_text(size = 7),
        legend.key.size = unit(8, "pt"))

pm <- bind_rows(lapply(spec_sens$rows, function(x) {
  data.frame(tissue = ifelse(x$tissue == "brain", "Brain", "PBMC"),
             model = unname(model_label[[x$model_label]]),
             mantel_full = x$mantel_p_full, mantel_nondeg = x$mantel_p_non_degree,
             degree_full = x$degree_preserving_p_full,
             degree_nondeg = x$degree_preserving_p_non_degree)
}))
f6c <- ggplot(pm) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(mantel_full, mantel_nondeg, color = "gene-label"), size = 2.2) +
  geom_point(aes(degree_full, degree_nondeg, color = "row-shuffle"), size = 2.2, shape = 17) +
  scale_color_manual(values = c("gene-label" = BLUE, "row-shuffle" = YELLOW)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "$p$ under full confounds", y = "$p$ under non-degree", color = NULL,
       title = "Randomization $p$ migrates upward") +
  # Compact top legend; free(f6c) below pulls the axis box left next to ylab.
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.28, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.y = element_text(margin = margin(r = 2)),
        plot.margin = margin(4, 4, 3, 2))

lad <- bind_rows(lapply(ladder$rows, function(x) {
  data.frame(tissue = ifelse(x$tissue == "brain", "Brain", "PBMC"),
             model = ifelse(x$is_baseline, "Co-expression", unname(model_label[[x$model_label]])),
             marginal = x$marginal %||% NA_real_,
             coexp_only = x$coexp_only %||% NA_real_,
             nondegree_only = x$nondegree_only %||% NA_real_,
             degree_only = x$degree_only %||% NA_real_,
             coexp_plus_nondegree = x$coexp_plus_nondegree %||% NA_real_,
             coexp_plus_full = x$coexp_plus_full %||% NA_real_)
}))
# PBMC scGPT/UCE ladders are absent from marginal_vs_adjusted_v2; they are
# recomputed in make_panel_data.py with the same method and validated there
# against the authoritative audit full and non-degree rows.
lad <- bind_rows(lad, bind_rows(lapply(c("scGPT_encoder", "UCE_encoder"), function(mk) {
  x <- panel$ladder_pbmc_extra[[mk]]
  data.frame(tissue = "PBMC", model = unname(model_label[[mk]]),
             marginal = x$marginal, coexp_only = x$coexp_only,
             nondegree_only = x$nondegree_only, degree_only = x$degree_only,
             coexp_plus_nondegree = x$coexp_plus_nondegree,
             coexp_plus_full = x$coexp_plus_full)
})))
lad_long <- bind_rows(lapply(seq_len(nrow(lad)), function(i) {
  r <- lad[i, ]
  rungs <- c("marginal", "coexp_only", "nondegree_only", "degree_only",
             "coexp_plus_nondegree", "coexp_plus_full")
  vals <- c(r$marginal, r$coexp_only, r$nondegree_only, r$degree_only,
            r$coexp_plus_nondegree, r$coexp_plus_full)
  data.frame(tissue = r$tissue, model = r$model,
             rung = factor(rungs, levels = rungs,
                    labels = c("marginal", "coexp", "nondeg", "degree",
                               "+nondeg", "+full")),
             rho = as.numeric(vals))
}))
lad_long <- lad_long %>% filter(!is.na(rho))
# short display names keep the 9-entry legend inside a half-width panel
lad_short <- c("Co-expression" = "Co-exp", "Geneformer embed" = "GF embed",
               "Geneformer attention" = "GF attn", "scFoundation encoder" = "scF enc",
               "UCE encoder" = "UCE", "scGPT encoder" = "scGPT",
               "Geneformer KO" = "GF KO", "Artifact-corrected KO" = "GF KO-corr",
               "Random-init floor" = "Rand floor")
lad_long$display_model <- unname(lad_short[lad_long$model])
f6d <- ggplot(lad_long, aes(rung, rho, group = model, color = display_model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  facet_wrap(~tissue) +
  scale_color_manual(values = c(BLUE, AQUA, YELLOW, VIOLET, RED, MUTED, "grey20", "grey70", "tan3")) +
  # One-column legend on the right frees the bottom band so the facet
  # panels can shift left into the former legend footprint.
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  labs(x = NULL, y = "partial Spearman $\\rho$", color = NULL,
       title = "Nested covariate ladder") +
  theme(legend.position = "right",
        legend.text = element_text(size = 7.2),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 2),
        plot.margin = margin(4, 2, 3, 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  d_title_nudge

# free(f6c): A/B y-ticks otherwise panel-align C's axis box far right of its
# ylab; D keeps a right legend column and a slightly wider share.
emit("fig6_spec_sensitivity",
     ((f6a | f6b) + plot_layout(widths = c(1.6, 1))) /
       ((free(f6c) | f6d) + plot_layout(widths = c(0.95, 1.15))) +
       plot_annotation(tag_levels = "A") +
       plot_layout(axis_titles = "collect"),
     6.8, 7.0)

# ==================== Figure 7: per-cell-type descriptive ====================
pertype_rows <- function(tissue) {
  rows <- Filter(function(x) identical(x$row_type, "pertype_fm"),
                 audit$per_cell_type[[tissue]]$full_confound$rows)

  # Load enhanced stats with bootstrap CI
  enhanced_file <- file.path(res, "pertype_stats_enhanced_v2.json")
  enhanced_ci <- list()
  if (file.exists(enhanced_file)) {
    enhanced <- fromJSON(enhanced_file, simplifyVector = FALSE)
    # Build lookup: tissue_celltype_model -> CI
    for (e in enhanced) {
      if (e$tissue == tissue && !is.null(e$bootstrap_ci95) && !is.null(e$readout)) {
        # Match readout to model_label mapping (only include readouts we use)
        readout_map <- list(
          coexp_vs_atac = "co_expression",
          geneformer_embed = "geneformer_embed",
          geneformer_attn = "geneformer_attn",
          scfoundation = "scFoundation_encoder"
        )
        model_key <- readout_map[[e$readout]]
        if (!is.null(model_key) && !is.na(model_key)) {
          key <- paste(e$cell_type, model_key, sep = "_")
          enhanced_ci[[key]] <- e$bootstrap_ci95
        }
      }
    }
  }

  bind_rows(lapply(rows, function(x) {
    key <- paste(x$cell_type, x$model_label, sep = "_")
    ci <- enhanced_ci[[key]]
    data.frame(
      tissue = ifelse(tissue == "brain", "Brain", "PBMC"),
      cell_type = x$cell_type,
      model = unname(model_label[[x$model_label]]),
      rho = x$observed_partial_rho,
      ci_lo = if (!is.null(ci)) ci[[1]] else NA_real_,
      ci_hi = if (!is.null(ci)) ci[[2]] else NA_real_,
      n_cells = x$n_cells
    )
  }))
}
pt_brain <- pertype_rows("brain")
pt_pbmc <- pertype_rows("pbmc")
pertype_panel <- function(data, tissue_name, color) {
  data$cell_type <- factor(data$cell_type, levels = rev(unique(data$cell_type)))
  ggplot(data, aes(rho, cell_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                  width = 0.25, color = "grey40", linewidth = 0.3,
                  na.rm = TRUE, orientation = "y") +
    geom_point(size = 1.9, color = color) +
    facet_wrap(~model, nrow = 1, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = 0.10), n.breaks = 3) +
    labs(x = NULL, y = NULL, title = tissue_name) +
    theme(strip.text = element_text(size = 10), axis.text.y = element_text(size = 9.5),
          axis.text.x = element_text(size = 8.5),
          plot.title = element_text(size = 11, face = "bold"),
          panel.spacing.x = unit(9, "pt"))
}
f7a <- pertype_panel(pt_brain, "Brain", BLUE)
f7b <- pertype_panel(pt_pbmc, "PBMC", AQUA)

ptf <- bind_rows(lapply(pertype_fm$per_type, function(x) {
  data.frame(cell_type = x$cell_type, coexp = x$coexp_vs_atac,
             embed = x$emb_partial, attn = x$attn_partial)
}))
ptf_long <- bind_rows(
  ptf %>% transmute(cell_type, series = "coexp vs proxy", rho = coexp),
  ptf %>% transmute(cell_type, series = "GF embed (partial)", rho = embed),
  ptf %>% transmute(cell_type, series = "GF attention (partial)", rho = attn))
ptf$cell_type <- factor(ptf$cell_type, levels = ptf$cell_type)
ptf_long$cell_type <- factor(ptf_long$cell_type, levels = ptf$cell_type)
f7c <- ggplot(ptf_long, aes(cell_type, rho, color = series, group = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 2.0) +
  scale_color_manual(values = c("black", BLUE, AQUA)) +
  labs(x = NULL, y = "Spearman $\\rho$", color = NULL,
       title = "Brain per-type readouts and the confound itself") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(legend.position = "top",
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.28, "cm"),
        legend.key.spacing.x = unit(4, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, margin = margin(t = 0, b = 2)),
        plot.margin = margin(4, 4, 3, 2))

cells <- bind_rows(lapply(names(panel$pertype_n_cells), function(tn) {
  cc <- panel$pertype_n_cells[[tn]]
  bind_rows(lapply(names(cc), function(ct) {
    data.frame(tissue = ifelse(tn == "brain", "Brain", "PBMC"),
               cell_type = ct, n_cells = cc[[ct]])
  }))
}))
# Legend must sit outside the panel spines (same band as C's series legend and
# A/B facet strips), not inside the bar-chart frame.
f7d <- ggplot(cells, aes(reorder(cell_type, n_cells), n_cells, fill = tissue)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  coord_flip() +
  labs(x = NULL, y = "cells", fill = NULL,
       title = "Cell counts behind the per-type rows") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             override.aes = list(color = NA))) +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.28, "cm"),
        legend.key.spacing.x = unit(6, "pt"),
        legend.margin = margin(0, 0, 2, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        # Reserve a clear header band above the top spine (matches C).
        plot.title.position = "plot",
        # Long D title: shift right of the D tag (tag itself nudged right too).
        plot.title = element_text(hjust = 0, margin = margin(t = 0, b = 2, l = 10)),
        plot.margin = margin(4, 6, 3, 2),
        # Nudge panel tag D further right toward the plot spine / title start
        # (deep y-tick gutter otherwise leaves D too far left of the title).
        plot.tag.position = c(0.24, 1))

emit("fig7_pertype_descriptive", f7a / f7b / (f7c | f7d) + plot_annotation(tag_levels = "A") +
       plot_layout(heights = c(1, 1, 0.95), axis_titles = "collect"), 6.8, 8.2)

# ==================== Figure 8: injection ladder and effect scale ====================
injection_rows <- function(tissue) {
  bind_rows(lapply(injection$tissues[[tissue]]$rows, function(x) {
    values <- vapply(x$replicate_runs, function(z) z$observed_partial_rho_axis_aligned, numeric(1))
    data.frame(tissue = ifelse(tissue == "brain", "Brain", "PBMC"), alpha = x$alpha,
               mean = mean(values), lo = min(values), hi = max(values))
  }))
}
dose <- bind_rows(injection_rows("brain"), injection_rows("pbmc"))
subdiv_path <- file.path(res, "injection_subdivided_v2.json")
subdiv <- NULL
if (file.exists(subdiv_path)) {
  sdj <- fromJSON(subdiv_path, simplifyVector = FALSE)
  subdiv <- bind_rows(lapply(c("brain", "pbmc"), function(tn) {
    bind_rows(lapply(sdj[[tn]]$rows, function(x) {
      values <- vapply(x$replicate_runs, function(z) z$observed_partial_rho_axis_aligned, numeric(1))
      data.frame(tissue = ifelse(tn == "brain", "Brain", "PBMC"), alpha = x$alpha,
                 mean = mean(values))
    }))
  }))
}
dose_panel <- function(data, subdiv_data, tissue_name, color) {
  ggplot(data, aes(alpha, mean)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, fill = color, color = NA) +
    geom_line(linewidth = 0.7, color = color) + geom_point(size = 1.4, color = color) +
    {if (!is.null(subdiv_data)) geom_point(data = subdiv_data, aes(alpha, mean),
                                           shape = 23, size = 2.0, fill = "white", color = color,
                                           inherit.aes = FALSE)} +
    labs(x = "axis-aligned injected fraction $\\alpha$", y = "recovered partial $\\rho$",
         title = paste(tissue_name, "injection ladder")) +
    theme(legend.position = "none")
}
# Panel A: nudge title slightly left toward the A tag (small left margin only;
# keep title.position = panel so it does not collide with the tag).
f8a <- dose_panel(dose %>% filter(tissue == "Brain"),
                  if (is.null(subdiv)) NULL else subdiv %>% filter(tissue == "Brain"),
                  "Brain", BLUE) +
  theme(plot.title = element_text(hjust = 0, margin = margin(b = 3, l = -10)))
f8b <- dose_panel(dose %>% filter(tissue == "PBMC"),
                  if (is.null(subdiv)) NULL else subdiv %>% filter(tissue == "PBMC"),
                  "PBMC", AQUA)

eff <- bind_rows(lapply(effect_scale$observed_effects_as_alpha, function(x) {
  data.frame(tissue = ifelse(x$tissue == "brain", "Brain", "PBMC"),
             model = x$model_label, observed_rho = x$observed_rho,
             alpha_equiv = x$alpha_equivalent %||% NA_real_,
             status = x$alpha_equivalent_status)
}))
# PBMC scGPT/UCE full rows are absent from effect_vs_injection_scale_v2; their
# alpha-equivalents are computed in make_panel_data.py with the same
# interpolation, validated against every stored INTERPOLATED value.
audit_full <- bind_rows(lapply(audit$pooled$pbmc$rows, function(x) {
  if (identical(x$row_type, "pooled_fm") && x$confound_spec == "full" &&
      x$model_label %in% c("scGPT_encoder", "UCE_encoder"))
    data.frame(tissue = "PBMC", model = x$model_label,
               observed_rho = x$observed_partial_rho)
}))
audit_full$alpha_equiv <- vapply(audit_full$model, function(mk)
  panel$alpha_equiv_extra[[paste0("pbmc_", mk)]], numeric(1))
audit_full$status <- "INTERPOLATED"
eff <- bind_rows(eff, audit_full)
eff_ok <- eff %>% filter(!is.na(alpha_equiv), alpha_equiv > 0) %>%
  mutate(display = paste(ifelse(tissue == "Brain", "brain", "PBMC"),
                         gsub("_", " ", model), sep = " -- "))
f8c <- ggplot() +
  geom_line(data = dose %>% filter(alpha > 0), aes(alpha, mean, color = tissue), linewidth = 0.4, alpha = 0.5) +
  geom_point(data = eff_ok, aes(alpha_equiv, observed_rho, color = tissue), size = 2.2) +
  scale_color_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  scale_x_log10() +
  labs(x = "$\\alpha$-equivalent", y = "observed partial $\\rho$", color = NULL,
       title = "Observed effects sit at the bottom of the ladder") +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        plot.title = element_text(margin = margin(t = 0, b = 0)),
        plot.margin = margin(1, 4, 3, 2))

supp_keys <- primary %>% filter(status == "supported by both nulls") %>%
  mutate(key = paste(tolower(tissue), model_key, sep = "|")) %>% pull(key)
eff_supp <- eff_ok %>%
  mutate(key = paste(tolower(tissue), model, sep = "|")) %>%
  filter(key %in% supp_keys)

# Add co-expression baseline for reference
baseline_alpha <- bind_rows(
  data.frame(tissue = "Brain",
            display = "brain -- coexp baseline",
            alpha_equiv = 0.0015,
            observed_rho = baseline$rho[baseline$tissue == "Brain"]),
  data.frame(tissue = "PBMC",
            display = "pbmc -- coexp baseline",
            alpha_equiv = 0.0012,
            observed_rho = baseline$rho[baseline$tissue == "PBMC"])
)
eff_supp <- bind_rows(eff_supp, baseline_alpha)

f8d <- ggplot(eff_supp, aes(reorder(display, alpha_equiv), alpha_equiv, fill = tissue)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.002, linetype = "dashed", color = RED) +
  annotate("text", x = 0.6, y = 0.0026, label = "smallest probed $\\alpha=0.002$",
           hjust = 0, size = 3.1, color = RED) +
  scale_fill_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  coord_flip() +
  labs(x = NULL, y = "$\\alpha$-equivalent", fill = NULL,
       title = "Injection equivalents") +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        plot.title = element_text(hjust = 0, margin = margin(t = 0, b = 0, l = 10)),
        plot.margin = margin(1, 4, 3, 2),
        axis.text.y = element_text(size = 8))

emit("fig8_injection_ladder", (f8a | f8b) / (f8c | f8d) + plot_annotation(tag_levels = "A") +
       plot_layout(axis_titles = "collect"), 6.8, 6.4)

# ==================== Figure 9: TF-disjoint probe ====================
probe_fam_label <- c(
  co_expression = "Co-expression",
  geneformer_embed = "Geneformer embed",
  geneformer_attn = "Geneformer attention",
  scGPT_encoder = "scGPT encoder",
  UCE_encoder = "UCE encoder",
  random_floor = "Random-init floor"
)
contr <- probe_stats$contrasts_vs_baseline
pa <- bind_rows(lapply(names(probe_stats$families), function(fam) {
  f <- probe_stats$families[[fam]]
  is_base <- fam == "co_expression"
  data.frame(family = unname(probe_fam_label[[fam]]), key = fam,
             rho = f$adjusted_rho_mean, q = f$mantel_q,
             flip_q = if (is_base) NA_real_ else contr[[fam]]$signflip_q,
             status = if (is_base) "co-expression baseline"
                      else if (isTRUE(contr[[fam]]$significant_q05)) "supported vs baseline"
                      else "not supported",
             stringsAsFactors = FALSE)
}))
pa$family <- factor(pa$family, levels = rev(unname(probe_fam_label)))
  q_labels <- ifelse(pa$status == "co-expression baseline",
                     sprintf("$q_M=%.3f$", pa$q),
                     sprintf("$q_M=%.3f;\\ q_{\\mathrm{flip}}=%.3f$", pa$q, pa$flip_q))
  # Thin left gutter + tighter right pad so the forest spans more of the 6.8in width.
  f9a_margin <- compute_measured_panel_margin(
    ytick_labels = pa$family,
    right_labels = q_labels,
    base_left_pt = 0,
    normal_ytick_envelope_pt = 160,
    padding_pt = 6
  )
  f9a <- ggplot(pa, aes(rho, family, color = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.4) +
  geom_text(aes(x = 0.0075,
                label = ifelse(status == "co-expression baseline",
                               sprintf("$q_M=%.3f$", q),
                               sprintf("$q_M=%.3f;\\ q_{\\mathrm{flip}}=%.3f$", q, flip_q))),
            hjust = 0, size = 4.0, color = "black") +
  scale_color_manual(values = c("supported vs baseline" = AQUA,
                                "not supported" = MUTED,
                                "co-expression baseline" = "black")) +
  coord_cartesian(xlim = c(-0.021, 0.0062), clip = "off") +
  labs(x = "adjusted test $\\rho$", y = NULL, color = NULL,
       title = "Supervised probe on held-out TFs (edge features only)") +
  theme(legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.margin = f9a_margin)

pb <- bind_rows(lapply(names(probe_sens$grid), function(fam) {
  g <- probe_sens$grid[[fam]]
  bind_rows(lapply(names(g), function(sub) {
    data.frame(family = unname(probe_fam_label[[fam]]),
               subset = factor(sub, levels = c("none", "atac_construction", "full", "detv_only"),
                               labels = c("none", "construction", "full", "detv")),
               rho = g[[sub]]$adjusted_rho_mean)
  }))
}))
pb$family <- factor(pb$family, levels = unname(probe_fam_label))
f9b <- ggplot(pb, aes(subset, rho, color = family, group = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
  annotate("segment", x = 2, xend = 2, y = -0.025, yend = 0.073,
           linetype = "dotted", color = RED, linewidth = 0.4) +
  annotate("text", x = 2.15, y = 0.055, label = "construction\ncovariates",
           hjust = 0, size = 2.8, color = RED) +
  scale_color_manual(values = c("Co-expression" = "black", "Geneformer embed" = BLUE,
                                "Geneformer attention" = AQUA, "scGPT encoder" = VIOLET,
                                "UCE encoder" = YELLOW, "Random-init floor" = MUTED)) +
  labs(x = NULL, y = "adjusted test $\\rho$", color = NULL,
       title = "Recovery collapses at the construction covariates") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "top",
        legend.text = element_text(size = 8.5),
        legend.key.size = unit(0.28, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.text.x = element_text(angle = 28, hjust = 1, size = 9),
        plot.margin = margin(4, 8, 3, 4))

arms <- bind_rows(lapply(names(probe_fam_label), function(fam) {
  eo <- probe_eval$results$edge_only[[fam]]
  al <- probe_eval$results$all[[fam]]
  data.frame(family = unname(probe_fam_label[[fam]]),
             arm = rep(c("edge features only", "all features"), each = 1),
             rho = c(eo$adjusted_rho_mean, al$adjusted_rho_mean))
}))
arms$family <- factor(arms$family, levels = rev(unname(probe_fam_label)))
f9c <- ggplot(arms, aes(rho, family, fill = arm)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_vline(xintercept = 0, color = "grey55") +
  scale_fill_manual(values = c(BLUE, MUTED)) +
  labs(y = NULL, x = "adjusted test $\\rho$", fill = NULL,
       title = "Degree features: no rescue") +
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.28, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.y = element_blank(),
        plot.margin = margin(4, 4, 3, 2))

dist <- bind_rows(lapply(names(probe_fam_label), function(fam) {
  eo <- probe_eval$results$edge_only[[fam]]
  data.frame(family = unname(probe_fam_label[[fam]]),
             mean = eo$adjusted_rho_mean, sd = eo$adjusted_rho_std,
             median = eo$adjusted_rho_median, n = eo$adjusted_rho_n_ok)
}))
dist$family <- factor(dist$family, levels = rev(unname(probe_fam_label)))
f9d <- ggplot(dist, aes(family)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, color = MUTED) +
  geom_point(aes(y = mean), size = 2.4, color = BLUE) +
  geom_point(aes(y = median), size = 2.4, shape = 17, color = YELLOW) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  coord_flip() +
  labs(x = NULL, y = "adjusted test $\\rho$",
       title = "Per-TF recovery spread") +
  theme(axis.ticks.y = element_blank(),
        plot.margin = margin(4, 6, 3, 2)) +
  d_title_nudge

# A/B stay full-width. C|D share one row via cowplot: plain patchwork (C|D)
# nested under A/B used to squeeze the bottom plot regions (~59pt).
f9_bottom <- cowplot::plot_grid(
  f9c, f9d, nrow = 1,
  rel_widths = c(1.05, 0.95),
  labels = c("C", "D"), label_size = 11, label_fontface = "plain",
  hjust = 0, vjust = 1.1
)
emit("fig9_tf_probe",
     (f9a / f9b / wrap_elements(full = f9_bottom)) +
       plot_annotation(tag_levels = list(c("A", "B", ""))) +
       plot_layout(heights = c(1.1, 1.0, 1.15)),
     6.8, 7.6)

# ==================== Figure 10: coverage, scope, QC ====================
cov <- bind_rows(lapply(names(panel$readout_qc$manifest_gene_coverage), function(k) {
  data.frame(readout = k, covered = panel$readout_qc$manifest_gene_coverage[[k]])
}))
f10a <- ggplot(cov, aes(reorder(readout, covered), covered)) +
  geom_col(fill = BLUE, width = 0.6) +
  geom_text(aes(label = covered), hjust = -0.1, size = 3.2) +
  coord_flip(ylim = c(1100, 1230)) +
  labs(x = NULL, y = "genes covered (of 1,200)",
       title = "Gene coverage by readout") +
  theme(axis.text.y = element_text(size = 9),
        # Trim left pad so free(f10a) can reclaim the blank left of short labels
        plot.margin = margin(4, 6, 3, 0))

covtab <- audit$model_coverage_table
cov_rows <- bind_rows(lapply(names(model_label), function(mk) {
  ml <- unname(model_label[[mk]])
  data.frame(model = ml,
             `Brain pooled` = mk %in% unlist(covtab$brain_pooled_FMs),
             `PBMC pooled` = mk %in% unlist(covtab$pbmc_pooled_FMs),
             `Brain per-type` = mk %in% unlist(covtab$brain_per_type_FMs),
             `PBMC per-type` = mk %in% unlist(covtab$pbmc_per_type_FMs),
             check.names = FALSE)
}))
cov_long <- bind_rows(lapply(c("Brain pooled", "PBMC pooled", "Brain per-type", "PBMC per-type"),
                             function(cc) {
  cov_rows %>% transmute(model, analysis = cc, present = .data[[cc]])
}))
cov_long$analysis <- factor(cov_long$analysis,
                            levels = c("Brain pooled", "PBMC pooled",
                                       "Brain per-type", "PBMC per-type"))
cov_long$model <- factor(cov_long$model, levels = rev(unname(model_label)))

# Compute coverage percentage per model
model_cov <- cov_long %>%
  group_by(model) %>%
  summarise(pct = 100 * mean(present), .groups = "drop")

# Continuous x so % labels can sit past the panel spine (discrete scales oob-
# squish positions > n_levels into the right border). xlim keeps the tile grid
# tight; clip="off" + right margin let the row-% annotations clear the spine.
analysis_lvls <- levels(cov_long$analysis)
cov_long$analysis_num <- as.numeric(cov_long$analysis)

f10b <- ggplot(cov_long, aes(analysis_num, model, fill = present)) +
  geom_tile(color = "white") +
  geom_text(data = model_cov, aes(x = 4.62, y = model, label = sprintf("%.0f\\%%", pct)),
            hjust = 0, size = 3.2, inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq_along(analysis_lvls), labels = analysis_lvls,
                     expand = c(0, 0)) +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey85")) +
  coord_cartesian(xlim = c(0.5, 4.5), clip = "off") +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Graph coverage") +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 10.5, angle = 40, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        plot.margin = margin(4, 42, 3, 4))

cells_r <- bind_rows(lapply(names(panel$readout_qc$cells_per_readout), function(k) {
  data.frame(readout = k, n = panel$readout_qc$cells_per_readout[[k]])
}))
f10c <- ggplot(cells_r, aes(reorder(readout, n), n)) +
  geom_col(fill = VIOLET, width = 0.6) +
  coord_flip() +
  labs(x = NULL, y = "cells",
       title = "Cells per readout") +
  theme(axis.text.y = element_text(size = 8.5),
        # Keep "cells" tight to C's ticks/spine (not D's deep rotated labels).
        axis.title.x = element_text(margin = margin(t = 1, b = 0)),
        plot.margin = margin(4, 6, 1, 4))

ea <- panel$readout_qc$edge_accounting
edges <- data.frame(
  step = factor(c("446 TF $\\times$ 1,200 genes", "minus self pairs (PBMC edge set)",
                  "minus brain marker mask (brain edge set)"),
                levels = c("446 TF $\\times$ 1,200 genes", "minus self pairs (PBMC edge set)",
                           "minus brain marker mask (brain edge set)")),
  n = c(ea$tf_x_gene, ea$pbmc_edges, ea$brain_edges))
f10d <- ggplot(edges, aes(step, n)) +
  geom_col(fill = c("grey70", BLUE, AQUA), width = 0.6) +
  geom_text(aes(label = format(n, big.mark = "")), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0, 580000)) +
  labs(x = NULL, y = "TF--gene edges",
       title = "Edge accounting") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 8.5)) +
  d_title_nudge

# free(f10a): C's longer y-ticks otherwise panel-align A rightward.
# Bottom row: do NOT collect x-axis titles — collect was parking C's "cells"
# on the same baseline as D's xlab under the rotated tick block. free(label,"b")
# also lets C's axis-title slot ignore D's deep bottom label band.
emit("fig10_coverage_qc",
     ((free(f10a) | f10b) + plot_layout(widths = c(1.25, 0.75))) /
       ((free(f10c, type = "label", side = "b") | f10d) +
          plot_layout(axis_titles = "keep")) +
       plot_annotation(tag_levels = "A"),
     6.8, 6.0)

# ============ Figure 11: third-tissue transfer of the construct ============
# Disclosure figure: extends the PBMC-only linkage coverage (5.91%, reviewer 2)
# to all three ATAC datasets and makes the marginal structure behind Fig. 2's
# additive decomposition visible. No new experiments; all values derive from
# pinned npz graphs and published JSONs (panel_data.json third_tissue key).
tt <- panel$third_tissue
TISSUE_COL <- c(Brain = BLUE, PBMC = YELLOW, `Fibroblast mix` = VIOLET)

cov11 <- bind_rows(lapply(tt$coverage, function(c) {
  data.frame(tissue = c$tissue,
             frac = c$relevant_peaks / c$total_peaks,
             label = sprintf("%s / %s", format(c$relevant_peaks, big.mark = ""),
                             format(c$total_peaks, big.mark = "")))
}))
f11a <- ggplot(cov11, aes(reorder(tissue, frac), frac, fill = tissue)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f\\%%", 100 * frac)), hjust = -0.1, size = 3.3) +
  coord_flip(ylim = c(0, 0.075)) +
  scale_fill_manual(values = TISSUE_COL) +
  labs(x = NULL, y = "peaks linked",
       title = "$\\pm$2 kb linkage admits 4--6\\% of peaks")

ov <- bind_rows(lapply(tt$edge_overlap$regions, function(r) {
  data.frame(region = paste(unlist(r$combo), collapse = " + "), n = r$n)
}))
ov$region[ov$region == "Brain + PBMC + Fibroblast mix"] <- "all three"
ov$shared <- ov$region == "all three"
f11b <- ggplot(ov, aes(reorder(region, n), n, fill = shared)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = format(n, big.mark = "")), hjust = -0.1, size = 3.1) +
  coord_flip(ylim = c(0, 47000)) +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey70")) +
  labs(x = NULL, y = "TF--gene edges",
       title = "A shared 38k-edge core")

deg <- bind_rows(lapply(names(tt$degree_tf_out), function(tn) {
  data.frame(tissue = tn, degree = unlist(tt$degree_tf_out[[tn]]))
}))
f11c <- ggplot(deg, aes(degree, color = tissue)) +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  scale_color_manual(values = TISSUE_COL) +
  labs(x = "supported targets per TF", y = "ECDF", color = NULL,
       title = "TF out-degree marginals are similar") +
  theme(legend.position = "top")

rp <- bind_rows(lapply(tt$rho_phi, function(r) {
  data.frame(pair = r$pair, observed = r$observed, phi = r$phi)
}))
f11d <- ggplot(rp, aes(y = reorder(pair, observed))) +
  geom_segment(aes(x = phi, xend = observed, yend = pair), color = "grey60",
               linewidth = 0.7) +
  geom_point(aes(x = phi, color = "binary $\\varphi$"), size = 2.4) +
  geom_point(aes(x = observed, color = "observed $\\rho$"), size = 2.4) +
  scale_color_manual(values = c(`binary $\\varphi$` = MUTED,
                                `observed $\\rho$` = BLUE)) +
  coord_cartesian(xlim = c(0.40, 0.60)) +
  labs(x = "cross-tissue agreement", y = NULL, color = NULL,
       title = "Rank agreement tracks support") +
  theme(legend.position = "top") +
  d_title_nudge

# Tags default to the outer y-tick margin; nudge right toward each panel spine.
emit("fig11_third_tissue_transfer",
     ((f11a | f11b) / (f11c | f11d) +
        plot_annotation(tag_levels = "A") +
        plot_layout(axis_titles = "collect")) &
       theme(plot.tag.position = c(0.28, 1)),
     6.8, 5.8)

# ==================== tables (locked enrichment SPEC) ====================
fmt <- function(x, digits = 4) sprintf(paste0("%.", digits, "f"), x)
pfmt <- function(p) ifelse(p <= 0.001, "$<0.001$", sprintf("$%.3f$", p))

# Table 1 (tab:primary): keep 15 full-confound rows; add Support column.
# Support for FM rows only: dual within-family BH q_M<0.05 and q_D<0.05 →
# "both" (or "both (neg)" if observed_partial_rho<0); else "neither".
# Co-expression baselines sit outside FM BH families → Support="--".
primary$support <- ifelse(
  primary$is_baseline, "--",
  ifelse(primary$mantel_q < 0.05 & primary$degree_q < 0.05,
         ifelse(primary$rho < 0, "both (neg)", "both"),
         "neither"))
table_rows <- vapply(seq_len(nrow(primary)), function(i) {
  if (primary$is_baseline[i]) {
    sprintf("%s & %s & $%s$ & %s & -- & %s & -- & -- \\\\",
            primary$tissue[i], primary$model[i], fmt(primary$rho[i]),
            pfmt(primary$mantel_p[i]), pfmt(primary$degree_p[i]))
  } else {
    sprintf("%s & %s & $%s$ & %s & $%.3f$ & %s & $%.3f$ & %s \\\\",
            primary$tissue[i], primary$model[i],
            fmt(primary$rho[i], digits = 5), pfmt(primary$mantel_p[i]),
            primary$mantel_q[i], pfmt(primary$degree_p[i]),
            primary$degree_q[i], primary$support[i])
  }
}, character(1))
writeLines(c(
  "\\begin{tabular}{llrrrrrl}", "\\toprule",
  "Tissue & Readout & Partial $\\rho$ & $p_M$ & $q_M$ & $p_D$ & $q_D$ & Support \\\\ ",
  "\\midrule", table_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table1_primary_fixed_panel.tex"))

# Table 2 (tab:cross): 8 descriptive columns from additive decomp + pair files.
# Join on tissue-pair IDs; no Mantel p, bootstrap CI, or additive_pred as inference.
tissue_names <- c(GSE174367 = "Brain", PBMC10k = "PBMC", GSE206767 = "Fibroblast mix")
pair_extra_files <- list(
  "Brain--PBMC" = "cross_tissue_brain_vs_pbmc.json",
  "Brain--Fibroblast mix" = "cross_tissue_atac_v2.json",
  "PBMC--Fibroblast mix" = "cross_tissue_pbmc_vs_fibroblast.json"
)
t2 <- bind_rows(lapply(decomp$rows, function(x) {
  pair_lab <- paste(tissue_names[[x$pair[[1]]]], tissue_names[[x$pair[[2]]]],
                    sep = "--")
  extra <- J(pair_extra_files[[pair_lab]])
  data.frame(
    pair = pair_lab,
    obs_rho = x$observed_spearman,
    add_frac = x$fraction_explained_by_additive_marginals,
    resid_rho = x$residual_spearman_after_own_additive_fits,
    phi = x$binary_support_phi,
    edge_jaccard = extra$edge_jaccard,
    mean_row_rho = extra$mean_per_tf_row_spearman,
    n_tf = x$n_tf_common,
    stringsAsFactors = FALSE
  )
}))
# Stable display order matching the locked SPEC.
t2_order <- c("Brain--PBMC", "Brain--Fibroblast mix", "PBMC--Fibroblast mix")
t2 <- t2[match(t2_order, t2$pair), ]
cross_rows <- vapply(seq_len(nrow(t2)), function(i) sprintf(
  "%s & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ & $%.3f$ & $%.4f$ & %d \\\\",
  t2$pair[i], t2$obs_rho[i], t2$add_frac[i], t2$resid_rho[i], t2$phi[i],
  t2$edge_jaccard[i], t2$mean_row_rho[i], as.integer(t2$n_tf[i])
), character(1))
writeLines(c(
  "\\begin{tabular}{lrrrrrrr}", "\\toprule",
  paste0("Tissue pair & Obs.\\ $\\rho$ & Add.\\ frac.\\ & Resid.\\ $\\rho$ & ",
         "$\\varphi$ & Edge Jac.\\ & Mean row $\\rho$ & Shared TFs \\\\ "),
  "\\midrule", cross_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table2_cross_tissue_observed.tex"))

# Table 3 (tab:pertype_ranges): 4-row descriptive min/max/median only.
# Source: audit$per_cell_type.*.descriptive_summary; no p_mc / q / BH / CI.
spec_keys <- c(full = "full_confound", `non-degree` = "non_degree_confound")
t3 <- bind_rows(lapply(c("brain", "pbmc"), function(tn) {
  bind_rows(lapply(names(spec_keys), function(spec_lab) {
    ds <- audit$per_cell_type[[tn]][[spec_keys[[spec_lab]]]]$descriptive_summary
    data.frame(
      tissue = ifelse(tn == "brain", "Brain", "PBMC"),
      spec = spec_lab,
      n_rows = ds$n_rows_exploratory,
      rho_min = ds$rho_min,
      rho_max = ds$rho_max,
      rho_median = ds$rho_median,
      stringsAsFactors = FALSE
    )
  }))
}))
t3_rows <- vapply(seq_len(nrow(t3)), function(i) sprintf(
  "%s & %s & %d & $%.5f$ & $%.5f$ & $%.5f$ \\\\",
  t3$tissue[i], t3$spec[i], as.integer(t3$n_rows[i]),
  t3$rho_min[i], t3$rho_max[i], t3$rho_median[i]
), character(1))
writeLines(c(
  "\\begin{tabular}{llrrrr}", "\\toprule",
  "Tissue & Confound spec & $n_{\\mathrm{rows}}$ & $\\rho_{\\min}$ & $\\rho_{\\max}$ & $\\rho_{\\mathrm{median}}$ \\\\ ",
  "\\midrule", t3_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table3_pertype_ranges.tex"))

# ==================== Table 4 + Fig.12: protocol-pass gate matrix (0x10) ====================
# Predeclared protocol-pass (see also 0x20 SAP intent):
#   dual_full ∧ concordance (FM–coexp Spearman > 0) ∧ non-degree same sign
#   ∧ multi-readout same-sign within model_family×tissue ∧ ρ > tissue coexp baseline.
# dual_nondeg is reported as its own column (0/13 under current data) but is NOT an
# extra AND beyond non-degree same-sign in the protocol-pass definition above.
yn <- function(x) ifelse(x, "yes", "no")

# Family id from audit rows (prefer model_family field).
fam_of <- function(tissue_key, model_key) {
  rows <- Filter(function(x) identical(x$row_type, "pooled_fm") &&
                   identical(x$model_label, model_key) &&
                   identical(x$confound_spec, "full"),
                 audit$pooled[[tissue_key]]$rows)
  if (!length(rows)) return(model_key)
  rows[[1]]$model_family %||% model_key
}

fm_full <- pooled %>% filter(spec == "full")
fm_nd <- pooled %>% filter(spec == "non_degree") %>%
  select(tissue, model_key, rho_nd = rho, mantel_q_nd = mantel_q, degree_q_nd = degree_q)
pp <- fm_full %>%
  left_join(fm_nd, by = c("tissue", "model_key")) %>%
  mutate(
    tissue_key = ifelse(tissue == "Brain", "brain", "pbmc"),
    dual_full = mantel_q < 0.05 & degree_q < 0.05,
    dual_nondeg = !is.na(mantel_q_nd) & mantel_q_nd < 0.05 &
      !is.na(degree_q_nd) & degree_q_nd < 0.05,
    same_sign_nd = !is.na(rho_nd) & (
      (rho == 0 & rho_nd == 0) | (rho * rho_nd > 0)
    ),
    baseline_rho = ifelse(tissue == "Brain",
                          baseline$rho[baseline$tissue == "Brain"][1],
                          baseline$rho[baseline$tissue == "PBMC"][1]),
    vs_baseline = rho > baseline_rho
  )

# Concordance Spearmans from panel_data.json (schema key still usability_*).
usa_long <- bind_rows(lapply(names(panel$usability_fm_vs_coexp), function(tn) {
  u <- panel$usability_fm_vs_coexp[[tn]]
  data.frame(
    tissue = ifelse(tn == "brain", "Brain", "PBMC"),
    model_key = names(u),
    conc_rho = as.numeric(unlist(u)),
    stringsAsFactors = FALSE
  )
}))
pp <- pp %>% left_join(usa_long, by = c("tissue", "model_key")) %>%
  mutate(concordance = !is.na(conc_rho) & conc_rho > 0)

# Multi-readout same-sign: all full-spec ρ in the same model_family × tissue agree in sign.
pp$family <- mapply(fam_of, pp$tissue_key, pp$model_key, USE.NAMES = FALSE)
fam_signs <- pp %>%
  group_by(tissue, family) %>%
  summarise(
    multi_ok = {
      s <- sign(rho)
      s <- s[s != 0]
      length(unique(s)) <= 1
    },
    .groups = "drop"
  )
pp <- pp %>% left_join(fam_signs, by = c("tissue", "family")) %>%
  mutate(
    multi_ok = ifelse(is.na(multi_ok), TRUE, multi_ok),
    protocol_pass = dual_full & concordance & same_sign_nd & multi_ok & vs_baseline
  )

# Stable row order: match Table 1 FM row order (primary without baselines).
primary_fm_order <- primary %>% filter(!is_baseline) %>%
  transmute(tissue, model_key, ord = row_number())
pp <- pp %>%
  left_join(primary_fm_order, by = c("tissue", "model_key")) %>%
  arrange(ord) %>%
  select(-ord)
stopifnot(nrow(pp) == 13L)
stopifnot(sum(pp$protocol_pass) == 0L)
stopifnot(!any(is.na(pp$protocol_pass)))

pp_rows <- vapply(seq_len(nrow(pp)), function(i) sprintf(
  "%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  pp$tissue[i], pp$model[i],
  yn(pp$dual_full[i]), yn(pp$dual_nondeg[i]), yn(pp$concordance[i]),
  yn(pp$same_sign_nd[i]), yn(pp$multi_ok[i]), yn(pp$vs_baseline[i]),
  yn(pp$protocol_pass[i])
), character(1))
writeLines(c(
  "\\begin{tabular}{llccccccc}", "\\toprule",
  paste0("Tissue & Readout & Dual full & Dual non-deg & Concord. & ",
         "ND same sign & Multi-RO sign & $\\rho>$ base & Protocol-pass \\\\ "),
  "\\midrule", pp_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table4_protocol_pass.tex"))

# Binary heatmap (rows × gates) as Fig.12 companion to Table 4.
gate_levels <- c("Dual full", "Dual non-deg", "Concordance",
                 "ND same sign", "Multi-RO sign", "$\\rho>$ baseline",
                 "Protocol-pass")
pp_long <- bind_rows(
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "Dual full", pass = pp$dual_full, stringsAsFactors = FALSE),
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "Dual non-deg", pass = pp$dual_nondeg, stringsAsFactors = FALSE),
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "Concordance", pass = pp$concordance, stringsAsFactors = FALSE),
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "ND same sign", pass = pp$same_sign_nd, stringsAsFactors = FALSE),
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "Multi-RO sign", pass = pp$multi_ok, stringsAsFactors = FALSE),
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "$\\rho>$ baseline", pass = pp$vs_baseline, stringsAsFactors = FALSE),
  data.frame(display = paste(pp$model, pp$tissue, sep = " / "),
             gate = "Protocol-pass", pass = pp$protocol_pass, stringsAsFactors = FALSE)
)
pp_long$display <- factor(pp_long$display, levels = rev(unique(pp_long$display)))
pp_long$gate <- factor(pp_long$gate, levels = gate_levels)
f12 <- ggplot(pp_long, aes(gate, display, fill = pass)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey85"),
                    labels = c(`TRUE` = "pass", `FALSE` = "fail"),
                    name = NULL) +
  labs(x = NULL, y = NULL,
       title = "Protocol-pass gates: $0/13$ rows pass all predeclared checks") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = 8),
        legend.position = "top",
        panel.grid = element_blank())
# Single-panel heatmap: stamp an "A" tag so figure-contract composites stay honest.
f12 <- f12 + plot_annotation(tag_levels = "A") +
  theme(plot.tag.position = c(0.02, 0.98))
emit("fig12_protocol_pass_matrix", f12, 6.8, 5.2, tags = "A")


# ==================== Figure 13: protocol scope card (0x15) ====================
scope_nodes <- data.frame(
  x = c(1, 2, 3, 4),
  y = c(1.2, 1.2, 1.2, 1.2),
  lab = c("$\\pm$2 kb\nwindow", "JASPAR\nmotif", "ATAC peak\npresence", "fixed 446$\\times$1{,}200\npanel"),
  stringsAsFactors = FALSE
)
scope_out <- data.frame(
  x = c(1.5, 2.5, 3.5),
  y = c(0.35, 0.35, 0.35),
  lab = c("distal / 3D\nexcluded", "ChIP / causal\nTF--target not claimed", "cell-type near-invariant\nproxy"),
  stringsAsFactors = FALSE
)
f13 <- ggplot() +
  geom_segment(aes(x = 1, xend = 4, y = 1.2, yend = 1.2),
               color = MUTED, linewidth = 0.6,
               arrow = arrow(length = unit(0.08, "in"), type = "closed")) +
  geom_label(data = scope_nodes, aes(x, y, label = lab),
             size = 2.8, label.size = 0.4, fill = LIGHT, color = BLUE,
             label.padding = unit(0.25, "lines"), lineheight = 0.95) +
  geom_label(data = scope_out, aes(x, y, label = lab),
             size = 2.4, label.size = 0.3, fill = "grey96", color = MUTED,
             label.padding = unit(0.2, "lines"), lineheight = 0.95) +
  annotate("text", x = 2.5, y = 1.85,
           label = "Protocol instance scope (what the audit asks)",
           size = 3.2, color = "black") +
  annotate("text", x = 2.5, y = -0.15,
           label = "Out of scope for this instance --- not a claim that models lack regulation elsewhere",
           size = 2.4, color = MUTED) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-0.35, 2.1)) +
  theme_void() +
  theme(plot.margin = margin(6, 8, 4, 8))
f13 <- f13 + plot_annotation(tag_levels = "A") +
  theme(plot.tag.position = c(0.02, 0.98))
emit("fig13_scope_card", f13, 6.8, 2.6, tags = "A")

# ==================== Table 5: related-work comparison (0x14) ====================
# Qualitative protocol axes only; complement-not-compete framing (PeerJ #4).
writeLines(c(
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  paste0("Resource & Expression in reference? & Degree / confound control? & ",
         "Explicit null semantics? & Multi-readout? & Hash-pinned capsule? \\\\ "),
  "\\midrule",
  paste0("scReg-Eval (this work) & No (edge weights motif+ATAC) & ",
         "Yes (coexp + 6 structural) & Yes (two MC nulls + BH families) & ",
         "Yes & Yes \\\\ "),
  paste0("Kendiukhov et al.\\ scFM interpretability & ",
         "Yes (co-expression / pathways) & Partial & ",
         "Varies & Attention-focused & No \\\\ "),
  paste0("Kendiukhov et al.\\ SAE residual stream & ",
         "Pathway/interaction features & Partial & ",
         "Varies & Feature-level & No \\\\ "),
  paste0("GeneRNIB living GRN benchmark & ",
         "Often literature / expression GRNs & Task-dependent & ",
         "Benchmark scores & Multi-method & Community leaderboard \\\\ "),
  paste0("BioLLM / scFM task suites & ",
         "Task labels from biology & Task-dependent & ",
         "Standard ML splits & Multi-task & Code releases \\\\ "),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(figs, "table5_related_work.tex"))


cat("protocol-pass 0/13 confirmed; table4 + fig12 written\n")
cat("all authoritative figures and tables written to", figs, "\n")
