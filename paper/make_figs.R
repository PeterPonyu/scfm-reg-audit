#!/usr/bin/env Rscript
# Authoritative figures for the fixed-panel audit. Run from paper/:
#   Rscript make_figs.R
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(jsonlite)
  library(tikzDevice)
})
options(tikzDefaultEngine = "pdftex")
FONTPKGS <- c("\\usepackage{amsmath}", "\\usepackage{amssymb}",
              "\\usepackage{helvet}", "\\renewcommand{\\familydefault}{\\sfdefault}")
# Text metrics must be measured at the manuscript's own base size. manuscript.tex
# is 11pt; tikzDevice defaults to 10pt, which under-reserves every string by ~10%
# and pushes left-anchored axis labels outside the bounding box, where the
# device's own clip path shears off their first characters.
options(tikzDocumentDeclaration = "\\documentclass[11pt]{article}",
        tikzMetricPackages = c(getOption("tikzMetricPackages"), FONTPKGS),
        tikzLatexPackages = c(getOption("tikzLatexPackages"), FONTPKGS),
        tikzMetricsDictionary = ".tikz_metrics_pdftex")

base <- Sys.getenv("SCFM_BASE", "..")
res <- file.path(base, "results", "v2")
figs <- file.path(base, "paper", "figs")
dir.create(figs, showWarnings = FALSE, recursive = TRUE)
audit <- fromJSON(file.path(res, "fixed_panel_audit_v2.json"), simplifyVector = FALSE)
injection <- fromJSON(file.path(res, "fixed_panel_signal_injection_v2.json"), simplifyVector = FALSE)
probe_stats <- fromJSON(file.path(res, "tf_probe_pair_stats_v2.json"), simplifyVector = FALSE)
probe_sens <- fromJSON(file.path(res, "tf_probe_pair_sensitivity_v2.json"), simplifyVector = FALSE)

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

emit <- function(name, plot, width, height) {
  tikz(file.path(figs, paste0(name, ".tex")), width = width, height = height,
       standAlone = FALSE, sanitize = FALSE)
  print(plot)
  dev.off()
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

# Figure 1: construct schematic and observed cross-tissue reproducibility.
# The flow is laid out as an S (two boxes per row) so the labels stay large
# enough to read after \fitfig scaling; a single row of four boxes at this
# font size is wider than the half-page panel.
flow <- data.frame(x = c(1, 3, 3, 1), y = c(2, 2, 1, 1),
                   label = c("accessible peak", "motif match", "TF-target weight", "fixed panel"))
f1a <- ggplot(flow, aes(x, y)) +
  geom_label(aes(label = label), size = 3.6, linewidth = 0.3, fill = "white") +
  geom_segment(aes(x = 1.62, xend = 2.38, y = 2, yend = 2),
               arrow = arrow(length = unit(0.08, "in")), color = MUTED) +
  geom_segment(aes(x = 3, xend = 3, y = 1.74, yend = 1.26),
               arrow = arrow(length = unit(0.08, "in")), color = MUTED) +
  geom_segment(aes(x = 2.38, xend = 1.62, y = 1, yend = 1),
               arrow = arrow(length = unit(0.08, "in")), color = MUTED) +
  annotate("text", x = 2, y = 0.52,
           label = "regulatory-potential proxy; not causal ground truth", size = 3.2, color = RED) +
  coord_cartesian(xlim = c(0.3, 3.7), ylim = c(0.35, 2.35), clip = "off") +
  labs(title = "Accessibility and motif evidence define the audited proxy") +
  theme_void(base_size = 11) + theme(plot.title = element_text(size = 11, hjust = 0))

cross <- bind_rows(lapply(audit$cross_tissue_construct_reproducibility$rows, function(x) {
  names <- c(GSE174367 = "Brain", PBMC10k = "PBMC", GSE206767 = "Fibroblast mix")
  data.frame(pair = paste(names[[x$pair[[1]]]], names[[x$pair[[2]]]], sep = "--"),
             rho = x$observed_spearman)
}))
f1b <- ggplot(cross, aes(reorder(pair, rho), rho)) +
  geom_col(fill = BLUE, width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.5, size = 3.3) +
  coord_cartesian(ylim = c(0, 0.58)) +
  labs(x = NULL, y = "observed Spearman $\\rho$", title = "Proxy structure is reproducible across tissues") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
emit("fig1_truth_construct", (f1a | f1b) + plot_annotation(tag_levels = "A"), 6.8, 2.8)

# Figure 2: primary full-confound fixed-panel effects and both randomization decisions.
primary <- pooled %>% filter(spec == "full")
# Interleave the co-expression baseline after each tissue block; it is a reference,
# not an FM row, and carries uncorrected p-values only (no BH family of its own).
primary <- bind_rows(primary, baseline)
primary$is_baseline <- primary$model_key == "co_expression"
primary$status <- ifelse(primary$is_baseline, "co-expression baseline", primary$status)
tissue_order <- unique(primary$tissue)
primary <- primary[order(match(primary$tissue, tissue_order), primary$is_baseline), ]
primary$display <- paste(primary$model, primary$tissue, sep = " -- ")
primary$display <- factor(primary$display, levels = rev(primary$display))
# Forest-plot layout: the q/p labels live in a fixed right-hand margin column
# (clip = "off"), so no label ever collides with the zero line, the y-axis row
# labels, or the panel border.
f2 <- ggplot(primary, aes(rho, display, color = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(shape = is_baseline), size = 2.5) +
  geom_text(aes(x = 0.0165,
                label = ifelse(is_baseline,
                               sprintf("$p_M=%.3f;\\ p_D=%.3f$", mantel_p, degree_p),
                               sprintf("$q_M=%.3f;\\ q_D=%.3f$", mantel_q, degree_q))),
            hjust = 0, size = 3.3, color = "black") +
  scale_color_manual(values = c("supported by both nulls" = AQUA,
                                "supported by one null" = YELLOW,
                                "not supported" = MUTED,
                                "co-expression baseline" = "black")) +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16), guide = "none") +
  coord_cartesian(xlim = c(-0.006, 0.0145), clip = "off") +
  labs(x = "partial Spearman $\\rho$", y = NULL, color = NULL,
       title = "Primary fixed-panel audit: small, readout-specific alignments") +
  theme(legend.position = "top", axis.ticks.y = element_blank(),
        plot.margin = margin(4, 118, 3, 4))
emit("fig2_decisive_result", f2, 6.8, 5.2)

# Figure 3: full versus non-degree confound specifications.
spec <- pooled %>% select(tissue, model, model_key, spec, rho, status)
spec$spec_label <- ifelse(spec$spec == "full", "Full confounds", "Non-degree sensitivity")
spec$display <- paste(spec$model, spec$tissue, sep = " -- ")
spec$display <- factor(spec$display, levels = rev(unique(spec$display)))
f3 <- ggplot(spec, aes(rho, display, group = display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(color = "grey65", linewidth = 0.45) +
  geom_point(aes(fill = spec_label), shape = 21, size = 2.2, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c("Full confounds" = BLUE, "Non-degree sensitivity" = YELLOW)) +
  labs(x = "partial Spearman $\\rho$", y = NULL, fill = NULL,
       title = "Effect direction and support depend on the confound specification") +
  theme(legend.position = "top", axis.ticks.y = element_blank())
emit("fig3_spec_sensitivity", f3, 6.8, 4.0)

# Figure 4: axis-aligned pipeline sensitivity diagnostic, explicitly non-inferential.
injection_rows <- function(tissue) {
  bind_rows(lapply(injection$tissues[[tissue]]$rows, function(x) {
    values <- vapply(x$replicate_runs, function(z) z$observed_partial_rho_axis_aligned, numeric(1))
    data.frame(tissue = ifelse(tissue == "brain", "Brain", "PBMC"), alpha = x$alpha,
               mean = mean(values), lo = min(values), hi = max(values))
  }))
}
dose <- bind_rows(injection_rows("brain"), injection_rows("pbmc"))
# Subdivided calibration points (alpha in {0.002, 0.005, 0.01}) from
# injection_subdivided_v2.json, if present -- these place the observed effects
# on the ladder at their own magnitude.
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
f4 <- ggplot(dose, aes(alpha, mean, color = tissue, fill = tissue)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  {if (!is.null(subdiv)) geom_point(data = subdiv, aes(alpha, mean, color = tissue),
                                    shape = 23, size = 2.0, fill = "white", inherit.aes = FALSE)} +
  scale_color_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  scale_fill_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  labs(x = "axis-aligned injected fraction $\\alpha$", y = "recovered partial $\\rho$",
       color = NULL, fill = NULL,
       title = "Pipeline sensitivity diagnostic (not a power or exclusion analysis)") +
  theme(legend.position = "top")
emit("fig4_pipeline_sensitivity", f4, 6.5, 3.2)

# Figure 5: descriptive per-cell-type full-confound effects.
pertype_rows <- function(tissue) {
  rows <- Filter(function(x) identical(x$row_type, "pertype_fm"),
                 audit$per_cell_type[[tissue]]$full_confound$rows)
  bind_rows(lapply(rows, function(x) data.frame(
    tissue = ifelse(tissue == "brain", "Brain", "PBMC"), cell_type = x$cell_type,
    model = unname(model_label[[x$model_label]]), rho = x$observed_partial_rho,
    n_cells = x$n_cells
  )))
}
pt_brain <- pertype_rows("brain")
pt_pbmc <- pertype_rows("pbmc")
pertype_panel <- function(data, tissue_name, color) {
  data$cell_type <- factor(data$cell_type, levels = rev(unique(data$cell_type)))
  ggplot(data, aes(rho, cell_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(size = 1.9, color = color) +
    facet_wrap(~model, nrow = 1, scales = "free_x") +
    # expand: keep points off the panel border, which clips them at the axis
    # limits. n.breaks caps the tick count so adjacent free_x panels do not
    # collide into runs like "0.0000.000" at their shared edge.
    scale_x_continuous(expand = expansion(mult = 0.10), n.breaks = 3) +
    labs(x = NULL, y = NULL, title = tissue_name) +
    theme(strip.text = element_text(size = 10), axis.text.y = element_text(size = 9.5),
          axis.text.x = element_text(size = 8.5),
          plot.title = element_text(size = 11, face = "bold"),
          panel.spacing.x = unit(9, "pt"))
}
f5 <- pertype_panel(pt_brain, "Brain", BLUE) /
      pertype_panel(pt_pbmc, "PBMC", AQUA) +
      plot_annotation(
        title = "Per-cell-type estimates are descriptive robustness checks",
        theme = theme(plot.title = element_text(size = 11.5, hjust = 0))) &
      labs(x = "full-confound partial $\\rho$ (descriptive)")
emit("fig5_pertype_descriptive", f5, 6.8, 5.2)

# Table 1: primary full-confound fixed-panel results (baseline rows carry raw p,
# not BH q — they sit outside the FM null families).
fmt <- function(x, digits = 4) sprintf(paste0("%.", digits, "f"), x)
pfmt <- function(p) ifelse(p <= 0.001, "$<0.001$", sprintf("$%.3f$", p))
table_rows <- vapply(seq_len(nrow(primary)), function(i) {
  if (primary$is_baseline[i]) {
    sprintf("%s & %s & $%s$ & %s & -- & %s & -- \\\\",
            primary$tissue[i], primary$model[i], fmt(primary$rho[i]),
            pfmt(primary$mantel_p[i]), pfmt(primary$degree_p[i]))
  } else {
    sprintf("%s & %s & $%s$ & %s & $%.3f$ & %s & $%.3f$ \\\\",
            primary$tissue[i], primary$model[i], fmt(primary$rho[i]), pfmt(primary$mantel_p[i]),
            primary$mantel_q[i], pfmt(primary$degree_p[i]), primary$degree_q[i])
  }
}, character(1))
writeLines(c(
  "\\begin{tabular}{llrrrrr}", "\\toprule",
  "Tissue & Readout & Partial $\\rho$ & $p_M$ & $q_M$ & $p_D$ & $q_D$ \\\\ ",
  "\\midrule", table_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table1_primary_fixed_panel.tex"))

# Table 2: observed-only cross-tissue construct reproducibility.
cross_rows <- vapply(seq_len(nrow(cross)), function(i) sprintf(
  "%s & $%.4f$ & 446 \\\\", cross$pair[i], cross$rho[i]), character(1))
writeLines(c(
  "\\begin{tabular}{lrr}", "\\toprule",
  "Tissue pair & Observed Spearman $\\rho$ & Shared TFs \\\\ ",
  "\\midrule", cross_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table2_cross_tissue_observed.tex"))

# Figure 6: TF-disjoint supervised probe. Panel A: confound-adjusted test rho per
# family with BOTH q values that gate the paper's conclusions -- the gene-label
# permutation q_M (is the family's own recovery nonzero) and the paired sign-flip
# q_flip against the co-expression baseline (the headline contrast). Colour
# follows the headline decision (supported against the baseline), so the panel
# cannot be read as rewarding the co-expression baseline's own nonzero q_M.
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
f6a <- ggplot(pa, aes(rho, family, color = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.4) +
  geom_text(aes(x = 0.0068,
                label = ifelse(status == "co-expression baseline",
                               sprintf("$q_M=%.3f$", q),
                               sprintf("$q_M=%.3f;\\ q_{\\mathrm{flip}}=%.3f$", q, flip_q))),
            hjust = 0, size = 3.0, color = "black") +
  scale_color_manual(values = c("supported vs baseline" = AQUA,
                                "not supported" = MUTED,
                                "co-expression baseline" = "black")) +
  coord_cartesian(xlim = c(-0.021, 0.0055), clip = "off") +
  labs(x = "adjusted test Spearman $\\rho$", y = NULL, color = NULL,
       title = "A. Supervised probe on held-out TFs (edge features only)") +
  theme(legend.position = "top", plot.margin = margin(4, 132, 3, 4))
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
f6b <- ggplot(pb, aes(subset, rho, color = family, group = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
  scale_color_manual(values = c("Co-expression" = "black", "Geneformer embed" = BLUE,
                                "Geneformer attention" = AQUA, "scGPT encoder" = VIOLET,
                                "UCE encoder" = YELLOW, "Random-init floor" = MUTED)) +
  labs(x = NULL, y = "adjusted test Spearman $\\rho$", color = NULL,
       title = "B. Effect of confound subset on probe recovery") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 28, hjust = 1, size = 9),
        legend.text = element_text(size = 9))
# stacked, not side-by-side: at 11pt base the two panel titles overlap when the
# panels share a 6.8in row
emit("fig6_tf_probe", f6a / f6b, 6.8, 5.6)

cat("all authoritative figures and tables written to", figs, "\n")
