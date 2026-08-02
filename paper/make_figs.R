#!/usr/bin/env Rscript
# Authoritative figures for the fixed-panel audit. Run from paper/:
#   python3 make_panel_data.py && Rscript make_figs.R
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
# A: flow schematic laid out as an S (two boxes per row) so labels stay large.
flow <- data.frame(x = c(1, 3, 3, 1), y = c(2, 2, 1, 1),
                   label = c("accessible peak", "motif match", "TF-target weight", "fixed panel"))
f1a <- ggplot(flow, aes(x, y)) +
  geom_label(aes(label = label), size = 3.6, linewidth = 0.3, fill = "white") +
  annotate("segment", x = 1.62, xend = 2.38, y = 2, yend = 2,
           arrow = arrow(length = unit(0.08, "in")), color = MUTED) +
  annotate("segment", x = 3, xend = 3, y = 1.74, yend = 1.26,
           arrow = arrow(length = unit(0.08, "in")), color = MUTED) +
  annotate("segment", x = 2.38, xend = 1.62, y = 1, yend = 1,
           arrow = arrow(length = unit(0.08, "in")), color = MUTED) +
  annotate("text", x = 2, y = 0.52,
           label = "regulatory-potential proxy; not causal ground truth", size = 3.2, color = RED) +
  coord_cartesian(xlim = c(0.3, 3.7), ylim = c(0.35, 2.35), clip = "off") +
  labs(title = "Accessibility and motif evidence define the audited proxy") +
  theme_void(base_size = 11) + theme(plot.title = element_text(size = 11, hjust = 0))

f1b <- ggplot(cross, aes(reorder(pair, rho), rho)) +
  geom_col(fill = BLUE, width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.5, size = 3.3) +
  coord_cartesian(ylim = c(0, 0.58)) +
  labs(x = NULL, y = "observed Spearman $\\rho$", title = "Proxy structure is reproducible across tissues") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

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
  annotate("text", x = 0.62, y = 6.3, label = "expected random (5.5)", size = 3.1,
           color = RED, hjust = 0) +
  labs(x = NULL, y = "motif hits per accessible peak",
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
  theme(legend.position = "top")

emit("fig1_truth_construct", (f1a | f1b) / (f1c | f1d) + plot_annotation(tag_levels = "A"),
     6.8, 5.8)

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
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))

f2b <- ggplot(dec, aes(reorder(pair, frac), frac)) +
  geom_col(fill = BLUE, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f\\%%", 100 * frac)), vjust = -0.5, size = 3.3) +
  coord_cartesian(ylim = c(0, 0.9)) +
  labs(x = NULL, y = "fraction of observed agreement",
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
  labs(x = NULL, y = "pairwise cell-type consensus $\\rho$",
       title = "Proxy is near cell-type-invariant") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

emit("fig2_cross_tissue_decomp", (f2a | f2b) / (f2c | f2d) + plot_annotation(tag_levels = "A"),
     6.8, 5.8)

# ==================== Figure 3: primary audit ====================
mk_forest <- function(df, title) {
  df$display <- factor(df$display, levels = rev(df$display))
  x_label <- if (df$is_baseline[1] %in% TRUE) 0.0165 else 0.0165
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
    coord_cartesian(xlim = c(-0.006, 0.0145), clip = "off") +
    labs(x = "partial Spearman $\\rho$", y = NULL, color = NULL, title = title) +
    theme(legend.position = "top", axis.ticks.y = element_blank(),
          plot.margin = margin(4, 118, 3, 4))
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
f3b <- mk_forest(nondeg, "Non-degree co-primary specification (no row supported by both nulls)")

f3c <- ggplot(primary %>% filter(!is_baseline), aes(mantel_q, degree_q, color = status)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.4) +
  scale_color_manual(values = support_colors) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "gene-label $q_M$ (log scale)", y = "row-shuffle $q_D$ (log scale)",
       color = NULL, title = "Both nulls below 0.05") +
  theme(legend.position = "top")

supp <- primary %>% filter(status == "supported by both nulls" & rho > 0)
base_cmp <- bind_rows(lapply(seq_len(nrow(supp)), function(i) {
  b <- baseline[baseline$tissue == supp$tissue[i], ]
  bind_rows(supp[i, ] %>% transmute(label = paste(model, tissue, sep = "\n"),
                                    kind = "FM row", rho = rho),
            data.frame(label = paste(supp$model[i], supp$tissue[i], sep = "\n"),
                       kind = "co-expression baseline", rho = b$rho[1]))
}))
f3d <- ggplot(base_cmp, aes(label, rho, fill = kind)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c(AQUA, "black")) +
  labs(x = NULL, y = "partial Spearman $\\rho$", fill = NULL,
       title = "Supported vs baseline") +
  theme(legend.position = "top", axis.text.x = element_text(size = 8.5))

emit("fig3_primary_audit", f3a / f3b / (f3c | f3d) + plot_annotation(tag_levels = "A") +
       plot_layout(heights = c(1.15, 1.0, 0.75)), 6.8, 9.2)

# ==================== Figure 4: usability check ====================
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
                           select(model_key, tissue, status, primary_rho = rho),
                         by = c("model_key", "tissue"))
usa$display <- factor(usa$display, levels = rev(usa$display))
f4a <- ggplot(usa, aes(fm_vs_coexp, display, color = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.5) +
  scale_color_manual(values = support_colors) +
  labs(x = "FM--co-expression Spearman $\\rho$ (usability check)", y = NULL, color = NULL,
       title = "Positives fail usability") +
  theme(legend.position = "top", axis.ticks.y = element_blank())

tile <- usa %>% transmute(display,
                          `supported positive` =
                            status == "supported by both nulls" & primary_rho > 0,
                          `passes usability ($\\rho>0$)` = fm_vs_coexp > 0)
tile_long <- bind_rows(
  tile %>% transmute(display, check = "supported positive", pass = `supported positive`),
  tile %>% transmute(display, check = "passes usability ($\\rho>0$)",
                     pass = `passes usability ($\\rho>0$)`))
f4b <- ggplot(tile_long, aes(check, display, fill = pass)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey80")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Support vs usability") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 12, hjust = 1),
        axis.ticks = element_blank())

ko_df <- data.frame(
  readout = rep(c("KO raw", "KO artifact-corrected"), each = 3),
  metric = rep(c("vs proxy (marginal)", "partial $|$ co-expression", "vs co-expression"), 2),
  rho = c(ko_stat$ko_vs_atac, ko_stat$ko_partial_given_coexp, ko_stat$ko_vs_coexp,
          ko_stat$ko_ctrl_vs_atac, ko_stat$ko_ctrl_partial_given_coexp, NA_real_))
ko_df$metric <- factor(ko_df$metric,
                       levels = c("vs proxy (marginal)", "partial $|$ co-expression",
                                  "vs co-expression"))
f4c <- ggplot(ko_df, aes(metric, rho, fill = readout)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = "grey55") +
  scale_fill_manual(values = c(MUTED, BLUE)) +
  labs(x = NULL, y = "Spearman $\\rho$", fill = NULL,
       title = "Knockout readout and its position-shift control") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 12, hjust = 1))

gf <- primary %>% filter(!is_baseline, model %in% c("Geneformer embed", "Geneformer attention"))
f4d <- ggplot(gf, aes(rho, model, group = tissue, color = tissue)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(color = "grey65") +
  geom_point(size = 2.6) +
  scale_color_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  labs(x = "partial Spearman $\\rho$ (full confounds)", y = NULL, color = NULL,
       title = "Two readouts of one model disagree in sign") +
  theme(legend.position = "top", axis.ticks.y = element_blank())

emit("fig4_usability_check", (f4a | f4b) / (f4c | f4d) + plot_annotation(tag_levels = "A"),
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
  coord_flip() +
  labs(x = NULL, y = "partial Spearman $\\rho$",
       title = "Row-shuffle null 95\\% band vs observed") +
  theme(axis.ticks.y = element_blank())

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
  scale_fill_manual(values = c(BLUE, YELLOW)) +
  labs(x = NULL, y = "Spearman $\\rho$", fill = NULL,
       title = "Attention alignment negative") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 18, hjust = 1, size = 8.5))

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
  geom_text(aes(y = observed, label = sprintf("$p=%.3f$", p)), hjust = -0.15, size = 3.1) +
  coord_flip() +
  labs(x = NULL, y = "partial $\\rho$ $|$ co-expression",
       title = "KO and control both align after co-expression alone") +
  theme(axis.ticks.y = element_blank())

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
  om %>% transmute(model, analysis = "authoritative audit", supported = authoritative),
  om %>% transmute(model, analysis = "attention-omission guard", supported = guard))
om_long$model <- factor(om_long$model, levels = rev(unique(om_long$model)))
f5d <- ggplot(om_long, aes(analysis, model, fill = supported)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey85")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Support decisions agree") +
  theme(legend.position = "top", axis.ticks = element_blank())

emit("fig5_null_diagnostics", (f5a | f5b) / (f5c | f5d) + plot_annotation(tag_levels = "A"),
     6.8, 6.0)

# ==================== Figure 6: specification dependence ====================
spec <- pooled %>% select(tissue, model, model_key, spec, rho, status)
spec$spec_label <- ifelse(spec$spec == "full", "Full", "Non-degree")
spec$display <- paste(spec$model, spec$tissue, sep = " -- ")
spec$display <- factor(spec$display, levels = rev(unique(spec$display)))
f6a <- ggplot(spec, aes(rho, display, group = display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(color = "grey65", linewidth = 0.45) +
  geom_point(aes(fill = spec_label), shape = 21, size = 2.2, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c("Full" = BLUE, "Non-degree" = YELLOW)) +
  labs(x = "partial Spearman $\\rho$", y = NULL, fill = NULL,
       title = "Specification changes effect direction") +
  theme(legend.position = "top", axis.ticks.y = element_blank())

sign_tile <- spec %>% mutate(sign = ifelse(rho > 0, "positive", "negative"))
f6b <- ggplot(sign_tile, aes(spec_label, display, fill = sign)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(positive = AQUA, negative = VIOLET)) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Sign-flip map") +
  theme(legend.position = "top", axis.ticks = element_blank())

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
  labs(x = "$p$ under full confounds (log)", y = "$p$ under non-degree (log)", color = NULL,
       title = "Randomization $p$ migrates upward") +
  theme(legend.position = "top")

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
                           labels = c("marginal", "coexp", "nondegree", "degree",
                                      "coexp+nondegree", "coexp+full")),
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
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(x = NULL, y = "partial Spearman $\\rho$", color = NULL,
       title = "Nested covariate ladder") +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 28, hjust = 1, size = 8))

emit("fig6_spec_sensitivity", (f6a | f6b) / (f6c | f6d) + plot_annotation(tag_levels = "A"),
     6.8, 7.0)

# ==================== Figure 7: per-cell-type descriptive ====================
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
  ptf %>% transmute(cell_type, series = "co-expression vs proxy", rho = coexp),
  ptf %>% transmute(cell_type, series = "Geneformer embed (partial)", rho = embed),
  ptf %>% transmute(cell_type, series = "Geneformer attention (partial)", rho = attn))
ptf$cell_type <- factor(ptf$cell_type, levels = ptf$cell_type)
ptf_long$cell_type <- factor(ptf_long$cell_type, levels = ptf$cell_type)
f7c <- ggplot(ptf_long, aes(cell_type, rho, color = series, group = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 2.0) +
  scale_color_manual(values = c("black", BLUE, AQUA)) +
  labs(x = NULL, y = "Spearman $\\rho$ (descriptive)", color = NULL,
       title = "Brain per-type readouts and the confound itself") +
  theme(legend.position = "top")

cells <- bind_rows(lapply(names(panel$pertype_n_cells), function(tn) {
  cc <- panel$pertype_n_cells[[tn]]
  bind_rows(lapply(names(cc), function(ct) {
    data.frame(tissue = ifelse(tn == "brain", "Brain", "PBMC"),
               cell_type = ct, n_cells = cc[[ct]])
  }))
}))
f7d <- ggplot(cells, aes(reorder(cell_type, n_cells), n_cells, fill = tissue)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  coord_flip() +
  labs(x = NULL, y = "cells", fill = NULL,
       title = "Cell counts behind the per-type rows") +
  theme(legend.position = "top")

emit("fig7_pertype_descriptive", f7a / f7b / (f7c | f7d) + plot_annotation(tag_levels = "A") +
       plot_layout(heights = c(1, 1, 0.9)), 6.8, 8.2)

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
f8a <- dose_panel(dose %>% filter(tissue == "Brain"),
                  if (is.null(subdiv)) NULL else subdiv %>% filter(tissue == "Brain"),
                  "Brain", BLUE)
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
  labs(x = "$\\alpha$-equivalent (log scale)", y = "observed partial $\\rho$", color = NULL,
       title = "Observed effects sit at the bottom of the ladder") +
  theme(legend.position = "top")

supp_keys <- primary %>% filter(status == "supported by both nulls") %>%
  mutate(key = paste(tolower(tissue), model_key, sep = "|")) %>% pull(key)
eff_supp <- eff_ok %>%
  mutate(key = paste(tolower(tissue), model, sep = "|")) %>%
  filter(key %in% supp_keys)
f8d <- ggplot(eff_supp, aes(reorder(display, alpha_equiv), alpha_equiv, fill = tissue)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.002, linetype = "dashed", color = RED) +
  annotate("text", x = 0.6, y = 0.0026, label = "smallest probed $\\alpha=0.002$",
           hjust = 0, size = 3.1, color = RED) +
  scale_fill_manual(values = c(Brain = BLUE, PBMC = AQUA)) +
  coord_flip() +
  labs(x = NULL, y = "$\\alpha$-equivalent of the observed effect", fill = NULL,
       title = "Injection equivalents") +
  theme(legend.position = "top", axis.text.y = element_text(size = 8))

emit("fig8_injection_ladder", (f8a | f8b) / (f8c | f8d) + plot_annotation(tag_levels = "A"),
     6.8, 6.4)

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
f9a <- ggplot(pa, aes(rho, family, color = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.4) +
  geom_text(aes(x = 0.0068,
                label = ifelse(status == "co-expression baseline",
                               sprintf("$q_M=%.3f$", q),
                               sprintf("$q_M=%.3f;\\ q_{\\mathrm{flip}}=%.3f$", q, flip_q))),
            hjust = 0, size = 4.0, color = "black") +
  scale_color_manual(values = c("supported vs baseline" = AQUA,
                                "not supported" = MUTED,
                                "co-expression baseline" = "black")) +
  coord_cartesian(xlim = c(-0.021, 0.0055), clip = "off") +
  labs(x = "adjusted test Spearman $\\rho$", y = NULL, color = NULL,
       title = "Supervised probe on held-out TFs (edge features only)") +
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
f9b <- ggplot(pb, aes(subset, rho, color = family, group = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
  scale_color_manual(values = c("Co-expression" = "black", "Geneformer embed" = BLUE,
                                "Geneformer attention" = AQUA, "scGPT encoder" = VIOLET,
                                "UCE encoder" = YELLOW, "Random-init floor" = MUTED)) +
  labs(x = NULL, y = "adjusted test Spearman $\\rho$", color = NULL,
       title = "Recovery collapses at the construction covariates") +
  # legend on top, two rows: a right-side legend eats half the row width and,
  # through patchwork's cross-row width alignment, squeezes the C/D row's plot
  # regions to zero width
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 28, hjust = 1, size = 9),
        legend.text = element_text(size = 9))

arms <- bind_rows(lapply(names(probe_fam_label), function(fam) {
  eo <- probe_eval$results$edge_only[[fam]]
  al <- probe_eval$results$all[[fam]]
  data.frame(family = unname(probe_fam_label[[fam]]),
             arm = rep(c("edge features only", "all features"), each = 1),
             rho = c(eo$adjusted_rho_mean, al$adjusted_rho_mean))
}))
arms$family <- factor(arms$family, levels = rev(unname(probe_fam_label)))
# horizontal bars: rotated x labels for the long family names collapse the plot
# region to a sliver when this panel shares a stacked patchwork row
f9c <- ggplot(arms, aes(rho, family, fill = arm)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_vline(xintercept = 0, color = "grey55") +
  scale_fill_manual(values = c(BLUE, MUTED)) +
  labs(y = NULL, x = "adjusted test Spearman $\\rho$", fill = NULL,
       title = "Degree features: no rescue") +
  theme(legend.position = "top", axis.ticks.y = element_blank())

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
  labs(x = NULL, y = "adjusted test Spearman $\\rho$ across TFs",
       title = "Per-TF recovery spread") +
  theme(axis.ticks.y = element_blank())

# four stacked full-width rows: in a three-row stack with (C | D) nested,
# patchwork's cross-row plot-region alignment squeezes the bottom row's plot
# regions to ~zero width (measured 59pt); full-width rows keep every plot region
# at the shared ~230pt column width
emit("fig9_tf_probe", f9a / f9b / f9c / f9d + plot_annotation(tag_levels = "A") +
       plot_layout(heights = c(1.05, 1.0, 0.85, 0.85)), 6.8, 8.4)

# ==================== Figure 10: coverage, scope, QC ====================
cov <- bind_rows(lapply(names(panel$readout_qc$manifest_gene_coverage), function(k) {
  data.frame(readout = k, covered = panel$readout_qc$manifest_gene_coverage[[k]])
}))
f10a <- ggplot(cov, aes(reorder(readout, covered), covered)) +
  geom_col(fill = BLUE, width = 0.6) +
  geom_text(aes(label = covered), hjust = -0.1, size = 3.2) +
  coord_flip(ylim = c(1100, 1230)) +
  labs(x = NULL, y = "manifest genes covered (of 1,200)",
       title = "Gene coverage by readout") +
  theme(axis.text.y = element_text(size = 9))

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
f10b <- ggplot(cov_long, aes(analysis, model, fill = present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey85")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Graph coverage") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1),
        axis.ticks = element_blank())

cells_r <- bind_rows(lapply(names(panel$readout_qc$cells_per_readout), function(k) {
  data.frame(readout = k, n = panel$readout_qc$cells_per_readout[[k]])
}))
f10c <- ggplot(cells_r, aes(reorder(readout, n), n)) +
  geom_col(fill = VIOLET, width = 0.6) +
  coord_flip() +
  labs(x = NULL, y = "cells",
       title = "Cells per readout") +
  theme(axis.text.y = element_text(size = 8.5))

ea <- panel$readout_qc$edge_accounting
edges <- data.frame(
  step = factor(c("446 TF $\\times$ 1,200 genes", "minus self pairs (PBMC edge set)",
                  "minus brain marker mask (brain edge set)"),
                levels = c("446 TF $\\times$ 1,200 genes", "minus self pairs (PBMC edge set)",
                           "minus brain marker mask (brain edge set)")),
  n = c(ea$tf_x_gene, ea$pbmc_edges, ea$brain_edges))
f10d <- ggplot(edges, aes(step, n)) +
  geom_col(fill = c("grey70", BLUE, AQUA), width = 0.6) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0, 580000)) +
  labs(x = NULL, y = "TF--gene edges",
       title = "Edge accounting") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 8.5))

emit("fig10_coverage_qc", (f10a | f10b) / (f10c | f10d) + plot_annotation(tag_levels = "A"),
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
             label = sprintf("%s / %s", format(c$relevant_peaks, big.mark = ","),
                             format(c$total_peaks, big.mark = ",")))
}))
f11a <- ggplot(cov11, aes(reorder(tissue, frac), frac, fill = tissue)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f\\%%", 100 * frac)), hjust = -0.1, size = 3.3) +
  coord_flip(ylim = c(0, 0.075)) +
  scale_fill_manual(values = TISSUE_COL) +
  labs(x = NULL, y = "fraction of all peaks linked",
       title = "$\\pm$2 kb linkage admits 4--6\\% of peaks")

ov <- bind_rows(lapply(tt$edge_overlap$regions, function(r) {
  data.frame(region = paste(unlist(r$combo), collapse = " + "), n = r$n)
}))
ov$region[ov$region == "Brain + PBMC + Fibroblast mix"] <- "all three"
ov$shared <- ov$region == "all three"
f11b <- ggplot(ov, aes(reorder(region, n), n, fill = shared)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3.1) +
  coord_flip(ylim = c(0, 47000)) +
  scale_fill_manual(values = c(`TRUE` = AQUA, `FALSE` = "grey70")) +
  labs(x = NULL, y = "non-self TF--gene edges",
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
  coord_cartesian(xlim = c(0.42, 0.56)) +
  labs(x = "cross-tissue agreement", y = NULL, color = NULL,
       title = "Rank agreement tracks support") +
  theme(legend.position = "top")

emit("fig11_third_tissue_transfer", (f11a | f11b) / (f11c | f11d) +
       plot_annotation(tag_levels = "A"), 6.8, 5.8)

# ==================== tables (unchanged) ====================
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

cross_rows <- vapply(seq_len(nrow(cross)), function(i) sprintf(
  "%s & $%.4f$ & 446 \\\\", cross$pair[i], cross$rho[i]), character(1))
writeLines(c(
  "\\begin{tabular}{lrr}", "\\toprule",
  "Tissue pair & Observed Spearman $\\rho$ & Shared TFs \\\\ ",
  "\\midrule", cross_rows, "\\bottomrule", "\\end{tabular}"
), file.path(figs, "table2_cross_tissue_observed.tex"))

cat("all authoritative figures and tables written to", figs, "\n")
