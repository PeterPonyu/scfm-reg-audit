# Figure typography and layout helpers for make_figs.R
# Unifies font packages, margins, and layout measurements with build_peerj_package.py

# Typography contract (matches figure_typography.py)
BASE_FONT_SIZE <- 11
MANUSCRIPT_FONT_SIZE <- 11
TARGET_WIDTH_INCHES <- 6.8
PEERJ_TEXT_BLOCK_INCHES <- 6.5

# Font packages: newtxtext/newtxmath (Times-compatible, matches PeerJ NimbusRomNo9L)
FONT_PACKAGES <- c("\\usepackage{amsmath}",
                   "\\usepackage{newtxtext}",
                   "\\usepackage{newtxmath}")

# tikzDevice options: 11pt base, font packages for metrics and rendering
setup_tikz_options <- function() {
  options(
    tikzDocumentDeclaration = sprintf("\\documentclass[%dpt]{article}", BASE_FONT_SIZE),
    tikzMetricPackages = c(getOption("tikzMetricPackages"), FONT_PACKAGES),
    tikzLatexPackages = c(getOption("tikzLatexPackages"), FONT_PACKAGES),
    tikzMetricsDictionary = ".tikz_metrics_pdftex"
  )
}

# Compute right margin for forest plots with q-value text columns
# Default accommodates "$q_M=0.000;\\ q_D=0.000$" at 11pt
compute_forest_margin <- function(qvalue_text_width_pt = 72, padding_pt = 12) {
  # Legacy reference values (empirical):
  # - fig3 panel A: margin(4, 118, 3, 4) -> ~118pt right
  # - fig9 panel A: margin(4, 132, 3, 4) -> ~132pt right
  #
  # Measured approach: qvalue_text_width_pt accounts for longest string,
  # padding_pt keeps text clear of plot region
  margin(4, qvalue_text_width_pt + padding_pt, 3, 4)
}

# Measure rendered label width using the active graphics font metrics.
# This intentionally measures the labels themselves rather than the axis region.
measure_label_width_pt <- function(labels, fontsize = BASE_FONT_SIZE - 1.5) {
  labels <- as.character(labels)
  if (!length(labels)) return(0)
  widths <- vapply(labels, function(label) {
    grob <- grid::textGrob(label, gp = grid::gpar(fontsize = fontsize))
    grid::convertWidth(grid::grobWidth(grob), unitTo = "pt", valueOnly = TRUE)
  }, numeric(1))
  max(widths, na.rm = TRUE)
}

# Compute margins for panels whose y tick labels and right-side annotations vary.
# ggplot2 reserves the measured y tick extent in the axis grob itself; only the
# amount beyond the normal 120pt label envelope needs an extra outer margin.
compute_measured_panel_margin <- function(ytick_labels = character(),
                                          right_labels = character(),
                                          fontsize = BASE_FONT_SIZE - 1.5,
                                          padding_pt = 12,
                                          base_left_pt = 4,
                                          normal_ytick_envelope_pt = 120) {
  left_width <- measure_label_width_pt(ytick_labels, fontsize = fontsize)
  right_width <- measure_label_width_pt(right_labels, fontsize = fontsize)
  ggplot2::margin(
    4,
    right_width + padding_pt,
    3,
    base_left_pt + max(0, left_width - normal_ytick_envelope_pt)
  )
}

# Compute margin adjustments for long y-tick labels
# Returns a ggplot2 margin() call with measured left padding
compute_ytick_margin <- function(max_label_width_pt = 60, base_left_pt = 4) {
  # For panels with long y-axis labels (e.g., model names, cell types),
  # measure the widest label and reserve space
  #
  # Standard margin: margin(4, 6, 3, 4) -> (top, right, bottom, left)
  # Adjust left to accommodate labels
  left_pt <- base_left_pt + max(0, max_label_width_pt - 40)  # 40pt baseline
  margin(4, 6, 3, left_pt)
}

# Standard theme with explicit typography contract
theme_screg <- function(base_size = BASE_FONT_SIZE) {
  theme_bw(base_size = base_size) + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = base_size, face = "plain", hjust = 0,
                              margin = margin(b = 3)),
    strip.text = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size - 1.5),
    legend.text = element_text(size = base_size - 1.5),
    axis.title = element_text(size = base_size - 1),
    axis.text = element_text(size = base_size - 1.5),
    plot.margin = margin(4, 6, 3, 4)  # default, override per panel as needed
  )
}

# Validate panel tags in emitted .tex file
# Prevents silent truncation (e.g., unescaped % turning rest of line into comment)
validate_panel_tags <- function(tex_path, expected_tags = LETTERS[1:4]) {
  lines <- readLines(tex_path, warn = FALSE)
  missing <- expected_tags[!vapply(expected_tags, function(t) {
    any(grepl(paste0("{", t, "}"), lines, fixed = TRUE))
  }, logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("emitted figure %s lost panel tags: %s",
                 basename(tex_path), paste(missing, collapse = ",")))
  }
}

# Unified emit function with tag validation
emit_figure <- function(name, plot, width, height, tags = LETTERS[1:4]) {
  base <- Sys.getenv("SCFM_BASE", "..")
  figs <- file.path(base, "paper", "figs")
  path <- file.path(figs, paste0(name, ".tex"))

  tikz(path, width = width, height = height, standAlone = FALSE, sanitize = FALSE)
  print(plot)
  dev.off()

  validate_panel_tags(path, tags)
  cat(name, "ok\n")
}
