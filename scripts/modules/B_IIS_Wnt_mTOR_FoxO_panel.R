# Draw the final Fig. 2C pathway schematic.
# The layout follows the manually redrawn reference while correcting pathway
# directionality and separating established edges from indirect cross-talk.
# This module is sourced by B_differential_expression_pipeline.R.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(grid)
  library(svglite)
  library(openxlsx)
  library(officer)
  library(rvg)
})

if (!exists("dir_base", inherits = TRUE)) {
  if (!exists("PROJECT_ROOT", inherits = TRUE)) {
    stop("PROJECT_ROOT or dir_base must be defined before sourcing this module.")
  }
  dir_base <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
}
if (!exists("output_dir", inherits = TRUE)) output_dir <- file.path(dir_base, "output")
if (!exists("figure_dir", inherits = TRUE)) figure_dir <- file.path(output_dir, "figures")
cross_stage_output_dir <- file.path(output_dir, "cross_stage_concordance")
invisible(lapply(
  c(output_dir, figure_dir, cross_stage_output_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

rebuild_openxlsx_compat <- function(path) {
  script <- file.path(
    dir_base, "scripts", "utilities", "rebuild_openxlsx_workbook.py"
  )
  python_candidates <- c(
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    Sys.which("python"), Sys.which("python3"),
    file.path(
      Sys.getenv("USERPROFILE"), ".cache", "codex-runtimes",
      "codex-primary-runtime", "dependencies", "python", "python.exe"
    )
  )
  python_candidates <- unique(
    python_candidates[nzchar(python_candidates) & file.exists(python_candidates)]
  )
  if (!length(python_candidates) || !file.exists(script)) {
    warning("Could not rebuild the openxlsx workbook for wider XLSX compatibility")
    return(invisible(FALSE))
  }
  rebuilt <- paste0(path, ".compat.xlsx")
  status <- system2(
    python_candidates[[1]],
    c(shQuote(script), shQuote(path), shQuote(rebuilt))
  )
  if (!identical(status, 0L) || !file.exists(rebuilt)) {
    warning("XLSX compatibility rebuild failed for ", path)
    return(invisible(FALSE))
  }
  file.copy(rebuilt, path, overwrite = TRUE)
  unlink(rebuilt)
  invisible(TRUE)
}

key <- readRDS(file.path(
  dir_base, "input_differential_expression", "workspaces",
  "B_differential_expression_key_results.rds"
))

status_code <- function(effect, padj) {
  case_when(
    !is.finite(effect) | !is.finite(padj) ~ NA_integer_,
    padj < 0.05 & effect > 0 ~ 1L,
    padj < 0.05 & effect < 0 ~ -1L,
    TRUE ~ 0L
  )
}

de <- key$hog_de_shared %>% mutate(stage = as.character(stage))
vv_l2 <- de %>%
  filter(stage == "L2") %>%
  transmute(
    HOG, vv_log2FC = log2FC_ash_vv, vv_padj = padj_global_vv,
    vv_status = status_code(log2FC_ash_vv, padj_global_vv)
  )
pd_focal <- de %>%
  filter(stage %in% c("L1", "L4", "L5")) %>%
  transmute(
    HOG, pd_stage = stage, pd_log2FC = log2FC_ash_pd,
    pd_padj = padj_global_pd,
    pd_status = status_code(log2FC_ash_pd, padj_global_pd)
  )
strict_overlap <- inner_join(pd_focal, vv_l2, by = "HOG") %>%
  filter(pd_status != 0L, vv_status != 0L, pd_status == vv_status) %>%
  mutate(direction = if_else(pd_status > 0L, "concordant up", "concordant down"))

highlighted_hogs <- tribble(
  ~HOG, ~node, ~display_label, ~expected_pd_stages,
  "N13.HOG0004507", "cbp", "nej/CBP", "L1",
  "N13.HOG0005758", "wdr59", "Wdr59\n(GATOR2)", "L1",
  "N13.HOG0006185", "gsk3", "sgg/GSK3beta", "L4/L5",
  "N13.HOG0006961", "pepck", "Pepck", "L4",
  "N13.HOG0007424", "pdk1", "Pdk1", "L4",
  "N13.HOG0008912", "glut", "sut1-3 /\nSLC2-like", "L4",
  "N13.HOG0009263", "nlk", "nmo/NLK", "L4",
  "N13.HOG0010043", "gbs", "Gbs-76A", "L4"
)

observed_highlights <- strict_overlap %>%
  filter(HOG %in% highlighted_hogs$HOG) %>%
  group_by(HOG) %>%
  summarise(
    observed_pd_stages = paste(
      c("L1", "L4", "L5")[c("L1", "L4", "L5") %in% unique(pd_stage)],
      collapse = "/"
    ),
    direction = if_else(all(direction == "concordant up"), "up", "mixed"),
    .groups = "drop"
  )

highlight_check <- highlighted_hogs %>%
  left_join(observed_highlights, by = "HOG") %>%
  mutate(
    valid = !is.na(observed_pd_stages) &
      observed_pd_stages == expected_pd_stages & direction == "up"
  )
if (!all(highlight_check$valid)) {
  stop(
    "Fig. 2C highlight validation failed for: ",
    paste(highlight_check$HOG[!highlight_check$valid], collapse = ", ")
  )
}

colours <- list(
  ink = "#202020",
  line = "#454545",
  red = "#D62828",
  red_border = "#B71C1C",
  stage = "#C82127",
  wnt_fill = "#EAF5FC",
  wnt_border = "#CFE7F7",
  iis_fill = "#EDF8EE",
  iis_border = "#D2EBD4",
  mtor_fill = "#FFF9E9",
  mtor_border = "#F2E4B5",
  foxo_fill = "#FFF2F4",
  foxo_border = "#F6DADF",
  output_fill = "#EEE9FA",
  output_border = "#9B8BD1",
  foxo_node = "#F9C56A",
  foxo_border_node = "#E99337",
  white = "#FFFFFF"
)

nodes <- tribble(
  ~node, ~label, ~x, ~y, ~w, ~h, ~highlight, ~stage_label, ~shape,
  "wnt", "Wnt/wg", 0.065, 0.795, 0.066, 0.040, FALSE, "", "pill",
  "fz", "fz/fz2 + arr", 0.165, 0.760, 0.087, 0.045, FALSE, "", "box",
  "dsh", "dsh", 0.275, 0.760, 0.052, 0.045, FALSE, "", "box",
  "gsk3", "sgg/GSK3beta", 0.390, 0.760, 0.105, 0.052, TRUE,
  "Pd L4/L5 -> Vv L2", "box",
  "nlk", "nmo/NLK", 0.505, 0.760, 0.075, 0.052, TRUE,
  "Pd L4 -> Vv L2", "box",
  "ilp", "Ilp", 0.055, 0.485, 0.045, 0.038, FALSE, "", "pill",
  "inr", "InR", 0.095, 0.405, 0.045, 0.052, FALSE, "", "box",
  "chico", "chico/IRS", 0.190, 0.405, 0.080, 0.052, FALSE, "", "box",
  "pi3k", "Pi3K21B/92E", 0.305, 0.405, 0.100, 0.052, FALSE, "", "box",
  "pdk1", "Pdk1", 0.415, 0.405, 0.062, 0.052, TRUE,
  "Pd L4 -> Vv L2", "box",
  "akt", "Akt1", 0.505, 0.405, 0.058, 0.052, FALSE, "", "box",
  "tsc", "TSC1/2\nTsc1 + gig", 0.605, 0.700, 0.070, 0.068, FALSE, "", "box",
  "rheb", "Rheb", 0.680, 0.700, 0.052, 0.048, FALSE, "", "box",
  "tor", "Tor", 0.752, 0.700, 0.052, 0.048, FALSE, "", "box",
  "s6k", "S6k", 0.825, 0.700, 0.052, 0.048, FALSE, "", "box",
  "growth", "Translation\nand growth", 0.947, 0.700, 0.086, 0.075, FALSE, "", "output",
  "thor", "Thor/4E-BP", 0.770, 0.585, 0.086, 0.050, FALSE, "", "box",
  "eif4e", "eIF4E1", 0.770, 0.505, 0.065, 0.050, FALSE, "", "box",
  "wdr59", "Wdr59\n(GATOR2)", 0.680, 0.855, 0.080, 0.060, TRUE,
  "Pd L1 -> Vv L2", "box",
  "foxo", "foxo", 0.690, 0.230, 0.085, 0.060, FALSE, "", "foxo",
  "jnk", "bsk/JNK", 0.600, 0.165, 0.066, 0.047, FALSE, "", "box",
  "sir2", "Sir2", 0.625, 0.100, 0.055, 0.047, FALSE, "", "box",
  "cbp", "nej/CBP", 0.715, 0.105, 0.080, 0.052, TRUE,
  "Pd L1 -> Vv L2", "box",
  "pepck", "Pepck", 0.805, 0.300, 0.070, 0.052, TRUE,
  "Pd L4 -> Vv L2", "box",
  "bnip3", "BNIP3", 0.805, 0.215, 0.070, 0.048, FALSE, "", "box",
  "gbs", "Gbs-76A", 0.805, 0.145, 0.075, 0.052, TRUE,
  "Pd L4 -> Vv L2", "box",
  "glut", "sut1-3 /\nSLC2-like", 0.915, 0.295, 0.082, 0.060, TRUE,
  "Pd L4 -> Vv L2", "box",
  "carb", "Carbohydrate uptake,\nmetabolism and storage", 0.940, 0.145,
  0.100, 0.078, FALSE, "", "output"
)

node_lookup <- nodes %>% select(node, x, y, w, h)

node_edge <- function(node, side = c("left", "right", "top", "bottom")) {
  side <- match.arg(side)
  row <- node_lookup %>% filter(.data$node == .env$node)
  stopifnot(nrow(row) == 1)
  switch(
    side,
    left = c(row$x - row$w / 2, row$y),
    right = c(row$x + row$w / 2, row$y),
    top = c(row$x, row$y + row$h / 2),
    bottom = c(row$x, row$y - row$h / 2)
  )
}

draw_panel <- function(xmin, xmax, ymin, ymax, fill, border, title, title_colour) {
  grid.roundrect(
    x = unit((xmin + xmax) / 2, "npc"),
    y = unit((ymin + ymax) / 2, "npc"),
    width = unit(xmax - xmin, "npc"),
    height = unit(ymax - ymin, "npc"),
    r = unit(0.12, "in"),
    gp = gpar(fill = fill, col = border, lwd = 1.4)
  )
  grid.text(
    title, x = unit(xmin + 0.014, "npc"), y = unit(ymax - 0.025, "npc"),
    just = c("left", "top"),
    gp = gpar(fontfamily = "Arial", fontsize = 17, fontface = "bold", col = title_colour)
  )
}

draw_arrow <- function(x1, y1, x2, y2, dashed = FALSE, colour = colours$line,
                       lwd = 1.25, arrowhead = TRUE) {
  grid.lines(
    x = unit(c(x1, x2), "npc"), y = unit(c(y1, y2), "npc"),
    gp = gpar(col = colour, fill = colour, lwd = lwd, lty = if (dashed) 2 else 1),
    arrow = if (arrowhead) arrow(type = "closed", length = unit(0.085, "in")) else NULL
  )
}

draw_curve_arrow <- function(x1, y1, x2, y2, curvature = 0.15,
                             dashed = TRUE, colour = colours$line) {
  grid.curve(
    x1 = unit(x1, "npc"), y1 = unit(y1, "npc"),
    x2 = unit(x2, "npc"), y2 = unit(y2, "npc"),
    curvature = curvature, square = FALSE,
    gp = gpar(col = colour, fill = colour, lwd = 1.0, lty = if (dashed) 2 else 1),
    arrow = arrow(type = "closed", length = unit(0.075, "in"))
  )
}

draw_inhibition <- function(x1, y1, x2, y2, colour = colours$line, lwd = 1.35) {
  dx <- x2 - x1
  dy <- y2 - y1
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len
  uy <- dy / len
  bar_half <- 0.008
  line_end_x <- x2 - ux * 0.006
  line_end_y <- y2 - uy * 0.006
  grid.lines(
    x = unit(c(x1, line_end_x), "npc"), y = unit(c(y1, line_end_y), "npc"),
    gp = gpar(col = colour, lwd = lwd)
  )
  grid.lines(
    x = unit(c(x2 - uy * bar_half, x2 + uy * bar_half), "npc"),
    y = unit(c(y2 + ux * bar_half, y2 - ux * bar_half), "npc"),
    gp = gpar(col = colour, lwd = lwd + 0.35)
  )
}

draw_node <- function(row) {
  fill <- if (row$highlight) colours$red else colours$white
  border <- if (row$highlight) colours$red_border else colours$ink
  text_colour <- if (row$highlight) colours$white else colours$ink
  if (row$shape == "output") {
    fill <- colours$output_fill
    border <- colours$output_border
  }
  if (row$shape == "foxo") {
    fill <- colours$foxo_node
    border <- colours$foxo_border_node
  }
  radius <- if (row$shape %in% c("pill", "foxo")) 0.22 else 0.07
  grid.roundrect(
    x = unit(row$x, "npc"), y = unit(row$y, "npc"),
    width = unit(row$w, "npc"), height = unit(row$h, "npc"),
    r = unit(radius, "in"),
    gp = gpar(fill = fill, col = border, lwd = if (row$highlight) 1.3 else 1.15)
  )
  label_size <- case_when(
    row$node == "foxo" ~ 17,
    row$node %in% c("growth", "carb") ~ 10.2,
    str_detect(row$label, "\\n") ~ 10.2,
    nchar(row$label) > 12 ~ 10.2,
    TRUE ~ 11.3
  )
  grid.text(
    row$label, x = unit(row$x, "npc"), y = unit(row$y, "npc"),
    gp = gpar(
      fontfamily = "Arial", fontsize = label_size,
      fontface = if (row$highlight || row$node == "foxo") "bold" else "plain",
      col = text_colour
    )
  )
  if (nzchar(row$stage_label)) {
    grid.text(
      row$stage_label,
      x = unit(row$x, "npc"), y = unit(row$y - row$h / 2 - 0.018, "npc"),
      just = c("centre", "top"),
      gp = gpar(fontfamily = "Arial", fontsize = 8.8, fontface = "bold", col = colours$stage)
    )
  }
}

draw_membrane <- function() {
  # Minimal membrane cues avoid implying a vertebrate-specific receptor anatomy.
  for (yy in c(0.385, 0.425)) {
    grid.lines(
      x = unit(c(0.020, 0.075), "npc"), y = unit(c(yy, yy), "npc"),
      gp = gpar(col = "#91A2B5", lwd = 2.0)
    )
  }
  grid.lines(
    x = unit(c(0.895, 0.895), "npc"), y = unit(c(0.265, 0.325), "npc"),
    gp = gpar(col = "#91A2B5", lwd = 2.2)
  )
  grid.lines(
    x = unit(c(0.935, 0.935), "npc"), y = unit(c(0.265, 0.325), "npc"),
    gp = gpar(col = "#91A2B5", lwd = 2.2)
  )
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  draw_panel(0.012, 0.535, 0.630, 0.965, colours$wnt_fill, colours$wnt_border,
             "Wnt signalling", "#355D9A")
  draw_panel(0.552, 0.990, 0.435, 0.965, colours$mtor_fill, colours$mtor_border,
             "mTOR signalling", "#654018")
  draw_panel(0.012, 0.535, 0.285, 0.600, colours$iis_fill, colours$iis_border,
             "Insulin / insulin-like signalling", "#153D28")
  draw_panel(0.580, 0.990, 0.035, 0.405, colours$foxo_fill, colours$foxo_border,
             "FoxO signalling and metabolic processes", "#654018")

  grid.text(
    "C", x = unit(0.006, "npc"), y = unit(0.992, "npc"),
    just = c("left", "top"),
    gp = gpar(fontfamily = "Arial", fontsize = 24, fontface = "bold", col = "black")
  )

  draw_membrane()

  # Canonical Wnt sequence: Dishevelled inhibits GSK3beta; NLK is shown as a
  # separate Wnt-associated kinase because a direct GSK3beta -> NLK step is not
  # supported by the pathway topology.
  draw_arrow(0.098, 0.795, 0.121, 0.770)
  draw_arrow(0.209, 0.760, 0.249, 0.760)
  draw_inhibition(0.301, 0.760, 0.335, 0.760)
  draw_curve_arrow(0.300, 0.782, 0.468, 0.782, curvature = 0.22, dashed = TRUE)

  # IIS cascade and its two major downstream branches.
  draw_arrow(0.055, 0.466, 0.083, 0.432)
  draw_arrow(0.118, 0.405, 0.150, 0.405)
  draw_arrow(0.230, 0.405, 0.255, 0.405)
  draw_arrow(0.355, 0.405, 0.383, 0.405)
  draw_arrow(0.446, 0.405, 0.476, 0.405)
  draw_inhibition(0.505, 0.432, 0.585, 0.666)
  draw_inhibition(0.534, 0.395, 0.646, 0.246)

  # mTORC1 regulation and outputs.
  draw_inhibition(0.640, 0.700, 0.649, 0.700)
  draw_arrow(0.706, 0.700, 0.726, 0.700)
  draw_arrow(0.778, 0.700, 0.799, 0.700)
  draw_arrow(0.851, 0.700, 0.904, 0.700)
  draw_inhibition(0.752, 0.676, 0.758, 0.618)
  draw_inhibition(0.770, 0.560, 0.770, 0.538)
  draw_curve_arrow(0.800, 0.505, 0.904, 0.675, curvature = -0.20, dashed = FALSE)
  draw_curve_arrow(0.700, 0.827, 0.748, 0.724, curvature = 0.08, dashed = TRUE)

  # Conservatively represented pathway cross-talk.
  draw_curve_arrow(0.440, 0.744, 0.570, 0.690, curvature = -0.12, dashed = TRUE)
  draw_curve_arrow(0.430, 0.735, 0.662, 0.260, curvature = 0.18, dashed = TRUE)
  draw_curve_arrow(0.505, 0.734, 0.690, 0.260, curvature = -0.14, dashed = TRUE)
  draw_arrow(0.633, 0.177, 0.650, 0.205)
  draw_curve_arrow(0.648, 0.120, 0.667, 0.205, curvature = -0.08, dashed = TRUE)
  draw_curve_arrow(0.715, 0.132, 0.699, 0.199, curvature = 0.08, dashed = TRUE)

  # FoxO outputs and independent carbohydrate entry through an SLC2-like
  # transporter. No FoxO -> transporter edge is drawn.
  draw_arrow(0.731, 0.246, 0.770, 0.284)
  draw_arrow(0.733, 0.230, 0.770, 0.215)
  grid.text(
    "extracellular sugar", x = unit(0.915, "npc"), y = unit(0.355, "npc"),
    gp = gpar(fontfamily = "Arial", fontsize = 8.5, col = colours$line)
  )
  draw_arrow(0.915, 0.342, 0.915, 0.328)
  draw_arrow(0.915, 0.265, 0.930, 0.190)
  draw_arrow(0.843, 0.145, 0.887, 0.145)
  draw_curve_arrow(0.823, 0.275, 0.900, 0.178, curvature = -0.12, dashed = FALSE)

  nodes %>% split(seq_len(nrow(.))) %>% walk(~draw_node(.x))

  # Legend and analysis criterion.
  draw_arrow(0.020, 0.205, 0.060, 0.205)
  grid.text(
    "Activation", x = unit(0.070, "npc"), y = unit(0.205, "npc"),
    just = c("left", "centre"),
    gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = colours$ink)
  )
  draw_inhibition(0.020, 0.172, 0.060, 0.172)
  grid.text(
    "Inhibition", x = unit(0.070, "npc"), y = unit(0.172, "npc"),
    just = c("left", "centre"),
    gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = colours$ink)
  )
  draw_arrow(0.020, 0.139, 0.060, 0.139, dashed = TRUE)
  grid.text(
    "Cross-talk or indirect modulation", x = unit(0.070, "npc"), y = unit(0.139, "npc"),
    just = c("left", "centre"),
    gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = colours$ink)
  )
  grid.roundrect(
    x = unit(0.088, "npc"), y = unit(0.090, "npc"),
    width = unit(0.140, "npc"), height = unit(0.042, "npc"),
    r = unit(0.06, "in"),
    gp = gpar(fill = colours$red, col = colours$red_border, lwd = 1.0)
  )
  grid.text(
    "Concordantly upregulated", x = unit(0.088, "npc"), y = unit(0.090, "npc"),
    gp = gpar(fontfamily = "Arial", fontsize = 10.5, fontface = "bold", col = "white")
  )
  grid.text(
    "Global BH FDR < 0.05 in both contrasts; matching ASH-shrunken effect directions",
    x = unit(0.275, "npc"), y = unit(0.090, "npc"),
    just = c("left", "centre"),
    gp = gpar(fontfamily = "Arial", fontsize = 8.8, col = colours$line)
  )
}

png_file <- file.path(figure_dir, "Fig2C_source.png")
svg_file <- file.path(figure_dir, "Fig2C_source.svg")
pdf_file <- file.path(figure_dir, "Fig2C_source.pdf")
pptx_file <- file.path(figure_dir, "Fig2C_source.pptx")

png(
  png_file, width = 5000, height = 3500, res = 350,
  type = "cairo-png", bg = "white"
)
draw_figure()
dev.off()

svglite(svg_file, width = 14.286, height = 10, bg = "white")
draw_figure()
dev.off()

cairo_pdf(pdf_file, width = 14.286, height = 10, family = "Arial", bg = "white")
draw_figure()
dev.off()

write_tsv(
  highlight_check %>%
    mutate(across(where(is.character), ~ str_replace_all(.x, "[\\r\\n]+", " "))),
  file.path(cross_stage_output_dir, "Fig2C_highlighted_HOGs.tsv")
)
write_tsv(
  nodes %>%
    mutate(across(where(is.character), ~ str_replace_all(.x, "[\\r\\n]+", " "))),
  file.path(cross_stage_output_dir, "Fig2C_nodes.tsv")
)
saveRDS(
  list(
    highlighted_hogs = highlight_check,
    nodes = nodes,
    source_overlap = strict_overlap,
    output_files = c(PNG = png_file, SVG = svg_file, PDF = pdf_file, PPTX = pptx_file)
  ),
  file.path(cross_stage_output_dir, "Fig2C_final_data.rds")
)

wb <- createWorkbook()
addWorksheet(wb, "Highlighted_HOGs")
writeDataTable(wb, "Highlighted_HOGs", highlight_check, tableStyle = "TableStyleMedium2")
freezePane(wb, "Highlighted_HOGs", firstRow = TRUE)
addWorksheet(wb, "Strict_concordant_HOGs")
writeDataTable(wb, "Strict_concordant_HOGs", strict_overlap, tableStyle = "TableStyleMedium2")
freezePane(wb, "Strict_concordant_HOGs", firstRow = TRUE)
saveWorkbook(
  wb,
  file.path(cross_stage_output_dir, "IIS_Wnt_mTOR_FoxO_pathway_HOGs.xlsx"),
  overwrite = TRUE
)
rebuild_openxlsx_compat(
  file.path(cross_stage_output_dir, "IIS_Wnt_mTOR_FoxO_pathway_HOGs.xlsx")
)

doc <- read_pptx()
presentation_xml <- doc$presentation$get()
slide_size_node <- xml2::xml_find_first(
  presentation_xml, ".//p:sldSz", xml2::xml_ns(presentation_xml)
)
xml2::xml_set_attr(slide_size_node, "cx", as.character(round(13.33 * 914400)))
xml2::xml_set_attr(slide_size_node, "cy", as.character(round(7.5 * 914400)))
xml2::xml_set_attr(slide_size_node, "type", "screen16x9")
doc <- add_slide(doc, layout = "Blank", master = "Office Theme")
doc <- ph_with(
  doc, dml(code = draw_figure()),
  location = ph_location(left = 0, top = 0, width = 13.33, height = 7.5)
)
print(doc, target = pptx_file)

writeLines(
  c(
    paste("Completed", Sys.time()),
    "Eight highlighted HOGs validated against the final strict overlap table.",
    "All highlighted HOGs are concordantly upregulated.",
    "sut1-3/SLC2-like is drawn as an independent carbohydrate transporter, not as a FoxO target.",
    "Dishevelled inhibits GSK3beta; no direct GSK3beta-to-NLK activation is shown.",
    "AKT inhibition of TSC1/2 and FoxO is shown explicitly.",
    "Dashed arrows denote cross-talk or indirect modulation, not demonstrated regulation in these wasps."
  ),
  file.path(cross_stage_output_dir, "IIS_Wnt_mTOR_FoxO_pathway_VALIDATION.txt")
)

message("Final Fig. 2C written to: ", figure_dir)
