# Cross-stage directional-concordance analyses for main Fig. 2 and Fig. S4.
#
# This module is sourced by B_differential_expression_pipeline.R after the
# final N13-HOG differential-expression results have been loaded. It can also
# be sourced independently when PROJECT_ROOT points to the repository root.

options(stringsAsFactors = FALSE)
set.seed(20260908)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(ordinal)
  library(topGO)
  library(clusterProfiler)
  library(glmnet)
  library(openxlsx)
  library(officer)
  library(rvg)
  library(patchwork)
})

if (!exists("dir_base", inherits = TRUE)) {
  if (!exists("PROJECT_ROOT", inherits = TRUE)) {
    stop("PROJECT_ROOT or dir_base must be defined before sourcing this module.")
  }
  dir_base <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
}
if (!exists("output_dir", inherits = TRUE)) {
  output_dir <- file.path(dir_base, "output")
}
if (!exists("figure_dir", inherits = TRUE)) {
  figure_dir <- file.path(output_dir, "figures")
}
if (!exists("supplemental_figure_dir", inherits = TRUE)) {
  supplemental_figure_dir <- file.path(output_dir, "supplemental_figures")
}
if (!exists("supplemental_table_dir", inherits = TRUE)) {
  supplemental_table_dir <- file.path(output_dir, "supplemental_tables")
}
cross_stage_output_dir <- file.path(output_dir, "cross_stage_concordance")
invisible(lapply(
  c(output_dir, figure_dir, supplemental_figure_dir, supplemental_table_dir,
    cross_stage_output_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

stages <- c("L1", "L2", "L3", "L4", "L5", "P")
bootstrap_B <- as.integer(Sys.getenv("B_CROSS_STAGE_BOOTSTRAP_B", "1000"))
selected_pairs <- tribble(
  ~comparison, ~pd_stage, ~vv_stage,
  "L1_to_L2", "L1", "L2",
  "L4_to_L2", "L4", "L2",
  "L5_to_L2", "L5", "L2",
  "P_to_P", "P", "P"
)

key_results_file <- file.path(
  dir_base, "input_differential_expression", "workspaces",
  "B_differential_expression_key_results.rds"
)
key <- if (file.exists(key_results_file)) readRDS(key_results_file) else list()
go_long <- readRDS(file.path(
  dir_base, "output_annotation", "key_dataframes_rds",
  "n13_hog_go_final_long.rds"
))
kegg_map <- readRDS(file.path(
  dir_base, "output_annotation", "key_dataframes_rds",
  "kegg_pathway_n13_hog_conserved_term2gene.rds"
))
kegg_names <- readRDS(file.path(
  dir_base, "output_annotation", "key_dataframes_rds",
  "kegg_pathway_term2name.rds"
))

hog_de_input <- if (exists("hog_de_shared", inherits = TRUE)) {
  get("hog_de_shared", inherits = TRUE)
} else {
  key$hog_de_shared
}
ortholog_input <- if (exists("orthologs_pls", inherits = TRUE)) {
  get("orthologs_pls", inherits = TRUE)
} else {
  key$orthologs_pls
}
if (is.null(hog_de_input) || is.null(ortholog_input)) {
  stop("The cached key results do not contain hog_de_shared and orthologs_pls.")
}

write_cross_stage_tsv <- function(x, filename) {
  readr::write_tsv(x, file.path(cross_stage_output_dir, filename), na = "")
}

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

save_editable_pptx <- function(plot, filename, width, height) {
  doc <- read_pptx()
  presentation_xml <- doc$presentation$get()
  slide_size_node <- xml2::xml_find_first(
    presentation_xml, ".//p:sldSz", xml2::xml_ns(presentation_xml)
  )
  xml2::xml_set_attr(
    slide_size_node, "cx", as.character(round(width * 914400))
  )
  xml2::xml_set_attr(
    slide_size_node, "cy", as.character(round(height * 914400))
  )
  if (abs(width / height - 16 / 9) < 0.02) {
    xml2::xml_set_attr(slide_size_node, "type", "screen16x9")
  }
  doc <- add_slide(doc, layout = "Blank", master = "Office Theme")
  doc <- ph_with(
    doc, dml(ggobj = plot),
    location = ph_location(left = 0, top = 0, width = width, height = height)
  )
  target <- if (identical(dirname(filename), ".")) {
    file.path(cross_stage_output_dir, filename)
  } else {
    filename
  }
  tryCatch(
    print(doc, target = target),
    error = function(e) {
      revised <- sub("\\.pptx$", "_revised.pptx", target, ignore.case = TRUE)
      warning(paste("Requested PPTX is locked; writing", basename(revised)))
      print(doc, target = revised)
    }
  )
}

parse_topgo_p <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  lt <- str_detect(x, "^<")
  out[lt] <- suppressWarnings(as.numeric(str_remove(x[lt], "^<\\s*")))
  out
}

de_long <- as_tibble(hog_de_input) %>%
  mutate(stage = as.character(stage)) %>%
  filter(stage %in% stages) %>%
  select(
    HOG, stage,
    log2FC_ash_pd, lfcSE_raw_pd, padj_global_pd,
    log2FC_ash_vv, lfcSE_raw_vv, padj_global_vv
  )

status_code <- function(effect, padj) {
  case_when(
    !is.finite(effect) | !is.finite(padj) ~ NA_integer_,
    padj < 0.05 & effect > 0 ~ 1L,
    padj < 0.05 & effect < 0 ~ -1L,
    TRUE ~ 0L
  )
}

status_long <- de_long %>%
  transmute(
    HOG, stage,
    pd_status = status_code(log2FC_ash_pd, padj_global_pd),
    vv_status = status_code(log2FC_ash_vv, padj_global_vv)
  )

get_pair <- function(pd_stage, vv_stage) {
  pd <- de_long %>%
    filter(stage == pd_stage) %>%
    transmute(
      HOG, pd_stage = pd_stage,
      pd_log2FC = log2FC_ash_pd, pd_SE = lfcSE_raw_pd,
      pd_padj = padj_global_pd,
      pd_status = status_code(log2FC_ash_pd, padj_global_pd)
    )
  vv <- de_long %>%
    filter(stage == vv_stage) %>%
    transmute(
      HOG, vv_stage = vv_stage,
      vv_log2FC = log2FC_ash_vv, vv_SE = lfcSE_raw_vv,
      vv_padj = padj_global_vv,
      vv_status = status_code(log2FC_ash_vv, padj_global_vv)
    )
  inner_join(pd, vv, by = "HOG") %>%
    filter(!is.na(pd_status), !is.na(vv_status)) %>%
    mutate(
      pd_DE = pd_status != 0,
      vv_DE = vv_status != 0,
      either_DE = pd_DE | vv_DE,
      both_DE = pd_DE & vv_DE,
      concordant_DE = both_DE & pd_status == vv_status,
      discordant_DE = both_DE & pd_status != vv_status,
      direction = case_when(
        concordant_DE & pd_status == 1 ~ "concordant up",
        concordant_DE & pd_status == -1 ~ "concordant down",
        discordant_DE ~ "discordant",
        TRUE ~ "not jointly significant"
      ),
      concordance_z = sign(pd_log2FC * vv_log2FC) *
        pmin(abs(pd_log2FC / pd_SE), abs(vv_log2FC / vv_SE))
    )
}

# ---------------------------------------------------------------------------
# A. Pairwise overlap and directional-concordance tests
# ---------------------------------------------------------------------------

stage_pairs <- expand_grid(pd_stage = stages, vv_stage = stages)
pair_tables <- pmap(stage_pairs, get_pair)

fit_pair_ordinal <- function(dat) {
  dat <- dat %>%
    mutate(vv_ordered = ordered(vv_status, levels = c(-1, 0, 1)))
  fit <- try(
    clm(vv_ordered ~ pd_status, data = dat, link = "logit", Hess = TRUE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(tibble(
      ordinal_log_odds = NA_real_, ordinal_SE = NA_real_,
      ordinal_z = NA_real_, ordinal_p = NA_real_,
      proportional_odds_p = NA_real_
    ))
  }
  co <- coef(summary(fit))
  row <- co[rownames(co) == "pd_status", , drop = FALSE]
  nominal <- try(nominal_test(fit), silent = TRUE)
  nominal_p <- NA_real_
  if (!inherits(nominal, "try-error") && "Pr(>Chi)" %in% colnames(nominal)) {
    rn <- rownames(nominal)
    hit <- which(rn == "pd_status")
    if (length(hit)) nominal_p <- nominal[hit[1], "Pr(>Chi)"]
  }
  tibble(
    ordinal_log_odds = row[1, "Estimate"],
    ordinal_SE = row[1, "Std. Error"],
    ordinal_z = row[1, "z value"],
    ordinal_p = row[1, "Pr(>|z|)"],
    proportional_odds_p = nominal_p
  )
}

pairwise_status_results <- map_dfr(pair_tables, function(dat) {
  N <- nrow(dat)
  pd_up <- sum(dat$pd_status == 1)
  pd_down <- sum(dat$pd_status == -1)
  vv_up <- sum(dat$vv_status == 1)
  vv_down <- sum(dat$vv_status == -1)
  expected_concordant <- (pd_up * vv_up + pd_down * vv_down) / N
  bind_cols(
    tibble(
      pd_stage = unique(dat$pd_stage),
      vv_stage = unique(dat$vv_stage),
      universe_n = N,
      pd_DE_n = sum(dat$pd_DE),
      vv_DE_n = sum(dat$vv_DE),
      either_DE_n = sum(dat$either_DE),
      jointly_DE_n = sum(dat$both_DE),
      concordant_DE_n = sum(dat$concordant_DE),
      discordant_DE_n = sum(dat$discordant_DE),
      expected_concordant_DE_n = expected_concordant,
      concordant_DE_enrichment = sum(dat$concordant_DE) / expected_concordant
    ),
    fit_pair_ordinal(dat)
  )
}) %>%
  mutate(
    ordinal_p_BH_36 = p.adjust(ordinal_p, "BH"),
    significance_FDR = case_when(
      ordinal_p_BH_36 < 0.001 ~ "***",
      ordinal_p_BH_36 < 0.01 ~ "**",
      ordinal_p_BH_36 < 0.05 ~ "*",
      ordinal_p_BH_36 < 0.10 ~ ".",
      TRUE ~ ""
    ),
    cell_label = paste0(
      concordant_DE_n, "/", jointly_DE_n,
      ifelse(significance_FDR == "", "",
             paste0("\n", significance_FDR))
    ),
    pd_stage = factor(pd_stage, levels = stages),
    vv_stage = factor(vv_stage, levels = stages)
  )

summarise_directional_concordance <- function(dat) {
  both <- dat %>% filter(both_DE)
  tab <- table(
    factor(both$pd_status, levels = c(-1, 1)),
    factor(both$vv_status, levels = c(-1, 1))
  )
  names(dimnames(tab)) <- c("Polistes", "Vespula")
  down_down <- unname(tab[1, 1])
  down_up <- unname(tab[1, 2])
  up_down <- unname(tab[2, 1])
  up_up <- unname(tab[2, 2])
  jointly_de <- sum(tab)
  concordant <- down_down + up_up
  discordant <- down_up + up_down
  pd_up <- up_down + up_up
  pd_down <- down_down + down_up
  vv_up <- down_up + up_up
  vv_down <- down_down + up_down
  expected_fraction <- if (jointly_de > 0) {
    (pd_up * vv_up + pd_down * vv_down) / jointly_de^2
  } else {
    NA_real_
  }
  fisher_fit <- if (jointly_de > 0) {
    fisher.test(tab, alternative = "greater")
  } else {
    NULL
  }
  log2_or_ha <- log2(
    ((up_up + 0.5) * (down_down + 0.5)) /
      ((up_down + 0.5) * (down_up + 0.5))
  )
  tibble(
    pd_stage = unique(dat$pd_stage),
    vv_stage = unique(dat$vv_stage),
    universe_n = nrow(dat),
    pd_DE_n = sum(dat$pd_DE),
    vv_DE_n = sum(dat$vv_DE),
    jointly_DE_n = jointly_de,
    concordant_DE_n = concordant,
    discordant_DE_n = discordant,
    up_up_n = up_up,
    down_down_n = down_down,
    up_down_n = up_down,
    down_up_n = down_up,
    concordance_fraction = concordant / jointly_de,
    expected_concordance_fraction_from_margins = expected_fraction,
    excess_concordance = concordant / jointly_de - expected_fraction,
    concordance_log2_odds_ratio_HA = log2_or_ha,
    positive_concordance_log2_odds_ratio = pmax(log2_or_ha, 0),
    fisher_one_sided_p_greater = if (is.null(fisher_fit)) {
      NA_real_
    } else {
      fisher_fit$p.value
    },
    fisher_conditional_odds_ratio = if (is.null(fisher_fit)) {
      NA_real_
    } else {
      unname(fisher_fit$estimate)
    }
  )
}

pairwise_concordance_results <- map_dfr(
  pair_tables, summarise_directional_concordance
) %>%
  mutate(
    fisher_p_FDR_36 = p.adjust(fisher_one_sided_p_greater, method = "BH"),
    significance_FDR = case_when(
      fisher_p_FDR_36 < 0.001 ~ "***",
      fisher_p_FDR_36 < 0.01 ~ "**",
      fisher_p_FDR_36 < 0.05 ~ "*",
      fisher_p_FDR_36 < 0.10 ~ ".",
      TRUE ~ ""
    ),
    cell_label = paste0(
      concordant_DE_n, "/", jointly_DE_n,
      if_else(significance_FDR == "", "",
              paste0("\n", significance_FDR))
    ),
    pd_stage = factor(pd_stage, levels = stages),
    vv_stage = factor(vv_stage, levels = stages)
  )

pairwise_fill_max <- max(
  pairwise_concordance_results$positive_concordance_log2_odds_ratio,
  na.rm = TRUE
)
p_pairwise_concordance <- ggplot(
  pairwise_concordance_results,
  aes(pd_stage, vv_stage, fill = positive_concordance_log2_odds_ratio)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(
    aes(label = cell_label), family = "Arial", size = 3.55,
    lineheight = 0.82
  ) +
  scale_fill_gradient(
    low = "white", high = "#D7191C",
    limits = c(0, pairwise_fill_max), oob = squish
  ) +
  coord_fixed() +
  labs(
    x = expression(italic("Polistes dominula")*" stage (late vs early)"),
    y = expression(italic("Vespula vulgaris")*" stage (queen vs worker)"),
    fill = expression("Concordance log"[2]*" odds ratio")
  ) +
  guides(fill = guide_colourbar(title.position = "top")) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(size = 12.5),
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(0.9, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

p_pairwise_ordinal <- ggplot(
  pairwise_status_results,
  aes(pd_stage, vv_stage, fill = ordinal_log_odds)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = cell_label), family = "Arial", size = 3.0,
            lineheight = 0.82) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    limits = c(-2.5, 2.5), oob = squish
  ) +
  coord_fixed() +
  labs(
    x = expression(italic("Polistes dominula")*" stage (late vs early)"),
    y = expression(italic("Vespula vulgaris")*" stage (queen vs worker)"),
    fill = "Ordinal\nlog odds"
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# ---------------------------------------------------------------------------
# B. Multivariable cumulative-link models: all Polistes stages together
# ---------------------------------------------------------------------------

status_wide <- status_long %>%
  pivot_wider(
    names_from = stage, values_from = c(pd_status, vv_status),
    names_sep = "_"
  )

fit_full_clm <- function(vv_stage) {
  response <- paste0("vv_status_", vv_stage)
  predictors <- paste0("pd_status_", stages)
  dat <- status_wide %>%
    select(HOG, all_of(response), all_of(predictors)) %>%
    filter(if_all(all_of(c(response, predictors)), ~ !is.na(.x))) %>%
    rename(y = all_of(response)) %>%
    mutate(y = ordered(y, levels = c(-1, 0, 1)))
  form <- reformulate(predictors, response = "y")
  fit <- clm(
    form, data = dat, link = "logit", Hess = TRUE,
    control = clm.control(maxIter = 200, gradTol = 1e-6)
  )
  max_gradient <- if (length(fit$maxGradient)) fit$maxGradient[[1]] else NA_real_
  condition_number <- if (length(fit$cond.H)) fit$cond.H[[1]] else NA_real_
  co <- as.data.frame(coef(summary(fit))) %>%
    rownames_to_column("term") %>%
    filter(term %in% predictors) %>%
    transmute(
      vv_stage = vv_stage,
      pd_stage = str_remove(term, "^pd_status_"),
      estimate = Estimate,
      SE = `Std. Error`,
      z = `z value`,
      p_value = `Pr(>|z|)`,
      n = nrow(dat),
      max_gradient = max_gradient,
      condition_number = condition_number
    )
  nominal <- try(as.data.frame(nominal_test(fit)) %>%
                   rownames_to_column("term"), silent = TRUE)
  if (!inherits(nominal, "try-error")) {
    pcol <- intersect(c("Pr(>Chi)", "Pr(>Chisq)"), names(nominal))
    if (length(pcol)) {
      nominal <- nominal %>%
        filter(term %in% predictors) %>%
        transmute(
          pd_stage = str_remove(term, "^pd_status_"),
          proportional_odds_p = .data[[pcol[1]]]
        )
      co <- left_join(co, nominal, by = "pd_stage")
    }
  }
  co
}

full_clm_results <- map_dfr(stages, fit_full_clm) %>%
  mutate(
    p_FDR_36 = p.adjust(p_value, "BH"),
    pd_stage = factor(pd_stage, levels = stages),
    vv_stage = factor(vv_stage, levels = stages),
    significance = case_when(
      p_FDR_36 < 0.001 ~ "***",
      p_FDR_36 < 0.01 ~ "**",
      p_FDR_36 < 0.05 ~ "*",
      p_FDR_36 < 0.10 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  left_join(
    pairwise_status_results %>%
      transmute(
        pd_stage = as.character(pd_stage),
        vv_stage = as.character(vv_stage),
        jointly_DE_n, concordant_DE_n,
        count_label = paste0(concordant_DE_n, "/", jointly_DE_n)
      ),
    by = c("pd_stage", "vv_stage")
  ) %>%
  mutate(
    pd_stage = factor(pd_stage, levels = stages),
    vv_stage = factor(vv_stage, levels = stages),
    panel_label = if_else(
      significance == "", count_label,
      paste0(count_label, "\n", significance)
    )
  )

p_full_clm <- ggplot(full_clm_results, aes(pd_stage, vv_stage, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = panel_label), family = "Arial", size = 3.15,
            lineheight = 0.82) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    oob = squish
  ) +
  coord_fixed() +
  labs(
    x = expression(italic("Polistes dominula")*" stage predictor"),
    y = expression(italic("Vespula vulgaris")*" stage response"),
    fill = "Partial ordinal\nlog odds"
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text = element_text(colour = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# ---------------------------------------------------------------------------
# C. Continuous-effect nonnegative ridge with all Polistes stages
# ---------------------------------------------------------------------------

vv_de_any <- de_long %>%
  group_by(HOG) %>%
  summarise(keep = any(padj_global_vv < 0.05, na.rm = TRUE), .groups = "drop") %>%
  filter(keep) %>%
  pull(HOG)

effect_wide <- de_long %>%
  filter(HOG %in% vv_de_any) %>%
  select(HOG, stage, log2FC_ash_pd, log2FC_ash_vv, lfcSE_raw_vv) %>%
  pivot_wider(
    names_from = stage,
    values_from = c(log2FC_ash_pd, log2FC_ash_vv, lfcSE_raw_vv),
    names_sep = "_"
  )

fit_all_stage_ridge <- function(vv_stage) {
  x_cols <- paste0("log2FC_ash_pd_", stages)
  y_col <- paste0("log2FC_ash_vv_", vv_stage)
  se_col <- paste0("lfcSE_raw_vv_", vv_stage)
  dat <- effect_wide %>%
    select(HOG, all_of(c(x_cols, y_col, se_col))) %>%
    filter(if_all(all_of(c(x_cols, y_col, se_col)), is.finite))
  X <- as.matrix(dat[, x_cols])
  colnames(X) <- stages
  y <- dat[[y_col]]
  w <- 1 / (dat[[se_col]]^2 + 1e-8)
  set.seed(20260907 + match(vv_stage, stages))
  foldid <- sample(rep(seq_len(10), length.out = nrow(dat)))
  fit <- cv.glmnet(
    X, y, weights = w, family = "gaussian", alpha = 0,
    lower.limits = 0, standardize = FALSE, intercept = TRUE,
    foldid = foldid, type.measure = "mse"
  )
  b <- as.matrix(coef(fit, s = "lambda.min"))
  beta_hat <- as.numeric(b[stages, 1])

  # Pairs bootstrap at the selected penalty. This is deliberately the same
  # stability summary used for the constrained coefficients elsewhere in the
  # exploratory analysis; it is not an unconstrained Wald test.
  boot <- matrix(NA_real_, nrow = bootstrap_B, ncol = length(stages))
  colnames(boot) <- stages
  for (i in seq_len(bootstrap_B)) {
    idx <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    fb <- try(
      glmnet(
        X[idx, , drop = FALSE], y[idx], weights = w[idx],
        family = "gaussian", alpha = 0, lambda = fit$lambda.min,
        lower.limits = 0, standardize = FALSE, intercept = TRUE
      ),
      silent = TRUE
    )
    if (!inherits(fb, "try-error")) {
      bb <- as.matrix(coef(fb, s = fit$lambda.min))
      boot[i, ] <- as.numeric(bb[stages, 1])
    }
  }
  boot_se <- apply(boot, 2, sd, na.rm = TRUE)
  boot_low <- apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  boot_high <- apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  selection_probability <- colMeans(boot > 1e-10, na.rm = TRUE)
  bootstrap_valid_n <- colSums(is.finite(boot))
  bootstrap_boundary_n <- colSums(boot <= 1e-10, na.rm = TRUE)
  bootstrap_p_one_sided <-
    (bootstrap_boundary_n + 1) / (bootstrap_valid_n + 1)
  tibble(
    vv_stage = vv_stage,
    pd_stage = stages,
    estimate = beta_hat,
    bootstrap_SE = boot_se,
    bootstrap_CI_low = boot_low,
    bootstrap_CI_high = boot_high,
    bootstrap_selection_probability = selection_probability,
    bootstrap_boundary_n = bootstrap_boundary_n,
    bootstrap_valid_n = bootstrap_valid_n,
    bootstrap_p_one_sided = bootstrap_p_one_sided,
    lambda_min = fit$lambda.min,
    cv_mse_min = min(fit$cvm),
    n = nrow(dat)
  )
}

all_stage_ridge_results <- map_dfr(stages, fit_all_stage_ridge) %>%
  mutate(
    bootstrap_p_FDR_36 = p.adjust(bootstrap_p_one_sided, "BH"),
    pd_stage = factor(pd_stage, levels = stages),
    vv_stage = factor(vv_stage, levels = stages),
    significance = case_when(
      bootstrap_p_FDR_36 < 0.001 ~ "***",
      bootstrap_p_FDR_36 < 0.01 ~ "**",
      bootstrap_p_FDR_36 < 0.05 ~ "*",
      bootstrap_p_FDR_36 < 0.10 ~ ".",
      TRUE ~ ""
    )
  )

ridge_fill_max <- max(all_stage_ridge_results$estimate, na.rm = TRUE)

p_all_stage_ridge <- ggplot(
  all_stage_ridge_results, aes(pd_stage, vv_stage, fill = estimate)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = significance), family = "Arial", size = 4.8) +
  scale_fill_gradient(
    low = "white", high = "#D7191C", limits = c(0, ridge_fill_max),
    oob = squish
  ) +
  coord_fixed() +
  labs(
    x = expression(italic("Polistes dominula")*" stage predictor"),
    y = expression(italic("Vespula vulgaris")*" stage response"),
    fill = "Nonnegative\npartial slope"
  ) +
  guides(fill = guide_colourbar(title.position = "top")) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text = element_text(colour = "black"),
    axis.title = element_text(size = 12.5),
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.5),
    legend.position = "top", legend.direction = "horizontal",
    legend.key.width = unit(2.0, "cm"),
    legend.title = element_text(size = 9), legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

tryCatch(
  ggsave(
    file.path(cross_stage_output_dir, "supplement_all_stage_nonnegative_ridge.pdf"),
    p_all_stage_ridge, width = 7.0, height = 6.4, device = cairo_pdf
  ),
  error = function(e) {
    warning("Base-name ridge PDF is locked; writing the revised PDF with a suffix")
    ggsave(
      file.path(cross_stage_output_dir, "supplement_all_stage_nonnegative_ridge_bootstrap.pdf"),
      p_all_stage_ridge, width = 7.0, height = 6.4, device = cairo_pdf
    )
  }
)
ggsave(file.path(cross_stage_output_dir, "supplement_all_stage_nonnegative_ridge.png"),
       p_all_stage_ridge, width = 7.0, height = 6.4, dpi = 400, bg = "white")
save_editable_pptx(
  p_all_stage_ridge, "supplement_all_stage_nonnegative_ridge.pptx", 7.0, 6.4
)

# ---------------------------------------------------------------------------
# C2. Exploratory nonnegative cumulative-logit model with pairs bootstrap
# ---------------------------------------------------------------------------

fit_nonnegative_clm <- function(X, y, start = NULL) {
  y <- as.integer(y)
  keep <- is.finite(y) & rowSums(!is.finite(X)) == 0
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  p <- ncol(X)

  nll <- function(par) {
    alpha1 <- par[1]
    alpha2 <- alpha1 + exp(par[2])
    beta <- par[2 + seq_len(p)]
    eta <- drop(X %*% beta)
    c1 <- plogis(alpha1 - eta)
    c2 <- plogis(alpha2 - eta)
    probs <- cbind(c1, c2 - c1, 1 - c2)
    probs <- pmax(probs, 1e-14)
    -sum(log(probs[cbind(seq_along(y), y)]))
  }

  if (is.null(start)) {
    p1 <- mean(y == 1)
    p12 <- mean(y <= 2)
    alpha1 <- qlogis(pmin(pmax(p1, 1e-4), 1 - 1e-4))
    alpha2 <- qlogis(pmin(pmax(p12, 1e-4), 1 - 1e-4))
    start <- c(alpha1, log(max(alpha2 - alpha1, 0.1)), rep(0.05, p))
  }

  fit <- nlminb(
    start = start,
    objective = nll,
    lower = c(-20, log(1e-4), rep(0, p)),
    upper = c(20, log(40), rep(10, p)),
    control = list(iter.max = 300, eval.max = 600, rel.tol = 1e-9)
  )
  list(
    par = fit$par,
    beta = setNames(fit$par[2 + seq_len(p)], colnames(X)),
    logLik = -fit$objective,
    convergence = fit$convergence,
    message = fit$message
  )
}

fit_bootstrap_nonnegative_clm <- function(vv_stage) {
  response <- paste0("vv_status_", vv_stage)
  predictors <- paste0("pd_status_", stages)
  dat <- status_wide %>%
    select(HOG, all_of(response), all_of(predictors)) %>%
    filter(if_all(all_of(c(response, predictors)), ~ !is.na(.x)))
  X <- as.matrix(dat[, predictors, drop = FALSE])
  colnames(X) <- stages
  y <- match(dat[[response]], c(-1, 0, 1))

  fit <- fit_nonnegative_clm(X, y)
  set.seed(20261900 + match(vv_stage, stages))
  boot <- matrix(NA_real_, nrow = bootstrap_B, ncol = length(stages))
  colnames(boot) <- stages
  boot_converged <- logical(bootstrap_B)
  for (i in seq_len(bootstrap_B)) {
    idx <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    fb <- try(
      fit_nonnegative_clm(X[idx, , drop = FALSE], y[idx], start = fit$par),
      silent = TRUE
    )
    if (!inherits(fb, "try-error") && fb$convergence == 0) {
      boot[i, ] <- fb$beta
      boot_converged[i] <- TRUE
    }
  }

  tibble(
    vv_stage = vv_stage,
    pd_stage = stages,
    estimate = as.numeric(fit$beta[stages]),
    bootstrap_SE = apply(boot, 2, sd, na.rm = TRUE),
    bootstrap_CI_low = apply(
      boot, 2, quantile, probs = 0.025, na.rm = TRUE
    ),
    bootstrap_CI_high = apply(
      boot, 2, quantile, probs = 0.975, na.rm = TRUE
    ),
    bootstrap_selection_probability = colMeans(boot > 1e-8, na.rm = TRUE),
    bootstrap_boundary_n = colSums(boot <= 1e-8, na.rm = TRUE),
    bootstrap_valid_n = colSums(is.finite(boot)),
    bootstrap_p_one_sided =
      (colSums(boot <= 1e-8, na.rm = TRUE) + 1) /
      (colSums(is.finite(boot)) + 1),
    model_logLik = fit$logLik,
    model_convergence = fit$convergence,
    bootstrap_converged = sum(boot_converged),
    bootstrap_B = bootstrap_B,
    n = nrow(dat)
  )
}

nonnegative_clm_results <- map_dfr(
  stages, fit_bootstrap_nonnegative_clm
) %>%
  mutate(
    bootstrap_p_FDR_36 = p.adjust(bootstrap_p_one_sided, "BH")
  ) %>%
  left_join(
    pairwise_status_results %>%
      transmute(
        pd_stage = as.character(pd_stage),
        vv_stage = as.character(vv_stage),
        jointly_DE_n, concordant_DE_n,
        count_label = paste0(concordant_DE_n, "/", jointly_DE_n)
      ),
    by = c("pd_stage", "vv_stage")
  ) %>%
  mutate(
    pd_stage = factor(pd_stage, levels = stages),
    vv_stage = factor(vv_stage, levels = stages),
    significance = case_when(
      bootstrap_p_FDR_36 < 0.001 ~ "***",
      bootstrap_p_FDR_36 < 0.01 ~ "**",
      bootstrap_p_FDR_36 < 0.05 ~ "*",
      bootstrap_p_FDR_36 < 0.10 ~ ".",
      TRUE ~ ""
    ),
    panel_label = if_else(
      significance == "", count_label,
      paste0(count_label, "\n", significance)
    )
  )

p_nonnegative_clm <- ggplot(
  nonnegative_clm_results, aes(pd_stage, vv_stage, fill = estimate)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = panel_label), family = "Arial", size = 3.75,
            lineheight = 0.82) +
  scale_fill_gradient(low = "white", high = "#D7191C", oob = squish) +
  coord_fixed() +
  labs(
    x = expression(italic("Polistes dominula")*" stage predictor"),
    y = expression(italic("Vespula vulgaris")*" stage response"),
    fill = "Nonnegative partial\nordinal log odds"
  ) +
  guides(fill = guide_colourbar(title.position = "top")) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text = element_text(colour = "black"),
    axis.title = element_text(size = 12.5),
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.5),
    legend.position = "top", legend.direction = "horizontal",
    legend.key.width = unit(2.0, "cm"),
    legend.title = element_text(size = 9), legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

figure2_file <- file.path(figure_dir, "Fig2.pptx")
ggsave(
  file.path(cross_stage_output_dir, "panelA_nonnegative_CLM_bootstrap.pdf"),
  p_nonnegative_clm, width = 7.0, height = 6.4, device = cairo_pdf
)
ggsave(
  file.path(cross_stage_output_dir, "panelA_nonnegative_CLM_bootstrap.png"),
  p_nonnegative_clm, width = 7.0, height = 6.4, dpi = 400, bg = "white"
)
save_editable_pptx(
  p_nonnegative_clm, "panelA_nonnegative_CLM_bootstrap.pptx", 7.0, 6.4
)

# ---------------------------------------------------------------------------
# D. Formal topGO weight01 and clusterProfiler-like KEGG enrichment
# ---------------------------------------------------------------------------

go_map <- go_long %>%
  filter(aspect == "biological_process") %>%
  transmute(HOG, GO = go_terms, term_name) %>%
  filter(!is.na(HOG), HOG != "", !is.na(GO), GO != "") %>%
  distinct()
gene2go <- split(go_map$GO, go_map$HOG)

# Accommodate the two column conventions used by saved mapping objects.
if (!"HOG" %in% names(kegg_map)) {
  kegg_term2gene <- kegg_map %>% transmute(term, gene = gene) %>% distinct()
} else {
  kegg_term2gene <- kegg_map %>% transmute(term, gene = HOG) %>% distinct()
}
kegg_term2name <- kegg_names %>%
  transmute(term = annotation_id, name = annotation_name) %>%
  distinct()

run_topgo <- function(pair_dat, comparison, set_name) {
  background <- intersect(pair_dat$HOG, names(gene2go))
  foreground <- switch(
    set_name,
    concordant_all = pair_dat$HOG[pair_dat$concordant_DE],
    concordant_up = pair_dat$HOG[pair_dat$direction == "concordant up"],
    concordant_down = pair_dat$HOG[pair_dat$direction == "concordant down"]
  ) %>% unique() %>% intersect(background)
  gene_list <- factor(as.integer(background %in% foreground), levels = c(0, 1))
  names(gene_list) <- background
  go_data <- new(
    "topGOdata", ontology = "BP", allGenes = gene_list,
    annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize = 20
  )
  result <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  tab <- GenTable(
    go_data, weight01 = result,
    orderBy = "weight01", topNodes = length(score(result))
  ) %>%
    as_tibble() %>%
    rename(
      GO = GO.ID, term_name = Term, background_n = Annotated,
      foreground_n = Significant, expected_n = Expected
    ) %>%
    mutate(
      weight01_p = parse_topgo_p(weight01),
      fold_enrichment = foreground_n / expected_n,
      comparison = comparison,
      set = set_name,
      foreground_total = length(foreground),
      background_total = length(background),
      foreground_HOGs = map_chr(
        GO,
        ~ paste(
          sort(intersect(genesInTerm(go_data, .x)[[1]], foreground)),
          collapse = "; "
        )
      )
    ) %>%
    arrange(weight01_p, desc(fold_enrichment))
  list(table = tab, data = go_data, foreground = foreground)
}

run_kegg <- function(pair_dat, comparison, set_name) {
  background <- intersect(pair_dat$HOG, unique(kegg_term2gene$gene))
  foreground <- switch(
    set_name,
    concordant_all = pair_dat$HOG[pair_dat$concordant_DE],
    concordant_up = pair_dat$HOG[pair_dat$direction == "concordant up"],
    concordant_down = pair_dat$HOG[pair_dat$direction == "concordant down"]
  ) %>% unique() %>% intersect(background)
  fit <- enricher(
    gene = foreground, universe = background,
    TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500
  )
  if (is.null(fit) || nrow(as.data.frame(fit)) == 0) return(tibble())
  as.data.frame(fit) %>%
    as_tibble() %>%
    transmute(
      comparison = comparison, set = set_name,
      pathway = ID, pathway_name = Description,
      foreground_ratio = GeneRatio, background_ratio = BgRatio,
      p_value = pvalue, p_BH = p.adjust, q_value = qvalue,
      foreground_HOGs = geneID, foreground_n = Count,
      foreground_total = length(foreground),
      background_total = length(background)
    ) %>%
    arrange(p_BH, p_value)
}

selected_pair_tables <- pmap(
  selected_pairs,
  function(comparison, pd_stage, vv_stage) {
    get_pair(pd_stage, vv_stage) %>% mutate(comparison = comparison)
  }
)
names(selected_pair_tables) <- selected_pairs$comparison

sets <- c("concordant_all", "concordant_up", "concordant_down")
topgo_objects <- list()
topgo_results <- list()
kegg_results <- list()
for (comparison in names(selected_pair_tables)) {
  for (set_name in sets) {
    key_name <- paste(comparison, set_name, sep = "__")
    topgo_objects[[key_name]] <- run_topgo(
      selected_pair_tables[[comparison]], comparison, set_name
    )
    topgo_results[[key_name]] <- topgo_objects[[key_name]]$table
    kegg_results[[key_name]] <- run_kegg(
      selected_pair_tables[[comparison]], comparison, set_name
    )
  }
}
topgo_results_all <- bind_rows(topgo_results)
kegg_results_all <- bind_rows(kegg_results)

theme_patterns <- tribble(
  ~functional_theme, ~pattern,
  "IIS/Wnt/mTOR/FoxO", "insulin|foxo|target of rapamycin|tor signaling|wnt|glycogen synthase kinase|nutrient sensing",
  "Ecdysone/development", "ecdys|ecdysone|juvenile hormone|nuclear receptor|steroid hormone receptor|molting|moulting|metamorph|instar|post-embryonic development",
  "Reproduction", "reproduct|oogen|gonad|ovary|ovarian|follic|vitell|germline|egg",
  "Lifespan/diapause/stress", "lifespan|aging|ageing|diapause|dormancy|starvation|hypoxia|oxidative stress|autophag",
  "Nutrient/lipid metabolism", "lipid|fatty acid|carbohydrate|glucose|glycogen|gluconeo|amino acid|energy reserve|nutrient"
)

classify_theme <- function(term) {
  out <- rep(NA_character_, length(term))
  for (i in seq_len(nrow(theme_patterns))) {
    hit <- is.na(out) & str_detect(
      str_to_lower(term), regex(theme_patterns$pattern[i], ignore_case = TRUE)
    )
    out[hit] <- theme_patterns$functional_theme[i]
  }
  out
}

eligible_topgo <- topgo_results_all %>%
  mutate(functional_theme = classify_theme(term_name)) %>%
  filter(
    weight01_p <= 0.05,
    foreground_n >= 2,
    fold_enrichment >= 1.5,
    !is.na(functional_theme)
  )

# Prefer direction-specific enrichments. If none exist for a comparison/theme,
# retain the strongest combined enrichment.
selected_go_terms <- eligible_topgo %>%
  mutate(set_priority = ifelse(set == "concordant_all", 2L, 1L)) %>%
  group_by(comparison, functional_theme) %>%
  arrange(set_priority, weight01_p, desc(fold_enrichment), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

term_hogs <- pmap_dfr(
  selected_go_terms,
  function(GO, term_name, background_n, foreground_n, expected_n,
           weight01, weight01_p, fold_enrichment, comparison, set,
           foreground_total, background_total, functional_theme, ...) {
    key_name <- paste(comparison, set, sep = "__")
    obj <- topgo_objects[[key_name]]
    hogs <- intersect(genesInTerm(obj$data, GO)[[1]], obj$foreground)
    tibble(
      comparison, set, functional_theme, GO, term_name,
      weight01_p, fold_enrichment, HOG = hogs
    )
  }
)

orth_annotation <- as_tibble(ortholog_input) %>%
  transmute(
    HOG = orthogroup, relationship_class,
    gene_pd, symbol_pd, description_pd,
    gene_vv, symbol_vv, description_vv,
    gene_dmel, symbol_dmel, description_dmel
  ) %>%
  distinct(HOG, .keep_all = TRUE)

first_symbol <- function(x) {
  x <- coalesce(x, "")
  str_trim(str_split_fixed(x, "\\s*\\|\\s*", 2)[, 1])
}

candidate_hogs <- term_hogs %>%
  inner_join(
    bind_rows(selected_pair_tables) %>%
      select(
        comparison, HOG, pd_stage, vv_stage,
        pd_log2FC, pd_SE, pd_padj,
        vv_log2FC, vv_SE, vv_padj,
        direction, concordance_z
      ),
    by = c("comparison", "HOG")
  ) %>%
  left_join(orth_annotation, by = "HOG") %>%
  mutate(
    plot_label = case_when(
      first_symbol(symbol_pd) != "" &
        !str_detect(first_symbol(symbol_pd), "^LOC[0-9]+$") ~ first_symbol(symbol_pd),
      first_symbol(symbol_dmel) != "" &
        !str_detect(first_symbol(symbol_dmel), "^CG[0-9]+$") ~ first_symbol(symbol_dmel),
      TRUE ~ HOG
    )
  ) %>%
  group_by(comparison, functional_theme) %>%
  arrange(desc(concordance_z), HOG, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

# A HOG can occur in several enriched GO terms. Assign each plotted HOG once,
# to the strongest supported term, so rows are not duplicated across facets.
automatic_selected_hogs <- candidate_hogs %>%
  arrange(weight01_p, desc(concordance_z), HOG) %>%
  group_by(HOG) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    HOG, plot_label, functional_theme,
    selected_from = comparison,
    selected_GO = GO,
    GO_term = term_name,
    topGO_weight01_p = weight01_p,
    max_concordance_z = concordance_z
  )

# Mechanistically interpretable representatives from the formally enriched GO
# or BH-significant KEGG sets. The Wang et al. column records direct matches to
# proteins shown in their IIS/FoxO/TOR schematic; it is not used as a test.
mechanistic_priority <- tribble(
  ~HOG, ~plot_label, ~functional_theme, ~selected_from, ~Wang_2021_component,
  "N13.HOG0005758", "Wdr59 (GATOR2)", "IIS/Wnt/mTOR/FoxO", "L1_to_L2", NA_character_,
  "N13.HOG0004507", "nej/CBP", "IIS/Wnt/mTOR/FoxO", "L1_to_L2", "CBP",
  "N13.HOG0008912", "sut1-3 / SLC2-like", "IIS/Wnt/mTOR/FoxO", "L4_to_L2", NA_character_,
  "N13.HOG0006185", "sgg/GSK3beta", "IIS/Wnt/mTOR/FoxO", "L4_to_L2; L5_to_L2", "GSK3beta",
  "N13.HOG0009263", "nmo/NLK", "IIS/Wnt/mTOR/FoxO", "L4_to_L2", NA_character_,
  "N13.HOG0010043", "Gbs-76A", "IIS/Wnt/mTOR/FoxO", "L4_to_L2", NA_character_,
  "N13.HOG0006961", "Pepck", "IIS/Wnt/mTOR/FoxO", "L4_to_L2", "PEPCK",
  "N13.HOG0007424", "Pdk1", "IIS/Wnt/mTOR/FoxO", "L4_to_L2", "PDK1",
  "N13.HOG0003585", "br", "Ecdysone/development", "L4_to_L2", NA_character_,
  "N13.HOG0004754", "EcR", "Ecdysone/development", "L4_to_L2", NA_character_,
  "N13.HOG0005004", "ss", "Reproduction", "L4_to_L2", NA_character_,
  "N13.HOG0006201", "Actbeta", "Reproduction", "L4_to_L2", NA_character_,
  "N13.HOG0000655", "bab1/2", "Reproduction", "L5_to_L2", NA_character_,
  "N13.HOG0006858", "Atg17", "Lifespan/diapause/stress", "L4_to_L2", NA_character_
)

go_hog_evidence <- eligible_topgo %>%
  separate_rows(foreground_HOGs, sep = ";\\s*") %>%
  rename(HOG = foreground_HOGs) %>%
  filter(HOG != "") %>%
  group_by(HOG) %>%
  summarise(
    enriched_GO_evidence = paste(
      unique(paste0(comparison, ": ", GO, " ", term_name)),
      collapse = "; "
    ),
    minimum_topGO_weight01_p = min(weight01_p),
    .groups = "drop"
  )

kegg_hog_evidence <- kegg_results_all %>%
  filter(p_BH < 0.05) %>%
  separate_rows(foreground_HOGs, sep = "/") %>%
  rename(HOG = foreground_HOGs) %>%
  filter(HOG != "") %>%
  group_by(HOG) %>%
  summarise(
    BH_significant_KEGG_evidence = paste(
      unique(paste0(comparison, ": ", pathway, " ", pathway_name)),
      collapse = "; "
    ),
    minimum_KEGG_BH_p = min(p_BH),
    .groups = "drop"
  )

concordance_summary <- bind_rows(selected_pair_tables) %>%
  filter(concordant_DE) %>%
  group_by(HOG) %>%
  summarise(max_concordance_z = max(concordance_z), .groups = "drop")

selected_hogs <- mechanistic_priority %>%
  left_join(go_hog_evidence, by = "HOG") %>%
  left_join(kegg_hog_evidence, by = "HOG") %>%
  left_join(concordance_summary, by = "HOG") %>%
  filter(!is.na(max_concordance_z))

heatmap_values_pd <- de_long %>%
  filter(stage %in% c("L1", "L4", "L5"), HOG %in% selected_hogs$HOG) %>%
  transmute(
    HOG, column = paste0("Polistes ", stage),
    log2FC = log2FC_ash_pd, padj_global = padj_global_pd
  )
heatmap_values_vv <- de_long %>%
  filter(stage == "L2", HOG %in% selected_hogs$HOG) %>%
  transmute(
    HOG, column = "Vespula L2",
    log2FC = log2FC_ash_vv, padj_global = padj_global_vv
  )

heatmap_data <- bind_rows(heatmap_values_pd, heatmap_values_vv) %>%
  left_join(selected_hogs, by = "HOG") %>%
  mutate(
    column = factor(
      column,
      levels = c("Polistes L1", "Polistes L4", "Polistes L5", "Vespula L2"),
      labels = c(
        "Polistes L1\nlate vs early", "Polistes L4\nlate vs early",
        "Polistes L5\nlate vs early", "Vespula L2\nqueen vs worker"
      )
    ),
    functional_theme = factor(
      functional_theme,
      levels = theme_patterns$functional_theme
    ),
    de_significance = case_when(
      padj_global < 0.001 ~ "***",
      padj_global < 0.01 ~ "**",
      padj_global < 0.05 ~ "*",
      TRUE ~ ""
    ),
    plot_label = factor(
      plot_label,
      levels = rev(unique(selected_hogs %>%
                            arrange(functional_theme, desc(max_concordance_z)) %>%
                            pull(plot_label)))
    )
  )

p_functional <- ggplot(heatmap_data, aes(column, plot_label, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = de_significance), family = "Arial", fontface = "bold",
    size = 3.35, colour = "black"
  ) +
  facet_grid(
    functional_theme ~ ., scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-3, 3), oob = squish
  ) +
  labs(x = NULL, y = NULL, fill = expression("ASH log"[2]*"FC")) +
  theme_minimal(base_family = "Arial", base_size = 8.2) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.1),
    axis.text.x = element_text(size = 8.7, colour = "black"),
    axis.text.y = element_text(size = 8.6, colour = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(0.9, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(4, 5, 4, 4)
  )

height_functional <- max(5.0, 1.8 + 0.28 * n_distinct(selected_hogs$HOG))
ggsave(file.path(cross_stage_output_dir, "panelB_GO_supported_functional_HOGs.pdf"),
       p_functional, width = 7.3, height = height_functional, device = cairo_pdf)
ggsave(file.path(cross_stage_output_dir, "panelB_GO_supported_functional_HOGs.png"),
       p_functional, width = 7.3, height = height_functional,
       dpi = 450, bg = "white")
save_editable_pptx(
  p_functional, "panelB_GO_supported_functional_HOGs.pptx",
  7.3, height_functional
)

# Main-text layout for the pairwise directional-concordance result and the
# functionally supported HOG heatmap.
figure2_ab_plot <- p_pairwise_concordance + p_functional +
  plot_layout(widths = c(0.95, 1.35)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 24))

# Final main-text Fig. 2 contains the pairwise directional-concordance matrix
# and the functionally supported HOG heatmap. The mechanistic pathway diagram
# is exported independently as Fig2C_source.* for manual redrawing and final
# assembly in PowerPoint.
ggsave(
  file.path(figure_dir, "Fig2.pdf"), figure2_ab_plot,
  width = 13.33, height = 7.5, device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "Fig2.png"), figure2_ab_plot,
  width = 13.33, height = 7.5, dpi = 450, bg = "white"
)
save_editable_pptx(
  figure2_ab_plot, figure2_file, 13.33, 7.5
)

# Supplemental Fig. S4 presents the two all-stage constrained models using
# a common 6 x 6 geometry and a single BH false-discovery-rate correction
# over all 36 coefficients in each panel.
figure_s4_plot <- p_nonnegative_clm + p_all_stage_ridge +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 24))

figure_s4_file <- file.path(supplemental_figure_dir, "FigS4.pptx")
ggsave(
  file.path(supplemental_figure_dir, "FigS4.pdf"), figure_s4_plot,
  width = 13.33, height = 7.5, device = cairo_pdf
)
ggsave(
  file.path(supplemental_figure_dir, "FigS4.png"), figure_s4_plot,
  width = 13.33, height = 7.5, dpi = 450, bg = "white"
)
save_editable_pptx(
  figure_s4_plot, figure_s4_file,
  13.33, 7.5
)

# ---------------------------------------------------------------------------
# E. Annotated workbooks and machine-readable outputs
# ---------------------------------------------------------------------------

go_annotation <- go_map %>%
  group_by(HOG) %>%
  summarise(
    GO_BP_IDs = paste(sort(unique(GO)), collapse = "; "),
    GO_BP_names = paste(sort(unique(term_name)), collapse = "; "),
    .groups = "drop"
  )
kegg_annotation <- kegg_term2gene %>%
  left_join(kegg_term2name, by = "term") %>%
  group_by(HOG = gene) %>%
  summarise(
    KEGG_IDs = paste(sort(unique(term)), collapse = "; "),
    KEGG_names = paste(sort(unique(name)), collapse = "; "),
    .groups = "drop"
  )

annotated_concordant <- imap(selected_pair_tables, function(dat, comparison) {
  dat %>%
    filter(concordant_DE) %>%
    left_join(orth_annotation, by = "HOG") %>%
    left_join(go_annotation, by = "HOG") %>%
    left_join(kegg_annotation, by = "HOG") %>%
    arrange(desc(concordance_z), HOG)
})

wb <- createWorkbook()
workbook_tables <- c(
  annotated_concordant,
  list(
    Pairwise_Fisher = pairwise_concordance_results %>%
      mutate(across(where(is.factor), as.character)),
    Nonnegative_CLM = nonnegative_clm_results %>%
      mutate(across(where(is.factor), as.character)),
    All_stage_ridge = all_stage_ridge_results %>%
      mutate(across(where(is.factor), as.character)),
    topGO_weight01 = topgo_results_all,
    topGO_shortlist = selected_go_terms,
    KEGG_clusterProfiler = kegg_results_all,
    PanelB_auto_selection = automatic_selected_hogs,
    PanelB_HOG_selection = selected_hogs,
    PanelB_plot_data = heatmap_data %>%
      mutate(across(where(is.factor), as.character))
  )
)
for (nm in names(workbook_tables)) {
  sheet <- substr(nm, 1, 31)
  addWorksheet(wb, sheet)
  writeDataTable(wb, sheet, workbook_tables[[nm]], tableStyle = "TableStyleMedium2")
  freezePane(wb, sheet, firstRow = TRUE)
}
concordant_workbook <- file.path(
  cross_stage_output_dir, "concordant_DE_HOGs_annotated.xlsx"
)
saveWorkbook(wb, concordant_workbook, overwrite = TRUE)
rebuild_openxlsx_compat(concordant_workbook)
file.copy(
  concordant_workbook,
  file.path(cross_stage_output_dir, "cross_stage_concordance_results.xlsx"),
  overwrite = TRUE
)

write_cross_stage_tsv(
  pairwise_concordance_results %>%
    mutate(
      across(where(is.factor), as.character),
      cell_label = str_replace_all(cell_label, "\\n", " ")
    ),
  "pairwise_concordance_one_sided_Fisher.tsv"
)
write_cross_stage_tsv(
  all_stage_ridge_results %>% mutate(across(where(is.factor), as.character)),
  "all_stage_nonnegative_ridge.tsv"
)
write_cross_stage_tsv(
  nonnegative_clm_results %>%
    mutate(
      across(where(is.factor), as.character),
      panel_label = str_replace_all(panel_label, "\\n", " ")
    ),
  "nonnegative_CLM_bootstrap_all_Polistes_stages.tsv"
)
write_cross_stage_tsv(topgo_results_all, "topGO_weight01_all_results.tsv")
write_cross_stage_tsv(selected_go_terms, "topGO_weight01_functional_shortlist.tsv")
write_cross_stage_tsv(kegg_results_all, "KEGG_clusterProfiler_enrichment.tsv")
write_cross_stage_tsv(selected_hogs, "panelB_GO_supported_HOG_selection.tsv")

save(
  pairwise_concordance_results, nonnegative_clm_results,
  all_stage_ridge_results,
  topgo_results_all, selected_go_terms, kegg_results_all,
  automatic_selected_hogs, selected_hogs, heatmap_data,
  file = file.path(cross_stage_output_dir, "cross_stage_concordance_workspace.RData")
)

writeLines(
  c(
    paste("Completed", Sys.time()),
    paste("Pairwise Fisher stage cells", nrow(pairwise_concordance_results)),
    paste("Nonnegative CLM coefficients", nrow(nonnegative_clm_results)),
    paste("Nonnegative CLM converged stage fits",
          sum(nonnegative_clm_results$model_convergence == 0) / length(stages)),
    paste("Bootstrap resamples per constrained model", bootstrap_B),
    paste("topGO result rows", nrow(topgo_results_all)),
    paste("KEGG result rows", nrow(kegg_results_all)),
    paste("Selected GO terms", nrow(selected_go_terms)),
    paste("Selected panel-B HOGs", nrow(selected_hogs)),
    R.version.string
  ),
  file.path(cross_stage_output_dir, "VALIDATION.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(cross_stage_output_dir, "sessionInfo.txt")
)
