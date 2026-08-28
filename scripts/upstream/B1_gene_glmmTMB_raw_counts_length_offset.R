# GLMMTMB DIFFERENTIAL EXPRESSION: ORIGINAL SALMON COUNTS + TXIMPORT OFFSET ####
#
# Separate alternative to B1_analysis_differential_expression_glmmSeq_glmmTMB_version.R.
# It fits the same stage-specific negative-binomial mixed model, but uses
# original Salmon estimated counts and a gene-by-sample offset that combines
# sample-specific average transcript length with library-size normalisation.

if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  stop("Set PROJECT_ROOT to the repository root before sourcing this script.")
}
dir_base <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
input_root <- file.path(dir_base, "input_differential_expression", "raw")
output_root <- file.path(dir_base, "input_differential_expression", "checkpoints")
output_dir <- file.path(output_root, "gene")
local_library <- file.path(dir_base, ".R_library")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (dir.exists(local_library)) {
  Sys.setenv(R_LIBS_USER = normalizePath(local_library, winslash = "/"))
  .libPaths(c(normalizePath(local_library, winslash = "/"), .libPaths()))
}

required_packages <- c("glmmTMB", "DESeq2", "readr", "dplyr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(glmmTMB)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

options(contrasts = c("contr.treatment", "contr.poly"))
stages <- c("L1", "L2", "L3", "L4", "L5", "P")
global_BH <- TRUE
n_cores <- suppressWarnings(as.integer(Sys.getenv("B1_N_CORES", "4")))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

read_gene_matrix <- function(path) {
  tab <- readr::read_tsv(path, show_col_types = FALSE)
  required <- c("gene_id", "gene_name")
  if (!all(required %in% names(tab))) {
    stop("Expected gene_id and gene_name columns in: ", path)
  }
  if (anyDuplicated(tab$gene_id)) stop("Duplicated gene IDs in: ", path)
  gene_name <- tab$gene_name
  mat <- tab |>
    dplyr::select(-gene_name) |>
    tibble::column_to_rownames("gene_id") |>
    as.matrix()
  storage.mode(mat) <- "numeric"
  colnames(mat) <- tolower(colnames(mat))
  list(matrix = mat, gene_name = stats::setNames(gene_name, rownames(mat)))
}

metadata_path <- function(species_code) {
  path <- file.path(input_root, species_code, paste0("metadata_", species_code, ".csv"))
  if (!file.exists(path)) stop("Metadata file not found for ", species_code, ": ", path)
  path
}

read_species_data <- function(species_code, low_level, high_level) {
  input_dir <- file.path(input_root, species_code)
  counts_obj <- read_gene_matrix(file.path(input_dir, "salmon.merged.gene_counts.tsv"))
  lengths_obj <- read_gene_matrix(file.path(input_dir, "salmon.merged.gene_lengths.tsv"))
  tpm_obj <- read_gene_matrix(file.path(input_dir, "salmon.merged.gene_tpm.tsv"))

  counts <- counts_obj$matrix
  lengths <- lengths_obj$matrix
  tpm <- tpm_obj$matrix

  if (!setequal(rownames(counts), rownames(lengths)) ||
      !setequal(rownames(counts), rownames(tpm))) {
    stop("Count, length and TPM matrices have different gene sets for ", species_code)
  }
  if (!setequal(colnames(counts), colnames(lengths)) ||
      !setequal(colnames(counts), colnames(tpm))) {
    stop("Count, length and TPM matrices have different sample sets for ", species_code)
  }
  lengths <- lengths[rownames(counts), colnames(counts), drop = FALSE]
  tpm <- tpm[rownames(counts), colnames(counts), drop = FALSE]

  metadata <- readr::read_csv(metadata_path(species_code), show_col_types = FALSE) |>
    dplyr::mutate(
      sampleNumber = tolower(sampleNumber),
      larvalStage = factor(larvalStage, levels = stages),
      lifeStage = factor(lifeStage, levels = c(low_level, high_level)),
      colonyNumber = factor(colonyNumber),
      colonyNested = interaction(lifeStage, colonyNumber, drop = TRUE)
    ) |>
    dplyr::distinct(sampleNumber, .keep_all = TRUE) |>
    dplyr::filter(sampleNumber %in% colnames(counts)) |>
    dplyr::arrange(match(sampleNumber, colnames(counts))) |>
    as.data.frame()

  stopifnot(
    nrow(metadata) == ncol(counts),
    identical(metadata$sampleNumber, colnames(counts)),
    !anyNA(metadata$larvalStage),
    !anyNA(metadata$lifeStage),
    !anyNA(metadata$colonyNested),
    all(is.finite(counts)),
    all(counts >= 0),
    all(is.finite(tpm)),
    all(tpm >= 0)
  )
  rownames(metadata) <- metadata$sampleNumber

  valid_length <- apply(lengths, 1, function(x) all(is.finite(x) & x > 0))
  if (!all(valid_length)) {
    warning(sum(!valid_length), " genes have invalid lengths and will be excluded for ", species_code)
  }

  list(
    counts = counts[valid_length, , drop = FALSE],
    lengths = lengths[valid_length, , drop = FALSE],
    tpm = tpm[valid_length, , drop = FALSE],
    gene_name = counts_obj$gene_name[rownames(counts)[valid_length]],
    metadata = metadata,
    n_invalid_length_genes = sum(!valid_length),
    low_level = low_level,
    high_level = high_level
  )
}

make_offset_data <- function(input) {
  txi <- list(
    abundance = input$tpm,
    counts = input$counts,
    length = input$lengths,
    countsFromAbundance = "no"
  )
  dds <- DESeq2::DESeqDataSetFromTximport(
    txi,
    colData = input$metadata,
    design = ~ 1
  )
  rounded_counts <- DESeq2::counts(dds)
  keep <- rowSums(rounded_counts >= 10) >= 5
  dds <- dds[keep, ]

  max_genes <- suppressWarnings(as.integer(Sys.getenv("B1_MAX_GENES", "0")))
  if (!is.na(max_genes) && max_genes > 0L && nrow(dds) > max_genes) {
    dds <- dds[seq_len(max_genes), ]
  }

  dds <- DESeq2::estimateSizeFactors(dds)
  normalisation_factors <- DESeq2::normalizationFactors(dds)
  if (is.null(normalisation_factors)) {
    stop("DESeq2 did not construct gene-by-sample normalisation factors")
  }
  offset_matrix <- log(normalisation_factors)
  if (any(!is.finite(offset_matrix))) stop("Non-finite values in offset matrix")

  list(
    dds = dds,
    counts = DESeq2::counts(dds),
    normalisation_factors = normalisation_factors,
    offset = offset_matrix,
    base_mean = rowMeans(DESeq2::counts(dds, normalized = TRUE)),
    keep = keep
  )
}

fit_one_gene <- function(i, counts, offsets, metadata, model_formula, control) {
  dat <- metadata
  dat$count <- as.numeric(counts[i, ])
  dat$offset_gs <- as.numeric(offsets[i, ])
  warning_messages <- character()

  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB::glmmTMB(
        model_formula,
        family = glmmTMB::nbinom2(link = "log"),
        data = dat,
        control = control
      ),
      warning = function(w) {
        warning_messages <<- unique(c(warning_messages, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(success = FALSE, error = conditionMessage(fit), warnings = paste(warning_messages, collapse = " | ")))
  }

  coef_table <- tryCatch(coef(summary(fit))$cond, error = function(e) NULL)
  fixed <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
  random_variance <- tryCatch(
    as.numeric(glmmTMB::VarCorr(fit)$cond[[1]][1, 1]),
    error = function(e) NA_real_
  )
  max_gradient <- tryCatch(max(abs(fit$fit$gradient), na.rm = TRUE), error = function(e) NA_real_)

  list(
    success = !is.null(coef_table) && !is.null(fixed),
    coefficients = fixed,
    standard_errors = if (is.null(coef_table)) NULL else coef_table[, "Std. Error"],
    p_values = if (is.null(coef_table)) NULL else coef_table[, "Pr(>|z|)"],
    dispersion = tryCatch(glmmTMB::sigma(fit), error = function(e) NA_real_),
    aic = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
    log_likelihood = tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_),
    random_intercept_variance = random_variance,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess),
    convergence_code = fit$fit$convergence,
    max_absolute_gradient = max_gradient,
    message = fit$fit$message,
    warnings = paste(warning_messages, collapse = " | "),
    error = ""
  )
}

fit_all_genes <- function(counts, offsets, metadata, model_formula, control, n_cores) {
  indices <- seq_len(nrow(counts))
  if (n_cores <= 1L) {
    return(lapply(indices, fit_one_gene, counts, offsets, metadata, model_formula, control))
  }

  cluster <- parallel::makePSOCKcluster(n_cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, suppressPackageStartupMessages(library(glmmTMB)))
  parallel::clusterExport(
    cluster,
    c("fit_one_gene", "counts", "offsets", "metadata", "model_formula", "control"),
    envir = environment()
  )
  parallel::parLapplyLB(cluster, indices, function(i) {
    fit_one_gene(i, counts, offsets, metadata, model_formula, control)
  })
}

run_species <- function(species_code, species_name, low_level, high_level, contrast_label) {
  message("Preparing original-count offsets for ", species_name, "...")
  input <- read_species_data(species_code, low_level, high_level)
  prepared <- make_offset_data(input)
  counts <- prepared$counts
  offsets <- prepared$offset
  expected_coefficients <- paste0("larvalStage", stages, ":lifeStage", high_level)

  model_formula <- count ~ 0 + larvalStage + larvalStage:lifeStage +
    offset(offset_gs) + (1 | colonyNested)
  control <- glmmTMB::glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000),
    rank_check = "warning"
  )

  message("Fitting ", nrow(counts), " genes for ", species_name, " using ", n_cores, " core(s)...")
  start_time <- Sys.time()
  fit_list <- fit_all_genes(
    counts, offsets, input$metadata, model_formula, control, n_cores
  )
  elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  names(fit_list) <- rownames(counts)

  all_coefficient_names <- unique(unlist(lapply(fit_list, function(x) names(x$coefficients))))
  coefficient_matrix <- matrix(
    NA_real_, nrow = nrow(counts), ncol = length(all_coefficient_names),
    dimnames = list(rownames(counts), all_coefficient_names)
  )
  se_matrix <- coefficient_matrix
  for (gene in rownames(counts)) {
    item <- fit_list[[gene]]
    if (isTRUE(item$success)) {
      coefficient_matrix[gene, names(item$coefficients)] <- item$coefficients
      se_matrix[gene, names(item$standard_errors)] <- item$standard_errors
    }
  }

  fit_diagnostics <- tibble::tibble(
    gene = rownames(counts),
    fit_success = vapply(fit_list, function(x) isTRUE(x$success), logical(1)),
    positive_definite_hessian = vapply(fit_list, function(x) isTRUE(x$positive_definite_hessian), logical(1)),
    convergence_code = vapply(fit_list, function(x) if (is.null(x$convergence_code)) NA_integer_ else as.integer(x$convergence_code), integer(1)),
    random_intercept_variance = vapply(fit_list, function(x) if (is.null(x$random_intercept_variance)) NA_real_ else x$random_intercept_variance, numeric(1)),
    dispersion = vapply(fit_list, function(x) if (is.null(x$dispersion)) NA_real_ else x$dispersion, numeric(1)),
    aic = vapply(fit_list, function(x) if (is.null(x$aic)) NA_real_ else x$aic, numeric(1)),
    max_absolute_gradient = vapply(fit_list, function(x) if (is.null(x$max_absolute_gradient)) NA_real_ else x$max_absolute_gradient, numeric(1)),
    warnings = vapply(fit_list, function(x) if (is.null(x$warnings)) "" else x$warnings, character(1)),
    error = vapply(fit_list, function(x) if (is.null(x$error)) "" else x$error, character(1))
  )

  finite_targets <- apply(
    coefficient_matrix[, expected_coefficients, drop = FALSE], 1,
    function(x) all(is.finite(x))
  ) & apply(
    se_matrix[, expected_coefficients, drop = FALSE], 1,
    function(x) all(is.finite(x) & x > 0)
  )
  valid_fit <- fit_diagnostics$fit_success &
    fit_diagnostics$positive_definite_hessian & finite_targets
  valid_fit[is.na(valid_fit)] <- FALSE
  fit_diagnostics$retained_for_inference <- valid_fit
  valid_genes <- fit_diagnostics$gene[valid_fit]

  results <- dplyr::bind_rows(lapply(stages, function(stage) {
    coefficient <- paste0("larvalStage", stage, ":lifeStage", high_level)
    estimate_loge <- coefficient_matrix[valid_genes, coefficient]
    se_loge <- se_matrix[valid_genes, coefficient]
    z_value <- estimate_loge / se_loge
    p_value <- 2 * stats::pnorm(-abs(z_value))
    tibble::tibble(
      stage = stage,
      gene = valid_genes,
      gene_name = unname(input$gene_name[valid_genes]),
      baseMean = unname(prepared$base_mean[valid_genes]),
      log2FC_raw = estimate_loge / log(2),
      lfcSE_raw = se_loge / log(2),
      pvalue = p_value,
      padj_stage_wise = p.adjust(p_value, method = "BH"),
      contrast = contrast_label,
      log2FC_shrunk = estimate_loge / log(2),
      lfcSE_shrunk = se_loge / log(2),
      z = z_value,
      method = "glmmTMB_raw_counts_gene_offset",
      lfc_shrinkage = "none"
    )
  }))
  results$padj_global <- p.adjust(results$pvalue, method = "BH")
  results$padj <- if (global_BH) results$padj_global else results$padj_stage_wise

  readr::write_tsv(
    results,
    file.path(output_dir, paste0("glmmTMB_raw_counts_gene_offset_", species_code, "_results.tsv")),
    na = "NA"
  )
  readr::write_tsv(
    fit_diagnostics,
    file.path(output_dir, paste0("glmmTMB_raw_counts_gene_offset_", species_code, "_fit_diagnostics.tsv")),
    na = "NA"
  )

  sample_offset_summary <- tibble::tibble(
    sample = colnames(counts),
    raw_library_size = colSums(counts),
    size_factor = DESeq2::sizeFactors(prepared$dds),
    geometric_mean_normalisation_factor = exp(colMeans(log(prepared$normalisation_factors))),
    minimum_normalisation_factor = apply(prepared$normalisation_factors, 2, min),
    maximum_normalisation_factor = apply(prepared$normalisation_factors, 2, max)
  )
  readr::write_tsv(
    sample_offset_summary,
    file.path(output_dir, paste0("glmmTMB_raw_counts_gene_offset_", species_code, "_sample_offsets.tsv"))
  )

  run_diagnostics <- tibble::tibble(
    species = species_name,
    method = "direct glmmTMB; glmmSeq-compatible model",
    model = paste(deparse(model_formula), collapse = " "),
    distribution = "negative binomial 2",
    test = "coefficient-wise asymptotic Wald z",
    normalisation = "original Salmon estimated counts plus DESeq2/tximport gene-by-sample average-transcript-length and library-size offset",
    n_samples = ncol(counts),
    n_colonies_nested_in_condition = nlevels(input$metadata$colonyNested),
    n_genes_input = nrow(input$counts),
    n_genes_invalid_length = input$n_invalid_length_genes,
    n_genes_after_prefilter = nrow(counts),
    n_genes_fitted_and_valid = length(valid_genes),
    n_fit_failures_or_invalid = nrow(counts) - length(valid_genes),
    n_near_zero_random_variance = sum(fit_diagnostics$random_intercept_variance < 1e-8, na.rm = TRUE),
    n_nonpositive_hessian = sum(!fit_diagnostics$positive_definite_hessian, na.rm = TRUE),
    n_stage_specific_tests = nrow(results),
    n_global_BH_0_05 = sum(results$padj_global < 0.05, na.rm = TRUE),
    elapsed_minutes = elapsed_minutes
  )
  readr::write_tsv(
    run_diagnostics,
    file.path(output_dir, paste0("glmmTMB_raw_counts_gene_offset_", species_code, "_run_diagnostics.tsv"))
  )

  saveRDS(
    list(
      species_code = species_code,
      species_name = species_name,
      model_formula = model_formula,
      coefficients = coefficient_matrix,
      standard_errors = se_matrix,
      fit_diagnostics = fit_diagnostics,
      results = results,
      counts = counts,
      normalisation_factors = prepared$normalisation_factors,
      offsets = prepared$offset,
      metadata = input$metadata,
      run_diagnostics = run_diagnostics
    ),
    file.path(output_dir, paste0("glmmTMB_raw_counts_gene_offset_", species_code, "_analysis_objects.rds")),
    compress = "xz"
  )

  results
}

species_to_run <- tolower(Sys.getenv("B1_SPECIES", "both"))
if (!species_to_run %in% c("both", "pd", "vv")) {
  stop("B1_SPECIES must be one of: both, pd, vv")
}
if (species_to_run %in% c("both", "pd")) {
  raw_offset_pd <- run_species("pd", "Polistes dominula", "E", "L", "L vs E")
}
if (species_to_run %in% c("both", "vv")) {
  raw_offset_vv <- run_species("vv", "Vespula vulgaris", "W", "Q", "Q vs W")
}

analysis_metadata <- tibble::tibble(
  key = c(
    "script", "analysis_date", "global_BH", "model", "random_effect",
    "distribution", "normalisation", "response", "test", "cores",
    "species_run", "DESeq2_version", "glmmTMB_version", "glmmSeq_note", "R_version"
  ),
  value = c(
    "B1_analysis_differential_expression_glmmTMB_raw_counts_gene_offset.R",
    as.character(Sys.Date()), as.character(global_BH),
    "count ~ 0 + larvalStage + larvalStage:lifeStage + offset(offset_gs) + (1 | colonyNested)",
    "colony nested within lifeStage/condition",
    "glmmTMB nbinom2", "DESeq2/tximport gene-by-sample normalisation factors",
    "rounded original Salmon estimated gene counts",
    "coefficient-wise asymptotic Wald z", as.character(n_cores), species_to_run,
    as.character(utils::packageVersion("DESeq2")),
    as.character(utils::packageVersion("glmmTMB")),
    "glmmSeq 0.5.7 was bypassed because it accepts only a sample-level sizeFactors vector, not a gene-by-sample offset matrix; the same glmmTMB model engine and contrasts were fitted directly",
    R.version.string
  )
)
metadata_filename <- if (species_to_run == "both") {
  "analysis_metadata.tsv"
} else {
  paste0("analysis_metadata_", species_to_run, ".tsv")
}
readr::write_tsv(analysis_metadata, file.path(output_dir, metadata_filename))
message("Raw-count gene-offset glmmTMB outputs written to: ", output_dir)
