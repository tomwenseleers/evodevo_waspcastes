# N13 HOG-LEVEL GLMMTMB DIFFERENTIAL EXPRESSION ####
#
# Gene-level Salmon counts, TPMs and effective lengths are summarised to the
# OrthoFinder N13 hierarchical orthogroup level with tximport::summarizeToGene.
# Each HOG is then fitted with the same negative-binomial mixed model used for
# the gene-level analysis, including a HOG-by-sample length/library offset and
# a colony-within-condition random intercept.

if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  stop("Set PROJECT_ROOT to the repository root before sourcing this script.")
}
dir_base <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
de_input_dir <- file.path(dir_base, "input_differential_expression")
metadata_dir <- file.path(de_input_dir, "raw")
rnaseq_results_dir <- file.path(
  dir_base, "nextflow_runs", "1_rnaseq", "output_files"
)
orthofinder_results_dir <- file.path(
  dir_base, "nextflow_runs", "2_EXCON",
  "2_EXCON_orthofinder_eggnogmapper_run", "results_EXCON", "orthofinder"
)
output_root <- file.path(de_input_dir, "checkpoints")
output_dir <- file.path(output_root, "HOG")
local_library <- file.path(dir_base, ".R_library")
hog_path <- file.path(
  orthofinder_results_dir, "Phylogenetic_Hierarchical_Orthogroups", "N13.tsv"
)
gene_result_dir <- file.path(output_root, "gene")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (dir.exists(local_library)) {
  Sys.setenv(R_LIBS_USER = normalizePath(local_library, winslash = "/"))
  .libPaths(c(normalizePath(local_library, winslash = "/"), .libPaths()))
}

required_packages <- c(
  "glmmTMB", "DESeq2", "tximport", "ashr", "readr", "dplyr", "tidyr",
  "purrr", "tibble", "openxlsx"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(glmmTMB)
  library(DESeq2)
  library(tximport)
  library(ashr)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

options(contrasts = c("contr.treatment", "contr.poly"))
stages <- c("L1", "L2", "L3", "L4", "L5", "P")
global_BH <- TRUE
n_cores <- suppressWarnings(as.integer(Sys.getenv("B1_N_CORES", "1")))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
if (n_cores > 1L) {
  warning("For memory safety, this workflow is intended to use B1_N_CORES=1 per species.")
}

species_info <- tibble::tribble(
  ~species_code, ~species_name,        ~hog_column,                   ~low_level, ~high_level, ~contrast,
  "pd",          "Polistes dominula",  "Polistes_dominula.clean",     "E",        "L",         "L vs E",
  "vv",          "Vespula vulgaris",   "Vespula_vulgaris.clean",      "W",        "Q",         "Q vs W"
)

collapse_ids <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  paste(x, collapse = ", ")
}

split_ids <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(character())
  unique(trimws(unlist(strsplit(x, "\\s*,\\s*"))))
}

read_n13_mapping <- function() {
  if (!file.exists(hog_path)) stop("N13 HOG table not found: ", hog_path)
  hog <- readr::read_tsv(hog_path, show_col_types = FALSE)
  required <- c("HOG", "OG", species_info$hog_column)
  if (!all(required %in% names(hog))) {
    stop("N13 table lacks required columns: ", paste(setdiff(required, names(hog)), collapse = ", "))
  }
  if (anyDuplicated(hog$HOG)) stop("Duplicated HOG identifiers in N13 table")

  membership <- purrr::map2_dfr(
    species_info$species_code,
    species_info$hog_column,
    function(code, column) {
      hog |>
        dplyr::select(HOG, OG, member_string = dplyr::all_of(column)) |>
        dplyr::mutate(gene_id = purrr::map(member_string, split_ids)) |>
        tidyr::unnest_longer(gene_id, values_to = "gene_id", keep_empty = FALSE) |>
        dplyr::transmute(HOG, OG, species_code = code, gene_id)
    }
  ) |>
    dplyr::left_join(species_info |> dplyr::select(species_code, species_name), by = "species_code") |>
    dplyr::select(HOG, OG, species_code, species_name, gene_id) |>
    dplyr::arrange(HOG, species_code, gene_id)

  if (anyDuplicated(membership[c("species_code", "gene_id")])) {
    stop("At least one gene maps to more than one N13 HOG")
  }

  wide <- membership |>
    dplyr::group_by(HOG, OG, species_code) |>
    dplyr::summarise(
      genomic_member_count = dplyr::n_distinct(gene_id),
      genomic_member_ids = collapse_ids(gene_id),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = species_code,
      values_from = c(genomic_member_count, genomic_member_ids),
      values_fill = list(genomic_member_count = 0L, genomic_member_ids = "")
    ) |>
    dplyr::right_join(hog |> dplyr::select(HOG, OG), by = c("HOG", "OG")) |>
    dplyr::mutate(
      genomic_member_count_pd = tidyr::replace_na(genomic_member_count_pd, 0L),
      genomic_member_count_vv = tidyr::replace_na(genomic_member_count_vv, 0L),
      genomic_member_ids_pd = tidyr::replace_na(genomic_member_ids_pd, ""),
      genomic_member_ids_vv = tidyr::replace_na(genomic_member_ids_vv, ""),
      present_in_both_species = genomic_member_count_pd > 0L & genomic_member_count_vv > 0L,
      relationship_class = dplyr::case_when(
        !present_in_both_species ~ "not_shared",
        genomic_member_count_pd == 1L & genomic_member_count_vv == 1L ~ "one_to_one",
        genomic_member_count_pd == 1L & genomic_member_count_vv > 1L ~ "one_to_many",
        genomic_member_count_pd > 1L & genomic_member_count_vv == 1L ~ "many_to_one",
        TRUE ~ "many_to_many"
      )
    ) |>
    dplyr::arrange(HOG)

  list(hog = hog, membership = membership, wide = wide)
}

write_mapping_outputs <- function(mapping) {
  readr::write_tsv(
    mapping$membership,
    file.path(output_dir, "N13_HOG_membership_long.tsv"),
    na = "NA"
  )
  readr::write_tsv(
    mapping$wide,
    file.path(output_dir, "N13_HOG_membership_summary.tsv"),
    na = "NA"
  )
  summary <- tibble::tibble(
    metric = c(
      "N13_HOGs_total", "N13_HOGs_shared", "one_to_one", "one_to_many",
      "many_to_one", "many_to_many", "Polistes_genes", "Vespula_genes"
    ),
    value = c(
      nrow(mapping$wide),
      sum(mapping$wide$present_in_both_species),
      sum(mapping$wide$relationship_class == "one_to_one"),
      sum(mapping$wide$relationship_class == "one_to_many"),
      sum(mapping$wide$relationship_class == "many_to_one"),
      sum(mapping$wide$relationship_class == "many_to_many"),
      sum(mapping$wide$genomic_member_count_pd),
      sum(mapping$wide$genomic_member_count_vv)
    )
  )
  readr::write_tsv(summary, file.path(output_dir, "N13_HOG_mapping_summary.tsv"))
}

read_gene_matrix <- function(path) {
  tab <- readr::read_tsv(path, show_col_types = FALSE)
  required <- c("gene_id", "gene_name")
  if (!all(required %in% names(tab))) {
    stop("Expected gene_id and gene_name columns in: ", path)
  }
  if (anyDuplicated(tab$gene_id)) stop("Duplicated gene IDs in: ", path)
  gene_name <- stats::setNames(tab$gene_name, tab$gene_id)
  mat <- tab |>
    dplyr::select(-gene_name) |>
    tibble::column_to_rownames("gene_id") |>
    as.matrix()
  storage.mode(mat) <- "numeric"
  colnames(mat) <- tolower(colnames(mat))
  list(matrix = mat, gene_name = gene_name)
}

metadata_path <- function(species_code) {
  path <- file.path(metadata_dir, paste0("metadata_", species_code, ".csv"))
  if (!file.exists(path)) stop("Metadata file not found for ", species_code, ": ", path)
  path
}

read_species_data <- function(info) {
  result_folder <- switch(
    info$species_code,
    pd = "Results_Polistes",
    vv = "Results_Vespula",
    stop("Unknown species code: ", info$species_code)
  )
  input_dir <- file.path(rnaseq_results_dir, result_folder, "star_salmon")
  counts_obj <- read_gene_matrix(file.path(input_dir, "salmon.merged.gene_counts.tsv"))
  lengths_obj <- read_gene_matrix(file.path(input_dir, "salmon.merged.gene_lengths.tsv"))
  tpm_obj <- read_gene_matrix(file.path(input_dir, "salmon.merged.gene_tpm.tsv"))
  counts <- counts_obj$matrix
  lengths <- lengths_obj$matrix[rownames(counts), colnames(counts), drop = FALSE]
  tpm <- tpm_obj$matrix[rownames(counts), colnames(counts), drop = FALSE]

  if (!setequal(rownames(counts), rownames(lengths)) ||
      !setequal(rownames(counts), rownames(tpm)) ||
      !setequal(colnames(counts), colnames(lengths)) ||
      !setequal(colnames(counts), colnames(tpm))) {
    stop("Count, length and TPM matrices do not align for ", info$species_code)
  }

  metadata <- readr::read_csv(metadata_path(info$species_code), show_col_types = FALSE) |>
    dplyr::mutate(
      sampleNumber = tolower(sampleNumber),
      larvalStage = factor(larvalStage, levels = stages),
      lifeStage = factor(lifeStage, levels = c(info$low_level, info$high_level)),
      colonyNumber = factor(colonyNumber),
      colonyNested = interaction(lifeStage, colonyNumber, drop = TRUE)
    ) |>
    dplyr::distinct(sampleNumber, .keep_all = TRUE) |>
    dplyr::filter(sampleNumber %in% colnames(counts)) |>
    dplyr::arrange(match(sampleNumber, colnames(counts))) |>
    as.data.frame()
  rownames(metadata) <- metadata$sampleNumber

  stopifnot(
    nrow(metadata) == ncol(counts),
    identical(metadata$sampleNumber, colnames(counts)),
    all(is.finite(counts) & counts >= 0),
    all(is.finite(tpm) & tpm >= 0)
  )

  valid_length <- apply(lengths, 1, function(x) all(is.finite(x) & x > 0))
  list(
    counts = counts[valid_length, , drop = FALSE],
    lengths = lengths[valid_length, , drop = FALSE],
    tpm = tpm[valid_length, , drop = FALSE],
    gene_name = counts_obj$gene_name[rownames(counts)[valid_length]],
    metadata = metadata,
    n_invalid_length_genes = sum(!valid_length)
  )
}

summarise_to_hog <- function(input, membership, species_code) {
  map <- membership |>
    dplyr::filter(species_code == .env$species_code) |>
    dplyr::filter(gene_id %in% rownames(input$counts)) |>
    dplyr::select(TXNAME = gene_id, GENEID = HOG)
  if (!nrow(map)) stop("No N13 member genes found for ", species_code)
  if (anyDuplicated(map$TXNAME)) stop("A gene maps to multiple HOGs for ", species_code)

  genes <- rownames(input$counts)[rownames(input$counts) %in% map$TXNAME]
  txi_gene <- list(
    abundance = input$tpm[genes, , drop = FALSE],
    counts = input$counts[genes, , drop = FALSE],
    length = input$lengths[genes, , drop = FALSE],
    countsFromAbundance = "no"
  )
  txi_hog <- tximport::summarizeToGene(
    txi_gene,
    tx2gene = map,
    countsFromAbundance = "no"
  )
  stopifnot(
    identical(rownames(txi_hog$counts), rownames(txi_hog$abundance)),
    identical(rownames(txi_hog$counts), rownames(txi_hog$length)),
    all(is.finite(txi_hog$length) & txi_hog$length > 0),
    isTRUE(all.equal(
      unname(colSums(txi_hog$counts)),
      unname(colSums(input$counts[genes, , drop = FALSE])),
      tolerance = 1e-8
    )),
    isTRUE(all.equal(
      unname(colSums(txi_hog$abundance)),
      unname(colSums(input$tpm[genes, , drop = FALSE])),
      tolerance = 1e-8
    ))
  )
  txi_hog
}

write_hog_matrix <- function(mat, species_code, suffix) {
  out <- tibble::as_tibble(mat, rownames = "HOG")
  readr::write_tsv(
    out,
    file.path(output_dir, paste0("N13_HOG_", species_code, "_", suffix, ".tsv")),
    na = "NA"
  )
}

make_offset_data <- function(txi_hog, metadata) {
  dds <- DESeq2::DESeqDataSetFromTximport(txi_hog, colData = metadata, design = ~ 1)
  rounded_counts <- DESeq2::counts(dds)
  keep <- rowSums(rounded_counts >= 10) >= 5
  dds <- dds[keep, ]

  # Estimate sample normalisation from the complete retained HOG set. Any
  # B1_MAX_HOGS limit is only a development/testing limit on model fitting.
  dds <- DESeq2::estimateSizeFactors(dds)

  max_hogs <- suppressWarnings(as.integer(Sys.getenv("B1_MAX_HOGS", "0")))
  if (!is.na(max_hogs) && max_hogs > 0L && nrow(dds) > max_hogs) {
    dds <- dds[seq_len(max_hogs), ]
  }

  normalisation_factors <- DESeq2::normalizationFactors(dds)
  if (is.null(normalisation_factors)) {
    stop("DESeq2 did not construct HOG-by-sample normalisation factors")
  }
  offset_matrix <- log(normalisation_factors)
  if (any(!is.finite(offset_matrix))) stop("Non-finite values in HOG offset matrix")

  list(
    dds = dds,
    counts = DESeq2::counts(dds),
    normalisation_factors = normalisation_factors,
    offset = offset_matrix,
    base_mean = rowMeans(DESeq2::counts(dds, normalized = TRUE)),
    keep = keep
  )
}

fit_one_hog <- function(i, counts, offsets, metadata, model_formula, control) {
  if (identical(tolower(Sys.getenv("B1_VERBOSE_FITS", "false")), "true")) {
    message("Starting fit ", i, ": ", rownames(counts)[i])
  }
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

  if (identical(tolower(Sys.getenv("B1_VERBOSE_FITS", "false")), "true")) {
    message("Finished fit ", i, ": ", rownames(counts)[i])
  }

  coef_table <- tryCatch(coef(summary(fit))$cond, error = function(e) NULL)
  fixed <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
  list(
    success = !is.null(coef_table) && !is.null(fixed),
    coefficients = fixed,
    standard_errors = if (is.null(coef_table)) NULL else coef_table[, "Std. Error"],
    dispersion = tryCatch(glmmTMB::sigma(fit), error = function(e) NA_real_),
    aic = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
    random_intercept_variance = tryCatch(as.numeric(glmmTMB::VarCorr(fit)$cond[[1]][1, 1]), error = function(e) NA_real_),
    positive_definite_hessian = isTRUE(fit$sdr$pdHess),
    convergence_code = fit$fit$convergence,
    max_absolute_gradient = tryCatch(max(abs(fit$fit$gradient), na.rm = TRUE), error = function(e) NA_real_),
    message = fit$fit$message,
    warnings = paste(warning_messages, collapse = " | "),
    error = ""
  )
}

fit_all_hogs <- function(counts, offsets, metadata, model_formula, control, n_cores) {
  indices <- seq_len(nrow(counts))
  if (n_cores <= 1L) {
    return(lapply(indices, fit_one_hog, counts, offsets, metadata, model_formula, control))
  }
  cluster <- parallel::makePSOCKcluster(n_cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, suppressPackageStartupMessages(library(glmmTMB)))
  parallel::clusterExport(
    cluster,
    c("fit_one_hog", "counts", "offsets", "metadata", "model_formula", "control"),
    envir = environment()
  )
  parallel::parLapplyLB(cluster, indices, function(i) {
    fit_one_hog(i, counts, offsets, metadata, model_formula, control)
  })
}

ash_shrink_by_stage <- function(results) {
  results |>
    dplyr::group_split(stage, .keep = TRUE) |>
    purrr::map_dfr(function(stage_results) {
      ok <- is.finite(stage_results$log2FC_raw) &
        is.finite(stage_results$lfcSE_raw) & stage_results$lfcSE_raw > 0
      stage_results$log2FC_ash <- NA_real_
      stage_results$lfcSE_ash <- NA_real_
      stage_results$ash_lfsr <- NA_real_
      stage_results$ash_qvalue <- NA_real_
      if (sum(ok) >= 2L) {
        ash_fit <- ashr::ash(
          betahat = stage_results$log2FC_raw[ok],
          sebetahat = stage_results$lfcSE_raw[ok],
          mixcompdist = "normal",
          method = "shrink"
        )
        stage_results$log2FC_ash[ok] <- ashr::get_pm(ash_fit)
        stage_results$lfcSE_ash[ok] <- ashr::get_psd(ash_fit)
        stage_results$ash_lfsr[ok] <- ashr::get_lfsr(ash_fit)
        stage_results$ash_qvalue[ok] <- ashr::get_qvalue(ash_fit)
      }
      stage_results
    }) |>
    dplyr::arrange(stage, HOG)
}

make_member_outputs <- function(info, membership, hog_results) {
  gene_result_path <- file.path(
    gene_result_dir,
    paste0("glmmTMB_raw_counts_gene_offset_", info$species_code, "_results.tsv")
  )
  if (!file.exists(gene_result_path)) {
    stop("Gene-level result file not found: ", gene_result_path)
  }
  gene_results <- readr::read_tsv(gene_result_path, show_col_types = FALSE) |>
    dplyr::transmute(
      stage,
      gene_id = gene,
      gene_name,
      gene_baseMean = baseMean,
      member_log2FC_raw = log2FC_raw,
      member_lfcSE_raw = lfcSE_raw,
      member_pvalue = pvalue,
      member_padj_global = padj_global
    ) |>
    dplyr::rename(HOG_join_gene = gene_id)

  gene_results <- gene_results |>
    dplyr::group_split(stage, .keep = TRUE) |>
    purrr::map_dfr(function(stage_results) {
      ok <- is.finite(stage_results$member_log2FC_raw) &
        is.finite(stage_results$member_lfcSE_raw) & stage_results$member_lfcSE_raw > 0
      stage_results$member_log2FC_ash <- NA_real_
      stage_results$member_lfcSE_ash <- NA_real_
      stage_results$member_ash_lfsr <- NA_real_
      if (sum(ok) >= 2L) {
        ash_fit <- ashr::ash(
          stage_results$member_log2FC_raw[ok],
          stage_results$member_lfcSE_raw[ok],
          mixcompdist = "normal",
          method = "shrink"
        )
        stage_results$member_log2FC_ash[ok] <- ashr::get_pm(ash_fit)
        stage_results$member_lfcSE_ash[ok] <- ashr::get_psd(ash_fit)
        stage_results$member_ash_lfsr[ok] <- ashr::get_lfsr(ash_fit)
      }
      stage_results
    })

  species_membership <- membership |>
    dplyr::filter(species_code == info$species_code) |>
    dplyr::rename(HOG_join_gene = gene_id)

  member_long <- gene_results |>
    dplyr::inner_join(species_membership, by = "HOG_join_gene") |>
    dplyr::rename(gene_id = HOG_join_gene) |>
    dplyr::select(
      species_code, species_name, HOG, OG, stage, gene_id, gene_name,
      gene_baseMean, member_log2FC_raw, member_lfcSE_raw,
      member_log2FC_ash, member_lfcSE_ash, member_ash_lfsr,
      member_pvalue, member_padj_global
    ) |>
    dplyr::arrange(HOG, stage, gene_id)

  hog_sign <- hog_results |>
    dplyr::select(HOG, stage, HOG_log2FC_raw = log2FC_raw, HOG_log2FC_ash = log2FC_ash)

  member_work <- member_long |>
    dplyr::filter(
      is.finite(member_log2FC_raw),
      is.finite(member_lfcSE_raw),
      member_lfcSE_raw > 0
    ) |>
    dplyr::mutate(member_precision = 1 / member_lfcSE_raw^2)

  weighted_centres <- member_work |>
    dplyr::group_by(HOG, OG, stage) |>
    dplyr::summarise(
      precision_weighted_member_log2FC_raw =
        sum(member_precision * member_log2FC_raw) / sum(member_precision),
      .groups = "drop"
    )

  diagnostics_observed <- member_work |>
    dplyr::left_join(weighted_centres, by = c("HOG", "OG", "stage")) |>
    dplyr::group_by(HOG, OG, stage) |>
    dplyr::summarise(
      analysed_member_count = dplyr::n(),
      analysed_member_ids = collapse_ids(gene_id),
      mean_member_log2FC_raw = mean(member_log2FC_raw),
      median_member_log2FC_raw = stats::median(member_log2FC_raw),
      mean_member_log2FC_ash = mean(member_log2FC_ash, na.rm = TRUE),
      median_member_log2FC_ash = stats::median(member_log2FC_ash, na.rm = TRUE),
      precision_weighted_member_log2FC_raw =
        dplyr::first(precision_weighted_member_log2FC_raw),
      precision_weighted_member_lfcSE_nominal = sqrt(1 / sum(member_precision)),
      minimum_SE_representative_gene_id = gene_id[which.min(member_lfcSE_raw)],
      minimum_SE_representative_gene_name = gene_name[which.min(member_lfcSE_raw)],
      minimum_SE_representative_baseMean = gene_baseMean[which.min(member_lfcSE_raw)],
      minimum_SE_representative_log2FC_raw = member_log2FC_raw[which.min(member_lfcSE_raw)],
      minimum_SE_representative_lfcSE_raw = member_lfcSE_raw[which.min(member_lfcSE_raw)],
      minimum_SE_representative_log2FC_ash = member_log2FC_ash[which.min(member_lfcSE_raw)],
      directional_member_count = sum(member_log2FC_raw != 0),
      positive_member_count = sum(member_log2FC_raw > 0),
      negative_member_count = sum(member_log2FC_raw < 0),
      member_log2FC_raw_minimum = min(member_log2FC_raw),
      member_log2FC_raw_maximum = max(member_log2FC_raw),
      member_log2FC_raw_range = diff(range(member_log2FC_raw)),
      member_log2FC_ash_minimum = min(member_log2FC_ash, na.rm = TRUE),
      member_log2FC_ash_maximum = max(member_log2FC_ash, na.rm = TRUE),
      member_log2FC_ash_range = diff(range(member_log2FC_ash, na.rm = TRUE)),
      cochran_Q_descriptive = sum(
        member_precision *
          (member_log2FC_raw - dplyr::first(precision_weighted_member_log2FC_raw))^2
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      proportion_members_in_dominant_direction = dplyr::if_else(
        directional_member_count > 0,
        pmax(positive_member_count, negative_member_count) / directional_member_count,
        NA_real_
      ),
      precision_weighted_member_SD = sqrt(
        cochran_Q_descriptive /
          (1 / precision_weighted_member_lfcSE_nominal^2)
      ),
      heterogeneity_df = analysed_member_count - 1L,
      heterogeneity_I2_descriptive = dplyr::if_else(
        analysed_member_count > 1L & cochran_Q_descriptive > 0,
        pmax(0, (cochran_Q_descriptive - heterogeneity_df) / cochran_Q_descriptive),
        0
      )
    )

  genomic_counts <- species_membership |>
    dplyr::group_by(HOG, OG) |>
    dplyr::summarise(
      genomic_member_count = dplyr::n_distinct(HOG_join_gene),
      genomic_member_ids = collapse_ids(HOG_join_gene),
      .groups = "drop"
    )

  diagnostics <- tidyr::crossing(
    genomic_counts,
    stage = stages
  ) |>
    dplyr::left_join(diagnostics_observed, by = c("HOG", "OG", "stage")) |>
    dplyr::left_join(hog_sign, by = c("HOG", "stage")) |>
    dplyr::mutate(
      species_code = info$species_code,
      species_name = info$species_name,
      analysed_member_count = tidyr::replace_na(analysed_member_count, 0L),
      proportion_members_agreeing_with_HOG_direction = dplyr::case_when(
        !is.finite(HOG_log2FC_raw) | HOG_log2FC_raw == 0 |
          is.na(directional_member_count) | directional_member_count == 0 ~ NA_real_,
        HOG_log2FC_raw > 0 ~ positive_member_count / directional_member_count,
        HOG_log2FC_raw < 0 ~ negative_member_count / directional_member_count,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::relocate(species_code, species_name, HOG, OG, stage) |>
    dplyr::arrange(HOG, stage)

  list(long = member_long, diagnostics = diagnostics)
}

run_species <- function(info, mapping) {
  info <- as.list(info)
  message("Preparing N13 HOG inputs for ", info$species_name, "...")
  input <- read_species_data(info)
  txi_hog <- summarise_to_hog(input, mapping$membership, info$species_code)
  write_hog_matrix(txi_hog$counts, info$species_code, "counts_original_salmon")
  write_hog_matrix(txi_hog$abundance, info$species_code, "TPM")
  write_hog_matrix(txi_hog$length, info$species_code, "effective_lengths")

  prepared <- make_offset_data(txi_hog, input$metadata)
  counts <- prepared$counts
  offsets <- prepared$offset
  expected_coefficients <- paste0("larvalStage", stages, ":lifeStage", info$high_level)
  model_formula <- count ~ 0 + larvalStage + larvalStage:lifeStage +
    offset(offset_gs) + (1 | colonyNested)
  control <- glmmTMB::glmmTMBControl(
    optCtrl = list(iter.max = 1000, eval.max = 1000),
    rank_check = "warning"
  )

  message("Fitting ", nrow(counts), " N13 HOGs for ", info$species_name,
          " using ", n_cores, " core(s)...")
  start_time <- Sys.time()
  fit_list <- fit_all_hogs(counts, offsets, input$metadata, model_formula, control, n_cores)
  elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  names(fit_list) <- rownames(counts)

  all_coefficient_names <- unique(unlist(lapply(fit_list, function(x) names(x$coefficients))))
  coefficient_matrix <- matrix(
    NA_real_, nrow = nrow(counts), ncol = length(all_coefficient_names),
    dimnames = list(rownames(counts), all_coefficient_names)
  )
  se_matrix <- coefficient_matrix
  for (hog in rownames(counts)) {
    item <- fit_list[[hog]]
    if (isTRUE(item$success)) {
      coefficient_matrix[hog, names(item$coefficients)] <- item$coefficients
      se_matrix[hog, names(item$standard_errors)] <- item$standard_errors
    }
  }

  fit_diagnostics <- tibble::tibble(
    HOG = rownames(counts),
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
  valid_hogs <- fit_diagnostics$HOG[valid_fit]

  results <- dplyr::bind_rows(lapply(stages, function(stage) {
    coefficient <- paste0("larvalStage", stage, ":lifeStage", info$high_level)
    estimate_loge <- coefficient_matrix[valid_hogs, coefficient]
    se_loge <- se_matrix[valid_hogs, coefficient]
    z_value <- estimate_loge / se_loge
    p_value <- 2 * stats::pnorm(-abs(z_value))
    tibble::tibble(
      species_code = info$species_code,
      species_name = info$species_name,
      stage = stage,
      HOG = valid_hogs,
      baseMean = unname(prepared$base_mean[valid_hogs]),
      log2FC_raw = estimate_loge / log(2),
      lfcSE_raw = se_loge / log(2),
      z = z_value,
      pvalue = p_value,
      padj_stage_wise = p.adjust(p_value, method = "BH"),
      contrast = info$contrast
    )
  })) |>
    dplyr::mutate(
      padj_global = p.adjust(pvalue, method = "BH"),
      padj = if (global_BH) padj_global else padj_stage_wise
    ) |>
    ash_shrink_by_stage() |>
    dplyr::left_join(
      mapping$wide |>
        dplyr::select(
          HOG, OG, present_in_both_species, relationship_class,
          genomic_member_count_pd, genomic_member_count_vv
        ),
      by = "HOG"
    ) |>
    dplyr::mutate(
      method = "direct_N13_HOG_glmmTMB_raw_counts_length_offset",
      lfc_shrinkage = "ashr_normal_mixture_separate_by_species_and_stage"
    ) |>
    dplyr::arrange(stage, HOG)

  member_outputs <- make_member_outputs(info, mapping$membership, results)
  diagnostics_joined <- member_outputs$diagnostics |>
    dplyr::left_join(
      results |>
        dplyr::select(
          HOG, stage, HOG_baseMean = baseMean,
          HOG_lfcSE_raw = lfcSE_raw, HOG_lfcSE_ash = lfcSE_ash, HOG_pvalue = pvalue,
          HOG_padj_global = padj_global, HOG_ash_lfsr = ash_lfsr
        ),
      by = c("HOG", "stage")
    )

  readr::write_tsv(
    results,
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_DE_results.tsv")),
    na = "NA"
  )
  readr::write_tsv(
    fit_diagnostics,
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_fit_diagnostics.tsv")),
    na = "NA"
  )
  readr::write_tsv(
    member_outputs$long,
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_member_effects_long.tsv")),
    na = "NA"
  )
  readr::write_tsv(
    diagnostics_joined,
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_member_diagnostics.tsv")),
    na = "NA"
  )

  sample_offset_summary <- tibble::tibble(
    sample = colnames(counts),
    raw_HOG_library_size = colSums(counts),
    size_factor = DESeq2::sizeFactors(prepared$dds),
    geometric_mean_normalisation_factor = exp(colMeans(log(prepared$normalisation_factors))),
    minimum_normalisation_factor = apply(prepared$normalisation_factors, 2, min),
    maximum_normalisation_factor = apply(prepared$normalisation_factors, 2, max)
  )
  readr::write_tsv(
    sample_offset_summary,
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_sample_offsets.tsv"))
  )

  run_diagnostics <- tibble::tibble(
    species = info$species_name,
    model = paste(deparse(model_formula), collapse = " "),
    distribution = "negative binomial 2",
    test = "coefficient-wise asymptotic Wald z",
    aggregation = "tximport::summarizeToGene using gene-to-N13-HOG mapping",
    normalisation = "original Salmon counts plus HOG-by-sample effective-length and library-size offset",
    n_samples = ncol(counts),
    n_colonies_nested_in_condition = nlevels(input$metadata$colonyNested),
    n_genes_in_N13_mapping = sum(mapping$membership$species_code == info$species_code),
    n_genes_invalid_length = input$n_invalid_length_genes,
    n_HOGs_summarised = nrow(txi_hog$counts),
    n_HOGs_after_prefilter = nrow(counts),
    n_HOGs_fitted_and_valid = length(valid_hogs),
    n_fit_failures_or_invalid = nrow(counts) - length(valid_hogs),
    n_nonpositive_hessian = sum(!fit_diagnostics$positive_definite_hessian, na.rm = TRUE),
    n_stage_specific_tests = nrow(results),
    n_global_BH_0_05 = sum(results$padj_global < 0.05, na.rm = TRUE),
    elapsed_minutes = elapsed_minutes
  )
  readr::write_tsv(
    run_diagnostics,
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_run_diagnostics.tsv"))
  )

  saveRDS(
    list(
      species_info = info,
      txi_HOG = txi_hog,
      counts_fitted = counts,
      normalisation_factors = prepared$normalisation_factors,
      offsets = prepared$offset,
      metadata = input$metadata,
      coefficients = coefficient_matrix,
      standard_errors = se_matrix,
      fit_diagnostics = fit_diagnostics,
      results = results,
      member_effects = member_outputs$long,
      member_diagnostics = diagnostics_joined,
      run_diagnostics = run_diagnostics
    ),
    file.path(output_dir, paste0("N13_HOG_", info$species_code, "_analysis_objects.rds")),
    compress = "xz"
  )
  invisible(results)
}

combine_species_outputs <- function(mapping) {
  result_paths <- file.path(output_dir, paste0("N13_HOG_", c("pd", "vv"), "_DE_results.tsv"))
  diagnostic_paths <- file.path(output_dir, paste0("N13_HOG_", c("pd", "vv"), "_member_diagnostics.tsv"))
  member_paths <- file.path(output_dir, paste0("N13_HOG_", c("pd", "vv"), "_member_effects_long.tsv"))
  if (!all(file.exists(c(result_paths, diagnostic_paths, member_paths)))) {
    stop("Both species outputs must exist before combining them")
  }
  results_long <- purrr::map_dfr(result_paths, readr::read_tsv, show_col_types = FALSE)
  diagnostics_long <- purrr::map_dfr(diagnostic_paths, readr::read_tsv, show_col_types = FALSE)
  members_long <- purrr::map_dfr(member_paths, readr::read_tsv, show_col_types = FALSE)

  results_wide <- results_long |>
    dplyr::select(
      species_code, HOG, OG, stage, baseMean, log2FC_raw, lfcSE_raw,
      log2FC_ash, lfcSE_ash, ash_lfsr, pvalue, padj_global
    ) |>
    tidyr::pivot_wider(
      names_from = species_code,
      values_from = c(
        baseMean, log2FC_raw, lfcSE_raw, log2FC_ash, lfcSE_ash,
        ash_lfsr, pvalue, padj_global
      ),
      names_glue = "{.value}_{species_code}"
    ) |>
    dplyr::left_join(mapping$wide, by = c("HOG", "OG")) |>
    dplyr::arrange(stage, HOG)

  shared_results <- results_wide |>
    dplyr::filter(present_in_both_species, is.finite(log2FC_ash_pd), is.finite(log2FC_ash_vv))

  readr::write_tsv(results_long, file.path(output_dir, "N13_HOG_DE_results_both_species_long.tsv"), na = "NA")
  readr::write_tsv(results_wide, file.path(output_dir, "N13_HOG_DE_results_both_species_wide.tsv"), na = "NA")
  readr::write_tsv(shared_results, file.path(output_dir, "N13_HOG_DE_results_shared_for_regression.tsv"), na = "NA")
  readr::write_tsv(diagnostics_long, file.path(output_dir, "N13_HOG_member_diagnostics_both_species.tsv"), na = "NA")
  readr::write_tsv(members_long, file.path(output_dir, "N13_HOG_member_effects_both_species_long.tsv"), na = "NA")

  # Machine-readable key dataframes for fast reload in later R sessions.
  saveRDS(results_long, file.path(output_dir, "N13_HOG_DE_results_both_species_long.rds"), compress = "xz")
  saveRDS(results_wide, file.path(output_dir, "N13_HOG_DE_results_both_species_wide.rds"), compress = "xz")
  saveRDS(shared_results, file.path(output_dir, "N13_HOG_DE_results_shared_for_regression.rds"), compress = "xz")
  saveRDS(diagnostics_long, file.path(output_dir, "N13_HOG_member_diagnostics_both_species.rds"), compress = "xz")
  saveRDS(members_long, file.path(output_dir, "N13_HOG_member_effects_both_species_long.rds"), compress = "xz")
  saveRDS(mapping$wide, file.path(output_dir, "N13_HOG_membership_summary.rds"), compress = "xz")
  saveRDS(mapping$membership, file.path(output_dir, "N13_HOG_membership_long.rds"), compress = "xz")

  HOG_membership <- mapping$wide
  HOG_membership_long <- mapping$membership
  save(
    results_long, results_wide, shared_results, diagnostics_long, members_long,
    HOG_membership, HOG_membership_long,
    file = file.path(output_dir, "N13_HOG_combined_analysis_workspace.RData"),
    compress = "xz"
  )

  workbook_path <- file.path(output_dir, "N13_HOG_differential_expression_and_member_diagnostics.xlsx")
  file_index <- tibble::tribble(
    ~content, ~file,
    "All HOG-level tests, long format", "N13_HOG_DE_results_both_species_long.tsv",
    "All HOG-level tests, wide format", "N13_HOG_DE_results_both_species_wide.tsv",
    "Shared HOG-level tests for cross-species regression", "N13_HOG_DE_results_shared_for_regression.tsv",
    "All member-gene effects", "N13_HOG_member_effects_both_species_long.tsv",
    "All member-level HOG diagnostics", "N13_HOG_member_diagnostics_both_species.tsv",
    "Combined reloadable R workspace", "N13_HOG_combined_analysis_workspace.RData",
    "Polistes complete analysis object", "N13_HOG_pd_analysis_objects.rds",
    "Vespula complete analysis object", "N13_HOG_vv_analysis_objects.rds"
  )
  openxlsx::write.xlsx(
    list(
      HOG_DE_shared = shared_results,
      HOG_membership = mapping$wide,
      file_index = file_index
    ),
    workbook_path,
    overwrite = TRUE
  )
  invisible(shared_results)
}

mapping <- read_n13_mapping()
write_mapping <- !identical(tolower(Sys.getenv("B1_WRITE_MAPPING", "true")), "false")
if (write_mapping) write_mapping_outputs(mapping)

prepare_only <- identical(tolower(Sys.getenv("B1_PREPARE_ONLY", "false")), "true")
combine_only <- identical(tolower(Sys.getenv("B1_COMBINE_ONLY", "false")), "true")
species_to_run <- tolower(Sys.getenv("B1_SPECIES", "both"))
if (!species_to_run %in% c("both", "pd", "vv")) {
  stop("B1_SPECIES must be one of: both, pd, vv")
}

if (prepare_only) {
  message("N13 HOG mapping outputs written to: ", output_dir)
} else if (combine_only) {
  combine_species_outputs(mapping)
  message("Combined N13 HOG outputs written to: ", output_dir)
} else {
  if (species_to_run %in% c("both", "pd")) {
    run_species(species_info |> dplyr::filter(species_code == "pd"), mapping)
  }
  if (species_to_run %in% c("both", "vv")) {
    run_species(species_info |> dplyr::filter(species_code == "vv"), mapping)
  }
  if (species_to_run == "both") combine_species_outputs(mapping)
}

analysis_metadata <- tibble::tibble(
  key = c(
    "script", "analysis_date", "HOG_node", "HOG_source", "global_BH",
    "model", "random_effect", "distribution", "aggregation", "normalisation",
    "response", "test", "HOG_ASH", "member_ASH", "cores", "species_run",
    "tximport_version", "ashr_version", "DESeq2_version", "glmmTMB_version", "R_version"
  ),
  value = c(
    "B1_analysis_differential_expression_N13_HOG_glmmTMB_raw_counts_gene_offset.R",
    as.character(Sys.Date()), "N13 (MRCA of sampled Polistinae and Vespinae)",
    "nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/results_EXCON/orthofinder/Phylogenetic_Hierarchical_Orthogroups/N13.tsv",
    as.character(global_BH),
    "count ~ 0 + larvalStage + larvalStage:lifeStage + offset(offset_gs) + (1 | colonyNested)",
    "colony nested within lifeStage/condition", "glmmTMB nbinom2",
    "tximport::summarizeToGene from genes to N13 HOGs",
    "DESeq2/tximport HOG-by-sample effective-length and library-size normalisation factors",
    "rounded summed original Salmon estimated counts", "coefficient-wise asymptotic Wald z",
    "ashr normal-mixture shrinkage separately by species and stage",
    "ashr normal-mixture shrinkage of member-gene coefficients separately by species and stage",
    as.character(n_cores), species_to_run,
    as.character(utils::packageVersion("tximport")),
    as.character(utils::packageVersion("ashr")),
    as.character(utils::packageVersion("DESeq2")),
    as.character(utils::packageVersion("glmmTMB")),
    R.version.string
  )
)
metadata_suffix <- if (prepare_only) {
  "mapping"
} else if (combine_only) {
  "combined"
} else {
  species_to_run
}
readr::write_tsv(
  analysis_metadata,
  file.path(output_dir, paste0("analysis_metadata_", metadata_suffix, ".tsv"))
)
