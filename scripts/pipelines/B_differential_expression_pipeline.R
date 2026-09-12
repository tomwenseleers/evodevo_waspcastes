# B_differential expression analysis.R
# Final transcriptomic analyses for the Vespula-Polistes caste evo-devo study.
#
# This script uses the final gene-level and N13 hierarchical-orthogroup (HOG)
# glmmTMB analyses. The computationally expensive negative-binomial mixed
# models are loaded from their saved outputs; set rerun_mixed_models to TRUE
# only when deliberately regenerating those upstream analyses.

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(mixOmics)
  library(export)
  library(ggthemes)
  library(RColorBrewer)
  library(scales)
  library(topGO)
  library(ontologyIndex)
  library(afex)
  library(estimatr)
  library(emmeans)
  library(marginaleffects)
  library(robustbase)
  library(glmnet)
  library(gtools)
  library(pheatmap)
  # topGO/AnnotationDbi export a dplyr::select() generic; keep dplyr verbs active.
  library(dplyr)
})

set.seed(20260826)
options(stringsAsFactors = FALSE)

if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  stop("PROJECT_ROOT must be set by B_differential expression analysis.R before sourcing this pipeline.")
}
dir_base <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
dir_annotations <- file.path(dir_base, "input_annotation")
dir_de <- file.path(dir_base, "input_differential_expression")
dir_cafe_inputs <- file.path(dir_base, "input_CAFE")
dir_go <- file.path(dir_annotations, "raw", "GO")
output_dir <- file.path(dir_base, "output")
figure_dir <- file.path(output_dir, "figures")
supplemental_figure_dir <- file.path(output_dir, "supplemental_figures")
supplemental_table_dir <- file.path(output_dir, "supplemental_tables")
invisible(lapply(
  c(output_dir, figure_dir, supplemental_figure_dir, supplemental_table_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

stages <- c("L1", "L2", "L3", "L4", "L5", "P")
alpha_de <- 0.05

f_genetable_pd <- file.path(dir_cafe_inputs, "annotation", "pd_genetable.csv")
f_genetable_vv <- file.path(dir_cafe_inputs, "annotation", "vv_genetable.csv")
f_obo <- file.path(dir_go, "go-basic.obo")
go_annot_orthologs <- file.path(
  dir_base, "output_annotation", "N13_HOG_GO_final_long.tsv.gz"
)
n13_hog_file <- file.path(
  dir_base, "nextflow_runs", "2_EXCON",
  "2_EXCON_orthofinder_eggnogmapper_run", "results_EXCON", "orthofinder",
  "Phylogenetic_Hierarchical_Orthogroups", "N13.tsv"
)

gene_model_script <- file.path(dir_base, "scripts", "upstream", "B1_gene_glmmTMB_raw_counts_length_offset.R")
hog_model_script <- file.path(dir_base, "scripts", "upstream", "B1_N13_HOG_glmmTMB_raw_counts_length_offset.R")

gene_result_dir <- file.path(dir_de, "checkpoints", "gene")
hog_result_dir <- file.path(dir_de, "checkpoints", "HOG")

gene_object_pd_file <- file.path(gene_result_dir, "glmmTMB_raw_counts_gene_offset_pd_analysis_objects.rds")
gene_object_vv_file <- file.path(gene_result_dir, "glmmTMB_raw_counts_gene_offset_vv_analysis_objects.rds")
hog_object_pd_file <- file.path(hog_result_dir, "N13_HOG_pd_analysis_objects.rds")
hog_object_vv_file <- file.path(hog_result_dir, "N13_HOG_vv_analysis_objects.rds")
hog_shared_file <- file.path(hog_result_dir, "N13_HOG_DE_results_shared_for_regression.rds")
hog_membership_file <- file.path(hog_result_dir, "N13_HOG_membership_long.rds")
hog_annotation_file <- file.path(
  dir_annotations, "reference", "pd_vv_N13_HOGs_annotated.rds"
)

required_inputs <- c(
  gene_model_script, hog_model_script,
  gene_object_pd_file, gene_object_vv_file,
  hog_object_pd_file, hog_object_vv_file,
  hog_shared_file, hog_membership_file, hog_annotation_file,
  f_genetable_pd, f_genetable_vv, f_obo, go_annot_orthologs, n13_hog_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required input files:\n", paste(missing_inputs, collapse = "\n"))
}


# 1. DIFFERENTIAL EXPRESSION ANALYSIS ####

## 1.1 Final gene-level glmmTMB analysis ####

# Upstream model fitted separately to each gene and species:
#   count_gs ~ 0 + larvalStage + larvalStage:lifeStage +
#              offset(offset_gs) + (1 | colonyNested)
# using glmmTMB::nbinom2. The response comprises rounded original Salmon
# estimated counts. offset_gs is the log of a gene-by-sample normalisation
# factor combining sample-specific effective transcript length and library
# size. colonyNested is colony nested within condition (E/L or W/Q).
# Coefficient-wise Wald tests compare L versus E in Polistes and Q versus W in
# Vespula within each larval stage. BH correction is applied globally across
# all gene-by-stage tests within species.

## 1.2 Load the saved final gene- and HOG-level analyses ####

if (!exists("rerun_mixed_models", inherits = TRUE)) rerun_mixed_models <- FALSE
if (isTRUE(rerun_mixed_models)) {
  stop(
    "The mixed-model fits are intentionally not rerun inside B_differential expression analysis.R. ",
    "Run these one-worker scripts independently, then restart this script:\n",
    gene_model_script, "\n", hog_model_script
  )
}

gene_obj_pd <- readRDS(gene_object_pd_file)
gene_obj_vv <- readRDS(gene_object_vv_file)
hog_obj_pd <- readRDS(hog_object_pd_file)
hog_obj_vv <- readRDS(hog_object_vv_file)

gene_de_pd <- as_tibble(gene_obj_pd$results) %>%
  mutate(stage = factor(stage, levels = stages))
gene_de_vv <- as_tibble(gene_obj_vv$results) %>%
  mutate(stage = factor(stage, levels = stages))

# Compatibility aliases for code and interactive inspection. These are final
# glmmTMB results, not DESeq2 results.
deseq_pd <- gene_de_pd
deseq_vv <- gene_de_vv

metadata_pd <- as.data.frame(hog_obj_pd$metadata)
metadata_vv <- as.data.frame(hog_obj_vv$metadata)
metadata_pd$larvalStage <- factor(metadata_pd$larvalStage, levels = stages)
metadata_vv$larvalStage <- factor(metadata_vv$larvalStage, levels = stages)

# Direct N13 HOG-level Salmon summaries, analogous to
# tximport::summarizeToGene(), are the PLS input matrices.
cts_pd <- hog_obj_pd$txi_HOG$counts
cts_vv <- hog_obj_vv$txi_HOG$counts
storage.mode(cts_pd) <- "double"
storage.mode(cts_vv) <- "double"

hog_de_shared <- readRDS(hog_shared_file) %>%
  mutate(stage = factor(stage, levels = stages))
hog_membership <- readRDS(hog_membership_file)
orthologs <- readRDS(hog_annotation_file) %>%
  mutate(orthogroup = as.character(orthogroup))

genetable_pd <- read_csv(f_genetable_pd, show_col_types = FALSE)
genetable_vv <- read_csv(f_genetable_vv, show_col_types = FALSE)


## 1.2.1 N13 HOG orthology-composition summary ####

# Report the unfiltered N13 mapping used to construct the cross-species
# transcriptomic units. These counts precede expression and PLS filtering.
count_hog_members <- function(x) {
  x <- if_else(is.na(x), "", str_trim(x))
  if_else(x == "", 0L, as.integer(str_count(x, fixed(",")) + 1L))
}

n13_raw <- read_tsv(
  n13_hog_file, col_types = cols(.default = col_character()),
  show_col_types = FALSE, progress = FALSE,
  na = character()
)
pd_n13_column <- "Polistes_dominula.clean"
vv_n13_column <- "Vespula_vulgaris.clean"
if (!all(c(pd_n13_column, vv_n13_column) %in% names(n13_raw))) {
  stop("The focal-species columns are missing from N13.tsv.")
}

n13_composition <- n13_raw %>%
  transmute(
    HOG,
    Polistes_members = count_hog_members(.data[[pd_n13_column]]),
    Vespula_members = count_hog_members(.data[[vv_n13_column]])
  )

n13_total_hogs <- nrow(n13_composition)
n13_shared_hogs <- sum(
  n13_composition$Polistes_members > 0L &
    n13_composition$Vespula_members > 0L
)

make_n13_summary_rows <- function(summary_group, metric, count, denominator,
                                  definition) {
  tibble(
    analysis_scope = "Unfiltered N13.tsv mapping before expression and PLS filtering",
    summary_group = summary_group,
    metric = metric,
    count = as.integer(count),
    denominator = as.integer(denominator),
    percentage = if_else(
      is.na(denominator) | denominator == 0L,
      NA_real_, round(100 * count / denominator, 2)
    ),
    definition = definition
  )
}

n13_orthology_composition_summary <- bind_rows(
  make_n13_summary_rows(
    "HOG presence", "All N13 HOGs", n13_total_hogs, n13_total_hogs,
    "All nonoverlapping HOGs defined at OrthoFinder node N13."
  ),
  make_n13_summary_rows(
    "HOG presence", "HOGs containing Polistes dominula",
    sum(n13_composition$Polistes_members > 0L), n13_total_hogs,
    "N13 HOGs with at least one P. dominula gene."
  ),
  make_n13_summary_rows(
    "HOG presence", "HOGs containing Vespula vulgaris",
    sum(n13_composition$Vespula_members > 0L), n13_total_hogs,
    "N13 HOGs with at least one V. vulgaris gene."
  ),
  make_n13_summary_rows(
    "HOG presence", "HOGs represented in both focal species",
    n13_shared_hogs, n13_total_hogs,
    "N13 HOGs with at least one gene from each focal species."
  ),
  make_n13_summary_rows(
    "HOG presence", "HOGs represented only in Polistes among the focal species",
    sum(n13_composition$Polistes_members > 0L &
          n13_composition$Vespula_members == 0L), n13_total_hogs,
    "May also contain genes from nonfocal taxa in the OrthoFinder analysis."
  ),
  make_n13_summary_rows(
    "HOG presence", "HOGs represented only in Vespula among the focal species",
    sum(n13_composition$Polistes_members == 0L &
          n13_composition$Vespula_members > 0L), n13_total_hogs,
    "May also contain genes from nonfocal taxa in the OrthoFinder analysis."
  ),
  make_n13_summary_rows(
    "HOG presence", "HOGs containing neither focal species",
    sum(n13_composition$Polistes_members == 0L &
          n13_composition$Vespula_members == 0L), n13_total_hogs,
    "N13 HOGs represented only by other taxa in the OrthoFinder analysis."
  ),
  make_n13_summary_rows(
    "Shared-HOG composition", "One Polistes gene : one Vespula gene",
    sum(n13_composition$Polistes_members == 1L &
          n13_composition$Vespula_members == 1L), n13_shared_hogs,
    "One-to-one HOG composition between the two focal species."
  ),
  make_n13_summary_rows(
    "Shared-HOG composition", "One Polistes gene : multiple Vespula genes",
    sum(n13_composition$Polistes_members == 1L &
          n13_composition$Vespula_members > 1L), n13_shared_hogs,
    "One-to-many HOG composition in the Polistes-to-Vespula direction."
  ),
  make_n13_summary_rows(
    "Shared-HOG composition", "Multiple Polistes genes : one Vespula gene",
    sum(n13_composition$Polistes_members > 1L &
          n13_composition$Vespula_members == 1L), n13_shared_hogs,
    "Many-to-one HOG composition in the Polistes-to-Vespula direction."
  ),
  make_n13_summary_rows(
    "Shared-HOG composition", "Multiple genes in both focal species",
    sum(n13_composition$Polistes_members > 1L &
          n13_composition$Vespula_members > 1L), n13_shared_hogs,
    "Many-to-many HOG composition between the two focal species."
  ),
  make_n13_summary_rows(
    "Mapped genes", "Polistes dominula genes assigned to N13 HOGs",
    sum(n13_composition$Polistes_members), NA_integer_,
    "Total P. dominula gene memberships in N13.tsv."
  ),
  make_n13_summary_rows(
    "Mapped genes", "Vespula vulgaris genes assigned to N13 HOGs",
    sum(n13_composition$Vespula_members), NA_integer_,
    "Total V. vulgaris gene memberships in N13.tsv."
  )
)

if (sum(n13_orthology_composition_summary$count[
      n13_orthology_composition_summary$summary_group ==
        "Shared-HOG composition"
    ]) != n13_shared_hogs) {
  stop("Shared N13 HOG relationship classes do not sum to the shared-HOG total.")
}

n13_orthology_summary_file <- file.path(
  output_dir, "N13_HOG_orthology_composition_summary.tsv"
)
write_tsv(n13_orthology_composition_summary, n13_orthology_summary_file)


## 1.3 Complete gene- and N13 HOG-level differential-expression exports ####

# These repository tables contain every fitted stage contrast. The primary
# significance field uses the manuscript-wide correction: BH across all
# stage-by-gene or stage-by-HOG tests within each species.
collapse_unique_values <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) NA_character_ else paste(x, collapse = "; ")
}

gene_annotation <- bind_rows(
  genetable_pd %>%
    transmute(
      species_code = "pd", gene_id = gene_pd,
      gene_symbol = symbol_pd, gene_description = description_pd,
      gene_type = genetype_pd, n_transcripts = ntranscripts_pd
    ),
  genetable_vv %>%
    transmute(
      species_code = "vv", gene_id = gene_vv,
      gene_symbol = symbol_vv, gene_description = description_vv,
      gene_type = genetype_vv, n_transcripts = ntranscripts_vv
    )
) %>%
  distinct(species_code, gene_id, .keep_all = TRUE)

hog_metadata <- bind_rows(
  as_tibble(hog_obj_pd$results),
  as_tibble(hog_obj_vv$results)
) %>%
  transmute(
    HOG, OG, present_in_both_species, relationship_class,
    genomic_member_count_pd, genomic_member_count_vv
  ) %>%
  distinct(HOG, .keep_all = TRUE)

hog_member_annotations <- hog_membership %>%
  left_join(gene_annotation, by = c("species_code", "gene_id")) %>%
  group_by(HOG, OG) %>%
  summarise(
    Polistes_member_gene_ids = collapse_unique_values(
      gene_id[species_code == "pd"]
    ),
    Polistes_member_gene_symbols = collapse_unique_values(
      gene_symbol[species_code == "pd"]
    ),
    Polistes_member_gene_descriptions = collapse_unique_values(
      gene_description[species_code == "pd"]
    ),
    Vespula_member_gene_ids = collapse_unique_values(
      gene_id[species_code == "vv"]
    ),
    Vespula_member_gene_symbols = collapse_unique_values(
      gene_symbol[species_code == "vv"]
    ),
    Vespula_member_gene_descriptions = collapse_unique_values(
      gene_description[species_code == "vv"]
    ),
    .groups = "drop"
  )

prepare_gene_de_export <- function(dat, species_code_value, species_value,
                                   contrast_value, positive_label,
                                   negative_label) {
  species_hog_map <- hog_membership %>%
    filter(species_code == species_code_value) %>%
    transmute(gene = gene_id, HOG, OG) %>%
    distinct(gene, .keep_all = TRUE)

  dat %>%
    mutate(stage = as.character(stage)) %>%
    left_join(species_hog_map, by = "gene") %>%
    left_join(hog_metadata, by = c("HOG", "OG")) %>%
    left_join(
      gene_annotation %>% filter(species_code == species_code_value),
      by = c("gene" = "gene_id")
    ) %>%
    mutate(
      global_FDR_lt_0.05 = is.finite(padj_global) & padj_global < alpha_de,
      effect_direction = case_when(
        !is.finite(log2FC_raw) ~ NA_character_,
        log2FC_raw > 0 ~ positive_label,
        log2FC_raw < 0 ~ negative_label,
        TRUE ~ "No estimated difference"
      ),
      significant_effect = if_else(
        global_FDR_lt_0.05,
        paste("Significant:", effect_direction),
        "Not significant at global FDR < 0.05"
      )
    ) %>%
    transmute(
      species_code = species_code_value,
      species = species_value,
      contrast = contrast_value,
      stage,
      gene_id = gene,
      gene_symbol,
      gene_description,
      gene_type,
      n_transcripts,
      N13_HOG = HOG,
      orthogroup = OG,
      HOG_present_in_both_species = present_in_both_species,
      HOG_relationship_class = relationship_class,
      HOG_member_count_Polistes = genomic_member_count_pd,
      HOG_member_count_Vespula = genomic_member_count_vv,
      base_mean = baseMean,
      log2FC_raw,
      SE_raw = lfcSE_raw,
      Wald_z = z,
      P_value = pvalue,
      stage_wise_FDR = padj_stage_wise,
      global_FDR = padj_global,
      global_FDR_lt_0.05,
      effect_direction,
      significant_effect,
      log2FC_ASH_shrunk = log2FC_shrunk,
      SE_ASH_shrunk = lfcSE_shrunk,
      analysis_method = method,
      shrinkage_method = lfc_shrinkage
    )
}

prepare_hog_de_export <- function(dat, species_code_value, species_value,
                                  contrast_value, positive_label,
                                  negative_label) {
  dat %>%
    mutate(stage = as.character(stage)) %>%
    left_join(hog_member_annotations, by = c("HOG", "OG")) %>%
    mutate(
      global_FDR_lt_0.05 = is.finite(padj_global) & padj_global < alpha_de,
      effect_direction = case_when(
        !is.finite(log2FC_raw) ~ NA_character_,
        log2FC_raw > 0 ~ positive_label,
        log2FC_raw < 0 ~ negative_label,
        TRUE ~ "No estimated difference"
      ),
      significant_effect = if_else(
        global_FDR_lt_0.05,
        paste("Significant:", effect_direction),
        "Not significant at global FDR < 0.05"
      )
    ) %>%
    transmute(
      species_code = species_code_value,
      species = species_value,
      contrast = contrast_value,
      stage,
      N13_HOG = HOG,
      orthogroup = OG,
      present_in_both_species,
      relationship_class,
      HOG_member_count_Polistes = genomic_member_count_pd,
      HOG_member_count_Vespula = genomic_member_count_vv,
      Polistes_member_gene_ids,
      Polistes_member_gene_symbols,
      Polistes_member_gene_descriptions,
      Vespula_member_gene_ids,
      Vespula_member_gene_symbols,
      Vespula_member_gene_descriptions,
      base_mean = baseMean,
      log2FC_raw,
      SE_raw = lfcSE_raw,
      Wald_z = z,
      P_value = pvalue,
      stage_wise_FDR = padj_stage_wise,
      global_FDR = padj_global,
      global_FDR_lt_0.05,
      effect_direction,
      significant_effect,
      log2FC_ASH_shrunk = log2FC_ash,
      SE_ASH_shrunk = lfcSE_ash,
      ASH_local_false_sign_rate = ash_lfsr,
      ASH_q_value = ash_qvalue,
      analysis_method = method,
      shrinkage_method = lfc_shrinkage
    )
}

full_gene_de_results <- bind_rows(
  prepare_gene_de_export(
    gene_de_pd, "pd", "Polistes dominula",
    "Late versus early season", "Higher in late-season females",
    "Higher in early-season females"
  ),
  prepare_gene_de_export(
    gene_de_vv, "vv", "Vespula vulgaris",
    "Queen versus worker", "Higher in queen-destined females",
    "Higher in worker-destined females"
  )
) %>%
  arrange(factor(species_code, c("pd", "vv")), factor(stage, stages), gene_id)

full_hog_de_results <- bind_rows(
  prepare_hog_de_export(
    as_tibble(hog_obj_pd$results), "pd", "Polistes dominula",
    "Late versus early season", "Higher in late-season females",
    "Higher in early-season females"
  ),
  prepare_hog_de_export(
    as_tibble(hog_obj_vv$results), "vv", "Vespula vulgaris",
    "Queen versus worker", "Higher in queen-destined females",
    "Higher in worker-destined females"
  )
) %>%
  arrange(factor(species_code, c("pd", "vv")), factor(stage, stages), N13_HOG)

full_gene_de_file <- file.path(
  output_dir, "full_gene_level_differential_expression_results.tsv.gz"
)
full_hog_de_file <- file.path(
  output_dir, "full_N13_HOG_level_differential_expression_results.tsv.gz"
)
write_tsv(full_gene_de_results, full_gene_de_file)
write_tsv(full_hog_de_results, full_hog_de_file)

full_de_export_index <- bind_rows(
  tibble(
    analysis_unit = "Gene",
    file = basename(full_gene_de_file),
    rows = nrow(full_gene_de_results),
    significant_global_FDR_lt_0.05 = sum(
      full_gene_de_results$global_FDR_lt_0.05, na.rm = TRUE
    ),
    compressed_bytes = file.info(full_gene_de_file)$size
  ),
  tibble(
    analysis_unit = "N13 HOG",
    file = basename(full_hog_de_file),
    rows = nrow(full_hog_de_results),
    significant_global_FDR_lt_0.05 = sum(
      full_hog_de_results$global_FDR_lt_0.05, na.rm = TRUE
    ),
    compressed_bytes = file.info(full_hog_de_file)$size
  )
)
write_tsv(full_de_export_index, file.path(output_dir, "full_DE_results_file_index.tsv"))


## 1.4.1 N13 HOG-level cross-species differential-expression table ####

make_species_hog_wide <- function(dat, species = c("Pd", "Vv")) {
  species <- match.arg(species)
  if (species == "Pd") {
    long <- dat %>%
      transmute(
        orthogroup = HOG, stage,
        baseMean = baseMean_pd,
        log2FCraw = log2FC_raw_pd,
        lfcSEraw = lfcSE_raw_pd,
        log2FCshrunk = log2FC_ash_pd,
        lfcSEshrunk = lfcSE_ash_pd,
        pvalue = pvalue_pd,
        padj = padj_global_pd,
        ash_lfsr = ash_lfsr_pd
      )
    suffix <- "Pd_{stage}_LvsE"
  } else {
    long <- dat %>%
      transmute(
        orthogroup = HOG, stage,
        baseMean = baseMean_vv,
        log2FCraw = log2FC_raw_vv,
        lfcSEraw = lfcSE_raw_vv,
        log2FCshrunk = log2FC_ash_vv,
        lfcSEshrunk = lfcSE_ash_vv,
        pvalue = pvalue_vv,
        padj = padj_global_vv,
        ash_lfsr = ash_lfsr_vv
      )
    suffix <- "Vv_{stage}_QvsW"
  }

  long %>%
    mutate(stage = as.character(stage)) %>%
    pivot_wider(
      names_from = stage,
      values_from = c(baseMean, log2FCraw, lfcSEraw, log2FCshrunk,
                      lfcSEshrunk, pvalue, padj, ash_lfsr),
      names_glue = paste0("{.value}_", suffix)
    )
}

hog_de_pd_wide <- make_species_hog_wide(hog_de_shared, "Pd")
hog_de_vv_wide <- make_species_hog_wide(hog_de_shared, "Vv")

ortholog_de_wide <- orthologs %>%
  left_join(hog_de_pd_wide, by = "orthogroup") %>%
  left_join(hog_de_vv_wide, by = "orthogroup")

write_tsv(ortholog_de_wide, file.path(output_dir, "N13_HOG_DE_results_annotated_wide.tsv"))
saveRDS(ortholog_de_wide, file.path(output_dir, "N13_HOG_DE_results_annotated_wide.rds"))
write_tsv(hog_de_shared, file.path(output_dir, "N13_HOG_DE_results_shared_long.tsv"))


## 1.5 Stacked column plot of final gene-level DE results ####

pd_hog_map <- hog_membership %>%
  filter(species_code == "pd") %>%
  transmute(HOG, gene = gene_id) %>%
  distinct()
vv_hog_map <- hog_membership %>%
  filter(species_code == "vv") %>%
  transmute(HOG, gene = gene_id) %>%
  distinct()
shared_species_hogs <- hog_membership %>%
  distinct(HOG, species_code) %>%
  count(HOG, name = "n_species") %>%
  filter(n_species == 2L) %>%
  pull(HOG)

sig_pd_hog_stage <- gene_de_pd %>%
  filter(!is.na(padj_global), padj_global < alpha_de) %>%
  inner_join(pd_hog_map, by = "gene") %>%
  distinct(HOG, stage)
sig_vv_hog_stage <- gene_de_vv %>%
  filter(!is.na(padj_global), padj_global < alpha_de) %>%
  inner_join(vv_hog_map, by = "gene") %>%
  distinct(HOG, stage)

is_noncoding <- function(x) {
  str_detect(
    str_to_lower(coalesce(as.character(x), "")),
    "(^|[^a-z])(lncrna|ncrna|mirna|snorna|snrna|rrna|trna|non[- ]coding)([^a-z]|$)"
  )
}

classify_significant_genes <- function(de, gene_table, hog_map, other_sig,
                                       species, gene_type_col) {
  gt <- gene_table %>%
    transmute(gene = .data[[paste0("gene_", tolower(species))]],
              gene_type = .data[[gene_type_col]])

  de %>%
    filter(!is.na(padj_global), padj_global < alpha_de) %>%
    dplyr::select(stage, gene, log2FC_raw, log2FC_shrunk, lfcSE_raw,
           pvalue, padj_global) %>%
    left_join(gt, by = "gene") %>%
    left_join(hog_map, by = "gene") %>%
    left_join(other_sig %>% mutate(shared_de = TRUE), by = c("HOG", "stage")) %>%
    mutate(
      species = species,
      category = case_when(
        is_noncoding(gene_type) ~ "Non-coding RNA",
        is.na(HOG) ~ "No N13 HOG",
        coalesce(shared_de, FALSE) ~ "DE in both species",
        TRUE ~ "DE in focal species only"
      ),
      stage = factor(stage, levels = stages)
    )
}

de_gene_categories <- bind_rows(
  classify_significant_genes(
    gene_de_pd, genetable_pd, pd_hog_map, sig_vv_hog_stage,
    species = "Pd", gene_type_col = "genetype_pd"
  ),
  classify_significant_genes(
    gene_de_vv, genetable_vv, vv_hog_map, sig_pd_hog_stage,
    species = "Vv", gene_type_col = "genetype_vv"
  )
) %>%
  mutate(
    species = factor(species, levels = c("Pd", "Vv"),
                     labels = c("Polistes dominula", "Vespula vulgaris")),
    category = factor(
      category,
      levels = c("DE in both species", "DE in focal species only",
                 "No N13 HOG", "Non-coding RNA")
    )
  )

de_stacked_counts <- de_gene_categories %>%
  count(species, stage, category, name = "n_genes") %>%
  complete(species, stage, category, fill = list(n_genes = 0L))

# Deliberately crude visual summary: collapse shared versus focal-only DE genes
# because this panel describes annotation/orthology coverage, not effect-size or
# cross-species concordance. Those questions are analysed at N13 HOG level in
# sections 4 and 5.
de_stacked_counts_plot <- de_gene_categories %>%
  mutate(
    plot_category = case_when(
      category == "Non-coding RNA" ~ "ncRNA",
      HOG %in% shared_species_hogs ~ "with orthology",
      TRUE ~ "without orthology"
    ),
    plot_category = factor(
      plot_category,
      levels = c("ncRNA", "without orthology", "with orthology")
    )
  ) %>%
  count(species, stage, plot_category, name = "n_genes") %>%
  complete(species, stage, plot_category, fill = list(n_genes = 0L))

write_tsv(de_gene_categories, file.path(output_dir, "gene_level_significant_DE_classification.tsv"))
write_tsv(de_stacked_counts, file.path(output_dir, "gene_level_significant_DE_stacked_counts.tsv"))
write_tsv(de_stacked_counts_plot, file.path(output_dir, "gene_level_significant_DE_stacked_counts_plot.tsv"))

p_de_stacked <- ggplot(
  de_stacked_counts_plot,
  aes(stage, n_genes, fill = plot_category)
) +
  geom_col(width = 0.90, colour = NA) +
  facet_wrap(~ species, nrow = 1, scales = "fixed") +
  scale_fill_manual(
    values = c(
      "ncRNA" = "#D8E8F4",
      "without orthology" = "#97C4DD",
      "with orthology" = "#3787BD"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Stage", y = "Number of DE genes", fill = NULL) +
  theme_few(base_size = 15) +
  theme(
    strip.text = element_text(face = "italic", size = 16),
    strip.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 0.47),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.47),
    panel.spacing = grid::unit(0.15, "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.18, 0.88),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.height = grid::unit(0.55, "cm"),
    legend.key.width = grid::unit(0.55, "cm")
  )

print(p_de_stacked)
graph2png(x = ggplot2::last_plot(), file = file.path(supplemental_figure_dir, "FigS2.png"), width = 8, height = 6)
graph2pdf(x = ggplot2::last_plot(), file = file.path(supplemental_figure_dir, "FigS2.pdf"), width = 8, height = 6)
graph2ppt(x = ggplot2::last_plot(), file = file.path(supplemental_figure_dir, "FigS2.pptx"), width = 8, height = 6)
# 2. LOW DIMENSIONAL PROJECTION: PLS ON COMBINED Pd + Vv DATA WITH PROJECTION ON 3 Season/Caste, Stage & Species AXES  ####

#### --- packages
library(edgeR)
library(mixOmics)
library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)

#### ---------------------------- HELPERS ------------------------------------

# N13 HOG counts were summarised upstream with tximport::summarizeToGene().
# convert RNAseq counts to logCPM
tmm_logcpm <- function(counts_mat, design_keep = NULL, do_filter = TRUE) {
  y <- DGEList(counts = counts_mat)
  y <- calcNormFactors(y, method = "TMM")
  keep <- rep(TRUE, nrow(y$counts))
  if (do_filter) {
    keep <- if (!is.null(design_keep)) filterByExpr(y, design = design_keep) else rowSums(y$counts) > 0
  }
  list(
    logcpm = cpm(y, log = TRUE, prior.count = 1)[keep, , drop=FALSE],
    dge    = y[keep, , keep.lib.sizes = FALSE]
  )
}

# detect elbow in VIP via kneedle method
vip_elbow_kneedle <- function(vip, 
                              min_keep = 100, 
                              max_keep = round(length(vip)*1) ) {
  v <- sort(as.numeric(vip), decreasing = TRUE)
  v <- v[1:max_keep]
  n <- length(v); if (is.null(max_keep)) max_keep <- n - 5L
  i <- seq_len(n)
  x <- (i - 1) / (n - 1); y <- (v - min(v)) / (max(v) - min(v) + 1e-12)
  num <- abs( (y[n]-y[1]) * x - (1-0) * y + 1*y[1] - y[n]*0 )
  den <- sqrt( (y[n]-y[1])^2 + 1 )
  d <- num/den
  rng <- seq.int(max(2L, min_keep), min(n - 1L, max_keep))
  k <- rng[ which.max(d[rng]) ]
  list(index = k, threshold = v[k])
}

# compute convex hulls
compute_hulls <- function(df, x = "Axis1_SC", y = "Axis2_Stage", group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    filter(n() >= 3) %>%
    dplyr::slice(chull(.data[[x]], .data[[y]])) %>%
    ungroup()
}

## build outcome matrix Y (Stage, Species, Season/Caste; QR orthogonalised to correct for slight data imbalances)

build_targets_Y <- function(meta, scale = FALSE) {
  Stage   <- factor(meta$Stage,   levels = levels(meta$Stage))
  Species <- factor(meta$Species, levels = c("Pd","Vv"))
  SC_raw  <- factor(meta$SC_raw,  levels = c("E","L","W","Q"))
  
  sc_bin <- ifelse((Species=="Pd" & SC_raw=="L") | (Species=="Vv" & SC_raw=="Q"), +1, -1)
  stage_num <- as.numeric(Stage)
  sp_bin <- ifelse(Species=="Vv", +1, -1)
  
  C  <- cbind(SC = sc_bin, Stage = stage_num, Species = sp_bin)
  QR <- qr(C)
  Y  <- qr.Q(QR)
  if (scale) Y <- scale(Y)   # optional scaling/standardisation
  colnames(Y) <- c("SC","Stage","Species") 
  
  # orient correctly
  if (cor(Y[,"SC"], sc_bin)<0)       Y[,"SC"]      <- -Y[,"SC"]
  if (cor(Y[,"Stage"], stage_num)<0) Y[,"Stage"]   <- -Y[,"Stage"]
  if (cor(Y[,"Species"], sp_bin)<0)  Y[,"Species"] <- -Y[,"Species"]
  
  list(Y = Y, sc_bin = sc_bin, stage_num = stage_num, sp_bin = sp_bin)
}


## VIP selection on a given (X,Y)
vip_select <- function(X, Y, min_keep = 100) {
  pls_sc    <- mixOmics::pls(X, matrix(Y[, "SC"],     ncol=1), ncomp=1, scale=FALSE)
  pls_stage <- mixOmics::pls(X, matrix(Y[, "Stage"],  ncol=1), ncomp=1, scale=FALSE)
  pls_sp    <- mixOmics::pls(X, matrix(Y[, "Species"],ncol=1), ncomp=1, scale=FALSE)
  
  vip_sc    <- mixOmics::vip(pls_sc)[,1]
  vip_stage <- mixOmics::vip(pls_stage)[,1]
  vip_sp    <- mixOmics::vip(pls_sp)[,1]
  
  vip_elbow_sc    <- vip_elbow_kneedle(vip_sc, min_keep)
  vip_elbow_stage <- vip_elbow_kneedle(vip_stage, min_keep)
  vip_elbow_sp    <- vip_elbow_kneedle(vip_sp, min_keep)
  
  thr_sc    <- vip_elbow_sc$threshold
  thr_stage <- vip_elbow_stage$threshold
  thr_sp    <- vip_elbow_sp$threshold
  
  sel_sc    <- names(vip_sc   )[vip_sc    >= thr_sc]
  sel_stage <- names(vip_stage)[vip_stage >= thr_stage]
  sel_sp    <- names(vip_sp   )[vip_sp    >= thr_sp]
  
  list(vip = list("SC" = vip_sc, "Stage" = vip_stage, "Species" = vip_sp),
       vip_elbow = list("SC" = vip_elbow_sc, "Stage" = vip_elbow_stage, "Species" = vip_elbow_sp),
       sel = list("SC" = sel_sc, "Stage" = sel_stage, "Species" = sel_sp),
       sel_union = union(union(sel_sc, sel_stage), sel_sp))
}

## Greedy 1–1 mapping of components to {SC, Stage, Species}
greedy_map <- function(C) {
  ys <- colnames(C)[order(apply(abs(C),2,max), decreasing = TRUE)]
  avail <- seq_len(nrow(C))
  m <- integer(ncol(C)); names(m) <- colnames(C)
  for (y in ys) {
    j <- avail[ which.max(abs(C[avail, y])) ]
    m[y] <- j
    avail <- setdiff(avail, j)
  }
  m
}

## Fit on TRAIN, project TEST, orient using TRAIN info (used by LOOCV)
fit_project_once <- function(Xtr, Xte, meta_tr, min_keep = 100, vip_preselect = TRUE) {
  btr <- build_targets_Y(meta_tr)
  Ytr <- btr$Y
  
  if (vip_preselect) {
    fold_vip <- vip_select(Xtr, Ytr, min_keep = min_keep)
    sel <- fold_vip$sel_union
  } else {
    fold_vip <- NULL
    sel <- colnames(Xtr)
  }
  Xtr2 <- Xtr[, sel, drop=FALSE]
  Xte2 <- Xte[, sel, drop=FALSE]
  
  fit   <- mixOmics::pls(Xtr2, Ytr, ncomp = 3, mode = "regression", scale = FALSE)
  SCOtr <- as.data.frame(fit$variates$X)
  
  C <- sapply(1:3, function(j) stats::cor(as.numeric(SCOtr[[j]]), Ytr))
  rownames(C) <- paste0("comp", 1:3); colnames(C) <- colnames(Ytr)
  map <- greedy_map(C)
  
  Axis_SC_tr    <- as.numeric(SCOtr[[ map["SC"]    ]])
  Axis_Stage_tr <- as.numeric(SCOtr[[ map["Stage"] ]])
  Axis_Sp_tr    <- as.numeric(SCOtr[[ map["Species"]]])
  
  f1 <- sign(stats::cor(Axis_SC_tr,    btr$sc_bin));    if (is.na(f1)) f1 <- 1
  f2 <- sign(stats::cor(Axis_Stage_tr, btr$stage_num)); if (is.na(f2)) f2 <- 1
  f3 <- sign(stats::cor(Axis_Sp_tr,    btr$sp_bin));    if (is.na(f3)) f3 <- 1
  
  pr   <- predict(fit, newdata = Xte2)
  SCOt <- as.data.frame(pr$variates)  # robust for 1-row predictions
  
  list(
    Axis1_SC      = f1 * as.numeric(SCOt[[ map["SC"]    ]]),
    Axis2_Stage   = f2 * as.numeric(SCOt[[ map["Stage"] ]]),
    Axis3_Species = f3 * as.numeric(SCOt[[ map["Species"]]]),
    selected_features = sel,
    selected_by_target = if (is.null(fold_vip)) {
      c(SC = length(sel), Stage = length(sel), Species = length(sel))
    } else {
      vapply(fold_vip$sel, length, integer(1))
    }
  )
}

## Select the TMM reference using training libraries only.
choose_tmm_reference <- function(counts, train_columns) {
  train_counts <- counts[, train_columns, drop = FALSE]
  lib_size <- colSums(train_counts)
  f75 <- apply(train_counts, 2, stats::quantile, probs = 0.75) / lib_size
  if (stats::median(f75) < 1e-20) {
    train_columns[which.max(colSums(sqrt(train_counts)))]
  } else {
    train_columns[which.min(abs(f75 - mean(f75)))]
  }
}

## TMM-normalise all libraries relative to a training-derived reference. The
## factor scaling is also based only on training libraries. The held-out sample
## contributes only its own counts when its normalisation factor is estimated.
tmm_logcpm_from_training <- function(counts, train_columns) {
  ref_sample <- choose_tmm_reference(counts, train_columns)
  ref_index <- match(ref_sample, colnames(counts))
  train_index <- match(train_columns, colnames(counts))
  norm_factors <- edgeR::calcNormFactors(
    counts, method = "TMM", refColumn = ref_index
  )
  norm_factors <- norm_factors /
    exp(mean(log(norm_factors[train_index])))
  effective_lib_size <- colSums(counts) * norm_factors
  list(
    logcpm = edgeR::cpm(
      counts, lib.size = effective_lib_size,
      log = TRUE, prior.count = 1
    ),
    ref_sample = ref_sample
  )
}

## Fully nested leave-one-sample-out cross-validation. Expression filtering,
## TMM reference selection, HOG-wise centring/scaling, VIP/Kneedle selection,
## PLS fitting, component assignment and orientation use training samples only.
loo_pls_emm_nested <- function(pd_counts, vv_counts, metadata_pd, metadata_vv,
                               meta, min_keep = 100,
                               vip_preselect = TRUE) {
  stopifnot(
    identical(rownames(pd_counts), rownames(vv_counts)),
    nrow(meta) == ncol(pd_counts) + ncol(vv_counts)
  )
  n <- nrow(meta)
  
  meta <- meta %>%
    dplyr::mutate(
      Species = factor(Species, levels = c("Pd","Vv")),
      Stage   = factor(Stage,   levels = levels(Stage)),
      SC2 = factor(ifelse((Species=="Pd" & SC_raw=="L") |
                            (Species=="Vv" & SC_raw=="Q"), "L/Q", "E/W"),
                   levels = c("E/W","L/Q"))
    )
  
  out <- vector("list", n)
  diagnostics <- vector("list", n)
  for (i in seq_len(n)) {
    te <- i; tr <- setdiff(seq_len(n), te)
    meta_tr <- meta[tr, , drop = FALSE]
    pd_train <- intersect(colnames(pd_counts), meta_tr$sample)
    vv_train <- intersect(colnames(vv_counts), meta_tr$sample)

    pd_metadata_tr <- metadata_pd[
      match(pd_train, metadata_pd$sampleNumber), , drop = FALSE
    ]
    pd_design_tr <- stats::model.matrix(
      ~ lifeStage + larvalStage, data = pd_metadata_tr
    )
    pd_dge_tr <- edgeR::DGEList(
      counts = pd_counts[, pd_train, drop = FALSE]
    )
    pd_dge_tr <- edgeR::calcNormFactors(pd_dge_tr, method = "TMM")
    keep <- edgeR::filterByExpr(pd_dge_tr, design = pd_design_tr)
    fold_features <- rownames(pd_counts)[keep]

    pd_norm <- tmm_logcpm_from_training(
      pd_counts[fold_features, , drop = FALSE], pd_train
    )
    vv_norm <- tmm_logcpm_from_training(
      vv_counts[fold_features, , drop = FALSE], vv_train
    )
    X0 <- t(cbind(
      pd_norm$logcpm[, colnames(pd_counts), drop = FALSE],
      vv_norm$logcpm[, colnames(vv_counts), drop = FALSE]
    ))
    X0 <- X0[meta$sample, , drop = FALSE]

    fold_means <- colMeans(X0[tr, , drop = FALSE])
    fold_sds <- apply(X0[tr, , drop = FALSE], 2, stats::sd)
    nonconstant <- is.finite(fold_sds) & fold_sds > 0
    X0 <- X0[, nonconstant, drop = FALSE]
    fold_means <- fold_means[nonconstant]
    fold_sds <- fold_sds[nonconstant]
    X <- sweep(X0, 2, fold_means, "-")
    X <- sweep(X, 2, fold_sds, "/")

    pr <- fit_project_once(
      Xtr = X[tr, , drop = FALSE],
      Xte = X[te, , drop = FALSE],
      meta_tr = meta_tr,
      min_keep = min_keep,
      vip_preselect = vip_preselect
    )
    out[[i]] <- tibble(
      sample         = meta$sample[te],
      Axis1_SC       = pr$Axis1_SC,
      Axis2_Stage    = pr$Axis2_Stage,
      Axis3_Species  = pr$Axis3_Species
    )
    diagnostics[[i]] <- tibble(
      sample = meta$sample[te],
      candidate_HOGs = ncol(X),
      selected_HOGs = length(pr$selected_features),
      selected_SC = unname(pr$selected_by_target["SC"]),
      selected_Stage = unname(pr$selected_by_target["Stage"]),
      selected_Species = unname(pr$selected_by_target["Species"]),
      pd_TMM_reference = pd_norm$ref_sample,
      vv_TMM_reference = vv_norm$ref_sample
    )
  }
  list(
    scores = dplyr::bind_rows(out) %>%
      dplyr::left_join(
        meta[, c("sample", "Stage", "Species", "SC2")], by = "sample"
      ),
    fold_diagnostics = dplyr::bind_rows(diagnostics)
  )
}


## 2.1 PLS FIT ####

### 1) Use direct N13 HOG count matrices
common_og <- intersect(rownames(cts_pd), rownames(cts_vv))
pd_og <- cts_pd[common_og, metadata_pd$sampleNumber, drop = FALSE]
vv_og <- cts_vv[common_og, metadata_vv$sampleNumber, drop = FALSE]
storage.mode(pd_og) <- "double"
storage.mode(vv_og) <- "double"

### 2) Calculate TMM + logCPM
pd_design_keep <- model.matrix(~ lifeStage + larvalStage, data = metadata_pd)
pd_norm <- tmm_logcpm(pd_og, design_keep = pd_design_keep, do_filter = TRUE)
vv_norm <- tmm_logcpm(vv_og, design_keep = NULL, do_filter = FALSE)

features <- intersect(rownames(pd_norm$logcpm), rownames(vv_norm$logcpm))
X_pd0 <- pd_norm$logcpm[features, metadata_pd$sampleNumber, drop=FALSE]
X_vv0 <- vv_norm$logcpm[features, metadata_vv$sampleNumber, drop=FALSE]

### 3) Global standardisation (preserves Species signal)
X_all0 <- cbind(X_pd0, X_vv0)                    # genes × samples
mu_all <- rowMeans(X_all0)
sd_all <- apply(X_all0, 1, sd); sd_all[sd_all == 0] <- 1
Z_all  <- t( (X_all0 - mu_all) / sd_all )        # samples × genes

### Metadata (combined)
meta_pd <- tibble(
  sample  = metadata_pd$sampleNumber,
  Species = factor("Pd", levels = c("Pd","Vv")),
  Stage   = factor(metadata_pd$larvalStage, levels = stages),
  SC_raw  = factor(metadata_pd$lifeStage, levels = c("E","L"))
)
meta_vv <- tibble(
  sample  = metadata_vv$sampleNumber,
  Species = factor("Vv", levels = c("Pd","Vv")),
  Stage   = factor(metadata_vv$larvalStage, levels = stages),
  SC_raw  = factor(metadata_vv$lifeStage, levels = c("W","Q"))
)
meta_all <- bind_rows(meta_pd, meta_vv)
stopifnot(nrow(Z_all) == nrow(meta_all))
rownames(Z_all) <- meta_all$sample

### 4) Construct target matrix Y (SC, Stage, Species), QR orthonormalised
bt  <- build_targets_Y(meta_all, scale = FALSE)
Y   <- bt$Y
sc_bin   <- bt$sc_bin
stage_num <- bt$stage_num
sp_bin    <- bt$sp_bin

### quick sanity checks
cor(Y[, "SC"],      sc_bin)      # 1
cor(Y[, "Stage"],   stage_num)   # 0.9963794
cor(Y[, "Species"], sp_bin)      # 0.9995007

### 5) Optional VIP-based feature preselection

vip_preselect = TRUE

# Optional VIP-based feature preselection: 
if (vip_preselect) {
     vipsel   <- vip_select(Z_all, Y, min_keep = 100)
     vipsel$vip_elbow$Species$index # 2210
     vipsel$vip_elbow$Stage$index # 1205
     vipsel$vip_elbow$SC$index # 909
     length(vipsel$sel_union) # 4750 in the final saved analysis
     sel <- vipsel$sel_union } else { sel <- colnames(Z_all) }

Z_sub <- Z_all[, sel, drop = FALSE]


### 6) 3-component PLS on combined data
plsfit <- mixOmics::pls(X = Z_sub, Y = Y, ncomp = 3, mode = "regression", scale = FALSE) # also possible to use scale=T, now I use raw QR orthogonalised responses

### Variance captured in X by 3 latent axes
plsfit$prop_expl_var$X
# comp1      comp2      comp3 
# 0.35939164 0.22769332 0.07815831

# -> A 3-component PLS model (mode=regression) was fitted to the globally 
# standardised expression matrix X and three orthogonal contrasts in Y: 
# Season/Caste, Stage and Species. These three latent variables summarised 
# 7.8%, 22.8% and 36.0% of the variation in X (cumulative 67.2%).



### Map components to targets and orient axes
SCO  <- as.data.frame(plsfit$variates$X)
cors <- sapply(1:3, function(j) cor(SCO[[j]], Y))
map  <- apply(abs(cors), 2, which.max); names(map) <- c("SC","Stage","Species")

Axis_SC_raw    <- SCO[[ map["SC"]    ]]
Axis_Stage_raw <- SCO[[ map["Stage"] ]]
Axis_Sp_raw    <- SCO[[ map["Species"]]]

flip1 <- sign(cor(Axis_SC_raw, sc_bin))
flip2 <- sign(cor(Axis_Stage_raw, stage_num))
flip3 <- sign(cor(Axis_Sp_raw, sp_bin))

Axis1_SC    <- flip1 * Axis_SC_raw
Axis2_Stage <- flip2 * Axis_Stage_raw
Axis3_Sp    <- flip3 * Axis_Sp_raw

# get features contributing most to each PLS axis
vips <- mixOmics::vip(plsfit)
vip_sc    <- sort(vips[,map["SC"]], decreasing = TRUE)
vip_stage <- sort(vips[,map["Stage"]], decreasing = TRUE)
vip_sp    <- sort(vips[,map["Species"]], decreasing = TRUE)

vip_elbow_sc    <- vip_elbow_kneedle(vip_sc, min_keep = 100)
vip_elbow_stage <- vip_elbow_kneedle(vip_stage, min_keep = 100)
vip_elbow_sp    <- vip_elbow_kneedle(vip_sp, min_keep = 100)

thr_sc    <- vip_elbow_sc$threshold
thr_stage <- vip_elbow_stage$threshold
thr_sp    <- vip_elbow_sp$threshold

sel_sc    <- names(vip_sc   )[vip_sc    >= thr_sc]
sel_stage <- names(vip_stage)[vip_stage >= thr_stage]
sel_sp    <- names(vip_sp   )[vip_sp    >= thr_sp]
length(sel_sc)    # 246
length(sel_stage) # 3565
length(sel_sp)    # 2190


### 7) Calculate and plot raw PLS scores
scores_all <- meta_all %>%
  mutate(Axis1_SC = Axis1_SC,
         Axis2_Stage = Axis2_Stage,
         Axis3_Species = Axis3_Sp,
         SC2 = factor(ifelse(sc_bin==+1, "L/Q", "E/W"), levels = c("E/W","L/Q")))


# plot Season/Caste x Stage PLS scores
scores_all2 <- scores_all
scores_all2$Species <- factor(scores_all2$Species,
                              levels = c("Pd", "Vv"),
                              labels = c("Polistes dominula", "Vespula vulgaris"))
scores_all2$SC2 <- factor(scores_all2$SC2,
                          levels = c("E/W", "L/Q"),
                          labels = c("Early / Worker", "Late / Queen"))

hull_all <- compute_hulls(scores_all, group_vars = c("Species","Stage","SC2"))

hull_all2 <- hull_all
hull_all2$Species <- factor(hull_all2$Species,
                            levels = c("Pd", "Vv"),
                            labels = c("Polistes dominula", "Vespula vulgaris"))
hull_all2$SC2 <- factor(hull_all2$SC2,
                          levels = c("E/W", "L/Q"),
                          labels = c("Early / Worker", "Late / Queen"))

stage_palette <- setNames(
  rev(RColorBrewer::brewer.pal(6, "Spectral")),
  levels(scores_all2$Stage)
)
darken_colour <- function(colour, proportion = 0.72) {
  rgb_values <- grDevices::col2rgb(colour) * proportion
  grDevices::rgb(
    rgb_values[1, ], rgb_values[2, ], rgb_values[3, ],
    maxColorValue = 255
  )
}
arrow_palette <- setNames(
  darken_colour(stage_palette, proportion = 0.72),
  names(stage_palette)
)

centroid_arrows <- scores_all2 %>%
  group_by(Species, Stage, SC2) %>%
  summarise(
    centroid_x = mean(Axis1_SC),
    centroid_y = mean(Axis2_Stage),
    sample_n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = SC2,
    values_from = c(centroid_x, centroid_y, sample_n),
    names_sep = "__"
  ) %>%
  transmute(
    Species,
    Stage,
    x = `centroid_x__Early / Worker`,
    y = `centroid_y__Early / Worker`,
    xend = `centroid_x__Late / Queen`,
    yend = `centroid_y__Late / Queen`,
    source_n = `sample_n__Early / Worker`,
    target_n = `sample_n__Late / Queen`,
    displacement = sqrt((xend - x)^2 + (yend - y)^2),
    arrow_colour = unname(arrow_palette[as.character(Stage)])
  )

if (nrow(centroid_arrows) != 12L || any(!is.finite(unlist(
  centroid_arrows[c("x", "y", "xend", "yend", "displacement")]
)))) {
  stop("Expected complete centroid arrows for six stages in both species.")
}
write_tsv(
  centroid_arrows %>% dplyr::select(-arrow_colour),
  file.path(output_dir, "Fig1A_centroid_arrow_coordinates.tsv")
)

arrow_spec_halo <- grid::arrow(
  angle = 19, length = grid::unit(0.155, "inches"), type = "closed"
)
arrow_spec_colour <- grid::arrow(
  angle = 19, length = grid::unit(0.155, "inches"), type = "closed"
)

p_all <- ggplot(scores_all2,
                aes(Axis1_SC, Axis2_Stage, color = Stage, shape = SC2)) +
  facet_wrap(~ Species, ncol = 1) +
  geom_polygon(data = hull_all2,
               aes(group = interaction(Species, Stage, SC2),
                   fill = Stage, color = Stage, alpha = SC2),
               linewidth = 0.25, inherit.aes = TRUE) +
  geom_point(size = 2) +
  geom_segment(
    data = centroid_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 2.20,
    lineend = "round",
    arrow = arrow_spec_halo
  ) +
  geom_segment(
    data = centroid_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = centroid_arrows$arrow_colour,
    linewidth = 1.28,
    lineend = "round",
    arrow = arrow_spec_colour
  ) +
  scale_shape_manual("Season / Caste", values = c(17, 16)) +
  scale_alpha_manual("Season / Caste", values = c(0.22, 0.60), guide = "none") +
  scale_color_manual(values = stage_palette, guide = "none") +
  scale_fill_manual(values = stage_palette, guide = "none") +
  xlab("Axis 1: Season / Caste") +
  ylab("Axis 2: Stage") +
  # labs(title = "Combined PLS — Pd + Vv") +
  theme_few() +
  # With rvg, 0.47 ggplot linewidth units exports as approximately 1.00 pt.
  theme(panel.border = element_rect(fill = NA, linewidth = 0.47),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.47),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 10),
        aspect.ratio = 0.93,
        panel.spacing = grid::unit(0.18, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=9))
print(p_all)

## sanity checks
cor(Axis1_SC,   Y[, "SC"])       # 0.5763527
cor(Axis2_Stage,Y[, "Stage"])    # 0.8513001
cor(Axis3_Sp,   Y[, "Species"])  # 0.9836731


### 8) Export PLS scores, oriented loadings, VIPs & annotate orthologs

## helper: sign-specific ranks (1 = most positive; 1 = most negative)
rank_signed <- function(x) {
  rpos <- rep(NA_integer_, length(x))
  rneg <- rep(NA_integer_, length(x))
  if (any(x > 0)) rpos[x > 0] <- rank(-x[x > 0], ties.method = "min")  # larger positive → smaller rank
  if (any(x < 0)) rneg[x < 0] <- rank( x[x < 0], ties.method = "min")  # more negative (smaller) → rank 1
  list(pos = rpos, neg = rneg)
}


## a) Sample factor scores (already oriented)
scores_export <- scores_all %>%
  dplyr::select(sample, Species, Stage, SC2,
                Axis1_SC, Axis2_Stage, Axis3_Species)

## b) Oriented X-loadings per biological axis
# loadings for variables (orthogroups) × components
Lx <- as.matrix(plsfit$loadings$X)
# keep only features present in the fitted subset (useful if sPLS was used)
feat <- intersect(colnames(Z_sub), rownames(Lx))
Lx   <- Lx[feat, , drop = FALSE]

loadings_df <- tibble::tibble(
  orthogroup       = rownames(Lx),
  loading_SC       =  flip1 * as.numeric(Lx[, map["SC"]]),
  loading_Stage    =  flip2 * as.numeric(Lx[, map["Stage"]]),
  loading_Species  =  flip3 * as.numeric(Lx[, map["Species"]])
)

## indicate top features & rank by loading (separately per sign)
loadings_df <- loadings_df %>%
  dplyr::mutate(
    top_SC      = orthogroup %in% sel_sc,
    top_SC_pos      = top_SC      & is.finite(loading_SC)      & loading_SC      > 0,
    top_SC_neg      = top_SC      & is.finite(loading_SC)      & loading_SC      < 0,
    
    top_Stage   = orthogroup %in% sel_stage,
    top_Stage_pos   = top_Stage   & is.finite(loading_Stage)   & loading_Stage   > 0,
    top_Stage_neg   = top_Stage   & is.finite(loading_Stage)   & loading_Stage   < 0,
    
    top_Species = orthogroup %in% sel_sp,
    top_Species_pos = top_Species & is.finite(loading_Species) & loading_Species > 0,
    top_Species_neg = top_Species & is.finite(loading_Species) & loading_Species < 0
  ) %>%
  
  ## sign-specific ranks based on factor loadings
  dplyr::mutate(
    rank_SC_loading_pos      = rank_signed(loading_SC)$pos,
    rank_SC_loading_neg      = rank_signed(loading_SC)$neg,
    rank_Stage_loading_pos   = rank_signed(loading_Stage)$pos,
    rank_Stage_loading_neg   = rank_signed(loading_Stage)$neg,
    rank_Species_loading_pos = rank_signed(loading_Species)$pos,
    rank_Species_loading_neg = rank_signed(loading_Species)$neg
  )


## c) VIP dataframe per biological axis (VIP is nonnegative, orientation-invariant)
vip_mat <- vips[feat, , drop = FALSE]

vip_df <- tibble::tibble(
  orthogroup  = rownames(vip_mat),
  VIP_SC      = as.numeric(vip_mat[, map["SC"]]),
  VIP_Stage   = as.numeric(vip_mat[, map["Stage"]]),
  VIP_Species = as.numeric(vip_mat[, map["Species"]])
)

## d) Annotate orthologs table (left join preserves all rows)
orthologs_pls <- orthologs %>%
  dplyr::left_join(loadings_df, by = "orthogroup") %>%
  dplyr::left_join(vip_df,      by = "orthogroup")


## e) Export PLS scores & loadings to file
write.csv(scores_export,  file = file.path(output_dir, "pls_scores_by_sample.csv"), row.names = FALSE)
write.csv(orthologs_pls,  file = file.path(output_dir, "pls_N13_HOGs_with_loadings_vip.csv"), row.names = FALSE)



## 2.2 CALCULATE LOOCV PERFORMANCE ####

# Calculate fully nested leave-one-out sample scores from the unfiltered shared
# HOG count matrices. No filtering, scaling or feature-selection information
# from the held-out sample is used to fit its fold-specific PLS model.
oos_nested <- loo_pls_emm_nested(
  pd_counts = pd_og,
  vv_counts = vv_og,
  metadata_pd = metadata_pd,
  metadata_vv = metadata_vv,
  meta = meta_all,
  min_keep = 100,
  vip_preselect = vip_preselect
)
oos <- oos_nested$scores
pls_loocv_fold_diagnostics <- oos_nested$fold_diagnostics

# calculate significance of emmeans contrasts for season/caste, stage & species based on leave-one-out sample scores
# TODO change lm_robust to MASS::lmrob & use marginaleffects package for significance
# (not supported by emmeans) to be consistent with robust models later on?

afex::set_sum_contrasts()
library(estimatr)
fit_axis1_oos <- estimatr::lm_robust(Axis1_SC ~ (SC2 + Stage + Species)^2, data = oos) 
sc_contr_oos <- emmeans::contrast(
  emmeans::emmeans(fit_axis1_oos, ~ SC2 | Stage * Species),
  method = "revpairwise",
  adjust = "none"
)
sc_contr_oos_result <- summary(
  sc_contr_oos,
  side = ">",
  adjust = "none"
)

fit_axis2_oos <- estimatr::lm_robust(Axis2_Stage ~ (SC2 + Stage + Species)^2, data = oos)
st_contr_oos  <- emmeans::contrast(emmeans::emmeans(fit_axis2_oos, ~ Stage | Species),
                                   method = "consec")

fit_axis3_oos <- estimatr::lm_robust(Axis3_Species ~ (SC2 + Stage + Species)^2, data = oos)
sp_contr_oos  <- emmeans::contrast(emmeans::emmeans(fit_axis3_oos, ~ Species | Stage),
                                   method = "revpairwise")

# Season/Caste contrasts per stage & species
write.csv(
  data.frame(sc_contr_oos_result),
  file = file.path(output_dir, "pls_season caste emmeans contrasts LOOCV scores.csv"),
  row.names = FALSE
)
sc_contr_oos_result
# Planned right-tailed contrasts test (L/Q) - (E/W) > 0 without multiplicity
# adjustment. Exact current results are exported above and in TableS5.tsv.

# stage contrasts per species
write.csv(data.frame(st_contr_oos),  file = file.path(output_dir, "pls_stage emmeans contrasts LOOCV scores.csv"), row.names = FALSE)
st_contr_oos

# -> biggest shift from L4 to L5 and L5 to P in Polistes dominula (late in development),
# but in L2 to L3, L3 to L4 and L4 to L5 in Vespula vulgaris (shifted earlier in development)

# species contrasts per stage
write.csv(data.frame(sp_contr_oos),  file = file.path(output_dir, "pls_species emmeans contrasts LOOCV scores.csv"), row.names = FALSE)
readr::write_tsv(
  pls_loocv_fold_diagnostics,
  file.path(output_dir, "pls_fully_nested_LOOCV_fold_diagnostics.tsv")
)
sp_contr_oos

# -> significant shifts by species across all stages




# 3. N13 HOG-LEVEL GO ENRICHMENT OF FEATURES ASSOCIATED WITH THE PLS AXES ####

# GO ontology
library(ontologyIndex) 
go <- get_ontology(f_obo, extract_tags = "everything")

# Final GO annotation: direct EXCON/Galaxy EggNOG and InterProScan GO retained
# through the Polistes-Vespula conserved-HOG rule, plus experimentally supported
# FlyBase GO transferred through the N13 Drosophila orthologue.
go_orthologs <- read_tsv(go_annot_orthologs, show_col_types = FALSE)

## HELPERS ---------------------------------------------------------------
build_gene2go <- function(go_tbl, id_col,
                          keep_aspect = NULL) {
  if (length(id_col) != 1L || !id_col %in% names(go_tbl)) {
    stop("id_col must name one column in the GO annotation table")
  }
  g <- go_tbl %>%
    mutate(
      go_terms = str_replace_all(go_terms, "\\s+", ""),
      go_terms = strsplit(go_terms, "[,;]+")
    ) %>%
    unnest(go_terms) %>%
    filter(str_detect(go_terms, "^GO:"))
  
  if (!is.null(keep_aspect)) {
    g <- g %>%
      mutate(aspect = case_when(
        aspect %in% c("biological_process","BP") ~ "BP",
        aspect %in% c("molecular_function","MF") ~ "MF",
        aspect %in% c("cellular_component","CC") ~ "CC",
        TRUE ~ NA_character_
      )) %>% filter(aspect == keep_aspect)
  }
  
  g <- g %>% distinct(.data[[id_col]], go_terms)
  gene2go  <- split(g$go_terms, g[[id_col]])
  universe <- sort(unique(g[[id_col]]))
  list(universe = universe, gene2go = gene2go)
}

run_topgo <- function(gene_ids, universe, gene2go,
                      ontology = c("BP","MF","CC"),
                      algorithm = c("elim","weight01","classic"),
                      statistic = c("fisher","ks","t"),
                      nodeSize = 20,
                      adjust = c("none", "BH")) {
  ontology  <- match.arg(ontology)
  algorithm <- match.arg(algorithm)
  statistic <- match.arg(statistic)
  adjust <- match.arg(adjust)

  allGenes <- factor(as.integer(universe %in% gene_ids))
  names(allGenes) <- universe

  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = allGenes,
                geneSel  = function(x) x == 1,
                annot    = annFUN.gene2GO,
                gene2GO  = gene2go,
                nodeSize = nodeSize)

  res <- runTest(GOdata, algorithm = algorithm, statistic = statistic)

  tab <- GenTable(GOdata, p = res, topNodes = length(score(res)), numChar = 10000) %>% # use all nodes
    as_tibble() %>%
    mutate(FoldEnrichment = Significant / pmax(Expected, 1e-12),
           p = suppressWarnings(as.numeric(p)))
  if (adjust == "BH") tab$p_BH <- p.adjust(tab$p, method = "BH")
  tab <- dplyr::arrange(tab, dplyr::across(dplyr::any_of(c("p_BH","p")))) %>%
    relocate(FoldEnrichment, .before = 6)
  tab
}

# wrap run_topgo() so it adds members
run_topgo_with_members <- function(gene_ids, universe, gene2go,
                                   ontology = c("BP","MF","CC"),
                                   algorithm = c("weight01","elim","classic"),
                                   statistic = c("fisher","ks","t"),
                                   nodeSize = 20,
                                   adjust = c("none","BH")) {
  ontology  <- match.arg(ontology)
  algorithm <- match.arg(algorithm)
  statistic <- match.arg(statistic)
  adjust    <- match.arg(adjust)
  
  fg <- intersect(gene_ids, universe)
  if (length(fg) < 2L) return(tibble())
  
  allGenes <- factor(as.integer(universe %in% fg)); names(allGenes) <- universe
  
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = allGenes,
                geneSel  = function(x) x == 1,
                annot    = annFUN.gene2GO,
                gene2GO  = gene2go,
                nodeSize = nodeSize)
  
  res <- runTest(GOdata, algorithm = algorithm, statistic = statistic)
  
  tab <- GenTable(GOdata, p = res, topNodes = length(score(res)), numChar = 10000) |>
    as_tibble() |>
    mutate(
      p = suppressWarnings(as.numeric(p)),
      Expected = as.numeric(Expected),
      Significant = as.integer(Significant),
      FoldEnrichment = Significant / pmax(Expected, 1e-12),
      ontology = ontology,
      algorithm = algorithm
    )
  
  if (adjust == "BH") tab <- tab %>% mutate(p_BH = p.adjust(p, "BH"))
  
  ## --- Members from topGO graph (guarantees Members_n == Significant) ---
  sig_genes <- names(allGenes)[allGenes==1]
  term_genes_list <- genesInTerm(GOdata, tab$GO.ID)  # list: GO -> universe genes in term (with propagation)
  mems <- lapply(term_genes_list, function(gs) intersect(gs, sig_genes))
  
  tab %>%
    mutate(
      Members    = mems,                       # list-col of orthogroups
      Members_n  = lengths(Members),
      MembersStr = vapply(Members, paste, collapse=";", FUN.VALUE=character(1))
    ) %>%
    relocate(FoldEnrichment, .before = p) %>%
    arrange(p, desc(FoldEnrichment))
}

# helpers to do GO enrichment analysis along PLS axes
# Build OG-level GO map (once)
og_map_pls <- build_gene2go(go_orthologs, id_col = "HOG")

# Universe = orthogroups that actually entered the PLS model (and have GO)
univ_pls <- intersect(colnames(Z_sub), og_map_pls$universe)

# Restrict the gene2go list to the universe to speed up topGO
g2g_pls   <- og_map_pls$gene2go[names(og_map_pls$gene2go) %in% univ_pls]

# small helper
.run_pls_go <- function(foreground, universe = univ_pls, gene2go = g2g_pls,
                        algorithm = "weight01", 
                        statistic = "fisher", 
                        nodeSize = 20,
                        adjust = "none") {
  fg <- intersect(foreground, universe)
  if (length(fg) < 2L) {
    return(list(BP = tibble(), MF = tibble(), CC = tibble()))
  }
  list(
    BP = run_topgo_with_members(fg, universe, gene2go, ontology = "BP",
                   algorithm = algorithm, statistic = statistic,
                   nodeSize = nodeSize, adjust = adjust),
    MF = run_topgo_with_members(fg, universe, gene2go, ontology = "MF",
                   algorithm = algorithm, statistic = statistic,
                   nodeSize = nodeSize, adjust = adjust),
    CC = run_topgo_with_members(fg, universe, gene2go, ontology = "CC",
                   algorithm = algorithm, statistic = statistic,
                   nodeSize = nodeSize, adjust = adjust)
  )
}

## Helpers to build a flat gene × term table of p-values across analyses
## GO name/definition lookup from go_orthologs tibble
go_lookup <- go_orthologs %>%
  mutate(go_terms = gsub("\\s+", "", go_terms)) %>%
  dplyr::filter(grepl("^GO:", go_terms)) %>%
  group_by(go_terms) %>%
  summarise(
    term_name  = dplyr::first(na.omit(term_name)),
    definition = dplyr::first(na.omit(definition)),
    .groups = "drop"
  ) %>% dplyr::rename(GO.ID = go_terms)

## Flatten one ontology table; add a constant label column
flatten_one_table <- function(tbl, label) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tibble())
  tbl %>%
    mutate(label = label) %>%                              # <-- add label
    left_join(go_lookup, by = "GO.ID") %>%
    dplyr::select(label, ontology, GO.ID, Term, term_name, definition,
           FoldEnrichment, p, dplyr::any_of("p_BH"),
           Significant, Expected, Members) %>%
    unnest_longer(Members, values_to = "orthogroup") %>%
    dplyr::filter(!is.na(orthogroup))
}

## Flatten a named list with $BP/$MF/$CC
flatten_go_results <- function(named_res_list, label) {
  bind_rows(
    flatten_one_table(named_res_list$BP, label),
    flatten_one_table(named_res_list$MF, label),
    flatten_one_table(named_res_list$CC, label)
  )
}




## 3.0 ORTHOGROUP-LEVEL ENRICHMENT TESTS OF ORTHOGROUPS MOST PREDICTIVE OF SEASON/CASTE, STAGE OR SPECIES IN PLS ANALYSIS #######################

# define relevant foreground sets
# top genes predictive of each axis (season/caste, stage or species) based on loadings PLS analysis & VIP scores+kneedle elbow method
sel_sc_pos <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$top_SC_pos==TRUE]))
length(sel_sc_pos) # 129
sel_sc_neg <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$top_SC_neg==TRUE]))
length(sel_sc_neg) # 117
sel_st_pos <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$top_Stage_pos==TRUE]))
length(sel_st_pos) # 1802
sel_st_neg <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$top_Stage_neg==TRUE]))
length(sel_st_neg) # 1763
sel_sp_pos <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$top_Species_pos==TRUE]))
length(sel_sp_pos) # 1210
sel_sp_neg <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$top_Species_neg==TRUE]))
length(sel_sp_neg) # 980

# Final foreground definition: top 300 HOGs in each loading direction.
K <- 300
sel_sc_pos_topK <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$rank_SC_loading_pos<=K]))
sel_sc_neg_topK <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$rank_SC_loading_neg<=K]))
sel_st_pos_topK <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$rank_Stage_loading_pos<=K]))
sel_st_neg_topK <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$rank_Stage_loading_neg<=K]))
sel_sp_pos_topK <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$rank_Species_loading_pos<=K]))
sel_sp_neg_topK <- as.vector(na.omit(orthologs_pls$orthogroup[orthologs_pls$rank_Species_loading_neg<=K]))


### 3.0.1 ENRICHMENT TEST FOR ORTHOGROUPS MOST PREDICTIVE OF LATE/QUEEN SEASON/CASTE ####
pls_SC_pos_go <- .run_pls_go(sel_sc_pos_topK)
pls_SC_pos_bp_enriched <- pls_SC_pos_go$BP %>%
  filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5) %>%
  dplyr::select(
    GO.ID, Term, FoldEnrichment, p, Significant, Expected,
    Members_n, MembersStr
  )
print(pls_SC_pos_bp_enriched)

### Dot plot
top <- 50L

## 1. Patterns 
## GO terms to drop (vertebrate-specific)
vertebrate_drop <- c("neutrophil|mammary gland|placenta|erythrocyte",
                     "connective tissue|glial|cerebr|pregnancy")

## Regexes for categories
# categories
cat1 = "Differential feeding / nutrient metabolism"
cat2 = "Reproduction / caste-related development"
cat3 = "Hibernation / diapause / long-term survival"
cat4 = "Other"

# Feeding / nutrient metabolism
pat_feed <- paste0(
  # nutrient / endocrine sensing (feeding angle)
  "insulin|glucose homeostasis|regulation of glucose|",
  
  # central carbon + energy metabolism
  "gluconeogen|carbohydrate|monosaccharide|",
  "carboxylic acid|monocarboxylic acid|dicarboxylic acid|citrate|",
  "generation of precursor metabolites|energy derivation|",
  "regulation of ATP biosynthetic|mitochondrial membrane potential|",
  
  # lipids/fatty acids/triglycerides (catabolism + metabolism)
  "fatty acid|triglyceride|lipid catabolic|lipid metabolic|",
  "fatty acid derivative|",
  
  # amino acids
  "amino acid|aspartate family|glutamine family|",
  
  # vitamins/cofactors/nucleotides
  "\\bNAD\\b|nicotinamide|vitamin|",
  "nucleobase|ribonucleoside|ribonucleotide|pyrimidine|purine|",
  
  # diet-related detox
  "xenobiotic|antibiotic|aldehyde"
)

# Reproduction / caste-related development
pat_repr <- paste0(
  "reproductive|gonad|follic|",
  "ovar|oocyt|oogen|",
  "epidermal growth factor|",
  "response to hormone"
)

# Overwintering / diapause / long-term survival
pat_hib <- paste0(
  "lifespan|hibernat|diapause|starvation|",
  "hypoxia|oxidative stress|cellular response to heat|",
  "response to activity|", # or could be left as other
  "negative regulation of apoptotic process|",  # or could be left as other
  "ERAD pathway|", # proteostasis/ER stress resistance, or could be left as other
  
  # reserve build-up & mobilisation framing
  "lipid storage|regulation of lipid storage|lipid local|adipo",
  
  # broader survival / maintenance cues
  "negative regulation of apoptotic process|",
  "protein stabilization"
)

## 2. Colour mapping for categories 

category_cols <- c("#164BCE", "#E02020", "#159230", "#111111")
names(category_cols) <- c(cat1, cat2, cat3, cat4)

## 3. Build the final display table.
tbl_eligible <- pls_SC_pos_go$BP %>%
  filter(
    p <= 0.02,
    Significant >= 3,
    FoldEnrichment >= 1.5,
    !str_detect(
      Term,
      regex(paste(vertebrate_drop, collapse = "|"), ignore_case = TRUE)
    )
  ) %>%
  dplyr::select(
    GO.ID, Term, FoldEnrichment, p,
    Significant, Expected, Members_n, MembersStr
  ) %>%
  mutate(
    Category = case_when(
      str_detect(Term, regex(pat_repr, ignore_case = TRUE)) ~ cat2,
      str_detect(Term, regex(pat_hib, ignore_case = TRUE)) ~ cat3,
      str_detect(Term, regex(pat_feed, ignore_case = TRUE)) ~ cat1,
      TRUE ~ cat4
    ),
    Category = factor(Category, levels = c(cat1, cat3, cat2, cat4))
  )

tbl <- bind_rows(
  tbl_eligible %>% slice_min(p, n = top, with_ties = FALSE),
  tbl_eligible %>% filter(Category == cat2)
) %>%
  distinct(GO.ID, .keep_all = TRUE) %>%
  arrange(desc(p)) %>%
  mutate(
    Term_key = paste(Term, GO.ID, sep = "___"),
    Term_key = factor(Term_key, levels = unique(Term_key))
  )

term_to_cat <- setNames(as.character(tbl$Category), as.character(tbl$Term_key))
axis_cols <- unname(category_cols[term_to_cat[levels(tbl$Term_key)]])

## 4. Plot
limits <- c(-10, -2)
breaks_p <- 10^seq(min(limits), max(limits), by = 2)
p_go_sc_pos <- ggplot(
  tbl,
  aes(
    x = Term_key,
    y = p,
    size = Significant,
    colour = Category
  )
) +
  geom_point() +
  coord_flip() +
  labs(
    y = expression(italic(p)~" value"),
    x = NULL,
    title = NULL,
    colour = "Category"
  ) +
  theme_minimal(base_size = 12) +
  scale_x_discrete(
    labels = function(x) str_wrap(sub("___GO:.*$", "", x), width = 52)
  ) +
  scale_y_log10(
    breaks = breaks_p,
    labels = function(x) parse(text = paste0("10^", log10(x))),
    limits = c(1E-10, 0.02)
  ) +
  scale_colour_manual(
    values = category_cols,
    breaks = c(cat1, cat3, cat2, cat4),
    labels = c("Nutrient metabolism", "Hibernation", "Reproduction", "Other")
  ) +
  scale_size_continuous(name = "Nr. of features", range = c(1.2, 5.2)) +
  theme(
    axis.text.y = element_text(hjust = 1, colour = axis_cols, size = 7.2),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    size = guide_legend(order = 1),
    colour = guide_legend(order = 2)
  )
print(p_go_sc_pos)
write_tsv(
  tbl %>% mutate(Term_key = as.character(Term_key)),
  file.path(output_dir, "GO_dotplot_PLS_SC_positive_plot_data.tsv")
)
graph2png(
  x = p_go_sc_pos,
  file = file.path(figure_dir, "Fig1B_source.png"),
  width = 8, height = 11
)
graph2pdf(
  x = p_go_sc_pos,
  file = file.path(figure_dir, "Fig1B_source.pdf"),
  width = 8, height = 11
)
graph2ppt(
  x = p_go_sc_pos,
  file = file.path(figure_dir, "Fig1B_source.pptx"),
  width = 8, height = 11
)

# Main-text Figure 1: editable vector PLS and GO panels on a 4:3 slide.
stage_legend_data <- tibble(
  Species = factor(
    "Vespula vulgaris",
    levels = c("Polistes dominula", "Vespula vulgaris")
  ),
  Stage = factor(c("P", "L5", "L4", "L3", "L2", "L1"), levels = stages),
  Axis1_SC = -41.5,
  Axis2_Stage = seq(8, -32, length.out = 6),
  label = c("P", "L5", "L4", "L3", "L2", "L1")
)
sc_legend_data <- tibble(
  Species = factor(
    "Vespula vulgaris",
    levels = c("Polistes dominula", "Vespula vulgaris")
  ),
  SC2 = factor(
    c("Early / Worker", "Late / Queen"),
    levels = c("Early / Worker", "Late / Queen")
  ),
  Axis1_SC = -41.5,
  Axis2_Stage = c(-55, -67),
  label = c("Early / Worker", "Late / Queen")
)
direction_legend_data <- tibble(
  Species = factor(
    c("Polistes dominula", "Vespula vulgaris"),
    levels = c("Polistes dominula", "Vespula vulgaris")
  ),
  x = 16, xend = 32,
  y = -80,
  label_x = 24,
  label_y = -70,
  label = c("LATE SEASON", "QUEEN-BIASED")
)

p_pls_figure <- p_all +
  geom_text(
    data = tibble(
      Species = factor(
        "Vespula vulgaris",
        levels = c("Polistes dominula", "Vespula vulgaris")
      ),
      Axis1_SC = -43.5, Axis2_Stage = 20, label = "Stage"
    ),
    aes(Axis1_SC, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.1
  ) +
  geom_point(
    data = stage_legend_data,
    aes(Axis1_SC, Axis2_Stage, colour = Stage),
    inherit.aes = FALSE, shape = 15, size = 3.0
  ) +
  geom_text(
    data = stage_legend_data,
    aes(Axis1_SC + 4, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 2.7
  ) +
  geom_text(
    data = tibble(
      Species = factor(
        "Vespula vulgaris",
        levels = c("Polistes dominula", "Vespula vulgaris")
      ),
      Axis1_SC = -43.5, Axis2_Stage = -44, label = "Season / Caste"
    ),
    aes(Axis1_SC, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.1
  ) +
  geom_point(
    data = sc_legend_data,
    aes(Axis1_SC, Axis2_Stage, shape = SC2),
    inherit.aes = FALSE, colour = "black", size = 3.0
  ) +
  geom_text(
    data = sc_legend_data,
    aes(Axis1_SC + 4, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 2.7
  ) +
  geom_segment(
    data = direction_legend_data,
    aes(x = x, y = y, xend = xend, yend = y),
    inherit.aes = FALSE,
    colour = "#7F7F7F",
    linewidth = 0.80,
    lineend = "round",
    arrow = grid::arrow(
      angle = 19, length = grid::unit(0.12, "inches"), type = "closed"
    )
  ) +
  geom_text(
    data = direction_legend_data,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    colour = "#7F7F7F",
    fontface = "bold",
    size = 3.0
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.margin = margin(4, 4, 4, 4)
  )

graph2png(
  x = p_pls_figure,
  file = file.path(figure_dir, "Fig1A_source.png"),
  width = 4.5, height = 7.5
)
graph2pdf(
  x = p_pls_figure,
  file = file.path(figure_dir, "Fig1A_source.pdf"),
  width = 4.5, height = 7.5
)
graph2ppt(
  x = p_pls_figure,
  file = file.path(figure_dir, "Fig1A_source.pptx"),
  width = 4.5, height = 7.5
)

figure1_file <- file.path(figure_dir, "Fig1.pptx")
figure1_pptx <- officer::read_pptx()
figure1_pptx <- officer::add_slide(
  figure1_pptx, layout = "Blank", master = "Office Theme"
)
figure1_pptx <- officer::ph_with(
  figure1_pptx,
  value = rvg::dml(ggobj = p_pls_figure),
  location = officer::ph_location(left = 0.20, top = 0.18, width = 4.15, height = 7.12)
)
figure1_pptx <- officer::ph_with(
  figure1_pptx,
  value = rvg::dml(ggobj = p_go_sc_pos),
  location = officer::ph_location(left = 4.30, top = 0.18, width = 5.55, height = 7.12)
)
figure1_pptx <- officer::ph_with(
  figure1_pptx,
  value = officer::fpar(
    officer::ftext(
      "A",
      officer::fp_text(font.family = "Arial", font.size = 24, bold = TRUE)
    )
  ),
  location = officer::ph_location(left = 0.22, top = 0.08, width = 0.55, height = 0.48)
)
figure1_pptx <- officer::ph_with(
  figure1_pptx,
  value = officer::fpar(
    officer::ftext(
      "B",
      officer::fp_text(font.family = "Arial", font.size = 24, bold = TRUE)
    )
  ),
  location = officer::ph_location(left = 4.32, top = 0.08, width = 0.55, height = 0.48)
)
print(figure1_pptx, target = figure1_file)



### 3.0.2 ENRICHMENT TEST FOR ORTHOGROUPS MOST PREDICTIVE OF EARLY/WORKER SEASON/CASTE ####
pls_SC_neg_go <- .run_pls_go(sel_sc_neg_topK)
pls_SC_neg_go$BP %>% dplyr::filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5) %>% as.data.frame()
# -> 8 sign enriched terms with p<0.02, Significant>=3 & FoldEnrichment>1.5, incl.
# GO.ID                                                                     Term Annotated Significant Expected
# 1 GO:0007188 adenylate cyclase-modulating G protein-coupled receptor signaling        28           7     1.67
# 3 GO:0045187 regulation of circadian sleep/wake cycle                                 23           6     1.37


### 3.0.3 ENRICHMENT TEST FOR ORTHOGROUPS MOST PREDICTIVE OF OLDER STAGES ####
pls_Stage_pos_go <- .run_pls_go(sel_st_pos_topK)
pls_Stage_pos_go$BP %>% dplyr::filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5) %>% as.data.frame()
# -> sign enriched terms related to various developmental processes

### 3.0.4 ENRICHMENT TEST FOR ORTHOGROUPS MOST PREDICTIVE OF YOUNGER STAGES ####
pls_Stage_neg_go <- .run_pls_go(sel_st_neg_topK)
pls_Stage_neg_go$BP %>% dplyr::filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5) %>% as.data.frame()
# -> sign enriched terms related to protein transport, rRNA processing etc

### 3.0.5 ENRICHMENT TEST FOR ORTHOGROUPS MOST PREDICTIVE OF Vespula vulgaris ####
pls_Species_pos_go <- .run_pls_go(sel_sp_pos_topK)
pls_Species_pos_go$BP %>% dplyr::filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5) %>% as.data.frame()
# -> sign enriched for very few processes

### 3.0.6 ENRICHMENT TEST FOR ORTHOGROUPS MOST PREDICTIVE OF Polistes dominula ####
pls_Species_neg_go <- .run_pls_go(sel_sp_neg_topK)
pls_Species_neg_go$BP %>% dplyr::filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5) %>% as.data.frame()
# -> 21 sign enriched terms, e.g. 
# GO.ID                                        Term Annotated Significant Expected
# 8  GO:0009411                              response to UV        49           9     3.29
# 9  GO:0009593              detection of chemical stimulus        24           6     1.61
# 13 GO:0010212              response to ionizing radiation        57           7     3.83
# 17 GO:0007608                 sensory perception of smell        38           7     2.55


flat_all <- bind_rows(
  flatten_go_results(pls_SC_pos_go, "SC_pos"),
  flatten_go_results(pls_SC_neg_go, "SC_neg"),
  flatten_go_results(pls_Stage_pos_go, "Stage_pos"),
  flatten_go_results(pls_Stage_neg_go, "Stage_neg"),
  flatten_go_results(pls_Species_pos_go, "Species_pos"),
  flatten_go_results(pls_Species_neg_go, "Species_neg")
)


# Export complete, unabridged topGO results and contributing N13 HOGs.
pls_go_lists <- list(
  SC_pos = pls_SC_pos_go, SC_neg = pls_SC_neg_go,
  Stage_pos = pls_Stage_pos_go, Stage_neg = pls_Stage_neg_go,
  Species_pos = pls_Species_pos_go, Species_neg = pls_Species_neg_go
)
go_term_results_all <- purrr::imap_dfr(pls_go_lists, function(res, label) {
  purrr::imap_dfr(res, function(tbl_go, ontology_name) {
    if (is.null(tbl_go) || !nrow(tbl_go)) return(tibble())
    tbl_go %>%
      mutate(label = label, ontology = ontology_name) %>%
      dplyr::select(label, ontology, everything(), -any_of("Members"))
  })
}) %>%
  left_join(go_lookup, by = "GO.ID") %>%
  mutate(
    Term = coalesce(term_name, Term),
    annotation_set = "final_direct_conserved_plus_Dmel_experimental",
    foreground_top_k = K,
    topGO_algorithm = "weight01",
    topGO_statistic = "fisher",
    nodeSize = 20L,
    multiple_testing_correction = "none"
  ) %>%
  dplyr::select(
    annotation_set, foreground_top_k, topGO_algorithm,
    topGO_statistic, nodeSize, multiple_testing_correction,
    everything(), -term_name
  )

go_term_results_enriched <- go_term_results_all %>%
  filter(p <= 0.02, Significant >= 3, FoldEnrichment >= 1.5)

pls_foregrounds <- list(
  SC_pos = sel_sc_pos_topK,
  SC_neg = sel_sc_neg_topK,
  Stage_pos = sel_st_pos_topK,
  Stage_neg = sel_st_neg_topK,
  Species_pos = sel_sp_pos_topK,
  Species_neg = sel_sp_neg_topK
) %>%
  purrr::imap_dfr(~ tibble(label = .y, HOG = .x)) %>%
  mutate(foreground_top_k = K, annotation_set = "final_direct_conserved_plus_Dmel_experimental")

write_tsv(go_term_results_all, file.path(output_dir, "topGO_PLS_all_terms.tsv"))
write_tsv(
  go_term_results_enriched,
  file.path(output_dir, "topGO_PLS_enriched_terms.tsv")
)
write_tsv(flat_all, file.path(output_dir, "topGO_PLS_contributing_N13_HOGs.tsv"))
write_tsv(pls_foregrounds, file.path(output_dir, "topGO_PLS_foreground_N13_HOGs.tsv"))
saveRDS(pls_go_lists, file.path(output_dir, "topGO_PLS_complete_results.rds"))

# One best row per orthogroup across all tests
flat_best <- flat_all %>%
  filter(is.finite(p)) %>%                                  # drop NA/Inf p's
  arrange(orthogroup, p, desc(FoldEnrichment), desc(Significant)) %>% 
  distinct(orthogroup, .keep_all = TRUE) %>%                # keep first per orthogroup
  ungroup()
length(unique(flat_best$orthogroup)) # 1677 out of 1800 annotated
sum(sel_sc_pos %in% unique(flat_best$orthogroup)) # 103 out of 129 in sel_sc_pos annotated
# Best per axis sign (e.g., SC_pos vs SC_neg)
flat_best_by_label <- flat_all %>%
  filter(is.finite(p)) %>%
  arrange(label, orthogroup, p, desc(FoldEnrichment), desc(Significant)) %>%
  distinct(label, orthogroup, .keep_all = TRUE) %>%
  ungroup()
length(unique(flat_best_by_label$orthogroup)) # 1677
# Best per ontology within each orthogroup
flat_best_by_ont <- flat_all %>%
  filter(is.finite(p)) %>%
  arrange(orthogroup, ontology, p, desc(FoldEnrichment), desc(Significant)) %>%
  distinct(orthogroup, ontology, .keep_all = TRUE) %>%
  ungroup()
length(unique(flat_best_by_ont$orthogroup)) # 1677


## Keep only contributing HOG rows for terms that pass the display criteria.
flat_all_filt <- flat_all %>%
  group_by(label, ontology, GO.ID) %>%
  mutate(keep_flag = any(!is.na(p) & p <= 0.02 &
                           !is.na(Significant) & Significant >= 3 &
                           !is.na(FoldEnrichment) & FoldEnrichment >= 1.5)) %>%
  ungroup() %>%
  filter(keep_flag) %>%
  dplyr::select(-keep_flag)


# 4. CROSS-STAGE DIRECTIONAL CONCORDANCE AND HETEROCHRONY ####

source(file.path(
  dir_base, "scripts", "modules", "B_cross_stage_concordance.R"
), local = TRUE)
source(file.path(
  dir_base, "scripts", "modules", "B_IIS_Wnt_mTOR_FoxO_panel.R"
), local = TRUE)

# 5. HEATMAPS OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED N13 HOGS IN BOTH SPECIES ####

library(dplyr)
library(tidyr)
library(stringr)
library(gtools)
library(pheatmap)

## -------------------- settings
top <- 300 # (all in sel_sc_pos_topK)
lfc_choice <- "shrunk"    # "raw" or "shrunk"
stages_use <- stages      # c("L1","L2","L3","L4","L5","P")
padj_cutoff <- 0.01

# NEW: Vv effect-size threshold used jointly with padj
# >0 enforces positive Vv DE, <0 enforces negative Vv DE
log2FC_cutoff <- 1      # adjust as desired

# Choose which PLS set to plot
sel_vec <- sel_sc_pos_topK
# sel_vec <- sel_sc_pos
# sel_vec <- sel_sc_neg
# sel_vec <- sel_sc_neg_topK

# Visual
cap <- 3
show_gaps <- TRUE         # now uses FAKE gaps (spacer rows)
clip_atzero <- TRUE      # TRUE to show only sign-consistent values

## -------------------- 1) select top PLS orthogroups in desired order
sel_ids <- sel_vec[1:min(top, length(sel_vec))]

wide_sel <- tibble(orthogroup = sel_ids) %>%
  left_join(ortholog_de_wide, by = "orthogroup") %>%
  mutate(orthogroup = factor(orthogroup, levels = sel_ids)) %>%
  arrange(orthogroup)

if (nrow(wide_sel) == 0) stop("No matching orthogroups found in ortholog_de_wide.")

## -------------------- 2) enforce stronger Vv filter:
## keep only genes with >=1 Vv stage where padj < padj_cutoff
## AND log2FC exceeds log2FC_cutoff in the same stage

vv_padj_cols <- sprintf("padj_Vv_%s_QvsW", stages_use)
vv_lfc_cols  <- sprintf("log2FC%s_Vv_%s_QvsW", lfc_choice, stages_use)

vv_padj_cols <- intersect(vv_padj_cols, names(wide_sel))
vv_lfc_cols  <- intersect(vv_lfc_cols,  names(wide_sel))

if (length(vv_padj_cols) != length(stages_use) ||
    length(vv_lfc_cols)  != length(stages_use)) {
  stop("Not all expected Vv padj/LFC columns found. Check column naming.")
}

P <- as.matrix(wide_sel %>% dplyr::select(all_of(vv_padj_cols)))
L <- as.matrix(wide_sel %>% dplyr::select(all_of(vv_lfc_cols)))

# stage-wise logical filter
sig_ok <- !is.na(P) & (P < padj_cutoff)

lfc_ok <- if (log2FC_cutoff >= 0) {
  !is.na(L) & (L > log2FC_cutoff)
} else {
  !is.na(L) & (L < log2FC_cutoff)
}

good_stage_mat <- sig_ok & lfc_ok

wide_sel <- wide_sel %>%
  mutate(anyDE_Vv_thresh = rowSums(good_stage_mat) > 0) %>%
  filter(anyDE_Vv_thresh)

if (nrow(wide_sel) == 0) {
  stop("After Vv padj+log2FC threshold filtering, no rows remain.")
}

## -------------------- 3) extract Pd and Vv LFC matrices
pd_cols <- sprintf("log2FC%s_Pd_%s_LvsE", lfc_choice, stages_use)
vv_cols <- sprintf("log2FC%s_Vv_%s_QvsW", lfc_choice, stages_use)

pd_cols <- intersect(pd_cols, names(wide_sel))
vv_cols <- intersect(vv_cols, names(wide_sel))

if (length(pd_cols) != length(stages_use) ||
    length(vv_cols) != length(stages_use)) {
  stop("Not all expected Pd/Vv LFC columns found. Check column naming.")
}

mat_pd <- wide_sel %>% dplyr::select(all_of(pd_cols)) %>% as.matrix()
mat_vv <- wide_sel %>% dplyr::select(all_of(vv_cols)) %>% as.matrix()

colnames(mat_pd) <- stages_use
colnames(mat_vv) <- stages_use

## -------------------- 4) build Pd/Vv star matrices based on padj
pd_padj_cols <- sprintf("padj_Pd_%s_LvsE", stages_use)
pd_padj_cols <- intersect(pd_padj_cols, names(wide_sel))

if (length(pd_padj_cols) != length(stages_use)) {
  stop("Not all expected Pd padj columns found. Check column naming.")
}

mat_pd_padj <- wide_sel %>% dplyr::select(all_of(pd_padj_cols)) %>% as.matrix()
mat_vv_padj <- wide_sel %>% dplyr::select(all_of(vv_padj_cols)) %>% as.matrix()

colnames(mat_pd_padj) <- stages_use
colnames(mat_vv_padj) <- stages_use

# Use ordinary ASCII significance symbols for compatibility with the editable
# PowerPoint exporter. Their vertical positions are corrected in the grid grob
# below because periods and asterisks sit at different heights on a text line.
heatmap_significance_symbol <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  if (p < 0.1) return(".")
  ""
}

pd_stars <- apply(mat_pd_padj, c(1,2), heatmap_significance_symbol)
vv_stars <- apply(mat_vv_padj, c(1,2), heatmap_significance_symbol)

## -------------------- 5) row labels from description_pd (fallbacks)
base_labels <- if ("description_pd" %in% names(wide_sel)) {
  as.character(wide_sel$description_pd)
} else {
  rep(NA_character_, nrow(wide_sel))
}

if ("gene_pd" %in% names(wide_sel)) {
  idx <- is.na(base_labels) | !nzchar(base_labels)
  base_labels[idx] <- wide_sel$gene_pd[idx]
}
idx <- is.na(base_labels) | !nzchar(base_labels)
base_labels[idx] <- as.character(wide_sel$orthogroup[idx])

base_labels <- make.unique(base_labels)

## -------------------- 6) pre-sort genes:
## by stage where mean(Pd,Vv) is maximal
mean_mat <- sapply(seq_along(stages_use), function(j) {
  rowMeans(cbind(mat_pd[, j], mat_vv[, j]), na.rm = TRUE)
})
colnames(mean_mat) <- stages_use

stage_max_idx <- apply(mean_mat, 1, which.max)
stage_max_val <- apply(mean_mat, 1, max, na.rm = TRUE)

ord <- order(stage_max_idx, -stage_max_val)

mat_pd <- mat_pd[ord, , drop = FALSE]
mat_vv <- mat_vv[ord, , drop = FALSE]
pd_stars <- pd_stars[ord, , drop = FALSE]
vv_stars <- vv_stars[ord, , drop = FALSE]
base_labels <- base_labels[ord]

n_genes <- length(base_labels)

## -------------------- 7) build interleaved matrix with optional FAKE gaps
rows_list <- list()
stars_list <- list()
labels_row <- character()
species_vec <- character()
rn_internal <- character()

for (i in seq_len(n_genes)) {
  # Pd row
  rows_list[[length(rows_list) + 1]] <- mat_pd[i, ]
  stars_list[[length(stars_list) + 1]] <- pd_stars[i, ]
  labels_row <- c(labels_row, base_labels[i])
  species_vec <- c(species_vec, "Pd")
  rn_internal <- c(rn_internal, paste0(base_labels[i], "__Pd"))
  
  # Vv row
  rows_list[[length(rows_list) + 1]] <- mat_vv[i, ]
  stars_list[[length(stars_list) + 1]] <- vv_stars[i, ]
  labels_row <- c(labels_row, "")
  species_vec <- c(species_vec, "Vv")
  rn_internal <- c(rn_internal, paste0(base_labels[i], "__Vv"))
  
  # Fake gap row after each pair except last
  if (show_gaps && i < n_genes) {
    rows_list[[length(rows_list) + 1]] <- rep(NA, length(stages_use))
    stars_list[[length(stars_list) + 1]] <- rep("", length(stages_use))
    labels_row <- c(labels_row, "")
    species_vec <- c(species_vec, " ")
    rn_internal <- c(rn_internal, paste0("gap_", i))
  }
}

labels_row <- gsub(", mitochondrial-like", "", labels_row)
labels_row <- gsub(", mitochondrial", "", labels_row)
mat_plot <- do.call(rbind, rows_list)
mat_stars_plot <- do.call(rbind, stars_list)

colnames(mat_plot) <- stages_use
colnames(mat_stars_plot) <- stages_use
rownames(mat_plot) <- rn_internal
rownames(mat_stars_plot) <- rn_internal

## -------------------- 8) row annotation (species; blank for gap rows)
row_anno <- data.frame(species = factor(species_vec, levels = c("Pd","Vv","")))
rownames(row_anno) <- rn_internal

anno_cols <- list(
  species = c(Pd = "grey90", Vv = "grey70", ` ` = "white")
)

## -------------------- 9) optional clipping at zero + capping
mat_disp <- mat_plot

if (clip_atzero) {
  if (log2FC_cutoff >= 0) {
    mat_disp <- pmax(mat_disp, 0)
    mat_disp <- pmin(mat_disp, cap)
  } else {
    mat_disp <- pmin(mat_disp, 0)
    mat_disp <- pmax(mat_disp, -cap)
  }
} else {
  mat_disp <- pmax(pmin(mat_disp, cap), -cap)
}

## -------------------- 10) choose color scale
if (clip_atzero) {
  if (log2FC_cutoff >= 0) {
    cols <- colorRampPalette(c("yellow", "red"))(200)
    breaks <- seq(0, cap, length.out = 201)
  } else {
    cols <- colorRampPalette(c("blue", "yellow"))(200)
    breaks <- seq(-cap, 0, length.out = 201)
  }
} else {
  cols <- colorRampPalette(c("blue", "yellow", "red"))(200)
  breaks <- seq(-cap, cap, length.out = 201)
}

## -------------------- 11) plot (no clustering), with stars
p_heatmap <- pheatmap(
  mat_disp,
  color = cols,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_row = labels_row,
  labels_col = stages_use,
  # annotation_row = row_anno,
  annotation_colors = anno_cols,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 11,
  angle_col = 0,
  display_numbers = mat_stars_plot,
  number_color = "black",
  na_col = "grey95",
  fontsize_number = 7 # ,
  # main = paste0(
  #   "Temporal DE dynamics (top=", min(top, length(sel_ids)), ") | ",
  #   "Vv filter: padj<", padj_cutoff, " & ",
  #   if (log2FC_cutoff >= 0) paste0("log2FC>", log2FC_cutoff) else paste0("log2FC<", log2FC_cutoff),
  #   if (clip_atzero) " | clipped at 0" else ""
  # )
)

# PowerPoint and PDF exporters calculate different text-box metrics for ".",
# "*", "**", and "***". Replace these text glyphs with vector geometry so all
# significance marks have exactly the same vertical centre at every zoom level.
draw_pheatmap_significance_symbols <- function(pheatmap_object,
                                                star_radius_pt = 0.65,
                                                star_spacing_pt = 1.9,
                                                dot_size_pt = 0.6,
                                                line_width = 0.45) {
  matrix_index <- which(pheatmap_object$gtable$layout$name == "matrix")
  if (length(matrix_index) != 1L) {
    stop("Could not identify the pheatmap matrix grob.")
  }

  matrix_grob <- pheatmap_object$gtable$grobs[[matrix_index]]
  text_index <- grep("text", names(matrix_grob$children))
  if (length(text_index) != 1L) {
    stop("Could not identify the pheatmap significance-label grob.")
  }

  number_grob <- matrix_grob$children[[text_index]]
  labels <- as.character(number_grob$label)
  x <- as.numeric(number_grob$x)
  y <- as.numeric(number_grob$y)

  dot_index <- which(labels == ".")
  dot_grob <- if (length(dot_index)) {
    grid::pointsGrob(
      x = grid::unit(x[dot_index], "npc"),
      y = grid::unit(y[dot_index], "npc"),
      pch = 16,
      size = grid::unit(dot_size_pt, "pt"),
      gp = grid::gpar(col = "black", fill = "black")
    )
  } else {
    grid::nullGrob()
  }

  star_index <- which(grepl("^\\*+$", labels))
  if (length(star_index)) {
    star_x_npc <- numeric()
    star_x_pt <- numeric()
    star_y_npc <- numeric()

    for (i in star_index) {
      n_stars <- nchar(labels[i])
      horizontal_offsets <- (
        seq_len(n_stars) - (n_stars + 1) / 2
      ) * star_spacing_pt
      star_x_npc <- c(star_x_npc, rep(x[i], n_stars))
      star_x_pt <- c(star_x_pt, horizontal_offsets)
      star_y_npc <- c(star_y_npc, rep(y[i], n_stars))
    }

    # Each asterisk is three centred line segments at 90, 30 and -30 degrees.
    angles <- rep(c(pi / 2, pi / 6, -pi / 6), times = length(star_x_npc))
    centre_x_npc <- rep(star_x_npc, each = 3L)
    centre_x_pt <- rep(star_x_pt, each = 3L)
    centre_y_npc <- rep(star_y_npc, each = 3L)
    dx <- star_radius_pt * cos(angles)
    dy <- star_radius_pt * sin(angles)

    star_grob <- grid::segmentsGrob(
      x0 = grid::unit(centre_x_npc, "npc") +
        grid::unit(centre_x_pt - dx, "pt"),
      y0 = grid::unit(centre_y_npc, "npc") - grid::unit(dy, "pt"),
      x1 = grid::unit(centre_x_npc, "npc") +
        grid::unit(centre_x_pt + dx, "pt"),
      y1 = grid::unit(centre_y_npc, "npc") + grid::unit(dy, "pt"),
      gp = grid::gpar(col = "black", lwd = line_width, lineend = "round")
    )
  } else {
    star_grob <- grid::nullGrob()
  }

  matrix_grob$children[[text_index]] <- grid::grobTree(dot_grob, star_grob)
  pheatmap_object$gtable$grobs[[matrix_index]] <- matrix_grob
  pheatmap_object
}

p_heatmap <- draw_pheatmap_significance_symbols(p_heatmap)

graph2png(fun = function() grid::grid.draw(p_heatmap$gtable), file=file.path(supplemental_figure_dir, "FigS3.png"), width=8, height=6)
graph2pdf(fun = function() grid::grid.draw(p_heatmap$gtable), file=file.path(supplemental_figure_dir, "FigS3.pdf"), width=8, height=6)
graph2ppt(fun = function() grid::grid.draw(p_heatmap$gtable), file=file.path(supplemental_figure_dir, "FigS3.pptx"), width=8, height=6)



# 6. MANUSCRIPT SUPPLEMENTARY TABLES ####

# Compact Word-compatible summaries underlying Figs. 1 and 2. Complete
# gene/HOG-level outputs remain available in the analysis output directory and
# are intentionally not reproduced as very large tables in the supplement.
supplementary_xlsx <- file.path(
  supplemental_table_dir, "TablesS4-S9.xlsx"
)

summarise_de_tests <- function(dat, species, unit, contrast,
                               effect_col, padj_col = "padj_global") {
  dat %>%
    mutate(stage = factor(stage, levels = stages)) %>%
    group_by(stage) %>%
    summarise(
      tested = sum(is.finite(pvalue)),
      significant = sum(.data[[padj_col]] < alpha_de, na.rm = TRUE),
      higher_in_late_or_queen = sum(
        .data[[padj_col]] < alpha_de & .data[[effect_col]] > 0,
        na.rm = TRUE
      ),
      higher_in_early_or_worker = sum(
        .data[[padj_col]] < alpha_de & .data[[effect_col]] < 0,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Species = species,
      Analysis_unit = unit,
      Contrast = contrast,
      Stage = as.character(stage),
      .before = 1
    ) %>%
    dplyr::select(
      Species, Analysis_unit, Stage, Contrast, tested, significant,
      higher_in_late_or_queen, higher_in_early_or_worker
    )
}

supp_table_s4 <- bind_rows(
  summarise_de_tests(
    gene_de_pd, "Polistes dominula", "Gene", "Late versus early",
    "log2FC_shrunk"
  ),
  summarise_de_tests(
    gene_de_vv, "Vespula vulgaris", "Gene", "Queen versus worker",
    "log2FC_shrunk"
  ),
  summarise_de_tests(
    as_tibble(hog_obj_pd$results), "Polistes dominula", "N13 HOG",
    "Late versus early", "log2FC_ash"
  ),
  summarise_de_tests(
    as_tibble(hog_obj_vv$results), "Vespula vulgaris", "N13 HOG",
    "Queen versus worker", "log2FC_ash"
  )
) %>%
  arrange(factor(Species, c("Polistes dominula", "Vespula vulgaris")),
          factor(Analysis_unit, c("Gene", "N13 HOG")),
          factor(Stage, stages))

axis_targets <- c("Season/caste", "Developmental stage", "Species")
axis_map <- unname(map[c("SC", "Stage", "Species")])
axis_variance <- as.numeric(plsfit$prop_expl_var$X[axis_map]) * 100
supp_table_s5a <- tibble(
  Axis = paste("Axis", 1:3),
  Biological_contrast = axis_targets,
  PLS_component = paste("Component", axis_map),
  X_variance_explained_percent = round(axis_variance, 2),
  HOGs_in_final_model = ncol(Z_sub)
)

supp_table_s5 <- as_tibble(as.data.frame(sc_contr_oos_result)) %>%
  transmute(
    Species = recode(as.character(Species),
                     Pd = "Polistes dominula", Vv = "Vespula vulgaris"),
    Stage = as.character(Stage),
    Contrast = as.character(contrast),
    Estimate = signif(estimate, 4),
    SE = signif(SE, 4),
    df = round(df, 1),
    t_ratio = signif(t.ratio, 4),
    One_sided_P = signif(p.value, 3)
  ) %>%
  arrange(factor(Species, c("Polistes dominula", "Vespula vulgaris")),
          factor(Stage, stages))

category_labels <- c(
  "Differential feeding / nutrient metabolism" = "Nutrient metabolism",
  "Hibernation / diapause / long-term survival" = "Hibernation",
  "Reproduction / caste-related development" = "Reproduction",
  "Other" = "Other"
)

# Build a focal-species consensus name for each N13 HOG. Symbols shared by
# Polistes and Vespula are preferred; informative focal symbols and then focal
# descriptions provide deterministic fallbacks when annotations differ.
split_focal_annotation <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(str_trim(x))) {
    return(character())
  }
  values <- str_trim(unlist(str_split(x, fixed("|"))))
  unique(values[nzchar(values)])
}

informative_focal_symbols <- function(x) {
  values <- split_focal_annotation(x)
  values[!str_detect(
    values,
    regex("^(LOC[0-9]+|ENS[A-Z0-9._-]+|CG[0-9]+|NA|N/A)$", ignore_case = TRUE)
  )]
}

informative_focal_descriptions <- function(x) {
  values <- split_focal_annotation(x) %>%
    str_remove(regex("^(LOW QUALITY PROTEIN:|PREDICTED:)\\s*", ignore_case = TRUE))
  values[!str_detect(
    values,
    regex("^(uncharacterized|hypothetical protein|predicted protein|NA|N/A)(\\s|$)",
          ignore_case = TRUE)
  )]
}

normalise_focal_name <- function(x) {
  x %>%
    str_remove(regex("-like$", ignore_case = TRUE)) %>%
    str_squish() %>%
    str_to_lower()
}

consensus_focal_hog_name <- function(gene_pd, symbol_pd, description_pd,
                                     gene_vv, symbol_vv, description_vv) {
  pd_symbols <- informative_focal_symbols(symbol_pd)
  vv_symbols <- informative_focal_symbols(symbol_vv)
  if (length(pd_symbols) > 0L && length(vv_symbols) > 0L) {
    shared <- intersect(
      normalise_focal_name(pd_symbols), normalise_focal_name(vv_symbols)
    )
    if (length(shared) > 0L) {
      return(pd_symbols[match(shared[[1]], normalise_focal_name(pd_symbols))])
    }
  }

  symbols <- unique(c(pd_symbols, vv_symbols))
  if (length(symbols) > 0L) {
    return(paste(head(symbols, 3L), collapse = " / "))
  }

  pd_descriptions <- informative_focal_descriptions(description_pd)
  vv_descriptions <- informative_focal_descriptions(description_vv)
  if (length(pd_descriptions) > 0L && length(vv_descriptions) > 0L) {
    shared <- intersect(
      normalise_focal_name(pd_descriptions),
      normalise_focal_name(vv_descriptions)
    )
    if (length(shared) > 0L) {
      matched <- pd_descriptions[
        match(shared[[1]], normalise_focal_name(pd_descriptions))
      ]
      return(str_remove(matched, regex("-like$", ignore_case = TRUE)))
    }
  }

  descriptions <- unique(c(pd_descriptions, vv_descriptions))
  if (length(descriptions) > 0L) {
    return(paste(head(descriptions, 2L), collapse = " / "))
  }

  focal_gene_ids <- unique(c(
    str_trim(unlist(str_split(coalesce(gene_pd, ""), "[,|]"))),
    str_trim(unlist(str_split(coalesce(gene_vv, ""), "[,|]")))
  ))
  focal_gene_ids <- focal_gene_ids[nzchar(focal_gene_ids)]
  if (length(focal_gene_ids) > 0L) {
    return(paste0(
      "uncharacterized (", paste(head(focal_gene_ids, 2L), collapse = " / "), ")"
    ))
  }
  "uncharacterized HOG"
}

n13_consensus_names <- orth_annotation %>%
  transmute(
    HOG,
    consensus_HOG_name = pmap_chr(
      list(gene_pd, symbol_pd, description_pd,
           gene_vv, symbol_vv, description_vv),
      consensus_focal_hog_name
    ) %>%
      str_replace_all(fixed(";"), ",") %>%
      str_squish()
  ) %>%
  distinct(HOG, .keep_all = TRUE)

consensus_name_lookup <- setNames(
  n13_consensus_names$consensus_HOG_name,
  n13_consensus_names$HOG
)

collapse_contributing_hog_names <- function(hog_string) {
  if (length(hog_string) == 0L || is.na(hog_string) || !nzchar(str_trim(hog_string))) {
    return(NA_character_)
  }
  hogs <- str_trim(unlist(str_split(hog_string, fixed(";"))))
  hogs <- hogs[nzchar(hogs)]
  consensus_names <- unname(consensus_name_lookup[hogs])
  consensus_names[is.na(consensus_names) | !nzchar(consensus_names)] <-
    "uncharacterized HOG"
  paste(consensus_names, collapse = "; ")
}

supp_table_s6 <- tbl %>%
  left_join(
    go_term_results_all %>%
      dplyr::select(GO.ID, Annotated) %>%
      distinct(),
    by = "GO.ID"
  ) %>%
  transmute(
    Category = unname(category_labels[as.character(Category)]),
    GO_ID = GO.ID,
    GO_term = as.character(Term),
    Annotated_HOGs = Annotated,
    Foreground_HOGs = Significant,
    Expected_HOGs = round(Expected, 2),
    Fold_enrichment = signif(FoldEnrichment, 3),
    weight01_P = signif(p, 3),
    `Contributing HOGs` = MembersStr,
    `Contributing consensus HOG names` = map_chr(
      MembersStr, collapse_contributing_hog_names
    )
  ) %>%
  arrange(factor(Category, c("Nutrient metabolism", "Hibernation",
                             "Reproduction", "Other")), weight01_P)

supp_table_s7 <- pairwise_concordance_results %>%
  transmute(
    Polistes_stage = as.character(pd_stage),
    Vespula_stage = as.character(vv_stage),
    Jointly_DE_HOGs = jointly_DE_n,
    Concordant_HOGs = concordant_DE_n,
    Concordant_fraction = signif(concordance_fraction, 3),
    Expected_fraction_from_margins =
      signif(expected_concordance_fraction_from_margins, 3),
    Concordance_log2_odds_ratio = signif(concordance_log2_odds_ratio_HA, 3),
    One_sided_Fisher_P = signif(fisher_one_sided_p_greater, 3),
    FDR_P_36_tests = signif(fisher_p_FDR_36, 3)
  ) %>%
  arrange(factor(Vespula_stage, stages), factor(Polistes_stage, stages))

supp_table_selected_hog_effects <- heatmap_data %>%
  mutate(
    column = as.character(column),
    effect = sprintf("%.3f (%s)", log2FC, if_else(de_significance == "", "ns", de_significance))
  ) %>%
  dplyr::select(HOG, column, effect) %>%
  pivot_wider(names_from = column, values_from = effect)

supp_table_selected_hogs <- selected_hogs %>%
  left_join(orth_annotation, by = "HOG") %>%
  left_join(n13_consensus_names, by = "HOG") %>%
  left_join(supp_table_selected_hog_effects, by = "HOG") %>%
  mutate(
    Selection_basis = case_when(
      !is.na(enriched_GO_evidence) & !is.na(BH_significant_KEGG_evidence) ~
        "Qualifying GO and BH-significant KEGG enrichment",
      !is.na(enriched_GO_evidence) ~ "Qualifying GO enrichment",
      !is.na(BH_significant_KEGG_evidence) ~
        "BH-significant KEGG enrichment",
      TRUE ~ "Connected component in the IIS/Wnt/mTOR/FoxO network"
    )
  ) %>%
  transmute(
    HOG,
    Displayed_gene = plot_label,
    Functional_group = functional_theme,
    Selected_comparison = selected_from,
    `Polistes L1 late vs early` = `Polistes L1\nlate vs early`,
    `Polistes L4 late vs early` = `Polistes L4\nlate vs early`,
    `Polistes L5 late vs early` = `Polistes L5\nlate vs early`,
    `Vespula L2 queen vs worker` = `Vespula L2\nqueen vs worker`,
    `Consensus HOG name` = consensus_HOG_name,
    `Polistes gene IDs` = gene_pd,
    `Polistes symbols` = symbol_pd,
    `Polistes descriptions` = description_pd,
    `Vespula gene IDs` = gene_vv,
    `Vespula symbols` = symbol_vv,
    `Vespula descriptions` = description_vv,
    Selection_basis,
    `Supporting enriched GO terms` = enriched_GO_evidence,
    `Minimum topGO weight01 P` = signif(minimum_topGO_weight01_p, 3),
    `Supporting BH-significant KEGG pathways` = BH_significant_KEGG_evidence,
    `Minimum KEGG BH P` = signif(minimum_KEGG_BH_p, 3),
    `Wang et al. 2021 pathway component` = Wang_2021_component,
    `Maximum concordance z` = signif(max_concordance_z, 3)
  ) %>%
  arrange(
    factor(Functional_group, levels = theme_patterns$functional_theme),
    Displayed_gene
  )

concordance_comparison_labels <- c(
  L4_to_L2 = "Polistes L4 to Vespula L2",
  P_to_P = "Pupae to pupae"
)
concordance_direction_labels <- c(
  concordant_up = "Concordant up",
  concordant_down = "Concordant down"
)
supp_table_concordant_go <- topgo_results_all %>%
  filter(
    comparison %in% names(concordance_comparison_labels),
    set %in% names(concordance_direction_labels),
    weight01_p <= 0.02,
    foreground_n >= 3,
    fold_enrichment >= 1.5
  ) %>%
  transmute(
    Comparison = unname(concordance_comparison_labels[comparison]),
    Direction = unname(concordance_direction_labels[set]),
    `GO ID` = GO,
    `GO term` = term_name,
    `Annotated HOGs` = background_n,
    `Foreground HOGs` = foreground_n,
    `Expected HOGs` = signif(expected_n, 3),
    `Fold enrichment` = signif(fold_enrichment, 3),
    `weight01 P` = signif(weight01_p, 3),
    `Contributing HOGs` = foreground_HOGs,
    `Contributing consensus HOG names` = map_chr(
      foreground_HOGs, collapse_contributing_hog_names
    )
  ) %>%
  arrange(
    factor(Comparison, unname(concordance_comparison_labels)),
    factor(Direction, unname(concordance_direction_labels)),
    `weight01 P`, `GO ID`
  )

# Supplementary numbering follows the analytical sequence: enrichment results
# precede the smaller mechanistic set selected from those results for Fig. 2B.
supp_table_s8 <- supp_table_concordant_go
supp_table_s9 <- supp_table_selected_hogs

supp_table_index <- tibble(
  Table = c("S4", "S5", "S6", "S7", "S8", "S9"),
  Worksheet = c(
    "S4_DE_summary", "S5_PLS_LOOCV",
    "S6_Fig1_GO", "S7_Fig2A_concordance", "S8_concordant_GO",
    "S9_Fig2B_HOGs"
  ),
  Description = c(
    "Differential-expression counts by species, analysis unit and stage",
    "Planned one-sided season/caste contrasts from fully nested leave-one-sample-out PLS axis-1 scores",
    "GO terms displayed in Fig. 1B",
    "Pairwise directional-concordance tests underlying Fig. 2A",
    "GO enrichment of concordantly regulated HOGs for Polistes L4 to Vespula L2 and pupae to pupae",
    "Functionally supported concordant HOGs displayed in Fig. 2B"
  )
)

supplementary_tables <- list(
  Table_index = supp_table_index,
  S4_DE_summary = supp_table_s4,
  S5_PLS_LOOCV = supp_table_s5,
  S6_Fig1_GO = supp_table_s6,
  S7_Fig2A_concordance = supp_table_s7,
  S8_concordant_GO = supp_table_s8,
  S9_Fig2B_HOGs = supp_table_s9
)

supplementary_table_source_dir <- file.path(output_dir, "workbook_sources")
dir.create(supplementary_table_source_dir, recursive = TRUE, showWarnings = FALSE)
purrr::iwalk(
  supplementary_tables,
  ~ readr::write_tsv(.x, file.path(supplementary_table_source_dir, paste0(.y, ".tsv")))
)

# Manuscript-facing filenames follow the final supplementary numbering. The
# more descriptive source filenames above are retained for the XLSX builder.
numbered_supplementary_tables <- list(
  TableS4 = supp_table_s4,
  TableS5 = supp_table_s5,
  TableS6 = supp_table_s6,
  TableS7 = supp_table_s7,
  TableS8 = supp_table_s8,
  TableS9 = supp_table_s9
)
purrr::iwalk(
  numbered_supplementary_tables,
  ~ readr::write_tsv(.x, file.path(supplemental_table_dir, paste0(.y, ".tsv")))
)
readr::write_tsv(supp_table_s5a, file.path(output_dir, "pls_axis_summary.tsv"))
legacy_s5_files <- c(
  file.path(supplemental_table_dir, "TableS5A.tsv"),
  file.path(supplemental_table_dir, "TableS5B.tsv"),
  file.path(supplementary_table_source_dir, "S5A_PLS_axes.tsv"),
  file.path(supplementary_table_source_dir, "S5B_PLS_LOOCV.tsv")
)
invisible(file.remove(legacy_s5_files[file.exists(legacy_s5_files)]))

legacy_s8_s9_files <- c(
  file.path(supplementary_table_source_dir, "S8_Fig2B_HOGs.tsv"),
  file.path(supplementary_table_source_dir, "S8_Fig2B_ridge.tsv"),
  file.path(supplementary_table_source_dir, "S9_concordant_GO.tsv")
)
invisible(file.remove(legacy_s8_s9_files[file.exists(legacy_s8_s9_files)]))

find_python <- function() {
  bundled_python <- file.path(
    Sys.getenv("USERPROFILE"), ".cache", "codex-runtimes",
    "codex-primary-runtime", "dependencies", "python", "python.exe"
  )
  candidates <- c(
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    bundled_python,
    Sys.which("python"),
    Sys.which("python3")
  )
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  usable <- vapply(
    candidates,
    function(candidate) {
      identical(
        suppressWarnings(system2(
          candidate,
          args = c("-c", shQuote("import openpyxl")),
          stdout = FALSE, stderr = FALSE
        )),
        0L
      )
    },
    logical(1)
  )
  candidates <- candidates[usable]
  if (!length(candidates)) {
    stop("Python with openpyxl is required to export the supplementary XLSX workbook")
  }
  candidates[[1]]
}

supplementary_workbook_builder <- file.path(
  dir_base, "scripts", "utilities", "build_supplementary_transcriptomic_workbook.py"
)
if (!file.exists(supplementary_workbook_builder)) {
  stop("Missing workbook builder: ", supplementary_workbook_builder)
}
supplementary_workbook_status <- system2(
  find_python(),
  args = c(
    shQuote(supplementary_workbook_builder),
    shQuote(supplementary_table_source_dir),
    shQuote(supplementary_xlsx)
  )
)
if (!identical(supplementary_workbook_status, 0L)) {
  stop("Supplementary workbook export failed with status ", supplementary_workbook_status)
}

# openxlsx can leave stale drawing relationships in workbooks that contain no
# drawings under some Windows installations. Rebuild these detailed workbooks
# with openpyxl so Excel and non-Microsoft readers open them without repair.
workbook_rebuilder <- file.path(
  dir_base, "scripts", "utilities", "rebuild_openxlsx_workbook.py"
)
if (!file.exists(workbook_rebuilder)) {
  stop("Missing workbook rebuilder: ", workbook_rebuilder)
}
detailed_workbooks <- file.path(
  cross_stage_output_dir,
  c(
    "concordant_DE_HOGs_annotated.xlsx",
    "cross_stage_concordance_results.xlsx",
    "IIS_Wnt_mTOR_FoxO_pathway_HOGs.xlsx"
  )
)
for (workbook_path in detailed_workbooks[file.exists(detailed_workbooks)]) {
  rebuild_status <- system2(
    find_python(),
    args = c(
      shQuote(workbook_rebuilder),
      shQuote(workbook_path),
      shQuote(workbook_path)
    )
  )
  if (!identical(rebuild_status, 0L)) {
    stop("Detailed workbook rebuild failed for ", workbook_path,
         " with status ", rebuild_status)
  }
}




# 7. SAVE FINAL ANALYSIS OBJECTS ####
repository_relative_path <- function(path) {
  path_abs <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_abs <- normalizePath(dir_base, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root_abs, "/")
  if (startsWith(tolower(path_abs), tolower(prefix))) {
    return(substring(path_abs, nchar(prefix) + 1L))
  }
  basename(path_abs)
}

analysis_metadata_final <- tibble(
  key = c(
    "script", "analysis_date", "gene_model", "HOG_model",
    "orthology_unit", "multiple_testing", "DEU_included",
    "PLS_features", "PLS_cross_validation", "PLS_score_contrasts",
    "PLS_X_variance_by_component",
    "GO_annotation", "GO_foreground", "GO_test", "GO_multiple_testing",
    "cross_stage_tests", "cross_stage_multiple_testing",
    "N13_orthology_composition",
    "full_gene_DE_results", "full_N13_HOG_DE_results",
    "Figure1", "Figure2", "Supplementary_tables"
  ),
  value = c(
    "B_differential expression analysis.R", as.character(Sys.Date()),
    "glmmTMB nbinom2 with gene-by-sample length/library offset and colonyNested random intercept",
    "same model fitted directly to tximport-summarised N13 HOG counts; ASH shrinkage by species and stage",
    "OrthoFinder hierarchical orthogroups at node N13",
    "global BH within species for differential-expression tests",
    "FALSE",
    as.character(length(sel)),
    paste(
      "fully nested leave-one-sample-out; expression filtering, TMM reference",
      "selection, HOG-wise centring/scaling, VIP/Kneedle selection, PLS fitting,",
      "component assignment and orientation estimated from training samples only"
    ),
    paste(
      "planned right-tailed late/queen minus early/worker contrasts within each",
      "stage and species; no multiplicity adjustment"
    ),
    paste(signif(plsfit$prop_expl_var$X, 6), collapse = ";"),
    "direct EXCON/Galaxy EggNOG and InterProScan conserved across wasps plus experimental FlyBase GO transferred through N13 Drosophila orthologues",
    "top 300 N13 HOGs in each oriented PLS loading direction",
    "topGO weight01 Fisher; nodeSize 20",
    "none; nominal topGO p values exported",
    paste0(
      "pairwise one-sided Fisher tests of directional concordance; all-stage ",
      "nonnegative cumulative-logit and inverse-variance-weighted ridge models; ",
      bootstrap_B, " pairs-bootstrap replicates with plus-one correction"
    ),
    "Benjamini-Hochberg FDR correction across all 36 stage pairs or coefficients per analysis",
    repository_relative_path(n13_orthology_summary_file),
    repository_relative_path(full_gene_de_file),
    repository_relative_path(full_hog_de_file),
    repository_relative_path(figure1_file),
    repository_relative_path(figure2_file),
    repository_relative_path(supplementary_xlsx)
  )
)
write_tsv(analysis_metadata_final, file.path(output_dir, "analysis_metadata.tsv"))
write_tsv(
  pairwise_concordance_results %>%
    mutate(cell_label = str_replace_all(cell_label, "\\n", " ")),
  file.path(output_dir, "pairwise_directional_concordance_results.tsv")
)
write_tsv(
  nonnegative_clm_results %>%
    mutate(panel_label = str_replace_all(panel_label, "\\n", " ")),
  file.path(output_dir, "all_stage_nonnegative_CLM_results.tsv")
)
write_tsv(
  all_stage_ridge_results,
  file.path(output_dir, "all_stage_nonnegative_ridge_results.tsv")
)
key_object_names <- intersect(
  c(
    "gene_de_pd", "gene_de_vv", "hog_de_shared",
    "n13_composition", "n13_orthology_composition_summary",
    "ortholog_de_wide", "de_gene_categories", "de_stacked_counts",
    "de_stacked_counts_plot",
    "scores_export", "centroid_arrows", "orthologs_pls", "pls_loocv_fold_diagnostics", "pls_go_lists",
    "go_term_results_all", "go_term_results_enriched", "flat_all",
    "pls_foregrounds", "tbl", "p_go_sc_pos", "p_pls_figure",
    "pairwise_concordance_results", "p_pairwise_concordance",
    "topgo_results_all", "kegg_results_all",
    "nonnegative_clm_results", "p_nonnegative_clm",
    "all_stage_ridge_results", "p_all_stage_ridge",
    "selected_hogs", "heatmap_data", "figure2_ab_plot", "figure_s4_plot",
    "supp_table_s4", "supp_table_s5a", "supp_table_s5",
    "n13_consensus_names",
    "supp_table_s6", "supp_table_s7", "supp_table_s8", "supp_table_s9",
    "wide_sel", "mat_pd", "mat_vv"
  ),
  ls(envir = .GlobalEnv)
)
key_results <- mget(key_object_names, envir = .GlobalEnv)
saveRDS(
  key_results,
  file.path(output_dir, "B_differential_expression_key_results.rds"),
  compress = "xz"
)

validation_checks <- c(
  fully_nested_PLS_LOOCV =
    nrow(pls_loocv_fold_diagnostics) == nrow(meta_all) &&
    all(pls_loocv_fold_diagnostics$candidate_HOGs > 0L) &&
    all(pls_loocv_fold_diagnostics$selected_HOGs > 0L),
  final_GO_input = identical(
    basename(go_annot_orthologs), "N13_HOG_GO_final_long.tsv.gz"
  ),
  top_300_foreground = identical(K, 300),
  full_GO_names = !any(str_detect(go_orthologs$term_name, fixed("..."))),
  no_duplicate_HOG_GO_pairs = !anyDuplicated(
    go_orthologs[c("HOG", "go_terms")]
  ),
  reproduction_terms_plotted = all(c(
    "female gonad development",
    "morphogenesis of follicular epithelium"
  ) %in% tbl$Term),
  complete_topGO_export = file.exists(file.path(
    output_dir, "topGO_PLS_all_terms.tsv"
  )),
  figure1_exported = file.exists(figure1_file),
  figure1_panel_labels = all(c("A", "B") %in% officer::pptx_summary(
    officer::read_pptx(figure1_file)
  )$text),
  figure1_category_order = identical(
    levels(tbl$Category), c(cat1, cat3, cat2, cat4)
  ),
  figure2_exported = file.exists(figure2_file),
  figure2_panel_labels = all(c("A", "B") %in% officer::pptx_summary(
    officer::read_pptx(figure2_file)
  )$text),
  finite_nonzero_bootstrap_p = all(c(
    nonnegative_clm_results$bootstrap_p_one_sided,
    all_stage_ridge_results$bootstrap_p_one_sided
  ) > 0 & c(
    nonnegative_clm_results$bootstrap_p_one_sided,
    all_stage_ridge_results$bootstrap_p_one_sided
  ) <= 1),
  global_36_test_FDR = all(c(
    pairwise_concordance_results$fisher_p_FDR_36,
    nonnegative_clm_results$bootstrap_p_FDR_36,
    all_stage_ridge_results$bootstrap_p_FDR_36
  ) >= 0 & c(
    pairwise_concordance_results$fisher_p_FDR_36,
    nonnegative_clm_results$bootstrap_p_FDR_36,
    all_stage_ridge_results$bootstrap_p_FDR_36
  ) <= 1),
  n13_orthology_summary_exported = file.exists(n13_orthology_summary_file),
  n13_shared_relationship_classes_complete = sum(
    n13_orthology_composition_summary$count[
      n13_orthology_composition_summary$summary_group ==
        "Shared-HOG composition"
    ]
  ) == n13_shared_hogs,
  full_DE_exports_exist = all(file.exists(c(
    full_gene_de_file, full_hog_de_file
  ))),
  full_DE_exports_under_100_MB = all(file.info(c(
    full_gene_de_file, full_hog_de_file
  ))$size < 100 * 1024^2),
  full_DE_global_FDR_flags_consistent = all(c(
    identical(
      full_gene_de_results$global_FDR_lt_0.05,
      is.finite(full_gene_de_results$global_FDR) &
        full_gene_de_results$global_FDR < alpha_de
    ),
    identical(
      full_hog_de_results$global_FDR_lt_0.05,
      is.finite(full_hog_de_results$global_FDR) &
        full_hog_de_results$global_FDR < alpha_de
    )
  )),
  supplementary_workbook_exported = file.exists(supplementary_xlsx),
  supplementary_GO_names_complete = !any(c(
    str_detect(supp_table_s6$GO_term, fixed("...")),
    str_detect(supp_table_s8$`GO term`, fixed("..."))
  )),
  supplementary_contributing_HOG_names_aligned = all(c(
    map2_lgl(
      supp_table_s6$`Contributing HOGs`,
      supp_table_s6$`Contributing consensus HOG names`,
      ~ length(str_split_1(.x, fixed(";"))) ==
        length(str_split_1(.y, fixed(";")))
    ),
    map2_lgl(
      supp_table_s8$`Contributing HOGs`,
      supp_table_s8$`Contributing consensus HOG names`,
      ~ length(str_split_1(.x, fixed(";"))) ==
        length(str_split_1(.y, fixed(";")))
    )
  )),
  supplementary_contributing_HOG_names_complete = all(c(
    !is.na(supp_table_s6$`Contributing consensus HOG names`) &
      nzchar(supp_table_s6$`Contributing consensus HOG names`),
    !is.na(supp_table_s8$`Contributing consensus HOG names`) &
      nzchar(supp_table_s8$`Contributing consensus HOG names`)
  )),
  supplementary_selected_HOG_provenance_complete =
    !anyDuplicated(supp_table_s9$HOG) &&
    setequal(supp_table_s9$HOG, selected_hogs$HOG) &&
    all(!is.na(supp_table_s9$`Consensus HOG name`)) &&
    all(nzchar(supp_table_s9$Selection_basis)),
  supplementary_tables_nonempty = all(vapply(
    supplementary_tables[-1], nrow, integer(1)
  ) > 0L)
)
if (!all(validation_checks)) {
  stop(
    "Final analysis validation failed: ",
    paste(names(validation_checks)[!validation_checks], collapse = ", ")
  )
}
writeLines(
  c(
    "VALIDATION_OK",
    paste(names(validation_checks), validation_checks, sep = "\t")
  ),
  file.path(output_dir, "VALIDATION.txt")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("B_differential expression analysis completed. Outputs: ", output_dir)

