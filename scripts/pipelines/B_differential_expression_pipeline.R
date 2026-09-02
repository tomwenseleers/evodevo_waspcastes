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
    panel.border = element_rect(fill = NA, linewidth = 0.5),
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
  
  if (vip_preselect) { sel  <- vip_select(Xtr, Ytr, min_keep = min_keep)$sel_union } else { 
    sel <- colnames(Xtr) }
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
    Axis3_Species = f3 * as.numeric(SCOt[[ map["Species"]]])
  )
}

## Leave-one-out driver producing out-of-sample axis scores
loo_pls_emm <- function(X, meta, min_keep = 100, vip_preselect = TRUE) {
  stopifnot(nrow(X) == nrow(meta))
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
  for (i in seq_len(n)) {
    te <- i; tr <- setdiff(seq_len(n), te)
    pr <- fit_project_once(Xtr = X[tr,,drop=FALSE], 
                           Xte = X[te,,drop=FALSE],
                           meta_tr = meta[tr,,drop=FALSE], 
                           min_keep = min_keep, 
                           vip_preselect = vip_preselect)
    out[[i]] <- tibble(
      sample         = meta$sample[te],
      Axis1_SC       = pr$Axis1_SC,
      Axis2_Stage    = pr$Axis2_Stage,
      Axis3_Species  = pr$Axis3_Species
    )
  }
  dplyr::bind_rows(out) %>%
    dplyr::left_join(meta[, c("sample","Stage","Species","SC2")], by = "sample")
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
p_all <- ggplot(scores_all2,
                aes(Axis1_SC, Axis2_Stage, color = Stage, shape = SC2)) +
  facet_wrap(~ Species, ncol = 1) +
  geom_polygon(data = hull_all2,
               aes(group = interaction(Species, Stage, SC2),
                   fill = Stage, color = Stage, alpha = SC2),
               linewidth = 0.25, inherit.aes = TRUE) +
  geom_point(size = 2) +
  scale_shape_manual("Season / Caste", values = c(17, 16)) +
  scale_alpha_manual("Season / Caste", values = c(0.22, 0.60)) +
  scale_color_manual(values = rev(brewer.pal(6, "Spectral")),
                     guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = rev(brewer.pal(6, "Spectral")),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Axis 1: Season / Caste") +
  ylab("Axis 2: Stage") +
  # labs(title = "Combined PLS — Pd + Vv") +
  theme_few() +
  theme(panel.border = element_rect(fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "italic", size = 10),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=9))
print(p_all)
graph2png(x = ggplot2::last_plot(), file =file.path(figure_dir, "Fig1A_source.png"), width=4.5, height=5.5) 
graph2pdf(x = ggplot2::last_plot(), file =file.path(figure_dir, "Fig1A_source.pdf"), width=4.5, height=5.5)
graph2ppt(x = ggplot2::last_plot(), file =file.path(figure_dir, "Fig1A_source.pptx"), width=4.5, height=5.5)

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

# calculate leave-one-out sample scores
oos <- loo_pls_emm(X = Z_sub, meta = meta_all, min_keep = 100, 
                   vip_preselect = vip_preselect)

# calculate significance of emmeans contrasts for season/caste, stage & species based on leave-one-out sample scores
# TODO change lm_robust to MASS::lmrob & use marginaleffects package for significance
# (not supported by emmeans) to be consistent with robust models later on?

afex::set_sum_contrasts()
library(estimatr)
fit_axis1_oos <- estimatr::lm_robust(Axis1_SC ~ (SC2 + Stage + Species)^2, data = oos) 
sc_contr_oos  <- emmeans::contrast(emmeans::emmeans(fit_axis1_oos, ~ SC2 | Stage*Species),
                                   method = "revpairwise")

fit_axis2_oos <- estimatr::lm_robust(Axis2_Stage ~ (SC2 + Stage + Species)^2, data = oos)
st_contr_oos  <- emmeans::contrast(emmeans::emmeans(fit_axis2_oos, ~ Stage | Species),
                                   method = "consec")

fit_axis3_oos <- estimatr::lm_robust(Axis3_Species ~ (SC2 + Stage + Species)^2, data = oos)
sp_contr_oos  <- emmeans::contrast(emmeans::emmeans(fit_axis3_oos, ~ Species | Stage),
                                   method = "revpairwise")

# Season/Caste contrasts per stage & species
write.csv(data.frame(sc_contr_oos),  file = file.path(output_dir, "pls_season caste emmeans contrasts LOOCV scores.csv"), row.names = FALSE)
sc_contr_oos
# Stage = L1, Species = Pd:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    15.16 3.68 124   4.125  0.0001
# 
# Stage = L2, Species = Pd:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    18.61 2.71 124   6.857  <.0001
# 
# Stage = L3, Species = Pd:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    14.40 3.06 124   4.713  <.0001
# 
# Stage = L4, Species = Pd:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    16.63 2.30 124   7.237  <.0001
# 
# Stage = L5, Species = Pd:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    12.18 2.62 124   4.649  <.0001
# 
# Stage = P, Species = Pd:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    15.74 3.41 124   4.617  <.0001
# 
# Stage = L1, Species = Vv:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)     7.31 3.11 124   2.352  0.0203
# 
# Stage = L2, Species = Vv:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)    10.76 2.23 124   4.818  <.0001
# 
# Stage = L3, Species = Vv:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)     6.55 2.60 124   2.522  0.0129
# 
# Stage = L4, Species = Vv:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)     8.77 2.41 124   3.644  0.0004
# 
# Stage = L5, Species = Vv:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)     4.32 2.18 124   1.984  0.0495
# 
# Stage = P, Species = Vv:
#   contrast      estimate   SE  df t.ratio p.value
# (L/Q) - (E/W)     7.89 3.06 124   2.581  0.0110

# -> significant shift in expected direction across all stages for both species

# stage contrasts per species
write.csv(data.frame(st_contr_oos),  file = file.path(output_dir, "pls_stage emmeans contrasts LOOCV scores.csv"), row.names = FALSE)
st_contr_oos
# Species = Pd:
#   contrast estimate   SE  df t.ratio p.value
# L2 - L1      8.00 9.52 124   0.840  0.8833
# L3 - L2     13.10 6.51 124   2.013  0.1855
# L4 - L3      7.38 4.05 124   1.822  0.2692
# L5 - L4     23.14 4.17 124   5.554  <.0001
# P - L5      11.34 4.05 124   2.798  0.0268
# 
# Species = Vv:
#   contrast estimate   SE  df t.ratio p.value
# L2 - L1      8.56 4.75 124   1.802  0.2745
# L3 - L2     24.53 6.51 124   3.767  0.0011
# L4 - L3     21.15 5.44 124   3.890  0.0006
# L5 - L4     20.49 3.79 124   5.410  <.0001
# P - L5       4.37 3.22 124   1.359  0.5438
# 
# Results are averaged over the levels of: SC2 
# P value adjustment: mvt method for 5 tests 

# -> biggest shift from L4 to L5 and L5 to P in Polistes dominula (late in development),
# but in L2 to L3, L3 to L4 and L4 to L5 in Vespula vulgaris (shifted earlier in development)

# species contrasts per stage
write.csv(data.frame(sp_contr_oos),  file = file.path(output_dir, "pls_species emmeans contrasts LOOCV scores.csv"), row.names = FALSE)
sp_contr_oos
# Stage = L1:
#   contrast estimate    SE  df t.ratio p.value
# Vv - Pd      84.8 2.250 124  37.744  <.0001
# 
# Stage = L2:
#   contrast estimate    SE  df t.ratio p.value
# Vv - Pd      81.3 1.850 124  44.001  <.0001
# 
# Stage = L3:
#   contrast estimate    SE  df t.ratio p.value
# Vv - Pd      76.4 2.260 124  33.748  <.0001
# 
# Stage = L4:
#   contrast estimate    SE  df t.ratio p.value
# Vv - Pd      71.9 0.930 124  77.230  <.0001
# 
# Stage = L5:
#   contrast estimate    SE  df t.ratio p.value
# Vv - Pd      67.8 1.270 124  53.316  <.0001
# 
# Stage = P:
#   contrast estimate    SE  df t.ratio p.value
# Vv - Pd      61.0 0.747 124  81.643  <.0001
# 
# Results are averaged over the levels of: SC2

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
    "Polistes dominula",
    levels = c("Polistes dominula", "Vespula vulgaris")
  ),
  Stage = factor(c("P", "L5", "L4", "L3", "L2", "L1"), levels = stages),
  Axis1_SC = -39,
  Axis2_Stage = seq(-7, -57, length.out = 6),
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
  Axis1_SC = -39,
  Axis2_Stage = c(-46, -58),
  label = c("Early / Worker", "Late / Queen")
)

p_pls_figure <- p_all +
  geom_text(
    data = tibble(
      Species = factor(
        "Polistes dominula",
        levels = c("Polistes dominula", "Vespula vulgaris")
      ),
      Axis1_SC = -41, Axis2_Stage = 4, label = "Stage"
    ),
    aes(Axis1_SC, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.7
  ) +
  geom_point(
    data = stage_legend_data,
    aes(Axis1_SC, Axis2_Stage, colour = Stage),
    inherit.aes = FALSE, shape = 15, size = 3.3
  ) +
  geom_text(
    data = stage_legend_data,
    aes(Axis1_SC + 4, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.2
  ) +
  geom_text(
    data = tibble(
      Species = factor(
        "Vespula vulgaris",
        levels = c("Polistes dominula", "Vespula vulgaris")
      ),
      Axis1_SC = -41, Axis2_Stage = -34, label = "Season / Caste"
    ),
    aes(Axis1_SC, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.7
  ) +
  geom_point(
    data = sc_legend_data,
    aes(Axis1_SC, Axis2_Stage, shape = SC2),
    inherit.aes = FALSE, colour = "black", size = 3.3
  ) +
  geom_text(
    data = sc_legend_data,
    aes(Axis1_SC + 4, Axis2_Stage, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3.2
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.margin = margin(4, 4, 4, 4)
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


# 4. ARE PAPER WASP GENE EXPRESSION PROFILES PREDICTIVE OF COMMON WASP CASTE DIFFERENCES? ####

### 4.1 SCATTERPLOT OF UNIVARIATE ROBUST REGRESSION MODELS OF DE IN Vv vs THAT IN SAME STAGE IN Pd ####

# robust regression models of Log2FC in Q vs W in Vespula vulgaris
# in different larval stages in function of Log2FC in L vs E in Polistes dominula
# with genes in both matched by orthogroup
# to test if gene regulatory modules overexpressed in late-season female larvae
# in primitively eusocial without morphologically distinct castes were co-opted 
# in the evolution of morphologically specialised castes

## ---- packages

library(dplyr)
library(rlang)
library(marginaleffects)
library(robustbase)
library(forcats)
# assumes `stages`, `ortholog_de_wide`, `sel_sc_pos_topK`, `sel_sc_neg_topK` exist

## ---- helpers

# choose rows (orthogroups) for the requested "universe"
pick_universe <- function(wide, universe = c("all_orthologs","de_either","de_both", "vv_de_any","pls_top"),
                          alpha = 0.05, pls_ids = NULL) {
  universe <- match.arg(universe)
  out <- wide
  
  # fast "any DE" flags from padj columns already inside `wide`
  vv_p <- as.matrix(dplyr::select(out, dplyr::starts_with("padj_Vv_")))
  pd_p <- as.matrix(dplyr::select(out, dplyr::starts_with("padj_Pd_")))
  anyDE_Vv <- rowSums((vv_p < alpha) & !is.na(vv_p), na.rm = TRUE) > 0
  anyDE_Pd <- rowSums((pd_p < alpha) & !is.na(pd_p), na.rm = TRUE) > 0
  
  keep <- switch(universe,
                 all_orthologs = rep(TRUE, nrow(out)),
                 de_either     = anyDE_Vv | anyDE_Pd,
                 de_both       = anyDE_Vv & anyDE_Pd,
                 vv_de_any     = anyDE_Vv,
                 pls_top       = out$orthogroup %in% pls_ids
  )
  out[keep, , drop = FALSE]
}


# build per-stage data with toggle for raw vs shrunk LFCs/SEs
make_stage_df <- function(stage, wide, lfc = c("raw","shrunk"), se = c("raw","auto","shrunk")) {
  lfc <- match.arg(lfc)
  se  <- match.arg(se)
  
  # dynamic names for LFCs
  xcol  <- paste0("log2FC", lfc, "_Pd_", stage, "_LvsE")
  ycol  <- paste0("log2FC", lfc, "_Vv_", stage, "_QvsW")
  
  # pick which SE to use
  se_eff <- if (se == "auto") lfc else se
  
  # expected SE column name
  wycol_try <- paste0("lfcSE", se_eff, "_Vv_", stage, "_QvsW")
  
  # graceful fallback if that column doesn't exist
  wycol <- if (wycol_try %in% names(wide)) {
    wycol_try
  } else {
    paste0("lfcSEraw_Vv_", stage, "_QvsW")
  }
  
  # other dynamic names
  padjP <- paste0("padj_Pd_", stage, "_LvsE")
  padjV <- paste0("padj_Vv_", stage, "_QvsW")
  bmP   <- paste0("baseMean_Pd_", stage, "_LvsE")
  bmV   <- paste0("baseMean_Vv_", stage, "_QvsW")
  
  out <- wide %>%
    transmute(
      orthogroup,
      stage = factor(stage, levels = stages),
      x  = !!sym(xcol),
      y  = !!sym(ycol),
      wy = !!sym(wycol),
      padj_pd = !!sym(padjP),
      padj_vv = !!sym(padjV),
      expr_gm_stage_best = sqrt(pmax(!!sym(bmP), 1e-8) * pmax(!!sym(bmV), 1e-8))
    ) %>%
    mutate(
      sig_both = (padj_pd < 0.05) & (padj_vv < 0.05),
      sig_lab  = ifelse(sig_both, "FDR < 0.05 in both", "NS in one/both"),
      w        = 1 / (wy^2 + 1e-8)
    ) %>%
    filter(is.finite(x), is.finite(y), is.finite(w))
  
  out
}

# Fit robust simple regressions y ~ x per stage 
fit_one_stage <- function(dat) {
  f <- try(lmrob(y ~ x, data = dat, weights = w,
                 init = "S",
                 control = lmrob.control(max.it = 5000, refine.tol = 1e-5, setting = "KS2014")),
           silent = TRUE)
  if (inherits(f, "try-error")) lm(y ~ x, data = dat, weights = w) else f
}

## assemble modeling data
# choose gene universe & options here:
# "all_orthologs" (primary), 
# "de_either" (DE in at least 1 stage in either species),
# "de_both" (DE in at least 1 stage in both species),
# "vv_de_any" (DE in at least 1 stage in Vv - to test predictive power for caste biased genes in Vv), 
#  or "pls_top" (top PLS selected features along season/caste axis)
universe_choice <- "vv_de_any" 
alpha_universe  <- 0.05
lfc_choice      <- "shrunk"
se_choice       <- "raw"        # or "auto"/"shrunk"
padj_method     <- "bonferroni" # or "BH" for FDR

wide_use <- pick_universe(
  ortholog_de_wide,
  universe = universe_choice,
  alpha    = alpha_universe,
  pls_ids  = c(sel_sc_pos_topK, sel_sc_neg_topK)
)

df <- dplyr::bind_rows(lapply(
  stages,
  make_stage_df,
  wide = wide_use,
  lfc  = lfc_choice,
  se   = se_choice
))

## ---- stage-wise fits & summaries
slopes_by_stage <- df %>%
  group_by(stage) %>%
  group_modify(~{
    fit <- fit_one_stage(.x)
    s   <- slopes(fit, variables = "x", newdata = datagrid(), conf_level = 0.95)
    tibble(
      estimate  = s$estimate[1],
      std.error = s$std.error[1],
      conf.low  = s$conf.low[1],
      conf.high = s$conf.high[1],
      p.value   = s$p.value[1],
      n         = nrow(.x)
    )
  }) %>%
  ungroup() %>%
  mutate(p.value.adj = p.adjust(p.value, padj_method)) %>%
  arrange(stage)

print(slopes_by_stage)
dim(slopes_by_stage) # 6 8 - correct

## ---- predictions for plotting
pred <- df %>%
  group_by(stage) %>%
  do({
    fit <- fit_one_stage(.)
    xr  <- range(.$x, finite = TRUE)
    grid <- tibble(x = seq(xr[1], xr[2], length.out = 200))
    p <- as.data.frame(marginaleffects::predictions(fit, newdata = grid, conf_level = 0.95))
    tibble(stage = unique(.$stage), x = p$x, estimate = p$estimate,
           conf.low = p$conf.low, conf.high = p$conf.high)
  }) %>% ungroup()


# --- Plot
p_stagewise_regression <- ggplot() +
  geom_ribbon(data = pred, aes(x, ymin = conf.low, ymax = conf.high, fill = stage),
              alpha = 0.15, colour = NA) +
  geom_point(data = df, aes(x, y, colour = stage, alpha = sig_lab), shape = 16) +
  scale_alpha_manual(values = c("FDR < 0.05 in both" = 1, "NS in one/both" = 0.05),
                     name = "", guide = "none") +
  geom_line(data = pred, aes(x, y = estimate, colour = stage), linewidth = 1) +
  scale_colour_brewer(palette = "Spectral", direction = -1, guide = guide_legend(reverse = TRUE)) +
  scale_fill_brewer(palette = "Spectral", direction = -1, guide = guide_legend(reverse = TRUE)) +
  labs(x = expression(italic("Polistes dominula")~log[2]*"FC (L vs E)"),
       y = expression(italic("Vespula vulgaris")~log[2]*"FC (Q vs W)")#,
       # subtitle = paste("Universe:", universe_choice, "| LFC:", lfc_choice),
       #title = "Stage-wise robust regressions (Vv ~ Pd)"
       ) +
  theme_few(base_size = 12) +
  geom_hline(yintercept=0, colour=alpha("black", 0.1)) +
  geom_vline(xintercept=0, colour=alpha("black", 0.1)) +
  theme(axis.line = element_line(colour = NA),
        axis.ticks = element_line(linewidth = 0.3, colour="black"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-7, 7))
graph2png(x = p_stagewise_regression, file =file.path(figure_dir, "Fig2A_source.png"), width=6, height=5)
graph2pdf(x = p_stagewise_regression, file =file.path(figure_dir, "Fig2A_source.pdf"), width=6, height=5)
graph2ppt(x = p_stagewise_regression, file =file.path(figure_dir, "Fig2A_source.pptx"), width=6, height=5)



### 4.2 TEST FOR HETEROCHRONY USING MULTIVARIATE DISTRIBUTED LAG ROBUST OR NONNEGATIVE RIDGE REGRESSIONS OF DE IN EACH STAGE OF Vv ~ SAME OR OTHER STAGES OF Pd ####

library(glmnet)
library(purrr)
library(tidyr)

# Build per-stage design matrix with flexible universe/LFC/SE/predictors
# predictors = "all"           -> use all Pd stages as X
#            = "older_or_equal"-> use Pd stages whose index >= y_stage index
#            = "" -> 
# LFC choice: lfc = "raw" | "shrunk"
# SE choice : se  = "auto" (match LFC if present, else fall back to raw/shrunk)
# scaling   : scale_X = TRUE to standardize columns of X

build_mult_data2 <- function(wide, stages, y_stage,
                             lfc = c("raw","shrunk"),
                             se  = c("auto","raw","shrunk"),
                             predictors = c("all","older_or_equal"),
                             scale_X = TRUE) {
  lfc <- match.arg(lfc)
  se  <- match.arg(se)
  predictors <- match.arg(predictors)
  
  # response and weights (Vv)
  ycol <- sprintf("log2FC%s_Vv_%s_QvsW", lfc, y_stage)
  
  if (se == "auto") {
    candidates <- c(sprintf("lfcSE%s_Vv_%s_QvsW", lfc, y_stage),
                    sprintf("lfcSEraw_Vv_%s_QvsW",    y_stage),
                    sprintf("lfcSEshrunk_Vv_%s_QvsW", y_stage))
    wcol <- candidates[candidates %in% names(wide)][1]
  } else {
    wcol <- sprintf("lfcSE%s_Vv_%s_QvsW", se, y_stage)
  }
  if (!wcol %in% names(wide)) stop("Weight SE column not found: ", wcol)
  
  # choose predictor stages
  if (predictors == "all") {
    use_stages <- stages
  } else {  # "older_or_equal": index >= y_stage
    idx_y <- match(y_stage, stages)
    if (is.na(idx_y)) stop("y_stage not found in 'stages'.")
    use_stages <- stages[seq.int(idx_y, length(stages))]
  }
  
  src_cols <- sprintf("log2FC%s_Pd_%s_LvsE", lfc, use_stages)
  new_cols <- paste0("x_", use_stages)
  
  dat <- wide %>%
    transmute(
      y = .data[[ycol]],
      w = 1 / ((.data[[wcol]])^2 + 1e-8),
      !!! rlang::set_names(rlang::syms(src_cols), new_cols)
    ) %>%
    dplyr::filter(is.finite(y), is.finite(w)) %>%
    tidyr::drop_na(dplyr::all_of(new_cols))
  
  if (nrow(dat) > 0L && scale_X) {
    dat[new_cols] <- lapply(dat[new_cols], function(z) as.numeric(scale(z)))
  }
  attr(dat, "x_names") <- new_cols
  dat
}

# Weighted least squares with intercept and optional slope constraint (single predictor)
wls_single <- function(dat, x_name, lower = -Inf) {
  x <- dat[[x_name]]
  y <- dat$y
  w <- dat$w
  
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  
  if (length(x) < 3) return(list(intercept = NA_real_, slope = NA_real_))
  
  wsum <- sum(w)
  wx <- sum(w * x) / wsum
  wy <- sum(w * y) / wsum
  
  num <- sum(w * (x - wx) * (y - wy))
  den <- sum(w * (x - wx)^2)
  
  b_uc <- if (den > 0) num / den else NA_real_
  a_uc <- wy - b_uc * wx
  
  # Exact constrained solution for slope >= 0:
  # if unconstrained slope < 0, optimum is at boundary b=0
  if (is.finite(lower) && lower == 0 && is.finite(b_uc) && b_uc < 0) {
    list(intercept = wy, slope = 0)
  } else {
    list(intercept = a_uc, slope = b_uc)
  }
}

bootstrap_wls_single <- function(dat, x_name, B = 200, lower = -Inf, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(dat)
  out <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    fb  <- wls_single(dat[idx, , drop = FALSE], x_name, lower = lower)
    out[b] <- fb$slope
  }
  out
}

# Fit weighted ridge with CV, optional coefficient constraints
fit_ridge_cv <- function(dat, x_names,
                         lower.limits = -Inf, upper.limits = Inf,
                         nfolds = 10, standardize = FALSE, ...) {
  X <- as.matrix(dat[, x_names, drop = FALSE])
  y <- dat$y
  w <- dat$w
  
  # glmnet expects lower/upper limits length = p or 1
  p <- ncol(X)
  if (length(lower.limits) == 1) lower.limits <- rep(lower.limits, p)
  if (length(upper.limits) == 1) upper.limits <- rep(upper.limits, p)
  
  cvfit <- cv.glmnet(
    x = X, y = y,
    alpha = 0,
    weights = w,
    nfolds = nfolds,
    standardize = standardize,
    intercept = TRUE,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    ...
  )
  
  # extract coefficients at lambda.min
  b <- as.matrix(coef(cvfit, s = "lambda.min"))
  # named vector for predictors only
  beta <- setNames(rep(NA_real_, length(x_names)), x_names)
  common <- intersect(x_names, rownames(b))
  beta[common] <- b[common, 1]
  
  list(cvfit = cvfit, beta = beta)
}


# Bootstrap ridge coefficients (pairs bootstrap)
bootstrap_ridge <- function(dat, x_names,
                            B = 200,
                            lower.limits = -Inf, upper.limits = Inf,
                            nfolds = 10, standardize = FALSE,
                            refit_cv = TRUE,
                            seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(dat)
  p <- length(x_names)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = p,
                     dimnames = list(NULL, x_names))
  
  # Option: if you prefer fixed lambda from original fit,
  # set refit_cv = FALSE and pass lambda_fixed externally.
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    d   <- dat[idx, , drop = FALSE]
    
    # If too small or degenerate, skip
    if (nrow(d) < 10) next
    
    fb <- try(
      fit_ridge_cv(
        d, x_names,
        lower.limits = lower.limits,
        upper.limits = upper.limits,
        nfolds = nfolds,
        standardize = standardize,
        ...
      ),
      silent = TRUE
    )
    if (inherits(fb, "try-error")) next
    
    boot_mat[b, ] <- fb$beta
  }
  
  boot_mat
}


# Finite two-sided test based on the centred bootstrap error distribution.
# The plus-one correction prevents impossible P = 0 values with finite B.
centred_bootstrap_p <- function(beta_hat, boot) {
  boot <- boot[is.finite(boot)]
  if (!length(boot) || !is.finite(beta_hat)) return(NA_real_)

  null_distance <- abs(beta_hat)
  bootstrap_error <- abs(boot - beta_hat)
  min(1, (sum(bootstrap_error >= null_distance) + 1) / (length(boot) + 1))
}

# Summarise bootstrap distribution into estimate/SE/CI/p
summarise_boot <- function(beta_hat, boot_mat) {
  # boot_mat: B x p
  se <- apply(boot_mat, 2, sd, na.rm = TRUE)
  ci <- t(apply(boot_mat, 2, function(z) {
    stats::quantile(z, c(0.025, 0.975), na.rm = TRUE)
  }))
  
  # Two-sided centred-bootstrap P value with finite-sample correction.
  pval <- vapply(seq_along(beta_hat), function(j) {
    centred_bootstrap_p(beta_hat[j], boot_mat[, j])
  }, numeric(1))
  
  tibble(
    estimate  = as.numeric(beta_hat),
    std.error = as.numeric(se),
    conf.low  = as.numeric(ci[, 1]),
    conf.high = as.numeric(ci[, 2]),
    p.value   = as.numeric(pval)
  )
}

fit_mult_stage2 <- function(dat, y_stage,
                            model = c("lmrob","ridge"),
                            # ridge options
                            lower.limits = -Inf,
                            upper.limits = Inf,
                            ridge_B = 200,
                            ridge_nfolds = 10,
                            ridge_standardize = FALSE,
                            ridge_refit_cv = TRUE,
                            ridge_seed = NULL) {
  model <- match.arg(model)
  
  x_names <- attr(dat, "x_names")
  
  if (length(x_names) == 0L || nrow(dat) < 10) {
    return(tibble(
      vv_stage = y_stage,
      pd_stage = character(),
      estimate = numeric(),
      std.error = numeric(),
      p.value = numeric(),
      conf.low = numeric(),
      conf.high = numeric(),
      n = nrow(dat)
    ))
  }
  
  if (model == "lmrob") {
    form <- as.formula(paste("y ~", paste(x_names, collapse = " + ")))
    fit <- try(
      robustbase::lmrob(
        form, data = dat, weights = w,
        init = "S",
        control = robustbase::lmrob.control(
          max.it = 5000, refine.tol = 1e-5, setting = "KS2014"
        )
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) fit <- lm(form, data = dat, weights = w)
    
    sm <- summary(fit)$coefficients
    sm <- sm[intersect(rownames(sm), x_names), , drop = FALSE]
    
    # lmrob/lm can supply SE and p directly;
    # CI not previously returned, but we now add it for consistency
    est <- sm[, "Estimate"]
    se  <- sm[, "Std. Error"]
    p   <- sm[, ncol(sm)]
    ci_low  <- est - 1.96 * se
    ci_high <- est + 1.96 * se
    
    tibble(
      vv_stage = y_stage,
      pd_stage = sub("^x_", "", rownames(sm)),
      estimate = as.numeric(est),
      std.error = as.numeric(se),
      conf.low = as.numeric(ci_low),
      conf.high = as.numeric(ci_high),
      p.value = as.numeric(p),
      n = nrow(dat)
    )
    
  } else {
    # ridge backend with bootstrap inference
    
    # ---- SPECIAL CASE: only 1 predictor -> use unregularised weighted fit
    if (length(x_names) == 1) {
      x1 <- x_names[1]
      
      # interpret scalar lower.limits for this single coefficient
      lower1 <- if (length(lower.limits) == 1) lower.limits else lower.limits[1]
      
      f1 <- wls_single(dat, x1, lower = lower1)
      boot <- bootstrap_wls_single(dat, x1, B = ridge_B, lower = lower1, seed = ridge_seed)
      
      se  <- sd(boot, na.rm = TRUE)
      ci  <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
      
      # Same centred two-sided bootstrap test used for multivariable fits.
      z <- boot[is.finite(boot)]
      p <- if (!length(z) || !is.finite(se) || se == 0) NA_real_
      else centred_bootstrap_p(f1$slope, z)
      
      return(tibble(
        vv_stage = y_stage,
        pd_stage = sub("^x_", "", x1),
        estimate = f1$slope,
        std.error = se,
        conf.low = as.numeric(ci[1]),
        conf.high = as.numeric(ci[2]),
        p.value = p,
        n = nrow(dat)
      ))
    }
    
    # ---- regular ridge path for >=2 predictors
    fr <- fit_ridge_cv(
      dat, x_names,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      nfolds = ridge_nfolds,
      standardize = ridge_standardize
    )
    
    boot_mat <- bootstrap_ridge(
      dat, x_names,
      B = ridge_B,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      nfolds = ridge_nfolds,
      standardize = ridge_standardize,
      refit_cv = ridge_refit_cv,
      seed = ridge_seed
    )
    
    summ <- summarise_boot(fr$beta, boot_mat)
    
    tibble(
      vv_stage = y_stage,
      pd_stage = sub("^x_", "", x_names)
    ) %>%
      bind_cols(summ) %>%
      mutate(n = nrow(dat))
  }
}

run_mult_grid <- function(wide, stages,
                          universe = c("all_orthologs","de_either","de_both","vv_de_any","pls_top"),
                          alpha = 0.05, pls_ids = NULL,
                          lfc = c("raw","shrunk"),
                          se  = c("auto","raw","shrunk"),
                          predictors = c("all","older_or_equal"),
                          scale_X = TRUE,
                          p_adj = c("bonferroni","BH"),
                          # NEW
                          model = c("lmrob","ridge"),
                          # ridge options
                          lower.limits = -Inf,
                          upper.limits = Inf,
                          ridge_B = 200,
                          ridge_nfolds = 10,
                          ridge_standardize = FALSE,
                          ridge_refit_cv = TRUE,
                          ridge_seed = NULL) {
  
  universe   <- match.arg(universe)
  lfc        <- match.arg(lfc)
  se         <- match.arg(se)
  predictors <- match.arg(predictors)
  p_adj      <- match.arg(p_adj)
  model      <- match.arg(model)
  
  wide_use <- pick_universe(wide, universe = universe, alpha = alpha, pls_ids = pls_ids)
  
  # run stage-wise models
  res <- purrr::map_dfr(stages, function(ys) {
    dat <- build_mult_data2(
      wide_use, stages = stages, y_stage = ys,
      lfc = lfc, se = se, predictors = predictors, scale_X = scale_X
    )
    
    fit_mult_stage2(
      dat, y_stage = ys,
      model = model,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      ridge_B = ridge_B,
      ridge_nfolds = ridge_nfolds,
      ridge_standardize = ridge_standardize,
      ridge_refit_cv = ridge_refit_cv,
      ridge_seed = ridge_seed
    )
  }) %>%
    mutate(
      vv_stage = factor(vv_stage, levels = stages),
      pd_stage = factor(pd_stage, levels = stages)
    )
  
  # FULL vv × pd template
  template <- tidyr::expand_grid(
    vv_stage = factor(stages, levels = stages),
    pd_stage = factor(stages, levels = stages)
  ) %>%
    mutate(
      allowed = if (predictors == "all") TRUE else
        as.integer(pd_stage) >= as.integer(vv_stage)
    )
  
  # join results onto full grid
  out <- template %>%
    left_join(res, by = c("vv_stage","pd_stage"))
  
  # set non-allowed cells explicitly to NA
  out <- out %>%
    mutate(
      estimate  = ifelse(allowed, estimate, NA_real_),
      std.error = ifelse(allowed, std.error, NA_real_),
      conf.low  = ifelse(allowed, conf.low, NA_real_),
      conf.high = ifelse(allowed, conf.high, NA_real_),
      p.value   = ifelse(allowed, p.value, NA_real_),
      n         = ifelse(allowed, n, NA_real_)
    )
  
  # p.adjust within each Vv stage, only for allowed + non-NA p's
  out <- out %>%
    group_by(vv_stage) %>%
    mutate(
      p.adj = {
        p <- p.value
        adj <- rep(NA_real_, length(p))
        idx <- which(allowed & !is.na(p))
        if (length(idx)) adj[idx] <- stats::p.adjust(p[idx], method = p_adj)
        adj
      }
    ) %>%
    ungroup() %>%
    mutate(
      lag = as.integer(vv_stage) - as.integer(pd_stage),
      universe = universe, lfc = lfc, se = se,
      predictors = predictors, model = model
    )
  
  out
}

# Run analysis
# choices: "all_orthologs","de_either","de_both","vv_de_any","pls_top"
universe_choice <- "vv_de_any"
alpha_universe  <- 0.05
lfc_choice      <- "shrunk"      # or "raw"
se_choice       <- "raw"         # or "auto"/"shrunk"
pred_policy     <- "older_or_equal"  # or "all"
padj_method     <- "bonferroni"
scale_X         <- FALSE

# # robust regression
# grid_mult <- run_mult_grid(
#   wide       = ortholog_de_wide,
#   stages     = stages,
#   universe   = universe_choice,
#   alpha      = alpha_universe,
#   pls_ids    = c(sel_sc_pos_topK, sel_sc_neg_topK),
#   lfc        = lfc_choice,
#   se         = se_choice,
#   predictors = pred_policy,
#   scale_X    = scale_X,
#   p_adj      = padj_method,
#   model      = "lmrob"
# )
# 
# # ridge regression
# grid_mult <- run_mult_grid(
#   wide       = ortholog_de_wide,
#   stages     = stages,
#   universe   = universe_choice,
#   alpha      = alpha_universe,
#   pls_ids    = c(sel_sc_pos_topK, sel_sc_neg_topK),
#   lfc        = lfc_choice,
#   se         = se_choice,
#   predictors = pred_policy,
#   scale_X    = scale_X,
#   p_adj      = padj_method,
#   model      = "ridge",
#   ridge_B    = 500,
#   ridge_nfolds = 10
# )

# ridge regression with nonnegativity constraints
grid_mult <- run_mult_grid(
  wide       = ortholog_de_wide,
  stages     = stages,
  universe   = universe_choice,
  alpha      = alpha_universe,
  pls_ids    = c(sel_sc_pos_topK, sel_sc_neg_topK),
  lfc        = lfc_choice,
  se         = se_choice,
  predictors = pred_policy,
  scale_X    = scale_X,
  p_adj      = padj_method,
  model      = "ridge",
  lower.limits = 0,          # <-- nonnegative coefficients
  ridge_B    = 500
)

# Heatmap of partial slopes
p_heterochrony <- ggplot(grid_mult, aes(x = pd_stage, y = vv_stage, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = gtools::stars.pval(p.adj)), size = 5) +
  scale_fill_gradient2("partial\nslope", low = "#92C2FFFF", mid = "white", high = "red2",
                       midpoint = 0, na.value = "grey90",
                       # oob     = scales::squish,
                       # limits = c(0, 0.7)
  ) +
  labs(x = expression(italic("Polistes dominula") ~ "stage (predictors)"),
       y = expression(italic("Vespula vulgaris") ~ "stage (response)"),
       fill = "Partial slope"
  ) +
  theme_minimal(base_size = 12)
graph2png(x = p_heterochrony, file =file.path(figure_dir, "Fig2B_source.png"), width=6, height=4)
graph2pdf(x = p_heterochrony, file =file.path(figure_dir, "Fig2B_source.pdf"), width=6, height=4)
graph2ppt(x = p_heterochrony, file =file.path(figure_dir, "Fig2B_source.pptx"), width=6, height=4)

# Main-text Figure 2: editable vector robust-regression and heterochrony panels.
# Wasp photographs are intentionally omitted so they can be placed manually.
p_stagewise_regression_figure <- p_stagewise_regression +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(12, 12, 8, 18)
  )

p_heterochrony_figure <- p_heterochrony +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 12, 12, 18)
  )

figure2_plot <- patchwork::wrap_plots(
  p_stagewise_regression_figure,
  p_heterochrony_figure,
  ncol = 1,
  heights = c(1.35, 1)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        # graph2ppt applies a 0.6 text scaling factor; 40 ggplot points gives
        # a 24-point PowerPoint text box. The OOXML correction below then
        # enforces the exact requested font properties in the exported deck.
        family = "Arial", face = "bold", size = 40, colour = "black"
      ),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

# graph2ppt/svglite currently drops the bold flag on patchwork tags. Correct
# only the A/B tag runs in the exported PPTX while preserving all editable
# vector elements and the slide master/layout.
enforce_pptx_panel_label_style <- function(
    file, labels = c("A", "B"), font_family = "Arial", font_size = 24) {
  stopifnot(file.exists(file))
  unpack_dir <- tempfile("figure2_pptx_")
  dir.create(unpack_dir)
  on.exit(unlink(unpack_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(file, exdir = unpack_dir)

  slide_files <- list.files(
    file.path(unpack_dir, "ppt", "slides"),
    pattern = "^slide[0-9]+\\.xml$", full.names = TRUE
  )
  matched <- character()

  for (slide_file in slide_files) {
    doc <- xml2::read_xml(slide_file)
    ns <- xml2::xml_ns(doc)
    text_nodes <- xml2::xml_find_all(doc, ".//a:t", ns)
    target_nodes <- text_nodes[xml2::xml_text(text_nodes) %in% labels]

    for (node in target_nodes) {
      label <- xml2::xml_text(node)
      run_properties <- xml2::xml_find_first(
        node, "parent::a:r/a:rPr", ns
      )
      if (inherits(run_properties, "xml_missing")) next
      xml2::xml_set_attrs(
        run_properties,
        c(sz = as.character(round(font_size * 100)), b = "1", i = "0")
      )
      for (font_node_name in c("a:latin", "a:cs")) {
        font_node <- xml2::xml_find_first(
          run_properties, paste0("./", font_node_name), ns
        )
        if (inherits(font_node, "xml_missing")) {
          font_node <- xml2::xml_add_child(
            run_properties, sub("^a:", "", font_node_name),
            .where = 1, .namespace = unname(ns[["a"]])
          )
        }
        xml2::xml_set_attr(font_node, "typeface", font_family)
      }
      matched <- c(matched, label)
    }
    xml2::write_xml(doc, slide_file)
  }

  if (!setequal(unique(matched), labels)) {
    stop(
      "Could not identify all requested PowerPoint panel labels: ",
      paste(setdiff(labels, unique(matched)), collapse = ", ")
    )
  }

  rebuilt <- tempfile(fileext = ".pptx")
  archive_files <- list.files(
    unpack_dir, recursive = TRUE, all.files = TRUE,
    include.dirs = FALSE, no.. = TRUE
  )
  zip::zipr(
    zipfile = rebuilt, files = archive_files, root = unpack_dir,
    include_directories = FALSE, mode = "mirror"
  )
  officer::read_pptx(rebuilt)
  if (!file.copy(rebuilt, file, overwrite = TRUE)) {
    stop("Could not replace PPTX after correcting panel-label typography.")
  }
  invisible(file)
}

figure2_file <- file.path(figure_dir, "Fig2.pptx")
graph2png(
  x = figure2_plot,
  file = file.path(figure_dir, "Fig2.png"),
  width = 7.5, height = 10
)
graph2pdf(
  x = figure2_plot,
  file = file.path(figure_dir, "Fig2.pdf"),
  width = 7.5, height = 10
)
graph2ppt(
  x = figure2_plot,
  file = figure2_file,
  width = 7.5, height = 10
)
enforce_pptx_panel_label_style(
  figure2_file, labels = c("A", "B"),
  font_family = "Arial", font_size = 24
)

# -> l2, L3 and L4 sign. predict differential expression in Q vs W in L2 in Vv 
# consistent with caste-biased gene expression in Vespula being 
# moved forward in development to the L2 stage, i.e. with heterochrony 
# (more specifically, predisplacement). 

# -> The positive P vs P fits the “shared maturation/overwintering” module being 
# reused (fat body buildup, stress tolerance, storage).



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
  supplemental_table_dir, "TablesS4-S8.xlsx"
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
  X_variance_explained_percent = axis_variance,
  HOGs_in_final_model = ncol(Z_sub)
)

supp_table_s5b <- as_tibble(as.data.frame(sc_contr_oos)) %>%
  transmute(
    Species = recode(as.character(Species),
                     Pd = "Polistes dominula", Vv = "Vespula vulgaris"),
    Stage = as.character(Stage),
    Contrast = as.character(contrast),
    Estimate = estimate,
    SE = SE,
    df = df,
    t_ratio = t.ratio,
    P_value = p.value
  ) %>%
  arrange(factor(Species, c("Polistes dominula", "Vespula vulgaris")),
          factor(Stage, stages))

category_labels <- c(
  "Differential feeding / nutrient metabolism" = "Nutrient metabolism",
  "Hibernation / diapause / long-term survival" = "Hibernation",
  "Reproduction / caste-related development" = "Reproduction",
  "Other" = "Other"
)
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
    Expected_HOGs = Expected,
    Fold_enrichment = FoldEnrichment,
    weight01_P = p,
    Contributing_HOGs = MembersStr
  ) %>%
  arrange(factor(Category, c("Nutrient metabolism", "Hibernation",
                             "Reproduction", "Other")), weight01_P)

supp_table_s7 <- slopes_by_stage %>%
  transmute(
    Stage = as.character(stage),
    HOGs = n,
    Robust_slope = estimate,
    SE = std.error,
    CI_95_lower = conf.low,
    CI_95_upper = conf.high,
    P_value = p.value,
    Bonferroni_P = p.value.adj
  ) %>%
  arrange(factor(Stage, stages))

supp_table_s8 <- grid_mult %>%
  filter(allowed) %>%
  transmute(
    Vespula_response_stage = as.character(vv_stage),
    Polistes_predictor_stage = as.character(pd_stage),
    Stage_lag = lag,
    HOGs = n,
    Partial_slope = estimate,
    Bootstrap_SE = std.error,
    CI_95_lower = conf.low,
    CI_95_upper = conf.high,
    Centred_bootstrap_P = p.value,
    Bonferroni_P = p.adj
  ) %>%
  arrange(factor(Vespula_response_stage, stages),
          factor(Polistes_predictor_stage, stages))

supp_table_index <- tibble(
  Table = c("S4", "S5A", "S5B", "S6", "S7", "S8"),
  Worksheet = c(
    "S4_DE_summary", "S5A_PLS_axes", "S5B_PLS_LOOCV",
    "S6_Fig1_GO", "S7_Fig2A_regressions", "S8_Fig2B_ridge"
  ),
  Description = c(
    "Differential-expression counts by species, analysis unit and stage",
    "PLS-axis mapping and variance explained",
    "Leave-one-sample-out season/caste contrasts on PLS axis 1",
    "GO terms displayed in Fig. 1B",
    "Stagewise robust-regression results underlying Fig. 2A",
    "Nonnegative-ridge partial slopes underlying Fig. 2B"
  )
)

supplementary_tables <- list(
  Table_index = supp_table_index,
  S4_DE_summary = supp_table_s4,
  S5A_PLS_axes = supp_table_s5a,
  S5B_PLS_LOOCV = supp_table_s5b,
  S6_Fig1_GO = supp_table_s6,
  S7_Fig2A_regressions = supp_table_s7,
  S8_Fig2B_ridge = supp_table_s8
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
  TableS5A = supp_table_s5a,
  TableS5B = supp_table_s5b,
  TableS6 = supp_table_s6,
  TableS7 = supp_table_s7,
  TableS8 = supp_table_s8
)
purrr::iwalk(
  numbered_supplementary_tables,
  ~ readr::write_tsv(.x, file.path(supplemental_table_dir, paste0(.y, ".tsv")))
)

find_python <- function() {
  candidates <- c(
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    Sys.which("python"),
    Sys.which("python3"),
    file.path(
      Sys.getenv("USERPROFILE"), ".cache", "codex-runtimes",
      "codex-primary-runtime", "dependencies", "python", "python.exe"
    )
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
    "PLS_features", "PLS_X_variance_by_component",
    "GO_annotation", "GO_foreground", "GO_test", "GO_multiple_testing",
    "ridge_bootstrap_test", "ridge_multiple_testing",
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
    paste(signif(plsfit$prop_expl_var$X, 6), collapse = ";"),
    "direct EXCON/Galaxy EggNOG and InterProScan conserved across wasps plus experimental FlyBase GO transferred through N13 Drosophila orthologues",
    "top 300 N13 HOGs in each oriented PLS loading direction",
    "topGO weight01 Fisher; nodeSize 20",
    "none; nominal topGO p values exported",
    "500 pairs-bootstrap replicates; two-sided centred-bootstrap P value with plus-one finite-sample correction",
    "Bonferroni correction across allowed Polistes predictors within each Vespula response stage",
    repository_relative_path(n13_orthology_summary_file),
    repository_relative_path(full_gene_de_file),
    repository_relative_path(full_hog_de_file),
    repository_relative_path(figure1_file),
    repository_relative_path(figure2_file),
    repository_relative_path(supplementary_xlsx)
  )
)
write_tsv(analysis_metadata_final, file.path(output_dir, "analysis_metadata.tsv"))
write_tsv(slopes_by_stage, file.path(output_dir, "stagewise_robust_regression_slopes.tsv"))
write_tsv(grid_mult, file.path(output_dir, "distributed_lag_nonnegative_ridge_results.tsv"))
key_object_names <- intersect(
  c(
    "gene_de_pd", "gene_de_vv", "hog_de_shared",
    "n13_composition", "n13_orthology_composition_summary",
    "ortholog_de_wide", "de_gene_categories", "de_stacked_counts",
    "de_stacked_counts_plot",
    "scores_export", "orthologs_pls", "pls_go_lists",
    "go_term_results_all", "go_term_results_enriched", "flat_all",
    "pls_foregrounds", "tbl", "p_go_sc_pos", "p_pls_figure",
    "slopes_by_stage", "p_stagewise_regression",
    "grid_mult", "p_heterochrony", "figure2_plot",
    "supp_table_s4", "supp_table_s5a", "supp_table_s5b",
    "supp_table_s6", "supp_table_s7", "supp_table_s8",
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
  finite_nonzero_bootstrap_p = all(
    is.na(grid_mult$p.value) |
      (grid_mult$p.value > 0 & grid_mult$p.value <= 1)
  ),
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
  supplementary_GO_names_complete = !any(str_detect(
    supp_table_s6$GO_term, fixed("...")
  )),
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

