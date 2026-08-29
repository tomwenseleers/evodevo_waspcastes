# Final differential-expression and cross-species transcriptomic workflow.

# FALSE loads the compact key-results RDS. TRUE reruns the downstream analyses
# from the versioned gene- and N13-HOG checkpoints. It does not refit glmmTMB.
rerun <- FALSE
rerun_mixed_models <- FALSE

resolve_repository_root <- function() {
  candidates <- character()
  command_file <- sub("^--file=", "", grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  ))
  if (length(command_file)) candidates <- c(candidates, dirname(command_file[[1L]]))
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile) && nzchar(frame$ofile)) {
      candidates <- c(candidates, dirname(frame$ofile))
    }
  }
  candidates <- c(candidates, getwd())
  for (candidate in unique(candidates)) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(candidate, "A_annotation.R")) &&
        dir.exists(file.path(candidate, "input_differential_expression")) &&
        dir.exists(file.path(candidate, "nextflow_runs"))) return(candidate)
  }
  stop("Could not locate the evodevo_waspcastes repository root.")
}

PROJECT_ROOT <- resolve_repository_root()
env_rerun <- Sys.getenv("EVODEVO_RERUN", unset = "")
if (nzchar(env_rerun)) rerun <- tolower(env_rerun) %in% c("1", "true", "yes")

if (isTRUE(rerun_mixed_models)) {
  stop(
    "The expensive glmmTMB fits are deliberately separate. Run ",
    "scripts/upstream/B1_gene_glmmTMB_raw_counts_length_offset.R and then ",
    "scripts/upstream/B1_N13_HOG_glmmTMB_raw_counts_length_offset.R with ",
    "PROJECT_ROOT set to this repository. Use one worker per species."
  )
}

if (isTRUE(rerun)) {
  source(file.path(
    PROJECT_ROOT, "scripts", "pipelines", "B_differential_expression_pipeline.R"
  ))
} else {
  workspace <- file.path(
    PROJECT_ROOT, "input_differential_expression", "workspaces",
    "B_differential_expression_key_results.rds"
  )
  if (!file.exists(workspace)) stop("Missing differential-expression workspace: ", workspace)
  list2env(readRDS(workspace), envir = .GlobalEnv)
  message(
    "Loaded the final differential-expression key objects. Set rerun <- TRUE ",
    "to regenerate the derived tables and figures from the saved model checkpoints."
  )
}
