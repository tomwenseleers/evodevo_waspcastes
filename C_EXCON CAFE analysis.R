# Final CAFE and node-wise GO-enrichment workflow for the wasp phylogeny.

# FALSE loads the compact key-results RDS. TRUE reruns the complete downstream
# CAFE summarisation, HOG annotation, TE filtering and node-wise topGO analyses.
rerun <- FALSE

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
        dir.exists(file.path(candidate, "input_CAFE")) &&
        dir.exists(file.path(candidate, "nextflow_runs"))) return(candidate)
  }
  stop("Could not locate the evodevo_waspcastes repository root.")
}

PROJECT_ROOT <- resolve_repository_root()
env_rerun <- Sys.getenv("EVODEVO_RERUN", unset = "")
if (nzchar(env_rerun)) rerun <- tolower(env_rerun) %in% c("1", "true", "yes")

if (isTRUE(rerun)) {
  options(evodevo.cafe.auto_run = TRUE)
  source(file.path(PROJECT_ROOT, "scripts", "pipelines", "C_EXCON_CAFE_pipeline.R"))
} else {
  workspace <- file.path(
    PROJECT_ROOT, "input_CAFE", "workspaces", "C_EXCON_CAFE_key_results.rds"
  )
  if (!file.exists(workspace)) stop("Missing CAFE workspace: ", workspace)
  list2env(readRDS(workspace), envir = .GlobalEnv)
  message(
    "Loaded the final CAFE/GO key objects. Set rerun <- TRUE to regenerate ",
    "the CAFE summaries, node-wise GO tables and Fig3."
  )
}
