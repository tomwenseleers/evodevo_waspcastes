# Final functional-annotation workflow for the wasp caste evo-devo study.

# FALSE loads the versioned workspace and existing compact annotation exports.
# TRUE reconstructs all annotation outputs from input_annotation/raw.
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
    if (file.exists(file.path(candidate, "A_annotation_FINAL.R")) &&
        dir.exists(file.path(candidate, "input_annotation"))) return(candidate)
  }
  stop("Could not locate the evodevo_waspcastes repository root.")
}

PROJECT_ROOT <- resolve_repository_root()
env_rerun <- Sys.getenv("EVODEVO_RERUN", unset = "")
if (nzchar(env_rerun)) rerun <- tolower(env_rerun) %in% c("1", "true", "yes")

if (isTRUE(rerun)) {
  source(file.path(PROJECT_ROOT, "scripts", "pipelines", "A_annotation_pipeline.R"))
} else {
  workspace <- file.path(
    PROJECT_ROOT, "input_annotation", "workspaces", "A_annotation_workspace.RData"
  )
  if (!file.exists(workspace)) stop("Missing annotation workspace: ", workspace)
  load(workspace, envir = .GlobalEnv)
  message(
    "Loaded the final annotation workspace. Set rerun <- TRUE to rebuild ",
    "output_annotation from the raw input files."
  )
}
