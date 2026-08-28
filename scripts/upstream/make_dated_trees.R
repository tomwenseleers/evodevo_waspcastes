# Date the EXCON/OrthoFinder topology for Vespidae and Aculeata.
#
# Branch lengths in the output trees are millions of years (Ma). The input
# topology is retained; only branch lengths are estimated. Run this file with
# source("scripts/upstream/make_dated_trees.R") or Rscript.

library(ape)

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
    for (up in 0:2) {
      probe <- candidate
      if (up > 0) probe <- normalizePath(
        file.path(candidate, paste(rep("..", up), collapse = "/")),
        winslash = "/", mustWork = FALSE
      )
      if (file.exists(file.path(probe, "A_annotation.R")) &&
          dir.exists(file.path(probe, "input_CAFE"))) return(probe)
    }
  }
  stop("Could not locate the evodevo_waspcastes repository root.")
}

if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  PROJECT_ROOT <- resolve_repository_root()
}

make_dated_trees <- function(
    analysis_dir = PROJECT_ROOT,
    input_tree = file.path(
      analysis_dir,
      "input_CAFE/tree_dating/orthofinder_pruned_tree_scaled.nwk"
    ),
    output_dir = file.path(analysis_dir, "input_CAFE", "tree_dating"),
    input_scale = 1000,
    chronos_model = "discrete",
    lambda = 1,
    rate_categories = 20) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(input_tree), input_scale > 0)

  repository_relative_path <- function(path, root = analysis_dir) {
    path_abs <- normalizePath(path, winslash = "/", mustWork = TRUE)
    root_abs <- normalizePath(root, winslash = "/", mustWork = TRUE)
    prefix <- paste0(root_abs, "/")
    if (startsWith(tolower(path_abs), tolower(prefix))) {
      return(substring(path_abs, nchar(prefix) + 1L))
    }
    basename(path_abs)
  }

  full_tree <- read.tree(input_tree)
  if (!is.rooted(full_tree) || is.null(full_tree$edge.length)) {
    stop("The input must be a rooted tree with branch lengths.")
  }
  if (any(!is.finite(full_tree$edge.length)) ||
      any(full_tree$edge.length <= 0)) {
    stop("All input branch lengths must be finite and positive.")
  }

  # EXCON's pruned_tree was written after multiplying lengths by 1000.
  full_tree$edge.length <- full_tree$edge.length / input_scale
  full_tree$node.label <- NULL

  clade_tips <- function(tree, tips) {
    tips <- unique(tips)
    missing <- setdiff(tips, tree$tip.label)
    if (length(missing)) {
      stop("Tips absent from tree: ", paste(missing, collapse = ", "))
    }
    extract.clade(tree, getMRCA(tree, tips))$tip.label
  }

  require_exact_clade <- function(tree, tips, label) {
    observed <- clade_tips(tree, tips)
    if (!setequal(observed, tips)) {
      stop(
        label, " is not monophyletic in the input tree. Extra tips: ",
        paste(setdiff(observed, tips), collapse = ", ")
      )
    }
    invisible(observed)
  }

  vespidae <- grep(
    "^(Ancistrocerus|Mischocyttarus|Polistes|Vespa|Vespula|Dolichovespula)_",
    full_tree$tip.label,
    value = TRUE
  )
  social_vespidae <- grep(
    "^(Mischocyttarus|Polistes|Vespa|Vespula|Dolichovespula)_",
    full_tree$tip.label,
    value = TRUE
  )
  vespidae_with_outgroups <- unique(c(
    vespidae,
    "Anoplius_nigerrimus",
    "Tiphia_femorata"
  ))
  formicidae <- grep(
    "^(Harpegnathos|Ooceraea|Linepithema|Atta|Acromyrmex|Monomorium|Solenopsis|Camponotus)_",
    full_tree$tip.label,
    value = TRUE
  )

  # The seed spans the two ends of Apoidea in this particular topology.
  apoidea <- clade_tips(
    full_tree,
    c("Nysson_spinosus", "Melipona_bicolor")
  )
  aculeata <- unique(c(
    vespidae, formicidae, apoidea,
    "Anoplius_nigerrimus", "Tiphia_femorata"
  ))

  require_exact_clade(full_tree, vespidae, "Sampled Vespidae")
  require_exact_clade(
    full_tree,
    vespidae_with_outgroups,
    "Sampled Vespidae + Anoplius + Tiphia"
  )
  require_exact_clade(full_tree, social_vespidae, "Polistinae + Vespinae")
  require_exact_clade(full_tree, formicidae, "Formicidae")
  require_exact_clade(full_tree, apoidea, "Apoidea")
  require_exact_clade(
    full_tree, c(formicidae, apoidea), "Formicidae + Apoidea"
  )
  require_exact_clade(full_tree, aculeata, "Aculeata")

  source_blaimer <- paste0(
    "Blaimer et al. 2023, Nature Communications 14:1212; ",
    "Supplementary Data 6, preferred topC-1 chronogram"
  )
  url_blaimer <- "https://doi.org/10.1038/s41467-023-36868-4"
  source_luo <- paste0(
    "Luo et al. 2022, BMC Ecology and Evolution 22:63; ",
    "Paleovespa (60-64 Ma), Symmorphus (90-94 Ma), ",
    "Polistinae-Vespinae estimate about 75 Ma"
  )
  url_luo <- "https://doi.org/10.1186/s12862-022-02017-6"

  # Point estimates are fixed in the primary chronograms so that shared nodes
  # have the same age in both trees. Blaimer intervals are published 95% HPDs.
  # The Vespidae ranges are conservative sensitivity bounds.
  specs <- list(
    Vespidae = data.frame(
      calibration = c("sampled_Vespidae", "Polistinae_plus_Vespinae"),
      fixed_age_Ma = c(110.0, 75.0),
      age_min_Ma = c(90.0, 60.0),
      age_max_Ma = c(123.9, 80.0),
      source = c(
        paste(source_luo, source_blaimer, sep = "; "),
        source_luo
      ),
      source_url = c(
        paste(url_luo, url_blaimer, sep = "; "),
        url_luo
      ),
      rationale = c(
        paste(
          "Fixed to Luo et al.'s Eumeninae-social Vespidae estimate;",
          "the sensitivity range combines the Eumeninae fossil and",
          "the upper 95% HPD for Vespoidea in Blaimer et al."
        ),
        paste(
          "Fixed to Luo et al.'s Polistinae-Vespinae estimate;",
          "the sensitivity range includes the oldest Vespinae fossil."
        )
      ),
      stringsAsFactors = FALSE
    ),
    Aculeata = data.frame(
      calibration = c(
        "Aculeata",
        "Formicidae_plus_Apoidea",
        "Apoidea",
        "sampled_Vespidae",
        "Polistinae_plus_Vespinae"
      ),
      fixed_age_Ma = c(142.30, 108.60, 102.13, 110.0, 75.0),
      age_min_Ma = c(129.1, 98.5, 92.1, 90.0, 60.0),
      age_max_Ma = c(157.6, 120.9, 113.8, 123.9, 80.0),
      source = c(
        source_blaimer, source_blaimer, source_blaimer,
        paste(source_luo, source_blaimer, sep = "; "), source_luo
      ),
      source_url = c(
        url_blaimer, url_blaimer, url_blaimer,
        paste(url_luo, url_blaimer, sep = "; "), url_luo
      ),
      rationale = c(
        "Fixed to the published median; range is its 95% HPD.",
        "Fixed to the published median; range is its 95% HPD.",
        "Fixed to the published median; range is its 95% HPD.",
        paste(
          "Fixed to Luo et al.'s Eumeninae-social Vespidae estimate;",
          "the sensitivity range combines the Eumeninae fossil and",
          "the upper 95% HPD for Vespoidea in Blaimer et al."
        ),
        paste(
          "Fixed to Luo et al.'s Polistinae-Vespinae estimate;",
          "the sensitivity range includes the oldest Vespinae fossil."
        )
      ),
      stringsAsFactors = FALSE
    )
  )

  clades <- list(
    Vespidae = vespidae,
    sampled_Vespidae = vespidae,
    Polistinae_plus_Vespinae = social_vespidae,
    Aculeata = aculeata,
    Formicidae_plus_Apoidea = c(formicidae, apoidea),
    Apoidea = apoidea
  )

  node_age <- function(tree, node) {
    depths <- node.depth.edgelength(tree)
    max(depths[seq_len(Ntip(tree))]) - depths[node]
  }

  date_one_tree <- function(tree_name, tips, spec, dated_source = NULL) {
    fit_warnings <- character()
    if (is.null(dated_source)) {
      tree <- keep.tip(full_tree, tips)
      tree$node.label <- NULL

      calibration_nodes <- vapply(
        spec$calibration,
        function(label) {
          getMRCA(tree, intersect(clades[[label]], tree$tip.label))
        },
        integer(1)
      )
      calibration <- makeChronosCalib(
        tree,
        node = calibration_nodes,
        age.min = spec$fixed_age_Ma,
        age.max = spec$fixed_age_Ma,
        soft.bounds = FALSE
      )

      fit <- withCallingHandlers(
        chronos(
          tree,
          lambda = lambda,
          model = chronos_model,
          calibration = calibration,
          control = chronos.control(
            nb.rate.cat = rate_categories,
            iter.max = 20000,
            eval.max = 20000,
            dual.iter.max = 50
          ),
          quiet = TRUE
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      dating_method <- "chronos"
      derived_from <- ""
    } else {
      missing <- setdiff(tips, dated_source$tip.label)
      if (length(missing)) {
        stop(
          tree_name, " tips absent from dated source: ",
          paste(missing, collapse = ", ")
        )
      }
      fit <- keep.tip(dated_source, tips)
      fit$node.label <- NULL
      dating_method <- "pruned_dated_tree"
      derived_from <- "Aculeata_dated_primary"
    }

    if (length(fit_warnings)) {
      stop(
        tree_name, " dating generated warnings: ",
        paste(unique(fit_warnings), collapse = " | ")
      )
    }
    if (!is.ultrametric(fit, tol = 1e-7)) {
      stop(tree_name, " result is not ultrametric.")
    }
    if (any(!is.finite(fit$edge.length)) || any(fit$edge.length <= 0)) {
      stop(tree_name, " result contains invalid branch lengths.")
    }

    fitted_nodes <- vapply(
      spec$calibration,
      function(label) getMRCA(fit, intersect(clades[[label]], fit$tip.label)),
      integer(1)
    )
    fitted_ages <- vapply(fitted_nodes, function(x) node_age(fit, x), numeric(1))
    tolerance <- 1e-5
    if (any(abs(fitted_ages - spec$fixed_age_Ma) > tolerance)) {
      stop(tree_name, " result does not reproduce a fixed calibration age.")
    }
    within_bounds <- fitted_ages >= spec$age_min_Ma - tolerance &
      fitted_ages <= spec$age_max_Ma + tolerance
    if (!all(within_bounds)) {
      stop(tree_name, " result violates one or more calibration bounds.")
    }

    prefix <- file.path(output_dir, paste0(tree_name, "_dated_primary"))
    write.tree(fit, paste0(prefix, ".nwk"))
    write.nexus(fit, file = paste0(prefix, ".nex"))

    pdf(
      paste0(prefix, ".pdf"),
      width = 12,
      height = max(8, 0.23 * Ntip(fit) + 2),
      onefile = FALSE
    )
    par(mar = c(4, 1, 3, 1))
    plot(
      ladderize(fit),
      direction = "rightwards",
      use.edge.length = TRUE,
      show.tip.label = TRUE,
      cex = if (Ntip(fit) > 30) 0.55 else 0.8,
      label.offset = max(node.depth.edgelength(fit)) * 0.01,
      no.margin = FALSE
    )
    axisPhylo(backward = TRUE)
    title(
      main = paste(gsub("_", " ", tree_name), "dated tree"),
      sub = "Branch lengths and axis are millions of years before present"
    )
    dev.off()

    depths <- node.depth.edgelength(fit)
    tip_depths <- depths[seq_len(Ntip(fit))]
    qc <- data.frame(
      tree = tree_name,
      input_file = repository_relative_path(input_tree),
      input_scale_removed = input_scale,
      n_tips = Ntip(fit),
      n_internal_nodes = fit$Nnode,
      dating_method = dating_method,
      derived_from = derived_from,
      model = chronos_model,
      lambda = lambda,
      rate_categories = rate_categories,
      root_age_Ma = max(tip_depths),
      min_root_to_tip_Ma = min(tip_depths),
      max_root_to_tip_Ma = max(tip_depths),
      min_edge_Ma = min(fit$edge.length),
      ultrametric = is.ultrametric(fit, tol = 1e-7),
      binary = is.binary(fit),
      warnings = "",
      stringsAsFactors = FALSE
    )
    fitted_spec <- transform(
      spec,
      tree = tree_name,
      node = fitted_nodes,
      fitted_age_Ma = fitted_ages,
      matches_fixed_age = abs(fitted_ages - fixed_age_Ma) <= tolerance,
      within_sensitivity_bounds = within_bounds
    )
    fitted_spec <- fitted_spec[, c(
      "tree", "calibration", "node", "fixed_age_Ma", "fitted_age_Ma",
      "matches_fixed_age", "age_min_Ma", "age_max_Ma",
      "within_sensitivity_bounds", "source", "source_url", "rationale"
    )]
    tip_table <- data.frame(
      tree = tree_name,
      tip = fit$tip.label,
      stringsAsFactors = FALSE
    )

    list(tree = fit, qc = qc, calibrations = fitted_spec, tips = tip_table)
  }

  aculeata_result <- date_one_tree(
    "Aculeata", aculeata, specs$Aculeata
  )
  results <- list(
    Vespidae = date_one_tree("Vespidae", vespidae, specs$Vespidae),
    Vespidae_with_outgroups = date_one_tree(
      "Vespidae_with_outgroups",
      vespidae_with_outgroups,
      specs$Vespidae,
      dated_source = aculeata_result$tree
    ),
    Aculeata = aculeata_result
  )

  write.table(
    do.call(rbind, lapply(results, `[[`, "calibrations")),
    file.path(output_dir, "dating_calibrations.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    do.call(rbind, lapply(results, `[[`, "qc")),
    file.path(output_dir, "dating_qc.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    do.call(rbind, lapply(results, `[[`, "tips")),
    file.path(output_dir, "dated_tree_tips.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  message("Wrote dated trees and metadata to: ", normalizePath(output_dir))
  invisible(results)
}

dated_tree_results <- make_dated_trees()
