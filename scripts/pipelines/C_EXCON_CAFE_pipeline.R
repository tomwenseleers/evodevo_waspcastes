# C_EXCON CAFE analysis_FINAL.R
#
# Self-contained final workflow for the wasp-focused EXCON/CAFE analysis.
# It summarises and annotates significant expansions and contractions at Nodes
# 20 and 25, excludes conservatively annotated transposable-element HOGs from
# downstream totals and interpretation, audits root-zero removals, performs
# node-wise topGO enrichment among significantly expanded non-TE HOGs, and
# exports the annotated tables and editable phylogeny PowerPoint.
#
# Final thresholds and enrichment settings:
#   CAFE family-wide p < 0.05 and focal Viterbi branch probability < 0.01;
#   topGO weight01/Fisher for BP, MF and CC, without p-value adjustment;
#   top_nodes_per_ontology = Inf.
#
# Source this complete file from RStudio to run the analysis. To load its
# functions without executing the workflow, first run:
#   options(evodevo.cafe.auto_run = FALSE)
#
# All inputs and outputs are resolved relative to this script's location in
# the repository root; the R working directory is not used to locate files.

resolve_repository_root <- function() {
  candidates <- character()

  command_file <- sub(
    "^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  )
  if (length(command_file)) {
    candidates <- c(candidates, dirname(command_file[[1L]]))
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    source_file <- frames[[i]]$ofile
    if (!is.null(source_file) && length(source_file) && nzchar(source_file)) {
      candidates <- c(candidates, dirname(source_file))
    }
  }
  candidates <- c(candidates, getwd())

  for (candidate in unique(candidates)) {
    candidate <- normalizePath(
      candidate, winslash = "/", mustWork = FALSE
    )
    if (
      file.exists(file.path(candidate, "A_annotation_FINAL.R")) &&
      dir.exists(file.path(candidate, "input_CAFE"))
    ) {
      return(candidate)
    }
  }
  stop(
    "Could not locate the evodevo_waspcastes repository root."
  )
}

if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  PROJECT_ROOT <- resolve_repository_root()
}
REPOSITORY_ROOT <- normalizePath(
  PROJECT_ROOT, winslash = "/", mustWork = TRUE
)

required_packages <- c(
  "data.table", "topGO", "GO.db", "AnnotationDbi", "ape", "export"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Required R packages are missing: ",
    paste(missing_packages, collapse = ", ")
  )
}
suppressPackageStartupMessages({
  library(data.table)
  library(topGO)
})

make_wasp_annotation_config <- function(project_dir, excon_dir, run_dir) {
  list(
    n0_file = file.path(run_dir, "cafe", "base", "N0.tsv"),
    vv_genetable_file = file.path(
      excon_dir, "annotation", "vv_genetable.csv"
    ),
    pd_genetable_file = file.path(
      excon_dir, "annotation", "pd_genetable.csv"
    ),
    dmel_genetable_file = file.path(
      excon_dir, "annotation", "dmel_genetable.csv"
    ),
    pd_vv_annotated_ortholog_file = file.path(
      excon_dir, "annotation", "pd_vv_orthofinder_orthologs_annotated.csv"
    ),
    pd_eggnog_annotation_file = file.path(
      excon_dir, "annotation", "Polistes_dominula_EggNOG.tsv"
    ),
    vv_eggnog_annotation_file = file.path(
      excon_dir, "annotation", "Vespula_vulgaris_EggNOG.tsv"
    ),
    ancistrocerus_ortholog_file = file.path(
      excon_dir, "orthology", "Ancistrocerus_nigricornis.clean.tsv"
    )
  )
}

# EMBEDDED GO-MAPPING AND topGO FUNCTIONS -----------------------------

assert_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path)
}

normalize_hog_go <- function(gene2go) {
  gene2go <- gene2go[lengths(gene2go) > 0L]
  lapply(gene2go, function(x) sort(unique(x[grepl("^GO:[0-9]+$", x)])))
}

load_og_go <- function(output_dir) {
  rds_path <- file.path(output_dir, "OG_GO_format_gene2GO.rds")
  tsv_path <- file.path(output_dir, "OG_GO_format.tsv")

  if (file.exists(rds_path)) {
    return(normalize_hog_go(readRDS(rds_path)))
  }
  assert_file(tsv_path, "OG_GO_format.tsv")
  message("Reading HOG-to-GO mapping: ", tsv_path)
  map <- fread(
    tsv_path,
    sep = "\t",
    header = FALSE,
    col.names = c("HOG", "GO"),
    colClasses = "character",
    quote = ""
  )
  normalize_hog_go(split(map$GO, map$HOG))
}

numeric_topgo_p <- function(x) {
  x <- trimws(as.character(x))
  suppressWarnings(as.numeric(sub("^<\\s*", "", x)))
}

write_go_member_table <- function(results, selected, path) {
  if (!nrow(results)) {
    empty <- data.table(
      HOG = character(),
      change = numeric(),
      branch_p = numeric(),
      GO.ID = character(),
      Term = character(),
      ontology = character(),
      p_raw = numeric(),
      topGOresult = character(),
      Significant = numeric(),
      Expected = numeric(),
      FoldEnrichment = numeric()
    )
    fwrite(empty, path, sep = "\t", quote = FALSE)
    return(invisible(empty))
  }

  member_rows <- results[
    !is.na(Members) & nzchar(Members),
    .(
      HOG = unlist(strsplit(Members, ";", fixed = TRUE), use.names = FALSE)
    ),
    by = .(
      GO.ID, Term, ontology, p_raw, topGOresult, Significant, Expected,
      FoldEnrichment
    )
  ]
  member_rows <- merge(
    member_rows,
    selected[, .(HOG, change, branch_p)],
    by = "HOG",
    all.x = TRUE
  )
  setcolorder(
    member_rows,
    c(
      "GO.ID", "Term", "ontology", "HOG", "change", "branch_p", "p_raw",
      "topGOresult", "Significant", "Expected", "FoldEnrichment"
    )
  )
  setorder(member_rows, ontology, p_raw, GO.ID, HOG)
  fwrite(member_rows, path, sep = "\t", quote = FALSE)
  invisible(member_rows)
}

run_topgo_test <- function(
    foreground, universe, gene2go, ontology = "BP",
    algorithm = "weight01", statistic = "fisher", node_size = 1L,
    top_nodes = 50L) {
  ontology <- match.arg(ontology, c("BP", "MF", "CC"))
  algorithm <- match.arg(
    algorithm, c("classic", "elim", "weight", "weight01", "lea", "parentchild")
  )
  if (statistic != "fisher") {
    stop("This binary foreground workflow currently supports statistic='fisher' only.")
  }

  universe <- sort(intersect(unique(universe), names(gene2go)))
  foreground <- sort(intersect(unique(foreground), universe))
  local_gene2go <- gene2go[universe]

  all_genes <- factor(as.integer(universe %in% foreground))
  names(all_genes) <- universe

  go_data <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = all_genes,
    geneSel = function(x) x == 1L,
    annot = annFUN.gene2GO,
    gene2GO = local_gene2go,
    nodeSize = as.integer(node_size)
  )
  result <- runTest(go_data, algorithm = algorithm, statistic = statistic)

  available <- length(score(result))
  n_to_report <- if (is.infinite(top_nodes)) available else min(as.integer(top_nodes), available)
  if (n_to_report < 1L) return(data.table())

  tab <- as.data.table(GenTable(
    go_data,
    topGOresult = result,
    orderBy = "topGOresult",
    ranksOf = "topGOresult",
    topNodes = n_to_report,
    numChar = 10000L
  ))
  full_go_terms <- vapply(tab$GO.ID, function(go_id) {
    term_object <- GO.db::GOTERM[[go_id]]
    if (is.null(term_object)) return(NA_character_)
    AnnotationDbi::Term(term_object)
  }, character(1L))
  tab[!is.na(full_go_terms), Term := full_go_terms[!is.na(full_go_terms)]]
  if (any(grepl("\\.\\.\\.$|…$", tab$Term))) {
    stop("At least one GO term remained ellipsis-truncated after GO.db lookup.")
  }
  tab[, p_raw := numeric_topgo_p(topGOresult)]

  adjustment_methods <- c(
    "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  )
  for (method in adjustment_methods) {
    tab[, (method) := p.adjust(p_raw, method = method)]
  }

  tab[, FoldEnrichment := fifelse(Expected > 0, Significant / Expected, Inf)]
  tab[, ontology := ontology]
  tab[, algorithm := algorithm]
  tab[, statistic := statistic]

  foreground_members <- names(all_genes)[all_genes == 1L]
  term_members <- genesInTerm(go_data, tab$GO.ID)
  members <- lapply(term_members, intersect, y = foreground_members)
  tab[, Members_n := lengths(members)]
  tab[, Members := vapply(members, paste, collapse = ";", FUN.VALUE = character(1))]

  setcolorder(
    tab,
    c(
      "GO.ID", "Term", "Annotated", "Significant", "Expected",
      "FoldEnrichment", "topGOresult", "p_raw", adjustment_methods,
      "ontology", "algorithm", "statistic", "Members_n", "Members"
    )
  )
  tab[]
}

# EMBEDDED HOG-ANNOTATION FUNCTIONS ----------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

split_ids <- function(x) {
  x <- x %||% ""
  x <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  unique(x[nzchar(x) & !is.na(x)])
}

collapse_values <- function(x, separator = "; ") {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) paste(x, collapse = separator) else ""
}

normalize_description <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^LOW QUALITY PROTEIN:\\s*", "", x, ignore.case = TRUE)
  x <- gsub("\\s+Add an annotation\\s*$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+isoform\\s+[A-Za-z0-9.-]+$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+", " ", x)
  bad <- is.na(x) |
    !nzchar(x) |
    grepl("^uncharacterized", x, ignore.case = TRUE) |
    grepl("^hypothetical protein", x, ignore.case = TRUE) |
    x %in% c("-", "NA", "n/a")
  x[bad] <- NA_character_
  x
}

read_eggnog_gene_annotations <- function(path, species, priority = 2L) {
  annotations <- fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    na.strings = c("", "-", "NA")
  )
  query_col <- names(annotations)[1L]
  query <- as.character(annotations[[query_col]])
  gene_id <- sub(".*\\|gene=([^|]+)\\|.*", "\\1", query)
  gene_id[gene_id == query] <- NA_character_

  data.table(
    species = species,
    gene_id = gene_id,
    gene_symbol = as.character(annotations[["Preferred_name"]]),
    gene_name = as.character(annotations[["Description"]]),
    annotation_source = "eggNOG-mapper exact-gene-ID fallback",
    annotation_priority = priority
  )[!is.na(gene_id) & nzchar(gene_id)]
}

first_informative <- function(x) {
  x <- as.character(x)
  hit <- which(!is.na(x) & nzchar(trimws(x)) & x != "-")
  if (length(hit)) x[hit[1L]] else NA_character_
}

read_direct_annotations <- function(config) {
  vv <- fread(config$vv_genetable_file, na.strings = c("", "NA"))
  vv <- vv[, .(
    species = "Vespula_vulgaris",
    gene_id = as.character(gene_vv),
    gene_symbol = as.character(symbol_vv),
    gene_name = as.character(description_vv),
    annotation_source = "NCBI: vv_genetable.csv",
    annotation_priority = 1L
  )]

  pd <- fread(config$pd_genetable_file, na.strings = c("", "NA"))
  pd <- pd[, .(
    species = "Polistes_dominula",
    gene_id = as.character(gene_pd),
    gene_symbol = as.character(symbol_pd),
    gene_name = as.character(description_pd),
    annotation_source = "NCBI: pd_genetable.csv",
    annotation_priority = 1L
  )]

  eggnog <- rbindlist(list(
    read_eggnog_gene_annotations(
      config$vv_eggnog_annotation_file,
      "Vespula_vulgaris"
    ),
    read_eggnog_gene_annotations(
      config$pd_eggnog_annotation_file,
      "Polistes_dominula"
    )
  ))

  # The orthology table is used only as a gene-ID annotation fallback. Its
  # run-specific OG numbers must not be joined to the current N0 HOG numbers.
  fallback_raw <- fread(
    config$pd_vv_annotated_ortholog_file,
    na.strings = c("", "NA")
  )
  fallback_rows <- vector("list", 2L * nrow(fallback_raw))
  k <- 0L
  for (i in seq_len(nrow(fallback_raw))) {
    for (prefix in c("vv", "pd")) {
      genes <- split_ids(fallback_raw[[paste0("gene_", prefix)]][i])
      if (!length(genes)) next
      k <- k + 1L
      fallback_rows[[k]] <- data.table(
        species = if (prefix == "vv") {
          "Vespula_vulgaris"
        } else {
          "Polistes_dominula"
        },
        gene_id = genes,
        gene_symbol = as.character(
          fallback_raw[[paste0("symbol_", prefix)]][i]
        ),
        gene_name = as.character(
          fallback_raw[[paste0("description_", prefix)]][i]
        ),
        annotation_source =
          "pd_vv_orthofinder_orthologs_annotated.csv fallback",
        annotation_priority = 3L
      )
    }
  }
  fallback_rows <- fallback_rows[seq_len(k)]
  fallback <- if (length(fallback_rows)) {
    rbindlist(fallback_rows, fill = TRUE)
  } else {
    data.table()
  }

  annotations <- rbindlist(list(vv, pd, eggnog, fallback), fill = TRUE)
  annotations[, gene_name := normalize_description(gene_name)]
  annotations[, gene_symbol := fifelse(
    is.na(gene_symbol) | gene_symbol %in% c("", "NA"),
    NA_character_,
    gene_symbol
  )]
  setorder(annotations, species, gene_id, annotation_priority)
  annotations[, {
    symbol_index <- which(
      !is.na(gene_symbol) & nzchar(gene_symbol) & gene_symbol != "-"
    )
    name_index <- which(!is.na(gene_name) & nzchar(gene_name))
    used_indices <- unique(c(
      head(symbol_index, 1L),
      head(name_index, 1L)
    ))
    list(
      gene_symbol = first_informative(gene_symbol),
      gene_name = first_informative(gene_name),
      annotation_source = collapse_values(
        annotation_source[used_indices]
      )
    )
  }, by = .(species, gene_id)]
}

read_dmel_annotations <- function(config) {
  dmel <- fread(config$dmel_genetable_file, na.strings = c("", "NA"))
  base <- dmel[, .(
    gene_id = as.character(gene_dmel),
    gene_symbol = as.character(symbol_dmel),
    gene_name = normalize_description(description_dmel)
  )]

  by_gene <- copy(base)
  by_gene[, lookup_id := gene_id]
  by_symbol <- copy(base)
  by_symbol[, lookup_id := gene_symbol]
  by_prefixed_symbol <- copy(base)
  by_prefixed_symbol[, lookup_id := paste0("Dmel_", gene_symbol)]

  result <- rbindlist(
    list(by_gene, by_symbol, by_prefixed_symbol),
    use.names = TRUE
  )
  result <- result[
    !is.na(lookup_id) & nzchar(lookup_id),
    .(lookup_id, gene_id, gene_symbol, gene_name)
  ]
  unique(result, by = "lookup_id")
}

make_member_long <- function(n0_selected, direct_annotations) {
  species <- c(
    "Vespula_vulgaris",
    "Polistes_dominula",
    "Ancistrocerus_nigricornis"
  )
  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(n0_selected))) {
    for (sp in species) {
      ids <- split_ids(n0_selected[[sp]][i])
      if (!length(ids)) next
      k <- k + 1L
      rows[[k]] <- data.table(
        HOG = n0_selected$HOG[i],
        OG = n0_selected$OG[i],
        species = sp,
        gene_id = ids
      )
    }
  }
  members <- rbindlist(rows[seq_len(k)], fill = TRUE)
  members <- merge(
    members,
    direct_annotations[, .(
      species, gene_id, gene_symbol, gene_name, annotation_source
    )],
    by = c("species", "gene_id"),
    all.x = TRUE,
    sort = FALSE
  )
  members[, annotation_quality := fifelse(
    !is.na(gene_name) | !is.na(gene_symbol),
    "direct NCBI-derived annotation",
    "not directly annotated"
  )]
  members[]
}

read_ancistrocerus_inferred_annotations <- function(
    config, ancistrocerus_gene_ids, direct_annotations, dmel_annotations) {
  if (!length(ancistrocerus_gene_ids)) return(data.table())

  orthologs <- fread(
    config$ancistrocerus_ortholog_file,
    na.strings = c("", "NA")
  )
  query_col <- "Ancistrocerus_nigricornis.clean"
  wanted_species <- c(
    "Vespula_vulgaris.clean",
    "Polistes_dominula.clean",
    "Drosophila_melanogaster.clean"
  )
  orthologs <- orthologs[Species %in% wanted_species]

  pair_rows <- list()
  k <- 0L
  for (i in seq_len(nrow(orthologs))) {
    query_ids <- intersect(
      split_ids(orthologs[[query_col]][i]),
      ancistrocerus_gene_ids
    )
    target_ids <- split_ids(orthologs$Orthologs[i])
    if (!length(query_ids) || !length(target_ids)) next
    k <- k + 1L
    pair_rows[[k]] <- CJ(
      gene_id = query_ids,
      ortholog_gene_id = target_ids,
      unique = TRUE
    )[, target_species := sub("\\.clean$", "", orthologs$Species[i])]
  }
  if (!k) return(data.table())
  pairs <- unique(rbindlist(pair_rows[seq_len(k)]))

  direct_target <- direct_annotations[, .(
    target_species = species,
    ortholog_gene_id = gene_id,
    transferred_symbol = gene_symbol,
    transferred_name = gene_name,
    transferred_from = species
  )]
  dmel_target <- dmel_annotations[, .(
    target_species = "Drosophila_melanogaster",
    ortholog_gene_id = lookup_id,
    transferred_symbol = gene_symbol,
    transferred_name = gene_name,
    transferred_from = "Drosophila_melanogaster"
  )]
  target_annotations <- rbindlist(
    list(direct_target, dmel_target),
    fill = TRUE
  )
  pairs <- merge(
    pairs,
    target_annotations,
    by = c("target_species", "ortholog_gene_id"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  pairs <- pairs[
    !is.na(transferred_name) | !is.na(transferred_symbol)
  ]
  if (!nrow(pairs)) return(data.table())

  pairs[, .(
    gene_symbol = collapse_values(transferred_symbol),
    gene_name = collapse_values(transferred_name),
    annotation_source = paste0(
      "OrthoFinder orthology transfer from: ",
      collapse_values(transferred_from)
    )
  ), by = gene_id]
}

annotate_ancistrocerus_members <- function(members, inferred) {
  if (!nrow(inferred)) return(members)
  idx <- members$species == "Ancistrocerus_nigricornis"
  members[, member_row_id := .I]
  anc <- merge(
    members[idx, .(member_row_id, gene_id)],
    inferred,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  matched <- !is.na(anc$gene_name) | !is.na(anc$gene_symbol)
  if (any(matched)) {
    target_rows <- anc$member_row_id[matched]
    members[target_rows, `:=`(
      gene_symbol = anc$gene_symbol[matched],
      gene_name = anc$gene_name[matched],
      annotation_source = anc$annotation_source[matched],
      annotation_quality = "orthology-inferred; not a direct NCBI name"
    )]
  }
  members[, member_row_id := NULL]
  members[]
}

build_dmel_hog_evidence <- function(n0_selected, dmel_annotations) {
  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(n0_selected))) {
    ids <- split_ids(n0_selected$Drosophila_melanogaster[i])
    if (!length(ids)) next
    k <- k + 1L
    rows[[k]] <- data.table(
      HOG = n0_selected$HOG[i],
      lookup_id = ids
    )
  }
  if (!k) return(data.table())
  result <- merge(
    rbindlist(rows[seq_len(k)]),
    dmel_annotations,
    by = "lookup_id",
    all.x = TRUE
  )
  result[, `:=`(
    species = "Drosophila_melanogaster",
    evidence_source = "FlyBase-derived dmel_genetable.csv",
    evidence_weight = 1
  )]
  result[]
}

consensus_rule <- function(text) {
  rules <- list(
    "odorant receptor family" =
      "odorant receptor|olfactory receptor",
    "gustatory receptor family" =
      "gustatory receptor",
    "odorant-binding protein family" =
      "odorant[- ]binding protein",
    "ionotropic/glutamate receptor family" =
      "ionotropic receptor|ionotropic glutamate receptor",
    "chemosensory protein family" =
      "chemosensory protein",
    "maltase/glucosidase family" =
      "maltase|alpha[- ]glucosidase|glucosidase",
    "histone H2B-like family" =
      "histone h2b",
    "fatty acid synthase" =
      "fatty acid synthase",
    "myelin transcription factor 1-like" =
      "myelin transcription factor 1|\\bmyt1\\b",
    "cytochrome P450 family" =
      "cytochrome p450",
    "cuticular protein family" =
      "cuticular protein|cuticle protein",
    "ribosomal protein family" =
      "ribosomal protein",
    "ubiquitin-related family" =
      "ubiquitin",
    "histone H3-like family" =
      "histone h3",
    "histone H4-like family" =
      "histone h4",
    "chitinase/chitin-binding family" =
      "chitinase|chitin[- ]binding",
    "serine protease family" =
      "serine protease|trypsin",
    "zinc-finger protein family" =
      "zinc[- ]finger",
    "ATP synthase family" =
      "atp synthase",
    "dynein family" =
      "dynein",
    "tubulin family" =
      "\\btubulin\\b",
    "actin family" =
      "\\bactin\\b",
    "transposable-element/reverse-transcriptase family" =
      "reverse transcriptase|transpos|retrotrans"
  )
  for (label in names(rules)) {
    if (grepl(rules[[label]], text, ignore.case = TRUE, perl = TRUE)) {
      return(label)
    }
  }
  NA_character_
}

choose_consensus <- function(names, symbols, sources) {
  names <- normalize_description(names)
  keep <- !is.na(names)
  names <- names[keep]
  sources <- sources[keep]
  symbols <- trimws(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  symbols <- symbols[!grepl(
    "^(LOC[0-9]+|ENS[A-Z0-9]*G[0-9]+)$",
    symbols,
    ignore.case = TRUE
  )]

  all_text <- tolower(paste(c(names, symbols), collapse = " | "))
  ruled <- consensus_rule(all_text)
  if (!is.na(ruled)) {
    return(list(
      consensus_gene_name = ruled,
      consensus_method = "conservative keyword family rule",
      consensus_evidence = collapse_values(head(unique(names), 8L))
    ))
  }

  if (length(names)) {
    weights <- ifelse(
      grepl("NCBI", sources, fixed = TRUE), 3,
      ifelse(grepl("FlyBase", sources, fixed = TRUE), 2, 1)
    )
    scored <- data.table(
      original = names,
      name_key = tolower(names),
      weight = weights
    )[, .(
      score = sum(weight),
      occurrences = .N,
      representative = original[which.min(nchar(original))]
    ), by = name_key]
    scored[, representative_length := nchar(representative)]
    setorder(scored, -score, -occurrences, representative_length)
    selected <- head(scored$representative, 2L)
    return(list(
      consensus_gene_name = paste(selected, collapse = " / "),
      consensus_method = "weighted recurring annotation",
      consensus_evidence = collapse_values(head(unique(names), 8L))
    ))
  }

  if (length(symbols)) {
    return(list(
      consensus_gene_name = paste(head(unique(symbols), 3L), collapse = " / "),
      consensus_method = "symbol fallback",
      consensus_evidence = collapse_values(head(unique(symbols), 8L))
    ))
  }

  list(
    consensus_gene_name = "unannotated HOG",
    consensus_method = "no usable local annotation",
    consensus_evidence = ""
  )
}

build_consensus_annotations <- function(members, dmel_evidence, all_hogs) {
  member_evidence <- members[
    !is.na(gene_name) | !is.na(gene_symbol),
    .(
      HOG,
      gene_name,
      gene_symbol,
      evidence_source = annotation_source
    )
  ]
  dmel_evidence <- dmel_evidence[, .(
    HOG,
    gene_name,
    gene_symbol,
    evidence_source
  )]
  evidence <- rbindlist(
    list(member_evidence, dmel_evidence),
    fill = TRUE
  )

  rbindlist(lapply(all_hogs, function(hog) {
    x <- evidence[HOG == hog]
    consensus <- choose_consensus(
      x$gene_name, x$gene_symbol, x$evidence_source
    )
    data.table(
      HOG = hog,
      consensus_gene_name = consensus$consensus_gene_name,
      consensus_method = consensus$consensus_method,
      consensus_evidence = consensus$consensus_evidence
    )
  }))
}

build_member_summary <- function(members, all_hogs) {
  summarized <- members[, .(
    member_count = uniqueN(gene_id),
    gene_ids = collapse_values(gene_id),
    gene_symbols = collapse_values(gene_symbol),
    gene_names = collapse_values(gene_name),
    annotation_sources = collapse_values(annotation_source),
    annotated_member_count = sum(
      !is.na(gene_name) | !is.na(gene_symbol)
    )
  ), by = .(HOG, species)]

  result <- data.table(HOG = all_hogs)
  for (sp in c(
    "Vespula_vulgaris",
    "Polistes_dominula",
    "Ancistrocerus_nigricornis"
  )) {
    part <- summarized[species == sp]
    setnames(
      part,
      c(
        "member_count", "gene_ids", "gene_symbols", "gene_names",
        "annotation_sources", "annotated_member_count"
      ),
      paste0(
        sp,
        c(
          "_member_count", "_gene_ids", "_gene_symbols", "_gene_names",
          "_annotation_sources", "_annotated_member_count"
        )
      )
    )
    part[, species := NULL]
    result <- merge(result, part, by = "HOG", all.x = TRUE)
  }
  count_cols <- grep(
    "(_member_count|_annotated_member_count)$",
    names(result),
    value = TRUE
  )
  for (column in count_cols) set(result, which(is.na(result[[column]])), column, 0L)
  text_cols <- setdiff(names(result), c("HOG", count_cols))
  for (column in text_cols) set(result, which(is.na(result[[column]])), column, "")
  result[]
}

add_interest_flags <- function(summary) {
  # Vectorized row-wise text is assembled separately to avoid cross-row mixing.
  text <- tolower(apply(summary[, .(
    consensus_gene_name,
    consensus_evidence,
    Vespula_vulgaris_gene_names,
    Polistes_dominula_gene_names,
    Ancistrocerus_nigricornis_gene_names
  )], 1L, paste, collapse = " | "))
  summary[, chemosensory_family_flag := grepl(
    paste(
      "odorant receptor|olfactory receptor|gustatory receptor|",
      "odorant[- ]binding protein|chemosensory protein|",
      "ionotropic receptor",
      sep = ""
    ),
    text,
    ignore.case = TRUE,
    perl = TRUE
  )]
  summary[]
}

# FINAL WASP CAFE ANALYSIS -------------------------------------------

# Summarise and annotate focal-node CAFE changes in EXCON run 14, then test
# GO enrichment among significantly expanded HOGs at nodes 20 and 25.
#
# This is a scripted R/RStudio workflow. Source the complete file from the
# repository root; no command-line arguments are parsed. Set
# options(evodevo.cafe.auto_run = FALSE) before sourcing to load the function
# without immediately executing it.

suppressPackageStartupMessages({
  library(data.table)
})

plot_wasp_tree <- function(
  tree_file,
  node_summary,
  figure_category_summary
) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required to plot the wasp phylogeny.")
  }

  tree <- ape::read.tree(tree_file)
  focal_nodes <- c("20", "25")
  focal_summary <- copy(node_summary)
  focal_summary[, node := as.character(node)]
  focal_summary <- focal_summary[node %in% focal_nodes]
  focal_summary[, node_order := match(node, focal_nodes)]
  setorder(focal_summary, node_order)
  focal_summary[, node_order := NULL]
  if (!identical(focal_summary$node, focal_nodes)) {
    stop("The node summary must contain wasp-tree nodes 20 and 25.")
  }

  expected_descendants <- list(
    `20` = c(
      "Mischocyttarus_mexicanus", "Polistes_dominula", "Polistes_fuscatus",
      "Polistes_exclamans", "Polistes_canadensis", "Vespula_germanica",
      "Vespula_pensylvanica", "Vespula_vulgaris", "Vespa_crabro",
      "Vespa_mandarinia", "Vespa_velutina", "Dolichovespula_media",
      "Dolichovespula_saxonica", "Dolichovespula_sylvestris"
    ),
    `25` = c(
      "Vespula_germanica", "Vespula_pensylvanica", "Vespula_vulgaris",
      "Vespa_crabro", "Vespa_mandarinia", "Vespa_velutina",
      "Dolichovespula_media", "Dolichovespula_saxonica",
      "Dolichovespula_sylvestris"
    )
  )
  for (node_id in focal_nodes) {
    observed <- sort(ape::extract.clade(tree, as.integer(node_id))$tip.label)
    if (!identical(observed, sort(expected_descendants[[node_id]]))) {
      stop("Unexpected descendant species for wasp-tree node ", node_id, ".")
    }
  }

  solitary <- c(
    "Anoplius_nigerrimus", "Tiphia_femorata", "Ancistrocerus_nigricornis"
  )
  polistinae <- c(
    "Mischocyttarus_mexicanus", "Polistes_dominula", "Polistes_fuscatus",
    "Polistes_exclamans", "Polistes_canadensis"
  )
  vespinae <- expected_descendants[["25"]]

  desired_bottom_to_top <- c(solitary, polistinae, vespinae)
  tree <- ape::rotateConstr(tree, desired_bottom_to_top)

  tip_colours <- ifelse(
    tree$tip.label %in% polistinae, "#306B80",
    ifelse(tree$tip.label %in% vespinae, "#B54737", "#505A63")
  )
  names(tip_colours) <- tree$tip.label
  display_tip_labels <- setNames(
    gsub("_", " ", tree$tip.label, fixed = TRUE),
    tree$tip.label
  )

  category_colours <- c(
    "Lipid / CHC" = "#F2AA00",
    "Chemoperception" = "#0AA36F",
    "Digestion / nutrition" = "#E56600",
    "Endocrine" = "#087DB6",
    "Behavioural signal" = "#087DB6",
    "Defence / venom" = "#A6498D",
    "Other / unresolved" = "#BAC0C5"
  )

  draw_vector_stroke <- function(x1, y1, x2, y2, colour, lwd = 1) {
    segments(x1, y1, x2, y2, col = colour, lwd = lwd)
  }

  draw_category_panel <- function(node_id, xleft, xright, border_colour) {
    rows <- copy(figure_category_summary[node == as.integer(node_id)])
    rows <- rows[n_HOGs > 0]
    rows[, category_order := match(figure_category, names(category_colours))]
    setorder(rows, category_order)
    total <- sum(rows$n_HOGs)
    if (!total) stop("No figure categories were supplied for node ", node_id, ".")

    panel_bottom <- 9.26
    panel_top <- 11.36
    rect(
      xleft, panel_bottom, xright, panel_top,
      col = "#FCFBFC", border = border_colour, lwd = 1
    )
    accent_colour <- if (node_id == "20") "#80639A" else "#B8615E"
    rect(
      xleft, panel_top - 0.05, xright, panel_top,
      col = accent_colour, border = NA
    )

    panel_title <- if (node_id == "20") {
      "STEM POLISTINAE + VESPINAE"
    } else {
      "STEM VESPINAE"
    }
    text(
      xleft + 0.39, panel_top - 0.38, panel_title,
      adj = c(0, 0.5), cex = 1.55, font = 2, col = "#263238"
    )
    text(
      xleft + 0.39, panel_top - 0.74,
      sprintf("Node %s | %d significantly expanded HOGs", node_id, total),
      adj = c(0, 0.5), cex = 1.15, col = "#5E6A73"
    )

    bar_left <- xleft + 0.28
    bar_right <- xright - 0.28
    bar_bottom <- panel_top - 1.18
    bar_top <- bar_bottom + 0.23
    cursor <- bar_left
    for (i in seq_len(nrow(rows))) {
      next_cursor <- cursor + (bar_right - bar_left) * rows$n_HOGs[[i]] / total
      rect(
        cursor, bar_bottom, next_cursor, bar_top,
        col = category_colours[[rows$figure_category[[i]]]], border = NA
      )
      cursor <- next_cursor
    }

    legend_cols <- 3L
    legend_y <- c(panel_bottom + 0.68, panel_bottom + 0.30)
    legend_x <- seq(xleft + 0.38, xright - 2.48, length.out = legend_cols)
    for (i in seq_len(nrow(rows))) {
      col_id <- (i - 1L) %% legend_cols + 1L
      row_id <- (i - 1L) %/% legend_cols + 1L
      if (row_id > length(legend_y)) next
      x <- legend_x[[col_id]]
      y <- legend_y[[row_id]]
      rect(
        x - 0.10, y - 0.08, x + 0.06, y + 0.08,
        col = category_colours[[rows$figure_category[[i]]]], border = NA
      )
      text(
        x + 0.26, y,
        paste0(rows$figure_category[[i]], "  ", rows$n_HOGs[[i]]),
        adj = c(0, 0.5), cex = 1.00, col = "#263238"
      )
    }
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = rep(0, 4), xpd = NA, family = "sans", bg = "white")
  plot.new()
  plot.window(xlim = c(0, 16.5), ylim = c(0, 11.7), xaxs = "i", yaxs = "i")
  n_tip <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  root_node <- setdiff(unique(tree$edge[, 1L]), tree$edge[, 2L])[[1L]]
  children_by_node <- split(tree$edge[, 2L], tree$edge[, 1L])

  ordered_tips <- function(node) {
    if (node <= n_tip) return(node)
    unlist(lapply(children_by_node[[as.character(node)]], ordered_tips))
  }
  displayed_tip_order <- ordered_tips(root_node)

  tree_bottom <- 1.78
  tree_top <- 8.90
  tip_y_values <- seq(tree_bottom, tree_top, length.out = n_tip)
  y_coord <- numeric(n_tip + n_node)
  y_coord[displayed_tip_order] <- tip_y_values
  calculate_node_y <- function(node) {
    if (node <= n_tip) return(y_coord[[node]])
    child_y <- vapply(
      children_by_node[[as.character(node)]], calculate_node_y, numeric(1L)
    )
    y_coord[[node]] <<- mean(range(child_y))
    y_coord[[node]]
  }
  calculate_node_y(root_node)

  node_depth <- ape::node.depth.edgelength(tree)
  tree_height <- max(node_depth[seq_len(n_tip)])
  age_scale <- ceiling(tree_height / 20) * 20
  tree_left <- 2.68
  tree_right <- 10.78
  age_to_x <- function(age) {
    tree_right - age / age_scale * (tree_right - tree_left)
  }
  x_coord <- age_to_x(tree_height - node_depth)
  tip_y <- setNames(y_coord[seq_len(n_tip)], tree$tip.label)

  group_ranges <- list(
    vespinae = range(tip_y[vespinae]) + c(-0.23, 0.23),
    polistinae = range(tip_y[polistinae]) + c(-0.23, 0.23),
    solitary = range(tip_y[solitary]) + c(-0.23, 0.23)
  )
  rect(
    0.35, group_ranges$vespinae[[1L]], 16.15, group_ranges$vespinae[[2L]],
    col = "#F7E9E7", border = NA
  )
  rect(
    0.35, group_ranges$polistinae[[1L]], 16.15, group_ranges$polistinae[[2L]],
    col = "#DCECF2", border = NA
  )
  rect(
    0.35, group_ranges$solitary[[1L]], 16.15, group_ranges$solitary[[2L]],
    col = "#F1F1F1", border = NA
  )

  draw_category_panel("20", 0.36, 8.17, "#8066A5")
  draw_category_panel("25", 8.35, 16.16, "#C5675B")

  edge_colour <- "#737373"
  phylogeny_lwd <- 1
  for (parent in as.integer(names(children_by_node))) {
    children <- children_by_node[[as.character(parent)]]
    child_y <- y_coord[children]
    draw_vector_stroke(
      x_coord[[parent]], min(child_y),
      x_coord[[parent]], max(child_y),
      colour = edge_colour, lwd = phylogeny_lwd
    )
  }
  for (edge_index in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[edge_index, 1L]
    child <- tree$edge[edge_index, 2L]
    draw_vector_stroke(
      x_coord[[parent]], y_coord[[child]],
      x_coord[[child]], y_coord[[child]],
      colour = edge_colour, lwd = phylogeny_lwd
    )
  }

  node_colours <- c("20" = "#80639A", "25" = "#B8615E")
  for (node_id in focal_nodes) {
    node_number <- as.integer(node_id)
    symbols(
      x_coord[[node_number]], y_coord[[node_number]],
      circles = 0.06, inches = FALSE, add = TRUE,
      fg = node_colours[[node_id]], bg = node_colours[[node_id]]
    )
  }

  tip_label_x <- tree_right + 0.10
  for (tip_index in seq_len(n_tip)) {
    tip_name <- tree$tip.label[[tip_index]]
    text(
      tip_label_x, y_coord[[tip_index]], display_tip_labels[[tip_name]],
      adj = c(0, 0.5), cex = 4 / 3, font = 3,
      col = tip_colours[[tip_name]]
    )
  }

  text(
    0.66, mean(group_ranges$vespinae), "ADVANCED\nEUSOCIAL",
    adj = c(0, 0.5), cex = 19 / 12, col = "#111111"
  )
  text(
    0.66, mean(group_ranges$polistinae), "PRIMITIVELY\nEUSOCIAL",
    adj = c(0, 0.5), cex = 19 / 12, col = "#111111"
  )
  text(
    0.66, mean(group_ranges$solitary), "SOLITARY",
    adj = c(0, 0.5), cex = 19 / 12, col = "#111111"
  )

  note_xleft <- 2.74
  note_xright <- 7.84
  note_ybottom <- 7.72
  note_ytop <- 8.51
  rect(
    note_xleft, note_ybottom, note_xright, note_ytop,
    col = "#EFECEC", border = "#263238", lwd = 1
  )
  text(
    note_xleft + 0.13, mean(c(note_ybottom, note_ytop)),
    "Significant expansions (+) / contractions (-)\nCAFE family p < 0.05; branch p < 0.01\nTE-derived HOGs excluded",
    adj = c(0, 0.5), cex = 1.24, col = "#59656E"
  )

  node_backgrounds <- c("20" = "#DED7E5", "25" = "#F3DDD8")
  node_offsets <- list(
    `20` = c(x = -0.59, y = 0.34),
    `25` = c(x = -0.655, y = 0.34)
  )
  node_box_widths <- c("20" = 0.89, "25" = 1.00)
  node_label_y_offsets <- c("20" = 0.175, "25" = 0.145)
  for (node_id in focal_nodes) {
    row <- focal_summary[node == node_id]
    node_number <- as.integer(node_id)
    node_text <- sprintf(
      "+%d / -%d",
      row$n_significant_expansions,
      row$n_significant_contractions
    )
    box_x <- x_coord[[node_number]] + node_offsets[[node_id]][["x"]]
    box_y <- y_coord[[node_number]] + node_offsets[[node_id]][["y"]]
    box_width <- node_box_widths[[node_id]]
    box_height <- 0.38
    rect(
      box_x - box_width / 2, box_y - box_height / 2,
      box_x + box_width / 2, box_y + box_height / 2,
      col = node_backgrounds[[node_id]], border = "#263238", lwd = 1
    )
    text(box_x, box_y, node_text, cex = 4 / 3, col = "#17242B")
    text(
      box_x, y_coord[[node_number]] - node_label_y_offsets[[node_id]],
      paste("Node", node_id), adj = c(0.5, 1), cex = 4 / 3,
      col = node_colours[[node_id]]
    )
  }

  axis_ages <- seq(0, age_scale, by = 20)
  axis_x <- age_to_x(axis_ages)
  axis_y <- 1.13
  draw_vector_stroke(
    min(axis_x), axis_y, max(axis_x), axis_y,
    colour = "#263238", lwd = 1
  )
  for (tick_x in axis_x) {
    draw_vector_stroke(
      tick_x, axis_y, tick_x, axis_y - 0.13,
      colour = "#263238", lwd = 1
    )
  }
  text(
    axis_x, axis_y - 0.32, labels = axis_ages,
    cex = 7 / 6, col = "#263238"
  )
  text(
    mean(range(axis_x)), axis_y - 0.69,
    "Millions of years before present",
    cex = 4 / 3, font = 1, col = "#263238"
  )

  invisible(tree)
}

run_wasp_cafe_analysis <- function(
  project_dir = REPOSITORY_ROOT,
  output_root = NULL,
  family_p_cutoff = 0.05,
  branch_p_cutoff = 0.01,
  go_p_cutoff = 0.05,
  algorithm = "weight01",
  statistic = "fisher",
  ontologies = c("BP", "MF", "CC"),
  top_nodes_per_ontology = Inf
) {
  project_dir <- normalizePath(
    project_dir, winslash = "/", mustWork = TRUE
  )

excon_dir <- file.path(project_dir, "input_CAFE")
mapping_dir <- file.path(excon_dir, "GO")
if (is.null(output_root)) {
  output_root <- file.path(project_dir, "output")
}
output_root <- normalizePath(
  output_root,
  winslash = "/",
  mustWork = FALSE
)

run_definitions <- list(
  wasp = list(
    run_label = "Vespidae_plus_outgroups",
    run_dir = file.path(excon_dir, "wasp_CAFE"),
    nodes = data.table(
      node = c("20", "25"),
      biological_label = c(
        "crown social Vespidae (Polistinae + Vespinae)",
        "crown Vespinae (advanced eusocial wasps)"
      ),
      transition_class = c(
        "origin of eusociality in sampled social Vespidae",
        "advanced eusociality"
      )
    )
  )
)

required_inputs <- c(
  file.path(mapping_dir, "OG_GO_format_gene2GO.rds"),
  file.path(run_definitions$wasp$run_dir, "cafe", "base", "N0.tsv"),
  file.path(excon_dir, "annotation", "000_BUSCO_annotations_summary.tsv"),
  file.path(excon_dir, "annotation", "expansion_functional_classification_curated.tsv")
)
for (definition in run_definitions) {
  required_inputs <- c(
    required_inputs,
    file.path(
      definition$run_dir, "cafe", "significant_families",
      "combined_changes_per_node.tsv"
    ),
    file.path(definition$run_dir, "cafe", "best", "Gamma_family_results.txt"),
    file.path(
      definition$run_dir, "cafe", "large_families", "Out_cafe_large",
      "Base_family_results.txt"
    ),
    file.path(definition$run_dir, "cafe", "base", "hog_gene_counts.tsv"),
    file.path(definition$run_dir, "cafe", "base", "hog_gene_counts_large.tsv"),
    file.path(definition$run_dir, "cafe", "best", "Gamma_change.tab"),
    file.path(
      definition$run_dir, "cafe", "large_families", "Out_cafe_large",
      "Base_change.tab"
    )
  )
}
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Required inputs are missing:\n", paste(missing_inputs, collapse = "\n"))
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
gene2go <- load_og_go(mapping_dir)

read_run_events <- function(definition, run_id) {
  event_file <- file.path(
    definition$run_dir, "cafe", "significant_families",
    "combined_changes_per_node.tsv"
  )
  events <- fread(event_file, colClasses = "character")
  read_family_results <- function(path, result_source) {
    x <- fread(path, colClasses = "character", check.names = FALSE)
    setnames(x, names(x)[1:2], c("HOG", "family_p"))
    x[, `:=`(
      HOG = as.character(HOG),
      family_p = as.numeric(family_p),
      family_result_source = result_source
    )]
    x[, .(HOG, family_p, family_result_source)]
  }
  family_results <- rbindlist(list(
    read_family_results(
      file.path(definition$run_dir, "cafe", "best", "Gamma_family_results.txt"),
      "main_model"
    ),
    read_family_results(
      file.path(
        definition$run_dir, "cafe", "large_families", "Out_cafe_large",
        "Base_family_results.txt"
      ),
      "large_family_individual_run"
    )
  ))
  events[, `:=`(
    node = as.character(Node),
    HOG = as.character(HOG),
    change = as.numeric(change),
    branch_p = as.numeric(pvalue),
    direction = as.character(direction),
    result_source = fifelse(
      as.character(source) == "large_family_own_lambda",
      "large_family_individual_run",
      as.character(source)
    ),
    run = run_id,
    run_label = definition$run_label
  )]
  events <- merge(events, family_results, by = "HOG", all.x = TRUE, sort = FALSE)
  events <- events[
    node %in% definition$nodes$node &
      !is.na(family_p) & family_p < family_p_cutoff &
      !is.na(branch_p) & branch_p < branch_p_cutoff &
      direction %in% c("expansion", "contraction")
  ]
  if (events[, any(result_source != family_result_source)]) {
    stop("CAFE event and family-result source mismatch in ", run_id, ".")
  }
  events <- merge(
    events[, .(
      run, run_label, node, HOG, change, family_p, branch_p, direction,
      result_source
    )],
    definition$nodes,
    by = "node",
    all.x = TRUE,
    sort = FALSE
  )
  events[, node_order := match(node, definition$nodes$node)]
  setorder(events, node_order, direction, branch_p, HOG)
  events[, node_order := NULL]
  events[]
}

read_tested_hogs <- function(definition) {
  main <- fread(
    file.path(definition$run_dir, "cafe", "best", "Gamma_change.tab"),
    select = "FamilyID",
    colClasses = "character"
  )$FamilyID
  large <- fread(
    file.path(
      definition$run_dir, "cafe", "large_families", "Out_cafe_large",
      "Base_change.tab"
    ),
    select = "FamilyID",
    colClasses = "character"
  )$FamilyID
  sort(unique(c(main, large)))
}

read_count_table <- function(path) {
  x <- fread(path, check.names = FALSE)
  x[, Desc := NULL]
  species <- setdiff(names(x), "HOG")
  for (column in species) set(x, j = column, value = as.numeric(x[[column]]))
  attr(x, "species") <- species
  x
}

summarise_root_excluded <- function(definition, run_id, track) {
  if (track == "main_model") {
    count_path <- file.path(definition$run_dir, "cafe", "base", "hog_gene_counts.tsv")
    result_path <- file.path(definition$run_dir, "cafe", "best", "Gamma_family_results.txt")
  } else {
    count_path <- file.path(
      definition$run_dir, "cafe", "base", "hog_gene_counts_large.tsv"
    )
    result_path <- file.path(
      definition$run_dir, "cafe", "large_families", "Out_cafe_large",
      "Base_family_results.txt"
    )
  }
  counts <- read_count_table(count_path)
  species <- attr(counts, "species")
  retained_hogs <- fread(result_path, colClasses = "character")[[1L]]
  excluded <- counts[!HOG %in% retained_hogs]
  if (!nrow(excluded)) return(data.table())

  polistinae <- intersect(c(
    "Mischocyttarus_mexicanus", "Polistes_dominula", "Polistes_fuscatus",
    "Polistes_exclamans", "Polistes_canadensis"
  ), species)
  vespinae <- intersect(c(
    "Vespula_germanica", "Vespula_pensylvanica", "Vespula_vulgaris",
    "Vespa_crabro", "Vespa_mandarinia", "Vespa_velutina",
    "Dolichovespula_media", "Dolichovespula_saxonica",
    "Dolichovespula_sylvestris"
  ), species)
  outgroups <- intersect(
    c("Ancistrocerus_nigricornis", "Anoplius_nigerrimus", "Tiphia_femorata"),
    species
  )
  count_matrix <- as.matrix(excluded[, ..species])
  storage.mode(count_matrix) <- "numeric"
  row_sums <- function(columns) {
    if (length(columns)) rowSums(excluded[, ..columns]) else rep(0, nrow(excluded))
  }
  tip_count <- function(species_name) {
    if (species_name %in% species) excluded[[species_name]] else rep(NA_real_, nrow(excluded))
  }
  data.table(
    run = run_id,
    run_label = definition$run_label,
    CAFE_track = track,
    HOG = excluded$HOG,
    n_species_in_count_table = length(species),
    n_species_present = rowSums(count_matrix > 0),
    max_tip_copy_number = apply(count_matrix, 1L, max),
    total_tip_copy_number = rowSums(count_matrix),
    polistinae_total_copies = row_sums(polistinae),
    vespinae_total_copies = row_sums(vespinae),
    social_wasp_total_copies = row_sums(c(polistinae, vespinae)),
    outgroup_total_copies = row_sums(outgroups),
    Vespula_vulgaris_copies = tip_count("Vespula_vulgaris"),
    Polistes_dominula_copies = tip_count("Polistes_dominula"),
    Ancistrocerus_nigricornis_copies = tip_count("Ancistrocerus_nigricornis"),
    inference = paste(
      "Prepared for CAFE but absent from family results; inferred root-size-zero",
      "filtering. This HOG was not tested by CAFE."
    )
  )
}

event_list <- Map(
  function(definition, run_id) read_run_events(definition, run_id),
  run_definitions,
  names(run_definitions)
)
all_events <- rbindlist(event_list, use.names = TRUE, fill = TRUE)
root_excluded <- rbindlist(lapply(c("main_model", "large_family_individual_run"), function(track) {
  summarise_root_excluded(run_definitions$wasp, "wasp_CAFE", track)
}), use.names = TRUE, fill = TRUE)
tested_hogs_all <- read_tested_hogs(run_definitions$wasp)
annotation_hogs <- sort(unique(c(all_events$HOG, tested_hogs_all, root_excluded$HOG)))

# Build the annotation table from the original broad-species N0 table so
# Drosophila evidence remains available even though it was not in the wasp tree.
annotation_config <- make_wasp_annotation_config(
  project_dir = project_dir,
  excon_dir = excon_dir,
  run_dir = run_definitions$wasp$run_dir
)

annotation_required <- c(
  annotation_config$n0_file,
  annotation_config$vv_genetable_file,
  annotation_config$pd_genetable_file,
  annotation_config$dmel_genetable_file,
  annotation_config$pd_vv_annotated_ortholog_file,
  annotation_config$pd_eggnog_annotation_file,
  annotation_config$vv_eggnog_annotation_file,
  annotation_config$ancistrocerus_ortholog_file
)
missing_annotation <- annotation_required[!file.exists(annotation_required)]
if (length(missing_annotation)) {
  stop(
    "Required annotation inputs are missing:\n",
    paste(missing_annotation, collapse = "\n")
  )
}

n0 <- fread(
  annotation_config$n0_file,
  select = c(
    "HOG", "OG", "Gene Tree Parent Clade",
    "Vespula_vulgaris", "Polistes_dominula",
    "Ancistrocerus_nigricornis", "Drosophila_melanogaster"
  ),
  colClasses = "character",
  na.strings = NULL,
  check.names = FALSE
)
n0_selected <- n0[HOG %in% annotation_hogs]
missing_hogs <- setdiff(annotation_hogs, n0_selected$HOG)
if (length(missing_hogs)) {
  stop("Focal HOGs absent from the wasp CAFE N0 table: ", paste(missing_hogs, collapse = ", "))
}

direct_annotations <- read_direct_annotations(annotation_config)
dmel_annotations <- read_dmel_annotations(annotation_config)
members <- make_member_long(n0_selected, direct_annotations)
anc_ids <- members[
  species == "Ancistrocerus_nigricornis",
  unique(gene_id)
]
anc_inferred <- read_ancistrocerus_inferred_annotations(
  annotation_config,
  anc_ids,
  direct_annotations,
  dmel_annotations
)
members <- annotate_ancistrocerus_members(members, anc_inferred)
dmel_evidence <- build_dmel_hog_evidence(n0_selected, dmel_annotations)
consensus <- build_consensus_annotations(members, dmel_evidence, annotation_hogs)
member_summary <- build_member_summary(members, annotation_hogs)

hog_annotations <- merge(
  n0_selected[, .(
    HOG,
    OG,
    gene_tree_parent_clade = `Gene Tree Parent Clade`
  )],
  consensus,
  by = "HOG",
  all.x = TRUE
)
hog_annotations <- merge(
  hog_annotations,
  member_summary,
  by = "HOG",
  all.x = TRUE
)

# Retain a previously curated label when the automated re-annotation is empty.
previous_file <- file.path(
  excon_dir, "annotation", "HOG_node_changes_and_annotations.tsv"
)
if (file.exists(previous_file)) {
  previous <- unique(fread(previous_file)[, .(
    HOG,
    previous_consensus_gene_name = consensus_gene_name,
    previous_consensus_method = consensus_method,
    previous_consensus_evidence = consensus_evidence
  )], by = "HOG")
  hog_annotations <- merge(hog_annotations, previous, by = "HOG", all.x = TRUE)
  use_previous <-
    !is.na(hog_annotations$previous_consensus_gene_name) &
    nzchar(hog_annotations$previous_consensus_gene_name) &
    hog_annotations$previous_consensus_gene_name != "unannotated HOG" &
    (is.na(hog_annotations$consensus_gene_name) |
       hog_annotations$consensus_gene_name == "unannotated HOG")
  hog_annotations[use_previous, `:=`(
    consensus_gene_name = previous_consensus_gene_name,
    consensus_method = paste0(
      previous_consensus_method,
      "; retained from previous curated table"
    ),
    consensus_evidence = previous_consensus_evidence
  )]
  hog_annotations[, c(
    "previous_consensus_gene_name", "previous_consensus_method",
    "previous_consensus_evidence"
  ) := NULL]
}
hog_annotations <- add_interest_flags(hog_annotations)
te_evidence_columns <- intersect(c(
  "consensus_gene_name", "consensus_evidence",
  "Vespula_vulgaris_gene_names", "Polistes_dominula_gene_names",
  "Ancistrocerus_nigricornis_gene_names", "Drosophila_melanogaster_gene_names"
), names(hog_annotations))
te_text <- do.call(
  paste,
  c(as.list(hog_annotations[, ..te_evidence_columns]), sep = " | ")
)
te_pattern <- paste(c(
  "transposable[- ]element", "retrotranspos", "retroviral",
  "reverse transcriptase", "rna-directed dna polymerase", "transposase",
  "piggybac", "tigger transposable", "mariner transpos", "helitron",
  "gypsy[- ]like", "copia[- ]like", "pao retrotranspos", "gag-pol polyprotein"
), collapse = "|")
hog_annotations[, TE_related_HOG_flag := grepl(
  te_pattern, te_text, ignore.case = TRUE, perl = TRUE
)]
hog_annotations[, TE_filter_reason := fifelse(
  TE_related_HOG_flag,
  "Conservative downstream TE exclusion based on direct EggNOG/member annotation or consensus name",
  ""
)]
setorder(hog_annotations, HOG)

all_annotated_events <- merge(
  all_events,
  hog_annotations,
  by = "HOG",
  all.x = TRUE,
  sort = FALSE
)
all_annotated_events[, run_order := match(run, names(run_definitions))]
all_annotated_events[, node_order := mapply(
  function(run_id, node_id) match(node_id, run_definitions[[run_id]]$nodes$node),
  run,
  node
)]
setorder(all_annotated_events, run_order, node_order, direction, branch_p, HOG)
all_annotated_events[, c("run_order", "node_order") := NULL]
te_excluded_events <- all_annotated_events[TE_related_HOG_flag %in% TRUE]
annotated_events <- all_annotated_events[!TE_related_HOG_flag %in% TRUE]
setcolorder(
  annotated_events,
  c(
    "run", "run_label", "node", "biological_label", "transition_class",
    "direction", "HOG", "OG", "gene_tree_parent_clade",
    "consensus_gene_name", "consensus_method", "consensus_evidence",
    "change", "family_p", "branch_p", "result_source",
    "TE_related_HOG_flag", "TE_filter_reason", "chemosensory_family_flag",
    setdiff(
      names(annotated_events),
      c(
        "run", "run_label", "node", "biological_label",
        "transition_class", "direction", "HOG", "OG",
        "gene_tree_parent_clade", "consensus_gene_name",
        "consensus_method", "consensus_evidence", "change", "family_p",
        "branch_p", "result_source", "TE_related_HOG_flag",
        "TE_filter_reason", "chemosensory_family_flag"
      )
    )
  )
)
setcolorder(te_excluded_events, names(annotated_events))

node_summary <- annotated_events[, .(
  n_significant_expansions = sum(direction == "expansion"),
  n_significant_contractions = sum(direction == "contraction"),
  n_expansions_main_model = sum(
    direction == "expansion" & result_source == "main_model"
  ),
  n_expansions_large_family = sum(
    direction == "expansion" & result_source != "main_model"
  ),
  n_contractions_main_model = sum(
    direction == "contraction" & result_source == "main_model"
  ),
  n_contractions_large_family = sum(
    direction == "contraction" & result_source != "main_model"
  ),
  n_unique_changed_HOGs = uniqueN(HOG),
  n_named_HOGs = uniqueN(HOG[consensus_gene_name != "unannotated HOG"]),
  n_chemosensory_HOGs = uniqueN(HOG[chemosensory_family_flag %in% TRUE])
), by = .(
  run, run_label, node, biological_label, transition_class
)]
te_excluded_summary <- te_excluded_events[, .(
  n_TE_expansions_excluded = sum(direction == "expansion"),
  n_TE_contractions_excluded = sum(direction == "contraction"),
  excluded_TE_HOGs = paste(sort(unique(HOG)), collapse = ";")
), by = .(run, node)]
node_summary <- merge(
  node_summary,
  te_excluded_summary,
  by = c("run", "node"),
  all.x = TRUE,
  sort = FALSE
)
node_summary[is.na(n_TE_expansions_excluded), n_TE_expansions_excluded := 0L]
node_summary[is.na(n_TE_contractions_excluded), n_TE_contractions_excluded := 0L]
node_summary[is.na(excluded_TE_HOGs), excluded_TE_HOGs := ""]
node_summary[, `:=`(
  family_wide_screen = paste0("p < ", family_p_cutoff),
  focal_branch_threshold = paste0("p < ", branch_p_cutoff),
  downstream_TE_policy = "excluded from totals and biological interpretation"
)]
node_summary[, run_order := match(run, names(run_definitions))]
node_summary[, node_order := mapply(
  function(run_id, node_id) match(node_id, run_definitions[[run_id]]$nodes$node),
  run,
  node
)]
setorder(node_summary, run_order, node_order)
node_summary[, c("run_order", "node_order") := NULL]

root_excluded <- merge(root_excluded, hog_annotations, by = "HOG", all.x = TRUE)
root_excluded[, present_in_social_wasps := social_wasp_total_copies > 0]
root_excluded[, biological_interest_class := fcase(
  chemosensory_family_flag %in% TRUE, "chemoperception",
  grepl(
    "cytochrome p450|desaturase|elongase|fatty acid|lipid",
    consensus_gene_name, ignore.case = TRUE
  ), "lipid metabolism / possible chemical communication",
  grepl(
    "venom|phospholipase|toxin|defen[cs]e",
    consensus_gene_name, ignore.case = TRUE
  ), "venom / defence",
  grepl(
    "insulin|juvenile hormone|ecdys|vitellogenin",
    consensus_gene_name, ignore.case = TRUE
  ), "endocrine / reproductive physiology",
  grepl(
    "protease|peptidase|maltase|glucosidase|amylase|lipase|digestion",
    consensus_gene_name, ignore.case = TRUE
  ), "digestion / nutrition",
  default = ""
)]
root_excluded[, interpretation_note := paste(
  "Not tested by CAFE because the family was removed after an inferred root",
  "size of zero; tip-copy counts are descriptive and no branch significance",
  "can be assigned from the wasp CAFE analysis."
)]

track_input_sizes <- rbindlist(lapply(c(
  "main_model", "large_family_individual_run"
), function(track) {
  count_file <- if (track == "main_model") {
    "hog_gene_counts.tsv"
  } else {
    "hog_gene_counts_large.tsv"
  }
  data.table(
    CAFE_track = track,
    n_prepared_for_CAFE = nrow(read_count_table(file.path(
      run_definitions$wasp$run_dir, "cafe", "base", count_file
    )))
  )
}))
root_zero_summary <- root_excluded[, .(
  n_removed_empty_at_root = .N,
  n_present_in_social_wasps = sum(present_in_social_wasps),
  n_annotated = sum(consensus_gene_name != "unannotated HOG", na.rm = TRUE),
  n_chemosensory = sum(chemosensory_family_flag %in% TRUE, na.rm = TRUE),
  n_other_biologically_relevant = sum(
    nzchar(biological_interest_class) &
      biological_interest_class != "chemoperception",
    na.rm = TRUE
  )
), by = .(run, run_label, CAFE_track)]
root_zero_summary <- merge(
  root_zero_summary, track_input_sizes,
  by = "CAFE_track", all.x = TRUE, sort = FALSE
)
root_zero_summary[, percent_of_prepared_families_removed :=
  100 * n_removed_empty_at_root / n_prepared_for_CAFE]
root_zero_summary[, interpretation := paste(
  "These HOGs were not tested by CAFE; no significant branch expansion or",
  "contraction can be inferred for them in the wasp CAFE analysis."
)]
root_zero_biologically_relevant <- root_excluded[
  present_in_social_wasps & nzchar(biological_interest_class)
]
setorder(root_zero_biologically_relevant, biological_interest_class, HOG)

figure_curation <- fread(file.path(
  excon_dir, "annotation",
  "expansion_functional_classification_curated.tsv"
))
figure_events <- merge(
  annotated_events[direction == "expansion", .(node, HOG)],
  figure_curation[, .(
    HOG, primary_category, vizueta_like_group, confidence
  )],
  by = "HOG",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(figure_events$vizueta_like_group)) {
  stop("At least one expanded focal HOG lacks a curated figure category.")
}
figure_events[, figure_category := fcase(
  vizueta_like_group == "Chemoperception", "Chemoperception",
  vizueta_like_group %chin% c(
    "CHC/lipid metabolism", "CHC/lipid metabolism candidate"
  ), "Lipid / CHC",
  vizueta_like_group %chin% c(
    "Digestion", "Digestion candidate", "Digestion/nutrient metabolism"
  ), "Digestion / nutrition",
  vizueta_like_group %chin% c(
    "Endocrine signalling", "Endocrine signalling candidate"
  ), "Endocrine",
  vizueta_like_group == "Endocrine/behavioural signalling", "Behavioural signal",
  vizueta_like_group %chin% c(
    "Venom", "Venom candidate", "Defence/detoxification candidate"
  ), "Defence / venom",
  default = "Other / unresolved"
)]
figure_category_levels <- c(
  "Lipid / CHC", "Chemoperception", "Digestion / nutrition",
  "Endocrine", "Behavioural signal", "Defence / venom",
  "Other / unresolved"
)
focal_transition_labels <- data.table(
  node = c("20", "25"),
  focal_transition = c(
    "Stem Polistinae+Vespinae; associated with the inferred origin of eusociality in this sampled wasp clade",
    "Stem Vespinae; associated with the origin of advanced eusociality and morphologically distinct castes"
  )
)
figure_category_summary <- figure_events[, .(
  n_HOGs = uniqueN(HOG),
  HOGs = paste(sort(unique(HOG)), collapse = ";")
), by = .(node, figure_category)]
figure_category_summary <- merge(
  figure_category_summary,
  focal_transition_labels,
  by = "node",
  all.x = TRUE,
  sort = FALSE
)
figure_category_summary[, category_order := match(
  figure_category, figure_category_levels
)]
figure_category_summary[, percent_of_expanded_HOGs :=
  100 * n_HOGs / sum(n_HOGs), by = node]
setcolorder(
  figure_category_summary,
  c(
    "node", "focal_transition", "figure_category", "n_HOGs",
    "HOGs", "category_order", "percent_of_expanded_HOGs"
  )
)
figure_category_summary[, node_order := match(node, c("20", "25"))]
setorder(figure_category_summary, node_order, category_order)
figure_category_summary[, node_order := NULL]
fwrite(
  figure_category_summary,
  file.path(
    output_root,
    "node20_node25_expansion_figure_category_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

fwrite(
  hog_annotations[HOG %in% annotated_events$HOG],
  file.path(output_root, "focal_significant_HOG_annotations_unique.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  annotated_events,
  file.path(output_root, "focal_node_significant_changes_annotated.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  te_excluded_events,
  file.path(output_root, "focal_TE_HOGs_excluded_from_interpretation.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  hog_annotations[
    HOG %in% tested_hogs_all & TE_related_HOG_flag %in% TRUE
  ],
  file.path(output_root, "tested_TE_HOGs_excluded_from_GO_and_interpretation.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  node_summary,
  file.path(output_root, "focal_node_significant_change_summary.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  root_excluded,
  file.path(output_root, "root_zero_filtered_HOGs_annotated.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  root_zero_summary,
  file.path(output_root, "root_zero_filter_summary.tsv"),
  sep = "\t",
  quote = FALSE
)
fwrite(
  root_zero_biologically_relevant,
  file.path(output_root, "root_zero_biologically_relevant_HOGs.tsv"),
  sep = "\t",
  quote = FALSE
)

member_evidence <- rbindlist(list(
  members[HOG %in% annotated_events$HOG, .(
    HOG, species, gene_id, gene_symbol, gene_name,
    annotation_source, annotation_quality
  )],
  dmel_evidence[HOG %in% annotated_events$HOG, .(
    HOG, species, gene_id, gene_symbol, gene_name,
    annotation_source = evidence_source,
    annotation_quality = "FlyBase-derived outgroup evidence"
  )]
), use.names = TRUE, fill = TRUE)
setorder(member_evidence, HOG, species, gene_id)
fwrite(
  member_evidence,
  file.path(output_root, "focal_significant_HOG_member_annotation_evidence.tsv"),
  sep = "\t",
  quote = FALSE
)

go_manifest <- list()
go_filtered_all <- list()
go_contributors_all <- list()

for (run_id in names(run_definitions)) {
  definition <- run_definitions[[run_id]]
  run_output <- file.path(output_root, run_id)
  dir.create(run_output, recursive = TRUE, showWarnings = FALSE)
  for (subdir in c("foregrounds", "results")) {
    generated_dir <- file.path(run_output, subdir)
    if (dir.exists(generated_dir)) {
      unlink(list.files(generated_dir, full.names = TRUE), recursive = TRUE)
    }
    dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir.create(file.path(run_output, "backgrounds"), recursive = TRUE, showWarnings = FALSE)

  tested_hogs_before_TE_filter <- read_tested_hogs(definition)
  tested_te_hogs <- hog_annotations[
    HOG %in% tested_hogs_before_TE_filter & TE_related_HOG_flag %in% TRUE,
    HOG
  ]
  tested_hogs <- setdiff(tested_hogs_before_TE_filter, tested_te_hogs)
  background <- sort(intersect(tested_hogs, names(gene2go)))
  writeLines(
    tested_hogs,
    file.path(run_output, "backgrounds", "tested_HOGs_all.txt")
  )
  writeLines(
    background,
    file.path(run_output, "backgrounds", "tested_GO_annotated_HOGs.txt")
  )
  writeLines(
    sort(tested_te_hogs),
    file.path(run_output, "backgrounds", "tested_TE_HOGs_excluded.txt")
  )

  for (i in seq_len(nrow(definition$nodes))) {
    node_id <- definition$nodes$node[[i]]
    node_label <- definition$nodes$biological_label[[i]]
    selected <- annotated_events[
      run == run_id & node == node_id & direction == "expansion",
      .(
        HOG, change, family_p, branch_p, result_source, consensus_gene_name,
        go_annotated = HOG %in% background
      )
    ]
    foreground <- intersect(selected$HOG, background)
    target_name <- paste0("Node_", node_id)

    fwrite(
      selected,
      file.path(
        run_output, "foregrounds",
        paste0(target_name, ".significant.expanded.metadata.tsv")
      ),
      sep = "\t",
      quote = FALSE
    )
    writeLines(
      selected$HOG,
      file.path(
        run_output, "foregrounds",
        paste0(target_name, ".significant.expanded.txt")
      )
    )

    result_dir <- file.path(
      run_output, "results", target_name, "significant_expanded",
      paste(algorithm, statistic, sep = "_")
    )
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

    combined <- lapply(ontologies, function(ontology) {
      message("Running ", run_id, " / ", target_name, " / ", ontology)
      run_topgo_test(
        foreground = foreground,
        universe = background,
        gene2go = gene2go,
        ontology = ontology,
        algorithm = algorithm,
        statistic = statistic,
        node_size = 1L,
        top_nodes = top_nodes_per_ontology
      )
    })
    ranked <- rbindlist(combined, use.names = TRUE, fill = TRUE)
    filtered <- ranked[!is.na(p_raw) & p_raw < go_p_cutoff][
      order(p_raw, -FoldEnrichment)
    ]

    fwrite(ranked, file.path(result_dir, "topGO_ranked.tsv"), sep = "\t")
    fwrite(filtered, file.path(result_dir, "topGO_filtered.tsv"), sep = "\t")
    ranked_members <- write_go_member_table(
      ranked,
      selected,
      file.path(result_dir, "topGO_ranked_contributing_HOGs.tsv")
    )
    filtered_members <- write_go_member_table(
      filtered,
      selected,
      file.path(result_dir, "topGO_filtered_contributing_HOGs.tsv")
    )

    metadata <- data.table(
      run = run_id,
      node = node_id,
      biological_label = node_label,
      n_significant_expanded_HOGs = nrow(selected),
      n_GO_annotated_foreground_HOGs = length(foreground),
      n_tested_HOGs_before_TE_filter = length(tested_hogs_before_TE_filter),
      n_TE_HOGs_excluded_from_background = length(tested_te_hogs),
      n_tested_HOGs = length(tested_hogs),
      n_GO_annotated_background_HOGs = length(background),
      n_nominally_significant_GO_terms = nrow(filtered),
      ontologies = paste(ontologies, collapse = ","),
      algorithm = algorithm,
      statistic = statistic,
      family_p_cutoff = family_p_cutoff,
      branch_p_cutoff = branch_p_cutoff,
      p_value_rule = paste0(
        "CAFE family p < ", family_p_cutoff,
        " and branch p < ", branch_p_cutoff
      ),
      downstream_TE_policy = paste(
        "Conservatively annotated TE HOGs excluded from foreground,",
        "background, totals and biological interpretation"
      ),
      GO_p_cutoff = go_p_cutoff,
      GO_p_adjustment = "none",
      top_nodes_per_ontology = "Inf"
    )
    fwrite(metadata, file.path(result_dir, "analysis_metadata.tsv"), sep = "\t")
    go_manifest[[paste(run_id, node_id, sep = "_")]] <- metadata

    if (nrow(filtered)) {
      filtered[, `:=`(
        run = run_id,
        node = node_id,
        biological_label = node_label
      )]
      go_filtered_all[[paste(run_id, node_id, sep = "_")]] <- filtered
    }
    if (nrow(filtered_members)) {
      filtered_members <- merge(
        filtered_members,
        hog_annotations[, .(HOG, consensus_gene_name)],
        by = "HOG",
        all.x = TRUE
      )
      filtered_members[, `:=`(
        run = run_id,
        node = node_id,
        biological_label = node_label
      )]
      go_contributors_all[[paste(run_id, node_id, sep = "_")]] <- filtered_members
    }
  }
}

go_manifest <- rbindlist(go_manifest, use.names = TRUE, fill = TRUE)
go_filtered <- rbindlist(go_filtered_all, use.names = TRUE, fill = TRUE)
go_contributors <- rbindlist(go_contributors_all, use.names = TRUE, fill = TRUE)

fwrite(go_manifest, file.path(output_root, "GO_analysis_manifest.tsv"), sep = "\t")
fwrite(
  go_filtered,
  file.path(output_root, "topGO_nominal_p_lt_0.05_all_focal_nodes.tsv"),
  sep = "\t"
)
fwrite(
  go_filtered[!is.na(p_raw) & p_raw < 0.01 & Members_n >= 2L],
  file.path(output_root, "topGO_p_lt_0.01_at_least_2_HOGs.tsv"),
  sep = "\t"
)
fwrite(
  go_contributors,
  file.path(output_root, "topGO_nominal_significant_contributing_HOGs.tsv"),
  sep = "\t"
)

# Add BUSCO context for every tip descending from the focal clades. The full
# table is retained because branch calls can be influenced by any included tip.
busco <- fread(
  file.path(excon_dir, "annotation", "000_BUSCO_annotations_summary.tsv")
)
read_tree_species <- function(definition) {
  tree <- paste(readLines(file.path(
    definition$run_dir, "cafe", "base", "SpeciesTree_rooted_ultra.txt"
  )), collapse = "")
  hits <- regmatches(
    tree,
    gregexpr("(?<=[(,])[^():,]+(?=:)", tree, perl = TRUE)
  )[[1L]]
  sort(unique(trimws(hits)))
}
wasp_species <- read_tree_species(run_definitions$wasp)
busco_public_columns <- intersect(c(
  "species", "source", "busco_version", "lineage", "mode",
  "busco_C_pct", "busco_S_pct", "busco_D_pct", "busco_F_pct",
  "busco_M_pct", "busco_n", "n_complete", "n_single", "n_duplicated",
  "n_fragmented", "n_missing"
), names(busco))
busco_context <- copy(busco[species %in% wasp_species, ..busco_public_columns])
busco_context[, analysis := "Vespidae_plus_outgroups"]
setcolorder(busco_context, c("analysis", setdiff(names(busco_context), "analysis")))
setorder(busco_context, analysis, species)
fwrite(
  busco_context,
  file.path(output_root, "BUSCO_context_included_species.tsv"),
  sep = "\t"
)

model_summary <- rbindlist(lapply(names(run_definitions), function(run_id) {
  definition <- run_definitions[[run_id]]
  model_selection <- fread(file.path(
    definition$run_dir, "cafe", "model_comparison", "model_selection.tsv"
  ))
  selected <- model_selection[Selected == "BEST"]
  results_lines <- readLines(file.path(
    definition$run_dir, "cafe", "best", "Gamma_results.txt"
  ))
  lambda_line <- grep("^Lambda:", results_lines, value = TRUE)
  data.table(
    run = run_id,
    run_label = definition$run_label,
    best_k = selected$k,
    n_standard_families = selected$n_families,
    neg_lnL = selected$neg_lnL,
    AIC = selected$AIC,
    lambda = as.numeric(sub("^Lambda:\\s*", "", lambda_line)),
    n_large_families = length(read_tested_hogs(definition)) - selected$n_families
  )
}))
fwrite(model_summary, file.path(output_root, "CAFE_model_summary.tsv"), sep = "\t")

# Compact manuscript-facing tables. Larger audit tables remain in output/ and
# can be regenerated from the versioned inputs and key-results RDS.
supplemental_table_dir <- file.path(output_root, "supplemental_tables")
dir.create(supplemental_table_dir, recursive = TRUE, showWarnings = FALSE)

gamma_lines <- readLines(file.path(
  run_definitions$wasp$run_dir, "cafe", "best", "Gamma_results.txt"
))
extract_gamma_value <- function(prefix) {
  line <- grep(paste0("^", prefix, ":"), gamma_lines, value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(trimws(sub(paste0("^", prefix, ":"), "", line[[1L]])))
}
standard_prepared <- track_input_sizes[
  CAFE_track == "main_model", n_prepared_for_CAFE
]
large_prepared <- track_input_sizes[
  CAFE_track == "large_family_individual_run", n_prepared_for_CAFE
]
standard_root_zero <- root_zero_summary[
  CAFE_track == "main_model", n_removed_empty_at_root
]
large_root_zero <- root_zero_summary[
  CAFE_track == "large_family_individual_run", n_removed_empty_at_root
]
node20_counts <- node_summary[node == "20"]
node25_counts <- node_summary[node == "25"]

table_s10 <- data.table(
  Quantity = c(
    "Included species", "Prepared standard-track HOGs",
    "Prepared high-differential HOGs", "Root-zero standard HOGs",
    "Root-zero high-differential HOGs", "Standard HOGs with results",
    "High-differential HOGs with results", "Gamma categories",
    "Root-family-size distribution", "Mean lambda", "Gamma alpha",
    "Symmetric error", "Focal family screen", "Focal branch threshold",
    "Node 20 non-TE changes", "Node 25 non-TE changes"
  ),
  Value = c(
    length(wasp_species), standard_prepared, large_prepared,
    standard_root_zero, large_root_zero,
    model_summary$n_standard_families, model_summary$n_large_families,
    model_summary$best_k, "Poisson", signif(model_summary$lambda, 6),
    signif(extract_gamma_value("Alpha"), 6),
    signif(extract_gamma_value("Epsilon"), 6),
    paste0("P < ", family_p_cutoff), paste0("Viterbi probability < ", branch_p_cutoff),
    paste0("+", node20_counts$n_significant_expansions, " / -", node20_counts$n_significant_contractions),
    paste0("+", node25_counts$n_significant_expansions, " / -", node25_counts$n_significant_contractions)
  ),
  Interpretation = c(
    "15 Vespidae plus Tiphia femorata and Anoplius nigerrimus outgroups",
    "Before CAFE root-size filtering", "Maximum among-species differential >20",
    "No family-wide or branch-specific result assigned",
    "No family-wide or branch-specific result assigned",
    "Joint Gamma-rate model", "Fitted as individual high-differential families",
    "Selected by AIC", "Favoured over uniform", "Per gene per million years",
    "Gamma-rate shape parameter", "Gene-count error model",
    "Nominal family-wide screen; no FDR correction",
    "Focal-branch Viterbi-probability threshold; no FDR correction",
    "Displayed and interpreted after TE exclusion",
    "Displayed and interpreted after TE exclusion"
  )
)
fwrite(table_s10, file.path(supplemental_table_dir, "TableS10.tsv"), sep = "\t")

event_figure_categories <- unique(figure_events[, .(node, HOG, figure_category)])
table_s11 <- merge(
  annotated_events,
  event_figure_categories,
  by = c("node", "HOG"), all.x = TRUE, sort = FALSE
)
table_s11[is.na(figure_category), figure_category := "Not classified"]
table_s11[, Track := fifelse(
  result_source == "main_model", "Standard", "Large-family"
)]
table_s11 <- table_s11[, .(
  Node = node,
  Direction = direction,
  HOG,
  Change = change,
  `Family P` = family_p,
  `Viterbi branch probability` = branch_p,
  Track,
  `Consensus HOG name` = consensus_gene_name,
  `Fig. 3 class` = figure_category,
  `Vespula vulgaris members` = Vespula_vulgaris_gene_names,
  `Polistes dominula members` = Polistes_dominula_gene_names,
  `Ancistrocerus nigricornis members` = fifelse(
    !is.na(Ancistrocerus_nigricornis_gene_names) &
      nzchar(Ancistrocerus_nigricornis_gene_names),
    Ancistrocerus_nigricornis_gene_names,
    Ancistrocerus_nigricornis_gene_ids
  )
)]
setorder(table_s11, Node, Direction, HOG)
fwrite(table_s11, file.path(supplemental_table_dir, "TableS11.tsv"), sep = "\t")

table_s12 <- go_filtered[, .(
  Node = node,
  Ontology = ontology,
  `GO ID` = GO.ID,
  `GO term` = Term,
  `Total HOGs (n)` = Annotated,
  `Foreground HOGs (n)` = Significant,
  Expected,
  `Fold enrichment` = FoldEnrichment,
  `weight01 P` = p_raw,
  `Contributing HOGs` = Members
)]
setorder(table_s12, Node, Ontology, `weight01 P`, `GO ID`)
fwrite(table_s12, file.path(supplemental_table_dir, "TableS12.tsv"), sep = "\t")

table_s13 <- root_zero_biologically_relevant[, .(
  HOG,
  `Interest class` = biological_interest_class,
  `Consensus HOG name` = consensus_gene_name,
  `Vespula vulgaris copies` = Vespula_vulgaris_copies,
  `Polistes dominula copies` = Polistes_dominula_copies,
  `Ancistrocerus nigricornis copies` = Ancistrocerus_nigricornis_copies,
  `Copies across all solitary/outgroup taxa` = outgroup_total_copies,
  `Total member copies across all social wasps` = social_wasp_total_copies
)]
setorder(table_s13, `Interest class`, HOG)
fwrite(table_s13, file.path(supplemental_table_dir, "TableS13.tsv"), sep = "\t")

annotation_qc <- data.table(
  item = c(
    "focal_event_rows", "unique_focal_HOGs", "named_unique_HOGs",
    "unannotated_unique_HOGs", "chemosensory_unique_HOGs",
    "focal_TE_HOG_rows_excluded", "tested_TE_HOGs_excluded_from_GO_universe",
    "root_zero_HOGs_not_tested_by_CAFE"
  ),
  value = c(
    nrow(annotated_events),
    uniqueN(annotated_events$HOG),
    uniqueN(annotated_events[consensus_gene_name != "unannotated HOG", HOG]),
    uniqueN(annotated_events[consensus_gene_name == "unannotated HOG", HOG]),
    uniqueN(annotated_events[chemosensory_family_flag %in% TRUE, HOG]),
    nrow(te_excluded_events),
    uniqueN(hog_annotations[
      HOG %in% tested_hogs_all & TE_related_HOG_flag %in% TRUE, HOG
    ]),
    nrow(root_excluded)
  )
)
fwrite(annotation_qc, file.path(output_root, "annotation_QC.tsv"), sep = "\t")

analysis_metadata <- data.table(
  item = c(
    "primary_analysis", "family_wide_screen", "focal_branch_threshold",
    "threshold_logic", "multiple_testing", "TE_policy",
    "TE_scope_caution", "root_zero_policy", "GO_method", "GO_ontologies",
    "text_truncation_policy"
  ),
  value = c(
    "EXCON/CAFE analysis of Vespidae plus Ancistrocerus, Anoplius and Tiphia outgroups",
    paste0("CAFE family-wide p < ", family_p_cutoff),
    paste0("Focal-branch Viterbi probability < ", branch_p_cutoff),
    "Both strict inequalities must be satisfied for a focal event",
    "No FDR correction applied to CAFE or topGO p-values",
    paste(
      "Conservatively annotated TE HOGs are excluded downstream from biological",
      "totals, functional classification, GO foregrounds and GO backgrounds"
    ),
    paste(
      "TEs were retained in the fitted CAFE model; this is a downstream",
      "interpretive exclusion, not a refitted TE-free CAFE analysis"
    ),
    paste(
      "Families absent from CAFE family results after inferred root-size-zero",
      "filtering are reported separately and cannot receive branch significance"
    ),
    paste0("topGO ", algorithm, " / ", statistic, "; nominal p < ", go_p_cutoff),
    paste(ontologies, collapse = ","),
    paste(
      "GO names are recovered in full from GO.db by GO identifier; annotation",
      "text is exported in full and is never shortened with an ellipsis"
    )
  )
)
fwrite(
  analysis_metadata,
  file.path(output_root, "CAFE_analysis_metadata.tsv"),
  sep = "\t",
  quote = FALSE
)

tree_powerpoint <- file.path(
  output_root, "figures", "Fig3.pptx"
)
dir.create(dirname(tree_powerpoint), recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("export", quietly = TRUE)) {
  stop("Package 'export' is required to create the tree PowerPoint.")
}
if (!nzchar(Sys.getenv("GDTOOLS_CACHE_DIR"))) {
  gdtools_cache <- file.path(tempdir(), "gdtools-font-cache")
  dir.create(gdtools_cache, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(GDTOOLS_CACHE_DIR = gdtools_cache)
}
export::graph2ppt(
  fun = function() {
    plot_wasp_tree(
      tree_file = file.path(
        run_definitions$wasp$run_dir,
        "cafe", "base", "SpeciesTree_rooted_ultra.txt"
      ),
      node_summary = node_summary,
      figure_category_summary = figure_category_summary
    )
  },
  file = tree_powerpoint,
  width = 16.5,
  height = 11.7,
  margins = c(top = 0, right = 0, bottom = 0, left = 0),
  vector.graphic = TRUE
)

CAFE_key_results <- list(
  node_summary = node_summary,
  annotated_events = annotated_events,
  focal_HOG_annotations = hog_annotations[HOG %in% unique(annotated_events$HOG)],
  root_zero_filtered_HOGs = root_excluded,
  root_zero_biologically_relevant_HOGs = root_zero_biologically_relevant,
  GO_manifest = go_manifest,
  GO_terms_nominal_p_lt_0.05 = go_filtered,
  GO_contributing_HOGs = go_contributors,
  model_summary = model_summary,
  BUSCO_context = busco_context,
  TableS10 = table_s10,
  TableS11 = table_s11,
  TableS12 = table_s12,
  TableS13 = table_s13,
  Fig3 = file.path("output", "figures", "Fig3.pptx")
)
saveRDS(
  CAFE_key_results,
  file.path(output_root, "C_EXCON_CAFE_key_results.rds"),
  compress = "xz"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "CAFE_sessionInfo.txt")
)
message("Wasp focal-node analysis complete: ", output_root)
invisible(list(
  node_summary = node_summary,
  annotated_events = annotated_events,
  GO_manifest = go_manifest,
  GO_terms_p_lt_0.05 = go_filtered,
  tree_powerpoint = tree_powerpoint,
  output_root = output_root
))
}

if (isTRUE(getOption("evodevo.cafe.auto_run", TRUE))) {
  wasp_analysis <- run_wasp_cafe_analysis(project_dir = REPOSITORY_ROOT)
}
