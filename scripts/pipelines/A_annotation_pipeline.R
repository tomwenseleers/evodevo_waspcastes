# Final functional annotation for Polistes dominula, Vespula vulgaris and
# their shared N13 hierarchical orthogroups. It incorporates:
#   - Galaxy InterProScan GO and pathway annotations;
#   - Galaxy EggNOG-mapper GO and KEGG annotations;
#   - EXCON EggNOG-mapper GO annotations;
#   - GO transfer from shared Drosophila melanogaster orthologues using
#     release-matched, experimentally supported FlyBase GO evidence.
#
# The final HOG annotation combines direct GO terms conserved across the two
# wasps with experimental Drosophila-transferred GO terms. Alternative direct,
# Drosophila-only and union sensitivity tables are not exported. KEGG and other
# EggNOG functional annotations are retained as before.

suppressPackageStartupMessages({
  library(dplyr)
  library(KEGGREST)
  library(ontologyIndex)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  stop("PROJECT_ROOT must be set by A_annotation_FINAL.R before sourcing this pipeline.")
}
base_dir <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
annotation_dir <- file.path(base_dir, "input_annotation")
raw_dir <- file.path(annotation_dir, "raw")
galaxy_dir <- file.path(raw_dir, "GalaxyOutput")
excon_dir <- file.path(raw_dir, "excon_go_files")
flybase_dir <- file.path(raw_dir, "FlyBase")
output_dir <- file.path(base_dir, "output_annotation")
legacy_output_dir <- output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list(
  pd_interpro = file.path(galaxy_dir, "Polistes_dominula_InterProScan.tsv.gz"),
  vv_interpro = file.path(galaxy_dir, "Vespula_vulgaris_InterProScan.tsv.gz"),
  pd_eggnog = file.path(galaxy_dir, "Polistes_dominula_EggNOG.tsv"),
  vv_eggnog = file.path(galaxy_dir, "Vespula_vulgaris_EggNOG.tsv"),
  pd_excon_go = file.path(excon_dir, "Polistes_dominula.go.txt"),
  vv_excon_go = file.path(excon_dir, "Vespula_vulgaris.go.txt"),
  dmel_flybase_go = file.path(flybase_dir, "gene_association.fb.gz"),
  hogs = file.path(annotation_dir, "reference", "pd_vv_N13_HOGs_annotated.rds"),
  vv_lookup = file.path(annotation_dir, "reference", "vv_refseq_to_ensemblrapid_full.csv"),
  obo = file.path(raw_dir, "GO", "go-basic.obo")
)

missing_files <- names(files)[!file.exists(unlist(files))]
if (length(missing_files)) {
  stop("Missing required inputs: ", paste(missing_files, collapse = ", "))
}

message("Reading GO ontology")
go <- get_ontology(files$obo, extract_tags = "everything")
ROOT_GO <- c("GO:0005575", "GO:0003674", "GO:0008150")
primary_ids <- go$id
obsolete_vec <- setNames(as.logical(go$obsolete), primary_ids)

as_char <- function(x) {
  if (is.list(x)) unlist(x, use.names = FALSE) else as.character(x)
}

tidy_definition <- function(x) {
  x <- as_char(x)
  x <- sub("\\s*\\[.*\\]$", "", x)
  gsub('^"|"$', "", x)
}

go_lookup <- tibble(
  go_terms = go$id,
  aspect = as_char(go$namespace),
  term_name = as_char(go$name),
  definition = if ("def" %in% names(go)) tidy_definition(go$def) else NA_character_
)

alt_to_primary <- {
  alt_lengths <- lengths(go$alt_id)
  if (sum(alt_lengths) > 0L) {
    setNames(rep(primary_ids, alt_lengths), as_char(go$alt_id))
  } else {
    character()
  }
}

replacement_target <- vapply(seq_along(primary_ids), function(i) {
  replacement <- as_char(go$replaced_by[[i]])
  replacement <- replacement[replacement %in% primary_ids]
  replacement <- replacement[!obsolete_vec[replacement]]
  if (length(replacement)) replacement[[1]] else NA_character_
}, FUN.VALUE = character(1))
replacement_target <- setNames(replacement_target, primary_ids)

canonicalise_go_ids <- function(ids, max_iter = 6L) {
  ids <- unique(as.character(ids))
  ids <- ids[str_detect(ids, "^GO:\\d{7}$")]
  if (!length(ids)) return(character())

  is_alt <- ids %in% names(alt_to_primary)
  if (any(is_alt)) ids[is_alt] <- alt_to_primary[ids[is_alt]]

  for (iteration in seq_len(max_iter)) {
    obsolete <- obsolete_vec[ids]
    obsolete[is.na(obsolete)] <- FALSE
    if (!any(obsolete)) break
    replacement <- replacement_target[ids]
    replaceable <- obsolete & !is.na(replacement)
    if (!any(replaceable)) break
    ids[replaceable] <- replacement[replaceable]
  }

  obsolete <- obsolete_vec[ids]
  obsolete[is.na(obsolete)] <- FALSE
  sort(unique(setdiff(ids[!obsolete], ROOT_GO)))
}

parse_go <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(character())
  tokens <- unlist(str_split(as.character(x), "[\\s,;|]+"), use.names = FALSE)
  canonicalise_go_ids(tokens)
}

most_specific <- function(terms) {
  terms <- intersect(unique(terms), go$id)
  if (!length(terms)) return(character())
  # EggNOG outputs ancestor-rich GO sets. A term is redundant whenever it is
  # the direct parent of another retained term. Removing all such parents in
  # one graph pass is equivalent to scanning every term's children, but linear
  # in the number of represented ontology edges rather than near-quadratic.
  represented_parents <- intersect(
    unique(unlist(go$parents[terms], use.names = FALSE)),
    terms
  )
  setdiff(terms, represented_parents)
}

go_closure_no_roots <- function(terms) {
  terms <- intersect(unique(terms), go$id)
  if (!length(terms)) return(character())
  setdiff(
    unique(c(terms, unlist(go$ancestors[terms], use.names = FALSE))),
    ROOT_GO
  )
}

empty_character <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) character() else unique(x[!is.na(x)])
}

interpro_columns <- c(
  "protein_ID", "sequence_MD5_digest", "sequence_length", "analysis",
  "signature_accession", "signature_description", "start_loc", "stop_loc",
  "score", "status", "date", "interproscan_annotations",
  "interproscan_description", "GO_annotations", "pathway_annotations"
)

read_interpro <- function(path, species, skip = 0L) {
  # Rows with only 13 fields legitimately lack the optional GO/pathway fields.
  raw <- suppressWarnings(read_tsv(
      path,
      col_names = interpro_columns,
      skip = skip,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE,
      progress = FALSE
    )) %>%
    transmute(
      species = .env$species,
      gene_id = str_extract(protein_ID, "LOC\\d+"),
      GO_annotations
    ) %>%
    filter(!is.na(gene_id), gene_id != "")

  go_long <- raw %>%
    select(species, gene_id, GO_annotations) %>%
    filter(!is.na(GO_annotations), GO_annotations != "", GO_annotations != "-") %>%
    distinct() %>%
    mutate(go_terms = map(GO_annotations, parse_go)) %>%
    select(-GO_annotations) %>%
    unnest_longer(go_terms) %>%
    filter(!is.na(go_terms), go_terms != "") %>%
    distinct(species, gene_id, go_terms) %>%
    mutate(source = "Galaxy_InterProScan")

  list(go = go_long)
}

read_eggnog <- function(path, species) {
  read_tsv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(
      species = species,
      gene_id = str_extract(.data[["#query"]], "LOC\\d+")
    ) %>%
    filter(!is.na(gene_id), gene_id != "")
}

eggnog_go_long <- function(eggnog) {
  eggnog %>%
    transmute(species, gene_id, go_terms = map(GOs, parse_go)) %>%
    unnest_longer(go_terms) %>%
    filter(!is.na(go_terms), go_terms != "") %>%
    distinct(species, gene_id, go_terms) %>%
    mutate(source = "Galaxy_EggNOG_mapper_2.1.8")
}

read_excon_go <- function(path, species, vv_lookup = NULL) {
  raw <- read_tsv(
    path,
    col_names = c("source_gene_id", "go_terms"),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )

  # EXCON stores one gene-GO pair per row. Canonicalise each distinct GO value
  # once, rather than reparsing the same identifier across >1 million rows.
  go_id_lookup <- tibble(raw_go_terms = unique(raw$go_terms)) %>%
    mutate(go_terms = map(raw_go_terms, parse_go)) %>%
    unnest_longer(go_terms) %>%
    filter(!is.na(go_terms), go_terms != "")

  raw <- raw %>%
    rename(raw_go_terms = go_terms) %>%
    left_join(go_id_lookup, by = "raw_go_terms") %>%
    select(-raw_go_terms)

  if (!is.null(vv_lookup)) {
    if (all(str_detect(raw$source_gene_id, "^LOC\\d+$"))) {
      raw <- raw %>% mutate(gene_id = source_gene_id)
    } else {
      raw <- raw %>%
        left_join(vv_lookup, by = c("source_gene_id" = "Ensembl_gene")) %>%
        mutate(gene_id = case_when(
          str_detect(source_gene_id, "^LOC\\d+$") ~ source_gene_id,
          str_detect(source_gene_id, "^ENSVVLG") ~ vv_gene,
          TRUE ~ NA_character_
        ))
    }
  } else {
    raw <- raw %>% mutate(gene_id = source_gene_id)
  }

  raw %>%
    transmute(
      species = species,
      gene_id,
      go_terms
    ) %>%
    filter(!is.na(gene_id), gene_id != "", !is.na(go_terms), go_terms != "") %>%
    distinct(species, gene_id, go_terms) %>%
    mutate(source = "EXCON_EggNOG_mapper")
}

message("Reading direct wasp annotation sources")
pd_interpro <- read_interpro(files$pd_interpro, "pd")
# The first four Vespula rows have 13 fields and no GO/pathway information.
vv_interpro <- read_interpro(files$vv_interpro, "vv", skip = 4L)
message("  InterProScan annotations read")
pd_eggnog <- read_eggnog(files$pd_eggnog, "pd")
vv_eggnog <- read_eggnog(files$vv_eggnog, "vv")
message("  Galaxy EggNOG annotations read")
vv_lookup <- read_csv(files$vv_lookup, show_col_types = FALSE)

pd_go_sources <- bind_rows(
  pd_interpro$go,
  eggnog_go_long(pd_eggnog),
  read_excon_go(files$pd_excon_go, "pd")
) %>% distinct()

vv_go_sources <- bind_rows(
  vv_interpro$go,
  eggnog_go_long(vv_eggnog),
  read_excon_go(files$vv_excon_go, "vv", vv_lookup)
) %>% distinct()
message("  EXCON GO annotations read and source sets merged")

go_sources <- bind_rows(pd_go_sources, vv_go_sources) %>%
  left_join(go_lookup, by = "go_terms")

summarise_gene_go <- function(source_long) {
  source_long %>%
    group_by(gene_id) %>%
    summarise(
      GO_union_list = list(sort(unique(go_terms))),
      GO_sources = list(sort(unique(source))),
      n_GO_sources = n_distinct(source),
      .groups = "drop"
    ) %>%
    mutate(GO_pruned_list = map(GO_union_list, most_specific))
}

pd_gene_go <- summarise_gene_go(pd_go_sources)
vv_gene_go <- summarise_gene_go(vv_go_sources)

gene_go_long <- function(gene_go, species) {
  gene_go %>%
    transmute(species, gene_id, go_terms = GO_pruned_list) %>%
    unnest_longer(go_terms) %>%
    filter(!is.na(go_terms), go_terms != "") %>%
    distinct() %>%
    left_join(go_lookup, by = "go_terms")
}

pd_go_long <- gene_go_long(pd_gene_go, "pd")
vv_go_long <- gene_go_long(vv_gene_go, "vv")

source_audit <- go_sources %>%
  group_by(species, source) %>%
  summarise(
    annotated_genes = n_distinct(gene_id),
    gene_GO_pairs = n(),
    GO_IDs = n_distinct(go_terms),
    .groups = "drop"
  )

pair_source_counts <- go_sources %>%
  distinct(species, gene_id, go_terms, source) %>%
  count(species, gene_id, go_terms, name = "n_sources")

source_unique_audit <- go_sources %>%
  distinct(species, gene_id, go_terms, source) %>%
  left_join(pair_source_counts, by = c("species", "gene_id", "go_terms")) %>%
  group_by(species, source) %>%
  summarise(
    source_unique_gene_GO_pairs = sum(n_sources == 1L),
    shared_gene_GO_pairs = sum(n_sources > 1L),
    .groups = "drop"
  )

source_audit <- source_audit %>%
  left_join(source_unique_audit, by = c("species", "source"))

message("Constructing direct N13 HOG GO annotations")
hogs <- readRDS(files$hogs) %>%
  transmute(
    HOG = ortholog,
    hog_n13,
    broad_orthogroup,
    gene_pd,
    gene_vv,
    gene_dmel,
    dmel_mapping_confidence
  )

collapse_hog_go <- function(hogs, gene_col, gene_go, prefix) {
  union_name <- paste0("GO_", prefix, "_union")
  pruned_name <- paste0("GO_", prefix, "_pruned")
  genes_name <- paste0("contributing_genes_", prefix)

  hogs %>%
    select(HOG, all_of(gene_col)) %>%
    separate_rows(all_of(gene_col), sep = "\\s*,\\s*") %>%
    mutate(!!gene_col := str_trim(.data[[gene_col]])) %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    distinct(HOG, .data[[gene_col]]) %>%
    left_join(gene_go, by = setNames("gene_id", gene_col)) %>%
    group_by(HOG) %>%
    summarise(
      !!union_name := list(sort(unique(unlist(GO_union_list)))),
      !!pruned_name := list(sort(unique(unlist(GO_pruned_list)))),
      !!genes_name := paste(sort(unique(.data[[gene_col]])), collapse = ";"),
      .groups = "drop"
    )
}

pd_hog_go <- collapse_hog_go(hogs, "gene_pd", pd_gene_go, "pd")
vv_hog_go <- collapse_hog_go(hogs, "gene_vv", vv_gene_go, "vv")

hog_go <- hogs %>%
  left_join(pd_hog_go, by = "HOG") %>%
  left_join(vv_hog_go, by = "HOG") %>%
  mutate(
    across(starts_with("GO_"), ~ map(.x, empty_character)),
    GO_pair_union = map2(
      GO_pd_union, GO_vv_union,
      ~ sort(unique(c(.x, .y)))
    ),
    GO_pair_conserved = map2(GO_pd_pruned, GO_vv_pruned, function(pd_terms, vv_terms) {
      if (!length(pd_terms) || !length(vv_terms)) return(character())
      most_specific(intersect(
        go_closure_no_roots(pd_terms),
        go_closure_no_roots(vv_terms)
      ))
    })
  )

hog_go_long <- function(hog_table, list_col, annotation_set) {
  hog_table %>%
    transmute(
      HOG, hog_n13, broad_orthogroup, gene_pd, gene_vv,
      annotation_set = annotation_set,
      go_terms = .data[[list_col]]
    ) %>%
    unnest_longer(go_terms) %>%
    filter(!is.na(go_terms), go_terms != "") %>%
    distinct(HOG, go_terms, .keep_all = TRUE) %>%
    left_join(go_lookup, by = "go_terms")
}

hog_go_union_long <- hog_go_long(hog_go, "GO_pair_union", "direct_wasp_union")
hog_go_conserved_long <- hog_go_long(
  hog_go, "GO_pair_conserved", "direct_wasp_conserved"
)

message("Reading release-matched FlyBase GO associations")
gaf_columns <- c(
  "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
  "DB_Reference", "Evidence_Code", "With_From", "Aspect",
  "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon",
  "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
)

dmel_go <- read_tsv(
  files$dmel_flybase_go,
  comment = "!",
  col_names = gaf_columns,
  col_types = cols(.default = col_character()),
  quote = "",
  show_col_types = FALSE,
  progress = FALSE
) %>%
  filter(
    str_detect(DB_Object_ID, "^FBgn\\d+$"),
    !str_detect(coalesce(Qualifier, ""), "(^|\\|)NOT($|\\|)"),
    Evidence_Code != "ND"
  ) %>%
  transmute(
    gene_dmel = DB_Object_ID,
    go_terms = map(GO_ID, parse_go),
    evidence_code = Evidence_Code,
    assigned_by = Assigned_By,
    reference = DB_Reference
  ) %>%
  unnest_longer(go_terms) %>%
  filter(!is.na(go_terms), go_terms != "") %>%
  distinct()

experimental_codes <- c(
  "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
  "HTP", "HDA", "HMP", "HGI", "HEP"
)

collapse_dmel_hog_go <- function(hogs, dmel_go, evidence_filter = NULL) {
  annotations <- dmel_go
  if (!is.null(evidence_filter)) {
    annotations <- annotations %>% filter(evidence_code %in% evidence_filter)
  }

  hogs %>%
    select(HOG, gene_dmel) %>%
    separate_rows(gene_dmel, sep = "\\s*,\\s*") %>%
    filter(!is.na(gene_dmel), gene_dmel != "") %>%
    distinct() %>%
    inner_join(annotations, by = "gene_dmel", relationship = "many-to-many") %>%
    group_by(HOG) %>%
    summarise(
      GO_dmel = list(most_specific(sort(unique(go_terms)))),
      contributing_dmel_genes = paste(sort(unique(gene_dmel)), collapse = ";"),
      evidence_codes = paste(sort(unique(evidence_code)), collapse = ";"),
      .groups = "drop"
    )
}

dmel_hog_experimental <- collapse_dmel_hog_go(
  hogs, dmel_go, experimental_codes
) %>%
  rename(
    GO_dmel_experimental = GO_dmel,
    contributing_dmel_genes_experimental = contributing_dmel_genes,
    dmel_evidence_codes_experimental = evidence_codes
  )

hog_go <- hog_go %>%
  left_join(dmel_hog_experimental, by = "HOG") %>%
  mutate(
    GO_dmel_experimental = map(GO_dmel_experimental, empty_character),
    GO_final = map2(
      GO_pair_conserved, GO_dmel_experimental,
      ~ most_specific(c(.x, .y))
    )
  )

hog_go_final_wide <- hog_go %>%
  transmute(
    HOG, hog_n13, broad_orthogroup,
    gene_pd, gene_vv, gene_dmel, dmel_mapping_confidence,
    contributing_genes_pd, contributing_genes_vv,
    contributing_dmel_genes = contributing_dmel_genes_experimental,
    dmel_evidence_codes = dmel_evidence_codes_experimental,
    GO_final
  )

hog_go_final_long <- hog_go %>%
  transmute(
    HOG, hog_n13, broad_orthogroup,
    gene_pd, gene_vv, gene_dmel, dmel_mapping_confidence,
    annotation_set = "final_direct_conserved_plus_Dmel_experimental",
    go_terms = GO_final,
    GO_direct_conserved = GO_pair_conserved,
    GO_dmel_experimental,
    contributing_dmel_genes = contributing_dmel_genes_experimental,
    dmel_evidence_codes = dmel_evidence_codes_experimental
  ) %>%
  unnest_longer(go_terms) %>%
  filter(!is.na(go_terms), go_terms != "") %>%
  mutate(
    direct_wasp_supported = map2_lgl(
      go_terms, GO_direct_conserved, ~ .x %in% .y
    ),
    dmel_experimental_supported = map2_lgl(
      go_terms, GO_dmel_experimental, ~ .x %in% .y
    )
  ) %>%
  select(-GO_direct_conserved, -GO_dmel_experimental) %>%
  distinct(HOG, go_terms, .keep_all = TRUE) %>%
  left_join(go_lookup, by = "go_terms")

go_coverage_summary <- tibble(
  annotation_set = "N13 final GO: direct conserved plus Dmel experimental",
  annotated_features = n_distinct(hog_go_final_long$HOG),
  feature_term_pairs = nrow(hog_go_final_long),
  distinct_terms = n_distinct(hog_go_final_long$go_terms)
)

message("Writing final GO outputs and provenance audit")
saveRDS(
  go_sources,
  file.path(output_dir, "wasp_direct_GO_provenance_long.rds"),
  compress = "xz"
)
write_tsv(source_audit, file.path(output_dir, "GO_source_audit.tsv"))
write_tsv(hog_go_final_long, file.path(output_dir, "N13_HOG_GO_final_long.tsv.gz"))
saveRDS(hog_go_final_wide, file.path(output_dir, "N13_HOG_GO_final_wide.rds"))

rm(
  pd_go_sources, vv_go_sources, go_sources, pair_source_counts,
  source_unique_audit, pd_gene_go, vv_gene_go, pd_hog_go, vv_hog_go,
  hog_go, hog_go_union_long, hog_go_conserved_long, dmel_go,
  dmel_hog_experimental, hog_go_final_wide,
  pd_go_long, vv_go_long, pd_interpro, vv_interpro, go, go_lookup
)
invisible(gc())

message("Parsing Galaxy EggNOG KEGG and related functional annotations")
split_field <- function(x, pattern) {
  if (is.na(x) || x == "" || x == "-") return(character())
  sort(unique(str_extract_all(x, pattern)[[1]]))
}

eggnog_function_long <- function(eggnog) {
  bind_rows(
    eggnog %>% transmute(species, gene_id, annotation_type = "KEGG_KO", annotation_id = map(KEGG_ko, ~ split_field(.x, "K\\d{5}"))),
    eggnog %>% transmute(species, gene_id, annotation_type = "KEGG_pathway", annotation_id = map(KEGG_Pathway, function(x) {
      ids <- split_field(x, "(?:map|ko)\\d{5}")
      sort(unique(str_replace(ids, "^ko", "map")))
    })),
    eggnog %>% transmute(species, gene_id, annotation_type = "KEGG_module", annotation_id = map(KEGG_Module, ~ split_field(.x, "M\\d{5}"))),
    eggnog %>% transmute(species, gene_id, annotation_type = "KEGG_reaction", annotation_id = map(KEGG_Reaction, ~ split_field(.x, "R\\d{5}"))),
    eggnog %>% transmute(species, gene_id, annotation_type = "EC", annotation_id = map(EC, function(x) {
      if (is.na(x) || x == "" || x == "-") character() else sort(unique(str_split(x, ",")[[1]]))
    })),
    eggnog %>% transmute(species, gene_id, annotation_type = "BRITE", annotation_id = map(BRITE, ~ split_field(.x, "ko\\d{5}"))),
    eggnog %>% transmute(species, gene_id, annotation_type = "CAZy", annotation_id = map(CAZy, function(x) {
      if (is.na(x) || x == "" || x == "-") character() else sort(unique(str_split(x, ",")[[1]]))
    }))
  ) %>%
    unnest_longer(annotation_id) %>%
    filter(!is.na(annotation_id), annotation_id != "", annotation_id != "-") %>%
    distinct(species, gene_id, annotation_type, annotation_id) %>%
    mutate(source = "Galaxy_EggNOG_mapper_2.1.8")
}

functional_gene_long <- bind_rows(
  eggnog_function_long(pd_eggnog),
  eggnog_function_long(vv_eggnog)
)

kegg_pathway_lookup_file <- file.path(output_dir, "KEGG_pathway_TERM2NAME.tsv")
if (file.exists(kegg_pathway_lookup_file)) {
  kegg_pathway_lookup <- read_tsv(kegg_pathway_lookup_file, show_col_types = FALSE)
} else if (file.exists(file.path(legacy_output_dir, "KEGG_pathway_TERM2NAME.tsv"))) {
  kegg_pathway_lookup <- read_tsv(
    file.path(legacy_output_dir, "KEGG_pathway_TERM2NAME.tsv"),
    show_col_types = FALSE
  )
} else {
  kegg_pathway_lookup <- tryCatch({
    pathway_names <- KEGGREST::keggList("pathway")
    tibble(
      annotation_id = names(pathway_names),
      annotation_name = unname(pathway_names)
    )
  }, error = function(error) {
    warning("Could not retrieve KEGG pathway names: ", conditionMessage(error))
    tibble(annotation_id = character(), annotation_name = character())
  })
  write_tsv(kegg_pathway_lookup, kegg_pathway_lookup_file)
}

# EggNOG-mapper 2.1.8 can report pathway identifiers retired from the current
# KEGG list endpoint. Preserve their last established names for complete labels.
legacy_kegg_pathway_names <- tribble(
  ~annotation_id, ~annotation_name,
  "map01130", "Biosynthesis of antibiotics",
  "map00072", "Synthesis and degradation of ketone bodies",
  "map00471", "D-Glutamine and D-glutamate metabolism",
  "map00472", "D-Arginine and D-ornithine metabolism",
  "map00281", "Geraniol degradation"
)
kegg_pathway_lookup <- bind_rows(
  kegg_pathway_lookup,
  anti_join(legacy_kegg_pathway_names, kegg_pathway_lookup, by = "annotation_id")
) %>%
  arrange(annotation_id)
write_tsv(kegg_pathway_lookup, kegg_pathway_lookup_file)

functional_gene_long <- functional_gene_long %>%
  left_join(kegg_pathway_lookup, by = "annotation_id")

hog_gene_map <- bind_rows(
  hogs %>% select(HOG, gene_id = gene_pd) %>% mutate(species = "pd"),
  hogs %>% select(HOG, gene_id = gene_vv) %>% mutate(species = "vv")
) %>%
  separate_rows(gene_id, sep = "\\s*,\\s*") %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  distinct(HOG, species, gene_id)

hog_function_long <- hog_gene_map %>%
  inner_join(functional_gene_long, by = c("species", "gene_id")) %>%
  group_by(HOG, annotation_type, annotation_id, annotation_name, source, species) %>%
  summarise(
    contributing_genes = paste(sort(unique(gene_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  group_by(HOG, annotation_type, annotation_id, annotation_name, source) %>%
  summarise(
    species_support = paste(sort(unique(species)), collapse = ";"),
    contributing_genes_pd = paste(contributing_genes[species == "pd"], collapse = ";"),
    contributing_genes_vv = paste(contributing_genes[species == "vv"], collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    support_class = if_else(species_support == "pd;vv", "both_species", "one_species")
  )

hog_function_conserved <- hog_function_long %>%
  filter(support_class == "both_species")

kegg_pathway_hog_union <- hog_function_long %>%
  filter(annotation_type == "KEGG_pathway")
kegg_pathway_hog_conserved <- hog_function_conserved %>%
  filter(annotation_type == "KEGG_pathway")

message("Writing final KEGG and related functional annotation outputs")
write_tsv(
  pd_eggnog %>% select(species, gene_id, everything()),
  file.path(output_dir, "pd_Galaxy_EggNOG_full.tsv.gz")
)
write_tsv(
  vv_eggnog %>% select(species, gene_id, everything()),
  file.path(output_dir, "vv_Galaxy_EggNOG_full.tsv.gz")
)
write_tsv(functional_gene_long, file.path(output_dir, "wasp_functional_annotation_gene_long.tsv.gz"))
write_tsv(hog_function_long, file.path(output_dir, "N13_HOG_functional_annotation_union_long.tsv.gz"))
write_tsv(hog_function_conserved, file.path(output_dir, "N13_HOG_functional_annotation_conserved_long.tsv.gz"))

write_tsv(
  kegg_pathway_hog_union %>% select(term = annotation_id, gene = HOG) %>% distinct(),
  file.path(output_dir, "KEGG_pathway_N13_HOG_union_TERM2GENE.tsv")
)
write_tsv(
  kegg_pathway_hog_conserved %>% select(term = annotation_id, gene = HOG) %>% distinct(),
  file.path(output_dir, "KEGG_pathway_N13_HOG_conserved_TERM2GENE.tsv")
)
write_tsv(
  functional_gene_long %>%
    filter(annotation_type == "KEGG_pathway") %>%
    select(species, term = annotation_id, gene = gene_id) %>%
    distinct(),
  file.path(output_dir, "KEGG_pathway_wasp_gene_TERM2GENE.tsv")
)

coverage_summary <- bind_rows(
  go_coverage_summary,
  tibble(
    annotation_set = c(
      "N13 KEGG pathway union", "N13 KEGG pathway conserved"
    ),
    annotated_features = c(
      n_distinct(kegg_pathway_hog_union$HOG),
      n_distinct(kegg_pathway_hog_conserved$HOG)
    ),
    feature_term_pairs = c(
      nrow(kegg_pathway_hog_union), nrow(kegg_pathway_hog_conserved)
    ),
    distinct_terms = c(
      n_distinct(kegg_pathway_hog_union$annotation_id),
      n_distinct(kegg_pathway_hog_conserved$annotation_id)
    )
  )
)
write_tsv(coverage_summary, file.path(output_dir, "annotation_coverage_summary.tsv"))

message("Saving reloadable workspace and key data frames")
key_rds_dir <- file.path(output_dir, "key_dataframes_rds")
dir.create(key_rds_dir, recursive = TRUE, showWarnings = FALSE)

key_dataframes <- list(
  go_source_audit = source_audit,
  annotation_coverage_summary = coverage_summary,
  n13_hog_go_final_long = hog_go_final_long,
  wasp_functional_annotation_gene_long = functional_gene_long,
  n13_hog_functional_annotation_union_long = hog_function_long,
  n13_hog_functional_annotation_conserved_long = hog_function_conserved,
  kegg_pathway_n13_hog_union_term2gene = kegg_pathway_hog_union %>%
    select(term = annotation_id, gene = HOG) %>%
    distinct(),
  kegg_pathway_n13_hog_conserved_term2gene = kegg_pathway_hog_conserved %>%
    select(term = annotation_id, gene = HOG) %>%
    distinct(),
  kegg_pathway_term2name = kegg_pathway_lookup
)

walk2(
  key_dataframes,
  names(key_dataframes),
  ~ saveRDS(.x, file.path(key_rds_dir, paste0(.y, ".rds")), compress = "gzip")
)

workspace_env <- list2env(key_dataframes, parent = emptyenv())
workspace_file <- file.path(output_dir, "final_GO_KEGG_analysis_workspace.RData")
save(
  list = names(key_dataframes),
  file = workspace_file,
  envir = workspace_env,
  compress = "gzip"
)

workspace_index <- tibble(
  object_name = names(key_dataframes),
  rows = map_int(key_dataframes, nrow),
  columns = map_int(key_dataframes, ncol),
  individual_rds = file.path(
    "key_dataframes_rds", paste0(names(key_dataframes), ".rds")
  )
)
write_tsv(workspace_index, file.path(output_dir, "workspace_file_index.tsv"))

writeLines(
  c(
    "Final GO and KEGG annotation",
    "",
    "Direct wasp GO combines Galaxy InterProScan, Galaxy EggNOG-mapper and",
    "EXCON EggNOG-mapper by set union at gene level. For each N13 HOG, direct",
    "terms are retained through the established Polistes-Vespula conserved-GO",
    "rule and combined with experimentally supported GO terms transferred from",
    "the mapped Drosophila melanogaster orthologue.",
    "",
    "The final GO table is N13_HOG_GO_final_long.tsv.gz. It uses EXP, IDA, IPI,",
    "IMP, IGI, IEP and their high-throughput equivalents from FlyBase FB2025_03,",
    "excluding NOT and ND annotations. Alternative GO sensitivity tables are",
    "not exported by this final script. Source provenance remains auditable.",
    "",
    "For KEGG pathway enrichment with clusterProfiler::enricher(), use one of",
    "the KEGG_pathway_*_TERM2GENE.tsv files plus KEGG_pathway_TERM2NAME.tsv.",
    "The union table maximises coverage; the conserved table requires the same",
    "pathway annotation in both wasp species represented by an N13 HOG.",
    "InterProScan contributes GO annotations only in this build; its separate",
    "MetaCyc/Reactome pathway field is not treated as KEGG annotation.",
    "InterPro rows with 13 rather than 15 fields lack optional GO/pathway fields;",
    "the parser audit found no shifted or overfilled records.",
    "Reload all key tables with load('final_GO_KEGG_analysis_workspace.RData')."
  ),
  file.path(output_dir, "README.txt")
)

validation_checks <- c(
  final_GO_nonempty = nrow(hog_go_final_long) > 0L,
  no_duplicate_HOG_GO_pairs = !anyDuplicated(
    hog_go_final_long[c("HOG", "go_terms")]
  ),
  canonical_GO_IDs = all(str_detect(
    hog_go_final_long$go_terms, "^GO:\\d{7}$"
  )),
  full_GO_names = !any(str_detect(hog_go_final_long$term_name, fixed("..."))),
  final_GO_exported = file.exists(file.path(
    output_dir, "N13_HOG_GO_final_long.tsv.gz"
  )),
  KEGG_exported = file.exists(file.path(
    output_dir, "KEGG_pathway_N13_HOG_union_TERM2GENE.tsv"
  )),
  workspace_exported = file.exists(workspace_file)
)
if (!all(validation_checks)) {
  stop(
    "Final annotation validation failed: ",
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

message("Completed: ", output_dir)
