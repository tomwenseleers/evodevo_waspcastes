# Supplementary R code associated with the article "Social insect castes evolved by re-timing ancient, seasonally expressed pathways"
## Anneline Mattens¹*, Viviana Di Pietro¹*, Cintia Akemi Oi², Ricardo Caliari Oliveira³, Chris Wyatt², Fabio Santos do Nascimento⁴ & Tom Wenseleers¹
### ¹ KU Leuven, Leuven, Belgium; ² The Sumner Lab, University College London, London, UK; ³ Universitat Autònoma de Barcelona, Barcelona, Spain; ⁴ Universidade de São Paulo, Ribeirão Preto, Brazil; *These authors contributed equally.; Correspondence: anneline.mattens@kuleuven.be; annelinemattens1@gmail.com; tom.wenseleers@kuleuven.be


This repository contains the code to reproduce the results presented in the article "Social insect castes evolved by re-timing ancient, seasonally expressed pathways", in which we show that queen-worker caste differentiation in vespine wasps recapitulates ancestral seasonally expressed pathways. Included are the final R workflows, compact reloadable analysis objects, required local inputs, and manuscript outputs (figures, supplementary figures and supplementary tables) related to functional annotation and the comparative transcriptomic and gene-family analyses.

All paths are defined relative to the repository root.

## Quick start

The three root scripts are the public entry points:

1. `A_annotation_FINAL.R`
2. `B_differential expression analysis_FINAL.R`
3. `C_EXCON CAFE analysis_FINAL.R`

Each begins with an editable `rerun <- FALSE` flag. With this default, the
script loads the corresponding compact RDS/RData objects without recomputing
the analysis. Set `rerun <- TRUE` to regenerate the derived outputs. The same
choice can be set non-interactively with the `EVODEVO_RERUN` environment
variable, for example in PowerShell:

```powershell
$env:EVODEVO_RERUN = "true"
Rscript "B_differential expression analysis_FINAL.R"
```

The B workflow reruns the downstream PLS, regression, GO, table, and figure
analyses from saved gene- and HOG-level model checkpoints. It deliberately does
not refit thousands of gene-wise mixed models. Full refits are available as the
two one-worker scripts in `scripts/upstream/` and should be run sequentially.

## Directory structure

```text
.
|-- A_annotation_FINAL.R
|-- B_differential expression analysis_FINAL.R
|-- C_EXCON CAFE analysis_FINAL.R
|-- input_annotation/
|   |-- raw/                 # InterProScan, EggNOG, FlyBase and GO inputs
|   |-- reference/           # N13 HOG and identifier reference objects
|   `-- workspaces/          # compact final annotation workspace
|-- output_annotation/       # final GO, KEGG and functional annotations
|-- input_differential_expression/
|   |-- raw/                 # Salmon counts, TPMs, lengths, metadata and N13.tsv
|   |-- checkpoints/         # fitted gene- and N13-HOG glmmTMB objects
|   `-- workspaces/          # compact downstream key-results RDS
|-- input_CAFE/
|   |-- wasp_CAFE/           # CAFE input/count/result files used downstream
|   |-- annotation/          # gene tables, curated HOG classes and BUSCO table
|   |-- GO/                  # HOG-to-GO mapping
|   |-- orthology/           # Ancistrocerus orthology annotation evidence
|   |-- tree_dating/         # dated wasp tree and calibration summary
|   `-- workspaces/          # compact CAFE/GO key-results RDS
|-- scripts/
|   |-- pipelines/           # complete analysis implementations
|   |-- upstream/            # costly model refits and dated-tree reconstruction
|   `-- utilities/           # supplementary-table workbook/export utilities
`-- output/
    |-- figures/             # Fig1, Fig2 and Fig3
    |-- supplemental_figures/# FigS1, FigS2 and FigS3
    |-- supplemental_tables/ # TableS1-TableS13 (S5 has A and B parts)
    `-- ...                  # complete derived analysis and audit tables
```

Large TSV files are not included in this repository, but can be recreated from the compact R objects by rerunning the relevant downstream script. 

## A. Functional annotation

`A_annotation_FINAL.R` integrates direct *Polistes dominula* and *Vespula
vulgaris* annotations from two EggNOG-mapper analyses and InterProScan. It also
transfers release-matched, experimentally supported FlyBase GO annotations
through *Drosophila melanogaster* members of the shared N13 hierarchical
orthogroups. The final N13-HOG GO set combines direct GO terms conserved across
the two wasps with these orthology-transferred FlyBase terms. KEGG and related
functional fields are imported from the Galaxy EggNOG-mapper tables.

The compact final mapping consumed by B is
`output_annotation/N13_HOG_GO_final_long.tsv.gz`. Full source-level provenance
is stored efficiently as `output_annotation/wasp_direct_GO_provenance_long.rds`.

## B. Differential expression and cross-species analyses

Raw Salmon gene counts, sample-specific effective gene lengths, TPMs, and
sample metadata are under `input_differential_expression/raw/`. At both gene
and N13-HOG levels, the fitted model was a negative-binomial (`nbinom2`)
`glmmTMB` model with stage-specific condition effects, a colony-within-condition
random intercept, and a gene/HOG-by-sample offset combining library-size and
effective-length normalisation. Original Salmon estimated counts were rounded
only at the count-model interface. Wald contrasts tested late versus early
season within each *P. dominula* stage and queen versus worker within each
*V. vulgaris* stage. BH FDR correction was applied globally across all stage-by-gene
or stage-by-HOG tests within each species, and coefficients were shrunk with
ASH separately by species and stage.

For cross-species analyses, genes were summarised directly to OrthoFinder N13
HOGs before model fitting. The final PLS, robust stage-wise regressions, and
nonnegative ridge heterochrony model use ASH-shrunken N13-HOG log2 fold changes.
Figure 1 GO foregrounds comprise the top 300 HOGs in each oriented PLS loading
direction. GO enrichment uses topGO's `weight01` algorithm with Fisher's exact
test and no post-hoc multiple-testing correction to the topology-adjusted
topGO values.

## C. EXCON, OrthoFinder, EggNOG-mapper and CAFE

The sequence-level steps were run externally with the Eco-Flow EXCON Nextflow
pipeline. EXCON's OrthoFinder stage generated the hierarchical orthogroups, and
the same orthology output was reused throughout this study. N13 and N0 identify
different hierarchical levels in that output: N13 HOGs are used for the
two-species transcriptomic comparison, whereas root-level N0 HOGs and resolved
gene trees are used for CAFE. HOG identifiers are analysis-specific and should
not be treated as stable identifiers across independent OrthoFinder reruns.

The orthology/annotation stage was run in the following form (machine-specific
paths omitted):

```bash
nextflow run main.nf -resume -profile local,docker \
  --input input.csv \
  --outdir results_orthology \
  --run_eggnog \
  --eggnog_data_dir eggnog_data \
  --max_cpus 24 --max_memory 96.GB
```

The wasp tree was pruned from the OrthoFinder species tree, re-estimated using
IQ-TREE 3, and dated with `scripts/upstream/make_dated_trees.R`. The dating
workflow uses `ape::chronos` with a discrete-rate model and literature-based
Aculeata, Apoidea, Vespidae, and Polistinae+Vespinae calibrations. The dated
Vespidae tree contains 15 Vespidae plus *Tiphia femorata* and *Anoplius
nigerrimus* as external outgroups and was passed to EXCON without rescaling.

The CAFE stage was run in the following form:

```bash
nextflow run main.nf -profile singularity --custom_config cluster.config \
  --input_tree Vespidae_with_outgroups_dated_primary.nwk \
  --input_tree_is_dated \
  --input_orthogroups N0.tsv \
  --cafe_focus_clades \
    'Vespa_crabro,Vespula_vulgaris|Vespa_crabro,Polistes_dominula|Polistes_dominula,Mischocyttarus_mexicanus' \
  --orthofinder_genetree_dir Resolved_Gene_Trees \
  --outdir wasp_CAFE \
  --cafe_max_differential 20 --cafe_filter_first -bg -resume
```

Standard families were fitted jointly. High-differential families were fitted
individually using the same fixed global lambda, Gamma settings, and error
model. The final downstream screen requires both a CAFE family-wide
`P < 0.05` and focal-branch Viterbi probability `< 0.01`. These nominal values
were not FDR-adjusted. Conservatively identified transposable-element HOGs were
retained in the fitted CAFE model but excluded from reported biological totals,
functional classification, GO foregrounds, and GO backgrounds. Families
removed after an inferred root size of zero are reported separately because no
family-wide or branch-specific CAFE result can be assigned to them.

Node-wise GO enrichment tests significantly expanded, non-TE HOGs at Node 20
(stem Polistinae+Vespinae) and Node 25 (stem Vespinae) against all GO-annotated,
tested non-TE HOGs. BP, MF, and CC are tested with topGO `weight01`/Fisher,
without FDR correction and without truncating the ranked ontology output.

## Manuscript outputs

- `output/figures/Fig1.*`: PLS developmental/season-caste structure and GO enrichment.
- `output/figures/Fig2.*`: cross-species stage-wise regressions and heterochrony model.
- `output/figures/Fig3.*`: dated wasp phylogeny with focal CAFE changes and expansion classes.
- `output/supplemental_figures/FigS1.png`: sampling and developmental-stage overview.
- `output/supplemental_figures/FigS2.*`: numbers of differentially expressed genes by stage.
- `output/supplemental_figures/FigS3.*`: stage-specific absolute shrunken log2 fold-change heatmap.
- `output/supplemental_tables/TableS1.tsv` to `TableS13.tsv`: final supplementary tables; Table S5 is split into `TableS5A.tsv` and `TableS5B.tsv`.

Tables S1-S3 and S9 summarise sampling, microsatellites, genome quality, and
tree calibrations. B generates Tables S4-S8. C generates Tables S10-S13.

## Software requirements

The analyses were validated with R 4.5.x. Principal R dependencies are:
`dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `KEGGREST`, `ontologyIndex`,
`edgeR`, `glmmTMB`, `DESeq2`, `tximport`, `ashr`, `mixOmics`, `topGO`, `GO.db`,
`AnnotationDbi`, `estimatr`, `robustbase`, `glmnet`, `ape`, `export`, `officer`,
`pheatmap`, and `openxlsx`. The supplementary XLSX builder additionally uses
Python 3 with `openpyxl`.

`FILE_MANIFEST.tsv` records the relative path, size, SHA-256 checksum, and role
of the packaged files (excluding the manifest itself). Because model
checkpoints and compressed raw annotation files are large, allow sufficient
disk space and clone time.

## Use of artificial intelligence

OpenAI Codex/GPT 5.6 was used as an AI-assisted coding and editing tool to help refactor and document R scripts, organise this repository, troubleshoot software execution, and prepare preliminary table and figure exports. The authors specified the study design, statistical models and tools to use, significance criteria, and biological questions; reviewed and edited all AI-assisted code and text; and reran and validated the reported workflows. Generative AI was not used to create or alter the underlying research data, and it did not independently determine the biological interpretations or conclusions. The authors take full responsibility for the accuracy and reproducibility of all analyses and repository contents.
