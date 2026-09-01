# Supplementary R code associated with the article “Social insect castes evolved by re-timing ancient, seasonally expressed pathways”

**Anneline Mattens<sup>1,\*</sup>, Viviana Di Pietro<sup>1,\*</sup>, Marthe Van Marsenille<sup>1,\*</sup>, Chris Wyatt<sup>2</sup>, Ricardo Caliari Oliveira<sup>3</sup>, Cintia Akemi Oi<sup>2</sup>, Fabio Santos do Nascimento<sup>4</sup> and Tom Wenseleers<sup>1</sup>**

<sup>1</sup> KU Leuven, Leuven, Belgium  
<sup>2</sup> The Sumner Lab, University College London, London, UK  
<sup>3</sup> Universitat Autònoma de Barcelona, Barcelona, Spain  
<sup>4</sup> Universidade de São Paulo, Ribeirão Preto, Brazil  

<sup>\*</sup> These authors contributed equally.

**Correspondence:** [anneline.mattens@kuleuven.be](mailto:anneline.mattens@kuleuven.be), [annelinemattens1@gmail.com](mailto:annelinemattens1@gmail.com) or [tom.wenseleers@kuleuven.be](mailto:tom.wenseleers@kuleuven.be)

This repository contains the code to reproduce the results presented in the article "Social insect castes evolved by re-timing ancient, seasonally expressed pathways", in which we show that queen-worker caste differentiation in vespine wasps recapitulates ancestral seasonally expressed pathways. Included are the final R workflows, compact reloadable analysis objects, required local inputs, and manuscript outputs (figures, supplementary figures and supplementary tables) related to functional annotation and the comparative transcriptomic and gene-family analyses.

All paths are defined relative to the repository root.

## Quick start

The three root scripts are the public entry points:

1. `A_annotation.R`
2. `B_differential expression analysis.R`
3. `C_EXCON CAFE analysis.R`

Each begins with an editable `rerun <- FALSE` flag. With this default, the
script loads the corresponding compact RDS/RData objects without recomputing
the analysis. Set `rerun <- TRUE` to regenerate the derived outputs. The same
choice can be set non-interactively with the `EVODEVO_RERUN` environment
variable, for example in PowerShell:

```powershell
$env:EVODEVO_RERUN = "true"
Rscript "B_differential expression analysis.R"
```

The B workflow reruns the downstream PLS, regression, GO, table, and figure
analyses from saved gene- and HOG-level model checkpoints. It deliberately does
not refit thousands of gene-wise mixed models. Full refits are available as the
two one-worker scripts in `scripts/upstream/` and should be run sequentially.

## Directory structure

```text
.
|-- A_annotation.R
|-- B_differential expression analysis.R
|-- C_EXCON CAFE analysis.R
|-- nextflow_runs/
|   |-- 1_rnaseq/            # nf-core/rnaseq inputs, launcher and retained outputs
|   `-- 2_EXCON/
|       |-- 1_download_genomes/                    # genome preparation and BUSCO summary
|       |-- 2_EXCON_orthofinder_eggnogmapper_run/ # EXCON orthology/annotation run
|       `-- 3_EXCON_CAFE_run/                     # wasp EXCON/CAFE run
|-- input_annotation/
|   |-- raw/                 # Galaxy, InterProScan, FlyBase and GO inputs
|   |-- reference/           # N13 HOG and identifier reference objects
|   `-- workspaces/          # compact final annotation workspace
|-- output_annotation/       # final GO, KEGG and functional annotations
|-- input_differential_expression/
|   |-- raw/                 # study-specific sample metadata
|   |-- checkpoints/         # fitted gene- and N13-HOG glmmTMB objects
|   `-- workspaces/          # compact downstream key-results RDS
|-- input_CAFE/
|   |-- annotation/          # gene tables and curated HOG classifications
|   |-- GO/                  # HOG-to-GO mapping
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

The retained Nextflow directories preserve substantially more of the original
pipeline structure than is strictly required by the R workflows. The R scripts
read Salmon matrices, EXCON EggNOG-mapper GO files, OrthoFinder HOG tables,
BUSCO summaries and CAFE results directly from these directories, avoiding
duplicate copies under the downstream `input_*` folders. Compact R objects are
also included so that the final analyses can be inspected without rerunning the
computationally expensive upstream pipelines.

## A. Functional annotation

`A_annotation.R` integrates direct *Polistes dominula* and *Vespula
vulgaris* annotations from two EggNOG-mapper analyses and InterProScan. It also
transfers release-matched, experimentally supported FlyBase GO annotations
through *Drosophila melanogaster* members of the shared N13 hierarchical
orthogroups. The final N13-HOG GO set combines direct GO terms conserved across
the two wasps with these orthology-transferred FlyBase terms. KEGG and related
functional fields are imported from the Galaxy EggNOG-mapper tables.

The EXCON EggNOG-mapper GO files are read directly from
`nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/results_EXCON/eggnogmapper/go_files/`;
the non-EXCON annotation sources remain under `input_annotation/raw/`.

The compact final mapping consumed by B is
`output_annotation/N13_HOG_GO_final_long.tsv.gz`. Full source-level provenance
is stored efficiently as `output_annotation/wasp_direct_GO_provenance_long.rds`.

## B. Differential expression and cross-species analyses

RNA-seq reads from each species were processed in separate nf-core/rnaseq
v3.14.0 runs under Nextflow v24.10.1. The pipeline used its `star_salmon`
workflow, with STAR v2.7.9a for alignment and Salmon v1.10.1 for transcript
quantification. Paired-end library strandedness was inferred from the reads,
and `--featurecounts_group_type transcript_id` was supplied in both runs. The
complete samplesheets, parameter files, software-version records, execution
reports and retained outputs are under `nextflow_runs/1_rnaseq/`. The two runs
can be launched in the same form with:

```bash
bash nextflow_runs/1_rnaseq/run_rnaseq_both_species.sh
```

The launcher resolves its working directories relative to its own location and
pins nf-core/rnaseq v3.14.0 by default. The archived samplesheets retain the
original FASTQ paths, while reference FASTA and GTF files are not duplicated in
this repository; these paths must therefore be supplied for a new installation
as described in the launcher comments.

Raw Salmon gene counts, sample-specific effective gene lengths and TPMs are
read directly from the species-specific `star_salmon/` directories under
`nextflow_runs/1_rnaseq/output_files/`. Study-specific sample metadata are
under `input_differential_expression/raw/`. At both gene
and N13-HOG levels, the fitted model was a negative-binomial (`nbinom2`)
`glmmTMB` model with stage-specific condition effects, a colony-within-condition
random intercept, and a gene/HOG-by-sample offset combining library-size and
effective-length normalisation. Original Salmon estimated counts were rounded
only at the count-model interface. Wald contrasts tested late versus early
season within each *P. dominula* stage and queen versus worker within each
*V. vulgaris* stage. BH correction was applied globally across all stage-by-gene
or stage-by-HOG tests within each species, and coefficients were shrunk with
ASH separately by species and stage.

Complete machine-readable results for every fitted contrast, rather than only
the significant subsets, are exported as
`output/full_gene_level_differential_expression_results.tsv.gz` and
`output/full_N13_HOG_level_differential_expression_results.tsv.gz`. Both tables
contain raw and ASH-shrunken effect estimates, Wald statistics, nominal p
values, stage-wise BH values, the primary global BH value, and an explicit
`global_FDR_lt_0.05` indicator. The gene-level table additionally contains gene
symbols, descriptions, gene types and N13-HOG assignments; the HOG-level table
contains member genes and aggregated gene annotations for both wasps. The
compressed TSV files can be read directly with `readr::read_tsv()`.

For cross-species analyses, genes were summarised directly to OrthoFinder N13
HOGs before model fitting. The final PLS, robust stage-wise regressions, and
nonnegative ridge heterochrony model use ASH-shrunken N13-HOG log2 fold changes.
Figure 1 GO foregrounds comprise the top 300 HOGs in each oriented PLS loading
direction. GO enrichment uses topGO's `weight01` algorithm with Fisher's exact
test. Given the dependencies in the GO ontology tree, and as recommended by the topGO package authors, these p values are not FDR corrected.

## C. EXCON, OrthoFinder, EggNOG-mapper and CAFE

The sequence-level steps were run externally with the Eco-Flow EXCON Nextflow
pipeline. EXCON's OrthoFinder stage generated the hierarchical orthogroups, and
the same orthology output was reused throughout this study. N13 and N0 identify
different hierarchical levels in that output: N13 HOGs are used for the
two-species transcriptomic comparison, whereas root-level N0 HOGs and resolved
gene trees are used for CAFE. HOG identifiers are analysis-specific and should
not be treated as stable identifiers across independent OrthoFinder reruns.

The EXCON orthology/EggNOG-mapper run is under
`nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/`. The wasp CAFE
run is under `nextflow_runs/2_EXCON/3_EXCON_CAFE_run/`. The downstream CAFE R
workflow reads the retained `results_EXCON/` output directly; the N13
transcriptomic workflow similarly reads `N13.tsv` directly from the retained
OrthoFinder hierarchical-orthogroup directory. The consolidated
*Ancistrocerus nigricornis* orthology table used as an annotation fallback is
retained in the same run under `results_EXCON/orthofinder/Orthologues/`.

The archived orthology/annotation run report records Nextflow v25.10.0,
OrthoFinder v3.1.3, eggNOG-mapper v2.1.13, AGAT v1.6.1 and gffread v0.12.7. The CAFE execution report records EXCON v2.3.2 (script ID
`5862d5fe137d0a3c587c84128848e3e0`) under Nextflow v26.04.6, with CAFE v4.2.1,
cafeplotter v0.2.0 and R v4.3.1. The retained `pipeline_info/` directories are
the authoritative provenance records for these upstream runs.

The orthology/annotation stage was run in the following form (machine-specific
paths omitted):

```bash
nextflow run main.nf -resume -profile local,docker \
  --input input.csv \
  --outdir results_EXCON \
  --run_eggnog \
  --eggnog_data_dir eggnog_data \
  --max_cpus 24 --max_memory 96.GB
```

The repository-relative launchers are
`nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/2.run_excon_pipeline.sh`
and `nextflow_runs/2_EXCON/3_EXCON_CAFE_run/2.run_excon_CAFE_pipeline.sh`.
Genome downloads, longest-isoform extraction, BUSCO assessment and regeneration
of `input_excon.csv` are likewise scripted under
`nextflow_runs/2_EXCON/1_download_genomes/` and the orthology-run directory.
Machine-specific executables, container profiles, data caches and cluster
configuration files can be supplied through the environment variables
documented at the start of each launcher.

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
  --outdir results_EXCON \
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
tested non-TE HOGs. BP, MF, and CC are tested with topGO `weight01`/Fisher. Given the dependencies in the GO ontology tree, and as recommended by the topGO package authors, these p values are not FDR corrected.

## Manuscript outputs

- `output/figures/Fig1.*`: PLS developmental/season-caste structure and GO enrichment.
- `output/figures/Fig2.*`: cross-species stage-wise regressions and heterochrony model.
- `output/figures/Fig3.*`: dated wasp phylogeny with focal CAFE changes and expansion classes.
- `output/supplemental_figures/FigS1.png`: sampling and developmental-stage overview.
- `output/supplemental_figures/FigS2.*`: numbers of differentially expressed genes by stage.
- `output/supplemental_figures/FigS3.*`: stage-specific absolute shrunken log2 fold-change heatmap.
- `output/supplemental_tables/TableS1.tsv` to `TableS13.tsv`: final supplementary tables; Table S5 is split into `TableS5A.tsv` and `TableS5B.tsv`.
- `output/full_gene_level_differential_expression_results.tsv.gz`: complete annotated gene-by-stage differential-expression results for both species.
- `output/full_N13_HOG_level_differential_expression_results.tsv.gz`: complete annotated N13-HOG-by-stage differential-expression results for both species.
- `output/full_DE_results_file_index.tsv`: row counts, significant-test counts and compressed sizes for the two complete DE tables.

Tables S1-S3 and S9 summarise sampling, microsatellites, genome quality, and
tree calibrations. B generates Tables S4-S8. C generates Tables S10-S13.

## Software requirements

The analyses were validated with R 4.5.x. Principal R dependencies are:
`dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `KEGGREST`, `ontologyIndex`,
`edgeR`, `glmmTMB`, `DESeq2`, `tximport`, `ashr`, `mixOmics`, `topGO`, `GO.db`,
`AnnotationDbi`, `estimatr`, `robustbase`, `glmnet`, `ape`, `export`, `officer`,
`pheatmap`, and `openxlsx`. The genome-download helper uses Python 3 with
`requests` and `beautifulsoup4`; the supplementary XLSX builder additionally
uses Python 3 with `openpyxl`.

Because model checkpoints and compressed raw annotation files are large, allow
sufficient disk space and clone time.

## Use of artificial intelligence

OpenAI Codex/GPT 5.6 was used as an AI-assisted coding and editing tool to help refactor and document R scripts, organise this repository, troubleshoot software execution, and prepare preliminary table and figure exports. The authors specified the study design, statistical models and tools to use, significance criteria, and biological questions; reviewed and edited all AI-assisted code and text; and reran and validated the reported workflows. Generative AI was not used to create or alter the underlying research data, and it did not independently determine the biological interpretations or conclusions. The authors take full responsibility for the accuracy and reproducibility of all analyses and repository contents.
