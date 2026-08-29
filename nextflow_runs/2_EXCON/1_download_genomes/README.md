# Download_genomes – genome & BUSCO prep pipeline

This folder contains all the scripts to:

1. Download genomes + annotations from Ensembl and/or NCBI.
2. Build non-redundant protein sets (one longest isoform per gene).
3. Compare Ensembl vs NCBI genomes with BUSCO for species where you have both.
4. Run BUSCO on the final chosen proteomes for all species.
5. Aggregate BUSCO outputs into a single summary table.

The end result is a curated set of BUSCO-scored proteomes that you can safely feed into Excon.

---

## Folder layout (key paths)

- `Download_genomes/`  ← **this folder**
  - `1.download_genomes_all.py`
  - `1.download_genomes_all.sh`
  - `Genomes_selection_ensembl.csv`
  - `Genomes_selection_ncbi.csv`
  - `Genomes_selection_both.csv`
  - `2.build_proteins_longest_all.sh`
  - `3.run_busco_all.sh`
  - `4.busco_annotations_aggregate.sh`
  - `4.busco_annotations_aggregate.py`
  - `genomes_out/`
    - `ensembl/<Species>/...`
    - `ncbi/<Species>/...`
    - `both/<Species>/...`
  - (later) `busco_annotation_all/` with all BUSCO runs
  - (later) `download_summary_all.tsv`
  - (later) `000_BUSCO_annotations_summary.tsv`

All commands below assume you are in:

```bash
/mnt/DATA/Genome_data_wasps/4.Excon_analysis/Download_genomes
```

---

## Overview of the workflow

1. **Step 1 – Download genomes/annotations**  
   Use `1.download_genomes_all.py` (via `1.download_genomes_all.sh`) to download genomes and annotations from Ensembl, NCBI, or both, based on the three `Genomes_selection_*.csv` tables. Output goes to `genomes_out/ensembl`, `genomes_out/ncbi`, and `genomes_out/both`.

2. **Step 2 – Build protein sets & keep longest isoform**  
   Run `2.build_proteins_longest_all.sh` to create protein FASTA files for each genome and collapse isoforms to a single “longest isoform per gene” proteome.

3. **Step 3 – Compare Ensembl vs NCBI with BUSCO (species with both)**  
   In `genomes_out/`, run `run_busco_both.sh`. This runs BUSCO on the Ensembl- and NCBI-derived proteomes for species in `genomes_out/both/` and helps decide which source to keep.

4. **Step 4 – Run BUSCO for the final choice for *all* species**  
   Back in the main `Download_genomes` directory, run `3.run_busco_all.sh`. This runs BUSCO on the chosen proteome per species (Ensembl *or* NCBI) and stores all runs in `busco_annotation_all/`.

5. **Step 5 – Aggregate BUSCO results**  
   Finally, run `4.busco_annotations_aggregate.sh`, which calls `4.busco_annotations_aggregate.py` and aggregates all BUSCO runs into `000_BUSCO_annotations_summary.tsv`.

Below, each step is explained in more detail: **what it does, why it exists, and how to run it.**

---

## Step 1 – Download genomes + annotations

### Files

- `1.download_genomes_all.py` – Python script that actually talks to Ensembl/NCBI.
- `1.download_genomes_all.sh` – small shell wrapper that calls the Python script with the correct CSVs and output folder.
- `Genomes_selection_ensembl.csv`
- `Genomes_selection_ncbi.csv`
- `Genomes_selection_both.csv`

### What & why

**Goal:**  
Create a *clean, unified* local collection of genomes + annotations, with a consistent folder structure, starting from your curated selection tables.

You may want some species **only from Ensembl**, some **only from NCBI**, and some **from both** sources (for later comparison). The script hides all the logic of Ensembl/NCBI FTP / dataset APIs, compressed files, and inconsistent file names.

`1.download_genomes_all.py`:

- Reads the three selection tables:
  - `Genomes_selection_ensembl.csv` – species to fetch from Ensembl.
  - `Genomes_selection_ncbi.csv` – species to fetch from NCBI (requires GCF/GCA columns).
  - `Genomes_selection_both.csv` – species where you want **both Ensembl + NCBI**.
- For each row, it:
  - Normalises the species name and builds a safe folder name (underscores).
  - For Ensembl:
    - Crawls the Ensembl organisms FTP for that species.
    - Tries to prefer **soft-masked** genome FASTA (if available).
    - Downloads a genome FASTA and a GFF3/GTF annotation and uncompresses them.
  - For NCBI:
    - Uses either the `datasets` CLI (if available) or the NCBI REST API.
    - Downloads a zip containing genome FASTA + GFF/GTF.
    - Extracts and renames files to something stable like  
      `GCF_xxx_yyy_genomic.fna` and `GCF_xxx_yyy_genomic.gff/gtf`.
- Fills a **summary table**: `genomes_out/download_summary_all.tsv`.

Directory layout (simplified):

```text
genomes_out/
  ensembl/
    Species_1/
      <Ensembl FASTA & GFF/GTF as downloaded>
  ncbi/
    Species_2/
      GCF_..._genomic.fna
      GCF_..._genomic.gff or .gtf
  both/
    Species_3/
      ensembl_softmasked.fa OR ensembl_genome.fa
      ensembl_annotation.gff3 or .gtf
      ncbi_GCF_..._genomic.fna
      ncbi_GCF_..._genomic.gff or .gtf
  download_summary_all.tsv
```

### How to run

Typical commands (from `Download_genomes/`):

```bash
# Ensembl only
python3 1.download_genomes_all.py   --ensembl_csv Genomes_selection_ensembl.csv   --out genomes_out

# NCBI only
python3 1.download_genomes_all.py   --ncbi_csv Genomes_selection_ncbi.csv   --out genomes_out

# Species with both
python3 1.download_genomes_all.py   --both_csv Genomes_selection_both.csv   --out genomes_out

# Or everything in one go
python3 1.download_genomes_all.py   --ensembl_csv Genomes_selection_ensembl.csv   --ncbi_csv Genomes_selection_ncbi.csv   --both_csv Genomes_selection_both.csv   --out genomes_out
```

Usually, `1.download_genomes_all.sh` wraps one of the “everything in one go” calls so you just do:

```bash
bash 1.download_genomes_all.sh
```

---

## Step 2 – Build protein files & keep the longest isoform

### File

- `2.build_proteins_longest_all.sh`

### What & why

**Goal:**  
For each genome + annotation pair downloaded in Step 1, generate a **clean protein FASTA** where each gene is represented by **one protein: the longest isoform**.

Why this is needed:

- BUSCO (and Excon) work best with a *non-redundant* proteome:
  - Without this step, genes with many isoforms would dominate statistics.
  - Keeping only the longest isoform per gene is a common standard for comparative analyses.
- It also ensures that Ensembl and NCBI proteomes are treated in a consistent way.

What the script does conceptually:

1. Loops over all species subfolders in `genomes_out/ensembl`, `genomes_out/ncbi` and `genomes_out/both`.
2. For each species and source, it:
   - Uses the genome FASTA + annotation (GFF3/GTF) to generate protein sequences (or uses existing protein FASTA if already present, depending on how you wrote it).
   - Groups transcripts by gene and **keeps the longest protein isoform** per gene.
3. Writes out something along the lines of:

```text
genomes_out/
  ensembl/Species/
    ENSEMBL_prot_long.faa        # one protein per gene
  ncbi/Species/
    NCBI_prot_long.faa           # one protein per gene
  both/Species/
    ENSEMBL_prot_long.faa
    NCBI_prot_long.faa
```

(Exact file names may differ slightly; see the header of `2.build_proteins_longest_all.sh` for the precise naming.)

These `*_prot_long` files are the **input for all BUSCO runs** in the next steps.

### How to run

From `Download_genomes/`:

```bash
bash 2.build_proteins_longest_all.sh
```

---

## Step 3 – BUSCO comparison for species with Ensembl + NCBI

### File

- `run_busco_both.sh` (run from inside `genomes_out/`)

### What & why

**Goal:**  
For species where you downloaded **both** Ensembl and NCBI genomes, decide **which source is better** based on BUSCO scores.

Why:

- Some species have multiple assemblies, and Ensembl vs NCBI may differ in quality.
- For Excon you want a *single*, best-quality proteome per species.
- BUSCO gives a standard completeness metric (C%, F%, M%), so comparing Ensembl vs NCBI at the **protein level** is an objective way to choose.

Conceptual behaviour of `run_busco_both.sh`:

1. You run it from `genomes_out/` so it can see the `both/` subfolder.
2. For each species in `genomes_out/both/`:
   - Runs BUSCO (mode `protein`) on the **Ensembl longest-isoform proteome**.
   - Runs BUSCO (mode `protein`) on the **NCBI longest-isoform proteome**.
3. Organises BUSCO run directories (one per combination of species × source).
4. Produces some summary table (name depends on your implementation) that reports, for each species:
   - BUSCO completeness (%C) for Ensembl and NCBI.
   - Fragmented and missing BUSCOs.
   - A “winner” or at least the information you use to manually choose which source to keep.

You then use this information to decide, per species in `both/`, whether you continue with the Ensembl or NCBI proteome.

### How to run

From inside `genomes_out/`:

```bash
cd genomes_out
bash run_busco_both.sh
cd ..
```

After this step, you know for each “both” species which source you will use later in `3.run_busco_all.sh`.

---

## Step 4 – BUSCO on the final proteome set for all species

### File

- `3.run_busco_all.sh`

### What & why

**Goal:**  
Run BUSCO on the **final chosen proteome** for every species (one proteome per species: Ensembl *or* NCBI), and store all BUSCO runs in a single, organised directory.

Why:

- Excon will later use these BUSCO scores to:
  - Filter out very poor assemblies.
  - Provide genome quality covariates.
- Having all BUSCO runs in a consistent layout makes downstream parsing and QC much easier.

Conceptual behaviour of `3.run_busco_all.sh`:

1. Uses the results/choices from `run_busco_both.sh` (for “both” species) plus:
   - Ensembl-only species from `genomes_out/ensembl/`.
   - NCBI-only species from `genomes_out/ncbi/`.
2. For each species:
   - Identifies the correct “longest-isoform” proteome to use.
   - Runs BUSCO in `protein` mode with your chosen lineage dataset.
   - Writes each run to a subdirectory in `busco_annotation_all/`, with names that encode species + source, e.g.:

```text
busco_annotation_all/
  Species1_ENSEMBL_prot_long/
    run_*/short_summary_*.txt
    run_*/full_table_*.tsv
  Species2_NCBI_prot_long/
    run_*/short_summary_*.txt
    run_*/full_table_*.tsv
  ...
```

3. When this script finishes, **all BUSCO runs** you care about are in `busco_annotation_all/`. This is exactly the directory used by the aggregator in Step 5.

### How to run

From `Download_genomes/`:

```bash
bash 3.run_busco_all.sh
```

---

## Step 5 – Aggregate BUSCO results into a single table

### Files

- `4.busco_annotations_aggregate.py` – Python script that parses BUSCO outputs.
- `4.busco_annotations_aggregate.sh` – shell wrapper that calls the Python script with appropriate arguments.

### What & why

**Goal:**  
Turn all the BUSCO run directories in `busco_annotation_all/` into a **single TSV summary** file that you can use for QC, plotting, filtering, and Excon inputs.

Why:

- BUSCO creates one directory per run, each with its own `short_summary*.txt` and `full_table*.tsv`.
- Manually extracting completeness, missing, fragmented counts, etc., across dozens of species is tedious and error-prone.
- A single table is easy to:
  - Join with your species metadata.
  - Use for filtering and plotting completeness distributions.

What `4.busco_annotations_aggregate.py` does:

1. Takes a BUSCO directory (by default `busco_annotation_all/`) and an output file name (by default `000_BUSCO_annotations_summary.tsv`).
2. Loops over all subdirectories (each BUSCO run).
3. For every run:
   - Infers **species** and **source** (Ensembl vs NCBI) from the folder name (e.g. `_ENSEMBL_prot_long` vs `_NCBI_prot_long`).
   - Finds:
     - `short_summary_*.txt`
     - Optionally `full_table_*.tsv`
   - Parses the `short_summary` to extract:
     - BUSCO version
     - Lineage dataset
     - Mode (protein)
     - Input file name
     - Percentages: C, S, D, F, M
     - Counts: n (total BUSCOs), #Complete, #Single, #Duplicated, #Fragmented, #Missing
4. Stores one row per run and writes a tab-separated file with columns such as:

   - `run_name`
   - `species`
   - `source` (ENSEMBL / NCBI)
   - `busco_version`
   - `lineage`
   - `mode`
   - `input_file`
   - `busco_C_pct`, `busco_S_pct`, `busco_D_pct`, `busco_F_pct`, `busco_M_pct`
   - `busco_n`
   - `n_complete`, `n_single`, `n_duplicated`, `n_fragmented`, `n_missing`
   - `short_summary`
   - `full_table`

   into `000_BUSCO_annotations_summary.tsv`.

### How to run

From `Download_genomes/`:

```bash
bash 4.busco_annotations_aggregate.sh
```

or directly (if you want to override defaults):

```bash
python3 4.busco_annotations_aggregate.py   --busco_dir busco_annotation_all   --out 000_BUSCO_annotations_summary.tsv
```

---

## Quick “one-shot” recap

Typical full pipeline:

```bash
# 1. Download genomes & annotations
bash 1.download_genomes_all.sh

# 2. Build protein sets (longest isoform per gene)
bash 2.build_proteins_longest_all.sh

# 3. Compare Ensembl vs NCBI for species in 'both'
cd genomes_out
bash run_busco_both.sh
cd ..

# (Use the results of step 3 to decide for each "both" species.)

# 4. Run BUSCO on the final proteome set for all species
bash 3.run_busco_all.sh

# 5. Aggregate BUSCO results
bash 4.busco_annotations_aggregate.sh
```

You can now use `download_summary_all.tsv` together with  
`000_BUSCO_annotations_summary.tsv` as the master overview of all genomes and their BUSCO quality for the Excon analyses.
