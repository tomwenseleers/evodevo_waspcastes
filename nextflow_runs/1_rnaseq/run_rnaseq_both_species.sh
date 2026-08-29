#!/usr/bin/env bash
# Run the two archived nf-core/rnaseq analyses. Defaults are relative to this
# script. Override paths with environment variables when FASTQs or references
# are stored elsewhere, for example:
#   GENOME_DIR=/data/references INPUT_DIR=/data/samplesheets bash "$0"
# The FASTQ paths inside each samplesheet must also exist on the current host.
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE="${RNASEQ_BASE:-${SCRIPT_DIR}}"
GENOME_DIR="${GENOME_DIR:-${BASE}/genoms}"
INPUT_DIR="${INPUT_DIR:-${BASE}/input_files}"
OUTPUT_DIR="${OUTPUT_DIR:-${BASE}/output_files}"
REPORTS_DIR="${OUTPUT_DIR}/reports"
WORK_DIR="${WORK_DIR:-${BASE}/nf-work}"
TMP_DIR="${TMP_DIR:-${BASE}/tmp}"
PROFILE="${PROFILE:-docker}"
NFCORE_RNASEQ_VERSION="${NFCORE_RNASEQ_VERSION:-3.14.0}"
NXF_VER="${NXF_VER:-24.10.1}"
MAX_CPUS="${MAX_CPUS:-40}"

# 2) Ensure folders exist
mkdir -p "${GENOME_DIR}" "${INPUT_DIR}" "${OUTPUT_DIR}/Results_Vespula" "${OUTPUT_DIR}/Results_Polistes" "${REPORTS_DIR}" "${WORK_DIR}" "${TMP_DIR}"

# 3) Temp env for Nextflow
export TMPDIR="${TMP_DIR}"
export NXF_TEMP="${TMP_DIR}"

# 4) Inputs (samplesheets)
P_SHEET="${INPUT_DIR}/samplesheetp.csv"
V_SHEET="${INPUT_DIR}/samplesheetv.csv"
P_FASTA="${POLISTES_FASTA:-${GENOME_DIR}/GCF_001465965.1_Pdom_r1.2_genomic.fna}"
P_GTF="${POLISTES_GTF:-${GENOME_DIR}/gene_annotation_polistes.gtf}"
V_FASTA="${VESPULA_FASTA:-${GENOME_DIR}/GCF_905475345.1_iyVesVulg1.1_genomic_softmasked.fa}"
V_GTF="${VESPULA_GTF:-${GENOME_DIR}/gene_annotation_vespula.gtf}"

for required in "${P_SHEET}" "${V_SHEET}" "${P_FASTA}" "${P_GTF}" "${V_FASTA}" "${V_GTF}"; do
    if [[ ! -f "${required}" ]]; then
        echo "[ERROR] Required input not found: ${required}" >&2
        exit 1
    fi
done

# 5) Timestamp for reports
STAMP=$(date +%Y%m%d_%H%M%S)

# 6) Run the pipeline(s)

# --- Polistes dominula ---
NXF_VER="${NXF_VER}" nextflow run nf-core/rnaseq -r "${NFCORE_RNASEQ_VERSION}" -profile "${PROFILE}" \
-w "${WORK_DIR}" \
--outdir "${OUTPUT_DIR}/Results_Polistes" \
--input "${P_SHEET}" \
--fasta "${P_FASTA}" \
--gtf "${P_GTF}" \
--featurecounts_group_type transcript_id \
--publish_dir_mode link \
--max_cpus "${MAX_CPUS}" \
-with-report "${REPORTS_DIR}/Polistes_report_${STAMP}.html" \
-with-timeline "${REPORTS_DIR}/Polistes_timeline_${STAMP}.html" \
-with-trace "${REPORTS_DIR}/Polistes_trace_${STAMP}.tsv" \
-with-dag "${REPORTS_DIR}/Polistes_flowchart_${STAMP}.png"

# --- Vespula vulgaris ---
NXF_VER="${NXF_VER}" nextflow run nf-core/rnaseq -r "${NFCORE_RNASEQ_VERSION}" -profile "${PROFILE}" \
-w "${WORK_DIR}" \
--outdir "${OUTPUT_DIR}/Results_Vespula" \
--input "${V_SHEET}" \
--fasta "${V_FASTA}" \
--gtf "${V_GTF}" \
--featurecounts_group_type transcript_id \
--publish_dir_mode link \
--max_cpus "${MAX_CPUS}" \
-with-report "${REPORTS_DIR}/Vespula_report_${STAMP}.html" \
-with-timeline "${REPORTS_DIR}/Vespula_timeline_${STAMP}.html" \
-with-trace "${REPORTS_DIR}/Vespula_trace_${STAMP}.tsv" \
-with-dag "${REPORTS_DIR}/Vespula_flowchart_${STAMP}.png"
