#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
GENOMES_OUT="${GENOMES_OUT:-${SCRIPT_DIR}/genomes_out}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/1.download_genomes_all.py" \
  --ensembl_csv "${SCRIPT_DIR}/Genomes_selection_ensembl.csv" \
  --ncbi_csv "${SCRIPT_DIR}/Genomes_selection_ncbi.csv" \
  --both_csv "${SCRIPT_DIR}/Genomes_selection_both.csv" \
  --out "${GENOMES_OUT}"
