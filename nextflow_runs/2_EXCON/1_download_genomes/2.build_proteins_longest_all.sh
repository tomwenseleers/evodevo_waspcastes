#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
GENOMES_OUT="${GENOMES_OUT:-${SCRIPT_DIR}/genomes_out}"
GFFREAD_BIN="${GFFREAD_BIN:-gffread}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/2.build_proteins_longest_all.py" \
  --base "${GENOMES_OUT}" \
  --gffread "${GFFREAD_BIN}"

