#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
BUSCO_OUTDIR="${BUSCO_OUTDIR:-${SCRIPT_DIR}/busco_annotation_all}"
SUMMARY_OUT="${SUMMARY_OUT:-${SCRIPT_DIR}/000_BUSCO_annotations_summary.tsv}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/4.busco_annotations_aggregate.py" \
  --busco_dir "${BUSCO_OUTDIR}" \
  --out "${SUMMARY_OUT}"
