#!/usr/bin/env bash
set -Eeuo pipefail

# Defaults are repository-relative. EXCON_PIPELINE may be a local checkout or
# the remote pipeline name Eco-Flow/excon. Set EXCON_REVISION to pin a known
# tag/commit; the archived orthology run report did not record one.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
EXCON_PIPELINE="${EXCON_PIPELINE:-Eco-Flow/excon}"
EXCON_REVISION="${EXCON_REVISION:-}"
INPUT_CSV="${INPUT_CSV:-${SCRIPT_DIR}/input_excon.csv}"
OUTDIR="${OUTDIR:-${SCRIPT_DIR}/results_EXCON}"
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-${SCRIPT_DIR}/eggnog_data}"
PROFILE="${PROFILE:-local,docker}"
NXF_VER="${NXF_VER:-25.10.0}"
MAX_CPUS="${MAX_CPUS:-24}"
MAX_MEMORY="${MAX_MEMORY:-96.GB}"

if [[ ! -f "${INPUT_CSV}" ]]; then
    echo "[ERROR] EXCON input CSV not found: ${INPUT_CSV}" >&2
    echo "Run 1.make_excon_input.py first or set INPUT_CSV." >&2
    exit 1
fi
mkdir -p "${OUTDIR}" "${EGGNOG_DATA_DIR}"

revision_args=()
if [[ -n "${EXCON_REVISION}" ]]; then
    revision_args=(-r "${EXCON_REVISION}")
fi

NXF_VER="${NXF_VER}" nextflow run "${EXCON_PIPELINE}" \
  "${revision_args[@]}" \
  -resume \
  -profile "${PROFILE}" \
  --input "${INPUT_CSV}" \
  --outdir "${OUTDIR}" \
  --run_eggnog \
  --eggnog_data_dir "${EGGNOG_DATA_DIR}" \
  --max_cpus "${MAX_CPUS}" \
  --max_memory "${MAX_MEMORY}"
