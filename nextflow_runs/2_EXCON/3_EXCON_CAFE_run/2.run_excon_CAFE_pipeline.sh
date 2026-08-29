#!/usr/bin/env bash
set -Eeuo pipefail

# Reproduce the archived wasp CAFE run. Resolved_Gene_Trees is generated and
# retained by the preceding EXCON/OrthoFinder run at the default path below.
# Override it with ORTHOFINDER_GENETREE_DIR only if moved. Where needed, provide
# the site-specific cluster configuration with CUSTOM_CONFIG.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
ORTHO_RUN="${SCRIPT_DIR}/../2_EXCON_orthofinder_eggnogmapper_run"

EXCON_PIPELINE="${EXCON_PIPELINE:-Eco-Flow/excon}"
EXCON_REVISION="${EXCON_REVISION:-2.3.2}"
NXF_VER="${NXF_VER:-26.04.6}"
PROFILE="${PROFILE:-singularity}"
CUSTOM_CONFIG="${CUSTOM_CONFIG:-}"

INPUT_TREE="${INPUT_TREE:-${REPO_ROOT}/input_CAFE/tree_dating/Vespidae_with_outgroups_dated_primary.nwk}"
INPUT_ORTHOGROUPS="${INPUT_ORTHOGROUPS:-${SCRIPT_DIR}/results_EXCON/cafe/base/N0.tsv}"
ORTHOFINDER_GENETREE_DIR="${ORTHOFINDER_GENETREE_DIR:-${ORTHO_RUN}/results_EXCON/orthofinder/Resolved_Gene_Trees}"
OUTDIR="${OUTDIR:-${SCRIPT_DIR}/results_EXCON}"
CAFE_FOCUS_CLADES="${CAFE_FOCUS_CLADES:-Vespa_crabro,Vespula_vulgaris|Vespa_crabro,Polistes_dominula|Polistes_dominula,Mischocyttarus_mexicanus}"
CAFE_MAX_DIFFERENTIAL="${CAFE_MAX_DIFFERENTIAL:-20}"

for required in "${INPUT_TREE}" "${INPUT_ORTHOGROUPS}"; do
    if [[ ! -f "${required}" ]]; then
        echo "[ERROR] Required input not found: ${required}" >&2
        exit 1
    fi
done
if [[ ! -d "${ORTHOFINDER_GENETREE_DIR}" ]]; then
    echo "[ERROR] Resolved gene-tree directory not found: ${ORTHOFINDER_GENETREE_DIR}" >&2
    echo "Set ORTHOFINDER_GENETREE_DIR to the EXCON/OrthoFinder output." >&2
    exit 1
fi

config_args=()
if [[ -n "${CUSTOM_CONFIG}" ]]; then
    if [[ ! -f "${CUSTOM_CONFIG}" ]]; then
        echo "[ERROR] CUSTOM_CONFIG does not exist: ${CUSTOM_CONFIG}" >&2
        exit 1
    fi
    config_args=(--custom_config "${CUSTOM_CONFIG}")
fi

mkdir -p "${OUTDIR}"
NXF_VER="${NXF_VER}" nextflow run "${EXCON_PIPELINE}" \
  -r "${EXCON_REVISION}" \
  -profile "${PROFILE}" \
  "${config_args[@]}" \
  --input_tree "${INPUT_TREE}" \
  --input_tree_is_dated \
  --input_orthogroups "${INPUT_ORTHOGROUPS}" \
  --cafe_focus_clades "${CAFE_FOCUS_CLADES}" \
  --orthofinder_genetree_dir "${ORTHOFINDER_GENETREE_DIR}" \
  --outdir "${OUTDIR}" \
  --cafe_max_differential "${CAFE_MAX_DIFFERENTIAL}" \
  --cafe_filter_first \
  -bg \
  -resume
