#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Default lineage for insects
LINEAGE_DEFAULT="${LINEAGE_DEFAULT:-insecta_odb10}"

# Special lineage for Caenorhabditis elegans
LINEAGE_CELEGANS="${LINEAGE_CELEGANS:-nematoda_odb10}"

# Where BUSCO stores/reads lineage data
BUSCO_DB="${BUSCO_DB:-${SCRIPT_DIR}/busco_downloads}"

# Where to put BUSCO result folders
OUTDIR="${BUSCO_OUTDIR:-${SCRIPT_DIR}/busco_annotation_all}"
GENOMES_OUT="${GENOMES_OUT:-${SCRIPT_DIR}/genomes_out}"

# CPUs per BUSCO job
NCPU="${NCPU:-8}"

# busco binary
BUSCO_BIN="${BUSCO_BIN:-busco}"
BUSCO_OFFLINE="${BUSCO_OFFLINE:-true}"

# ---------------------------
# FUNCTIONS
# ---------------------------

choose_lineage() {
    local species="$1"
    if [[ "$species" == "Caenorhabditis_elegans" ]]; then
        echo "$LINEAGE_CELEGANS"
    else
        echo "$LINEAGE_DEFAULT"
    fi
}

run_busco() {
    local fa="$1"
    local species="$2"
    local run_name="$3"

    local lineage
    lineage="$(choose_lineage "$species")"

    echo "[BUSCO] $run_name"
    echo "        species : $species"
    echo "        input   : $fa"
    echo "        lineage : $lineage"
    echo

    local offline_args=()
    if [[ "${BUSCO_OFFLINE}" == "true" ]]; then
        offline_args+=(--offline)
    fi

    "${BUSCO_BIN}" \
        -i "$fa" \
        -l "$lineage" \
        -m proteins \
        -o "$run_name" \
        --out_path "$OUTDIR" \
        --download_path "$BUSCO_DB" \
        "${offline_args[@]}" \
        --cpu "$NCPU"
}

# ---------------------------
# MAIN
# ---------------------------

mkdir -p "${OUTDIR}"

echo "[INFO] Running BUSCO (proteins mode) for ALL species"
echo "[INFO] Default lineage : ${LINEAGE_DEFAULT}"
echo "[INFO] C. elegans     : ${LINEAGE_CELEGANS}"
echo "[INFO] Outdir         : ${OUTDIR}"
echo

# Build set of species that have BOTH Ensembl + NCBI in genomes_out/both/
declare -A in_both

shopt -s nullglob
for fa in "${GENOMES_OUT}"/both/00_proteins_out/*_NCBI_proteins_longest.fa; do
    base=$(basename "$fa")
    sp=${base%_NCBI_proteins_longest.fa}
    in_both["$sp"]=1
done
shopt -u nullglob

# -------- Ensembl-only species (ignore those with "both") --------
echo "[INFO] === Ensembl-only species (genomes_out/ensembl) ==="

shopt -s nullglob
ensembl_files=("${GENOMES_OUT}"/ensembl/00_proteins_out/*_ENSEMBL_proteins_longest.fa)
if (( ${#ensembl_files[@]} == 0 )); then
    echo "[WARN] No Ensembl *_ENSEMBL_proteins_longest.fa found."
else
    for fa in "${ensembl_files[@]}"; do
        base=$(basename "$fa")
        sp=${base%_ENSEMBL_proteins_longest.fa}

        # Skip if this species also exists in 'both' (we will use NCBI there)
        if [[ -n "${in_both[$sp]:-}" ]]; then
            echo "[SKIP] ${sp} (also in 'both'; will use NCBI from genomes_out/both)"
            continue
        fi

        run_busco "$fa" "$sp" "${sp}_ENSEMBL_prot_long"
    done
fi
shopt -u nullglob

# -------- NCBI-only species (ignore those with 'both') --------
echo
echo "[INFO] === NCBI-only species (genomes_out/ncbi) ==="

shopt -s nullglob
ncbi_files=("${GENOMES_OUT}"/ncbi/00_proteins_out/*_NCBI_proteins_longest.fa)
if (( ${#ncbi_files[@]} == 0 )); then
    echo "[WARN] No NCBI *_NCBI_proteins_longest.fa found."
else
    for fa in "${ncbi_files[@]}"; do
        base=$(basename "$fa")
        sp=${base%_NCBI_proteins_longest.fa}

        # Skip if this species also exists in 'both' (we will use NCBI from 'both')
        if [[ -n "${in_both[$sp]:-}" ]]; then
            echo "[SKIP] ${sp} (also in 'both'; will use NCBI from genomes_out/both)"
            continue
        fi

        run_busco "$fa" "$sp" "${sp}_NCBI_prot_long"
    done
fi
shopt -u nullglob

# -------- Species with both Ensembl + NCBI (use NCBI only) --------
echo
echo "[INFO] === Species with both Ensembl + NCBI (using NCBI proteins from genomes_out/both) ==="

shopt -s nullglob
both_ncbi_files=("${GENOMES_OUT}"/both/00_proteins_out/*_NCBI_proteins_longest.fa)
if (( ${#both_ncbi_files[@]} == 0 )); then
    echo "[WARN] No *_NCBI_proteins_longest.fa found in genomes_out/both/00_proteins_out."
else
    for fa in "${both_ncbi_files[@]}"; do
        base=$(basename "$fa")
        sp=${base%_NCBI_proteins_longest.fa}

        run_busco "$fa" "$sp" "${sp}_NCBI_prot_long"
    done
fi
shopt -u nullglob

echo
echo "[OK] BUSCO proteins runs for all species finished."
echo "    Results in: ${OUTDIR}/<species>_ENSEMBL_prot_long/ and/or ${OUTDIR}/<species>_NCBI_prot_long/"
