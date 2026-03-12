#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly PIPELINE_WRAPPER="bin/run_pipeline_test.sh"
readonly VALIDATOR="bin/validate_hpc_matrix.py"
readonly MEDIUM_INPUT_BUILDER="bin/build_real_run_inputs.py"
readonly VALID_MODES="prepare db-full db-reuse db-matrix p1 p2 all"

show_usage() {
    cat <<'EOF'
Usage:
  bin/run_oist_hpc_matrix.sh [--dry-run] [--resume] --hpc-root PATH <prepare|db-full|db-reuse|db-matrix|p1|p2|all> [options]

Run the OIST HPC validation campaign through the medium real-data case.

Required:
  --hpc-root PATH             Campaign root on HPC-native storage.

Options:
  --dry-run                   Print commands and filesystem actions without executing them.
  --resume                    Pass -resume to Nextflow and wrapper runs where supported.
  --medium-candidates-tsv PATH
                              Candidate TSV for generating medium-cohort inputs.
  --medium-sample-csv PATH    Medium-cohort sample sheet. Used when candidates TSV is omitted.
  --medium-metadata PATH      Medium-cohort metadata TSV. Used when candidates TSV is omitted.
  --sample-count N            Expected medium-cohort size. Default: 20.
  --singularity-cache PATH    Override the Singularity cache root.
  -h, --help                  Show this help message.

Modes:
  prepare     Prepare the tracked 9-sample cohort on the login node.
  db-full     Run the full runtime-database prep gate on SLURM.
  db-reuse    Re-run the full runtime-database prep gate and expect reuse.
  db-matrix   Run disposable database existence-state cases.
  p1          Run the tracked 9-sample real pipeline gate.
  p2          Run the medium Mycoplasmatota/Bacillota real-data gate.
  all         Run prepare, db-full, db-reuse, db-matrix, p1, and p2 in order.

Notes:
  - Resource settings come from the coded OIST profile and process defaults.
  - PADLOC is not part of runtime database prep because its fixed image bundles the DB.
EOF
}

print_command() {
    local escaped
    local parts=()
    local arg
    for arg in "$@"; do
        printf -v escaped '%q' "${arg}"
        parts+=("${escaped}")
    done
    printf '%s\n' "${parts[*]}"
}

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

is_valid_mode() {
    case "$1" in
        prepare | db-full | db-reuse | db-matrix | p1 | p2 | all)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

run_or_print() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command "$@"
        return 0
    fi
    "$@"
}

run_expect_failure() {
    local expected_text="$1"
    local log_path="$2"
    shift 2

    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command "$@"
        return 0
    fi

    mkdir -p "$(dirname "${log_path}")"
    if "$@" >"${log_path}" 2>&1; then
        fail "Command unexpectedly succeeded: $(print_command "$@")"
    fi
    if ! grep -Fq "${expected_text}" "${log_path}"; then
        fail "Failure log ${log_path} does not contain expected text: ${expected_text}"
    fi
}

ensure_medium_inputs() {
    if [[ -n "${MEDIUM_CANDIDATES_TSV}" ]]; then
        generate_medium_inputs
        return 0
    fi
    if [[ -z "${MEDIUM_SAMPLE_CSV}" || -z "${MEDIUM_METADATA}" ]]; then
        fail "--medium-candidates-tsv or both --medium-sample-csv and --medium-metadata are required for ${MODE}"
    fi
}

generate_medium_inputs() {
    local generated_dir="${HPC_ROOT}/medium_inputs/generated"

    run_or_print \
        python3 "${MEDIUM_INPUT_BUILDER}" \
        --candidate-tsv "${MEDIUM_CANDIDATES_TSV}" \
        --outdir "${generated_dir}" \
        --sample-count "${SAMPLE_COUNT}"

    MEDIUM_SAMPLE_CSV="${generated_dir}/sample_sheet.csv"
    MEDIUM_METADATA="${generated_dir}/metadata.tsv"
}

require_path() {
    local path="$1"
    local description="$2"
    if [[ ! -e "${path}" ]]; then
        fail "Missing ${description}: ${path}"
    fi
}

require_golden_db_tree() {
    require_path "${TAXDUMP_DIR}/.nf_myco_ready.json" "taxdump ready marker"
    require_path "${CHECKM2_DIR}/.nf_myco_ready.json" "CheckM2 ready marker"
    require_path "${BUSCO_DIR}/.nf_myco_ready.json" "BUSCO ready marker"
    require_path "${EGGNOG_DIR}/.nf_myco_ready.json" "eggNOG ready marker"
}

require_tracked_cohort() {
    require_path "${ACCEPT_ROOT}/generated/sample_sheet.csv" "tracked cohort sample sheet"
    require_path "${ACCEPT_ROOT}/generated/metadata.tsv" "tracked cohort metadata"
}

run_prepare() {
    run_or_print "${PIPELINE_WRAPPER}" prepare --work-root "${ACCEPT_ROOT}"
}

run_dbprep_wrapper() {
    local expected_status="$1"
    local resume_args=()

    if [[ "${RESUME}" == "true" ]]; then
        resume_args+=(--resume)
    fi

    run_or_print \
        "${PIPELINE_WRAPPER}" dbprep-slurm \
        "${resume_args[@]}" \
        --dbprep-profile oist \
        --work-root "${ACCEPT_ROOT}" \
        --taxdump "${TAXDUMP_DIR}" \
        --checkm2-db "${CHECKM2_DIR}" \
        --busco-db "${BUSCO_DIR}" \
        --eggnog-db "${EGGNOG_DIR}" \
        --singularity-cache-dir "${SINGULARITY_CACHE}"

    run_or_print \
        python3 "${VALIDATOR}" dbprep \
        --results-dir "${ACCEPT_ROOT}/runs/dbprep-slurm/results" \
        --taxdump "${TAXDUMP_DIR}" \
        --checkm2-db "${CHECKM2_DIR}" \
        --busco-db "${BUSCO_DIR}" \
        --eggnog-db "${EGGNOG_DIR}" \
        --expected-component taxdump \
        --expected-component checkm2 \
        --expected-component busco_root \
        --expected-component eggnog \
        --expected-status "${expected_status}" \
        --expected-arg=--taxdump \
        --expected-arg=--checkm2_db \
        --expected-arg=--busco_db \
        --expected-arg=--eggnog_db
}

run_partial_prepare_case() {
    local work_dir="$1"
    local outdir="$2"
    shift 2
    local resume_args=()
    if [[ "${RESUME}" == "true" ]]; then
        resume_args+=(-resume)
    fi
    run_or_print \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${work_dir}" \
        "${resume_args[@]}" \
        "$@" \
        --outdir "${outdir}" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
}

seed_valid_without_marker_cases() {
    local root="$1"

    run_or_print mkdir -p "${root}/checkm2" "${root}/busco/bacillota_odb12" \
        "${root}/busco/mycoplasmatota_odb12" "${root}/eggnog"
    run_or_print cp "${CHECKM2_DIR}/"*.dmnd "${root}/checkm2/"
    run_or_print cp "${BUSCO_DIR}/bacillota_odb12/dataset.cfg" "${root}/busco/bacillota_odb12/"
    run_or_print cp "${BUSCO_DIR}/mycoplasmatota_odb12/dataset.cfg" "${root}/busco/mycoplasmatota_odb12/"
    run_or_print cp "${EGGNOG_DIR}/eggnog.db" "${root}/eggnog/"
    run_or_print cp "${EGGNOG_DIR}/eggnog_proteins.dmnd" "${root}/eggnog/"
}

run_db_matrix() {
    local case_root="${CASE_ROOT}"
    local fail_log

    if [[ "${DRY_RUN}" != "true" ]]; then
        require_golden_db_tree
    fi

    run_or_print mkdir -p "${case_root}"

    local valid_root="${case_root}/db3_valid_without_marker"
    seed_valid_without_marker_cases "${valid_root}"
    run_partial_prepare_case \
        "${valid_root}/work_checkm2" \
        "${valid_root}/results_checkm2" \
        --checkm2_db "${valid_root}/checkm2" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${valid_root}/results_checkm2" \
        --checkm2-db "${valid_root}/checkm2" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db
    run_partial_prepare_case \
        "${valid_root}/work_busco" \
        "${valid_root}/results_busco" \
        --busco_db "${valid_root}/busco" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${valid_root}/results_busco" \
        --busco-db "${valid_root}/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db
    run_partial_prepare_case \
        "${valid_root}/work_eggnog" \
        "${valid_root}/results_eggnog" \
        --eggnog_db "${valid_root}/eggnog" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${valid_root}/results_eggnog" \
        --eggnog-db "${valid_root}/eggnog" \
        --expected-component eggnog \
        --expected-status prepared \
        --expected-arg=--eggnog_db

    run_partial_prepare_case \
        "${case_root}/db4_taxdump/work" \
        "${case_root}/db4_taxdump/results" \
        --taxdump "${case_root}/db4_taxdump/db/ncbi_taxdump_20240914" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_taxdump/results" \
        --taxdump "${case_root}/db4_taxdump/db/ncbi_taxdump_20240914" \
        --expected-component taxdump \
        --expected-status prepared \
        --expected-arg=--taxdump
    run_partial_prepare_case \
        "${case_root}/db4_checkm2/work" \
        "${case_root}/db4_checkm2/results" \
        --checkm2_db "${case_root}/db4_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_checkm2/results" \
        --checkm2-db "${case_root}/db4_checkm2/db/checkm2/CheckM2_database" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db
    run_partial_prepare_case \
        "${case_root}/db4_busco/work" \
        "${case_root}/db4_busco/results" \
        --busco_db "${case_root}/db4_busco/db/busco" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_busco/results" \
        --busco-db "${case_root}/db4_busco/db/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db
    run_partial_prepare_case \
        "${case_root}/db4_eggnog/work" \
        "${case_root}/db4_eggnog/results" \
        --eggnog_db "${case_root}/db4_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_eggnog/results" \
        --eggnog-db "${case_root}/db4_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --expected-component eggnog \
        --expected-status prepared \
        --expected-arg=--eggnog_db

    fail_log="${case_root}/db5_taxdump_missing_no_download.log"
    run_expect_failure \
        "Taxdump destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_taxdump/work" \
        --taxdump "${case_root}/db5_taxdump/db/ncbi_taxdump_20240914" \
        --outdir "${case_root}/db5_taxdump/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db5_checkm2_missing_no_download.log"
    run_expect_failure \
        "CheckM2 destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_checkm2/work" \
        --checkm2_db "${case_root}/db5_checkm2/db/checkm2/CheckM2_database" \
        --outdir "${case_root}/db5_checkm2/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db5_busco_missing_no_download.log"
    run_expect_failure \
        "BUSCO destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_busco/work" \
        --busco_db "${case_root}/db5_busco/db/busco" \
        --outdir "${case_root}/db5_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db5_eggnog_missing_no_download.log"
    run_expect_failure \
        "eggNOG destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_eggnog/work" \
        --eggnog_db "${case_root}/db5_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --outdir "${case_root}/db5_eggnog/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"

    run_or_print mkdir -p "${case_root}/db6_taxdump/db" "${case_root}/db6_checkm2/db/checkm2" \
        "${case_root}/db6_busco/db" "${case_root}/db6_eggnog/db/Eggnog_db"
    run_or_print touch "${case_root}/db6_taxdump/db/ncbi_taxdump_20240914"
    run_or_print touch "${case_root}/db6_checkm2/db/checkm2/CheckM2_database"
    run_or_print touch "${case_root}/db6_busco/db/busco"
    run_or_print touch "${case_root}/db6_eggnog/db/Eggnog_db/Eggnog_Diamond_db"
    fail_log="${case_root}/db6_taxdump_file.log"
    run_expect_failure \
        "Taxdump destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_taxdump/work" \
        --taxdump "${case_root}/db6_taxdump/db/ncbi_taxdump_20240914" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_taxdump/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db6_checkm2_file.log"
    run_expect_failure \
        "CheckM2 destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_checkm2/work" \
        --checkm2_db "${case_root}/db6_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_checkm2/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db6_busco_file.log"
    run_expect_failure \
        "BUSCO destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_busco/work" \
        --busco_db "${case_root}/db6_busco/db/busco" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db6_eggnog_file.log"
    run_expect_failure \
        "eggNOG destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_eggnog/work" \
        --eggnog_db "${case_root}/db6_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_eggnog/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"

    run_or_print mkdir -p "${case_root}/db7_checkm2/db/checkm2/CheckM2_database"
    run_or_print mkdir -p "${case_root}/db7_busco/db/busco/bacillota_odb12"
    run_or_print mkdir -p "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db"
    run_or_print touch "${case_root}/db7_checkm2/db/checkm2/CheckM2_database/broken.txt"
    run_or_print touch "${case_root}/db7_busco/db/busco/bacillota_odb12/dataset.cfg"
    run_or_print touch "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db/broken.txt"
    fail_log="${case_root}/db7_checkm2_invalid.log"
    run_expect_failure \
        "CheckM2 destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_checkm2/work" \
        --checkm2_db "${case_root}/db7_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_checkm2/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db7_busco_invalid.log"
    run_expect_failure \
        "BUSCO destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_busco/work" \
        --busco_db "${case_root}/db7_busco/db/busco" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    fail_log="${case_root}/db7_eggnog_invalid.log"
    run_expect_failure \
        "eggNOG destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_eggnog/work" \
        --eggnog_db "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_eggnog/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"

    run_partial_prepare_case \
        "${case_root}/db7_checkm2/work_force" \
        "${case_root}/db7_checkm2/results_force" \
        --checkm2_db "${case_root}/db7_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_checkm2/results_force" \
        --checkm2-db "${case_root}/db7_checkm2/db/checkm2/CheckM2_database" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db
    run_partial_prepare_case \
        "${case_root}/db7_busco/work_force" \
        "${case_root}/db7_busco/results_force" \
        --busco_db "${case_root}/db7_busco/db/busco" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_busco/results_force" \
        --busco-db "${case_root}/db7_busco/db/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db
    run_partial_prepare_case \
        "${case_root}/db7_eggnog/work_force" \
        "${case_root}/db7_eggnog/results_force" \
        --eggnog_db "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_eggnog/results_force" \
        --eggnog-db "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --expected-component eggnog \
        --expected-status prepared \
        --expected-arg=--eggnog_db

    run_or_print mkdir -p "${case_root}/db9_busco/db/busco/bacillota_odb12"
    run_or_print touch "${case_root}/db9_busco/db/busco/bacillota_odb12/dataset.cfg"
    fail_log="${case_root}/db9_busco_missing_lineage.log"
    run_expect_failure \
        "BUSCO destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db9_busco/work" \
        --busco_db "${case_root}/db9_busco/db/busco" \
        --download_missing_databases true \
        --outdir "${case_root}/db9_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    run_partial_prepare_case \
        "${case_root}/db9_busco/work_force" \
        "${case_root}/db9_busco/results_force" \
        --busco_db "${case_root}/db9_busco/db/busco" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db9_busco/results_force" \
        --busco-db "${case_root}/db9_busco/db/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db

    run_or_print mkdir -p "${case_root}/db10_checkm2/db/checkm2/CheckM2_database/CheckM2_database"
    run_or_print touch "${case_root}/db10_checkm2/db/checkm2/CheckM2_database/CheckM2_database/CheckM2_database.dmnd"
    run_partial_prepare_case \
        "${case_root}/db10_checkm2/work" \
        "${case_root}/db10_checkm2/results" \
        --checkm2_db "${case_root}/db10_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db10_checkm2/results" \
        --checkm2-db "${case_root}/db10_checkm2/db/checkm2/CheckM2_database" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db

    fail_log="${case_root}/db11_readonly_checkout.log"
    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command test ! -e "${REPO_ROOT}/download_eggnog_data.py.patched"
    else
        if [[ -e "${REPO_ROOT}/download_eggnog_data.py.patched" ]]; then
            fail "Unexpected eggNOG patched downloader in project checkout"
        fi
    fi
}

run_real_case() {
    local work_dir="$1"
    local outdir="$2"
    local sample_csv="$3"
    local metadata_tsv="$4"
    local resume_args=()

    if [[ "${RESUME}" == "true" ]]; then
        resume_args+=(-resume)
    fi

    run_or_print \
        nextflow run . -profile oist \
        -work-dir "${work_dir}" \
        "${resume_args[@]}" \
        --sample_csv "${sample_csv}" \
        --metadata "${metadata_tsv}" \
        --taxdump "${TAXDUMP_DIR}" \
        --checkm2_db "${CHECKM2_DIR}" \
        --busco_db "${BUSCO_DIR}" \
        --eggnog_db "${EGGNOG_DIR}" \
        --singularity_cache_dir "${SINGULARITY_CACHE}" \
        --outdir "${outdir}"
}

run_p1() {
    if [[ "${DRY_RUN}" != "true" ]]; then
        require_tracked_cohort
        require_golden_db_tree
    fi
    run_real_case \
        "${RESULT_ROOT}/p1/work" \
        "${RESULT_ROOT}/p1/out" \
        "${ACCEPT_ROOT}/generated/sample_sheet.csv" \
        "${ACCEPT_ROOT}/generated/metadata.tsv"
    run_or_print \
        python3 "${VALIDATOR}" tracked-run \
        --outdir "${RESULT_ROOT}/p1/out" \
        --metadata-tsv "${ACCEPT_ROOT}/generated/metadata.tsv" \
        --db-root "${DB_ROOT}"
}

run_p2() {
    ensure_medium_inputs
    if [[ "${DRY_RUN}" != "true" ]]; then
        require_golden_db_tree
    fi
    run_real_case \
        "${RESULT_ROOT}/p2/work" \
        "${RESULT_ROOT}/p2/out" \
        "${MEDIUM_SAMPLE_CSV}" \
        "${MEDIUM_METADATA}"
    run_or_print \
        python3 "${VALIDATOR}" medium-run \
        --outdir "${RESULT_ROOT}/p2/out" \
        --db-root "${DB_ROOT}" \
        --sample-count "${SAMPLE_COUNT}" \
        --allowed-phylum Mycoplasmatota \
        --allowed-phylum Bacillota
}

main() {
    local mode=""

    DRY_RUN=false
    RESUME=false
    HPC_ROOT=""
    MEDIUM_SAMPLE_CSV=""
    MEDIUM_METADATA=""
    MEDIUM_CANDIDATES_TSV=""
    SAMPLE_COUNT=20
    SINGULARITY_CACHE=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --resume)
                RESUME=true
                shift
                ;;
            --hpc-root)
                HPC_ROOT="$2"
                shift 2
                ;;
            --medium-sample-csv)
                MEDIUM_SAMPLE_CSV="$2"
                shift 2
                ;;
            --medium-candidates-tsv)
                MEDIUM_CANDIDATES_TSV="$2"
                shift 2
                ;;
            --medium-metadata)
                MEDIUM_METADATA="$2"
                shift 2
                ;;
            --sample-count)
                SAMPLE_COUNT="$2"
                shift 2
                ;;
            --singularity-cache)
                SINGULARITY_CACHE="$2"
                shift 2
                ;;
            -h|--help)
                show_usage
                return 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                fail "Unknown option: $1"
                ;;
            *)
                if [[ -n "${mode}" ]]; then
                    fail "Multiple modes provided: ${mode} and $1"
                fi
                mode="$1"
                shift
                ;;
        esac
    done

    if [[ -z "${mode}" ]]; then
        show_usage >&2
        return 1
    fi
    if ! is_valid_mode "${mode}"; then
        fail "Unsupported mode ${mode}. Expected one of: ${VALID_MODES}"
    fi
    if [[ -z "${HPC_ROOT}" ]]; then
        fail "--hpc-root is required"
    fi

    MODE="${mode}"
    ACCEPT_ROOT="${HPC_ROOT}/acceptance"
    DB_ROOT="${HPC_ROOT}/db"
    RESULT_ROOT="${HPC_ROOT}/results"
    CASE_ROOT="${HPC_ROOT}/db_cases"
    if [[ -z "${SINGULARITY_CACHE}" ]]; then
        SINGULARITY_CACHE="${HPC_ROOT}/singularity_cache"
    fi
    TAXDUMP_DIR="${DB_ROOT}/ncbi_taxdump_20240914"
    CHECKM2_DIR="${DB_ROOT}/checkm2/CheckM2_database"
    BUSCO_DIR="${DB_ROOT}/busco"
    EGGNOG_DIR="${DB_ROOT}/Eggnog_db/Eggnog_Diamond_db"

    if [[ "${DRY_RUN}" != "true" ]]; then
        mkdir -p "${ACCEPT_ROOT}" "${DB_ROOT}" "${RESULT_ROOT}" "${CASE_ROOT}" "${SINGULARITY_CACHE}"
    fi

    case "${MODE}" in
        prepare)
            run_prepare
            ;;
        db-full)
            run_dbprep_wrapper prepared
            ;;
        db-reuse)
            run_dbprep_wrapper present
            ;;
        db-matrix)
            run_db_matrix
            ;;
        p1)
            run_p1
            ;;
        p2)
            run_p2
            ;;
        all)
            run_prepare
            run_dbprep_wrapper prepared
            run_dbprep_wrapper present
            run_db_matrix
            run_p1
            run_p2
            ;;
    esac
}

main "$@"
