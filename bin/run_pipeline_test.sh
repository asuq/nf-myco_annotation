#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly RUNNER="bin/run_acceptance_tests.py"
readonly VALID_MODES="prepare unit stub local slurm dbprep-slurm all"
readonly BASE_CONFIG="${REPO_ROOT}/conf/base.config"
readonly CURRENT_RUNTIME_DB_HELPER_IMAGE="quay.io/asuq1617/nf-myco_db:0.3"
readonly STALE_RUNTIME_DB_HELPER_IMAGE="quay.io/asuq1617/nf-myco_db:0.2"
readonly MINIMUM_HOST_PYTHON="3.12"

show_usage() {
    cat <<'EOF'
Usage:
  bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|dbprep-slurm|all> [args...]

Manual wrapper around bin/run_acceptance_tests.py.

Wrapper options:
  --dry-run     Print the delegated Python command and exit.
  -h, --help    Show this help message.

Notes:
  - Arguments after the mode are forwarded unchanged to bin/run_acceptance_tests.py.
  - Use "<mode> --help" to see the delegated harness help for that mode.
  - Host python3 must be >= 3.12 for the acceptance harness and validators.
  - Real-data modes still use params.ccfinder_container from pipeline config.
  - dbprep-slurm and all depend on the current runtime-db helper image and Codetta-aware helper CLIs.
EOF
}

is_valid_mode() {
    local mode="$1"
    case "${mode}" in
        prepare | unit | stub | local | slurm | dbprep-slurm | all)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
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

preflight_runtime_db_helper_image() {
    if grep -Fq "${STALE_RUNTIME_DB_HELPER_IMAGE}" "${BASE_CONFIG}"; then
        printf '%s\n' \
            "Error: ${BASE_CONFIG} still pins stale runtime-db helper image ${STALE_RUNTIME_DB_HELPER_IMAGE}. Update it to ${CURRENT_RUNTIME_DB_HELPER_IMAGE} before running dbprep-slurm or all." \
            >&2
        return 1
    fi
}

preflight_host_python() {
    local version_text

    if ! command -v python3 >/dev/null 2>&1; then
        printf '%s\n' \
            "Error: python3 was not found on PATH. Load Python ${MINIMUM_HOST_PYTHON} or adjust PATH before running ${RUNNER}." \
            >&2
        return 1
    fi

    version_text="$(python3 --version 2>&1 | tr -d '\n')"
    if ! python3 -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 12) else 1)' >/dev/null 2>&1; then
        printf '%s\n' \
            "Error: host python3 must be >= ${MINIMUM_HOST_PYTHON} for ${RUNNER} and the acceptance validators. Current interpreter: ${version_text}. Load Python ${MINIMUM_HOST_PYTHON} or adjust PATH before rerunning." \
            >&2
        return 1
    fi
}

main() {
    local dry_run=false
    local mode
    local command=()

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --dry-run)
                dry_run=true
                shift
                ;;
            -h | --help)
                show_usage
                return 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                printf 'Error: unknown wrapper option %s\n' "$1" >&2
                show_usage >&2
                return 1
                ;;
            *)
                break
                ;;
        esac
    done

    if [[ $# -eq 0 ]]; then
        show_usage >&2
        return 1
    fi

    mode="$1"
    shift

    if ! is_valid_mode "${mode}"; then
        printf 'Error: unsupported mode %s. Expected one of: %s\n' "${mode}" "${VALID_MODES}" >&2
        return 1
    fi

    if [[ "${dry_run}" != "true" ]]; then
        preflight_host_python
    fi

    case "${mode}" in
        dbprep-slurm | all)
            preflight_runtime_db_helper_image
            ;;
    esac

    command=("python3" "${RUNNER}" "${mode}" "$@")

    if [[ "${dry_run}" == "true" ]]; then
        cd "${REPO_ROOT}"
        print_command "${command[@]}"
        return 0
    fi

    cd "${REPO_ROOT}"
    exec "${command[@]}"
}

main "$@"
