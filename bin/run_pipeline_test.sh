#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly RUNNER="bin/run_acceptance_tests.py"
readonly VALID_MODES="prepare unit stub local slurm all"

show_usage() {
    cat <<'EOF'
Usage:
  bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|all> [args...]

Manual wrapper around bin/run_acceptance_tests.py.

Wrapper options:
  --dry-run     Print the delegated Python command and exit.
  -h, --help    Show this help message.

Notes:
  - Arguments after the mode are forwarded unchanged to bin/run_acceptance_tests.py.
  - Use "<mode> --help" to see the delegated harness help for that mode.
  - Real-data modes still use params.ccfinder_container from pipeline config.
EOF
}

is_valid_mode() {
    local mode="$1"
    case "${mode}" in
        prepare | unit | stub | local | slurm | all)
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
