#!/usr/bin/env bash
set -euo pipefail

real_mv=/usr/bin/mv
protected=(tool_bin ccfinder_run ccfinder_raw)

if (( $# <= 1 )); then
    exec "${real_mv}" "$@"
fi

destination="${!#}"
filtered_sources=()

for ((i = 1; i < $#; i++)); do
    arg="${!i}"
    source_name="$(basename "${arg}")"
    keep_source=true

    for protected_name in "${protected[@]}"; do
        if [[ "${source_name}" == "${protected_name}" ]]; then
            keep_source=false
            break
        fi
    done

    if [[ "${keep_source}" == true ]]; then
        filtered_sources+=("${arg}")
    fi
done

if (( ${#filtered_sources[@]} == 0 )); then
    exit 0
fi

exec "${real_mv}" "${filtered_sources[@]}" "${destination}"
