#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: calculate_assembly_stats.sh --staged-manifest PATH --output PATH

Compute in-house assembly statistics from staged FASTA files using seqtk.
The staged manifest must contain accession and staged_filename columns.
EOF
}

log_error() {
    printf 'ERROR: %s\n' "$1" >&2
}

find_column_index() {
    local -n header_ref=$1
    local target="$2"
    local index=0
    local normalised_target

    normalised_target="$(printf '%s' "$target" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')"
    for column in "${header_ref[@]}"; do
        if [[ "$(printf '%s' "$column" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')" == "${normalised_target}" ]]; then
            printf '%s\n' "${index}"
            return 0
        fi
        index=$((index + 1))
    done

    return 1
}

resolve_genome_path() {
    local manifest_dir="$1"
    local staged_filename="$2"

    if [[ "${staged_filename}" = /* ]]; then
        printf '%s\n' "${staged_filename}"
        return 0
    fi
    printf '%s/%s\n' "${manifest_dir}" "${staged_filename}"
}

compute_assembly_stats() {
    local fasta_path="$1"
    local lengths_file summary_file n50 target scaffolds genome_size

    lengths_file="$(mktemp)"
    summary_file="$(mktemp)"
    CLEANUP_FILES+=("${lengths_file}" "${summary_file}")

    if ! seqtk comp "${fasta_path}" | awk -v summary_file="${summary_file}" '
        BEGIN {
            count = 0
            total = 0
        }
        {
            if ($2 ~ /^[0-9]+$/) {
                print $2
                count += 1
                total += $2
            }
        }
        END {
            if (count == 0) {
                exit 1
            }
            printf "%d\t%d\n", count, total > summary_file
        }
    ' > "${lengths_file}"; then
        return 1
    fi

    IFS=$'\t' read -r scaffolds genome_size < "${summary_file}"
    target=$(( (genome_size + 1) / 2 ))
    n50="$(
        sort -nr "${lengths_file}" | awk -v target="${target}" '
            {
                cumulative += $1
                if (cumulative >= target) {
                    print $1
                    exit
                }
            }
        '
    )"

    if [[ -z "${n50}" ]]; then
        return 1
    fi

    printf '%s\t%s\t%s\n' "${n50}" "${scaffolds}" "${genome_size}"
}

main() {
    local staged_manifest="" output="" manifest_dir header_line duplicate_accessions=""
    local accession_index staged_filename_index required_field_count sort_key malformed_row_info

    while (($# > 0)); do
        case "$1" in
            --staged-manifest)
                staged_manifest="$2"
                shift 2
                ;;
            --output)
                output="$2"
                shift 2
                ;;
            -h|--help)
                usage
                return 0
                ;;
            *)
                log_error "Unknown argument: $1"
                usage >&2
                return 1
                ;;
        esac
    done

    if [[ -z "${staged_manifest}" || -z "${output}" ]]; then
        log_error "--staged-manifest and --output are required."
        usage >&2
        return 1
    fi
    if [[ ! -f "${staged_manifest}" ]]; then
        log_error "Missing staged manifest: ${staged_manifest}"
        return 1
    fi
    if ! command -v seqtk >/dev/null 2>&1; then
        log_error "seqtk is required to calculate in-house assembly statistics."
        return 1
    fi

    manifest_dir="$(cd "$(dirname "${staged_manifest}")" && pwd)"
    IFS=$'\t' read -r -a header_line < "${staged_manifest}"
    if ((${#header_line[@]} == 0)); then
        log_error "Staged manifest is empty: ${staged_manifest}"
        return 1
    fi

    if ! accession_index="$(find_column_index header_line "accession")"; then
        log_error "Staged manifest is missing the accession column."
        return 1
    fi
    if ! staged_filename_index="$(find_column_index header_line "staged_filename")"; then
        log_error "Staged manifest is missing the staged_filename column."
        return 1
    fi
    required_field_count=$((accession_index > staged_filename_index ? accession_index + 1 : staged_filename_index + 1))

    malformed_row_info="$(
        awk -F '\t' -v min_fields="${required_field_count}" '
            NR > 1 && NF < min_fields {
                printf "%d\t%d\n", NR, NF
                exit
            }
        ' "${staged_manifest}"
    )"
    if [[ -n "${malformed_row_info}" ]]; then
        local malformed_row_number observed_fields
        IFS=$'\t' read -r malformed_row_number observed_fields <<< "${malformed_row_info}"
        log_error "Staged manifest row ${malformed_row_number} is malformed: expected at least ${required_field_count} tab-delimited fields, found ${observed_fields}."
        return 1
    fi

    duplicate_accessions="$(
        awk -F '\t' -v accession_index="$((accession_index + 1))" 'NR > 1 { print $accession_index }' "${staged_manifest}" \
            | sed '/^[[:space:]]*$/d' \
            | sort \
            | uniq -d
    )"
    if [[ -n "${duplicate_accessions}" ]]; then
        log_error "Staged manifest contains duplicate accession values: ${duplicate_accessions}"
        return 1
    fi
    sort_key="$((accession_index + 1))"

    mkdir -p "$(dirname "${output}")"
    {
        printf 'accession\tn50\tscaffolds\tgenome_size\n'
        while IFS=$'\t' read -r -a fields; do
            local accession staged_filename fasta_path stats_line

            if ((${#fields[@]} <= staged_filename_index)) || ((${#fields[@]} <= accession_index)); then
                log_error "Staged manifest row is missing required fields."
                return 1
            fi

            accession="${fields[accession_index]}"
            staged_filename="${fields[staged_filename_index]}"
            if [[ -z "${accession}" ]]; then
                log_error "Staged manifest contains an empty accession."
                return 1
            fi
            if [[ -z "${staged_filename}" ]]; then
                log_error "Staged manifest contains an empty staged_filename for ${accession}."
                return 1
            fi

            fasta_path="$(resolve_genome_path "${manifest_dir}" "${staged_filename}")"
            if [[ ! -f "${fasta_path}" ]]; then
                log_error "Missing staged FASTA for ${accession}: ${fasta_path}"
                return 1
            fi

            if ! stats_line="$(compute_assembly_stats "${fasta_path}")"; then
                log_error "Failed to calculate assembly statistics for ${accession}."
                return 1
            fi

            printf '%s\t%s\n' "${accession}" "${stats_line}"
        done < <(
            tail -n +2 "${staged_manifest}" \
                | LC_ALL=C sort -t $'\t' -k"${sort_key},${sort_key}"
        )
    } > "${output}"
}

CLEANUP_FILES=()
trap 'if ((${#CLEANUP_FILES[@]} > 0)); then rm -f "${CLEANUP_FILES[@]}"; fi' EXIT

main "$@"
