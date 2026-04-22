#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: backfill_assembly_stats_gc.sh --results-dir PATH [--input PATH] [--output PATH] [--seqtk-binary CMD]

Backfill gc_content into an existing assembly_stats.tsv for one previous run.
The script preserves the existing table rows and columns, appending gc_content
derived from the published staged genomes under samples/<accession>/staged/.
EOF
}

cleanup_output_tmp() {
    if [[ -n "${OUTPUT_TMP}" ]]; then
        rm -f "${OUTPUT_TMP}"
    fi
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

resolve_staged_fasta() {
    local results_dir="$1"
    local accession="$2"
    local staged_dir fasta_candidates
    local -a fasta_paths

    staged_dir="${results_dir}/samples/${accession}/staged"
    if [[ ! -d "${staged_dir}" ]]; then
        log_error "Missing staged directory for ${accession}: ${staged_dir}"
        return 1
    fi

    shopt -s nullglob
    fasta_paths=("${staged_dir}"/*.fasta)
    shopt -u nullglob

    if ((${#fasta_paths[@]} == 0)); then
        log_error "Expected exactly one staged FASTA for ${accession}, found none in ${staged_dir}."
        return 1
    fi
    if ((${#fasta_paths[@]} > 1)); then
        fasta_candidates="$(printf '%s, ' "${fasta_paths[@]}")"
        fasta_candidates="${fasta_candidates%, }"
        log_error "Expected exactly one staged FASTA for ${accession}, found multiple: ${fasta_candidates}"
        return 1
    fi

    printf '%s\n' "${fasta_paths[0]}"
}

run_seqtk_comp() {
    local fasta_path="$1"

    SEQTK_BINARY="${seqtk_binary}" FASTA_PATH="${fasta_path}" \
        bash -lc 'eval "$SEQTK_BINARY" comp "$FASTA_PATH"'
}

compute_gc_content() {
    local fasta_path="$1"
    local gc_content

    gc_content="$(
        run_seqtk_comp "${fasta_path}" | awk '
            BEGIN {
                adenine = 0
                cytosine = 0
                guanine = 0
                thymine = 0
            }
            {
                if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ && $6 ~ /^[0-9]+$/) {
                    adenine += $3
                    cytosine += $4
                    guanine += $5
                    thymine += $6
                }
            }
            END {
                canonical = adenine + cytosine + guanine + thymine
                if (canonical == 0) {
                    exit 1
                }
                value = ((cytosine + guanine) / canonical) * 100
                formatted = sprintf("%.6f", value)
                sub(/0+$/, "", formatted)
                sub(/\.$/, "", formatted)
                print formatted
            }
        '
    )"

    if [[ -z "${gc_content}" ]]; then
        return 1
    fi
    printf '%s\n' "${gc_content}"
}

main() {
    local results_dir="" input="" output="" output_dir
    local accession_index required_column
    local duplicate_accessions="" gc_content accession staged_fasta
    local -a header_array fields

    seqtk_binary="seqtk"

    while (($# > 0)); do
        case "$1" in
            --results-dir)
                results_dir="$2"
                shift 2
                ;;
            --input)
                input="$2"
                shift 2
                ;;
            --output)
                output="$2"
                shift 2
                ;;
            --seqtk-binary)
                seqtk_binary="$2"
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

    if [[ -z "${results_dir}" ]]; then
        log_error "--results-dir is required."
        usage >&2
        return 1
    fi
    if [[ ! -d "${results_dir}" ]]; then
        log_error "Missing results directory: ${results_dir}"
        return 1
    fi

    input="${input:-${results_dir}/cohort/assembly_stats/assembly_stats.tsv}"
    output="${output:-${results_dir}/cohort/assembly_stats/assembly_stats.with_gc.tsv}"

    if [[ ! -f "${input}" ]]; then
        log_error "Missing input assembly stats table: ${input}"
        return 1
    fi

    IFS=$'\t' read -r -a header_array < <(head -n 1 "${input}" | tr -d '\r')
    if ((${#header_array[@]} == 0)); then
        log_error "Input assembly stats table is empty: ${input}"
        return 1
    fi

    for required_column in accession n50 scaffolds genome_size; do
        if ! find_column_index header_array "${required_column}" >/dev/null; then
            log_error "Input assembly stats table is missing the ${required_column} column."
            return 1
        fi
    done
    if find_column_index header_array "gc_content" >/dev/null 2>&1; then
        log_error "Input assembly stats table already contains gc_content: ${input}"
        return 1
    fi

    accession_index="$(find_column_index header_array "accession")"
    duplicate_accessions="$(
        awk -F '\t' -v accession_index="$((accession_index + 1))" 'NR > 1 { print $accession_index }' "${input}" \
            | sed '/^[[:space:]]*$/d' \
            | sort \
            | uniq -d
    )"
    if [[ -n "${duplicate_accessions}" ]]; then
        log_error "Input assembly stats table contains duplicate accession values: ${duplicate_accessions}"
        return 1
    fi

    if ! bash -lc 'command -v "$1" >/dev/null 2>&1' _ "${seqtk_binary%% *}"; then
        if [[ "${seqtk_binary}" != *" "* ]]; then
            log_error "seqtk is required to backfill gc_content."
            return 1
        fi
    fi

    output_dir="$(dirname "${output}")"
    mkdir -p "${output_dir}"
    OUTPUT_TMP="$(mktemp "${output_dir}/.assembly_stats_gc.XXXXXX")"
    trap cleanup_output_tmp EXIT

    {
        printf '%s\tgc_content\n' "$(IFS=$'\t'; printf '%s' "${header_array[*]}")"

        while IFS=$'\t' read -r -a fields; do
            if ((${#fields[@]} == 0)); then
                continue
            fi
            if ((${#fields[@]} <= accession_index)); then
                log_error "Assembly stats row is missing the accession field."
                return 1
            fi

            accession="${fields[accession_index]}"
            if [[ -z "${accession}" ]]; then
                log_error "Assembly stats row contains an empty accession."
                return 1
            fi

            staged_fasta="$(resolve_staged_fasta "${results_dir}" "${accession}")" || return 1
            if ! gc_content="$(compute_gc_content "${staged_fasta}")"; then
                log_error "Failed to calculate gc_content for ${accession} from ${staged_fasta}."
                return 1
            fi

            printf '%s\t%s\n' "$(IFS=$'\t'; printf '%s' "${fields[*]}")" "${gc_content}"
        done < <(tail -n +2 "${input}" | tr -d '\r')
    } > "${OUTPUT_TMP}"

    mv "${OUTPUT_TMP}" "${output}"
    OUTPUT_TMP=""
    trap - EXIT
}

seqtk_binary=""
OUTPUT_TMP=""
main "$@"
