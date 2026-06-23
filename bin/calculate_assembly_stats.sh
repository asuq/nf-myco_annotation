#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: calculate_assembly_stats.sh --staged-manifest PATH --output PATH [--jobs INT]

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

validate_jobs() {
    local jobs="$1"

    if ! [[ "${jobs}" =~ ^[1-9][0-9]*$ ]]; then
        log_error "--jobs must be a positive integer."
        return 1
    fi
}

resolve_script_path() {
    local invocation="$1"

    if [[ "${invocation}" == */* ]]; then
        printf '%s\n' "${invocation}"
        return 0
    fi

    command -v "${invocation}"
}

write_worker_failure() {
    local error_path="$1"
    local accession="$2"
    local fasta_path="$3"
    local reason="$4"
    local seqtk_stderr="${5:-}"
    local diagnostic_path="${6:-}"
    local first_nonempty file_size

    {
        printf 'Failed to calculate assembly statistics for %s.\n' "${accession}"
        printf 'FASTA: %s\n' "${fasta_path}"
        printf 'Reason: %s\n' "${reason}"
        if [[ -e "${fasta_path}" ]]; then
            file_size="$(wc -c < "${fasta_path}" 2>/dev/null || true)"
            if [[ -n "${file_size}" ]]; then
                printf 'FASTA size bytes: %s\n' "${file_size//[[:space:]]/}"
            fi
            first_nonempty="$(awk 'NF { print; exit }' "${fasta_path}" 2>/dev/null || true)"
            if [[ -z "${first_nonempty}" ]]; then
                printf 'First non-empty FASTA line: <none>\n'
            elif [[ "${first_nonempty}" != ">"* ]]; then
                printf 'First non-empty FASTA line: %s\n' "${first_nonempty:0:160}"
            fi
        fi
        if [[ -n "${diagnostic_path}" && -s "${diagnostic_path}" ]]; then
            printf 'Diagnostic: '
            cat "${diagnostic_path}"
        fi
        if [[ -n "${seqtk_stderr}" && -s "${seqtk_stderr}" ]]; then
            printf 'seqtk stderr:\n'
            cat "${seqtk_stderr}"
        fi
    } > "${error_path}"
}

compute_assembly_stats() {
    local accession="$1"
    local fasta_path="$2"
    local temp_dir="$3"
    local result_path="$4"
    local error_path="$5"
    local lengths_file summary_file seqtk_stderr diagnostic_file n50 target scaffolds genome_size
    local adenine_count cytosine_count guanine_count thymine_count canonical_bases gc_content
    local pipeline_statuses reason

    lengths_file="$(mktemp "${temp_dir}/lengths.XXXXXX")"
    summary_file="$(mktemp "${temp_dir}/summary.XXXXXX")"
    seqtk_stderr="$(mktemp "${temp_dir}/seqtk_stderr.XXXXXX")"
    diagnostic_file="$(mktemp "${temp_dir}/diagnostic.XXXXXX")"

    if [[ ! -s "${fasta_path}" ]]; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "FASTA is empty or unavailable."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    set +e
    seqtk comp "${fasta_path}" 2> "${seqtk_stderr}" | awk \
        -v summary_file="${summary_file}" \
        -v diagnostic_file="${diagnostic_file}" '
        BEGIN {
            count = 0
            total = 0
            adenine = 0
            cytosine = 0
            guanine = 0
            thymine = 0
        }
        {
            if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ && $6 ~ /^[0-9]+$/) {
                print $2
                count += 1
                total += $2
                adenine += $3
                cytosine += $4
                guanine += $5
                thymine += $6
            }
        }
        END {
            if (count == 0) {
                print "seqtk comp produced no valid contig rows." > diagnostic_file
                exit 1
            }
            if ((adenine + cytosine + guanine + thymine) == 0) {
                print "FASTA contains no canonical A/C/G/T bases." > diagnostic_file
                exit 1
            }
            printf "%d\t%d\t%d\t%d\t%d\t%d\n", \
                count, total, adenine, cytosine, guanine, thymine > summary_file
        }
    ' > "${lengths_file}"
    pipeline_statuses=("${PIPESTATUS[@]}")
    set -e

    if ((pipeline_statuses[0] != 0)); then
        reason="seqtk comp exited with status ${pipeline_statuses[0]}."
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "${reason}" \
            "${seqtk_stderr}" \
            "${diagnostic_file}"
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi
    if ((pipeline_statuses[1] != 0)); then
        reason="seqtk comp output did not contain usable assembly statistics."
        if [[ -s "${diagnostic_file}" ]]; then
            reason="$(cat "${diagnostic_file}")"
        fi
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "${reason}" \
            "${seqtk_stderr}" \
            "${diagnostic_file}"
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    if ! IFS=$'\t' read -r scaffolds genome_size adenine_count cytosine_count guanine_count thymine_count < "${summary_file}"; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "Could not read seqtk summary output."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi
    target=$(( (genome_size + 1) / 2 ))
    if ! n50="$(
        sort -nr "${lengths_file}" | awk -v target="${target}" '
            {
                cumulative += $1
                if (cumulative >= target) {
                    print $1
                    exit
                }
            }
        '
    )"; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "Could not calculate N50 from contig lengths."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    if [[ -z "${n50}" ]]; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "N50 calculation produced an empty value."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    canonical_bases=$((adenine_count + cytosine_count + guanine_count + thymine_count))
    if ((canonical_bases == 0)); then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "FASTA contains no canonical A/C/G/T bases."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    if ! gc_content="$(
        awk -v cytosine="${cytosine_count}" -v guanine="${guanine_count}" -v canonical="${canonical_bases}" '
            BEGIN {
                value = ((cytosine + guanine) / canonical) * 100
                formatted = sprintf("%.6f", value)
                sub(/0+$/, "", formatted)
                sub(/\.$/, "", formatted)
                print formatted
            }
        '
    )"; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "Could not calculate GC content."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi
    if [[ -z "${gc_content}" ]]; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "GC-content calculation produced an empty value."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    if ! printf '%s\t%s\t%s\t%s\t%s\n' \
        "${accession}" \
        "${n50}" \
        "${scaffolds}" \
        "${genome_size}" \
        "${gc_content}" \
        > "${result_path}"; then
        write_worker_failure \
            "${error_path}" \
            "${accession}" \
            "${fasta_path}" \
            "Could not write worker result row."
        rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
        return 1
    fi

    rm -f "${lengths_file}" "${summary_file}" "${seqtk_stderr}" "${diagnostic_file}"
}

run_worker() {
    local row_index="$1"
    local accession="$2"
    local fasta_path="$3"
    local worker_dir="$4"
    local result_path error_path temp_dir

    result_path="${worker_dir}/results/${row_index}.tsv"
    error_path="${worker_dir}/errors/${row_index}.log"
    temp_dir="${worker_dir}/tmp/${row_index}"
    rm -rf "${temp_dir}"
    mkdir -p "${temp_dir}"
    rm -f "${result_path}" "${error_path}"

    if ! compute_assembly_stats "${accession}" "${fasta_path}" "${temp_dir}" "${result_path}" "${error_path}"; then
        rm -rf "${temp_dir}"
        return 1
    fi

    rm -rf "${temp_dir}"
}

collect_worker_error() {
    local worker_dir="$1"
    local task_count="$2"
    local row_index error_path

    for ((row_index = 1; row_index <= task_count; row_index++)); do
        error_path="${worker_dir}/errors/${row_index}.log"
        if [[ -s "${error_path}" ]]; then
            cat "${error_path}"
            return 0
        fi
    done

    return 1
}

merge_worker_results() {
    local worker_dir="$1"
    local task_count="$2"
    local output_path="$3"
    local output_tmp row_index result_path

    output_tmp="${worker_dir}/assembly_stats.tsv"
    {
        printf 'accession\tn50\tscaffolds\tgenome_size\tgc_content\n'
        for ((row_index = 1; row_index <= task_count; row_index++)); do
            result_path="${worker_dir}/results/${row_index}.tsv"
            if [[ ! -f "${result_path}" ]]; then
                log_error "Assembly-stat worker did not produce row ${row_index}."
                return 1
            fi
            cat "${result_path}"
        done
    } > "${output_tmp}"

    mv "${output_tmp}" "${output_path}"
}

main() {
    local staged_manifest="" output="" jobs="1" manifest_dir header_line duplicate_accessions=""
    local accession_index staged_filename_index

    if [[ "${1:-}" == "--worker" ]]; then
        shift
        if (($# != 4)); then
            log_error "Worker mode expects row index, accession, FASTA path, and worker dir."
            return 1
        fi
        run_worker "$1" "$2" "$3" "$4"
        return $?
    fi

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
            --jobs)
                jobs="$2"
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
    if ! validate_jobs "${jobs}"; then
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

    local output_dir worker_dir task_file worker_script task_count
    output_dir="$(dirname "${output}")"
    mkdir -p "${output_dir}"
    output_dir="$(cd "${output_dir}" && pwd)"
    worker_dir="$(mktemp -d "${output_dir}/.calculate_assembly_stats.XXXXXX")"
    CLEANUP_DIRS+=("${worker_dir}")
    mkdir -p "${worker_dir}/results" "${worker_dir}/errors" "${worker_dir}/tmp"
    task_file="${worker_dir}/tasks.nul"
    worker_script="$(resolve_script_path "$0")"
    task_count=0

    while IFS=$'\t' read -r -a fields; do
        local accession staged_filename fasta_path

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

        task_count=$((task_count + 1))
        printf '%s\0%s\0%s\0%s\0' \
            "${task_count}" \
            "${accession}" \
            "${fasta_path}" \
            "${worker_dir}" \
            >> "${task_file}"
    done < <(tail -n +2 "${staged_manifest}" | tr -d '\r')

    if [[ -s "${task_file}" ]]; then
        if ! xargs -0 -n 4 -P "${jobs}" bash "${worker_script}" --worker < "${task_file}"; then
            local worker_error
            worker_error="$(collect_worker_error "${worker_dir}" "${task_count}" || true)"
            if [[ -n "${worker_error}" ]]; then
                log_error "${worker_error}"
            else
                log_error "Assembly-stat worker execution failed."
            fi
            return 1
        fi
    fi

    if ! merge_worker_results "${worker_dir}" "${task_count}" "${output}"; then
        return 1
    fi
}

CLEANUP_FILES=()
CLEANUP_DIRS=()
trap 'if ((${#CLEANUP_FILES[@]} > 0)); then rm -f "${CLEANUP_FILES[@]}"; fi; if ((${#CLEANUP_DIRS[@]} > 0)); then rm -rf "${CLEANUP_DIRS[@]}"; fi' EXIT

main "$@"
