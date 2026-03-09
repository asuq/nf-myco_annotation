/*
 * Run CRISPRCasFinder with the locked standalone invocation semantics.
 */
process CCFINDER {
    tag "${meta.accession}"
    label 'process_high'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/ccfinder" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(genome), val(gcode)

    output:
    tuple val(meta), path('ccfinder'), path('result.json'), path('ccfinder.log'), emit: results
    tuple val(meta), path('result.json'), emit: result_json
    path 'versions.yml', emit: versions

    script:
    def ccfinderRoot = '/usr/local/CRISPRCasFinder'
    def extraArgs = (params.ccfinder_extra_args ?: '').toString()
    """
    ccfinder_root='${ccfinderRoot}'
    if [[ ! -f "\${ccfinder_root}/CRISPRCasFinder.pl" ]]; then
        echo "CRISPRCasFinder.pl not found under \${ccfinder_root}." >&2
        exit 1
    fi
    real_muscle="\$(command -v muscle || true)"
    if [[ -z "\${real_muscle}" ]]; then
        echo "muscle not found in PATH." >&2
        exit 1
    fi
    real_mv="\$(command -v mv || true)"
    if [[ -z "\${real_mv}" ]]; then
        echo "mv not found in PATH." >&2
        exit 1
    fi

    task_root="\$PWD"
    tool_bin="\${task_root}/tool_bin"
    run_root="\${task_root}/ccfinder_run"
    tool_output_root="\${task_root}/ccfinder_raw"
    genome_path="\$(cd "\$(dirname "${genome}")" && pwd)/\$(basename "${genome}")"
    genome_name="\$(basename "${genome}")"

    mkdir -p "\${tool_bin}" "\${run_root}"
    {
        printf '%s\n' '#!/usr/bin/env bash'
        printf '%s\n' 'set -euo pipefail'
        printf 'real_muscle=%q\n' "\${real_muscle}"
        printf '%s\n' 'translated_args=()'
        printf '%s\n' 'while ((\$# > 0)); do'
        printf '%s\n' '    case "\$1" in'
        printf '%s\n' '        -align)'
        printf '%s\n' '            translated_args+=(-in "\$2")'
        printf '%s\n' '            shift 2'
        printf '%s\n' '            ;;'
        printf '%s\n' '        -output)'
        printf '%s\n' '            translated_args+=(-out "\$2")'
        printf '%s\n' '            shift 2'
        printf '%s\n' '            ;;'
        printf '%s\n' '        *)'
        printf '%s\n' '            translated_args+=("\$1")'
        printf '%s\n' '            shift'
        printf '%s\n' '            ;;'
        printf '%s\n' '    esac'
        printf '%s\n' 'done'
        printf '%s\n' 'exec "\$real_muscle" "\${translated_args[@]}"'
    } > "\${tool_bin}/muscle"
    chmod +x "\${tool_bin}/muscle"
    {
        printf '%s\n' '#!/usr/bin/env bash'
        printf '%s\n' 'set -euo pipefail'
        printf 'real_mv=%q\n' "\${real_mv}"
        printf 'protected=(%q %q %q)\n' "\$(basename "\${tool_bin}")" "\$(basename "\${run_root}")" "\$(basename "\${tool_output_root}")"
        printf '%s\n' 'if (( \$# <= 1 )); then'
        printf '%s\n' '    exec "\$real_mv" "\$@"'
        printf '%s\n' 'fi'
        printf '%s\n' 'destination="\${!#}"'
        printf '%s\n' 'filtered_sources=()'
        printf '%s\n' 'for ((i = 1; i < \$#; i++)); do'
        printf '%s\n' '    arg="\${!i}"'
        printf '%s\n' '    source_name="\$(basename "\$arg")"'
        printf '%s\n' '    keep_source=true'
        printf '%s\n' '    for protected_name in "\${protected[@]}"; do'
        printf '%s\n' '        if [[ "\$source_name" == "\$protected_name" ]]; then'
        printf '%s\n' '            keep_source=false'
        printf '%s\n' '            break'
        printf '%s\n' '        fi'
        printf '%s\n' '    done'
        printf '%s\n' '    if [[ "\$keep_source" == true ]]; then'
        printf '%s\n' '        filtered_sources+=("\$arg")'
        printf '%s\n' '    fi'
        printf '%s\n' 'done'
        printf '%s\n' 'if (( \${#filtered_sources[@]} == 0 )); then'
        printf '%s\n' '    exit 0'
        printf '%s\n' 'fi'
        printf '%s\n' 'exec "\$real_mv" "\${filtered_sources[@]}" "\$destination"'
    } > "\${tool_bin}/mv"
    chmod +x "\${tool_bin}/mv"
    export PATH="\${tool_bin}:\$PATH"

    awk -v output_root="\${task_root}" '
        /^>/ {
            if (output_path != "") {
                close(output_path)
            }
            contig_id = \$1
            sub(/^>/, "", contig_id)
            sub(/\\.[0-9]+\$/, "", contig_id)
            output_path = output_root "/" contig_id ".fna"
            print \$0 > output_path
            next
        }
        {
            sequence_line = \$0
            while (length(sequence_line) > 60) {
                print substr(sequence_line, 1, 60) >> output_path
                sequence_line = substr(sequence_line, 61)
            }
            if (length(sequence_line) > 0) {
                print sequence_line >> output_path
            }
        }
    ' "\${genome_path}"

    pushd "\${run_root}" >/dev/null
    cp "\${genome_path}" "\${genome_name}"
    : > index.html
    set +e
    perl "\${ccfinder_root}/CRISPRCasFinder.pl" -in "\${genome_name}" \
        -outdir "\${tool_output_root}" \
        -soFile "\${ccfinder_root}/sel392v2.so" \
        -DBcrispr "\${ccfinder_root}/supplementary_files/CRISPR_crisprdb.csv" \
        -repeats "\${ccfinder_root}/supplementary_files/Repeat_List.csv" \
        -DIRrepeat "\${ccfinder_root}/supplementary_files/repeatDirection.tsv" \
        -cpuMacSyFinder ${task.cpus} -cpuProkka ${task.cpus} \
        -log -html -levelMin 2 \
        -cas -ccvRep -getSummaryCasfinder -gcode "${gcode}" \
        ${extraArgs}
    exit_code=\$?
    set -e
    popd >/dev/null

    mkdir -p ccfinder
    result_json_path=\$(find "\${tool_output_root}" "\${run_root}" -type f -name 'result.json' | head -n 1 || true)
    internal_log_path=\$(find "\${tool_output_root}" "\${run_root}" -type f \\( -name 'ccfinder.log' -o -name 'logFile_*' \\) | head -n 1 || true)
    if [[ -n "\${result_json_path}" ]]; then
        perl -0pi -e 's/:\\s*(?=,|\\}|\\])/: null/g' "\${result_json_path}"
    fi
    if [[ -d "\${tool_output_root}" ]]; then
        cp -R "\${tool_output_root}/". ccfinder/
    fi

    if [[ -n "\${result_json_path}" ]]; then
        cp "\${result_json_path}" result.json
    else
        : > result.json
    fi

    if [[ -n "\${internal_log_path}" ]]; then
        cp "\${internal_log_path}" ccfinder.log
    else
        : > ccfinder.log
    fi
    printf 'exit_code=%s\n' "\$exit_code" >> ccfinder.log

    {
        printf '"%s":\n' "${task.process}"
        printf '  crisprcasfinder: "%s"\n' "\$(perl "\${ccfinder_root}/CRISPRCasFinder.pl" -v 2>&1 | head -n 1 || echo NA)"
    } > versions.yml
    """

    stub:
    '''
    mkdir -p ccfinder
    cat <<'EOF' > result.json
    {"Sequences": []}
    EOF
    cp result.json ccfinder/result.json
    : > ccfinder.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      crisprcasfinder: "stub"
    EOF
    '''
}
