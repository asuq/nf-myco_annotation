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

    set +e
    perl "\${ccfinder_root}/CRISPRCasFinder.pl" -in "${genome}" \
        -outdir "\$PWD/ccfinder" \
        -soFile "\${ccfinder_root}/sel392v2.so" \
        -DBcrispr "\${ccfinder_root}/supplementary_files/CRISPR_crisprdb.csv" \
        -repeats "\${ccfinder_root}/supplementary_files/Repeat_List.csv" \
        -DIRrepeat "\${ccfinder_root}/supplementary_files/repeatDirection.tsv" \
        -CASFinder "\${ccfinder_root}/CasFinder-2.0.3" \
        -cpuMacSyFinder ${task.cpus} -cpuProkka ${task.cpus} \
        -log -html -levelMin 2 \
        -cas -ccvRep -getSummaryCasfinder -gcode "${gcode}" \
        ${extraArgs} \
        > ccfinder.log 2>&1
    exit_code=\$?
    set -e

    mkdir -p ccfinder
    result_json_path=\$(find ccfinder -type f -name 'result.json' | head -n 1 || true)
    if [[ -n "\${result_json_path}" ]]; then
        cp "\${result_json_path}" result.json
    else
        : > result.json
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> ccfinder.log

    cat <<EOF > versions.yml
    "${task.process}":
      crisprcasfinder: "$(perl "\${ccfinder_root}/CRISPRCasFinder.pl" -v 2>&1 | head -n 1 || echo NA)"
    EOF
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
