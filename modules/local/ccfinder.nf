/*
 * Placeholder for CRISPRCasFinder execution.
 */
process CCFINDER {
    input:
    tuple val(meta), path(genome)
    val gcode

    output:
    tuple val(meta), path('ccfinder'), emit: outdir
    tuple val(meta), path('result.json'), emit: result_json
    path 'versions.yml', emit: versions

    script:
    '''
    echo "CCFINDER is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    mkdir -p ccfinder
    : > result.json
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
