/*
 * Placeholder for a single CheckM2 invocation.
 */
process CHECKM2 {
    input:
    tuple val(meta), path(genome)
    val translation_table

    output:
    tuple val(meta), val(translation_table), path('checkm2_output'), emit: outdir
    tuple val(meta), val(translation_table), path('quality_report.tsv'), emit: quality_report
    path 'versions.yml', emit: versions

    script:
    '''
    echo "CHECKM2 is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    mkdir -p checkm2_output
    : > quality_report.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
