/*
 * Placeholder for an offline BUSCO run.
 */
process BUSCO {
    input:
    tuple val(meta), path(genome)
    val lineage
    path dataset_dir

    output:
    tuple val(meta), val(lineage), path('busco_output'), emit: outdir
    tuple val(meta), val(lineage), path('short_summary.json'), emit: summary
    path 'versions.yml', emit: versions

    script:
    '''
    echo "BUSCO is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    mkdir -p busco_output
    : > short_summary.json
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
