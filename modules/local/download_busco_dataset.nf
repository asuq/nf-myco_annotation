/*
 * Placeholder for BUSCO dataset download preparation.
 */
process DOWNLOAD_BUSCO_DATASET {
    input:
    val lineage

    output:
    tuple val(lineage), path('busco_dataset'), emit: dataset
    tuple val(lineage), path('download.log'), emit: log
    path 'versions.yml', emit: versions

    script:
    '''
    echo "DOWNLOAD_BUSCO_DATASET is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    mkdir -p busco_dataset
    : > download.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
