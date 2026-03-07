/*
 * Placeholder for ANI clustering and representative selection.
 */
process CLUSTER_ANI {
    input:
    path ani_matrix
    path ani_metadata

    output:
    path 'cluster.tsv', emit: clusters
    path 'representatives.tsv', emit: representatives
    path 'versions.yml', emit: versions

    script:
    '''
    echo "CLUSTER_ANI is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > cluster.tsv
    : > representatives.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
