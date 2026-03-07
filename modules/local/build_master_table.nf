/*
 * Placeholder for final master-table assembly.
 */
process BUILD_MASTER_TABLE {
    input:
    path table_inputs

    output:
    path 'master_table.csv', emit: master_table
    path 'versions.yml', emit: versions

    stub:
    '''
    : > master_table.csv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "BUILD_MASTER_TABLE is a placeholder module." >&2
    exit 1
    '''
}
