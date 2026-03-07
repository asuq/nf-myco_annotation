/*
 * Placeholder for authoritative sample-status table writing.
 */
process WRITE_SAMPLE_STATUS {
    input:
    path status_inputs

    output:
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    stub:
    '''
    : > sample_status.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "WRITE_SAMPLE_STATUS is a placeholder module." >&2
    exit 1
    '''
}
