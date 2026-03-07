/*
 * Placeholder for all-vs-all FastANI.
 */
process FASTANI {
    input:
    path fastani_paths

    output:
    path 'fastani.matrix', emit: matrix
    path 'fastani.tsv', emit: raw_output
    path 'versions.yml', emit: versions

    stub:
    '''
    : > fastani.matrix
    : > fastani.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "FASTANI is a placeholder module." >&2
    exit 1
    '''
}
