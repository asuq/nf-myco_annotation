/*
 * Placeholder for PADLOC input cleanup from Prokka outputs.
 */
process PREPARE_PADLOC_INPUTS {
    input:
    path gff
    path faa

    output:
    path 'padloc_input.gff', emit: gff
    path 'padloc_input.faa', emit: faa
    path 'versions.yml', emit: versions

    script:
    '''
    echo "PREPARE_PADLOC_INPUTS is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > padloc_input.gff
    : > padloc_input.faa
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
