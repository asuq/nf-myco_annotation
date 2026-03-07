/*
 * Placeholder for PADLOC execution.
 */
process PADLOC {
    input:
    path clean_gff
    path faa

    output:
    path 'padloc', emit: outdir
    path 'versions.yml', emit: versions

    stub:
    '''
    mkdir -p padloc
    : > padloc/placeholder.txt
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "PADLOC is a placeholder module." >&2
    exit 1
    '''
}
