/*
 * Placeholder for Prokka annotation.
 */
process PROKKA {
    input:
    tuple val(meta), path(genome)
    val gcode

    output:
    tuple val(meta), path('prokka'), emit: outdir
    path 'versions.yml', emit: versions

    script:
    '''
    echo "PROKKA is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    mkdir -p prokka
    : > prokka/placeholder.txt
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
