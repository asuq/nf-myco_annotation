/*
 * Placeholder for Barrnap rRNA detection.
 */
process BARRNAP {
    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path('rrna.gff'), emit: gff
    tuple val(meta), path('rrna.fa'), emit: fasta
    path 'versions.yml', emit: versions

    script:
    '''
    echo "BARRNAP is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > rrna.gff
    : > rrna.fa
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
