/*
 * Placeholder for internal-ID-safe FASTA staging.
 */
process STAGE_INPUTS {
    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path('staged.fasta'), emit: staged_fasta
    tuple val(meta), path('staged.fasta.fai'), emit: fai, optional: true
    path 'versions.yml', emit: versions

    script:
    '''
    echo "STAGE_INPUTS is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > staged.fasta
    : > staged.fasta.fai
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
