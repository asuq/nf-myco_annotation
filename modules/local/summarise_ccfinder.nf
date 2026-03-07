/*
 * Placeholder for CRISPRCasFinder result parsing.
 */
process SUMMARISE_CCFINDER {
    input:
    path ccfinder_results

    output:
    path 'ccfinder_strains.tsv', emit: strains
    path 'ccfinder_contigs.tsv', emit: contigs
    path 'ccfinder_crisprs.tsv', emit: crisprs
    path 'versions.yml', emit: versions

    stub:
    '''
    : > ccfinder_strains.tsv
    : > ccfinder_contigs.tsv
    : > ccfinder_crisprs.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "SUMMARISE_CCFINDER is a placeholder module." >&2
    exit 1
    '''
}
