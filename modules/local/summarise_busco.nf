/*
 * Placeholder for BUSCO summary parsing.
 */
process SUMMARISE_BUSCO {
    input:
    path busco_summaries

    output:
    path 'busco_summary.tsv', emit: summary
    path 'versions.yml', emit: versions

    script:
    '''
    echo "SUMMARISE_BUSCO is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > busco_summary.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
