/*
 * Placeholder for taxonomy expansion from a pinned taxdump.
 */
process TAXONOMY_EXPAND {
    input:
    path validated_samples
    path metadata
    path taxdump

    output:
    path 'taxonomy_expanded.tsv', emit: taxonomy
    path 'versions.yml', emit: versions

    script:
    '''
    echo "TAXONOMY_EXPAND is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > taxonomy_expanded.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
