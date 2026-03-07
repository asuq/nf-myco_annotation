/*
 * Placeholder for sample and metadata validation.
 */
process VALIDATE_INPUTS {
    input:
    path sample_csv
    path metadata

    output:
    path 'validated_samples.tsv', emit: validated_samples
    path 'accession_map.tsv', emit: accession_map
    path 'validation_warnings.tsv', emit: validation_warnings
    path 'versions.yml', emit: versions

    stub:
    '''
    : > validated_samples.tsv
    : > accession_map.tsv
    : > validation_warnings.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "VALIDATE_INPUTS is a placeholder module." >&2
    exit 1
    '''
}
