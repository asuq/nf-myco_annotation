/*
 * Placeholder for ANI eligibility filtering and FastANI manifest building.
 */
process BUILD_FASTANI_INPUTS {
    input:
    path sample_tables

    output:
    path 'fastani_paths.txt', emit: paths
    path 'ani_metadata.tsv', emit: metadata
    path 'ani_exclusions.tsv', emit: exclusions
    path 'versions.yml', emit: versions

    script:
    '''
    echo "BUILD_FASTANI_INPUTS is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > fastani_paths.txt
    : > ani_metadata.tsv
    : > ani_exclusions.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
