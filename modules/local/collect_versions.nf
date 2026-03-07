/*
 * Placeholder for merged tool and database version reporting.
 */
process COLLECT_VERSIONS {
    input:
    path version_files

    output:
    path 'tool_and_db_versions.tsv', emit: versions_table

    stub:
    '''
    : > tool_and_db_versions.tsv
    '''

    script:
    '''
    echo "COLLECT_VERSIONS is a placeholder module." >&2
    exit 1
    '''
}
