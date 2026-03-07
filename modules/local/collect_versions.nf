/*
 * Placeholder for merged tool and database version reporting.
 */
process COLLECT_VERSIONS {
    input:
    path version_files

    output:
    path 'tool_and_db_versions.tsv', emit: versions_table

    script:
    '''
    echo "COLLECT_VERSIONS is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > tool_and_db_versions.tsv
    '''
}
