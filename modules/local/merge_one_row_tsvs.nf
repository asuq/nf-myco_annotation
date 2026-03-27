/*
 * Merge canonical one-row TSV fragments into one accession-sorted summary
 * table using a locked header contract.
 */
process MERGE_ONE_ROW_TSVS {
    tag "merge_tsv_rows"
    label 'process_single'
    cache 'deep'

    input:
    path header_asset, name: 'canonical_header.tsv'
    path row_files, name: 'rows/row??.tsv'

    output:
    path 'merged.tsv', emit: merged
    path 'versions.yml', emit: versions

    script:
    def rowFileList = row_files instanceof Collection ? row_files : [row_files]
    def rowArgs = rowFileList.collect { "--input \"${it}\"" }.join(' \\\n        ')
    """
    script_path="\$(command -v merge_one_row_tsvs.py)"
    python3 "\${script_path}" \
        --header "${header_asset}" \
        ${rowArgs} \
        --output merged.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/merge_one_row_tsvs.py"
    EOF
    """

    stub:
    def rowFileList = row_files instanceof Collection ? row_files : [row_files]
    def rowArgs = rowFileList.collect { "--input \"${it}\"" }.join(' \\\n        ')
    """
    script_path="\$(command -v merge_one_row_tsvs.py)"
    python3 "\${script_path}" \
        --header "${header_asset}" \
        ${rowArgs} \
        --output merged.tsv

    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/merge_one_row_tsvs.py"
    EOF
    """
}
