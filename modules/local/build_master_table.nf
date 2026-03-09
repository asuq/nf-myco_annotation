/*
 * Assemble the final master table from the validated manifest, metadata block,
 * and the stable derived-summary tables collected earlier in the workflow.
 */
process BUILD_MASTER_TABLE {
    tag "master_table"
    label 'process_single'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['master_table.tsv']
                ? filename
                : null
        },
    )

    input:
    path validated_samples
    path metadata
    path append_columns
    path taxonomy
    path checkm2
    path sixteen_s_status
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path ccfinder_strains
    path ani_summary
    path assembly_stats

    output:
    path 'master_table.tsv', emit: master_table
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    """
    build_master_table.py \
        --validated-samples "${validated_samples}" \
        --metadata "${metadata}" \
        --append-columns "${append_columns}" \
        --taxonomy "${taxonomy}" \
        --checkm2 "${checkm2}" \
        --16s-status "${sixteen_s_status}" \
        ${buscoArgs} \
        --ccfinder-strains "${ccfinder_strains}" \
        --ani "${ani_summary}" \
        --assembly-stats "${assembly_stats}" \
        --output master_table.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_master_table.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > master_table.tsv
    Accession	Tax_ID	16S	Gcode	Low_quality	Cluster_ID
    sample_a	123	Yes	4	false	cluster_1
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_master_table.py"
    EOF
    '''
}
