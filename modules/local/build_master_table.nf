/*
 * Assemble the final master table from the validated manifest, metadata block,
 * and the stable derived-summary tables collected earlier in the workflow.
 */
process BUILD_MASTER_TABLE {
    tag "master_table"
    label 'process_single'
    cache 'deep'
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
    val busco_lineages
    path taxonomy
    path checkm2
    path sixteen_s_status
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path codetta_summary
    path ccfinder_strains
    path ani_summary
    path assembly_stats

    output:
    path 'master_table.tsv', emit: master_table
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    def lineageArgs = (busco_lineages as List<String>).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n        ')
    """
    build_master_table.py \
        --validated-samples "${validated_samples}" \
        --metadata "${metadata}" \
        ${lineageArgs} \
        --taxonomy "${taxonomy}" \
        --checkm2 "${checkm2}" \
        --16s-status "${sixteen_s_status}" \
        ${buscoArgs} \
        --codetta-summary "${codetta_summary}" \
        --ccfinder-strains "${ccfinder_strains}" \
        --ani "${ani_summary}" \
        --assembly-stats "${assembly_stats}" \
        --output master_table.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_master_table.py"
    EOF
    """.stripIndent()

    stub:
    def appendColumns = [
        'is_new',
        'superkingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'Completeness_gcode4',
        'Completeness_gcode11',
        'Contamination_gcode4',
        'Contamination_gcode11',
        'Coding_Density_gcode4',
        'Coding_Density_gcode11',
        'Average_Gene_Length_gcode4',
        'Average_Gene_Length_gcode11',
        'Total_Coding_Sequences_gcode4',
        'Total_Coding_Sequences_gcode11',
        'GC_Content',
        'Gcode',
        'Codetta_Genetic_Code',
        'Codetta_NCBI_Table_Candidates',
        'Low_quality',
        '16S',
    ] + (busco_lineages as List<String>).collect { "BUSCO_${it}" } + [
        'CRISPRS',
        'SPACERS_SUM',
        'CRISPR_FRAC',
        'Cluster_ID',
        'Is_Representative',
        'ANI_to_Representative',
        'Score',
    ]
    def appendStubValues = [
        is_new: 'false',
        Gcode: '4',
        GC_Content: '50',
        Codetta_Genetic_Code: 'FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        Codetta_NCBI_Table_Candidates: '1;11',
        Low_quality: 'false',
        '16S': 'Yes',
        CRISPRS: '2',
        SPACERS_SUM: '7',
        CRISPR_FRAC: '0.1',
        Cluster_ID: 'cluster_1',
        Is_Representative: 'yes',
        ANI_to_Representative: '100',
        Score: '0.95',
    ]
    def appendRow = appendColumns.collect { column ->
        if (appendStubValues.containsKey(column)) {
            return appendStubValues[column]
        }
        return column.startsWith('BUSCO_')
            ? 'C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200'
            : 'NA'
    }.join('\t')
    """
    metadata_header="\$(head -n 1 "${metadata}")"
    metadata_row="\$(awk 'NR == 2 { print; exit }' "${metadata}")"
    printf '%s\t%s\n' "\${metadata_header}" "${appendColumns.join('\t')}" > master_table.tsv
    printf '%s\t%s\n' "\${metadata_row}" "${appendRow}" >> master_table.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_master_table.py"
    EOF
    """.stripIndent()
}
