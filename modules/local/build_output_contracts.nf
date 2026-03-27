/*
 * Build runtime reporting contracts from the configured BUSCO lineage list.
 */
process BUILD_OUTPUT_CONTRACTS {
    tag "output_contracts"
    label 'process_single'
    cache 'deep'

    input:
    val busco_lineages

    output:
    path 'master_table_append_columns.txt', emit: append_columns
    path 'sample_status_columns.txt', emit: sample_status_columns
    path 'versions.yml', emit: versions

    script:
    def lineageArgs = (busco_lineages as List<String>).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n        ')
    """script_path="\$(command -v build_output_contracts.py)"
python3 "\${script_path}" \
    ${lineageArgs} \
    --append-columns-output master_table_append_columns.txt \
    --sample-status-columns-output sample_status_columns.txt

cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/build_output_contracts.py"
EOF
"""

    stub:
    def appendColumns = [
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
        'Gcode',
        'Codetta_Genetic_Code',
        'Codetta_NCBI_Table_Candidates',
        'Low_quality',
        '16S',
        *((busco_lineages as List<String>).collect { "BUSCO_${it}" }),
        'CRISPRS',
        'SPACERS_SUM',
        'CRISPR_FRAC',
        'Cluster_ID',
        'Is_Representative',
        'ANI_to_Representative',
        'Score',
    ].join('\n')
    def sampleStatusColumns = [
        'accession',
        'internal_id',
        'is_new',
        'validation_status',
        'taxonomy_status',
        'barrnap_status',
        'checkm2_gcode4_status',
        'checkm2_gcode11_status',
        'gcode_status',
        'gcode',
        'low_quality',
        *((busco_lineages as List<String>).collect { "busco_${it}_status" }),
        'codetta_status',
        'prokka_status',
        'ccfinder_status',
        'padloc_status',
        'eggnog_status',
        'ani_included',
        'ani_exclusion_reason',
        'warnings',
        'notes',
    ].join('\n')
    """cat <<'EOF' > master_table_append_columns.txt
${appendColumns}
EOF
cat <<'EOF' > sample_status_columns.txt
${sampleStatusColumns}
EOF
cat <<'EOF' > versions.yml
"${task.process}":
  python: "stub"
  script: "bin/build_output_contracts.py"
EOF
"""
}
