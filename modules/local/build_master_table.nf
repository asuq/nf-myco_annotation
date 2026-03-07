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
            filename in ['master_table.tsv', 'sample_status.tsv']
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
    path busco_tables
    path ccfinder_strains
    path ani

    output:
    path 'master_table.tsv', emit: master_table
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    """
    python3 "${projectDir}/bin/build_master_table.py" \
        --validated-samples "${validated_samples}" \
        --metadata "${metadata}" \
        --append-columns "${append_columns}" \
        --taxonomy "${taxonomy}" \
        --checkm2 "${checkm2}" \
        --16s-status "${sixteen_s_status}" \
        ${buscoArgs} \
        --ccfinder-strains "${ccfinder_strains}" \
        --ani "${ani}" \
        --output master_table.tsv \
        --sample-status-output sample_status.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_master_table.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > master_table.tsv
    Accession	Tax_ID	16S	Gcode	Low_quality	Cluster_ID
    sample_a	123	Yes	4	false	cluster_1
    EOF
    cat <<'EOF' > sample_status.tsv
    accession	internal_id	is_new	validation_status	taxonomy_status	barrnap_status	checkm2_gcode4_status	checkm2_gcode11_status	gcode_status	gcode	low_quality	busco_bacillota_odb12_status	busco_mycoplasmatota_odb12_status	prokka_status	ccfinder_status	padloc_status	eggnog_status	ani_included	ani_exclusion_reason	warnings	notes
    sample_a	sample_a	false	done	done	done	done	done	done	4	false	done	done	na	done	na	na	true
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_master_table.py"
    EOF
    '''
}
