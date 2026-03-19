/*
 * Build the authoritative sample-status table by overlaying downstream module
 * outcomes onto the validation-time status seed.
 */
process WRITE_SAMPLE_STATUS {
    tag "sample_status"
    label 'process_single'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'sample_status.tsv' ? filename : null },
    )

    input:
    path validated_samples
    path initial_status
    path columns
    path metadata
    path taxonomy
    path checkm2
    path sixteen_s_status
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path codetta_summary
    path ccfinder_strains
    path prokka_manifest
    path padloc_manifest
    path eggnog_manifest
    path ani_summary
    path assembly_stats
    val primary_busco_column

    output:
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    """
    script_path="\$(command -v build_sample_status.py)"
    python3 "\${script_path}" \
        --validated-samples "${validated_samples}" \
        --initial-status "${initial_status}" \
        --metadata "${metadata}" \
        --taxonomy "${taxonomy}" \
        --checkm2 "${checkm2}" \
        --16s-status "${sixteen_s_status}" \
        ${buscoArgs} \
        --codetta-summary "${codetta_summary}" \
        --ccfinder-strains "${ccfinder_strains}" \
        --prokka-manifest "${prokka_manifest}" \
        --padloc-manifest "${padloc_manifest}" \
        --eggnog-manifest "${eggnog_manifest}" \
        --ani "${ani_summary}" \
        --assembly-stats "${assembly_stats}" \
        --primary-busco-column "${primary_busco_column}" \
        --columns "${columns}" \
        --output sample_status.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_sample_status.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > sample_status.tsv
    accession	internal_id	is_new	validation_status	taxonomy_status	barrnap_status	checkm2_gcode4_status	checkm2_gcode11_status	gcode_status	gcode	low_quality	busco_bacillota_odb12_status	busco_mycoplasmatota_odb12_status	codetta_status	prokka_status	ccfinder_status	padloc_status	eggnog_status	ani_included	ani_exclusion_reason	warnings	notes
    sample_a	sample_a	false	done	done	done	done	done	done	4	false	done	done	done	done	done	done	done	true		stub_warning	stub note
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_sample_status.py"
    EOF
    '''
}
