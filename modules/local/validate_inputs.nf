/*
 * Validate the sample manifest and metadata table before any per-sample work.
 */
process VALIDATE_INPUTS {
    tag "${sample_csv.baseName}"
    label 'process_single'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['validated_samples.tsv', 'accession_map.tsv', 'validation_warnings.tsv', 'sample_status.tsv']
                ? filename
                : null
        },
    )

    input:
    path sample_csv
    path metadata

    output:
    path 'validated_samples.tsv', emit: validated_samples
    path 'accession_map.tsv', emit: accession_map
    path 'validation_warnings.tsv', emit: validation_warnings
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    script:
    """
    validate_inputs.py \
        --sample-csv "${sample_csv}" \
        --metadata "${metadata}" \
        --outdir .

    cat <<EOF > versions.yml
    "${task.process}":
      python: "$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/validate_inputs.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > validated_samples.tsv
    accession	is_new	assembly_level	genome_fasta	internal_id
    sample_a	true	Scaffold	/work/sample_a.fna	sample_a
    EOF
    cat <<'EOF' > accession_map.tsv
    accession	internal_id	is_new	assembly_level	genome_fasta	metadata_present
    sample_a	sample_a	true	Scaffold	/work/sample_a.fna	true
    EOF
    cat <<'EOF' > validation_warnings.tsv
    accession	warning_code	message
    EOF
    cat <<'EOF' > sample_status.tsv
    accession	internal_id	is_new	validation_status	taxonomy_status	barrnap_status	checkm2_gcode4_status	checkm2_gcode11_status	gcode_status	gcode	low_quality	busco_bacillota_odb12_status	busco_mycoplasmatota_odb12_status	prokka_status	ccfinder_status	padloc_status	eggnog_status	ani_included	ani_exclusion_reason	warnings	notes
    sample_a	sample_a	true	done	na	na	na	na	na	NA	NA	na	na	na	na	na	na	na		internal_id_collision_resolved	stub warning
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/validate_inputs.py"
    EOF
    '''
}
