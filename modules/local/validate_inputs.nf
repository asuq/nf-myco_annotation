/*
 * Validate the sample manifest and metadata table before any per-sample work.
 */
process VALIDATE_INPUTS {
    tag "${sample_csv.baseName}"

    input:
    path sample_csv
    path metadata

    output:
    path 'validated_samples.tsv', emit: validated_samples
    path 'accession_map.tsv', emit: accession_map
    path 'validation_warnings.tsv', emit: validation_warnings
    path 'versions.yml', emit: versions

    script:
    """
    python3 "${projectDir}/bin/validate_inputs.py" \
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
    accession	is_new	assembly_level	genome_fasta	internal_id	metadata_present
    sample_a	true	Scaffold	/work/sample_a.fna	sample_a	true
    EOF
    cat <<'EOF' > accession_map.tsv
    accession	internal_id	is_new	assembly_level	genome_fasta	metadata_present
    sample_a	sample_a	true	Scaffold	/work/sample_a.fna	true
    EOF
    cat <<'EOF' > validation_warnings.tsv
    accession	warning_code	message
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/validate_inputs.py"
    EOF
    '''
}
