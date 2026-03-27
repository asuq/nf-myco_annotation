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
    path sample_status_columns

    output:
    path 'validated_samples.tsv', emit: validated_samples
    path 'accession_map.tsv', emit: accession_map
    path 'validation_warnings.tsv', emit: validation_warnings
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    script:
    """validate_inputs.py \
    --sample-csv "${sample_csv}" \
    --metadata "${metadata}" \
    --sample-status-columns "${sample_status_columns}" \
    --defer-genome-fasta-check \
    --outdir .

cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/validate_inputs.py"
EOF
"""

    stub:
    def stubGenome = file("${projectDir}/assets/testdata/stub/genomes/TEST_ACC.fasta").toString()
    """cat <<'EOF' > validated_samples.tsv
accession	is_new	assembly_level	genome_fasta	internal_id
TEST_ACC	false	NA	${stubGenome}	TEST_ACC
EOF
cat <<'EOF' > accession_map.tsv
accession	internal_id	is_new	assembly_level	genome_fasta	metadata_present
TEST_ACC	TEST_ACC	false	NA	${stubGenome}	true
EOF
cat <<'EOF' > validation_warnings.tsv
accession	warning_code	message
EOF
header="\$(paste -sd '\t' "${sample_status_columns}")"
printf '%s\n' "\${header}" > sample_status.tsv
status_values=()
while IFS= read -r column; do
    case "\${column}" in
        accession)
            status_values+=("TEST_ACC")
            ;;
        internal_id)
            status_values+=("TEST_ACC")
            ;;
        is_new)
            status_values+=("false")
            ;;
        validation_status)
            status_values+=("done")
            ;;
        warnings)
            status_values+=("stub_warning")
            ;;
        notes)
            status_values+=("stub warning")
            ;;
        *_status|ani_included)
            status_values+=("na")
            ;;
        gcode|low_quality)
            status_values+=("NA")
            ;;
        *)
            status_values+=("")
            ;;
    esac
done < "${sample_status_columns}"
tab_char="\$(printf '\t')"
printf '%s\n' "\$(IFS="\${tab_char}"; printf '%s' "\${status_values[*]}")" >> sample_status.tsv
cat <<'EOF' > versions.yml
"${task.process}":
  python: "stub"
  script: "bin/validate_inputs.py"
EOF
"""
}
