/*
 * Build the authoritative sample-status table by overlaying downstream module
 * outcomes onto the validation-time status seed.
 */
process WRITE_SAMPLE_STATUS {
    tag "sample_status"
    label 'process_single'
    cache 'deep'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'sample_status.tsv' ? filename : null },
    )

    input:
    path validated_samples
    path initial_status, name: 'initial_status.tsv'
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
    """.stripIndent()

    stub:
    """
    header="\$(paste -sd '\t' "${columns}")"
    printf '%s\n' "\${header}" > sample_status.tsv
    status_values=()
    while IFS= read -r column; do
        case "\${column}" in
            accession)
                status_values+=("sample_a")
                ;;
            internal_id)
                status_values+=("sample_a")
                ;;
            is_new)
                status_values+=("false")
                ;;
            gcode)
                status_values+=("4")
                ;;
            low_quality)
                status_values+=("false")
                ;;
            ani_included)
                status_values+=("true")
                ;;
            ani_exclusion_reason)
                status_values+=("")
                ;;
            warnings)
                status_values+=("stub_warning")
                ;;
            notes)
                status_values+=("stub note")
                ;;
            *_status)
                status_values+=("done")
                ;;
            *)
                status_values+=("")
                ;;
        esac
    done < "${columns}"
    tab_char="\$(printf '\t')"
    printf '%s\n' "\$(IFS="\${tab_char}"; printf '%s' "\${status_values[*]}")" >> sample_status.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_sample_status.py"
    EOF
    """.stripIndent()
}
