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
    path append_columns
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
    """
    build_master_table.py \
        --validated-samples "${validated_samples}" \
        --metadata "${metadata}" \
        --append-columns "${append_columns}" \
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
    """
    metadata_header="\$(head -n 1 "${metadata}")"
    metadata_row="\$(awk 'NR == 2 { print; exit }' "${metadata}")"
    append_header="\$(paste -sd '\t' "${append_columns}")"
    printf '%s\t%s\n' "\${metadata_header}" "\${append_header}" > master_table.tsv
    append_values=()
    while IFS= read -r column; do
        case "\${column}" in
            Gcode)
                append_values+=("4")
                ;;
            Codetta_Genetic_Code)
                append_values+=("FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
                ;;
            Codetta_NCBI_Table_Candidates)
                append_values+=("1;11")
                ;;
            Low_quality)
                append_values+=("false")
                ;;
            16S)
                append_values+=("Yes")
                ;;
            BUSCO_*)
                append_values+=("C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200")
                ;;
            CRISPRS)
                append_values+=("2")
                ;;
            SPACERS_SUM)
                append_values+=("7")
                ;;
            CRISPR_FRAC)
                append_values+=("0.1")
                ;;
            Cluster_ID)
                append_values+=("cluster_1")
                ;;
            Is_Representative)
                append_values+=("yes")
                ;;
            ANI_to_Representative)
                append_values+=("100")
                ;;
            Score)
                append_values+=("0.95")
                ;;
            *)
                append_values+=("NA")
                ;;
        esac
    done < "${append_columns}"
    tab_char="\$(printf '\t')"
    append_row="\$(IFS="\${tab_char}"; printf '%s' "\${append_values[*]}")"
    printf '%s\t%s\n' "\${metadata_row}" "\${append_row}" >> master_table.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_master_table.py"
    EOF
    """.stripIndent()
}
