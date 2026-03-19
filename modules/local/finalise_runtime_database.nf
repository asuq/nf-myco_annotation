/*
 * Validate one prepared runtime database directory and emit one report row.
 */
process FINALISE_RUNTIME_DATABASE {
    tag "${component}"
    label 'process_single'
    label 'finalise_runtime_database'

    input:
    tuple val(component), val(destination), path(mode_file), val(source_label), path(lineages_file)

    output:
    path("${component}_report.tsv"), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    script_path="/usr/local/bin/finalise_runtime_database.py"
    mode="\$(tr -d '\n' < "${mode_file}")"

    helper_args=(
        --component
        "${component}"
        --destination
        "${destination}"
        --report
        "${component}_report.tsv"
        --mode
        "\${mode}"
        --source
        "${source_label}"
    )

    while IFS= read -r lineage; do
        [[ -z "\${lineage}" ]] && continue
        helper_args+=(--busco-lineage "\${lineage}")
    done < "${lineages_file}"

    /usr/local/bin/python3 "\${script_path}" "\${helper_args[@]}"

    python_version="\$(/usr/local/bin/python3 --version 2>&1 | sed 's/^Python //')"
    {
        printf '"%s":\n' "${task.process}"
        printf '  python: "%s"\n' "\${python_version}"
        printf '  helper: "%s"\n' 'bin/finalise_runtime_database.py'
    } > versions.yml
    """

    stub:
    """
    mode="\$(tr -d '\n' < "${mode_file}")"
    details=""
    case "${component}" in
        checkm2)
            details='required_paths=CheckM2_database.dmnd;dmnd=CheckM2_database.dmnd'
            ;;
        codetta)
            details='required_paths=Pfam-A_enone.hmm,Pfam-A_enone.hmm.h3f,Pfam-A_enone.hmm.h3i,Pfam-A_enone.hmm.h3m,Pfam-A_enone.hmm.h3p;profile=Pfam-A_enone.hmm'
            ;;
        busco_root)
            required_paths=()
            while IFS= read -r lineage; do
                [[ -z "\${lineage}" ]] && continue
                required_paths+=("\${lineage}/dataset.cfg")
            done < "${lineages_file}"
            details="required_paths=\$(IFS=,; echo "\${required_paths[*]}");lineages=\$(paste -sd ';' "${lineages_file}")"
            ;;
        eggnog)
            details='required_paths=eggnog.db,eggnog_proteins.dmnd;files=2'
            ;;
    esac
    status='prepared'
    if [[ "\${mode}" == 'reuse' ]]; then
        status='present'
    fi

    {
        printf 'component\tstatus\tsource\tdestination\tdetails\n'
        printf '%s\t%s\t%s\t%s\t%s\n' "${component}" "\${status}" "${source_label}" "${destination}" "\${details}"
    } > "${component}_report.tsv"

    {
        printf '"%s":\n' "${task.process}"
        printf '  python: "%s"\n' 'stub'
        printf '  helper: "%s"\n' 'bin/finalise_runtime_database.py'
    } > versions.yml
    """
}
