/*
 * Prepare the configured BUSCO lineage datasets under the canonical db_root.
 */
process PREP_BUSCO_DATABASES {
    tag "busco"
    label 'process_medium'

    input:
    tuple path(db_root), path(source_root_input), val(source_enabled), val(download_enabled), val(version), val(lineages), val(link_mode), val(scratch_root), val(force)

    output:
    path("busco_report.tsv"), emit: report
    path 'versions.yml', emit: versions

    script:
    def renderedLineages = (lineages as List<String>).join('\n')
    """
    script_path="\$(command -v prepare_runtime_databases.py)"
    db_root_path="\$(cd "${db_root}" && pwd -P)"
    destination_root="\${db_root_path}/busco"

    cat <<'EOF' > busco_lineages.txt
    ${renderedLineages}
    EOF

    helper_args=(
        --busco-dest-root
        "\${destination_root}"
        --report
        busco_report.tsv
        --link-mode
        "${link_mode}"
    )

    if [[ "${download_enabled}" == "true" ]]; then
        helper_args+=(--download)
    fi

    if [[ "${force}" == "true" ]]; then
        helper_args+=(--force)
    fi

    if [[ -n "${version ?: ''}" ]]; then
        helper_args+=(--busco-version "${version}")
    fi

    if [[ -n "${scratch_root ?: ''}" ]]; then
        helper_args+=(--scratch-root "${scratch_root}")
    fi

    if [[ "${source_enabled}" == "true" ]]; then
        source_root_path="\$(cd "\$(dirname "${source_root_input}")" && pwd -P)/\$(basename "${source_root_input}")"
        while IFS= read -r lineage; do
            candidate_path="\${source_root_path}/\${lineage}"
            if [[ ! -e "\${candidate_path}" && -e "\${candidate_path}.tar.gz" ]]; then
                candidate_path="\${candidate_path}.tar.gz"
            elif [[ ! -e "\${candidate_path}" && -e "\${candidate_path}.zip" ]]; then
                candidate_path="\${candidate_path}.zip"
            fi
            helper_args+=(--busco-lineage-source "\${lineage}=\${candidate_path}")
        done < busco_lineages.txt
    fi

    python3 "\${script_path}" "\${helper_args[@]}"

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      helper: "bin/prepare_runtime_databases.py"
    EOF
    """

    stub:
    """
    db_root_path="\$(cd "${db_root}" && pwd -P)"
    destination_root="\${db_root_path}/busco"
    mkdir -p "\${destination_root}"

    printf '%s\n' ${(lineages as List<String>).collect { "'${it}'" }.join(' ')} > busco_lineages.txt

    {
        printf 'component\tstatus\tsource\tdestination\tdetails\n'
        while IFS= read -r lineage; do
        mkdir -p "\${destination_root}/\${lineage}"
        : > "\${destination_root}/\${lineage}/dataset.cfg"
            printf 'busco:%s\tprepared\tstub\t%s/%s\trequired_paths=dataset.cfg;source_mode=stub\n' "\${lineage}" "\${destination_root}" "\${lineage}"
        done < busco_lineages.txt
        printf 'busco_root\tprepared\tderived\t%s\trequired_paths=stub;source_mode=stub\n' "\${destination_root}"
    } > busco_report.tsv

    printf '%s\n' \
        '"${task.process}":' \
        '  python: "stub"' \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """
}
