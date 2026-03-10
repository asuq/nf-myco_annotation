/*
 * Prepare one runtime database component under the canonical db_root.
 */
process PREP_RUNTIME_DATABASE {
    tag "${component}"
    label 'process_medium'

    input:
    tuple val(component), path(db_root), path(source_input), val(source_enabled), val(download_enabled), val(version), val(link_mode), val(scratch_root), val(force)

    output:
    path("${component}_report.tsv"), emit: report
    path 'versions.yml', emit: versions

    script:
    def sourceFlag = [
        taxdump: '--taxdump-source',
        checkm2: '--checkm2-source',
        eggnog: '--eggnog-source',
        padloc: '--padloc-source',
    ][component]
    def destinationFlag = [
        taxdump: '--taxdump-dest',
        checkm2: '--checkm2-dest',
        eggnog: '--eggnog-dest',
        padloc: '--padloc-dest',
    ][component]
    def versionFlag = [
        taxdump: '--taxdump-version',
        checkm2: '--checkm2-version',
        eggnog: '--eggnog-version',
        padloc: '--padloc-version',
    ][component]
    """
    script_path="\$(command -v prepare_runtime_databases.py)"
    db_root_path="\$(cd "${db_root}" && pwd -P)"
    destination_path="\${db_root_path}/${component}"

    helper_args=(
        "${destinationFlag}"
        "\${destination_path}"
        --report
        "${component}_report.tsv"
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
        helper_args+=("${versionFlag}" "${version}")
    fi

    if [[ -n "${scratch_root ?: ''}" ]]; then
        helper_args+=(--scratch-root "${scratch_root}")
    fi

    if [[ "${source_enabled}" == "true" ]]; then
        source_path="\$(cd "\$(dirname "${source_input}")" && pwd -P)/\$(basename "${source_input}")"
        helper_args+=("${sourceFlag}" "\${source_path}")
    fi

    python3 "\${script_path}" "\${helper_args[@]}"

    python_version="\$(python3 --version 2>&1 | sed 's/^Python //')"
    printf '%s\n' \
        '"${task.process}":' \
        "  python: \"\${python_version}\"" \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """

    stub:
    """
    db_root_path="\$(cd "${db_root}" && pwd -P)"
    destination_path="\${db_root_path}/${component}"
    mkdir -p "\${destination_path}"

    case "${component}" in
        taxdump)
            : > "\${destination_path}/names.dmp"
            : > "\${destination_path}/nodes.dmp"
            ;;
        checkm2)
            : > "\${destination_path}/CheckM2_database.dmnd"
            ;;
        eggnog)
            : > "\${destination_path}/eggnog.db"
            : > "\${destination_path}/eggnog_proteins.dmnd"
            ;;
        padloc)
            mkdir -p "\${destination_path}/hmm"
            : > "\${destination_path}/hmm/padlocdb.hmm"
            ;;
    esac

    printf 'component\tstatus\tsource\tdestination\tdetails\n%s\tprepared\tstub\t%s\trequired_paths=stub;source_mode=stub\n' \
        "${component}" \
        "\${destination_path}" \
        > "${component}_report.tsv"

    printf '%s\n' \
        '"${task.process}":' \
        '  python: "stub"' \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """
}
