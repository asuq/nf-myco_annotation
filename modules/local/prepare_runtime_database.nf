/*
 * Prepare the curated taxdump in its final runtime directory.
 */
process PREP_TAXDUMP_DATABASE {
    tag "taxdump"
    label 'process_medium'

    input:
    tuple val(destination), val(download_enabled), val(version), val(scratch_root), val(force)

    output:
    path('taxdump_report.tsv'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    script_path="${projectDir}/bin/prepare_runtime_databases.py"

    helper_args=(
        --taxdump-dest
        "${destination}"
        --report
        taxdump_report.tsv
    )

    if [[ "${download_enabled}" == "true" ]]; then
        helper_args+=(--download)
    fi

    if [[ "${force}" == "true" ]]; then
        helper_args+=(--force)
    fi

    if [[ -n "${version ?: ''}" ]]; then
        helper_args+=(--taxdump-version "${version}")
    fi

    if [[ -n "${scratch_root ?: ''}" ]]; then
        helper_args+=(--scratch-root "${scratch_root}")
    fi

    "\${script_path}" "\${helper_args[@]}"

    python_version="\$(/usr/local/bin/python3 --version 2>&1 | sed 's/^Python //')"
    printf '%s\n' \
        '"${task.process}":' \
        "  python: \"\${python_version}\"" \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """

    stub:
    """
    mkdir -p "${destination}"
    : > "${destination}/names.dmp"
    : > "${destination}/nodes.dmp"

    {
        printf 'component\tstatus\tsource\tdestination\tdetails\n'
        printf 'taxdump\tprepared\tstub\t%s\trequired_paths=names.dmp,nodes.dmp;source_mode=stub\n' "${destination}"
    } > taxdump_report.tsv

    printf '%s\n' \
        '"${task.process}":' \
        '  python: "stub"' \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """
}
