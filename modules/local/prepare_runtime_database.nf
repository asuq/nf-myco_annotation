/*
 * Prepare the curated taxdump in its final runtime directory.
 */
process PREP_TAXDUMP_DATABASE {
    tag "taxdump"
    label 'process_medium'
    label 'prep_taxdump_database'
    stageInMode 'symlink'

    input:
    tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(version), path(scratch_root_path), val(force)

    output:
    path('taxdump_report.tsv'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination_parent}/${destination_name}"
    script_path="\$(command -v prepare_runtime_databases.py)"
    python_path="\$(command -v python3)"

    helper_args=(
        --taxdump-dest
        "\${destination_path}"
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

    helper_args+=(--scratch-root "${scratch_root_path}")

    "\${python_path}" "\${script_path}" "\${helper_args[@]}"

    python_version="\$("\${python_path}" --version 2>&1 | sed 's/^Python //')"
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

/*
 * Prepare the curated Codetta profile bundle in its final runtime directory.
 */
process PREP_CODETTA_DATABASE {
    tag "codetta"
    label 'process_medium'
    label 'prep_codetta_database'
    stageInMode 'symlink'

    input:
    tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), path(scratch_root_path), val(force)

    output:
    path('codetta_report.tsv'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination_parent}/${destination_name}"
    script_path="\$(command -v prepare_runtime_databases.py)"
    python_path="\$(command -v python3)"

    helper_args=(
        --codetta-dest
        "\${destination_path}"
        --report
        codetta_report.tsv
    )

    if [[ "${download_enabled}" == "true" ]]; then
        helper_args+=(--download)
    fi

    if [[ "${force}" == "true" ]]; then
        helper_args+=(--force)
    fi

    helper_args+=(--scratch-root "${scratch_root_path}")

    "\${python_path}" "\${script_path}" "\${helper_args[@]}"

    python_version="\$("\${python_path}" --version 2>&1 | sed 's/^Python //')"
    printf '%s\n' \
        '"${task.process}":' \
        "  python: \"\${python_version}\"" \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """

    stub:
    """
    mkdir -p "${destination}"
    : > "${destination}/Pfam-A_enone.hmm"
    : > "${destination}/Pfam-A_enone.hmm.h3f"
    : > "${destination}/Pfam-A_enone.hmm.h3i"
    : > "${destination}/Pfam-A_enone.hmm.h3m"
    : > "${destination}/Pfam-A_enone.hmm.h3p"
    : > "${destination}/.nf_myco_ready.json"

    {
        printf 'component\tstatus\tsource\tdestination\tdetails\n'
        printf 'codetta\tprepared\tstub\t%s\trequired_paths=Pfam-A_enone.hmm,Pfam-A_enone.hmm.h3f,Pfam-A_enone.hmm.h3i,Pfam-A_enone.hmm.h3m,Pfam-A_enone.hmm.h3p;source_mode=stub\n' "${destination}"
    } > codetta_report.tsv

    printf '%s\n' \
        '"${task.process}":' \
        '  python: "stub"' \
        '  helper: "bin/prepare_runtime_databases.py"' \
        > versions.yml
    """
}
