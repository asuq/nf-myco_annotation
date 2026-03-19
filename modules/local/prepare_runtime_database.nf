/*
 * Prepare the curated taxdump in its final runtime directory.
 */
process PREP_TAXDUMP_DATABASE {
    tag "taxdump"
    label 'process_medium'
    label 'prep_taxdump_database'

    input:
    tuple val(destination), val(download_enabled), val(version), val(scratch_root), val(force)

    output:
    path('taxdump_report.tsv'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    script_path="/usr/local/bin/prepare_runtime_databases.py"

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

    /usr/local/bin/python3 "\${script_path}" "\${helper_args[@]}"

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

/*
 * Prepare the curated Codetta profile bundle in its final runtime directory.
 */
process PREP_CODETTA_DATABASE {
    tag "codetta"
    label 'process_medium'
    label 'prep_codetta_database'

    input:
    tuple val(destination), val(download_enabled), val(scratch_root), val(force)

    output:
    path('codetta_report.tsv'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    script_path="/usr/local/bin/prepare_runtime_databases.py"

    helper_args=(
        --codetta-dest
        "${destination}"
        --report
        codetta_report.tsv
    )

    if [[ "${download_enabled}" == "true" ]]; then
        helper_args+=(--download)
    fi

    if [[ "${force}" == "true" ]]; then
        helper_args+=(--force)
    fi

    if [[ -n "${scratch_root ?: ''}" ]]; then
        helper_args+=(--scratch-root "${scratch_root}")
    fi

    /usr/local/bin/python3 "\${script_path}" "\${helper_args[@]}"

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
