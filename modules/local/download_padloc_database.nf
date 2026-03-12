/*
 * Reuse or download the PADLOC database into its final directory.
 */
process DOWNLOAD_PADLOC_DATABASE {
    tag "padloc"
    label 'process_medium'
    label 'download_padloc_database'

    input:
    tuple val(destination), val(download_enabled), val(force)

    output:
    tuple val('padloc'), val(destination), path('prep_mode.txt'), val('native_download'), path('lineages.txt'), emit: finalise_input
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination}"
    mode="prepared"

    mkdir -p padloc_bin padloc_bootstrap_data
    cp "\$(command -v padloc)" padloc_bin/padloc.real
    sed -i 's#mkdir -p "\\${SRC_DIR}/../data"#mkdir -p "\\${PADLOC_BOOTSTRAP_DATA}"#' padloc_bin/padloc.real
    chmod +x padloc_bin/padloc.real
    export PADLOC_WRAPPER_REAL="\$PWD/padloc_bin/padloc.real"
    cat <<'EOF' > padloc_bin/padloc
#!/usr/bin/env bash
set -euo pipefail
exec bash "\${PADLOC_WRAPPER_REAL}" "\$@"
EOF
    chmod +x padloc_bin/padloc
    export PADLOC_BOOTSTRAP_DATA="\$PWD/padloc_bootstrap_data"
    export PATH="\$PWD/padloc_bin:\$PATH"

    if [[ -e "\${destination_path}" && ! -d "\${destination_path}" ]]; then
        echo "PADLOC destination must be a directory: \${destination_path}" >&2
        exit 1
    fi

    if [[ -d "\${destination_path}" && -f "\${destination_path}/.nf_myco_ready.json" && -f "\${destination_path}/hmm/padlocdb.hmm" ]]; then
        mode="reuse"
    else
        if [[ -d "\${destination_path}" ]]; then
            if [[ "${force}" != "true" ]]; then
                echo "PADLOC destination exists but is not ready: \${destination_path}" >&2
                exit 1
            fi
            rm -rf "\${destination_path}"
        fi

        if [[ "${download_enabled}" != "true" ]]; then
            echo "PADLOC destination is missing and download is disabled: \${destination_path}" >&2
            exit 1
        fi

        mkdir -p "\${destination_path}"
        padloc --data "\${destination_path}" --db-update
    fi

    printf '%s\n' "\${mode}" > prep_mode.txt
    : > lineages.txt

    {
        printf '"%s":\n' "${task.process}"
        printf '  padloc: "%s"\n' "\$(padloc --version 2>&1 | awk 'NF { print; exit }' || echo NA)"
    } > versions.yml
    """

    stub:
    """
    mkdir -p "${destination}/hmm"
    : > "${destination}/hmm/padlocdb.hmm"
    : > "${destination}/.nf_myco_ready.json"
    printf '%s\n' 'reuse' > prep_mode.txt
    : > lineages.txt
    {
        printf '"%s":\n' "${task.process}"
        printf '  padloc: "%s"\n' 'stub'
    } > versions.yml
    """
}
