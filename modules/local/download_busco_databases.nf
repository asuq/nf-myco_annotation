/*
 * Reuse or download the configured BUSCO lineage datasets into one root.
 */
process DOWNLOAD_BUSCO_DATABASES {
    tag "busco"
    label 'process_medium'
    label 'download_busco_databases'

    input:
    tuple val(destination), val(download_enabled), val(lineages), val(force)

    output:
    tuple val('busco_root'), val(destination), path('prep_mode.txt'), val('native_download'), path('lineages.txt'), emit: finalise_input
    path 'versions.yml', emit: versions

    script:
    def renderedLineagesShell = (lineages as List<String>).collect { value ->
        "'${value.replace("'", "'\"'\"'")}'"
    }.join(' ')
    """
    destination_path="${destination}"
    mode="prepared"
    lineages=(${renderedLineagesShell})

    printf '%s\n' "\${lineages[@]}" > lineages.txt

    if [[ -e "\${destination_path}" && ! -d "\${destination_path}" ]]; then
        echo "BUSCO destination must be a directory: \${destination_path}" >&2
        exit 1
    fi

    is_ready=true
    if [[ ! -d "\${destination_path}" || ! -f "\${destination_path}/.nf_myco_ready.json" ]]; then
        is_ready=false
    else
        while IFS= read -r lineage; do
            [[ -z "\${lineage}" ]] && continue
            if [[ ! -f "\${destination_path}/\${lineage}/dataset.cfg" ]]; then
                is_ready=false
                break
            fi
        done < lineages.txt
    fi

    if [[ "\${is_ready}" == "true" ]]; then
        mode="reuse"
    else
        if [[ -d "\${destination_path}" ]]; then
            if [[ "${force}" != "true" ]]; then
                echo "BUSCO destination exists but is not ready: \${destination_path}" >&2
                exit 1
            fi
            rm -rf "\${destination_path}"
        fi

        if [[ "${download_enabled}" != "true" ]]; then
            echo "BUSCO destination is missing and download is disabled: \${destination_path}" >&2
            exit 1
        fi

        mkdir -p "\${destination_path}"
        while IFS= read -r lineage; do
            [[ -z "\${lineage}" ]] && continue
            busco --download_path "\${destination_path}" --download "\${lineage}"
            if [[ -d "\${destination_path}/lineages/\${lineage}" && ! -d "\${destination_path}/\${lineage}" ]]; then
                mv "\${destination_path}/lineages/\${lineage}" "\${destination_path}/\${lineage}"
            fi
        done < lineages.txt
    fi

    printf '%s\n' "\${mode}" > prep_mode.txt

    {
        printf '"%s":\n' "${task.process}"
        printf '  busco: "%s"\n' "\$(busco --version 2>&1 | head -n 1 || echo NA)"
    } > versions.yml
    """

    stub:
    def stubRenderedLineagesShell = (lineages as List<String>).collect { value ->
        "'${value.replace("'", "'\"'\"'")}'"
    }.join(' ')
    """
    mkdir -p "${destination}"
    lineages=(${stubRenderedLineagesShell})
    printf '%s\n' "\${lineages[@]}" > lineages.txt
    for lineage in "\${lineages[@]}"; do
        mkdir -p "${destination}/\${lineage}"
        : > "${destination}/\${lineage}/dataset.cfg"
    done
    : > "${destination}/.nf_myco_ready.json"
    printf '%s\n' 'reuse' > prep_mode.txt
    {
        printf '"%s":\n' "${task.process}"
        printf '  busco: "%s"\n' 'stub'
    } > versions.yml
    """
}
