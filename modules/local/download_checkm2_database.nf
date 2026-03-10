/*
 * Reuse or download the CheckM2 database into its final directory.
 */
process DOWNLOAD_CHECKM2_DATABASE {
    tag "checkm2"
    label 'process_medium'

    input:
    tuple val(destination), val(download_enabled), val(force)

    output:
    tuple val('checkm2'), val(destination), path('prep_mode.txt'), val('native_download'), path('lineages.txt'), emit: finalise_input
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination}"
    mode="prepared"

    if [[ -e "\${destination_path}" && ! -d "\${destination_path}" ]]; then
        echo "CheckM2 destination must be a directory: \${destination_path}" >&2
        exit 1
    fi

    is_ready=false
    if [[ -d "\${destination_path}" && -f "\${destination_path}/.nf_myco_ready.json" ]]; then
        shopt -s nullglob
        dmnd_candidates=("\${destination_path}"/*.dmnd)
        shopt -u nullglob
        if [[ \${#dmnd_candidates[@]} -eq 1 ]]; then
            is_ready=true
        fi
    fi

    if [[ "\${is_ready}" == "true" ]]; then
        mode="reuse"
    else
        if [[ -d "\${destination_path}" ]]; then
            if [[ "${force}" != "true" ]]; then
                echo "CheckM2 destination exists but is not ready: \${destination_path}" >&2
                exit 1
            fi
            rm -rf "\${destination_path}"
        fi

        if [[ "${download_enabled}" != "true" ]]; then
            echo "CheckM2 destination is missing and download is disabled: \${destination_path}" >&2
            exit 1
        fi

        mkdir -p "\${destination_path}"
        checkm2 database --download --path "\${destination_path}" --no_write_json_db
    fi

    printf '%s\n' "\${mode}" > prep_mode.txt
    : > lineages.txt

    {
        printf '"%s":\n' "${task.process}"
        printf '  checkm2: "%s"\n' "\$(checkm2 --version 2>&1 | head -n 1 || echo NA)"
    } > versions.yml
    """

    stub:
    """
    mkdir -p "${destination}"
    : > "${destination}/CheckM2_database.dmnd"
    : > "${destination}/.nf_myco_ready.json"
    printf '%s\n' 'reuse' > prep_mode.txt
    : > lineages.txt
    {
        printf '"%s":\n' "${task.process}"
        printf '  checkm2: "%s"\n' 'stub'
    } > versions.yml
    """
}
