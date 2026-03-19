/*
 * Reuse or download the CheckM2 database into its final directory.
 */
process DOWNLOAD_CHECKM2_DATABASE {
    tag "checkm2"
    label 'process_medium'
    label 'download_checkm2_database'
    stageInMode 'symlink'

    input:
    tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(force)

    output:
    tuple val('checkm2'), val(destination), path(destination_parent), val(destination_name), path('prep_mode.txt'), val('native_download'), path('lineages.txt'), emit: finalise_input
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination_parent}/${destination_name}"
    mode="prepared"

    count_top_level_dmnd() {
        shopt -s nullglob
        dmnd_candidates=("\$1"/*.dmnd)
        shopt -u nullglob
        printf '%s\n' "\${#dmnd_candidates[@]}"
    }

    normalise_download_layout() {
        nested_root="\$1/CheckM2_database"

        if [[ ! -d "\${nested_root}" ]]; then
            return
        fi

        if [[ "\$(count_top_level_dmnd "\$1")" != "0" ]]; then
            return
        fi

        shopt -s dotglob nullglob
        nested_entries=("\${nested_root}"/*)
        shopt -u dotglob nullglob

        if [[ \${#nested_entries[@]} -eq 0 ]]; then
            rmdir "\${nested_root}"
            return
        fi

        for entry in "\${nested_entries[@]}"; do
            mv "\${entry}" "\$1/"
        done
        rmdir "\${nested_root}"
    }

    if [[ -e "\${destination_path}" && ! -d "\${destination_path}" ]]; then
        echo "CheckM2 destination must be a directory: \${destination_path}" >&2
        exit 1
    fi

    if [[ -d "\${destination_path}" ]]; then
        normalise_download_layout "\${destination_path}"
    fi

    is_ready=false
    if [[ -d "\${destination_path}" && -f "\${destination_path}/.nf_myco_ready.json" ]]; then
        if [[ "\$(count_top_level_dmnd "\${destination_path}")" == "1" ]]; then
            is_ready=true
        fi
    fi

    if [[ "\${is_ready}" == "true" ]]; then
        mode="reuse"
    else
        if [[ -d "\${destination_path}" && "\$(count_top_level_dmnd "\${destination_path}")" == "1" ]]; then
            mode="prepared"
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
            normalise_download_layout "\${destination_path}"
        fi
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
