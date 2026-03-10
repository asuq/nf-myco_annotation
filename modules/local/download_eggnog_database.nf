/*
 * Reuse or download the eggNOG database into its final directory.
 */
process DOWNLOAD_EGGNOG_DATABASE {
    tag "eggnog"
    label 'process_high'

    input:
    tuple val(destination), val(download_enabled), val(force)

    output:
    tuple val('eggnog'), val(destination), path('prep_mode.txt'), val('native_download'), path('lineages.txt'), emit: finalise_input
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination}"
    mode="prepared"

    if [[ -e "\${destination_path}" && ! -d "\${destination_path}" ]]; then
        echo "eggNOG destination must be a directory: \${destination_path}" >&2
        exit 1
    fi

    if [[ -d "\${destination_path}" && -f "\${destination_path}/.nf_myco_ready.json" && -f "\${destination_path}/eggnog.db" && -f "\${destination_path}/eggnog_proteins.dmnd" ]]; then
        mode="reuse"
    else
        if [[ -d "\${destination_path}" ]]; then
            if [[ "${force}" != "true" ]]; then
                echo "eggNOG destination exists but is not ready: \${destination_path}" >&2
                exit 1
            fi
            rm -rf "\${destination_path}"
        fi

        if [[ "${download_enabled}" != "true" ]]; then
            echo "eggNOG destination is missing and download is disabled: \${destination_path}" >&2
            exit 1
        fi

        mkdir -p "\${destination_path}"
        download_eggnog_data.py --data_dir "\${destination_path}" -y ${force ? '-f' : ''}
    fi

    printf '%s\n' "\${mode}" > prep_mode.txt
    : > lineages.txt

    {
        printf '"%s":\n' "${task.process}"
        printf '  eggnog_mapper: "%s"\n' "\$(emapper.py --version 2>&1 | head -n 1 || echo NA)"
    } > versions.yml
    """

    stub:
    """
    mkdir -p "${destination}"
    : > "${destination}/eggnog.db"
    : > "${destination}/eggnog_proteins.dmnd"
    : > "${destination}/.nf_myco_ready.json"
    printf '%s\n' 'reuse' > prep_mode.txt
    : > lineages.txt
    {
        printf '"%s":\n' "${task.process}"
        printf '  eggnog_mapper: "%s"\n' 'stub'
    } > versions.yml
    """
}
