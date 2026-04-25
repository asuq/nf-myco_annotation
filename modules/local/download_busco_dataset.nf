/*
 * Download one BUSCO lineage dataset for later offline reuse.
 */
process DOWNLOAD_BUSCO_DATASET {
    tag "${lineage}"
    label 'process_single'

    input:
    tuple val(lineage), path(download_parent), val(download_name)

    output:
    tuple val(lineage), path("${lineage}"), emit: dataset
    tuple val(lineage), path('download.log'), emit: log
    path 'versions.yml', emit: versions

    script:
    """
    download_root="${download_parent}/${download_name}"
    mkdir -p "\${download_root}"

    set +e
    busco \
        --download_path "\${download_root}" \
        --download "${lineage}" \
        > download.log 2>&1
    exit_code=\$?
    set -e

    if [[ \$exit_code -ne 0 ]]; then
        printf 'exit_code=%s\n' "\$exit_code" >> download.log
        exit "\$exit_code"
    fi

    dataset_dir=""
    if [[ -d "\${download_root}/${lineage}" ]]; then
        dataset_dir="\${download_root}/${lineage}"
    elif [[ -d "\${download_root}/lineages/${lineage}" ]]; then
        dataset_dir="\${download_root}/lineages/${lineage}"
    fi

    if [[ -z "\${dataset_dir}" ]]; then
        echo "Could not locate downloaded BUSCO dataset directory for ${lineage}." >> download.log
        exit 1
    fi

    ln -s "\${dataset_dir}" "${lineage}"

    {
        printf '"%s":\n' "${task.process}"
        printf '  busco: "%s"\n' "\$(command -v busco >/dev/null 2>&1 && busco --version 2>&1 | head -n 1 || echo NA)"
    } > versions.yml
    """

    stub:
    """
    mkdir -p "${lineage}"
    : > "${lineage}/dataset.cfg"
    : > download.log
    {
        printf '"%s":\n' "${task.process}"
        printf '  busco: "%s"\n' 'stub'
    } > versions.yml
    """
}
