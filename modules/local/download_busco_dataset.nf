/*
 * Download one BUSCO lineage dataset for later offline reuse.
 */
process DOWNLOAD_BUSCO_DATASET {
    tag "${lineage}"
    label 'process_medium'

    input:
    val lineage

    output:
    tuple val(lineage), path('busco_dataset'), emit: dataset
    tuple val(lineage), path('download.log'), emit: log
    path 'versions.yml', emit: versions

    script:
    def downloadRoot = params.busco_download_dir ?: "${params.outdir}/resources/busco"
    """
    mkdir -p "${downloadRoot}"

    set +e
    busco \
        --download_path "${downloadRoot}" \
        --download "${lineage}" \
        > download.log 2>&1
    exit_code=\$?
    set -e

    if [[ \$exit_code -ne 0 ]]; then
        printf 'exit_code=%s\n' "\$exit_code" >> download.log
        exit "\$exit_code"
    fi

    dataset_dir=""
    if [[ -d "${downloadRoot}/${lineage}" ]]; then
        dataset_dir="${downloadRoot}/${lineage}"
    elif [[ -d "${downloadRoot}/lineages/${lineage}" ]]; then
        dataset_dir="${downloadRoot}/lineages/${lineage}"
    fi

    if [[ -z "\${dataset_dir}" ]]; then
        echo "Could not locate downloaded BUSCO dataset directory for ${lineage}." >> download.log
        exit 1
    fi

    ln -s "\${dataset_dir}" busco_dataset

    cat <<EOF > versions.yml
    "${task.process}":
      busco: "$(command -v busco >/dev/null 2>&1 && busco --version 2>&1 | head -n 1 || echo NA)"
    EOF
    """

    stub:
    """
    mkdir -p busco_dataset
    : > busco_dataset/dataset.cfg
    : > download.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      busco: "stub"
    EOF
    """
}
