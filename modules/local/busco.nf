/*
 * Run BUSCO in genome/offline mode for one sample-lineage pair. Failures are
 * converted into empty JSON summaries so downstream parsing stays sample-level.
 */
process BUSCO {
    tag "${meta.accession} / ${lineage}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/busco/${lineage}" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(genome), val(lineage), path(dataset_dir)

    output:
    tuple val(meta), val(lineage), path("busco_${lineage}"), path('short_summary.json'), path('busco.log'), emit: results
    tuple val(meta), val(lineage), path('short_summary.json'), emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    output_dir="busco_${lineage}"
    busco_download_root="busco_downloads"
    staged_lineage_dir="\${busco_download_root}/lineages/${lineage}"
    dataset_source="\$(cd "${dataset_dir}" && pwd)"

    mkdir -p "\$(dirname "\${staged_lineage_dir}")"
    ln -s "\${dataset_source}" "\${staged_lineage_dir}"

    set +e
    busco \
        --in "${genome}" \
        --download_path "\${busco_download_root}" \
        --lineage_dataset "${lineage}" \
        --mode genome \
        --offline \
        --cpu ${task.cpus} \
        --out "\${output_dir}" \
        --out_path . \
        > busco.log 2>&1
    exit_code=\$?
    set -e

    mkdir -p "\${output_dir}"
    summary_json=\$(find "\${output_dir}" -type f \\( -name 'short_summary*.json' -o -name 'short_summary.json' \\) | head -n 1 || true)
    if [[ -n "\${summary_json}" ]]; then
        cp "\${summary_json}" short_summary.json
    else
        : > short_summary.json
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> busco.log

    busco_version="\$(command -v busco >/dev/null 2>&1 && busco --version 2>&1 | awk 'NF { value=\$0 } END { if (value) print value }' || true)"
    busco_version="\${busco_version:-NA}"
    printf '"%s":\n  busco: "%s"\n' \
      "${task.process}" \
      "\${busco_version}" \
      > versions.yml
    """

    stub:
    """
    mkdir -p "busco_${lineage}"
    cat <<'EOF' > short_summary.json
    {"results": {"C": 98.0, "S": 98.0, "D": 0.0, "F": 1.0, "M": 1.0, "n": 200}}
    EOF
    : > busco.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      busco: "stub"
    EOF
    """
}
