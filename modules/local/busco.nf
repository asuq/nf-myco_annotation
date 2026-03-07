/*
 * Run BUSCO in genome/offline mode for one sample-lineage pair. Failures are
 * converted into empty JSON summaries so downstream parsing stays sample-level.
 */
process BUSCO {
    tag "${meta.accession} / ${lineage}"

    input:
    tuple val(meta), path(genome), val(lineage), path(dataset_dir)

    output:
    tuple val(meta), val(lineage), path("busco_${lineage}"), path('short_summary.json'), path('busco.log'), emit: results
    tuple val(meta), val(lineage), path('short_summary.json'), emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    output_dir="busco_${lineage}"

    set +e
    busco \
        --in "${genome}" \
        --lineage_dataset "${dataset_dir}" \
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

    cat <<EOF > versions.yml
    "${task.process}":
      busco: "$(command -v busco >/dev/null 2>&1 && busco --version 2>&1 | head -n 1 || echo NA)"
    EOF
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
