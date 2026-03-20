/*
 * Run eggNOG with the required FASTA-header cleanup and emit a comment-free
 * annotation table for downstream inspection.
 */
process EGGNOG {
    tag "${meta.accession}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/eggnog" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(faa), path(eggnog_db)

    output:
    tuple val(meta), path('eggnog'), path('eggnog_annotations.tsv'), path('eggnog.log'), emit: results
    tuple val(meta), path('eggnog_annotations.tsv'), emit: annotations
    path 'versions.yml', emit: versions

    script:
    def prefix = (meta.internal_id ?: meta.accession).toString()
    def extraArgs = (params.eggnog_extra_args ?: '').toString()
    """
    awk '
        /^>/ { gsub(/[[:space:]]+/, "_") }
        { print }
    ' "${faa}" > clean_input.faa

    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    tmp_dir="\$PWD/eggnog_tmp"
    : > eggnog.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> eggnog.log
        rm -rf eggnog "\${tmp_dir}"
        mkdir -p eggnog "\${tmp_dir}"
        set +e
        emapper.py \
            -i clean_input.faa \
            -m diamond \
            --cpu ${task.cpus} \
            --output_dir eggnog \
            -o "${prefix}" \
            --temp_dir "\${tmp_dir}" \
            --report_orthologs \
            --data_dir "${eggnog_db}" \
            --sensmode ultra-sensitive \
            --override \
            ${extraArgs} \
            >> eggnog.log 2>&1
        exit_code=\$?
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_eggnog=%s\n' "\${attempt}" >> eggnog.log
        (( attempt += 1 ))
    done

    ann="eggnog/${prefix}.emapper.annotations"
    if [[ "\${exit_code}" -eq 0 && -f "\${ann}" ]]; then
        grep -v '^##' "\${ann}" | sed -e '1s/^#//' > eggnog_annotations.tsv
    else
        : > eggnog_annotations.tsv
    fi

    rm -rf "\${tmp_dir}"
    printf 'exit_code=%s\n' "\$exit_code" >> eggnog.log

    eggnog_version="\$(python -c \"import importlib.metadata as m; print(m.version('eggnog-mapper'))\" 2>/dev/null || echo NA)"
    printf '"%s":\n  eggnog_mapper: "%s"\n' \
      "${task.process}" \
      "\${eggnog_version}" \
      > versions.yml
    """

    stub:
    '''
    mkdir -p eggnog
    cat <<'EOF' > eggnog_annotations.tsv
    query	seed_ortholog	evalue	score
    sample_a_1	ortholog_1	1e-20	200
    EOF
    : > eggnog.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      eggnog_mapper: "stub"
    EOF
    '''
}
