/*
 * Run eggNOG with the required FASTA-header cleanup and emit a comment-free
 * annotation table for downstream inspection.
 */
process EGGNOG {
    tag "${meta.accession}"
    label 'process_high'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/eggnog" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path('eggnog'), path('eggnog_annotations.tsv'), path('eggnog.log'), emit: results
    tuple val(meta), path('eggnog_annotations.tsv'), emit: annotations
    path 'versions.yml', emit: versions

    script:
    def eggnogDb = (params.eggnog_db ?: '').toString()
    def prefix = (meta.internal_id ?: meta.accession).toString()
    def extraArgs = (params.eggnog_extra_args ?: '').toString()
    """
    if [[ -z '${eggnogDb}' ]]; then
        echo "params.eggnog_db is required for EGGNOG." >&2
        exit 1
    fi

    awk '
        /^>/ { gsub(/[[:space:]]+/, "_") }
        { print }
    ' "${faa}" > clean_input.faa

    tmp_dir="\$PWD/eggnog_tmp"
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
        --data_dir "${eggnogDb}" \
        --sensmode ultra-sensitive \
        --override \
        ${extraArgs} \
        > eggnog.log 2>&1
    exit_code=\$?
    set -e

    ann="eggnog/${prefix}.emapper.annotations"
    if [[ -f "\${ann}" ]]; then
        grep -v '^##' "\${ann}" | sed -e '1s/^#//' > eggnog_annotations.tsv
    else
        : > eggnog_annotations.tsv
    fi

    rm -rf "\${tmp_dir}"
    printf 'exit_code=%s\n' "\$exit_code" >> eggnog.log

    cat <<EOF > versions.yml
    "${task.process}":
      eggnog_mapper: "$(emapper.py --version 2>&1 | head -n 1 || echo NA)"
    EOF
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
