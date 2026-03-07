/*
 * Run eggNOG with the required FASTA-header cleanup and emit a comment-free
 * annotation table for downstream inspection.
 */
process EGGNOG {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path('eggnog'), path('eggnog_annotations.tsv'), path('eggnog.log'), emit: results
    tuple val(meta), path('eggnog_annotations.tsv'), emit: annotations
    path 'versions.yml', emit: versions

    script:
    def defaultTemplate = params.eggnog_db ? 'emapper.py --data_dir {db_dir} --cpu {cpus} -i {faa} --itype proteins --output_dir {outdir} -o eggnog' : ''
    def commandTemplate = ((params.eggnog_command_template ?: defaultTemplate) as String).replace("'", "'\"'\"'")
    def eggnogDb = (params.eggnog_db ?: '').toString().replace("'", "'\"'\"'")
    def invocationMode = params.eggnog_command_template ? 'template' : 'default'
    """
    if [[ -z '${commandTemplate}' ]]; then
        echo "Either params.eggnog_command_template or params.eggnog_db is required for EGGNOG." >&2
        exit 1
    fi

    awk '
        /^>/ { gsub(/[[:space:]]+/, "_") }
        { print }
    ' "${faa}" > eggnog_input.faa

    faa_path="\$PWD/eggnog_input.faa"
    outdir_path="\$PWD/eggnog"
    prefix_path="eggnog"
    faa_escaped=\$(printf '%q' "\${faa_path}")
    outdir_escaped=\$(printf '%q' "\${outdir_path}")
    prefix_escaped=\$(printf '%q' "\${prefix_path}")
    cpus_escaped=\$(printf '%q' "${task.cpus}")
    db_dir_escaped=\$(printf '%q' '${eggnogDb}')
    command_template='${commandTemplate}'
    command_template="\${command_template//\\{faa\\}/\${faa_escaped}}"
    command_template="\${command_template//\\{outdir\\}/\${outdir_escaped}}"
    command_template="\${command_template//\\{prefix\\}/\${prefix_escaped}}"
    command_template="\${command_template//\\{cpus\\}/\${cpus_escaped}}"
    command_template="\${command_template//\\{db_dir\\}/\${db_dir_escaped}}"

    set +e
    eval "\${command_template}" > eggnog.log 2>&1
    exit_code=\$?
    set -e

    mkdir -p eggnog
    annotations_file=\$(find eggnog -type f \\( -name '*.annotations' -o -name '*.emapper.annotations' -o -name '*.tsv' \\) | head -n 1 || true)
    if [[ -n "\${annotations_file}" ]]; then
        grep -v '^#' "\${annotations_file}" > eggnog_annotations.tsv || true
    else
        : > eggnog_annotations.tsv
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> eggnog.log

    cat <<EOF > versions.yml
    "${task.process}":
      invocation_mode: "${invocationMode}"
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
      invocation_mode: "stub"
    EOF
    '''
}
