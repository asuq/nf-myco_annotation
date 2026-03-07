/*
 * Run PADLOC via a user-configurable command template to avoid hard-coding the
 * wrong invocation while still producing a real module wrapper.
 */
process PADLOC {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(clean_gff), path(faa)

    output:
    tuple val(meta), path('padloc'), path('padloc.log'), emit: results
    path 'versions.yml', emit: versions

    script:
    def commandTemplate = (params.padloc_command_template ?: '').replace("'", "'\"'\"'")
    """
    if [[ -z '${commandTemplate}' ]]; then
        echo "params.padloc_command_template is required for PADLOC." >&2
        exit 1
    fi

    gff_path=\$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "${clean_gff}")
    faa_path=\$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "${faa}")
    outdir_path="\$PWD/padloc"
    gff_escaped=\$(printf '%q' "\${gff_path}")
    faa_escaped=\$(printf '%q' "\${faa_path}")
    outdir_escaped=\$(printf '%q' "\${outdir_path}")
    command_template='${commandTemplate}'
    command_template="\${command_template//\\{gff\\}/\${gff_escaped}}"
    command_template="\${command_template//\\{faa\\}/\${faa_escaped}}"
    command_template="\${command_template//\\{outdir\\}/\${outdir_escaped}}"

    set +e
    eval "\${command_template}" > padloc.log 2>&1
    exit_code=\$?
    set -e

    mkdir -p padloc
    printf 'exit_code=%s\n' "\$exit_code" >> padloc.log

    cat <<EOF > versions.yml
    "${task.process}":
      invocation_mode: "template"
    EOF
    """

    stub:
    '''
    mkdir -p padloc
    : > padloc/placeholder.txt
    : > padloc.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      invocation_mode: "stub"
    EOF
    '''
}
