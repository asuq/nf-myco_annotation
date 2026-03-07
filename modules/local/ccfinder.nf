/*
 * Run CRISPRCasFinder through a user-configurable command template so the
 * legacy invocation can be preserved without hard-coding the wrong flags here.
 */
process CCFINDER {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(genome), val(gcode)

    output:
    tuple val(meta), path('ccfinder'), path('result.json'), path('ccfinder.log'), emit: results
    tuple val(meta), path('result.json'), emit: result_json
    path 'versions.yml', emit: versions

    script:
    def commandTemplate = (params.ccfinder_command_template ?: '').replace("'", "'\"'\"'")
    """
    if [[ -z '${commandTemplate}' ]]; then
        echo "params.ccfinder_command_template is required for CCFINDER." >&2
        exit 1
    fi

    genome_path=\$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "${genome}")
    outdir_path="\$PWD/ccfinder"
    genome_escaped=\$(printf '%q' "\${genome_path}")
    outdir_escaped=\$(printf '%q' "\${outdir_path}")
    gcode_escaped=\$(printf '%q' "${gcode}")
    command_template='${commandTemplate}'
    command_template="\${command_template//\\{genome\\}/\${genome_escaped}}"
    command_template="\${command_template//\\{outdir\\}/\${outdir_escaped}}"
    command_template="\${command_template//\\{gcode\\}/\${gcode_escaped}}"

    set +e
    eval "\${command_template}" > ccfinder.log 2>&1
    exit_code=\$?
    set -e

    mkdir -p ccfinder
    if [[ -f ccfinder/result.json ]]; then
        cp ccfinder/result.json result.json
    else
        : > result.json
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> ccfinder.log

    cat <<EOF > versions.yml
    "${task.process}":
      invocation_mode: "template"
    EOF
    """

    stub:
    '''
    mkdir -p ccfinder
    cat <<'EOF' > result.json
    {"Sequences": []}
    EOF
    cp result.json ccfinder/result.json
    : > ccfinder.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      invocation_mode: "stub"
    EOF
    '''
}
