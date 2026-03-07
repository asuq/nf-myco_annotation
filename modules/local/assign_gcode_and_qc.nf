/*
 * Placeholder for CheckM2 summarisation, gcode assignment, and QC status.
 */
process ASSIGN_GCODE_AND_QC {
    input:
    path checkm2_gcode4
    path checkm2_gcode11

    output:
    path 'gcode_qc.tsv', emit: summary
    path 'versions.yml', emit: versions

    stub:
    '''
    : > gcode_qc.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''

    script:
    '''
    echo "ASSIGN_GCODE_AND_QC is a placeholder module." >&2
    exit 1
    '''
}
