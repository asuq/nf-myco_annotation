/*
 * Merge paired CheckM2 reports into a single per-sample QC summary.
 */
process ASSIGN_GCODE_AND_QC {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/checkm2" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(checkm2_gcode4_report, name: 'checkm2_gcode4_report.tsv'), path(checkm2_gcode11_report, name: 'checkm2_gcode11_report.tsv')

    output:
    tuple val(meta), path('checkm2_summary.tsv'), emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    summarise_checkm2.py \
        --accession "${meta.accession}" \
        --gcode4-report "${checkm2_gcode4_report}" \
        --gcode11-report "${checkm2_gcode11_report}" \
        --gcode-rule "${params.gcode_rule}" \
        --output checkm2_summary.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/summarise_checkm2.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > checkm2_summary.tsv
    accession	Completeness_gcode4	Completeness_gcode11	Contamination_gcode4	Contamination_gcode11	Coding_Density_gcode4	Coding_Density_gcode11	Average_Gene_Length_gcode4	Average_Gene_Length_gcode11	Total_Coding_Sequences_gcode4	Total_Coding_Sequences_gcode11	Gcode	Low_quality	checkm2_status	warnings
    sample_a	95	82	2	1	0.9	0.8	900	850	800	780	4	false	done
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/summarise_checkm2.py"
    EOF
    '''
}
