/*
 * Merge runtime database reports and render a reusable Nextflow argument file.
 */
process MERGE_RUNTIME_DATABASE_REPORTS {
    tag "runtime-database-report"
    label 'process_single'
    publishDir(
        "${params.outdir}",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename in ['runtime_database_report.tsv', 'nextflow_args.txt'] ? filename : null },
    )

    input:
    path reports, name: 'reports/report??.tsv'
    val db_root

    output:
    path 'runtime_database_report.tsv', emit: report
    path 'nextflow_args.txt', emit: args
    path 'versions.yml', emit: versions

    script:
    """
    script_path="\$(command -v merge_runtime_database_reports.py)"
    report_args=()
    for report_file in reports/report*.tsv; do
        report_args+=(--report "\${report_file}")
    done
    python3 "\${script_path}" \
        --db-root "${db_root}" \
        --output-report runtime_database_report.tsv \
        --output-args nextflow_args.txt \
        "\${report_args[@]}"

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      helper: "bin/merge_runtime_database_reports.py"
    EOF
    """

    stub:
    """
    cat <<'EOF' > runtime_database_report.tsv
    component	status	source	destination	details
    taxdump	prepared	stub	${db_root}/taxdump	required_paths=names.dmp,nodes.dmp;source_mode=stub
    checkm2	prepared	stub	${db_root}/checkm2	required_paths=CheckM2_database.dmnd;source_mode=stub
    busco_root	prepared	derived	${db_root}/busco	required_paths=bacillota_odb12/dataset.cfg,mycoplasmatota_odb12/dataset.cfg;source_mode=stub
    eggnog	prepared	stub	${db_root}/eggnog	required_paths=eggnog.db,eggnog_proteins.dmnd;source_mode=stub
    padloc	prepared	stub	${db_root}/padloc	required_paths=hmm/padlocdb.hmm;source_mode=stub
    EOF

    cat <<'EOF' > nextflow_args.txt
    nextflow run . \\
      --taxdump ${db_root}/taxdump \\
      --checkm2_db ${db_root}/checkm2 \\
      --busco_download_dir ${db_root}/busco \\
      --eggnog_db ${db_root}/eggnog \\
      --padloc_db ${db_root}/padloc
    EOF

    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      helper: "bin/merge_runtime_database_reports.py"
    EOF
    """
}
