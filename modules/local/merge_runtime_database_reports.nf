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

    output:
    path 'runtime_database_report.tsv', emit: report
    path 'nextflow_args.txt', emit: args
    path 'versions.yml', emit: versions

    script:
    """
    script_path="${projectDir}/bin/merge_runtime_database_reports.py"
    report_args=()
    for report_file in reports/report*.tsv; do
        report_args+=(--report "\${report_file}")
    done
    "\${script_path}" \
        --output-report runtime_database_report.tsv \
        --output-args nextflow_args.txt \
        "\${report_args[@]}"

    {
        printf '"%s":\n' "${task.process}"
        printf '  python: "%s"\n' "\$(/usr/local/bin/python3 --version 2>&1 | sed 's/^Python //')"
        printf '  helper: "%s"\n' 'bin/merge_runtime_database_reports.py'
    } > versions.yml
    """

    stub:
    """
    script_path="${projectDir}/bin/merge_runtime_database_reports.py"
    report_args=()
    for report_file in reports/report*.tsv; do
        report_args+=(--report "\${report_file}")
    done
    "\${script_path}" \
        --output-report runtime_database_report.tsv \
        --output-args nextflow_args.txt \
        "\${report_args[@]}"

    {
        printf '"%s":\n' "${task.process}"
        printf '  python: "%s"\n' 'stub'
        printf '  helper: "%s"\n' 'bin/merge_runtime_database_reports.py'
    } > versions.yml
    """
}
