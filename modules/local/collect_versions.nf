/*
 * Merge per-process versions files with runtime/container/resource context into
 * one final TSV report.
 */
process COLLECT_VERSIONS {
    tag "versions"
    label 'process_single'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'tool_and_db_versions.tsv' ? filename : null },
    )

    input:
    path version_files, name: 'version_files/versions??.yml'
    val nextflow_version
    val pipeline_version
    val git_commit
    val container_engine

    output:
    path 'tool_and_db_versions.tsv', emit: versions_table

    script:
    def versionFileList = version_files instanceof Collection ? version_files : [version_files]
    def versionArgs = versionFileList.collect { "--version-file \"${it}\"" }.join(' \\\n        ')
    def lineageArgs = (params.busco_lineages as List<String>).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n        ')
    def containerRefs = [
        'python'  : params.python_container ?: 'NA',
        'seqtk'   : params.seqtk_container ?: 'NA',
        'barrnap' : params.barrnap_container ?: 'NA',
        'checkm2' : params.checkm2_container ?: 'NA',
        'busco'   : params.busco_container ?: 'NA',
        'prokka'  : params.prokka_container ?: 'NA',
        'ccfinder': params.ccfinder_container ?: 'NA',
        'padloc'  : params.padloc_container ?: 'NA',
        'eggnog'  : params.eggnog_container ?: 'NA',
        'fastani' : params.fastani_container ?: 'NA',
    ]
    def containerArgs = containerRefs.collect { name, value ->
        "--container-ref \"${name}=${value}\""
    }.join(' \\\n        ')
    """
    script_path="\$(command -v collect_versions.py)"
    python3 "\${script_path}" \
        ${versionArgs} \
        --nextflow-version "${nextflow_version}" \
        --pipeline-version "${pipeline_version}" \
        --git-commit "${git_commit}" \
        --container-engine "${container_engine}" \
        --use-biocontainers "${params.use_biocontainers}" \
        --checkm2-db "${params.checkm2_db ?: 'NA'}" \
        --checkm2-db-label "${params.checkm2_db_label ?: 'NA'}" \
        --taxdump "${params.taxdump ?: 'NA'}" \
        --taxdump-label "${params.taxdump_label ?: 'NA'}" \
        ${lineageArgs} \
        --busco-download-dir "${params.busco_download_dir ?: 'NA'}" \
        --eggnog-db "${params.eggnog_db ?: 'NA'}" \
        --eggnog-db-label "${params.eggnog_db_label ?: 'NA'}" \
        --padloc-db-label "${params.padloc_db_label ?: 'NA'}" \
        ${containerArgs} \
        --output tool_and_db_versions.tsv
    """

    stub:
    '''
    cat <<'EOF' > tool_and_db_versions.tsv
    component	kind	version	image_or_path	notes
    nextflow	runtime	stub	NA	workflow
    python	runtime	stub	NA	reported by VALIDATE_INPUTS
    EOF
    '''
}
