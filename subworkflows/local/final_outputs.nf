include { BUILD_MASTER_TABLE } from '../../modules/local/build_master_table'
include { COLLECT_VERSIONS } from '../../modules/local/collect_versions'
include { SELECT_ANI_REPRESENTATIVES } from '../../modules/local/select_ani_representatives'
include { WRITE_SAMPLE_STATUS } from '../../modules/local/write_sample_status'

/*
 * Build the final cohort tables from the keyed per-sample summaries and gather
 * one authoritative provenance report.
 */
workflow FINAL_OUTPUTS {
    take:
    validated_samples
    initial_status
    metadata
    taxonomy
    checkm2_summaries
    sixteen_s_summaries
    busco_summaries
    ccfinder_summaries
    prokka_results
    padloc_results
    eggnog_results
    ani_clusters
    ani_metadata
    ani_matrix
    version_files
    nextflow_version
    pipeline_version
    git_commit
    container_engine

    main:
    checkm2_seed = Channel.value(file("${projectDir}/assets/testdata/headers/checkm2_summary.tsv"))
    sixteen_s_seed = Channel.value(file("${projectDir}/assets/testdata/headers/16s_status.tsv"))
    ccfinder_seed = Channel.value(file("${projectDir}/assets/testdata/headers/ccfinder_strains.tsv"))
    append_columns = Channel.value(file("${projectDir}/assets/master_table_append_columns.txt"))
    sample_status_columns = Channel.value(file("${projectDir}/assets/sample_status_columns.txt"))
    primaryBuscoColumn = params.busco_primary_column ?: "BUSCO_${params.busco_lineages[0]}"

    extractExitCode = { logFile ->
        def file = logFile.toFile()
        if (!file.exists()) {
            return 'NA'
        }
        def matches = file.readLines().findAll { line -> line.startsWith('exit_code=') }
        return matches ? matches[-1].split('=', 2)[1].trim() : 'NA'
    }

    countTopLevelFiles = { directoryPath ->
        def files = directoryPath.toFile().listFiles()
        return files == null ? 0 : files.count { file -> file.isFile() }
    }

    combined_checkm2 = checkm2_seed
        .mix(checkm2_summaries.map { meta, summary -> summary })
        .collectFile(
            name: 'checkm2_summaries.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_16s = sixteen_s_seed
        .mix(sixteen_s_summaries.map { meta, best16s, status -> status })
        .collectFile(
            name: '16s_statuses.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_ccfinder = ccfinder_seed
        .mix(ccfinder_summaries.map { meta, strains, contigs, crisprs -> strains })
        .collectFile(
            name: 'ccfinder_strains.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    collected_busco = busco_summaries

    prokkaManifest = Channel
        .of('accession\texit_code\tgff_size\tfaa_size')
        .concat(
            prokka_results.map { meta, prokkaDir, gff, faa, log ->
                "${meta.accession}\t${extractExitCode.call(log)}\t${gff.toFile().length()}\t${faa.toFile().length()}"
            }
        )
        .collectFile(name: 'prokka_manifest.tsv', newLine: true, sort: false)

    padlocManifest = Channel
        .of('accession\texit_code\tresult_file_count')
        .concat(
            padloc_results.map { meta, padlocDir, log ->
                "${meta.accession}\t${extractExitCode.call(log)}\t${countTopLevelFiles.call(padlocDir)}"
            }
        )
        .collectFile(name: 'padloc_manifest.tsv', newLine: true, sort: false)

    eggnogManifest = Channel
        .of('accession\texit_code\tannotations_size\tresult_file_count')
        .concat(
            eggnog_results.map { meta, eggnogDir, annotations, log ->
                "${meta.accession}\t${extractExitCode.call(log)}\t${annotations.toFile().length()}\t${countTopLevelFiles.call(eggnogDir)}"
            }
        )
        .collectFile(name: 'eggnog_manifest.tsv', newLine: true, sort: false)

    SELECT_ANI_REPRESENTATIVES(
        ani_clusters,
        ani_metadata,
        ani_matrix,
    )

    BUILD_MASTER_TABLE(
        validated_samples,
        metadata,
        append_columns,
        taxonomy,
        combined_checkm2,
        combined_16s,
        collected_busco,
        combined_ccfinder,
        SELECT_ANI_REPRESENTATIVES.out.ani_summary,
    )

    WRITE_SAMPLE_STATUS(
        validated_samples,
        initial_status,
        sample_status_columns,
        metadata,
        taxonomy,
        combined_checkm2,
        combined_16s,
        collected_busco,
        combined_ccfinder,
        prokkaManifest,
        padlocManifest,
        eggnogManifest,
        SELECT_ANI_REPRESENTATIVES.out.ani_summary,
        primaryBuscoColumn,
    )

    final_versions = SELECT_ANI_REPRESENTATIVES.out.versions
        .mix(BUILD_MASTER_TABLE.out.versions)
        .mix(WRITE_SAMPLE_STATUS.out.versions)

    collected_versions = version_files
        .mix(final_versions)
        .collect()

    COLLECT_VERSIONS(
        collected_versions,
        nextflow_version,
        pipeline_version,
        git_commit,
        container_engine,
    )

    emit:
    master_table = BUILD_MASTER_TABLE.out.master_table
    sample_status = WRITE_SAMPLE_STATUS.out.sample_status
    ani_representatives = SELECT_ANI_REPRESENTATIVES.out.ani_representatives
    versions_table = COLLECT_VERSIONS.out.versions_table
    versions = final_versions
}
