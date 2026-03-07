include { BUILD_MASTER_TABLE } from '../../modules/local/build_master_table'
include { COLLECT_VERSIONS } from '../../modules/local/collect_versions'
include { SELECT_ANI_REPRESENTATIVES } from '../../modules/local/select_ani_representatives'

/*
 * Build the final cohort tables from the keyed per-sample summaries and gather
 * one authoritative provenance report.
 */
workflow FINAL_OUTPUTS {
    take:
    validated_samples
    metadata
    taxonomy
    checkm2_summaries
    sixteen_s_summaries
    busco_summaries
    ccfinder_summaries
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
        .map { meta, lineage, summary -> summary }
        .collect()

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

    final_versions = SELECT_ANI_REPRESENTATIVES.out.versions
        .mix(BUILD_MASTER_TABLE.out.versions)

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
    sample_status = BUILD_MASTER_TABLE.out.sample_status
    ani_representatives = SELECT_ANI_REPRESENTATIVES.out.ani_representatives
    versions_table = COLLECT_VERSIONS.out.versions_table
    versions = final_versions
}
