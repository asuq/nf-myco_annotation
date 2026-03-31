include { BUILD_COHORT_16S } from '../../modules/local/build_cohort_16s'

/*
 * Build the cohort-level intact and partial 16S FASTA sets from the per-sample
 * 16S summaries emitted by the QC branch.
 */
workflow COHORT_16S {
    take:
    sixteen_s_summaries
    metadata

    main:
    collected_best_16s = sixteen_s_summaries
        .map { meta, best16s, status -> best16s }
        .collect()

    collected_status_tables = sixteen_s_summaries
        .map { meta, best16s, status -> status }
        .collect()

    BUILD_COHORT_16S(
        collected_best_16s,
        collected_status_tables,
        metadata,
    )

    emit:
    all_best_16S = BUILD_COHORT_16S.out.best_fasta
    all_best_16S_manifest = BUILD_COHORT_16S.out.best_manifest
    all_partial_16S = BUILD_COHORT_16S.out.partial_fasta
    all_partial_16S_manifest = BUILD_COHORT_16S.out.partial_manifest
    versions = BUILD_COHORT_16S.out.versions
}
