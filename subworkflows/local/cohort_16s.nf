include { PUBLISH_COHORT_16S } from '../../modules/local/publish_cohort_16s'
include { SUMMARISE_16S } from '../../modules/local/summarise_16s'

/*
 * Convert per-sample Barrnap outputs into per-sample 16S summaries and
 * concatenate separate intact and partial cohort 16S FASTA sets.
 */
workflow COHORT_16S {
    take:
    barrnap_outputs

    main:
    SUMMARISE_16S(barrnap_outputs)

    cohortManifestHeader = "${projectDir}/assets/testdata/headers/16s_status.tsv"

    collected_all_best_16S = SUMMARISE_16S.out.intact_cohort_candidates
        .map { meta, cohortBest16S -> cohortBest16S }
        .collectFile(name: 'all_best_16S.fna')

    intact_manifest_header = Channel.fromPath(cohortManifestHeader, checkIfExists: true)
    intact_manifest_rows = SUMMARISE_16S.out.intact_manifest_rows
        .map { meta, manifestRow -> manifestRow }
        .filter { manifestRow -> manifestRow.toFile().length() > 0 }

    collected_all_best_16S_manifest = intact_manifest_header
        .concat(intact_manifest_rows)
        .collectFile(
            name: 'all_best_16S_manifest.tsv',
            keepHeader: true,
            skip: 1,
        )

    collected_all_partial_16S = SUMMARISE_16S.out.partial_cohort_candidates
        .map { meta, cohortPartial16S -> cohortPartial16S }
        .collectFile(name: 'all_partial_16S.fna')

    partial_manifest_header = Channel.fromPath(cohortManifestHeader, checkIfExists: true)
    partial_manifest_rows = SUMMARISE_16S.out.partial_manifest_rows
        .map { meta, manifestRow -> manifestRow }
        .filter { manifestRow -> manifestRow.toFile().length() > 0 }

    collected_all_partial_16S_manifest = partial_manifest_header
        .concat(partial_manifest_rows)
        .collectFile(
            name: 'all_partial_16S_manifest.tsv',
            keepHeader: true,
            skip: 1,
        )

    PUBLISH_COHORT_16S(
        collected_all_best_16S,
        collected_all_best_16S_manifest,
        collected_all_partial_16S,
        collected_all_partial_16S_manifest,
    )

    emit:
    sample_summaries = SUMMARISE_16S.out.sample_summaries
    all_best_16S = PUBLISH_COHORT_16S.out.best_fasta
    all_best_16S_manifest = PUBLISH_COHORT_16S.out.best_manifest
    all_partial_16S = PUBLISH_COHORT_16S.out.partial_fasta
    all_partial_16S_manifest = PUBLISH_COHORT_16S.out.partial_manifest
    versions = SUMMARISE_16S.out.versions
}
