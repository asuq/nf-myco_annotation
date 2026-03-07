include { SUMMARISE_16S } from '../../modules/local/summarise_16s'

/*
 * Convert per-sample Barrnap outputs into per-sample 16S summaries and
 * concatenate the subset eligible for the cohort FASTA.
 */
workflow COHORT_16S {
    take:
    barrnap_outputs

    main:
    SUMMARISE_16S(barrnap_outputs)

    all_best_16S = SUMMARISE_16S.out.cohort_candidates
        .map { meta, cohortBest16S -> cohortBest16S }
        .collectFile(name: 'all_best_16S.fna')

    all_best_16S_manifest = SUMMARISE_16S.out.cohort_manifest_rows
        .map { meta, manifestRow -> manifestRow }
        .collectFile(
            name: 'all_best_16S_manifest.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    emit:
    sample_summaries = SUMMARISE_16S.out.sample_summaries
    all_best_16S = all_best_16S
    all_best_16S_manifest = all_best_16S_manifest
    versions = SUMMARISE_16S.out.versions
}
