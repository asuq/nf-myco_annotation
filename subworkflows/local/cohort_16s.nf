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
    unpackTuple = { item, channelName, expectedSize ->
        if (!(item instanceof List)) {
            def actualType = item == null ? 'null' : item.getClass().getName()
            throw new IllegalArgumentException(
                "${channelName} expected a ${expectedSize}-value tuple, received ${actualType}."
            )
        }
        if (item.size() != expectedSize) {
            def preview = item.take(Math.min(item.size(), 3))
            throw new IllegalArgumentException(
                "${channelName} expected a ${expectedSize}-value tuple, " +
                    "received ${item.size()} value(s): ${preview}"
            )
        }
        return item
    }

    collected_best_16s = sixteen_s_summaries
        .map { item ->
            def values = unpackTuple.call(item, 'sixteen_s_summaries', 3)
            values[1]
        }
        .collect()

    collected_status_tables = sixteen_s_summaries
        .map { item ->
            def values = unpackTuple.call(item, 'sixteen_s_summaries', 3)
            values[2]
        }
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
