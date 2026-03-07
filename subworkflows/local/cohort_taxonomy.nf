include { TAXONOMY_EXPAND } from '../../modules/local/taxonomy_expand'

/*
 * Expand taxonomy for the requested sample set from a pinned user-supplied
 * taxdump directory.
 */
workflow COHORT_TAXONOMY {
    take:
    validated_samples
    metadata
    taxdump

    main:
    TAXONOMY_EXPAND(validated_samples, metadata, taxdump)

    emit:
    taxonomy = TAXONOMY_EXPAND.out.taxonomy
    versions = TAXONOMY_EXPAND.out.versions
}
