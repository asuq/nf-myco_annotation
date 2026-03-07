include { TAXONOMY_EXPAND } from '../../modules/local/taxonomy_expand'

/*
 * Placeholder for cohort-level taxonomy expansion.
 */
workflow COHORT_TAXONOMY {
    take:
    taxonomy_inputs

    main:
    results = taxonomy_inputs

    emit:
    results = results
}
