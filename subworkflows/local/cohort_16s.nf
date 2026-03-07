include { SUMMARISE_16S } from '../../modules/local/summarise_16s'

/*
 * Placeholder for cohort-level 16S aggregation.
 */
workflow COHORT_16S {
    take:
    barrnap_outputs

    main:
    results = barrnap_outputs

    emit:
    results = results
}
