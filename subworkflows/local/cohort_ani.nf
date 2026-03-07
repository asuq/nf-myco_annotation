include { BUILD_FASTANI_INPUTS } from '../../modules/local/build_fastani_inputs'
include { FASTANI } from '../../modules/local/fastani'
include { CLUSTER_ANI } from '../../modules/local/cluster_ani'

/*
 * Placeholder for cohort-level ANI preparation, FastANI, and clustering.
 */
workflow COHORT_ANI {
    take:
    ani_inputs

    main:
    results = ani_inputs

    emit:
    results = results
}
