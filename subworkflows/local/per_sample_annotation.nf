include { PROKKA } from '../../modules/local/prokka'
include { CCFINDER } from '../../modules/local/ccfinder'
include { SUMMARISE_CCFINDER } from '../../modules/local/summarise_ccfinder'
include { PREPARE_PADLOC_INPUTS } from '../../modules/local/prepare_padloc_inputs'
include { PADLOC } from '../../modules/local/padloc'
include { EGGNOG } from '../../modules/local/eggnog'

/*
 * Placeholder for the gcode-gated annotation branch.
 */
workflow PER_SAMPLE_ANNOTATION {
    take:
    annotation_inputs

    main:
    results = annotation_inputs

    emit:
    results = results
}
