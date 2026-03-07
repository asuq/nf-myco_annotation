include { BARRNAP } from '../../modules/local/barrnap'
include { CHECKM2 } from '../../modules/local/checkm2'
include { ASSIGN_GCODE_AND_QC } from '../../modules/local/assign_gcode_and_qc'
include { BUSCO } from '../../modules/local/busco'
include { SUMMARISE_BUSCO } from '../../modules/local/summarise_busco'

/*
 * Placeholder for per-sample QC and BUSCO processing.
 */
workflow PER_SAMPLE_QC {
    take:
    qc_inputs

    main:
    results = qc_inputs

    emit:
    results = results
}
