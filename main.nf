#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO_DATASET_PREP } from './subworkflows/local/busco_dataset_prep'
include { INPUT_VALIDATION_AND_STAGING } from './subworkflows/local/input_validation_and_staging'
include { PER_SAMPLE_QC } from './subworkflows/local/per_sample_qc'

workflow {
    if (!params.sample_csv) {
        error "params.sample_csv is required."
    }
    if (!params.metadata) {
        error "params.metadata is required."
    }
    if (!(params.busco_lineages instanceof List) || params.busco_lineages.isEmpty()) {
        error "params.busco_lineages must be a non-empty list."
    }

    log.warn 'This phase wires only input validation, staging, Barrnap, paired CheckM2, gcode/QC summarisation, BUSCO dataset prep, and BUSCO offline runs.'
    log.warn 'Prokka, CRISPRCasFinder, PADLOC, eggNOG, ANI clustering, and final table assembly remain intentionally unwired here.'

    sampleCsv = Channel.fromPath(params.sample_csv, checkIfExists: true)
    metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    buscoLineages = Channel.fromList(params.busco_lineages as List<String>)

    INPUT_VALIDATION_AND_STAGING(sampleCsv, metadata)
    BUSCO_DATASET_PREP(buscoLineages)
    PER_SAMPLE_QC(
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        BUSCO_DATASET_PREP.out.datasets,
    )

    versions = INPUT_VALIDATION_AND_STAGING.out.versions
        .mix(BUSCO_DATASET_PREP.out.versions)
        .mix(PER_SAMPLE_QC.out.versions)

    emit:
    validated_samples = INPUT_VALIDATION_AND_STAGING.out.validated_samples
    accession_map = INPUT_VALIDATION_AND_STAGING.out.accession_map
    validation_warnings = INPUT_VALIDATION_AND_STAGING.out.validation_warnings
    initial_sample_status = INPUT_VALIDATION_AND_STAGING.out.sample_status
    staged_genomes = INPUT_VALIDATION_AND_STAGING.out.staged_genomes
    barrnap = PER_SAMPLE_QC.out.barrnap
    checkm2_gcode4 = PER_SAMPLE_QC.out.checkm2_gcode4
    checkm2_gcode11 = PER_SAMPLE_QC.out.checkm2_gcode11
    gcode_qc = PER_SAMPLE_QC.out.gcode_qc
    busco_datasets = BUSCO_DATASET_PREP.out.datasets
    busco = PER_SAMPLE_QC.out.busco
    busco_summaries = PER_SAMPLE_QC.out.busco_summaries
    versions = versions
}
