#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO_DATASET_PREP } from './subworkflows/local/busco_dataset_prep'
include { COHORT_16S } from './subworkflows/local/cohort_16s'
include { COHORT_ANI } from './subworkflows/local/cohort_ani'
include { COHORT_TAXONOMY } from './subworkflows/local/cohort_taxonomy'
include { INPUT_VALIDATION_AND_STAGING } from './subworkflows/local/input_validation_and_staging'
include { PER_SAMPLE_QC } from './subworkflows/local/per_sample_qc'

workflow {
    if (!params.sample_csv) {
        error "params.sample_csv is required."
    }
    if (!params.metadata) {
        error "params.metadata is required."
    }
    if (!params.taxdump) {
        error "params.taxdump is required."
    }
    if (!(params.busco_lineages instanceof List) || params.busco_lineages.isEmpty()) {
        error "params.busco_lineages must be a non-empty list."
    }

    log.warn 'This phase wires validation, staging, taxonomy expansion, 16S summarisation, BUSCO parsing, and the ANI cohort branch.'
    log.warn 'Prokka, CRISPRCasFinder, PADLOC, eggNOG, and final table/status aggregation remain intentionally unwired here.'

    sampleCsv = Channel.fromPath(params.sample_csv, checkIfExists: true)
    metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    taxdump = Channel.fromPath(params.taxdump, checkIfExists: true)
    buscoLineages = Channel.fromList(params.busco_lineages as List<String>)

    INPUT_VALIDATION_AND_STAGING(sampleCsv, metadata)
    BUSCO_DATASET_PREP(buscoLineages)
    COHORT_TAXONOMY(INPUT_VALIDATION_AND_STAGING.out.validated_samples, metadata, taxdump)
    PER_SAMPLE_QC(
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        BUSCO_DATASET_PREP.out.datasets,
    )
    COHORT_16S(PER_SAMPLE_QC.out.barrnap)
    COHORT_ANI(
        INPUT_VALIDATION_AND_STAGING.out.validated_samples,
        metadata,
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        PER_SAMPLE_QC.out.gcode_qc,
        COHORT_16S.out.sample_summaries,
        PER_SAMPLE_QC.out.busco_summaries,
    )

    versions = INPUT_VALIDATION_AND_STAGING.out.versions
        .mix(BUSCO_DATASET_PREP.out.versions)
        .mix(COHORT_TAXONOMY.out.versions)
        .mix(PER_SAMPLE_QC.out.versions)
        .mix(COHORT_16S.out.versions)
        .mix(COHORT_ANI.out.versions)

    emit:
    validated_samples = INPUT_VALIDATION_AND_STAGING.out.validated_samples
    accession_map = INPUT_VALIDATION_AND_STAGING.out.accession_map
    validation_warnings = INPUT_VALIDATION_AND_STAGING.out.validation_warnings
    initial_sample_status = INPUT_VALIDATION_AND_STAGING.out.sample_status
    staged_genomes = INPUT_VALIDATION_AND_STAGING.out.staged_genomes
    taxonomy = COHORT_TAXONOMY.out.taxonomy
    barrnap = PER_SAMPLE_QC.out.barrnap
    checkm2_gcode4 = PER_SAMPLE_QC.out.checkm2_gcode4
    checkm2_gcode11 = PER_SAMPLE_QC.out.checkm2_gcode11
    gcode_qc = PER_SAMPLE_QC.out.gcode_qc
    sixteen_s = COHORT_16S.out.sample_summaries
    all_best_16S = COHORT_16S.out.all_best_16S
    all_best_16S_manifest = COHORT_16S.out.all_best_16S_manifest
    busco_datasets = BUSCO_DATASET_PREP.out.datasets
    busco = PER_SAMPLE_QC.out.busco
    busco_summaries = PER_SAMPLE_QC.out.busco_summaries
    parsed_busco_summaries = COHORT_ANI.out.parsed_busco
    fastani_inputs = COHORT_ANI.out.fastani_inputs
    fastani_paths = COHORT_ANI.out.fastani_paths
    ani_metadata = COHORT_ANI.out.ani_metadata
    ani_exclusions = COHORT_ANI.out.ani_exclusions
    fastani_matrix = COHORT_ANI.out.fastani_matrix
    fastani_raw = COHORT_ANI.out.fastani_raw
    fastani_log = COHORT_ANI.out.fastani_log
    ani_clusters = COHORT_ANI.out.clusters
    ani_representatives = COHORT_ANI.out.representatives
    versions = versions
}
