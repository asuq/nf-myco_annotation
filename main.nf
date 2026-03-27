#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUILD_OUTPUT_CONTRACTS } from './modules/local/build_output_contracts'
include { BUSCO_DATASET_PREP } from './subworkflows/local/busco_dataset_prep'
include { COHORT_16S } from './subworkflows/local/cohort_16s'
include { COHORT_ANI } from './subworkflows/local/cohort_ani'
include { COHORT_TAXONOMY } from './subworkflows/local/cohort_taxonomy'
include { FINAL_OUTPUTS } from './subworkflows/local/final_outputs'
include { INPUT_VALIDATION_AND_STAGING } from './subworkflows/local/input_validation_and_staging'
include { PER_SAMPLE_ANNOTATION } from './subworkflows/local/per_sample_annotation'
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
    if (!params.checkm2_db) {
        error "params.checkm2_db is required."
    }
    if (!params.codetta_db) {
        error "params.codetta_db is required."
    }
    if (!params.eggnog_db) {
        error "params.eggnog_db is required."
    }
    if (!(params.busco_lineages instanceof List) || params.busco_lineages.isEmpty()) {
        error "params.busco_lineages must be a non-empty list."
    }

    log.warn 'PADLOC and eggNOG outputs are retained in sample folders but are intentionally excluded from master_table.tsv.'

    sampleCsv = Channel.fromPath(params.sample_csv, checkIfExists: true)
    metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    taxdump = Channel.fromPath(params.taxdump, checkIfExists: true)
    checkm2Db = Channel.fromPath(params.checkm2_db, checkIfExists: true)
    codettaDb = Channel.fromPath(params.codetta_db, checkIfExists: true)
    eggnogDb = Channel.fromPath(params.eggnog_db, checkIfExists: true)
    buscoLineagesList = params.busco_lineages as List<String>
    buscoLineages = Channel.fromList(buscoLineagesList)

    BUILD_OUTPUT_CONTRACTS(Channel.value(buscoLineagesList))

    INPUT_VALIDATION_AND_STAGING(
        sampleCsv,
        metadata,
        BUILD_OUTPUT_CONTRACTS.out.sample_status_columns,
    )
    BUSCO_DATASET_PREP(buscoLineages)
    COHORT_TAXONOMY(INPUT_VALIDATION_AND_STAGING.out.validated_samples, metadata, taxdump)
    PER_SAMPLE_QC(
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        checkm2Db,
        BUSCO_DATASET_PREP.out.datasets,
    )
    COHORT_16S(PER_SAMPLE_QC.out.barrnap)
    PER_SAMPLE_ANNOTATION(
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        PER_SAMPLE_QC.out.gcode_qc,
        codettaDb,
        eggnogDb,
    )
    COHORT_ANI(
        INPUT_VALIDATION_AND_STAGING.out.validated_samples,
        metadata,
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        PER_SAMPLE_QC.out.gcode_qc,
        COHORT_16S.out.sample_summaries,
        PER_SAMPLE_QC.out.busco_summaries,
    )
    FINAL_OUTPUTS(
        INPUT_VALIDATION_AND_STAGING.out.validated_samples,
        INPUT_VALIDATION_AND_STAGING.out.sample_status,
        metadata,
        COHORT_TAXONOMY.out.taxonomy,
        PER_SAMPLE_QC.out.gcode_qc,
        COHORT_16S.out.sample_summaries,
        COHORT_ANI.out.parsed_busco,
        PER_SAMPLE_ANNOTATION.out.codetta_summary,
        PER_SAMPLE_ANNOTATION.out.ccfinder_summary,
        PER_SAMPLE_ANNOTATION.out.prokka,
        PER_SAMPLE_ANNOTATION.out.padloc,
        PER_SAMPLE_ANNOTATION.out.eggnog,
        COHORT_ANI.out.clusters,
        COHORT_ANI.out.ani_metadata,
        COHORT_ANI.out.assembly_stats,
        COHORT_ANI.out.fastani_matrix,
        BUILD_OUTPUT_CONTRACTS.out.append_columns,
        BUILD_OUTPUT_CONTRACTS.out.sample_status_columns,
        INPUT_VALIDATION_AND_STAGING.out.versions
            .mix(BUILD_OUTPUT_CONTRACTS.out.versions)
            .mix(BUSCO_DATASET_PREP.out.versions)
            .mix(COHORT_TAXONOMY.out.versions)
            .mix(PER_SAMPLE_QC.out.versions)
            .mix(COHORT_16S.out.versions)
            .mix(PER_SAMPLE_ANNOTATION.out.versions)
            .mix(COHORT_ANI.out.versions),
        workflow.nextflow.version.toString(),
        workflow.manifest.version ?: 'NA',
        workflow.commitId ?: 'NA',
        workflow.containerEngine ?: 'none',
    )
}
