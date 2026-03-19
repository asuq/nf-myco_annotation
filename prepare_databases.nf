#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Register the main workflow subworkflows so shared config selectors remain
 * valid without weakening selector validation globally.
 */
include { BUSCO_DATASET_PREP as UNUSED_BUSCO_DATASET_PREP } from './subworkflows/local/busco_dataset_prep'
include { COHORT_16S as UNUSED_COHORT_16S } from './subworkflows/local/cohort_16s'
include { COHORT_ANI as UNUSED_COHORT_ANI } from './subworkflows/local/cohort_ani'
include { COHORT_TAXONOMY as UNUSED_COHORT_TAXONOMY } from './subworkflows/local/cohort_taxonomy'
include { FINAL_OUTPUTS as UNUSED_FINAL_OUTPUTS } from './subworkflows/local/final_outputs'
include { INPUT_VALIDATION_AND_STAGING as UNUSED_INPUT_VALIDATION_AND_STAGING } from './subworkflows/local/input_validation_and_staging'
include { PER_SAMPLE_ANNOTATION as UNUSED_PER_SAMPLE_ANNOTATION } from './subworkflows/local/per_sample_annotation'
include { PER_SAMPLE_QC as UNUSED_PER_SAMPLE_QC } from './subworkflows/local/per_sample_qc'
include { RUNTIME_DATABASE_PREP } from './subworkflows/local/runtime_database_prep'

workflow {
    if (!(params.busco_lineages instanceof List) || params.busco_lineages.isEmpty()) {
        error "params.busco_lineages must be a non-empty list."
    }

    def parseBoolean = { value ->
        if (value instanceof Boolean) {
            return value
        }
        if (value == null) {
            return false
        }
        return value.toString().toBoolean()
    }

    def normaliseDestination = { rawValue ->
        rawValue ? new File(rawValue.toString()).absolutePath : null
    }

    def destinations = [
        taxdump   : normaliseDestination.call(params.taxdump),
        checkm2   : normaliseDestination.call(params.checkm2_db),
        busco_root: normaliseDestination.call(params.busco_db),
        codetta   : normaliseDestination.call(params.codetta_db),
        eggnog    : normaliseDestination.call(params.eggnog_db),
    ]
    if (!destinations.values().any { it != null }) {
        error "At least one database destination must be set for prepare_databases.nf."
    }

    def downloadEnabled = parseBoolean.call(params.download_missing_databases)
    def forceRebuild = parseBoolean.call(params.force_runtime_database_rebuild)
    def scratchRoot = params.runtime_db_scratch_root
        ? new File(params.runtime_db_scratch_root.toString()).absolutePath
        : null

    taxdumpRequest = destinations.taxdump
        ? Channel.of(
            tuple(
                destinations.taxdump,
                downloadEnabled,
                params.taxdump_version,
                scratchRoot,
                forceRebuild,
            )
        )
        : Channel.empty()
    checkm2Request = destinations.checkm2
        ? Channel.of(tuple(destinations.checkm2, downloadEnabled, forceRebuild))
        : Channel.empty()
    buscoRequest = destinations.busco_root
        ? Channel.of(
            tuple(
                destinations.busco_root,
                downloadEnabled,
                params.busco_lineages as List<String>,
                forceRebuild,
            )
        )
        : Channel.empty()
    codettaRequest = destinations.codetta
        ? Channel.of(
            tuple(
                destinations.codetta,
                downloadEnabled,
                scratchRoot,
                forceRebuild,
            )
        )
        : Channel.empty()
    eggnogRequest = destinations.eggnog
        ? Channel.of(tuple(destinations.eggnog, downloadEnabled, forceRebuild))
        : Channel.empty()
    RUNTIME_DATABASE_PREP(
        taxdumpRequest,
        checkm2Request,
        buscoRequest,
        codettaRequest,
        eggnogRequest,
    )
}
