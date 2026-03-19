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
        rawValue ? new File(rawValue.toString()).canonicalFile.absolutePath : null
    }
    def buildMountedDestination = { destination ->
        if (!destination) {
            return null
        }
        def destinationFile = new File(destination)
        def parentFile = destinationFile.getParentFile()
        if (parentFile == null) {
            error "Runtime database destination must have a parent directory: ${destination}"
        }
        parentFile.mkdirs()
        return [
            destination,
            file(parentFile.absolutePath),
            destinationFile.getName(),
        ]
    }
    def destinations = [
        taxdump   : normaliseDestination.call(params.taxdump),
        checkm2   : normaliseDestination.call(params.checkm2_db),
        busco_root: normaliseDestination.call(params.busco_db),
        codetta   : normaliseDestination.call(params.codetta_db),
        eggnog    : normaliseDestination.call(params.eggnog_db),
    ]
    def mountedDestinations = destinations.collectEntries { component, destination ->
        [(component): buildMountedDestination.call(destination)]
    }
    if (!destinations.values().any { it != null }) {
        error "At least one database destination must be set for prepare_databases.nf."
    }

    def downloadEnabled = parseBoolean.call(params.download_missing_databases)
    def forceRebuild = parseBoolean.call(params.force_runtime_database_rebuild)
    def scratchRoot = params.runtime_db_scratch_root
        ? new File(params.runtime_db_scratch_root.toString()).canonicalFile.absolutePath
        : null
    def buildMountedScratchRoot = { destinationParent ->
        if (scratchRoot != null) {
            def scratchFile = new File(scratchRoot)
            scratchFile.mkdirs()
            return file(scratchFile.absolutePath)
        }
        return destinationParent
    }

    taxdumpRequest = destinations.taxdump
        ? Channel.of(
            tuple(
                mountedDestinations.taxdump[0],
                mountedDestinations.taxdump[1],
                mountedDestinations.taxdump[2],
                downloadEnabled,
                params.taxdump_version,
                buildMountedScratchRoot.call(mountedDestinations.taxdump[1]),
                forceRebuild,
            )
        )
        : Channel.empty()
    checkm2Request = destinations.checkm2
        ? Channel.of(
            tuple(
                mountedDestinations.checkm2[0],
                mountedDestinations.checkm2[1],
                mountedDestinations.checkm2[2],
                downloadEnabled,
                forceRebuild,
            )
        )
        : Channel.empty()
    buscoRequest = destinations.busco_root
        ? Channel.of(
            tuple(
                mountedDestinations.busco_root[0],
                mountedDestinations.busco_root[1],
                mountedDestinations.busco_root[2],
                downloadEnabled,
                params.busco_lineages as List<String>,
                forceRebuild,
            )
        )
        : Channel.empty()
    codettaRequest = destinations.codetta
        ? Channel.of(
            tuple(
                mountedDestinations.codetta[0],
                mountedDestinations.codetta[1],
                mountedDestinations.codetta[2],
                downloadEnabled,
                buildMountedScratchRoot.call(mountedDestinations.codetta[1]),
                forceRebuild,
            )
        )
        : Channel.empty()
    eggnogRequest = destinations.eggnog
        ? Channel.of(
            tuple(
                mountedDestinations.eggnog[0],
                mountedDestinations.eggnog[1],
                mountedDestinations.eggnog[2],
                downloadEnabled,
                forceRebuild,
            )
        )
        : Channel.empty()
    RUNTIME_DATABASE_PREP(
        taxdumpRequest,
        checkm2Request,
        buscoRequest,
        codettaRequest,
        eggnogRequest,
    )
}
