#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUNTIME_DATABASE_PREP } from './subworkflows/local/runtime_database_prep'

workflow {
    if (!params.db_root) {
        error "params.db_root is required."
    }
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

    def downloadEnabled = parseBoolean.call(params.download_missing_databases)
    def forceRebuild = parseBoolean.call(params.force_runtime_database_rebuild)
    def linkMode = params.runtime_db_link_mode ?: 'copy'
    def scratchRoot = params.runtime_db_scratch_root
        ? new File(params.runtime_db_scratch_root.toString()).absolutePath
        : null

    def dbRootFile = new File(params.db_root.toString()).absoluteFile
    dbRootFile.mkdirs()
    def dbRoot = file(dbRootFile.toString(), checkIfExists: true)
    def placeholder = file("${projectDir}/assets/runtime_db_placeholder.txt", checkIfExists: true)

    def resolveLocalSource = { rawValue, componentName ->
        if (!rawValue) {
            return [placeholder, false]
        }
        def localPath = new File(rawValue.toString()).absoluteFile
        if (!localPath.exists()) {
            log.warn "Local ${componentName} source does not exist and will be ignored: ${localPath}"
            return [placeholder, false]
        }
        return [file(localPath.toString(), checkIfExists: true), true]
    }

    def taxdumpSource = resolveLocalSource.call(params.taxdump_source, 'taxdump')
    def checkm2Source = resolveLocalSource.call(params.checkm2_source, 'CheckM2')
    def buscoSource = resolveLocalSource.call(params.busco_source_root, 'BUSCO')
    def eggnogSource = resolveLocalSource.call(params.eggnog_source, 'eggNOG')
    def padlocSource = resolveLocalSource.call(params.padloc_source, 'PADLOC')

    taxdumpRequest = Channel.of(
        tuple(
            'taxdump',
            dbRoot,
            taxdumpSource[0],
            taxdumpSource[1],
            downloadEnabled,
            params.taxdump_version,
            linkMode,
            scratchRoot,
            forceRebuild,
        )
    )
    checkm2Request = Channel.of(
        tuple(
            'checkm2',
            dbRoot,
            checkm2Source[0],
            checkm2Source[1],
            downloadEnabled,
            params.checkm2_version,
            linkMode,
            scratchRoot,
            forceRebuild,
        )
    )
    buscoRequest = Channel.of(
        tuple(
            dbRoot,
            buscoSource[0],
            buscoSource[1],
            downloadEnabled,
            params.busco_version,
            params.busco_lineages as List<String>,
            linkMode,
            scratchRoot,
            forceRebuild,
        )
    )
    eggnogRequest = Channel.of(
        tuple(
            'eggnog',
            dbRoot,
            eggnogSource[0],
            eggnogSource[1],
            downloadEnabled,
            params.eggnog_version,
            linkMode,
            scratchRoot,
            forceRebuild,
        )
    )
    padlocRequest = Channel.of(
        tuple(
            'padloc',
            dbRoot,
            padlocSource[0],
            padlocSource[1],
            downloadEnabled,
            params.padloc_version,
            linkMode,
            scratchRoot,
            forceRebuild,
        )
    )

    RUNTIME_DATABASE_PREP(
        taxdumpRequest,
        checkm2Request,
        buscoRequest,
        eggnogRequest,
        padlocRequest,
        dbRootFile.toString(),
    )
}
