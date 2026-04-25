#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO_DATASET_PREP } from './subworkflows/local/busco_dataset_prep'
include { COHORT_16S } from './subworkflows/local/cohort_16s'
include { COHORT_ANI } from './subworkflows/local/cohort_ani'
include { COHORT_TAXONOMY } from './subworkflows/local/cohort_taxonomy'
include { FINAL_OUTPUTS } from './subworkflows/local/final_outputs'
include { INPUT_VALIDATION_AND_STAGING } from './subworkflows/local/input_validation_and_staging'
include { PER_SAMPLE_ANNOTATION } from './subworkflows/local/per_sample_annotation'
include { PER_SAMPLE_QC } from './subworkflows/local/per_sample_qc'

workflow {
    def normaliseBuscoLineages = { rawValue ->
        def rawItems = rawValue instanceof List ? rawValue : [rawValue]
        def lineages = rawItems
            .findAll { it != null }
            .collectMany { it.toString().split(',') as List }
            .collect { it.trim() }
            .findAll { it }
        if (lineages.isEmpty()) {
            error "params.busco_lineages must resolve to one or more lineage names."
        }
        def seen = [] as Set
        def duplicates = lineages.findAll { !seen.add(it) }.unique()
        if (!duplicates.isEmpty()) {
            error "params.busco_lineages must resolve to unique lineage names: ${duplicates.join(', ')}"
        }
        return lineages
    }

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
    buscoLineagesList = normaliseBuscoLineages.call(params.busco_lineages)
    primaryBuscoColumn = (params.busco_primary_column ?: "BUSCO_${buscoLineagesList[0]}").toString()

    log.warn 'PADLOC and eggNOG outputs are retained in sample folders but are intentionally excluded from master_table.tsv.'

    sampleCsv = Channel.fromPath(params.sample_csv, checkIfExists: true)
    metadata = Channel.value(file(params.metadata, checkIfExists: true))
    taxdump = Channel.fromPath(params.taxdump, checkIfExists: true)
    checkm2Db = Channel.fromPath(params.checkm2_db, checkIfExists: true)
    codettaDb = Channel.fromPath(params.codetta_db, checkIfExists: true)
    eggnogDb = Channel.fromPath(params.eggnog_db, checkIfExists: true)
    buscoLineages = Channel.fromList(buscoLineagesList)

    INPUT_VALIDATION_AND_STAGING(
        sampleCsv,
        metadata,
        Channel.value(buscoLineagesList),
    )
    BUSCO_DATASET_PREP(buscoLineages)
    COHORT_TAXONOMY(INPUT_VALIDATION_AND_STAGING.out.validated_samples, metadata, taxdump)
    PER_SAMPLE_QC(
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        checkm2Db,
        BUSCO_DATASET_PREP.out.datasets,
    )
    COHORT_16S(PER_SAMPLE_QC.out.sixteen_s_summaries, metadata)
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
        PER_SAMPLE_QC.out.sixteen_s_summaries,
        PER_SAMPLE_QC.out.busco_summaries,
        Channel.value(primaryBuscoColumn),
    )
    FINAL_OUTPUTS(
        INPUT_VALIDATION_AND_STAGING.out.validated_samples,
        INPUT_VALIDATION_AND_STAGING.out.sample_status,
        metadata,
        COHORT_TAXONOMY.out.taxonomy,
        PER_SAMPLE_QC.out.gcode_qc,
        PER_SAMPLE_QC.out.sixteen_s_summaries,
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
        Channel.value(buscoLineagesList),
        Channel.value(primaryBuscoColumn),
        INPUT_VALIDATION_AND_STAGING.out.versions
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
