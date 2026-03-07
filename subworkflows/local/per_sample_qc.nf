include { BARRNAP } from '../../modules/local/barrnap'
include { CHECKM2 as CHECKM2_GCODE4 } from '../../modules/local/checkm2'
include { CHECKM2 as CHECKM2_GCODE11 } from '../../modules/local/checkm2'
include { ASSIGN_GCODE_AND_QC } from '../../modules/local/assign_gcode_and_qc'
include { BUSCO } from '../../modules/local/busco'
include { SUMMARISE_BUSCO } from '../../modules/local/summarise_busco'

/*
 * Run the per-sample QC branch: Barrnap, paired CheckM2 runs, and BUSCO across
 * the configured lineage datasets.
 */
workflow PER_SAMPLE_QC {
    take:
    sample_genomes
    busco_datasets

    main:
    BARRNAP(sample_genomes)

    CHECKM2_GCODE4(sample_genomes, Channel.value(4))
    CHECKM2_GCODE11(sample_genomes, Channel.value(11))

    paired_checkm2_reports = CHECKM2_GCODE4.out.quality_report
        .mix(CHECKM2_GCODE11.out.quality_report)
        .map { meta, translationTable, report ->
            tuple(
                meta.accession,
                [
                    meta: meta,
                    translation_table: translationTable as Integer,
                    report: report,
                ],
            )
        }
        .groupTuple()
        .map { accession, groupedReports ->
            def reportByTable = groupedReports.collectEntries { entry ->
                [(entry.translation_table as Integer): entry.report]
            }
            tuple(groupedReports[0].meta, reportByTable[4], reportByTable[11])
        }

    ASSIGN_GCODE_AND_QC(paired_checkm2_reports)

    busco_jobs = sample_genomes
        .combine(busco_datasets)
        .map { meta, genome, lineage, datasetDir ->
            tuple(meta, genome, lineage, datasetDir)
        }

    BUSCO(busco_jobs)
    SUMMARISE_BUSCO(BUSCO.out.summary)

    versions = BARRNAP.out.versions
        .mix(CHECKM2_GCODE4.out.versions)
        .mix(CHECKM2_GCODE11.out.versions)
        .mix(ASSIGN_GCODE_AND_QC.out.versions)
        .mix(BUSCO.out.versions)
        .mix(SUMMARISE_BUSCO.out.versions)

    emit:
    barrnap = BARRNAP.out.results
    checkm2 = ASSIGN_GCODE_AND_QC.out.summary
    busco = SUMMARISE_BUSCO.out.summary
    versions = versions
}
