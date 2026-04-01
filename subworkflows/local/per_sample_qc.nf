include { BARRNAP } from '../../modules/local/barrnap'
include { SUMMARISE_16S } from '../../modules/local/summarise_16s'
include { CHECKM2 as CHECKM2_GCODE4 } from '../../modules/local/checkm2'
include { CHECKM2 as CHECKM2_GCODE11 } from '../../modules/local/checkm2'
include { ASSIGN_GCODE_AND_QC } from '../../modules/local/assign_gcode_and_qc'
include { BUSCO } from '../../modules/local/busco'

/*
 * Run the per-sample QC branch limited to Barrnap, per-sample 16S summary
 * generation, paired CheckM2 runs, gcode assignment, and raw BUSCO offline
 * executions across the configured lineages.
 */
workflow PER_SAMPLE_QC {
    take:
    sample_genomes
    checkm2_db
    busco_datasets

    main:
    BARRNAP(sample_genomes)
    SUMMARISE_16S(BARRNAP.out.results)

    CHECKM2_GCODE4(sample_genomes.combine(checkm2_db), Channel.value(4))
    CHECKM2_GCODE11(sample_genomes.combine(checkm2_db), Channel.value(11))

    checkm2_gcode4_reports = CHECKM2_GCODE4.out.quality_report
        .map { meta, translationTable, report ->
            tuple(meta.accession, meta, report)
        }

    checkm2_gcode11_reports = CHECKM2_GCODE11.out.quality_report
        .map { meta, translationTable, report ->
            tuple(meta.accession, report)
        }

    paired_checkm2_reports = checkm2_gcode4_reports
        .join(checkm2_gcode11_reports)
        .map { accession, meta, gcode4Report, gcode11Report ->
            tuple(meta, gcode4Report, gcode11Report)
        }

    ASSIGN_GCODE_AND_QC(paired_checkm2_reports)

    busco_jobs = sample_genomes
        .combine(busco_datasets)
        .map { meta, genome, lineage, datasetDir ->
            tuple(meta, genome, lineage, datasetDir)
        }

    BUSCO(busco_jobs)

    versions = BARRNAP.out.versions
        .mix(SUMMARISE_16S.out.versions)
        .mix(CHECKM2_GCODE4.out.versions)
        .mix(CHECKM2_GCODE11.out.versions)
        .mix(ASSIGN_GCODE_AND_QC.out.versions)
        .mix(BUSCO.out.versions)

    emit:
    barrnap = BARRNAP.out.results
    sixteen_s_summaries = SUMMARISE_16S.out.summaries
    checkm2_gcode4 = CHECKM2_GCODE4.out.results
    checkm2_gcode11 = CHECKM2_GCODE11.out.results
    gcode_qc = ASSIGN_GCODE_AND_QC.out.summary
    busco = BUSCO.out.results
    busco_summaries = BUSCO.out.summary
    versions = versions
}
