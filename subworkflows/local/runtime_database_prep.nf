include { PREP_TAXDUMP_DATABASE } from '../../modules/local/prepare_runtime_database'
include { PREP_CODETTA_DATABASE } from '../../modules/local/prepare_runtime_database'
include { DOWNLOAD_CHECKM2_DATABASE } from '../../modules/local/download_checkm2_database'
include { DOWNLOAD_BUSCO_DATABASES } from '../../modules/local/download_busco_databases'
include { DOWNLOAD_EGGNOG_DATABASE } from '../../modules/local/download_eggnog_database'
include { FINALISE_RUNTIME_DATABASE as FINALISE_CHECKM2_DATABASE } from '../../modules/local/finalise_runtime_database'
include { FINALISE_RUNTIME_DATABASE as FINALISE_BUSCO_DATABASE } from '../../modules/local/finalise_runtime_database'
include { FINALISE_RUNTIME_DATABASE as FINALISE_EGGNOG_DATABASE } from '../../modules/local/finalise_runtime_database'
include { MERGE_RUNTIME_DATABASE_REPORTS } from '../../modules/local/merge_runtime_database_reports'

/*
 * Prepare runtime databases in their final tool-consumable directories.
 */
workflow RUNTIME_DATABASE_PREP {
    take:
    taxdump_request
    checkm2_request
    busco_request
    codetta_request
    eggnog_request

    main:
    PREP_TAXDUMP_DATABASE(taxdump_request)

    DOWNLOAD_CHECKM2_DATABASE(checkm2_request)
    FINALISE_CHECKM2_DATABASE(DOWNLOAD_CHECKM2_DATABASE.out.finalise_input)

    DOWNLOAD_BUSCO_DATABASES(busco_request)
    FINALISE_BUSCO_DATABASE(DOWNLOAD_BUSCO_DATABASES.out.finalise_input)

    PREP_CODETTA_DATABASE(codetta_request)

    DOWNLOAD_EGGNOG_DATABASE(eggnog_request)
    FINALISE_EGGNOG_DATABASE(DOWNLOAD_EGGNOG_DATABASE.out.finalise_input)

    reports = PREP_TAXDUMP_DATABASE.out.report
        .mix(FINALISE_CHECKM2_DATABASE.out.report)
        .mix(FINALISE_BUSCO_DATABASE.out.report)
        .mix(PREP_CODETTA_DATABASE.out.report)
        .mix(FINALISE_EGGNOG_DATABASE.out.report)
        .collect()

    MERGE_RUNTIME_DATABASE_REPORTS(reports)

    emit:
    report = MERGE_RUNTIME_DATABASE_REPORTS.out.report
    args = MERGE_RUNTIME_DATABASE_REPORTS.out.args
    versions = PREP_TAXDUMP_DATABASE.out.versions
        .mix(DOWNLOAD_CHECKM2_DATABASE.out.versions)
        .mix(DOWNLOAD_BUSCO_DATABASES.out.versions)
        .mix(PREP_CODETTA_DATABASE.out.versions)
        .mix(DOWNLOAD_EGGNOG_DATABASE.out.versions)
        .mix(FINALISE_CHECKM2_DATABASE.out.versions)
        .mix(FINALISE_BUSCO_DATABASE.out.versions)
        .mix(FINALISE_EGGNOG_DATABASE.out.versions)
        .mix(MERGE_RUNTIME_DATABASE_REPORTS.out.versions)
}
