include { PREP_RUNTIME_DATABASE as PREP_TAXDUMP } from '../../modules/local/prepare_runtime_database'
include { PREP_RUNTIME_DATABASE as PREP_CHECKM2 } from '../../modules/local/prepare_runtime_database'
include { PREP_RUNTIME_DATABASE as PREP_EGGNOG } from '../../modules/local/prepare_runtime_database'
include { PREP_RUNTIME_DATABASE as PREP_PADLOC } from '../../modules/local/prepare_runtime_database'
include { PREP_BUSCO_DATABASES } from '../../modules/local/prepare_busco_databases'
include { MERGE_RUNTIME_DATABASE_REPORTS } from '../../modules/local/merge_runtime_database_reports'

/*
 * Prepare runtime databases under one canonical db_root.
 */
workflow RUNTIME_DATABASE_PREP {
    take:
    taxdump_request
    checkm2_request
    busco_request
    eggnog_request
    padloc_request
    db_root

    main:
    PREP_TAXDUMP(taxdump_request)
    PREP_CHECKM2(checkm2_request)
    PREP_BUSCO_DATABASES(busco_request)
    PREP_EGGNOG(eggnog_request)
    PREP_PADLOC(padloc_request)

    reports = PREP_TAXDUMP.out.report
        .mix(PREP_CHECKM2.out.report)
        .mix(PREP_BUSCO_DATABASES.out.report)
        .mix(PREP_EGGNOG.out.report)
        .mix(PREP_PADLOC.out.report)
        .collect()

    MERGE_RUNTIME_DATABASE_REPORTS(reports, db_root)

    emit:
    report = MERGE_RUNTIME_DATABASE_REPORTS.out.report
    args = MERGE_RUNTIME_DATABASE_REPORTS.out.args
    versions = PREP_TAXDUMP.out.versions
        .mix(PREP_CHECKM2.out.versions)
        .mix(PREP_BUSCO_DATABASES.out.versions)
        .mix(PREP_EGGNOG.out.versions)
        .mix(PREP_PADLOC.out.versions)
        .mix(MERGE_RUNTIME_DATABASE_REPORTS.out.versions)
}
