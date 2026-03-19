include { PROKKA } from '../../modules/local/prokka'
include { CCFINDER } from '../../modules/local/ccfinder'
include { CODETTA } from '../../modules/local/codetta'
include { SUMMARISE_CODETTA } from '../../modules/local/summarise_codetta'
include { SUMMARISE_CCFINDER } from '../../modules/local/summarise_ccfinder'
include { PADLOC } from '../../modules/local/padloc'
include { EGGNOG } from '../../modules/local/eggnog'

/*
 * Run the gcode-gated annotation branch by joining staged genomes to the
 * per-sample CheckM2 summary that assigned the translation table.
 */
workflow PER_SAMPLE_ANNOTATION {
    take:
    sample_genomes
    gcode_summaries
    codetta_db
    eggnog_db

    main:
    def eggnogOnlyAccessions = null
    if (params.eggnog_only_accessions != null) {
        def rawAccessions = params.eggnog_only_accessions instanceof Collection
            ? (params.eggnog_only_accessions as Collection)
            : params.eggnog_only_accessions.toString().split(',')
        eggnogOnlyAccessions = rawAccessions
            .collect { it.toString().trim() }
            .findAll { !it.isEmpty() } as Set
    }

    sample_by_accession = sample_genomes.map { meta, genome ->
        tuple(meta.accession, meta, genome)
    }

    gcode_by_accession = gcode_summaries
        .map { meta, summary -> summary }
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.accession.toString(), row.Gcode.toString())
        }

    annotation_candidates = sample_by_accession
        .join(gcode_by_accession)
        .filter { accession, meta, genome, gcode ->
            ['4', '11'].contains(gcode.toString())
        }
        .map { accession, meta, genome, gcode ->
            tuple(meta, genome, gcode.toString())
        }

    CODETTA(sample_genomes.combine(codetta_db))
    SUMMARISE_CODETTA(CODETTA.out.summary_input)

    PROKKA(annotation_candidates)
    CCFINDER(annotation_candidates)
    SUMMARISE_CCFINDER(CCFINDER.out.result_json)
    PADLOC(PROKKA.out.padloc_inputs)

    eggnog_inputs = PROKKA.out.eggnog_inputs
    if (eggnogOnlyAccessions != null && !eggnogOnlyAccessions.isEmpty()) {
        eggnog_inputs = eggnog_inputs.filter { meta, faa ->
            eggnogOnlyAccessions.contains(meta.accession.toString())
        }
    }

    EGGNOG(eggnog_inputs.combine(eggnog_db))

    versions = PROKKA.out.versions
        .mix(CODETTA.out.versions)
        .mix(SUMMARISE_CODETTA.out.versions)
        .mix(CCFINDER.out.versions)
        .mix(SUMMARISE_CCFINDER.out.versions)
        .mix(PADLOC.out.versions)
        .mix(EGGNOG.out.versions)

    emit:
    codetta = CODETTA.out.results
    codetta_summary = SUMMARISE_CODETTA.out.summary
    prokka = PROKKA.out.results
    ccfinder = CCFINDER.out.results
    ccfinder_summary = SUMMARISE_CCFINDER.out.summaries
    padloc = PADLOC.out.results
    eggnog = EGGNOG.out.results
    versions = versions
}
