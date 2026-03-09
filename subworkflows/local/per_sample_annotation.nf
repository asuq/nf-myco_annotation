include { PROKKA } from '../../modules/local/prokka'
include { CCFINDER } from '../../modules/local/ccfinder'
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

    main:
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

    PROKKA(annotation_candidates)
    CCFINDER(annotation_candidates)
    ccfinder_summary_inputs = CCFINDER.out.results.map { meta, ccfinder_dir, result_json, log ->
        tuple(meta, ccfinder_dir, result_json)
    }

    SUMMARISE_CCFINDER(ccfinder_summary_inputs)
    PADLOC(PROKKA.out.padloc_inputs)
    EGGNOG(PROKKA.out.eggnog_inputs)

    versions = PROKKA.out.versions
        .mix(CCFINDER.out.versions)
        .mix(SUMMARISE_CCFINDER.out.versions)
        .mix(PADLOC.out.versions)
        .mix(EGGNOG.out.versions)

    emit:
    prokka = PROKKA.out.results
    ccfinder = CCFINDER.out.results
    ccfinder_summary = SUMMARISE_CCFINDER.out.summaries
    padloc = PADLOC.out.results
    eggnog = EGGNOG.out.results
    versions = versions
}
