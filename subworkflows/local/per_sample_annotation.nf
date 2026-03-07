include { PROKKA } from '../../modules/local/prokka'
include { CCFINDER } from '../../modules/local/ccfinder'
include { SUMMARISE_CCFINDER } from '../../modules/local/summarise_ccfinder'
include { PADLOC } from '../../modules/local/padloc'
include { EGGNOG } from '../../modules/local/eggnog'

/*
 * Run the gcode-gated annotation branch. The caller may pass mixed samples; the
 * workflow keeps only gcode 4/11 entries for downstream annotation modules.
 */
workflow PER_SAMPLE_ANNOTATION {
    take:
    annotation_inputs

    main:
    annotation_candidates = annotation_inputs
        .filter { meta, genome, gcode ->
            ['4', '11'].contains(gcode.toString())
        }
        .map { meta, genome, gcode ->
            tuple(meta, genome, gcode.toString())
        }

    PROKKA(annotation_candidates)
    CCFINDER(annotation_candidates)
    SUMMARISE_CCFINDER(CCFINDER.out.result_json)
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
