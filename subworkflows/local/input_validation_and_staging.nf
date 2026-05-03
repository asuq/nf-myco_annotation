include { VALIDATE_INPUTS } from '../../modules/local/validate_inputs'
include { STAGE_INPUTS } from '../../modules/local/stage_inputs'

/*
 * Validate top-level inputs once, then stage every requested genome to a
 * canonical internal-ID-safe FASTA for downstream per-sample tools.
 */
workflow INPUT_VALIDATION_AND_STAGING {
    take:
    sample_csv
    metadata
    busco_lineages

    main:
    VALIDATE_INPUTS(sample_csv, metadata, busco_lineages)

    sample_genomes = VALIDATE_INPUTS.out.validated_samples
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = row.collectEntries { key, value -> [(key): value] }
            def genomePath = java.nio.file.Path.of(row.genome_fasta).toRealPath()
            tuple(meta, genomePath)
        }

    STAGE_INPUTS(sample_genomes)

    versions = VALIDATE_INPUTS.out.versions.mix(STAGE_INPUTS.out.versions)

    emit:
    validated_samples = VALIDATE_INPUTS.out.validated_samples
    accession_map = VALIDATE_INPUTS.out.accession_map
    validation_warnings = VALIDATE_INPUTS.out.validation_warnings
    sample_status = VALIDATE_INPUTS.out.sample_status
    staged_genomes = STAGE_INPUTS.out.staged_fasta
    staged_fai = STAGE_INPUTS.out.fai
    versions = versions
}
