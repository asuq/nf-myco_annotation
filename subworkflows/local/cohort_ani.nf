include { SUMMARISE_BUSCO } from '../../modules/local/summarise_busco'
include { BUILD_FASTANI_INPUTS } from '../../modules/local/build_fastani_inputs'
include { FASTANI } from '../../modules/local/fastani'
include { CLUSTER_ANI } from '../../modules/local/cluster_ani'

/*
 * Build the ANI-eligible cohort, run all-vs-all FastANI, and cluster the
 * eligible subset using the existing cluster_ani.py implementation.
 */
workflow COHORT_ANI {
    take:
    validated_samples
    metadata
    staged_genomes
    gcode_qc
    sixteen_s_summaries
    busco_summaries

    main:
    SUMMARISE_BUSCO(busco_summaries)

    staged_manifest = Channel
        .of('accession\tinternal_id\tstaged_filename')
        .concat(
            staged_genomes.map { meta, genome ->
                "${meta.accession}\t${meta.internal_id}\t${genome.getName()}"
            }
        )
        .collectFile(name: 'staged_genomes.tsv', newLine: true)

    staged_fasta_files = staged_genomes
        .map { meta, genome -> genome }
        .collect()

    combined_checkm2 = gcode_qc
        .map { meta, summary -> summary }
        .collectFile(
            name: 'checkm2_summaries.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_16s = sixteen_s_summaries
        .map { meta, best16s, status -> status }
        .collectFile(
            name: '16s_statuses.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_busco = SUMMARISE_BUSCO.out.summary
        .map { meta, lineage, summary -> summary }
        .collectFile(
            name: 'busco_summaries.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    BUILD_FASTANI_INPUTS(
        validated_samples,
        metadata,
        staged_manifest,
        staged_fasta_files,
        combined_checkm2,
        combined_16s,
        combined_busco,
    )
    FASTANI(BUILD_FASTANI_INPUTS.out.fastani_inputs, BUILD_FASTANI_INPUTS.out.paths)
    CLUSTER_ANI(FASTANI.out.matrix, BUILD_FASTANI_INPUTS.out.metadata)

    versions = SUMMARISE_BUSCO.out.versions
        .mix(BUILD_FASTANI_INPUTS.out.versions)
        .mix(FASTANI.out.versions)
        .mix(CLUSTER_ANI.out.versions)

    emit:
    parsed_busco = SUMMARISE_BUSCO.out.summary
    fastani_inputs = BUILD_FASTANI_INPUTS.out.fastani_inputs
    fastani_paths = BUILD_FASTANI_INPUTS.out.paths
    ani_metadata = BUILD_FASTANI_INPUTS.out.metadata
    ani_exclusions = BUILD_FASTANI_INPUTS.out.exclusions
    fastani_matrix = FASTANI.out.matrix
    fastani_raw = FASTANI.out.raw_output
    fastani_log = FASTANI.out.log
    clusters = CLUSTER_ANI.out.clusters
    representatives = CLUSTER_ANI.out.representatives
    versions = versions
}
