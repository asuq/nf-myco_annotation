/*
 * Publish the collected cohort 16S artefacts under a stable results path.
 */
process PUBLISH_COHORT_16S {
    tag "cohort_16s"
    label 'process_single'
    cache 'deep'
    publishDir(
        "${params.outdir}/cohort/16s",
        mode: 'copy',
        overwrite: true,
    )

    input:
    path cohort_fasta
    path cohort_manifest

    output:
    path 'all_best_16S.fna', emit: fasta
    path 'all_best_16S_manifest.tsv', emit: manifest

    script:
    """
    cp "${cohort_fasta}" all_best_16S.fna
    cp "${cohort_manifest}" all_best_16S_manifest.tsv
    """

    stub:
    '''
    cat <<'EOF' > all_best_16S.fna
    >sample_a 16S ribosomal RNA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    cat <<'EOF' > all_best_16S_manifest.tsv
    accession	16S	best_16S_header	best_16S_length	include_in_all_best_16S	warnings
    sample_a	Yes	sample_a 16S ribosomal RNA	80	true
    EOF
    '''
}
