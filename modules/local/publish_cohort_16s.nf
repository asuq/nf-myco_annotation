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
    path cohort_best_fasta
    path cohort_best_manifest
    path cohort_partial_fasta
    path cohort_partial_manifest

    output:
    path 'all_best_16S.fna', emit: best_fasta
    path 'all_best_16S_manifest.tsv', emit: best_manifest
    path 'all_partial_16S.fna', emit: partial_fasta
    path 'all_partial_16S_manifest.tsv', emit: partial_manifest

    script:
    """
    cp "${cohort_best_fasta}" all_best_16S.fna
    cp "${cohort_best_manifest}" all_best_16S_manifest.tsv
    cp "${cohort_partial_fasta}" all_partial_16S.fna
    cp "${cohort_partial_manifest}" all_partial_16S_manifest.tsv
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
    : > all_partial_16S.fna
    cat <<'EOF' > all_partial_16S_manifest.tsv
    accession	16S	best_16S_header	best_16S_length	include_in_all_best_16S	warnings
    EOF
    '''
}
