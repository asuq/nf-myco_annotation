/*
 * Summarise one sample's Barrnap results and emit separate intact and partial
 * cohort-candidate files without parsing TSV state in Groovy.
 */
process SUMMARISE_16S {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/16s" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['16S_status.tsv', 'best_16S.fna']
                ? filename
                : null
        },
    )

    input:
    tuple val(meta), path(rrna_gff), path(rrna_fasta), path(barrnap_log)

    output:
    tuple val(meta), path('best_16S.fna'), path('16S_status.tsv'), emit: sample_summaries
    tuple val(meta), path('cohort_best_16S.fna'), emit: intact_cohort_candidates
    tuple val(meta), path('cohort_intact_manifest_row.tsv'), emit: intact_manifest_rows
    tuple val(meta), path('cohort_partial_16S.fna'), emit: partial_cohort_candidates
    tuple val(meta), path('cohort_partial_manifest_row.tsv'), emit: partial_manifest_rows
    path 'versions.yml', emit: versions

    script:
    def isAtypical = (meta.is_atypical ?: false).toString()
    """
    summarise_16s.py \
        --accession "${meta.accession}" \
        --rrna-gff "${rrna_gff}" \
        --rrna-fasta "${rrna_fasta}" \
        --outdir . \
        --is-atypical "${isAtypical}"

    cohort_status="\$(awk -F '\t' 'NR == 2 { print \$2 }' 16S_status.tsv)"
    include_in_all="\$(awk -F '\t' 'NR == 2 { print \$5 }' 16S_status.tsv)"

    if [[ -s best_16S.fna ]] && [[ "\${include_in_all}" == "true" ]]; then
        cp best_16S.fna cohort_best_16S.fna
        cp 16S_status.tsv cohort_intact_manifest_row.tsv
    else
        : > cohort_best_16S.fna
        head -n 1 16S_status.tsv > cohort_intact_manifest_row.tsv
    fi

    if [[ -s best_16S.fna ]] && [[ "\${cohort_status}" == "partial" ]]; then
        cp best_16S.fna cohort_partial_16S.fna
        cp 16S_status.tsv cohort_partial_manifest_row.tsv
    else
        : > cohort_partial_16S.fna
        head -n 1 16S_status.tsv > cohort_partial_manifest_row.tsv
    fi

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/summarise_16s.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > best_16S.fna
    >sample_a 16S ribosomal RNA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    cat <<'EOF' > 16S_status.tsv
    accession	16S	best_16S_header	best_16S_length	include_in_all_best_16S	warnings
    sample_a	Yes	sample_a 16S ribosomal RNA	80	true
    EOF
    cp best_16S.fna cohort_best_16S.fna
    cp 16S_status.tsv cohort_intact_manifest_row.tsv
    : > cohort_partial_16S.fna
    head -n 1 16S_status.tsv > cohort_partial_manifest_row.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/summarise_16s.py"
    EOF
    '''
}
