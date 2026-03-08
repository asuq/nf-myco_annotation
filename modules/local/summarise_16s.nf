/*
 * Summarise one sample's Barrnap results and emit cohort-candidate files that
 * the subworkflow can concatenate without parsing TSV state in Groovy.
 */
process SUMMARISE_16S {
    tag "${meta.accession}"
    label 'process_single'

    input:
    tuple val(meta), path(rrna_gff), path(rrna_fasta), path(barrnap_log)

    output:
    tuple val(meta), path('best_16S.fna'), path('16S_status.tsv'), emit: sample_summaries
    tuple val(meta), path('cohort_best_16S.fna'), emit: cohort_candidates
    tuple val(meta), path('cohort_manifest_row.tsv'), emit: cohort_manifest_rows
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

    if [[ -s best_16S.fna ]] && awk -F '\t' 'NR == 2 && \$5 == "true" { found = 1 } END { exit(found ? 0 : 1) }' 16S_status.tsv; then
        cp best_16S.fna cohort_best_16S.fna
        cp 16S_status.tsv cohort_manifest_row.tsv
    else
        : > cohort_best_16S.fna
        head -n 1 16S_status.tsv > cohort_manifest_row.tsv
    fi

    cat <<EOF > versions.yml
    "${task.process}":
      python: "$(python3 --version 2>&1 | sed 's/^Python //')"
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
    cp 16S_status.tsv cohort_manifest_row.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/summarise_16s.py"
    EOF
    '''
}
