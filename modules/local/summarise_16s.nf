/*
 * Placeholder for per-sample 16S status calls and cohort FASTA assembly.
 */
process SUMMARISE_16S {
    input:
    path barrnap_results

    output:
    path 'best_16S.fna', emit: best_16s
    path '16S_status.tsv', emit: status_table
    path 'all_best_16S.fna', emit: cohort_fasta
    path 'all_best_16S_manifest.tsv', emit: cohort_manifest
    path 'versions.yml', emit: versions

    script:
    '''
    echo "SUMMARISE_16S is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    : > best_16S.fna
    : > 16S_status.tsv
    : > all_best_16S.fna
    : > all_best_16S_manifest.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
