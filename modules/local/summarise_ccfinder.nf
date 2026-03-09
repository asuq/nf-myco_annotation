/*
 * Summarise one CRISPRCasFinder JSON output into strain, contig, and CRISPR
 * TSVs. The Python CLI is already failure-tolerant for malformed sample input.
 */
process SUMMARISE_CCFINDER {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/ccfinder" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(ccfinder_dir), path(result_json)

    output:
    tuple val(meta), path('ccfinder_strains.tsv'), path('ccfinder_contigs.tsv'), path('ccfinder_crisprs.tsv'), emit: summaries
    tuple val(meta), path('ccfinder_strains.tsv'), emit: strains
    path 'versions.yml', emit: versions

    script:
    """
    summarise_ccfinder.py \
        --accession "${meta.accession}" \
        --ccfinder-dir "${ccfinder_dir}" \
        --result-json "${result_json}" \
        --outdir .

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/summarise_ccfinder.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > ccfinder_strains.tsv
    accession	CRISPRS	SPACERS_SUM	CRISPR_FRAC	ccfinder_status	warnings
    sample_a	2	7	0.1	done
    EOF
    cat <<'EOF' > ccfinder_contigs.tsv
    accession	contig_id	contig_length	CRISPRS	SPACERS_SUM	CRISPR_FRAC
    sample_a	contig_1	1000	2	7	0.1
    EOF
    cat <<'EOF' > ccfinder_crisprs.tsv
    accession	contig_id	crispr_id	evidence_level	spacer_count	start	end	crispr_length
    sample_a	contig_1	crispr_1	4	3	10	80	71
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/summarise_ccfinder.py"
    EOF
    '''
}
