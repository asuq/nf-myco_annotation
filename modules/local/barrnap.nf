/*
 * Run Barrnap per sample, publish the raw Barrnap outputs, and emit the
 * per-sample 16S summary artefacts used by downstream cohort and reporting
 * steps.
 */
process BARRNAP {
    tag "${meta.accession}"
    label 'process_small'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/barrnap" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )
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
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path('rrna.gff'), path('rrna.fa'), path('barrnap.log'), emit: results
    tuple val(meta), path("${meta.internal_id}_best_16S.fna"), path("${meta.internal_id}_16S_status.tsv"), emit: sixteen_s_summaries
    path 'versions.yml', emit: versions

    script:
    def barrnapKingdom = params.barrnap_kingdom ?: 'bac'
    """
    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    : > barrnap.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> barrnap.log
        rm -f rrna.gff rrna.fa
        set +e
        barrnap \
            --threads ${task.cpus} \
            --kingdom ${barrnapKingdom} \
            --outseq rrna.fa \
            "${genome}" \
            > rrna.gff 2>> barrnap.log
        exit_code=\$?
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_barrnap=%s\n' "\${attempt}" >> barrnap.log
        (( attempt += 1 ))
    done

    if [[ "\${exit_code}" -ne 0 ]]; then
        : > rrna.gff
        : > rrna.fa
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> barrnap.log

    summarise_script="\$(command -v summarise_16s.py)"
    python3 "\${summarise_script}" \
        --accession "${meta.accession}" \
        --rrna-gff rrna.gff \
        --rrna-fasta rrna.fa \
        --outdir .

    cp best_16S.fna "${meta.internal_id}_best_16S.fna"
    cp 16S_status.tsv "${meta.internal_id}_16S_status.tsv"

    barrnap_version="\$(command -v barrnap >/dev/null 2>&1 && barrnap --version 2>&1 | awk 'NF { value=\$0 } END { if (value) print value }' || true)"
    barrnap_version="\${barrnap_version:-NA}"
    printf '"%s":\n  barrnap: "%s"\n  python: "%s"\n  summarise_16s: "%s"\n' \
      "${task.process}" \
      "\${barrnap_version}" \
      "\$(python3 --version 2>&1 | sed 's/^Python //')" \
      "bin/summarise_16s.py" \
      > versions.yml
    """

    stub:
    """
    cat <<'EOF' > rrna.gff
    contig1	barrnap	rRNA	1	100	5.0	+	.	Name=16S_rRNA
    EOF
    cat <<'EOF' > rrna.fa
    >sample_a 16S ribosomal RNA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    : > barrnap.log
    cat <<'EOF' > best_16S.fna
    >sample_a 16S ribosomal RNA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    cat <<'EOF' > 16S_status.tsv
    accession	16S	best_16S_header	best_16S_length	warnings
    sample_a	Yes	sample_a 16S ribosomal RNA	80
    EOF
    cp best_16S.fna "${meta.internal_id}_best_16S.fna"
    cp 16S_status.tsv "${meta.internal_id}_16S_status.tsv"
    cat <<'EOF' > versions.yml
    "${task.process}":
      barrnap: "stub"
      python: "stub"
      summarise_16s: "bin/summarise_16s.py"
    EOF
    """
}
