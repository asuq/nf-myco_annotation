/*
 * Run Barrnap per sample and always emit parseable output files for downstream
 * summarisation, even when Barrnap itself fails for that sample.
 */
process BARRNAP {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path('rrna.gff'), path('rrna.fa'), path('barrnap.log'), emit: results
    path 'versions.yml', emit: versions

    script:
    def barrnapKingdom = params.barrnap_kingdom ?: 'bac'
    """
    set +e
    barrnap \
        --threads ${task.cpus} \
        --kingdom ${barrnapKingdom} \
        --outseq rrna.fa \
        "${genome}" \
        > rrna.gff 2> barrnap.log
    exit_code=\$?
    set -e

    if [[ \$exit_code -ne 0 ]]; then
        : > rrna.gff
        : > rrna.fa
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> barrnap.log

    cat <<EOF > versions.yml
    "${task.process}":
      barrnap: "$(command -v barrnap >/dev/null 2>&1 && barrnap --version 2>&1 | head -n 1 || echo NA)"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > rrna.gff
    contig1	barrnap	rRNA	1	100	5.0	+	.	Name=16S_rRNA
    EOF
    cat <<'EOF' > rrna.fa
    >sample_a 16S ribosomal RNA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    : > barrnap.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      barrnap: "stub"
    EOF
    '''
}
