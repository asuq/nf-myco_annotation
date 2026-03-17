/*
 * Run Prokka for gcode-qualified samples and emit stable file handles for the
 * downstream PADLOC and eggNOG wrappers.
 */
process PROKKA {
    tag "${meta.accession}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/prokka" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(genome), val(gcode)

    output:
    tuple val(meta), path('prokka'), path('prokka.gff'), path('prokka.faa'), path('prokka.log'), emit: results
    tuple val(meta), path('prokka.gff'), path('prokka.faa'), emit: padloc_inputs
    tuple val(meta), path('prokka.faa'), emit: eggnog_inputs
    path 'versions.yml', emit: versions

    script:
    def internalId = (meta.internal_id ?: meta.accession).toString()
    def rawLocustag = internalId.toUpperCase().replaceAll(/[^A-Z0-9]/, '')
    def locustag = rawLocustag ? rawLocustag.take(20) : 'PROKKA'
    if (!(locustag ==~ /^[A-Z].*/)) {
        locustag = "L${locustag}".take(20)
    }
    """
    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    : > prokka.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> prokka.log
        rm -rf prokka
        set +e
        prokka "${genome}" \
            --outdir prokka \
            --prefix "${internalId}" \
            --locustag "${locustag}" \
            --compliant \
            --gcode "${gcode}" \
            --cpus ${task.cpus} \
            --rfam \
            >> prokka.log 2>&1
        exit_code=\$?
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_prokka=%s\n' "\${attempt}" >> prokka.log
        (( attempt += 1 ))
    done

    mkdir -p prokka

    prokka_gff=\$(find prokka -maxdepth 1 -type f -name '*.gff' | head -n 1 || true)
    if [[ "\${exit_code}" -eq 0 && -n "\${prokka_gff}" ]]; then
        cp "\${prokka_gff}" prokka.gff
    else
        : > prokka.gff
    fi

    prokka_faa=\$(find prokka -maxdepth 1 -type f -name '*.faa' | head -n 1 || true)
    if [[ "\${exit_code}" -eq 0 && -n "\${prokka_faa}" ]]; then
        cp "\${prokka_faa}" prokka.faa
    else
        : > prokka.faa
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> prokka.log

    prokka_version="\$(command -v prokka >/dev/null 2>&1 && prokka --version 2>&1 | awk 'NF { value=\$0 } END { if (value) print value }' || true)"
    prokka_version="\${prokka_version:-NA}"
    printf '"%s":\n  prokka: "%s"\n' \
      "${task.process}" \
      "\${prokka_version}" \
      > versions.yml
    """

    stub:
    '''
    mkdir -p prokka
    cat <<'EOF' > prokka.gff
    ##gff-version 3
    contig1	Prokka	CDS	1	30	.	+	0	ID=gnl|Prokka|gene_1
    ##FASTA
    >contig1
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    cat <<'EOF' > prokka.faa
    >sample_a_1
    MAAAAAAAAA
    EOF
    cp prokka.gff prokka/sample_a.gff
    cp prokka.faa prokka/sample_a.faa
    : > prokka.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      prokka: "stub"
    EOF
    '''
}
