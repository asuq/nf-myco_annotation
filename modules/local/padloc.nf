/*
 * Run PADLOC after trimming Prokka-specific prefixes from the staged GFF.
 */
process PADLOC {
    tag "${meta.accession}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/padloc" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(gff), path(faa)

    output:
    tuple val(meta), path('padloc'), path('padloc.log'), emit: results
    path 'versions.yml', emit: versions

    script:
    def extraArgs = (params.padloc_extra_args ?: '').toString()
    """
    awk '
        BEGIN { in_fasta = 0 }
        /^##FASTA/ { in_fasta = 1 }
        in_fasta == 0 {
            gsub(/gnl\\|Prokka\\|/, "")
            print
        }
    ' "${gff}" > padloc_input.gff

    cp "${faa}" padloc_input.faa

    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    : > padloc.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> padloc.log
        rm -rf padloc
        mkdir -p padloc
        set +e
        padloc --faa padloc_input.faa --gff padloc_input.gff --cpu ${task.cpus} --outdir "\$PWD/padloc" ${extraArgs} \
            >> padloc.log 2>&1
        exit_code=\$?
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_padloc=%s\n' "\${attempt}" >> padloc.log
        (( attempt += 1 ))
    done

    printf 'exit_code=%s\n' "\$exit_code" >> padloc.log

    padloc_version="\$(padloc --version 2>&1 | awk 'NF { print; exit }' || echo NA)"
    printf '"%s":\n  transform: "%s"\n  padloc: "%s"\n' \
      "${task.process}" \
      'strip_prokka_prefix_and_fasta_tail' \
      "\${padloc_version}" \
      > versions.yml
    """

    stub:
    '''
    cat <<'EOF' > padloc_input.gff
    ##gff-version 3
    contig1	Prokka	CDS	1	30	.	+	0	ID=gene_1
    EOF
    cat <<'EOF' > padloc_input.faa
    >sample_a_1
    MAAAAAAAAA
    EOF
    mkdir -p padloc
    : > padloc/placeholder.txt
    : > padloc.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      transform: "stub"
      padloc: "stub"
    EOF
    '''
}
