/*
 * Run PADLOC with the locked inline GFF cleanup so the wrapper keeps a single
 * module boundary around both preparation and execution.
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

    mkdir -p padloc

    set +e
    padloc --faa padloc_input.faa --gff padloc_input.gff --cpu ${task.cpus} --outdir "\$PWD/padloc" ${extraArgs} \
        > padloc.log 2>&1
    exit_code=\$?
    set -e

    printf 'exit_code=%s\n' "\$exit_code" >> padloc.log

    cat <<EOF > versions.yml
    "${task.process}":
      transform: "strip_prokka_prefix_and_fasta_tail"
      padloc: "\$(padloc --help 2>&1 | head -n 1 || echo NA)"
    EOF
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
