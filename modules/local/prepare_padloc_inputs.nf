/*
 * Convert Prokka outputs into the PADLOC-safe GFF/FAA pair required by v1.
 */
process PREPARE_PADLOC_INPUTS {
    tag "${meta.accession}"

    input:
    tuple val(meta), path(gff), path(faa)

    output:
    tuple val(meta), path('padloc_input.gff'), path('padloc_input.faa'), emit: inputs
    path 'versions.yml', emit: versions

    script:
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

    cat <<EOF > versions.yml
    "${task.process}":
      transform: "strip_prokka_prefix_and_fasta_tail"
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
    cat <<'EOF' > versions.yml
    "${task.process}":
      transform: "stub"
    EOF
    '''
}
