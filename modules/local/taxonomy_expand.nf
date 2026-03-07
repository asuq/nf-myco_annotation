/*
 * Expand metadata Tax_ID values into lineage columns from a pinned user-supplied
 * taxdump. Missing or invalid Tax_ID values are left for downstream status
 * handling by omission from the keyed output table.
 */
process TAXONOMY_EXPAND {
    tag "taxonomy"
    label 'process_single'
    publishDir(
        { "${params.outdir}/cohort/taxonomy" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    path validated_samples
    path metadata
    path taxdump

    output:
    path 'taxonomy_expanded.tsv', emit: taxonomy
    path 'versions.yml', emit: versions

    script:
    """
    python3 "${projectDir}/bin/taxonomy_expand.py" \
        --validated-samples "${validated_samples}" \
        --metadata "${metadata}" \
        --taxdump "${taxdump}" \
        --output taxonomy_expanded.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/taxonomy_expand.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > taxonomy_expanded.tsv
    Tax_ID	superkingdom	phylum	class	order	family	genus	species
    123	Bacteria	Firmicutes	Mollicutes	Mycoplasmatales	Mycoplasmataceae	Mycoplasma	Mycoplasma_testus
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/taxonomy_expand.py"
    EOF
    '''
}
