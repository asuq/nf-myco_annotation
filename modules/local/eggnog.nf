/*
 * Placeholder for eggNOG annotation.
 */
process EGGNOG {
    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path('eggnog'), emit: outdir
    tuple val(meta), path('eggnog_annotations.tsv'), emit: annotations
    path 'versions.yml', emit: versions

    script:
    '''
    echo "EGGNOG is a placeholder module." >&2
    exit 1
    '''

    stub:
    '''
    mkdir -p eggnog
    : > eggnog_annotations.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      placeholder: "true"
    EOF
    '''
}
