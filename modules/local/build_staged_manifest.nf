/*
 * Build one deterministic staged-genome manifest TSV for cohort ANI helpers.
 */
process BUILD_STAGED_MANIFEST {
    tag "staged_manifest"
    label 'process_single'
    cache 'deep'

    input:
    val manifest_rows

    output:
    path 'staged_genomes.tsv', emit: manifest
    path 'versions.yml', emit: versions

    script:
    def manifestRowList = manifest_rows instanceof Collection ? manifest_rows : [manifest_rows]
    def manifestRowsText = manifestRowList.join('\n')
    """
    python_path="\$(command -v python3)"
    script_path="\$(command -v build_staged_manifest.py)"

    cat <<'EOF' > staged_manifest_rows.tsv
${manifestRowsText}
EOF

    "\${python_path}" "\${script_path}" \
        --input staged_manifest_rows.tsv \
        --output staged_genomes.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_staged_manifest.py"
    EOF
    """

    stub:
    def manifestRowList = manifest_rows instanceof Collection ? manifest_rows : [manifest_rows]
    def manifestRowsText = manifestRowList.join('\n')
    """
    python_path="\$(command -v python3)"
    script_path="\$(command -v build_staged_manifest.py)"

    cat <<'EOF' > staged_manifest_rows.tsv
${manifestRowsText}
EOF

    "\${python_path}" "\${script_path}" \
        --input staged_manifest_rows.tsv \
        --output staged_genomes.tsv

    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/build_staged_manifest.py"
    EOF
    """
}
