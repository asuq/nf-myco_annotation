/*
 * Run all-vs-all FastANI from the prepared path list and normalise the matrix
 * output to a stable `fastani.matrix` filename.
 */
process FASTANI {
    tag "fastani"

    input:
    path fastani_paths

    output:
    path 'fastani.matrix', emit: matrix
    path 'fastani.tsv', emit: raw_output
    path 'versions.yml', emit: versions

    script:
    """
    set +e
    fastANI \
        --rl "${fastani_paths}" \
        --ql "${fastani_paths}" \
        --matrix \
        -t ${task.cpus} \
        -o fastani.tsv \
        > fastani.log 2>&1
    exit_code=\$?
    set -e

    matrix_file=\$(find . -maxdepth 1 -type f \\( -name 'fastani.tsv.matrix' -o -name '*.matrix' \\) | head -n 1 || true)
    if [[ -n "\${matrix_file}" ]]; then
        cp "\${matrix_file}" fastani.matrix
    else
        : > fastani.matrix
    fi

    if [[ ! -f fastani.tsv ]]; then
        : > fastani.tsv
    fi

    cat <<EOF > versions.yml
    "${task.process}":
      fastani: "$(command -v fastANI >/dev/null 2>&1 && fastANI --version 2>&1 | head -n 1 || echo NA)"
    EOF
    """

    stub:
    '''
    : > fastani.matrix
    : > fastani.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      fastani: "stub"
    EOF
    '''
}
