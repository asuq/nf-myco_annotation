/*
 * Decide ANI eligibility, build canonical FastANI input files, and emit the
 * ANI-ready metadata table required by cluster_ani.py.
 */
process BUILD_FASTANI_INPUTS {
    tag "fastani_prep"
    label 'process_single'
    cache 'deep'
    publishDir(
        { "${params.outdir}/cohort/fastani" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['ani_metadata.tsv', 'ani_exclusions.tsv', 'fastani_paths.txt']
                ? filename
                : null
        },
    )

    input:
    path validated_samples
    path metadata
    path staged_manifest
    path staged_fastas
    path checkm2_summaries
    path sixteen_s_statuses
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path assembly_stats
    val primary_busco_column

    output:
    path 'fastani_inputs', emit: fastani_inputs
    path 'fastani_paths.txt', emit: paths
    path 'ani_metadata.tsv', emit: metadata
    path 'ani_exclusions.tsv', emit: exclusions
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    """build_fastani_inputs.py \
    --validated-samples "${validated_samples}" \
    --metadata "${metadata}" \
    --staged-manifest "${staged_manifest}" \
    --checkm2 "${checkm2_summaries}" \
    --16s-status "${sixteen_s_statuses}" \
    ${buscoArgs} \
    --primary-busco-column "${primary_busco_column}" \
    --assembly-stats "${assembly_stats}" \
    --outdir .

    cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/build_fastani_inputs.py"
EOF
"""

    stub:
    """mkdir -p fastani_inputs
cat <<'EOF' > fastani_inputs/sample_a.fasta
>sample_a
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF
cat <<'EOF' > fastani_paths.txt
fastani_inputs/sample_a.fasta
EOF
cat <<'EOF' > ani_metadata.tsv
accession	matrix_name	path	assembly_level	gcode	checkm2_completeness	checkm2_contamination	n50	scaffolds	genome_size	organism_name	${primary_busco_column}
sample_a	fastani_inputs/sample_a.fasta	fastani_inputs/sample_a.fasta	Scaffold	11	95	1	50000	3	800000	Sample_A	C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200
EOF
cat <<'EOF' > ani_exclusions.tsv
accession	internal_id	ani_included	ani_exclusion_reason
sample_a	sample_a	true
EOF
cat <<'EOF' > versions.yml
"${task.process}":
  python: "stub"
  script: "bin/build_fastani_inputs.py"
EOF
"""
}
