/*
 * Run ANI clustering and representative selection using the existing Python
 * implementation adapted for the pipeline metadata schema.
 */
process CLUSTER_ANI {
    tag "ani_cluster"

    input:
    path ani_matrix
    path ani_metadata

    output:
    path 'cluster.tsv', emit: clusters
    path 'representatives.tsv', emit: representatives
    path 'versions.yml', emit: versions

    script:
    def buscoColumn = params.busco_primary_column ?: "BUSCO_${params.busco_lineages[0]}"
    def aniThreshold = params.ani_threshold ?: 0.95
    """
    python3 "${projectDir}/bin/cluster_ani.py" \
        --ani-matrix "${ani_matrix}" \
        --ani-metadata "${ani_metadata}" \
        --busco-column "${buscoColumn}" \
        --threshold ${aniThreshold} \
        --threads ${task.cpus} \
        --outdir .

    cat <<EOF > versions.yml
    "${task.process}":
      python: "$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/cluster_ani.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > cluster.tsv
    Accession	Cluster_ID	Is_Representative	ANI_to_Representative	Score	Path
    sample_a	cluster_1	yes	100	1.0	/work/sample_a.fna
    EOF
    cat <<'EOF' > representatives.tsv
    Cluster_ID	Representative_Accession	Organism_Name	CheckM2_Completeness	CheckM2_Contamination	BUSCO	Assembly_Level	N50	Cluster_Size
    cluster_1	sample_a	Sample A	95	1	C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200	Scaffold	50000	1
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/cluster_ani.py"
    EOF
    '''
}
