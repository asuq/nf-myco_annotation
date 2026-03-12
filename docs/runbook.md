# Runbook

## Required inputs

The pipeline requires:

- `--sample_csv`: manifest CSV with the locked lower-case headers
- `--metadata`: metadata table with `Accession` or `accession`
- `--taxdump`: pinned NCBI taxdump directory containing `names.dmp` and `nodes.dmp`
- `--checkm2_db`: CheckM2 database directory containing one top-level `.dmnd`
- `--busco_db` or `--prepare_busco_datasets true`
- `--eggnog_db` for real eggNOG runs
- `--padloc_db` for real PADLOC runs

Optional labels used only in `tool_and_db_versions.tsv`:

- `--taxdump_label`
- `--checkm2_db_label`
- `--eggnog_db_label`
- `--padloc_db_label`

Shared helper image requirement:

- `params.python_container` is the single helper image for Python-based tasks
  such as validation, ANI clustering, ANI representative selection, final
  status assembly, and version collection
- it must provide `numpy` and `scipy`
- the repo-owned Dockerfile for this image lives under
  `docker/python_helper/Dockerfile`
- local Docker and HPC Singularity runs should reuse that one helper image to
  reduce container pulls

Runtime database prep helper image:

- `prepare_databases.nf` always uses the dedicated helper image
  `quay.io/asuq1617/nf-myco_db:0.2`
- the repo-owned Dockerfile for this image lives under
  `docker/runtime_db_helper/Dockerfile`

## Runtime database preparation

Use `nextflow run prepare_databases.nf` when one or more runtime databases are
missing on the target machine.

- It is a separate operator-facing Nextflow entry point, not part of the normal
  analysis workflow runtime.
- It uses the same database directory params as the main workflow.
- An existing valid database directory is reused in place.
- A missing database directory is populated there when
  `--download_missing_databases true` is set.
- Taxdump uses the shared helper downloader. CheckM2, BUSCO, eggNOG, and PADLOC
  use their own tool-native download or update commands.
- Scratch work uses `--runtime_db_scratch_root` when set and otherwise stays
  beside the destination rather than in system tmp.
- A prepared destination is reusable only when it validates and contains
  `.nf_myco_ready.json`.
- Concurrent preparation is blocked by a per-destination lock sidecar ending in
  `.nf_myco_prepare.lock`.
- The helper logic lives in `bin/prepare_runtime_databases.py`, but that script
  is now the implementation detail behind the prep workflow rather than the
  main operator interface.

Supported runtime databases:

- taxdump
- CheckM2
- BUSCO lineage directories
- eggNOG
- PADLOC

Example:

```bash
nextflow run prepare_databases.nf -profile oist \
  --taxdump /shared/db/ncbi_taxdump_20240914 \
  --checkm2_db /shared/db/checkm2/CheckM2_database \
  --busco_db /shared/db/busco \
  --eggnog_db /shared/db/Eggnog_db/Eggnog_Diamond_db \
  --padloc_db /shared/db/padloc \
  --download_missing_databases true \
  --runtime_db_scratch_root /shared/db/.scratch \
  --outdir /shared/db/runtime-prep
```

The prep workflow writes `runtime_database_report.tsv` and `nextflow_args.txt`
under `--outdir` on success.

## Profiles

- `debug`: composable behaviour profile that defaults eggNOG smoke runs to `GCA_000027325.1`
- `local`: local executor
- `slurm`: SLURM executor with optional `params.slurm_queue` and `params.slurm_cluster_options`
- `singularity`: Singularity execution with optional `params.singularity_cache_dir` and `params.singularity_run_options`
- `oist`: standalone OIST HPC profile with SLURM and Singularity enabled, using the submitting user account
- `test`: local fixture profile for `-stub-run`

## Minimal test path

Use the frozen fixtures under `assets/testdata/stub/`:

```bash
nextflow run . -profile test -stub-run
```

This validates the DSL wiring, channel contracts, publish locations, and final
table/version outputs without requiring external databases or tool containers.

## Manual wrapper

For a single operator-facing entrypoint, use `bin/run_pipeline_test.sh`.

See `docs/run_pipeline_test.md` for wrapper usage, prerequisites, examples, and
expected output locations.

## Acceptance harness

Use `bin/run_acceptance_tests.py` for the layered acceptance workflow:

- `prepare`: download or reuse the tracked acceptance source genomes and build a generated cohort under `assets/testdata/local/acceptance/`
- `unit`: run the Python unit-test layer for fine-grained and minor edge cases
- `stub`: run the full-pipeline `-stub-run` smoke test
- `local`: run the generated positive cohort with `-profile debug,local,docker`
- `slurm`: run the same cohort with `-profile debug,slurm,singularity` and compare its stable outputs against the latest successful local run
- `dbprep-slurm`: run `prepare_databases.nf` on SLURM, download the curated runtime databases, and validate the prepared database tree
- `all`: run `prepare`, `unit`, `stub`, `local`, and `slurm` in sequence

Tracked cohort descriptors live under `assets/testdata/acceptance/`. Large
downloads, generated manifests, and run artefacts stay under the ignored local
tree `assets/testdata/local/acceptance/`.

Real-data `local`, `slurm`, and `all` runs require:

- `--taxdump`
- `--checkm2-db`
- `--busco-db` or `--prepare-busco-datasets`
- `--eggnog-db`
- `--padloc-db`

These acceptance runs also assume the CRISPRCasFinder image is already set via
`params.ccfinder_container` in pipeline config. The harness does not accept a
separate CCFINDER container argument.

Acceptance runs use the `debug` profile by default, which sets
`params.eggnog_only_accessions = 'GCA_000027325.1'`. Normal raw pipeline runs
still execute eggNOG for every gcode-qualified sample unless you opt into
`debug` or set that parameter explicitly. If you override `--local-profile` or
`--slurm-profile` in the harness, include `debug` yourself if you still want
the smoke-only eggNOG behaviour.

For OIST or any other full-eggNOG HPC validation, use raw `nextflow run .`
instead of the default SLURM acceptance wrapper. The wrapper's SLURM mode is
still centred on the debug acceptance cohort and local-baseline comparison.

Example local acceptance run:

```bash
python3 bin/run_acceptance_tests.py local \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db
```

Example SLURM acceptance run:

```bash
python3 bin/run_acceptance_tests.py slurm \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db \
  --slurm-queue short
```

Example SLURM database-prep run:

```bash
python3 bin/run_acceptance_tests.py dbprep-slurm \
  --dbprep-profile oist \
  --work-root /path/to/work-root \
  --taxdump /path/to/db/ncbi_taxdump_20240914 \
  --checkm2-db /path/to/db/checkm2/CheckM2_database \
  --busco-db /path/to/db/busco \
  --eggnog-db /path/to/db/Eggnog_db/Eggnog_Diamond_db \
  --padloc-db /path/to/db/padloc \
  --slurm-queue short \
  --singularity-cache-dir /path/to/singularity-cache
```

## Example real runs

Local:

```bash
nextflow run . \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --outdir results
```

Local debug smoke variant:

```bash
nextflow run . -profile debug,local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --outdir results
```

Local debug run with an overridden eggNOG smoke accession:

```bash
nextflow run . -profile debug,local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --eggnog_only_accessions SOME_OTHER_ACCESSION \
  --outdir results
```

SLURM:

```bash
nextflow run . -profile slurm \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --outdir results
```

Singularity:

```bash
nextflow run . -profile singularity \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --outdir results
```

OIST:

```bash
nextflow run . -profile oist \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --singularity_cache_dir /path/to/singularity-cache \
  --outdir results
```

## OIST HPC testing

Use the OIST profile directly for full eggNOG runs:

```bash
nextflow run . -profile oist ...
```

Do not include `debug` and do not set `--eggnog_only_accessions` when the goal
is to validate full eggNOG execution on HPC.

Build everything on HPC-native storage. Do not copy the local database tree to
HPC. Use three stages in order:

1. prepare the tracked acceptance cohort on the HPC login node
2. prepare runtime databases on SLURM with `dbprep-slurm`
3. run the real pipeline with raw `nextflow run . -profile oist`

Step-by-step procedure:

1. Start a persistent shell on the login node and move into the repository.

```bash
tmux new -s nf_myco_hpc
cd /path/to/nf-myco_annotation
```

2. Define one clean HPC root.

```bash
export HPC_ROOT=/scratch/$USER/nf_myco_hpc_$(date +%Y%m%d)
export ACCEPT_ROOT=$HPC_ROOT/acceptance
export DB_ROOT=$HPC_ROOT/db
export RESULT_ROOT=$HPC_ROOT/results
export SINGULARITY_CACHE=$HPC_ROOT/singularity_cache

mkdir -p "$ACCEPT_ROOT" "$DB_ROOT" "$RESULT_ROOT" "$SINGULARITY_CACHE"
```

3. Run preflight checks.

```bash
python3 --version
nextflow -version
singularity --version
sbatch --version
nextflow config -profile oist >/dev/null
```

4. Prepare the tracked 9-sample cohort on the HPC login node.

```bash
bin/run_pipeline_test.sh prepare \
  --work-root "$ACCEPT_ROOT"
```

Expected result:

- `$ACCEPT_ROOT/generated/sample_sheet.csv`
- `$ACCEPT_ROOT/generated/metadata.tsv`
- `$ACCEPT_ROOT/download_checksums.tsv`
- downloaded FASTA files under `$ACCEPT_ROOT/downloads/`

If this step fails, bring back:

- `$ACCEPT_ROOT/prepare.log` if present
- `$ACCEPT_ROOT/generated/`
- `$ACCEPT_ROOT/download_checksums.tsv`

5. Prepare the runtime databases on SLURM with the wrapper gate.

Set the explicit HPC destination directories first:

```bash
export TAXDUMP_DIR=$DB_ROOT/ncbi_taxdump_20240914
export CHECKM2_DIR=$DB_ROOT/checkm2/CheckM2_database
export BUSCO_DIR=$DB_ROOT/busco
export EGGNOG_DIR=$DB_ROOT/Eggnog_db/Eggnog_Diamond_db
export PADLOC_DIR=$DB_ROOT/padloc
```

```bash
bin/run_pipeline_test.sh dbprep-slurm \
  --dbprep-profile oist \
  --work-root "$ACCEPT_ROOT" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2-db "$CHECKM2_DIR" \
  --busco-db "$BUSCO_DIR" \
  --eggnog-db "$EGGNOG_DIR" \
  --padloc-db "$PADLOC_DIR" \
  --singularity-cache-dir "$SINGULARITY_CACHE"
```

Expected result:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/runtime_database_report.tsv`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/nextflow_args.txt`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/trace.tsv`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/report.html`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/timeline.html`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/dag.html`

Validate the prepared DB directories:

- taxdump:
  - `$TAXDUMP_DIR/names.dmp`
  - `$TAXDUMP_DIR/nodes.dmp`
- CheckM2:
  - exactly one top-level `*.dmnd` in `$CHECKM2_DIR`
- BUSCO:
  - `$BUSCO_DIR/bacillota_odb12/dataset.cfg`
  - `$BUSCO_DIR/mycoplasmatota_odb12/dataset.cfg`
- eggNOG:
  - `$EGGNOG_DIR/eggnog.db`
  - `$EGGNOG_DIR/eggnog_proteins.dmnd`
- PADLOC:
  - `$PADLOC_DIR/hmm/padlocdb.hmm`
- ready markers:
  - `.nf_myco_ready.json` in each prepared DB root

If this step fails, bring back:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/`
- failed task dir under `$ACCEPT_ROOT/runs/dbprep-slurm/work/`
- `.nextflow.log`
- `.command.sh`, `.command.out`, `.command.err`, and `.exitcode` from the failed task

6. Run one optional structural smoke test.

```bash
nextflow run . -profile test -stub-run --outdir "$RESULT_ROOT/stub"
```

7. Run the tracked 9-sample cohort on OIST with full eggNOG.

```bash
nextflow run . -profile oist \
  -work-dir "$RESULT_ROOT/p1/work" \
  --sample_csv "$ACCEPT_ROOT/generated/sample_sheet.csv" \
  --metadata "$ACCEPT_ROOT/generated/metadata.tsv" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --padloc_db "$PADLOC_DIR" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULT_ROOT/p1/out" \
  --max_cpus 64 \
  --max_memory 256.GB \
  --max_time 72.h
```

Do not include `debug` and do not set `--eggnog_only_accessions`.

Monitor:

```bash
squeue -u "$USER"
tail -f .nextflow.log
```

Expected final outputs:

- `$RESULT_ROOT/p1/out/tables/master_table.tsv`
- `$RESULT_ROOT/p1/out/tables/sample_status.tsv`
- `$RESULT_ROOT/p1/out/tables/tool_and_db_versions.tsv`
- `$RESULT_ROOT/p1/out/pipeline_info/trace.tsv`
- `$RESULT_ROOT/p1/out/pipeline_info/report.html`
- `$RESULT_ROOT/p1/out/pipeline_info/timeline.html`
- `$RESULT_ROOT/p1/out/pipeline_info/dag.html`

Acceptance checks:

- `sample_status.tsv` has no failed samples
- `tool_and_db_versions.tsv` has no `unknown` rows
- `eggnog_mapper` is `2.1.13`
- PADLOC tool version is clean
- DB path rows point at `$DB_ROOT/...`
- for full eggNOG runs, `eggnog_status` is `done` for eligible samples and
  `skipped` only when `gcode = NA`

8. Prepare and run a medium real-data cohort of about 20 to 30 samples.

Use your own `sample_csv` and `metadata.tsv`, but make sure the cohort includes:

- gcode4 candidates
- gcode11 candidates
- at least one CRISPR-positive sample
- at least one CRISPR-negative sample
- at least one atypical-excluded sample
- at least one atypical-exception sample
- at least one ANI-near pair
- at least two `is_new=true` rows

Run it with the same command pattern, but use:

- `-work-dir "$RESULT_ROOT/p2/work"`
- `--outdir "$RESULT_ROOT/p2/out"`

9. Run the large or full real-data cohort.

Use the same command pattern again with:

- `-work-dir "$RESULT_ROOT/p3/work"`
- `--outdir "$RESULT_ROOT/p3/out"`

Only increase `--max_cpus`, `--max_memory`, or `--max_time` if the first two
HPC runs show that the defaults are too tight.

10. If a run fails, download the useful artefacts back to local.

For database-prep failures, bring back:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/`
- failed task dir under `$ACCEPT_ROOT/runs/dbprep-slurm/work/`
- `.nextflow.log`

For pipeline failures, bring back:

- `$RESULT_ROOT/<case>/out/tables/`
- `$RESULT_ROOT/<case>/out/pipeline_info/`
- failed task dir under `$RESULT_ROOT/<case>/work/`
- `.nextflow.log`

Minimum useful artefacts:

- `.command.sh`
- `.command.out`
- `.command.err`
- `.exitcode`
- `versions.yml`
- task-specific logs such as `fastani.log`, `ccfinder.log`, or `result.json`

## Final outputs

Published final tables are written under `results/tables/`:

- `master_table.tsv`
- `sample_status.tsv`
- `tool_and_db_versions.tsv`

The master table preserves the original metadata block, then appends derived
columns in the order locked by `assets/master_table_append_columns.txt`.

## Notes

- PADLOC and eggNOG outputs are retained per sample but are not merged into the final master table.
- Original accessions remain the published sample-folder names. Internal sanitized IDs are execution-only.
- The shared Python helper image now includes `numpy` and `scipy` so ANI clustering and representative selection reuse the same helper container as the other Python tasks.
