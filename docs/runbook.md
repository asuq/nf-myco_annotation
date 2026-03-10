# Runbook

## Required inputs

The pipeline requires:

- `--sample_csv`: manifest CSV with the locked lower-case headers
- `--metadata`: metadata table with `Accession` or `accession`
- `--taxdump`: pinned NCBI taxdump directory containing `names.dmp` and `nodes.dmp`
- `--checkm2_db`: CheckM2 database path
- `--busco_download_dir` or `--prepare_busco_datasets true`
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

- `params.runtime_db_helper_container` is the optional dedicated helper image
  for `prepare_databases.nf`
- the repo-owned Dockerfile for this image lives under
  `docker/runtime_db_helper/Dockerfile`
- if `params.runtime_db_helper_container` is unset, the prep workflow runs on
  the host and requires `aria2c` plus Python 3.12 on `PATH`

## Runtime database preparation

Use `nextflow run prepare_databases.nf` when one or more runtime databases are
missing on the target machine.

- It is a separate operator-facing Nextflow entry point, not part of the normal
  analysis workflow runtime.
- It prepares databases in place under one canonical `--db_root`.
- It can reuse local staged sources, download curated remote sources, or mix the
  two in one run.
- Archive extraction uses destination-local scratch by default, not system tmp.
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
  --db_root /shared/db/runtime \
  --taxdump_source /staged/db_sources/taxdump_20240914 \
  --checkm2_source /staged/db_sources/checkm2/CheckM2_database.dmnd \
  --busco_source_root /staged/db_sources/busco \
  --eggnog_source /staged/db_sources/eggnog_data \
  --padloc_source /staged/db_sources/padloc_data \
  --runtime_db_link_mode symlink \
  --outdir /shared/db/runtime-prep
```

The prep workflow writes `runtime_database_report.tsv` and `nextflow_args.txt`
under `--outdir` on success.

## Profiles

- `debug`: composable behaviour profile that defaults eggNOG smoke runs to `GCA_000027325.1`
- `local`: local executor
- `slurm`: SLURM executor with optional `params.slurm_queue`, `params.slurm_account`, and `params.slurm_cluster_options`
- `singularity`: Singularity execution with optional `params.singularity_cache_dir` and `params.singularity_run_options`
- `oist`: standalone OIST HPC profile with SLURM and Singularity enabled
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
- `all`: run `prepare`, `unit`, `stub`, `local`, and `slurm` in sequence

Tracked cohort descriptors live under `assets/testdata/acceptance/`. Large
downloads, generated manifests, and run artefacts stay under the ignored local
tree `assets/testdata/local/acceptance/`.

Real-data `local`, `slurm`, and `all` runs require:

- `--taxdump`
- `--checkm2-db`
- `--busco-download-dir` or `--prepare-busco-datasets`
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
  --busco-download-dir /path/to/busco-lineages \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db
```

Example SLURM acceptance run:

```bash
python3 bin/run_acceptance_tests.py slurm \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-download-dir /path/to/busco-lineages \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db \
  --slurm-queue short \
  --slurm-account my_account
```

## Example real runs

Local:

```bash
nextflow run . \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_download_dir /path/to/busco-lineages \
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
  --busco_download_dir /path/to/busco-lineages \
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
  --busco_download_dir /path/to/busco-lineages \
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
  --busco_download_dir /path/to/busco-lineages \
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
  --busco_download_dir /path/to/busco-lineages \
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
  --busco_download_dir /path/to/busco-lineages \
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

Step-by-step procedure:

1. Start a persistent shell on the login node and move into the repository.

```bash
tmux new -s nf_myco_hpc
cd /path/to/nf-myco_annotation
```

2. Define one clean test root and the shared working paths.

```bash
export HPC_ROOT=/scratch/$USER/nf_myco_hpc_$(date +%Y%m%d)
export DB_SRC_ROOT=/shared/staged_db_sources
export DB_TEST_ROOT=$HPC_ROOT/db_prep_tests
export DB_RUNTIME_ROOT=$HPC_ROOT/db_runtime
export WORK_ROOT=$HPC_ROOT/work
export RESULTS_ROOT=$HPC_ROOT/results
export SINGULARITY_CACHE=$HPC_ROOT/singularity_cache

mkdir -p "$DB_TEST_ROOT" "$DB_RUNTIME_ROOT" "$WORK_ROOT" "$RESULTS_ROOT" "$SINGULARITY_CACHE"
```

3. Run preflight checks.

```bash
python3 --version
nextflow -version
singularity --version
sbatch --version
nextflow config -profile oist >/dev/null
```

4. Test runtime database preparation from directory sources.

```bash
nextflow run prepare_databases.nf -profile oist \
  --db_root "$DB_TEST_ROOT/db1" \
  --taxdump_source "$DB_SRC_ROOT/taxdump_20240914" \
  --checkm2_source "$DB_SRC_ROOT/checkm2/CheckM2_database.dmnd" \
  --busco_source_root "$DB_SRC_ROOT/busco" \
  --eggnog_source "$DB_SRC_ROOT/eggnog_data" \
  --padloc_source "$DB_SRC_ROOT/padloc_data" \
  --runtime_db_link_mode symlink \
  --runtime_db_scratch_root "$DB_TEST_ROOT/db1/.scratch" \
  --outdir "$DB_TEST_ROOT/db1/out"
```

Expected result:

- `runtime_database_report.tsv` shows `prepared`
- every destination contains `.nf_myco_ready.json`
- `nextflow_args.txt` is ready to paste into the main workflow

5. Re-run the exact same command from step 4.

Expected result:

- report rows become `present`

6. Test archive-source handling in a disposable destination.

```bash
mkdir -p "$DB_TEST_ROOT/archives"
tar -C "$DB_SRC_ROOT" -czf "$DB_TEST_ROOT/archives/taxdump.tar.gz" taxdump_20240914
tar -C "$DB_SRC_ROOT/busco" -czf "$DB_TEST_ROOT/archives/bacillota_odb12.tar.gz" bacillota_odb12
tar -C "$DB_SRC_ROOT/busco" -czf "$DB_TEST_ROOT/archives/mycoplasmatota_odb12.tar.gz" mycoplasmatota_odb12
tar -C "$DB_SRC_ROOT" -czf "$DB_TEST_ROOT/archives/eggnog.tar.gz" eggnog_data
tar -C "$DB_SRC_ROOT" -czf "$DB_TEST_ROOT/archives/padloc.tar.gz" padloc_data

nextflow run prepare_databases.nf -profile oist \
  --db_root "$DB_TEST_ROOT/db2" \
  --taxdump_source "$DB_TEST_ROOT/archives/taxdump.tar.gz" \
  --checkm2_source "$DB_SRC_ROOT/checkm2/CheckM2_database.dmnd" \
  --busco_source_root "$DB_TEST_ROOT/archives" \
  --eggnog_source "$DB_TEST_ROOT/archives/eggnog.tar.gz" \
  --padloc_source "$DB_TEST_ROOT/archives/padloc.tar.gz" \
  --runtime_db_link_mode symlink \
  --runtime_db_scratch_root "$DB_TEST_ROOT/db2/.scratch" \
  --outdir "$DB_TEST_ROOT/db2/out"
```

Expected result:

- preparation succeeds
- no `.prepare-*` directories remain after success

7. Test invalid-destination failure and `--force` recovery.

```bash
mkdir -p "$DB_TEST_ROOT/db3/padloc"
printf 'broken\n' > "$DB_TEST_ROOT/db3/padloc/broken.txt"

nextflow run prepare_databases.nf -profile oist \
  --db_root "$DB_TEST_ROOT/db3" \
  --padloc_source "$DB_SRC_ROOT/padloc_data" \
  --outdir "$DB_TEST_ROOT/db3/out"
```

Expected result:

- non-zero exit because the destination is invalid

Then rerun with `--force`:

```bash
nextflow run prepare_databases.nf -profile oist \
  --db_root "$DB_TEST_ROOT/db3" \
  --padloc_source "$DB_SRC_ROOT/padloc_data" \
  --force_runtime_database_rebuild true \
  --outdir "$DB_TEST_ROOT/db3/out"
```

Expected result:

- success
- `broken.txt` is gone

8. Prepare the final runtime database root that all real HPC runs will share.

```bash
nextflow run prepare_databases.nf -profile oist \
  --db_root "$DB_RUNTIME_ROOT" \
  --taxdump_source "$DB_SRC_ROOT/taxdump_20240914" \
  --checkm2_source "$DB_SRC_ROOT/checkm2/CheckM2_database.dmnd" \
  --busco_source_root "$DB_SRC_ROOT/busco" \
  --eggnog_source "$DB_SRC_ROOT/eggnog_data" \
  --padloc_source "$DB_SRC_ROOT/padloc_data" \
  --runtime_db_link_mode symlink \
  --runtime_db_scratch_root "$DB_RUNTIME_ROOT/.scratch" \
  --outdir "$DB_RUNTIME_ROOT/out"
```

9. Run one optional structural smoke test.

```bash
nextflow run . -profile test -stub-run --outdir "$RESULTS_ROOT/stub"
```

10. Generate the tracked 9-sample cohort from the locked acceptance assets.

```bash
python3 bin/run_acceptance_tests.py prepare --work-root "$WORK_ROOT/p1"
```

This creates the sample sheet and metadata from the tracked cohort plan under
`assets/testdata/acceptance/`.

11. Run the tracked 9-sample cohort on OIST with full eggNOG.

```bash
nextflow run . -profile oist \
  -work-dir "$WORK_ROOT/p1/runs/full_eggnog/work" \
  --sample_csv "$WORK_ROOT/p1/generated/sample_sheet.csv" \
  --metadata "$WORK_ROOT/p1/generated/metadata.tsv" \
  --taxdump "$DB_RUNTIME_ROOT/taxdump_20240914" \
  --checkm2_db "$DB_RUNTIME_ROOT/checkm2" \
  --busco_download_dir "$DB_RUNTIME_ROOT/busco" \
  --eggnog_db "$DB_RUNTIME_ROOT/eggnog" \
  --padloc_db "$DB_RUNTIME_ROOT/padloc" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULTS_ROOT/p1" \
  --max_cpus 64 \
  --max_memory 256.GB \
  --max_time 72.h
```

12. Prepare and run a medium real-data cohort of about 20 to 30 samples.

Use your own `sample_csv` and `metadata.tsv`, but make sure the cohort includes:

- gcode4 candidates
- gcode11 candidates
- at least one CRISPR-positive sample
- at least one CRISPR-negative sample
- at least one atypical-excluded sample
- at least one atypical-exception sample
- at least one ANI-near pair
- at least two `is_new=true` rows

Run it with the same command pattern, but use a fresh `-work-dir` and
`--outdir`.

13. Run the large or full real-data cohort.

Use the same command pattern again with its own `-work-dir` and `--outdir`.
Increase `--max_time`, `--max_memory`, or `--slurm_queue` only if your site
policy or the earlier runs show that the defaults are too tight.

14. Monitor all real runs.

```bash
squeue -u "$USER"
tail -f .nextflow.log
```

Full-eggNOG success criteria:

- `results/tables/master_table.tsv`, `sample_status.tsv`, and `tool_and_db_versions.tsv` exist
- `results/pipeline_info/trace.tsv`, `report.html`, `timeline.html`, and `dag.html` exist
- `tool_and_db_versions.tsv` points at the prepared runtime database root
- `eggnog_status` is `done` for gcode-qualified samples and `skipped` only when `gcode = NA`
- eligible sample folders contain fresh eggNOG, PADLOC, and CCFINDER outputs

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
