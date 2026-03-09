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
- local Docker and HPC Apptainer runs should reuse that one helper image to
  reduce container pulls

## Profiles

- `local`: local executor
- `slurm`: SLURM executor with optional `params.slurm_queue`, `params.slurm_account`, and `params.slurm_cluster_options`
- `apptainer`: Apptainer execution with optional `params.apptainer_cache_dir` and `params.apptainer_run_options`
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
- `local`: run the generated positive cohort with `-profile local,docker`
- `slurm`: run the same cohort with `-profile slurm,apptainer` and compare its stable outputs against the latest successful local run
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

Acceptance runs also short-circuit eggNOG to one tracked smoke accession by
setting `params.eggnog_only_accessions` internally. Normal raw pipeline runs
still execute eggNOG for every gcode-qualified sample unless that param is set
explicitly.

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

Apptainer:

```bash
nextflow run . -profile apptainer \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --busco_download_dir /path/to/busco-lineages \
  --eggnog_db /path/to/eggnog-db \
  --padloc_db /path/to/padloc-db \
  --outdir results
```

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
