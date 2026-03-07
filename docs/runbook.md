# Runbook

## Required inputs

The pipeline requires:

- `--sample_csv`: manifest CSV with the locked lower-case headers
- `--metadata`: metadata table with `Accession` or `accession`
- `--taxdump`: pinned NCBI taxdump directory containing `names.dmp` and `nodes.dmp`
- `--checkm2_db`: CheckM2 database path
- `--busco_download_dir` or `--prepare_busco_datasets true`
- `--eggnog_db` for real eggNOG runs

Optional labels used only in `tool_and_db_versions.tsv`:

- `--taxdump_label`
- `--checkm2_db_label`
- `--eggnog_db_label`
- `--padloc_db_label`

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
- Custom process images are intentionally configured by placeholder params; manual image pinning remains user-controlled.
