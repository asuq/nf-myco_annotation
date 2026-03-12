# Manual Pipeline Test Wrapper

## Purpose

`bin/run_pipeline_test.sh` is a thin operator-facing wrapper around
`bin/run_acceptance_tests.py`.

Use it when you want one stable command for manual pipeline testing without
remembering the lower-level Python entrypoint.

The wrapper does not re-implement pipeline logic. It forwards work to the
acceptance harness, which remains the single source of truth for validation,
Nextflow command construction, and output checks.

For runtime database preparation on HPC, prefer the wrapper's `dbprep-slurm`
mode. The wrapper's `prepare` mode only prepares the tracked acceptance cohort
and its cached source genomes.

Use raw `nextflow run prepare_databases.nf` only when you need direct operator
control outside the standard SLURM validation gate.

For a scripted SLURM validation of runtime database preparation itself, use the
wrapper's `dbprep-slurm` mode. That mode runs `prepare_databases.nf`, downloads
the curated runtime databases, and validates the prepared database tree.

For OIST or any other full-eggNOG HPC validation, use this sequence:

1. `bin/run_pipeline_test.sh prepare` on the login node
2. `bin/run_pipeline_test.sh dbprep-slurm` on SLURM
3. raw `nextflow run . -profile oist` for the real pipeline run

For the broader OIST campaign through the medium real-data case, use the
dedicated wrapper:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root all \
  --medium-sample-csv /path/to/medium_sample_sheet.csv \
  --medium-metadata /path/to/medium_metadata.tsv
```

That wrapper keeps using the coded OIST resource defaults. It does not add
`--max_cpus`, `--max_memory`, or `--max_time` overrides.

The wrapper's normal `slurm` path remains centred on the debug acceptance
cohort and local-baseline comparison.

## Command

```bash
bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|dbprep-slurm|all> [args...]
```

Supported modes:

- `prepare`: download or reuse the tracked acceptance source genomes
- `unit`: run the Python unit-test layer
- `stub`: run the stub smoke test through the acceptance harness
- `local`: run the real-data acceptance cohort with the local profile
- `slurm`: run the real-data acceptance cohort with the SLURM profile and compare against the latest local baseline
- `dbprep-slurm`: run the runtime database prep workflow on SLURM and validate the prepared database tree
- `all`: run `prepare`, `unit`, `stub`, `local`, and `slurm` in sequence

Wrapper options:

- `--dry-run`: print the delegated `python3 bin/run_acceptance_tests.py ...` command and exit
- `-h`, `--help`: show wrapper help before a mode is chosen

Arguments after the mode are forwarded unchanged to
`bin/run_acceptance_tests.py`.

If you need mode-specific help, run:

```bash
bin/run_pipeline_test.sh <mode> --help
```

## Prerequisites

Always required:

- `bash`
- `python3`

Required by pipeline-running modes:

- `nextflow` for `stub`, `local`, `slurm`, `dbprep-slurm`, and `all`
- `sbatch` for `slurm`, `dbprep-slurm`, and `all`

Required by real-data modes:

- `--taxdump`
- `--checkm2-db`
- `--busco-db` or `--prepare-busco-datasets`
- `--eggnog-db`

CRISPRCasFinder still uses `params.ccfinder_container` from `nextflow.config`.
The wrapper does not accept a separate container override.

The wrapper does not prepare runtime databases such as taxdump, CheckM2,
BUSCO, or eggNOG. Prepare those separately before real-data runs when they are
not already available. PADLOC uses the fixed database bundled in the default
PADLOC image.

Acceptance-backed `local`, `slurm`, and `all` runs now use the composable
`debug` profile by default. That profile restricts eggNOG to the tracked smoke
accession `GCA_000027325.1`. Raw pipeline runs are unchanged and still execute
eggNOG for every gcode-qualified sample unless you opt into `debug` or set
`params.eggnog_only_accessions` explicitly. If you override the harness
profiles, include `debug` yourself if you still want the smoke-only eggNOG
behaviour.

## Examples

Quick stub smoke test:

```bash
bin/run_pipeline_test.sh stub
```

Prepare the cached real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh prepare
```

Preview a real local run without executing it:

```bash
bin/run_pipeline_test.sh --dry-run local \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

Run the local real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh local \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

Run the SLURM real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh slurm \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db \
  --slurm-queue short
```

Run the SLURM runtime database-prep validation:

```bash
bin/run_pipeline_test.sh dbprep-slurm \
  --dbprep-profile oist \
  --work-root /path/to/work-root \
  --taxdump /path/to/db/ncbi_taxdump_20240914 \
  --checkm2-db /path/to/db/checkm2/CheckM2_database \
  --busco-db /path/to/db/busco \
  --eggnog-db /path/to/db/Eggnog_db/Eggnog_Diamond_db \
  --slurm-queue short \
  --singularity-cache-dir /path/to/singularity-cache
```

Run the full layered workflow:

```bash
bin/run_pipeline_test.sh all \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

## Output locations

The wrapper uses the acceptance harness defaults unless you override them.

- Stub smoke-test results are written under `assets/testdata/local/acceptance/runs/stub/`
- Generated acceptance inputs, downloads, and run artefacts are written under `assets/testdata/local/acceptance/`
- Final published tables for successful runs are under each run directory's `results/tables/`

For the lower-level interface and the raw Nextflow examples, see
`docs/runbook.md`.
