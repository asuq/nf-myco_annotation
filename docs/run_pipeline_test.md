# Manual Pipeline Test Wrapper

## Purpose

`bin/run_pipeline_test.sh` is a thin operator-facing wrapper around
`bin/run_acceptance_tests.py`.

Use it when you want one stable command for manual pipeline testing without
remembering the lower-level Python entrypoint.

The wrapper does not re-implement pipeline logic. It forwards work to the
acceptance harness, which remains the single source of truth for validation,
Nextflow command construction, and output checks.

## Command

```bash
bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|all> [args...]
```

Supported modes:

- `prepare`: download or reuse the tracked acceptance source genomes
- `unit`: run the Python unit-test layer
- `stub`: run the stub smoke test through the acceptance harness
- `local`: run the real-data acceptance cohort with the local profile
- `slurm`: run the real-data acceptance cohort with the SLURM profile and compare against the latest local baseline
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

- `nextflow` for `stub`, `local`, `slurm`, and `all`
- `sbatch` for `slurm` and `all`

Required by real-data modes:

- `--taxdump`
- `--checkm2-db`
- `--busco-download-dir` or `--prepare-busco-datasets`
- `--eggnog-db`
- `--padloc-db`

CRISPRCasFinder still uses `params.ccfinder_container` from `nextflow.config`.
The wrapper does not accept a separate container override.

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
  --busco-download-dir /path/to/busco-lineages \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db
```

Run the local real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh local \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-download-dir /path/to/busco-lineages \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db
```

Run the SLURM real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh slurm \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-download-dir /path/to/busco-lineages \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db \
  --slurm-queue short \
  --slurm-account my_account
```

Run the full layered workflow:

```bash
bin/run_pipeline_test.sh all \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --busco-download-dir /path/to/busco-lineages \
  --eggnog-db /path/to/eggnog-db \
  --padloc-db /path/to/padloc-db
```

## Output locations

The wrapper uses the acceptance harness defaults unless you override them.

- Stub smoke-test results are written under `assets/testdata/local/acceptance/runs/stub/`
- Generated acceptance inputs, downloads, and run artefacts are written under `assets/testdata/local/acceptance/`
- Final published tables for successful runs are under each run directory's `results/tables/`

For the lower-level interface and the raw Nextflow examples, see
`docs/runbook.md`.
