# nf-myco_annotation

[![Nextflow](https://img.shields.io/badge/version-25.04.8-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

`nf-myco_annotation` is a Nextflow DSL2 pipeline for cohort-scale genome
annotation, QC, taxonomy expansion, and ANI-based reporting.

The implemented workflow validates the input manifest and metadata, stages
genomes to execution-safe internal IDs, runs Barrnap, paired CheckM2
translation-table predictions, BUSCO, Codetta, Prokka, CRISPRCasFinder,
PADLOC, eggNOG, and FastANI-based clustering, then publishes three final
reporting tables:

- `master_table.tsv`
- `sample_status.tsv`
- `tool_and_db_versions.tsv`

This README is the operator-facing landing page for normal setup and runs. For
full operator detail, acceptance workflows, and troubleshooting, see
[`docs/runbook.md`](docs/runbook.md). For the implemented architecture and
reporting contracts, see
[`docs/implemented_pipeline.md`](docs/implemented_pipeline.md).

## Table of contents

- [Project summary](#project-summary)
- [Entrypoints](#entrypoints)
- [Prerequisites](#prerequisites)
- [Required inputs and databases](#required-inputs-and-databases)
- [Quick-start workflow](#quick-start-workflow)
- [Main run examples](#main-run-examples)
- [Output summary](#output-summary)
- [Profiles](#profiles)
- [Documentation map](#documentation-map)


## Project summary

The main analysis workflow is designed around these steps:

- validate the sample manifest and metadata table
- stage each genome under a sanitized internal ID for tool-safe execution
- expand taxonomy from a pinned user-supplied taxdump
- run per-sample QC with Barrnap, paired CheckM2 runs, and BUSCO
- run Codetta for every sample and publish all Codetta artefacts under
  `samples/<accession>/codetta/`
- run gcode-gated annotation with Prokka, CRISPRCasFinder, PADLOC, and eggNOG
- build ANI inputs, run all-vs-all FastANI, cluster genomes, and select
  representatives
- publish final cohort tables and a combined versions report

PADLOC and eggNOG outputs are retained in per-sample folders but are not
merged into `master_table.tsv`. Codetta is merged into the final reporting
tables through `Codetta_Genetic_Code`, `Codetta_NCBI_Table_Candidates`, and
`codetta_status`. In `tool_and_db_versions.tsv`, Codetta is reported as
compatible with `v2.0`, the default container image is
`quay.io/asuq1617/codetta:2.0`, and the pinned upstream `main` commit used for
the image build is recorded explicitly.

## Entrypoints

The repository exposes two operator-facing Nextflow entrypoints.

### `main.nf`

Use `main.nf` for normal analysis runs. It expects the analysis inputs and
runtime databases to be ready before execution, except that BUSCO lineage
datasets can be downloaded through the workflow if
`--prepare_busco_datasets true` is used.

### `prepare_databases.nf`

Use `prepare_databases.nf` before analysis when one or more runtime databases
are missing on the target machine. It prepares taxdump, CheckM2, Codetta,
BUSCO, and eggNOG resources in place and writes
`runtime_database_report.tsv` plus `nextflow_args.txt` under its `--outdir`.

## Prerequisites

- Nextflow `>= 25.04.8`
- either Docker or Singularity for containerised runs
- access to the configured tool images
- a pinned taxdump directory
- a CheckM2 database directory
- a Codetta profile directory containing `Pfam-A_enone.hmm` and its pressed
  HMMER sidecars
- BUSCO lineage datasets under `--busco_db`, or permission to let the workflow
  download them with `--prepare_busco_datasets true`
- an eggNOG data directory for real eggNOG runs

Typical profile combinations are:

- local container run: `-profile local,docker`
- SLURM container run: `-profile slurm,singularity`
- OIST HPC run: `-profile oist`

## Required inputs and databases

The main workflow expects these inputs:

| Input | Required details |
| --- | --- |
| Sample manifest | CSV with the locked lower-case columns `accession`, `is_new`, `assembly_level`, and `genome_fasta` |
| Metadata table | Delimited table with `Accession` or `accession` as the join key |
| Taxdump | Directory containing at least `names.dmp` and `nodes.dmp` |
| CheckM2 database | Directory containing exactly one top-level `.dmnd` file |
| Codetta database | Directory containing `Pfam-A_enone.hmm` plus `.h3f`, `.h3i`, `.h3m`, and `.h3p` |
| BUSCO source | `--busco_db` with lineage directories, or `--prepare_busco_datasets true` |
| eggNOG database | Directory passed through `--eggnog_db` |

Useful defaults from the implementation:

- `--busco_lineages` defaults to `bacillota_odb12,mycoplasmatota_odb12`
- `--gcode_rule` defaults to `strict_delta`
- `--ani_threshold` defaults to `0.95`
- `--outdir` defaults to `results`

## Quick-start workflow

### 1. Optionally prepare runtime databases

Use the dedicated prep entrypoint when one or more databases do not already
exist in their final tool-consumable locations.

```bash
nextflow run prepare_databases.nf -profile local,docker \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --download_missing_databases true \
  --outdir runtime-prep
```

### 2. Run the main analysis workflow

```bash
nextflow run . -profile local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

### 3. Optionally run the stub smoke test

```bash
nextflow run . -profile test -stub-run
```

This validates the DSL wiring, channel contracts, publish locations, and final
report outputs without requiring real external databases or tool containers.

## Main run examples

Local Docker run:

```bash
nextflow run . -profile local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

SLURM plus Singularity run:

```bash
nextflow run . -profile slurm,singularity \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

OIST profile run:

```bash
nextflow run . -profile oist \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --singularity_cache_dir /path/to/singularity-cache \
  --outdir results
```

Main-workflow BUSCO dataset download variant:

```bash
nextflow run . -profile local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --prepare_busco_datasets true \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

## Output summary

By default the workflow publishes under `results/`:

```text
results/
  cohort/
  pipeline_info/
  samples/
  tables/
```

The most important outputs are:

- `results/tables/master_table.tsv`
- `results/tables/sample_status.tsv`
- `results/tables/tool_and_db_versions.tsv`
- `results/samples/<accession>/codetta/` for per-sample Codetta outputs and summaries
- `results/samples/<accession>/...` for the remaining per-sample tool outputs
- `results/cohort/fastani/` for FastANI inputs, metadata, exclusions, and matrix outputs
- `results/cohort/ani_clusters/` for clustering and representative tables
- `results/cohort/taxonomy/` for expanded taxonomy output
- `results/pipeline_info/` for Nextflow trace, report, timeline, and DAG files

For the full published layout and the exact reporting column contracts, see
[`docs/implemented_pipeline.md`](docs/implemented_pipeline.md).

## Profiles

- `local`: local executor
- `docker`: enables Docker execution
- `slurm`: enables the SLURM executor
- `singularity`: enables Singularity execution
- `oist`: OIST HPC profile with SLURM and Singularity enabled
- `test`: fixture-backed local profile for `-stub-run`
- `debug`: composable profile that sets
  `--eggnog_only_accessions GCA_000027325.1`

Use `debug` only when that eggNOG short-circuit is acceptable. For full eggNOG
validation on HPC, use the normal execution profile instead of `debug`.

## Documentation map

- `README.md`: onboarding, prerequisites, quick-start, and normal run examples
- [`docs/runbook.md`](docs/runbook.md): extended operator detail, runtime
  database preparation notes, acceptance workflows, and troubleshooting
- [`docs/implemented_pipeline.md`](docs/implemented_pipeline.md): implemented
  workflow structure, published layout, and final reporting contracts
- [`docs/run_pipeline_test.md`](docs/run_pipeline_test.md): wrapper-specific
  usage for `bin/run_pipeline_test.sh`
