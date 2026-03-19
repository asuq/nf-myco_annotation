# Runbook

## Required inputs

The pipeline requires:

- `--sample_csv`: manifest CSV with the locked lower-case headers
- `--metadata`: metadata table with `Accession` or `accession`
- `--taxdump`: pinned NCBI taxdump directory containing `names.dmp` and `nodes.dmp`
- `--checkm2_db`: CheckM2 database directory containing one top-level `.dmnd`
- `--codetta_db`: Codetta profile directory containing `Pfam-A_enone.hmm` and
  its `.h3f`, `.h3i`, `.h3m`, and `.h3p` sidecars
- `--busco_db` or `--prepare_busco_datasets true`
- `--eggnog_db` for real eggNOG runs

### Preparing metadata.tsv

`accessions.txt` should contain one accession per line. `datasets` and
`dataformat` must be installed and available on `PATH`. This prepares only
`metadata.tsv`; the sample manifest still needs to be prepared separately.

```bash
datasets summary genome accession \
  --inputfile accessions.txt \
  --assembly-source GenBank \
  --assembly-version latest \
  --as-json-lines | \
  dataformat tsv genome \
    --fields accession,assminfo-name,organism-name,assminfo-atypicalwarnings,assminfo-notes,assmstats-total-sequence-len,type_material-label,organism-tax-id,assminfo-level,assminfo-release-date,assminfo-sequencing-tech,assminfo-assembly-method,assmstats-number-of-contigs,assmstats-number-of-scaffolds,assmstats-scaffold-n50,assmstats-genome-coverage,source_database,assminfo-bioproject,assminfo-biosample-accession,assminfo-submitter \
  > metadata.raw.tsv
```

```bash
python3 - <<'PY'
import csv
from pathlib import Path

input_path = Path("metadata.raw.tsv")
output_path = Path("metadata.tsv")

column_map = [
    ("accession", "accession"),
    ("assminfo-name", "Assembly_Name"),
    ("organism-name", "Organism_Name"),
    ("assminfo-atypicalwarnings", "Atypical_Warnings"),
    ("assminfo-notes", "Notes"),
    ("assmstats-total-sequence-len", "Genome_Size"),
    ("type_material-label", "Type_Material"),
    ("organism-tax-id", "Tax_ID"),
    ("assminfo-level", "Assembly_Level"),
    ("assminfo-release-date", "Release_Date"),
    ("assminfo-sequencing-tech", "Sequencing_Tech"),
    ("assminfo-assembly-method", "Assembly_Method"),
    ("assmstats-number-of-contigs", "Contigs"),
    ("assmstats-number-of-scaffolds", "Scaffolds"),
    ("assmstats-scaffold-n50", "N50"),
    ("assmstats-genome-coverage", "Genome_Coverage"),
    ("source_database", "Source_Database"),
    ("assminfo-bioproject", "BioProject"),
    ("assminfo-biosample-accession", "BioSample"),
    ("assminfo-submitter", "Submitter"),
]

with input_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = [source for source, _target in column_map if source not in reader.fieldnames]
    if missing:
        raise SystemExit(
            "metadata.raw.tsv is missing expected column(s): " + ", ".join(missing)
        )

    with output_path.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=[target for _source, target in column_map],
            delimiter="\t",
        )
        writer.writeheader()
        for row in reader:
            writer.writerow({target: row[source] for source, target in column_map})
PY
```

The rename step matters because downstream consumers currently look for
headers such as `Tax_ID`, `Organism_Name`, `Assembly_Level`,
`Atypical_Warnings`, `Notes`, `Genome_Size`, `Scaffolds`, and `N50`.

Optional labels used only in `tool_and_db_versions.tsv`:

- `--taxdump_label`
- `--checkm2_db_label`
- `--codetta_db_label`
- `--eggnog_db_label`

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
  `quay.io/asuq1617/nf-myco_db:0.3`
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
- Taxdump uses the shared helper downloader. CheckM2, BUSCO, and eggNOG use
  their own tool-native download or update commands.
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
- Codetta
- BUSCO lineage directories
- eggNOG

Example:

```bash
nextflow run prepare_databases.nf -profile oist \
  --taxdump /shared/db/ncbi_taxdump_20240914 \
  --checkm2_db /shared/db/checkm2/CheckM2_database \
  --codetta_db /shared/db/codetta/Pfam-A_enone \
  --busco_db /shared/db/busco \
  --eggnog_db /shared/db/Eggnog_db/Eggnog_Diamond_db \
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
expected output locations. Its `dbprep-slurm` and `all` modes now preflight the
repo config so they fail fast if `conf/base.config` still pins the stale
runtime-db helper image tag.

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
- `--codetta-db`
- `--busco-db` or `--prepare-busco-datasets`
- `--eggnog-db`

These acceptance runs also assume the CRISPRCasFinder image is already set via
`params.ccfinder_container` in pipeline config. The harness does not accept a
separate CCFINDER container argument.

PADLOC uses the fixed database bundled in the default PADLOC image, so neither
the pipeline nor the acceptance harness accepts an external PADLOC database
path anymore.

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
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

Example SLURM acceptance run:

```bash
python3 bin/run_acceptance_tests.py slurm \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db \
  --slurm-queue short
```

Example SLURM database-prep run:

```bash
python3 bin/run_acceptance_tests.py dbprep-slurm \
  --dbprep-profile oist \
  --work-root /path/to/work-root \
  --taxdump /path/to/db/ncbi_taxdump_20240914 \
  --checkm2-db /path/to/db/checkm2/CheckM2_database \
  --codetta-db /path/to/db/codetta/Pfam-A_enone \
  --busco-db /path/to/db/busco \
  --eggnog-db /path/to/db/Eggnog_db/Eggnog_Diamond_db \
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
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

Local debug smoke variant:

```bash
nextflow run . -profile debug,local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

Local debug run with an overridden eggNOG smoke accession:

```bash
nextflow run . -profile debug,local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
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
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

Singularity:

```bash
nextflow run . -profile singularity \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

OIST:

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

For the scripted campaign through the medium Mycoplasmatota/Bacillota run, use
the dedicated wrapper:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root all
```

Optional:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root --gcode-rule delta_then_11 all
```

That wrapper:

- prepares the tracked 9-sample cohort on HPC-native storage
- prepares the fixed medium Mycoplasmatota/Bacillota cohort from the
  repo-tracked source catalog and cohort plan under `assets/testdata/medium/`
- runs the full runtime DB gate and the disposable DB existence-state cases
- runs the tracked 9-sample real pipeline gate
- runs the medium Mycoplasmatota/Bacillota real-data gate
- writes the medium run `sample_sheet.csv`, `metadata.tsv`, and
  `source_stats.tsv` under `"$HPC_ROOT/medium/generated"`

It intentionally uses the resource defaults already coded in the OIST profile
and process configuration. It does not pass explicit `--max_cpus`,
`--max_memory`, or `--max_time` overrides.
Nextflow runtime tasks retry up to 3 total attempts. Database download tasks
retry up to 2 total attempts. Sample-level soft-fail modules also retry their
wrapped tools internally up to `soft_fail_attempts` total attempts before
emitting fallback placeholder outputs.

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
which python3
python3 --version
nextflow -version
singularity --version
sbatch --version
nextflow config -profile oist >/dev/null
```

The host-side wrappers and validators require `python3 >= 3.12`. If the login
node still exposes an older interpreter such as `Python 3.6.8`, load a newer
module or adjust `PATH` before you use the wrappers:

```bash
module avail python
module load python/3.12
hash -r
which python3
python3 --version
```

If the login node cannot provide `python3 >= 3.12`, skip the wrappers and run
raw Nextflow directly instead. In that fallback path you lose the wrapper-level
acceptance checks, but the actual pipeline tasks still run inside containers.

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
export CODETTA_DIR=$DB_ROOT/codetta/Pfam-A_enone
export BUSCO_DIR=$DB_ROOT/busco
export EGGNOG_DIR=$DB_ROOT/Eggnog_db/Eggnog_Diamond_db
```

```bash
bin/run_pipeline_test.sh dbprep-slurm \
  --dbprep-profile oist \
  --work-root "$ACCEPT_ROOT" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2-db "$CHECKM2_DIR" \
  --codetta-db "$CODETTA_DIR" \
  --busco-db "$BUSCO_DIR" \
  --eggnog-db "$EGGNOG_DIR" \
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
- Codetta:
  - `$CODETTA_DIR/Pfam-A_enone.hmm`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3f`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3i`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3m`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3p`
- BUSCO:
  - `$BUSCO_DIR/bacillota_odb12/dataset.cfg`
  - `$BUSCO_DIR/mycoplasmatota_odb12/dataset.cfg`
- eggNOG:
  - `$EGGNOG_DIR/eggnog.db`
  - `$EGGNOG_DIR/eggnog_proteins.dmnd`
- ready markers:
  - `.nf_myco_ready.json` in each prepared DB root

If this step fails, bring back:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/`
- failed task dir under `$ACCEPT_ROOT/runs/dbprep-slurm/work/`
- `.nextflow.log`
- `.command.sh`, `.command.out`, `.command.err`, and `.exitcode` from the failed task

Fallback when host Python cannot be fixed:

```bash
nextflow run prepare_databases.nf -profile oist \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --download_missing_databases true \
  --runtime_db_scratch_root "$DB_ROOT/.scratch" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$ACCEPT_ROOT/dbprep_manual"
```

That path still prepares the databases, but it skips the wrapper-level
acceptance validation.

6. Run one optional structural smoke test.

```bash
nextflow run . -profile test -stub-run --outdir "$RESULT_ROOT/stub"
```

7. Run the tracked 9-sample cohort on OIST with full eggNOG.

Do not jump straight to `all` on a new HPC setup. Make this tracked run pass
first, then move on to the full scripted campaign.

```bash
nextflow run . -profile oist \
  -work-dir "$RESULT_ROOT/p1/work" \
  --sample_csv "$ACCEPT_ROOT/generated/sample_sheet.csv" \
  --metadata "$ACCEPT_ROOT/generated/metadata.tsv" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULT_ROOT/p1/out"
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

8. Prepare and run the fixed medium Mycoplasmatota/Bacillota cohort.

The wrapper prepares this cohort automatically from the repo-tracked medium
catalogue and plan before `p2` and `all`, or you can run that step directly:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root "$HPC_ROOT" medium-prepare
```

The generated medium inputs are written under:

- `"$HPC_ROOT/medium/generated/sample_sheet.csv"`
- `"$HPC_ROOT/medium/generated/metadata.tsv"`
- `"$HPC_ROOT/medium/generated/source_stats.tsv"`
- `"$HPC_ROOT/medium/download_checksums.tsv"`

Then run the medium case with:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root "$HPC_ROOT" p2
```

Only after `prepare`, `dbprep-slurm`, and the tracked `p1` run succeed should
you launch the full scripted campaign:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root "$HPC_ROOT" all
```

9. If a run fails, download the useful artefacts back to local.

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
Codetta contributes `Codetta_Genetic_Code` and
`Codetta_NCBI_Table_Candidates`, and `sample_status.tsv` includes
`codetta_status`.

## Notes

- PADLOC and eggNOG outputs are retained per sample but are not merged into the final master table.
- Original accessions remain the published sample-folder names. Internal sanitized IDs are execution-only.
- The shared Python helper image now includes `numpy` and `scipy` so ANI clustering and representative selection reuse the same helper container as the other Python tasks.
- Codetta provenance is split deliberately: `tool_and_db_versions.tsv` reports
  Codetta as `v2.0`, uses the default image
  `quay.io/asuq1617/codetta:2.0`, and records the pinned upstream `main`
  source commit used for that image build.
