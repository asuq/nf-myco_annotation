## Codetta Singularity Input Collision

### Failure

The OIST HPC `db-full` run failed in `PREP_CODETTA_DATABASE` with:

```text
Process input file name collision -- There are multiple input files for each
of the following file names: codetta
```

### Root Cause

`prepare_databases.nf` stages two path inputs into the Codetta helper process:

- the destination parent directory
- the helper scratch root

When `params.runtime_db_scratch_root` is unset, the workflow previously reused
the destination parent itself as the staged scratch root. For Codetta on HPC,
both staged paths therefore had the same basename, `codetta`, which is illegal
for Nextflow process input staging under Singularity.

### Patch

- Implementation commits:
  - `b111dea` `fix(dbprep): avoid codetta input name collisions`
  - `78fff8e` `test(dbprep): cover singularity scratch staging layout`
  - `6c3ffe1` `docs(debug): record codetta singularity collision fix`
- The default helper scratch root now lives under the destination parent as:
  `.nf_myco_runtime_db_scratch`
- If an explicit `runtime_db_scratch_root` resolves to the same directory as
  the destination parent, the workflow also falls back to that internal scratch
  directory.
- This preserves the same-filesystem requirement for rename-based finalisation
  while avoiding staged input basename collisions.

### Verification

- `python3 -m unittest tests.test_nextflow_module_syntax`
- `nextflow run prepare_databases.nf -profile test -stub-run --taxdump ... --codetta_db ...`
