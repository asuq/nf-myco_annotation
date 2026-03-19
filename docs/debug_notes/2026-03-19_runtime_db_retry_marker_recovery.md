## Runtime DB Retry Marker Recovery

### Failure

The OIST HPC `db-full` run retried the native download modules, but the second
attempt then failed with messages such as:

```text
CheckM2 destination exists but is not ready: checkm2/CheckM2_database
```

The same pattern appeared for BUSCO and eggNOG.

### Root Cause

The native download modules stage the destination parent into the task work
directory and write directly into the final host destination. If a first
download attempt fails after creating a partial directory, Nextflow retries the
process, but the second attempt previously treated that partial directory as a
user-owned invalid destination and failed unless
`--force_runtime_database_rebuild` was set.

That behaviour is correct for genuinely pre-existing invalid destinations, but
it is wrong for retryable download failures created by the same workflow run.

### Patch

- Implementation commits:
  - `444a6dc` `fix(dbprep): recover native downloads on retry`
  - `5217115` `test(dbprep): cover native download retry markers`
  - `0253b29` `docs(debug): record native download retry recovery`
- Added per-component retry markers for:
  - CheckM2
  - BUSCO
  - eggNOG
- The marker is created when a native download starts.
- If a later retry sees an invalid partial destination and the matching retry
  marker is present, the module now removes the partial destination and retries
  the download automatically.
- If the invalid destination existed before this workflow run and no retry
  marker is present, the module still fails without
  `--force_runtime_database_rebuild`, preserving the safety contract.

### Verification

- `python3 -m unittest tests.test_nextflow_module_syntax`
- `nextflow run prepare_databases.nf -profile test -stub-run --taxdump ... --checkm2_db ... --codetta_db ... --busco_db ... --eggnog_db ...`
