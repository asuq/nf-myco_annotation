# HPC `db-full` wrapper idempotence mismatch

## Failure

The OIST matrix run completed `prepare_databases.nf` successfully, but the
wrapper failed immediately afterwards in the `db-full` validation step with:

`Unexpected dbprep statuses for: busco_root, checkm2, eggnog, taxdump`

## Cause

`bin/run_oist_hpc_matrix.sh` always validated `db-full` as a fresh prep and
therefore always expected `prepared` rows from
`runtime_database_report.tsv`.

That assumption is only true on a brand-new `HPC_ROOT`. When the same root is
reused, the database-prep workflow correctly reports `present` for all four
database roots, but the wrapper still validated against `prepared`.

## Fix

Make `db-full` inspect the four canonical ready markers before running the
validator:

- no ready markers -> expect `prepared`
- all four ready markers -> expect `present`
- mixed state -> fail fast with a clear inconsistency error

`db-reuse` remains the explicit second-pass reuse check and still expects
`present`.

## Validation

- `python3 -m unittest tests.test_run_oist_hpc_matrix tests.test_validate_hpc_matrix`
- `python3 -m unittest discover -s tests -p 'test_*.py'`
- `bash bin/run_oist_hpc_matrix.sh --dry-run --hpc-root /tmp/hpc_matrix all`
