# HPC db-matrix stale case-root mismatch

## Failure

The OIST campaign rerun failed during `db-matrix` with:

`Unexpected dbprep statuses for: checkm2`

The failure appeared after the `db3_valid_without_marker` CheckM2 case completed
successfully.

## Cause

The wrapper reused the same disposable `db_cases/` tree across reruns and only
created directories with `mkdir -p`. That allowed earlier markers and outputs to
persist. On rerun, the `db3_valid_without_marker/checkm2` case was no longer a
true marker-free valid root, so the prep report returned `present` instead of
the wrapper's expected `prepared`.

This was a wrapper-state problem, not a runtime DB-prep failure.

## Fix

Reset the entire disposable `db_cases/` tree at the start of `db-matrix`.

Guard the reset so it only targets a resolved path that:

- is non-empty
- is not `/`
- ends in `/db_cases`

Shared reusable roots such as `acceptance/`, `db/`, `results/`, and
`singularity_cache/` remain untouched.

## Validation

- `python3 -m unittest tests.test_run_oist_hpc_matrix tests.test_validate_hpc_matrix`
- `python3 -m unittest discover -s tests -p 'test_*.py'`
- `bash bin/run_oist_hpc_matrix.sh --dry-run --hpc-root /tmp/hpc_matrix all`
