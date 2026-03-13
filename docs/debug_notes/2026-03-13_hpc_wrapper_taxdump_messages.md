# HPC taxdump negative-case wrapper mismatch

## Failure

The OIST matrix run reached the negative taxdump case
`db5_taxdump_missing_no_download` and then failed in the wrapper with:

`Failure log ... does not contain expected text: Taxdump destination is missing and download is disabled`

## Cause

Taxdump prep no longer emits the old module-local taxdump error strings. It now
uses the shared helper in `bin/prepare_runtime_databases.py`, which reports:

- `No local source was supplied for taxdump, and remote download is disabled.`
- `Destination must be a directory for taxdump: ...`

The wrapper still expected the earlier taxdump-specific wording, so the
negative-case validation failed even though the workflow behaviour itself was
correct.

## Fix

Align only the wrapper-side expected strings for the taxdump negative cases:

- missing-without-download
- file-instead-of-directory

No runtime DB-prep logic changed.

## Validation

- `python3 -m unittest tests.test_run_oist_hpc_matrix tests.test_validate_hpc_matrix`
- `python3 -m unittest discover -s tests -p 'test_*.py'`
- `bash bin/run_oist_hpc_matrix.sh --dry-run --hpc-root /tmp/hpc_matrix all`
