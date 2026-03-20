## HPC Versions DB-Root Symlink Mismatch

### Failure

The tracked OIST `p1` pipeline run completed successfully, but the wrapper
validator failed afterwards with:

```text
ERROR: Database paths do not point at the HPC DB root:
busco_datasets, checkm2_db, codetta_db, eggnog_db, taxdump
```

### Root Cause

`bin/validate_hpc_matrix.py` compared the resolved `--db-root` path against the
raw `image_or_path` strings from `tool_and_db_versions.tsv` using
`startswith(...)`.

On OIST, the operator-facing path may be a storage alias such as `/flash/...`,
while `Path.resolve()` can normalise it to a different real path. That made the
validator reject correct database rows purely because the same location was
spelled through different path aliases.

### Patch

- Normalised each reported database path before comparing it with the resolved
  HPC DB root.
- Switched the check from raw string prefix matching to resolved path ancestry.
- Added a regression that covers one symlinked DB-root alias.

### Verification

- `python3 -m unittest tests.test_validate_hpc_matrix`
