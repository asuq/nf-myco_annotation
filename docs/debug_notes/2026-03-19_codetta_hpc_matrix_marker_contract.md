## Codetta HPC Matrix Marker Contract

### Failure

The OIST HPC `db-matrix` run failed in the Codetta case that seeds a directory
containing:

- `Pfam-A_enone.hmm`
- `Pfam-A_enone.hmm.h3f`
- `Pfam-A_enone.hmm.h3i`
- `Pfam-A_enone.hmm.h3m`
- `Pfam-A_enone.hmm.h3p`

but no `.nf_myco_ready.json`.

The observed helper error was:

```text
Destination is valid but incomplete for codetta; missing .nf_myco_ready.json
```

### Root Cause

The HPC matrix assumed that Codetta should behave like CheckM2, BUSCO, and
eggNOG in the `db3_valid_without_marker` case. That assumption was wrong.

Codetta is currently prepared through the single-step helper path in
`bin/prepare_runtime_databases.py`, not through the native
`download -> finalise` pattern. The helper intentionally treats a valid Codetta
directory without `.nf_myco_ready.json` as incomplete rather than reusable.

So the failure was a contract mismatch between the OIST matrix wrapper and the
implemented Codetta helper semantics.

### Patch

- Implementation commit:
  - `1f85d0c` `fix(hpc): align codetta db-matrix marker contract`
- Test commit:
  - `4e7a948` `test(hpc): cover codetta marker contract`

Changes made:

- kept Codetta on the existing helper-prepared path
- changed `db3_valid_without_marker/codetta` in
  `bin/run_oist_hpc_matrix.sh` from an expected success to an expected failure
- kept the seeded Codetta files unchanged so the negative case still covers the
  marker-less profile directory
- added helper tests proving that:
  - a valid Codetta directory without `.nf_myco_ready.json` is rejected
  - the same directory succeeds when `--force` is used
- added wrapper and validator regressions so the HPC contract stays aligned

### Deferred Alternative

The broader refactor to split Codetta into `download -> finalise` was deferred.
That would change the runtime-database contract rather than just fixing the HPC
matrix mismatch.

### Verification

- `python3 -m unittest tests.test_prepare_runtime_databases tests.test_run_oist_hpc_matrix tests.test_validate_hpc_matrix`
- `bash -n bin/run_oist_hpc_matrix.sh`
- `bash bin/run_oist_hpc_matrix.sh --dry-run --hpc-root /tmp/nf_myco_hpc db-matrix`
