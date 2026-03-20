# 2026-03-20 medium HPC `gcode_na` validation mismatch

## Failure

The OIST `p2` medium campaign completed successfully in Nextflow, but the
wrapper failed afterwards with:

```text
ERROR: Found failed sample-status rows: ...
```

The affected rows were all Bacillota medium fixtures whose
`sample_status.tsv` entries had:

- `gcode_status=failed`
- `gcode=NA`
- `warnings` containing `gcode_na`
- downstream annotation statuses set to `skipped`, not `failed`

## Root cause

`bin/validate_hpc_matrix.py` used the same blanket `assert_no_failed_statuses()`
rule for both the tracked `p1` run and the fixed medium `p2` run. That rule
rejects any row containing the literal value `failed`.

This was too strict for the medium cohort. The fixed medium fixtures include
Bacillota cases that can remain unresolved under the pipeline's default
`strict_delta` gcode policy. Those rows are expected to stop at gcode
assignment, emit `gcode_na`, and skip the gcode-dependent annotation stages.

## Patches

- `468aacb` `fix(hpc): allow isolated medium gcode_na rows`
  - added a medium-specific sample-status validator
  - kept `p1` on the original strict rule
  - allowed only the isolated pattern where `gcode_status` is the sole failed
    column and `warnings` contains `gcode_na`
- `7ae6751` `test(hpc): cover medium gcode_na validation`
  - added a `medium-run` acceptance test for the allowed isolated `gcode_na`
    case
  - added rejection tests for extra failed statuses and missing `gcode_na`
    warnings
  - kept explicit coverage that the strict helper still rejects any failed row

## Verification

- `python3 -m unittest tests.test_validate_hpc_matrix`
  - `Ran 11 tests`
  - `OK`
- `p1` behaviour remains unchanged because `validate_tracked_run()` still calls
  `assert_no_failed_statuses()`
