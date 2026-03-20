# 2026-03-21 medium HPC secondary BUSCO validation mismatch

## Failure

The OIST `p2` medium campaign completed successfully in Nextflow, but the
wrapper failed afterwards with:

```text
ERROR: Found unexpected medium sample-status failures:
GCF_000009045.1_MED1 (busco_mycoplasmatota_odb12_status)
```

The affected row showed an isolated failure in the secondary BUSCO lineage for
a Bacillota sample.

## Root cause

`bin/validate_hpc_matrix.py` had already been relaxed for isolated
`gcode_status=failed` rows with `gcode_na`, but it still rejected every BUSCO
failure. In the medium campaign this was too strict for degraded failures in
the non-primary BUSCO lineage.

The pipeline uses the first configured BUSCO lineage, currently
`bacillota_odb12`, as the primary lineage for ANI eligibility and
representative scoring. A failure in the secondary lineage can still leave the
sample operationally valid for the medium acceptance gate as long as the
primary lineage succeeded.

## Patches

- `eb82f81` `fix(hpc): allow secondary busco failures in p2`
  - added BUSCO-lineage-aware logic to the medium validator
  - continued to reject any primary BUSCO lineage failure
  - allowed only isolated secondary BUSCO failures that record
    `busco_summary_failed`
- `6ff50af` `test(hpc): cover secondary busco medium validation`
  - added an accepted `medium-run` case for an isolated secondary BUSCO
    failure
  - added rejection coverage for a primary BUSCO failure
  - added rejection coverage for a secondary BUSCO failure without the expected
    warning token

## Verification

- `python3 -m unittest tests.test_validate_hpc_matrix`
  - `Ran 14 tests`
  - `OK`
- `p1` remains strict because tracked-run validation still uses
  `assert_no_failed_statuses()`
