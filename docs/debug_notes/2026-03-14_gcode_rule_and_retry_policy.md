## 2026-03-14: gcode rule selection and retry policy

### Failure
- The medium HPC run completed, but `sample_status.tsv` still contained failed
  rows.
- Two causes were confirmed:
  - close valid CheckM2 pairs in the fixed medium cohort produced `Gcode = NA`
    under the strict delta rule
  - one `CHECKM2_GCODE11` task for `GCA_965226525.1_MED1` hit a transient
    `OSError: [Errno 5] Input/output error` inside the CheckM2 container

### Decision
- Keep the original gcode resolution as rule `strict_delta`.
- Add rule `delta_then_11`, which applies the strict delta thresholds first and
  then falls back to `11` for valid close-score pairs.
- Expose the rule through `params.gcode_rule`, defaulting to `strict_delta`.
- Normal runtime tasks should make 3 total attempts.
- Database download tasks should make 2 total attempts.
- CheckM2 keeps its degraded sample-level behaviour, but now retries the tool
  call internally before emitting an empty summary report.

### Patch direction
- Add explicit `task_attempts` and `db_download_attempts` config params.
- Remove most hard-coded `maxRetries = 0` runtime overrides.
- Thread `gcode_rule` through the HPC wrapper and acceptance harness.
- Update `summarise_checkm2.py` and the relevant regression tests.
