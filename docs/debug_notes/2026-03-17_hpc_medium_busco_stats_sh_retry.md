## 2026-03-17: medium BUSCO run failed for one sample on `stats.sh --version`

### Failure
- The medium HPC run completed successfully, but post-run validation failed with:
  `Found failed sample-status rows: STREP_FILLER_1`
- The failing sample row showed:
  - `busco_bacillota_odb12_status = failed`
  - `warnings = busco_summary_failed`
  - `ani_exclusion_reason = missing_primary_busco`

### Root cause
- The published BUSCO artefacts for `STREP_FILLER_1` showed a transient BUSCO
  runtime failure before the real lineage analysis started.
- `busco.log` recorded:
  - `Command '['stats.sh', '--version']' returned non-zero exit status 134.`
- The sibling sample `STREP_FILLER_2` succeeded on the same lineage and source
  family, so this was not a deterministic sample-content failure.
- The BUSCO module was still running BUSCO only once, then converting any
  failure into an empty `short_summary.json`.
- That empty summary propagated into `busco_summary_failed` and a failed sample
  row, even though the run would likely have succeeded on retry.

### Fix
- Add an internal retry loop to the BUSCO module using the shared
  `params.task_attempts` total-attempt policy.
- Keep the existing degraded fallback after the final attempt:
  - preserve `busco.log`
  - emit an empty `short_summary.json`
  - let downstream summarisation record the sample-level failure without killing
    the whole pipeline
- This mirrors the existing hardening already applied to CheckM2.
