## 2026-03-17: medium validator rejected the deliberate missing-metadata sample

### Failure
- The medium HPC run completed successfully, but post-run validation failed with:
  `Observed disallowed phyla in medium run: NA`
- The offending row was `MYCO_NEW_1`.

### Root cause
- `MYCO_NEW_1` is the fixed medium cohort's deliberate `missing_metadata_case`.
- It is prepared with `is_new=true` and `include_metadata=false`, so the generated
  metadata table intentionally omits its metadata row.
- The pipeline therefore emits `Tax_ID = NA` and taxonomy columns, including
  `phylum`, as `NA`.
- This matches the design spec and the acceptance-harness role checks.
- The medium validator was too strict and treated every `phylum = NA` value as a
  disallowed phylum.

### Fix
- Keep the pipeline and fixed medium cohort unchanged.
- Allow `phylum = NA` only for accessions tagged `missing_metadata_case` in the
  fixed medium cohort plan.
- Keep rejecting any unexpected `phylum = NA` rows and any non-allowed,
  non-`NA` phyla.
- Make the fixed medium wrapper path pass the medium metadata, cohort plan, and
  source catalog into the validator so the stricter tracked-cohort contract can
  be applied.
