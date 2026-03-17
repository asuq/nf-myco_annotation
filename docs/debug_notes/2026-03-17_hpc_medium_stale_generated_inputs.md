## 2026-03-17: fixed medium cohort changes were hidden by stale generated inputs

### Failure
- The fixed medium ANI pair was corrected in the repo, but the same
  `missing_ani_cluster_pair` error still appeared on HPC.

### Root cause
- The OIST wrapper only regenerated the tracked medium inputs when
  `medium/generated/sample_sheet.csv` or `medium/generated/metadata.tsv` were
  missing.
- On an existing `HPC_ROOT`, `p2` therefore reused previously generated inputs
  that still pointed at the old invalid ANI pair, even though the repo-tracked
  cohort plan had already been fixed.

### Fix
- Keep operator-provided `--medium-sample-csv` and `--medium-metadata`
  overrides unchanged.
- For the tracked fixed medium cohort path, always rerun `medium-prepare`
  before `p2` and `all`.
- Preparation still reuses cached downloads, so this refresh is cheap but keeps
  the generated inputs aligned with the current tracked assets.
