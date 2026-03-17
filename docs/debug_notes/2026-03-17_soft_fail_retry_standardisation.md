## 2026-03-17: standardise internal retries for soft-fail runtime modules

### Background
- The pipeline already had global Nextflow retry control through
  `params.task_attempts` and `conf/base.config`.
- Several sample-level modules intentionally swallowed tool failures and
  emitted placeholder outputs so the whole pipeline could continue.
- Because those modules exited successfully after degrading their outputs,
  Nextflow never saw a failed task and its retry logic could not help.

### Problem
- BUSCO was hardened first after an HPC failure in which one sample hit a
  transient `stats.sh --version` crash.
- That fix worked, but it left BUSCO as a one-off special case even though
  Barrnap, CheckM2, Prokka, CRISPRCasFinder, PADLOC, and eggNOG had the same
  soft-fail structure.

### Fix
- Keep Nextflow retry control unchanged.
- Add a separate `params.soft_fail_attempts` knob for internal retries in the
  soft-fail sample-level modules.
- Standardise the same retry-before-fallback pattern across:
  - Barrnap
  - CheckM2
  - BUSCO
  - Prokka
  - CRISPRCasFinder
  - PADLOC
  - eggNOG

### Result
- Transient tool errors now get multiple attempts even when the module keeps
  the pipeline alive by degrading outputs on final failure.
- Hard-fail tasks such as FastANI and database download tasks remain on their
  existing retry paths.
