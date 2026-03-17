## 2026-03-17: fixed medium cohort used an invalid ANI candidate pair

### Failure
- The medium HPC run completed successfully, but post-run validation failed with
  `missing_ani_cluster_pair`.

### Root cause
- The fixed medium cohort validator requires at least two
  `ani_cluster_candidate` rows to share the same non-`NA` `Cluster_ID`.
- The tracked medium cohort plan defined:
  - `BACI-PAIR` from `GCF_000009045.1`
  - `BACI PAIR` from `GCF_000006885.1`
- Those sources are different organisms:
  - `GCF_000009045.1` is *Bacillus subtilis subsp. subtilis str. 168*
  - `GCF_000006885.1` is *Streptococcus pneumoniae TIGR4*
- They cannot satisfy the fixed ANI-pair contract at the pipeline's default ANI
  threshold, so this was a bad tracked cohort rather than a pipeline bug.

### Fix
- Keep the collision pair accessions unchanged:
  - `BACI-PAIR`
  - `BACI PAIR`
- Point both of those tracked ANI candidates at the same Bacillota source
  accession so the fixed medium cohort genuinely contains one reproducible ANI
  pair.
- Add a regression test that the tracked medium cohort's
  `ani_cluster_candidate` rows all come from one source accession.
