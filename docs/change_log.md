# Change Log

## v0.3.0 - 2026-08-01

- Added an ANI recovery CLI with detailed logging, reusable published assembly
  statistics, safer manifest handling, and parallelised rescue workers.
- Added sequence-derived GC content, `is_new` reporting, optional inclusion of
  genomes with incomplete 16S calls in ANI clustering, and separate
  low-quality 16S cohort outputs.
- Added reusable `nf-helper` site configuration, GWDG SCC and MPCDF Viper CPU
  profiles, and bounded Viper Slurm submission and status polling behaviour.
- Updated the workflow for the strict Nextflow 26 syntax parser and set
  Nextflow 26 as the minimum supported release line.
- Improved BUSCO dataset preparation, lineage parsing, output compression, and
  summary parsing while avoiding concurrent writes to shared databases.
- Improved assembly-statistics failure reporting, temporary-file handling,
  parallel GC backfilling, and protection against pipeline race conditions.
- Corrected per-sample publication layout, retained Prokka GenBank output, and
  hardened final-output channel handling.
- Added a locked Pixi test environment and expanded configuration, integration,
  HPC, reporting, and recovery regression coverage.

This release preserves the software, container, database, and `nf-helper`
revisions already pinned by the latest commit on `main`.

## 2026-03-07

- Scaffolded the repository skeleton from `docs/design_spec.md`.
- Added placeholder Nextflow entrypoints, config files, docs, and assets.
