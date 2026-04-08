# AGENTS.md

- Read `docs/design_spec.md` before proposing or making pipeline changes. It is the source of truth for v1 behaviour, and it wins if this file is less specific.
- Build and maintain the pipeline in Nextflow DSL2, targeting Nextflow `>= 25.04.8`.
- Keep orchestration in Nextflow. Put parsing, validation, summarisation, and merge logic in small CLIs under `bin/`.
- Preserve original metadata columns exactly as provided and in original order. Append derived columns only after the metadata block.
- Keep the final reporting contracts centered on the emitted master table, `sample_status.tsv`, and `tool_and_db_versions.tsv`.
- Use original accession values in the master table and in published output folder names. Use sanitized internal IDs only for execution-time filenames, prefixes, and other tool-safe internals.
- Prefer BioContainers. Custom process images are added manually; when no usable container is available, add only a clearly documented placeholder under `docker/<tool>/` unless explicitly asked to implement the image.
- Use nf-core-like retry and resource-handling patterns in Nextflow configuration and process definitions.
- Do not add PADLOC or eggNOG outputs to the master table unless `docs/design_spec.md` changes.
- Keep explicit Nextflow module interfaces with named outputs, and include `stub` sections for testability.
- Keep tracked project documentation in `docs/`.
- Keep all local working docs, debug notes, commit notes, and other repo-specific notes under `.untracked/`; do not use `/tmp` for repo-specific notes.
- Classify local notes under `.untracked/` by purpose:
  - `.untracked/debug/` for incident-specific debugging and recovery notes
  - `.untracked/ledger/` for current bug and state ledgers
  - `.untracked/notes/` for maintenance, planning, and change-history notes
- Keep the `.untracked/` top level limited to bucket directories plus a small README.
