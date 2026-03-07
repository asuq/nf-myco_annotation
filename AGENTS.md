# AGENTS.md

- Read `docs/design_spec.md` before proposing or making pipeline changes. It is the source of truth for v1 behaviour, and it wins if this file is less specific.
- Build and maintain the pipeline in Nextflow DSL2, targeting Nextflow `>= 25.04.8`.
- Keep orchestration in Nextflow. Put parsing, validation, summarisation, and merge logic in small CLIs under `bin/`.
- Preserve original metadata columns exactly as provided and in original order. Append derived columns only after the metadata block.
- Keep the final reporting contracts centered on `master_table.csv`, `sample_status.tsv`, and `tool_and_db_versions.tsv`.
- Use original accession values in the master table and in published output folder names. Use sanitized internal IDs only for execution-time filenames, prefixes, and other tool-safe internals.
- Prefer BioContainers. Add `docker/<tool>/Dockerfile` only when no usable container is available.
- Use nf-core-like retry and resource-handling patterns in Nextflow configuration and process definitions.
- Do not add PADLOC or eggNOG outputs to the master table unless `docs/design_spec.md` changes.
- Keep explicit Nextflow module interfaces with named outputs, and include `stub` sections for testability.
