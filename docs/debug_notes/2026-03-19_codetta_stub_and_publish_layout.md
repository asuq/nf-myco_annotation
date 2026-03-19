## 2026-03-19: Codetta stub output and publish-layout fixes

### Background
- Codetta was added as a new always-on per-sample annotation step.
- The change touched runtime database preparation, per-sample execution,
  final reporting, and stub-mode acceptance coverage.

### Failures
- The first local `main.nf -stub-run` failed because
  `SUMMARISE_CODETTA` declared `versions.yml` but its stub block did not
  create that file.
- The first successful stub run still produced an incomplete
  `tool_and_db_versions.tsv` because the `COLLECT_VERSIONS` stub did not emit
  Codetta tool, container, or database rows.
- The published sample layout placed raw Codetta files under
  `samples/<accession>/codetta/codetta/` because `CODETTA` published the
  output directory directly into an already tool-named target directory.

### Fix
- Add `versions.yml` generation to the `SUMMARISE_CODETTA` stub.
- Extend the `COLLECT_VERSIONS` stub output with `codetta`, `codetta`
  container, and `codetta_db` rows.
- Change the `CODETTA` `publishDir` target from
  `samples/<accession>/codetta` to `samples/<accession>` so the emitted
  `codetta/` directory lands exactly once.
- Add a regression assertion in `tests/test_nextflow_module_syntax.py` for the
  `CODETTA` publish target.

### Result
- `main.nf -profile test -stub-run` completes successfully.
- `prepare_databases.nf -profile test -stub-run` completes successfully with
  Codetta database preparation wired in.
- The final stub outputs now include Codetta columns and statuses in the
  reporting tables, Codetta provenance in `tool_and_db_versions.tsv`, and the
  expected sample-level publish layout under `samples/<accession>/codetta/`.
