## Codetta Container Missing `ps`

### Failure

The OIST HPC real pipeline run reached the `PER_SAMPLE_ANNOTATION:CODETTA`
process, but every Codetta task failed before Codetta itself ran. Nextflow
reported:

```text
Command 'ps' required by nextflow to collect task metrics cannot be found
```

The failing image was `quay.io/asuq1617/codetta:2.0`.

### Root Cause

The custom Codetta image is based on `python:3.12-slim`. That base image does
not include `ps`. Nextflow expects `ps` inside the task container so it can
collect process metrics while the task is running.

The Codetta image installed the Codetta runtime dependencies, but it did not
install `procps`, so the container contract was incomplete for Nextflow on HPC.

### Patch

- Added `procps` to `docker/codetta/Dockerfile`.
- Tightened `tests/test_codetta_container.py` so the runtime contract now
  requires `command -v ps`.
- Implementation commits:
  - `a9abf47` `fix(codetta): add procps to runtime image`
  - `62899b3` `test(codetta): require ps in container contract`

### Verification

- `python3 -m unittest tests.test_codetta_container`
- `export RUN_DOCKER_TESTS=1 && python3 -m unittest tests.test_codetta_container`
- rebuild the image and publish the refreshed `quay.io/asuq1617/codetta:2.0`
  tag before rerunning HPC
