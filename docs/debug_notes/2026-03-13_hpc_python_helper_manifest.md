## HPC python-scipy helper image manifest mismatch

### Failure
- OIST pipeline runs failed before `VALIDATE_INPUTS` could start.
- Singularity could not pull `docker://quay.io/asuq1617/python-scipy:3.12`.
- Error:
  - `no image satisfies requested platform: linux/amd64`

### Cause
- The repo already pointed `params.python_container` at `quay.io/asuq1617/python-scipy:3.12`.
- The local `docker/python_helper/Dockerfile` built and ran correctly for `linux/amd64`.
- The published tag on Quay was therefore out of sync with the tested local image and did not expose a usable amd64 manifest to Singularity on OIST.

### Fix
- Verify the helper image locally with:
  - `docker build --platform linux/amd64 -t codex-python-scipy-test:3.12 docker/python_helper`
  - `docker run --rm --platform linux/amd64 codex-python-scipy-test:3.12 python -c "import numpy, scipy; ..."`
- Republish `quay.io/asuq1617/python-scipy:3.12` with an explicit `linux/amd64` manifest.

### Follow-up
- Add a local container contract test for the shared Python helper image so the amd64 scientific stack is checked before HPC runs.
