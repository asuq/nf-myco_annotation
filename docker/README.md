# Docker Images

Module-specific Dockerfiles live under `docker/<tool>/Dockerfile` when
BioContainers do not cover a required runtime.

Current repo-owned images:

- `docker/ccfinder/Dockerfile`: CRISPRCasFinder runtime
- `docker/python_helper/Dockerfile`: shared Python helper runtime used by
  `params.python_container`, including the ANI scientific stack (`numpy` and
  `scipy`)
