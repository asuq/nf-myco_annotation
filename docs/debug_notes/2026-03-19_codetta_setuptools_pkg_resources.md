## 2026-03-19: Codetta container `pkg_resources` compatibility fix

### Failure
- The first Docker build of `docker/codetta/Dockerfile` completed, but the log
  showed `ModuleNotFoundError: No module named 'pkg_resources'` during
  `setup.sh`.
- Manual inspection of the built image confirmed the problem:
  `python3 -c "import pkg_resources"` failed inside the container.
- The build appeared successful only because the upstream `setup.sh` did not
  fail the shell after the internal Python error.

### Cause
- Current Setuptools releases have removed `pkg_resources`.
- Codetta `v2.0` still calls `pkg_resources` through its setup helper, so the
  image needs an older compatible Setuptools release.

### Fix
- Change the image dependency install step to:
  `python3 -m pip install --no-cache-dir numpy scipy "setuptools<82"`.
- Add an explicit `python3 -c "import pkg_resources"` check in the Dockerfile
  so the build fails immediately if the compatibility shim is missing.
- Extend the Codetta container contract test to assert that
  `pkg_resources` imports successfully at runtime.

### Result
- Rebuilding `quay.io/asuq1617/codetta:2.0` now completes with
  `Codetta setup complete!`.
- The runtime contract passes and the image exposes:
  - Codetta scripts on `PATH`
  - the bundled `resources/` directory
  - metadata files for `v2.0`, pinned source commit, and source branch
