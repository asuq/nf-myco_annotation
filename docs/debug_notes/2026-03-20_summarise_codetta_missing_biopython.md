## `SUMMARISE_CODETTA` Missing Biopython

### Failure

The OIST HPC real pipeline run reached `PER_SAMPLE_ANNOTATION:SUMMARISE_CODETTA`
and failed with:

```text
ModuleNotFoundError: No module named 'Bio'
```

The failing container was `quay.io/asuq1617/python-scipy:3.12`.

### Root Cause

`bin/summarise_codetta.py` imports `Bio.Data.CodonTable`, so the shared Python
helper image must provide Biopython. The current helper image only installed
NumPy and SciPy, which was enough for ANI helpers but not for Codetta summary
parsing.

This was a helper-image contract gap rather than a Codetta workflow bug.

### Patch

- Added `biopython==1.84` to `docker/python_helper/Dockerfile`.
- Tightened `tests/test_python_helper_container.py` so the helper image must
  import Biopython as well as NumPy and SciPy.
- Updated `docs/runbook.md` so the shared helper-image requirements mention
  Biopython explicitly.
- Implementation commits:
  - `12d5612` `fix(helper): add biopython to python image`
  - `e148ab8` `test(helper): require biopython in python image`

### Verification

- `python3 -m unittest tests.test_python_helper_container tests.test_summarise_codetta`
- `export RUN_DOCKER_TESTS=1 && python3 -m unittest tests.test_python_helper_container`
- rebuild and publish `quay.io/asuq1617/python-scipy:3.12`
