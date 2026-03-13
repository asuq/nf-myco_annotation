## HPC dbprep MD5 mismatch retry

- Date: 2026-03-13
- Context: OIST `dbprep-slurm` campaign
- Symptom:
  - repeated taxdump download failures with:
    - `Checksum mismatch for new_taxdump.tar.gz: expected ... observed ...`
- Cause:
  - the shared runtime DB helper treated one checksum mismatch as a hard failure
  - transient remote corruption or incomplete transfer could therefore fail the whole prep run immediately
- Fix:
  - add a bounded retry path in `prepare_runtime_databases.py`
  - retry only MD5 mismatches
  - maximum attempts: 3 total downloads
  - keep all other checksum failures as immediate errors
- Validation:
  - unit test proving persistent mismatch fails after 3 attempts
  - unit test proving a later successful download is accepted
  - full Python suite passed after the change
