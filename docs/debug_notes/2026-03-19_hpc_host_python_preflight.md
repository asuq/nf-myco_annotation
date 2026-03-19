## HPC Host Python Preflight

### Failure

The OIST HPC wrapper failed immediately on the login node with:

```text
SyntaxError: future feature annotations is not defined
```

The reported host interpreter was `Python 3.6.8`.

### Root Cause

The acceptance harness and HPC validators are invoked directly through
`python3` on the host shell before Nextflow submits work to SLURM. This repo
now assumes Python `3.12` for host-side helper scripts, but the wrappers did
not preflight the interpreter version. As a result, an old login-node
interpreter surfaced as a raw Python syntax error instead of a clear operator
message.

### Patch

- Implementation commits:
  - `a57ab74` `fix(hpc): preflight host python in operator wrappers`
  - `1327d11` `test(hpc): cover old host python wrapper failures`
- Added host `python3 >= 3.12` preflights to:
  - [bin/run_pipeline_test.sh](/Users/asuq/Documents/Lab/Coding/nf-myco_annotation/bin/run_pipeline_test.sh)
  - [bin/run_oist_hpc_matrix.sh](/Users/asuq/Documents/Lab/Coding/nf-myco_annotation/bin/run_oist_hpc_matrix.sh)
- Added regression tests that simulate `Python 3.6.8` on `PATH`.
- Updated the HPC-facing docs to:
  - call out the host Python requirement explicitly
  - recommend `prepare -> dbprep-slurm -> p1 -> all`
  - show the missing Codetta paths in the HPC runbook examples

### Expected Operator Outcome

When the login node exposes an older interpreter, the wrappers should now stop
immediately with a direct message telling the operator to load Python `3.12` or
adjust `PATH`.
