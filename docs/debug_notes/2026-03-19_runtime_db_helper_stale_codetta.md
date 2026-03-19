## 2026-03-19: runtime-db helper image skew and HPC Codetta wrapper gap

- Context: local `prepare_databases.nf` acceptance with `--taxdump` and
  `--codetta_db`, plus the OIST HPC dbprep and matrix wrappers.

### Failures seen

1. `prepare_databases.nf` failed in `PREP_CODETTA_DATABASE` because the
   published helper image `quay.io/asuq1617/nf-myco_db:0.2` bundled an older
   `prepare_runtime_databases.py` that did not recognise `--codetta-dest`.
2. The runtime-db Nextflow modules hardcoded helper script paths under
   `/usr/local/bin`, which made image/code skew harder to spot and recover
   from.
3. `bin/run_oist_hpc_matrix.sh` still treated Codetta as optional. Its dbprep
   wrapper, validator expectations, ready-root checks, real pipeline runs, and
   disposable DB matrix cases all omitted `--codetta_db`.
4. `bin/validate_hpc_matrix.py` could not validate a Codetta dbprep component
   or `--codetta_db` in `nextflow_args.txt`.

### Root cause

- Repo source already supported Codetta runtime database preparation, but the
  published helper image tag lagged behind the repo.
- The HPC wrapper layer was updated for taxdump, CheckM2, BUSCO, and eggNOG,
  but Codetta was never threaded through the parallel dbprep validation path.

### Fixes applied

- Moved runtime-db prep and finalise modules to resolve helper CLIs and
  `python3` from `PATH`, instead of assuming one fixed absolute image path.
- Bumped the dedicated helper image contract in `conf/base.config` to
  `quay.io/asuq1617/nf-myco_db:0.3`.
- Added a wrapper preflight in `bin/run_pipeline_test.sh` so `dbprep-slurm`
  and `all` fail fast if `conf/base.config` still pins the stale helper tag.
- Threaded Codetta through `bin/run_oist_hpc_matrix.sh`:
  - canonical HPC DB root `db/codetta/Pfam-A_enone`
  - ready-root checks
  - `db-full` and `db-reuse` validator expectations
  - real pipeline runs
  - disposable db-matrix cases for valid, missing, file, and partial states
- Extended `bin/validate_hpc_matrix.py` to accept `--codetta-db`, validate the
  `codetta` component, and treat `--codetta_db` as an expected dbprep flag.
- Strengthened the helper image container test so the built image must expose
  the Codetta dbprep flags and a manifest containing the `codetta` component.

### Recovery path

1. Rebuild the helper image from `docker/runtime_db_helper/Dockerfile` as
   `quay.io/asuq1617/nf-myco_db:0.3`.
2. Re-run the local `prepare_databases.nf` acceptance path to confirm Codetta
   and taxdump prepare successfully in one run.
3. Re-run the OIST HPC gates:
   - `bin/run_pipeline_test.sh dbprep-slurm ...`
   - `bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root db-full`
   - `bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root db-matrix`
