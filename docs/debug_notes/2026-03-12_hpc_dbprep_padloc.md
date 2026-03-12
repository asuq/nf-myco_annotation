## HPC dbprep PADLOC failure chain

- Date: 2026-03-12
- Context: `bin/run_pipeline_test.sh dbprep-slurm --dbprep-profile oist`

### Failures seen

1. Prep helper image `quay.io/asuq1617/nf-myco_db:0.1` lacked `ps`, so Nextflow task metrics failed before taxdump prep.
2. PADLOC BioContainer tried to create `/usr/local/data` before honouring `--data`, which failed under read-only Singularity images.
3. The task-local PADLOC wrapper used shell literals such as `"$@"` and `${SRC_DIR}` inside Nextflow script strings without escaping, so Groovy attempted to interpolate them during module compilation.
4. After the compile-time fixes, the upstream PADLOC launcher still failed at runtime because it initialised `DATA=$(normpath "${SRC_DIR}/../data")` before processing `--data`, and the BioContainer does not ship that bundled `../data` directory.

### Fixes applied

- Bumped helper image tag to `quay.io/asuq1617/nf-myco_db:0.2`.
- Patched PADLOC prep and analysis modules to:
  - copy the real PADLOC launcher into the task directory;
  - rewrite the bootstrap data path to a writable task-local directory;
  - rewrite the launcher `DATA=` initialiser to the same task-local bootstrap directory;
  - execute the patched wrapper through a task-local shim.
- Escaped all shell literals in the wrapper block so Nextflow passes them through unchanged.
- Exposed `--force-runtime-database-rebuild` through `dbprep-slurm` so partial DB trees can be rebuilt in place after failed HPC attempts.
- Revalidated the PADLOC path locally with:
  - targeted syntax/config tests;
  - `prepare_databases.nf -profile test -stub-run`;
  - a focused Docker-backed PADLOC launcher contract test against the original BioContainer.

### Expected recovery command

```bash
bin/run_pipeline_test.sh dbprep-slurm \
  --resume \
  --force-runtime-database-rebuild \
  --dbprep-profile oist \
  --work-root /path/to/acceptance \
  --taxdump /path/to/db/ncbi_taxdump_20240914 \
  --checkm2-db /path/to/db/checkm2/CheckM2_database \
  --busco-db /path/to/db/busco \
  --eggnog-db /path/to/db/Eggnog_db/Eggnog_Diamond_db \
  --padloc-db /path/to/db/padloc \
  --singularity-cache-dir /path/to/singularity_cache
```
