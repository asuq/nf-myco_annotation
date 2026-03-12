## HPC dbprep PADLOC failure chain

- Date: 2026-03-12
- Context: `bin/run_pipeline_test.sh dbprep-slurm --dbprep-profile oist`

### Failures seen

1. Prep helper image `quay.io/asuq1617/nf-myco_db:0.1` lacked `ps`, so Nextflow task metrics failed before taxdump prep.
2. PADLOC BioContainer tried to create `/usr/local/data` before honouring `--data`, which failed under read-only Singularity images.
3. The task-local PADLOC wrapper used shell literals such as `"$@"` and `${SRC_DIR}` inside Nextflow script strings without escaping, so Groovy attempted to interpolate them during module compilation.
4. After the compile-time fixes, the upstream PADLOC launcher still failed at runtime because it initialised `DATA=$(normpath "${SRC_DIR}/../data")` before processing `--data`, and the BioContainer does not ship that bundled `../data` directory.

### Final fix

- Bumped helper image tag to `quay.io/asuq1617/nf-myco_db:0.2`.
- Replaced PADLOC runtime patching with a fixed custom image:
  - `quay.io/asuq1617/padloc:2.0.0-nfmyco1`
  - based on `quay.io/biocontainers/padloc:2.0.0--hdfd78af_1`
  - launcher patched at image build time so it no longer initialises the missing bundled `../data` path.
- Exposed `--force-runtime-database-rebuild` through `dbprep-slurm` so partial DB trees can be rebuilt in place after failed HPC attempts.

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
