# HPC eggNOG DB prep read-only checkout fix

## Failure

`DOWNLOAD_EGGNOG_DATABASE` failed on OIST while trying to write a patched
copy of the bundled downloader into the software checkout.

The software checkout is readable on compute nodes but not writable.

## Cause

The original BioContainer ships `download_eggnog_data.py` with stale
`eggnogdb.embl.de` download URLs. The workflow was patching that script at
runtime, and an earlier interpolation bug placed the patched copy under the
read-only project checkout instead of the task work directory.

## Fix

The final fix is to use a custom image with the downloader patched at build
time. That removes the need to write any patched helper script during runtime,
so the software checkout can remain read-only on compute nodes.

## Local validation

- `python3 -m unittest tests.test_nextflow_module_syntax tests.test_nextflow_config_contracts tests.test_eggnog_container`
- `nextflow run prepare_databases.nf -profile test -stub-run --taxdump /tmp/taxdump --checkm2_db /tmp/checkm2 --busco_db /tmp/busco --eggnog_db /tmp/eggnog --outdir /tmp/runtime_db_prep_clean_env`

The fixed image contract and the prep workflow stub both passed locally.
