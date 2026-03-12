# HPC eggNOG DB prep read-only checkout fix

## Failure

`DOWNLOAD_EGGNOG_DATABASE` failed on OIST while trying to write:

- `/bucket/.../nf-myco_annotation/download_eggnog_data.py.patched`

The software checkout is readable on compute nodes but not writable.

## Cause

The eggNOG prep module used `${PWD}` inside a Groovy-interpolated Nextflow script
block:

- `patched_script="${PWD}/download_eggnog_data.py.patched"`

That resolved outside the task work directory, so the module tried to write the
patched downloader into the project checkout.

## Fix

Write the patched helper into the task work directory explicitly with shell PWD:

- `patched_script="$PWD/download_eggnog_data.py.patched"`

## Local validation

- `python3 -m unittest tests.test_nextflow_module_syntax tests.test_nextflow_config_contracts`
- `nextflow run prepare_databases.nf -profile test -stub-run --taxdump /tmp/taxdump --checkm2_db /tmp/checkm2 --busco_db /tmp/busco --eggnog_db /tmp/eggnog --outdir /tmp/runtime_db_prep_readonly_fix`

Both passed locally.
